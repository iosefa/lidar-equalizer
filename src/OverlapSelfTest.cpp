#include "OverlapSelfTest.hpp"

#include "OverlapFlagger.hpp"
#include "OverlapOps.hpp"

#include <pdal/Dimension.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>

#include <cmath>
#include <iostream>
#include <vector>

namespace overlap
{
namespace
{
struct TestPoint
{
    double x;
    double y;
    double z;
    std::uint16_t psid;
    std::uint8_t channel;
    double scanAngle;
    std::uint8_t returnNumber;
    std::uint16_t intensity;
    double gpsTime;
};

std::size_t countMarked(const OverlapMask& mask)
{
    std::size_t total = 0;
    for (auto v : mask.pointMask)
        total += v ? 1 : 0;
    return total;
}

} // namespace

bool runOverlapSelfTest()
{
    bool passed = true;

    pdal::PointTable table;
    auto layout = table.layout();
    layout->registerDim(pdal::Dimension::Id::X);
    layout->registerDim(pdal::Dimension::Id::Y);
    layout->registerDim(pdal::Dimension::Id::Z);
    layout->registerDim(pdal::Dimension::Id::PointSourceId);
    layout->registerDim(pdal::Dimension::Id::ScanChannel);
    layout->registerDim(pdal::Dimension::Id::ScanAngleRank);
    layout->registerDim(pdal::Dimension::Id::Intensity);
    layout->registerDim(pdal::Dimension::Id::ReturnNumber);
    layout->registerDim(pdal::Dimension::Id::GpsTime);
    layout->registerDim(pdal::Dimension::Id::Classification);
    table.finalize();

    auto view = std::make_shared<pdal::PointView>(table);

    const TestPoint pts[] = {
        {0.20, 0.20, 10.0, 1, 0, -4.0, 1, 80, 10.0},
        {0.30, 0.20, 10.1, 2, 1, 8.0, 1, 90, 11.0},
        {0.35, 0.25, 10.2, 2, 1, 2.0, 2, 120, 9.5},
        {1.50, 0.20, 9.5, 1, 0, 5.0, 1, 70, 12.0},
        {0.25, 1.20, 9.8, 3, 0, -1.0, 1, 60, 13.0},
    };

    for (std::size_t i = 0; i < std::size(pts); ++i)
    {
        const TestPoint& p = pts[i];
        const auto pid = static_cast<pdal::PointId>(i);
        view->setField(pdal::Dimension::Id::X, pid, p.x);
        view->setField(pdal::Dimension::Id::Y, pid, p.y);
        view->setField(pdal::Dimension::Id::Z, pid, p.z);
        view->setField(pdal::Dimension::Id::PointSourceId, pid, p.psid);
        view->setField(pdal::Dimension::Id::ScanChannel, pid, p.channel);
        view->setField(pdal::Dimension::Id::ScanAngleRank, pid, static_cast<std::int16_t>(p.scanAngle));
        view->setField(pdal::Dimension::Id::Intensity, pid, p.intensity);
        view->setField(pdal::Dimension::Id::ReturnNumber, pid, p.returnNumber);
        view->setField(pdal::Dimension::Id::GpsTime, pid, p.gpsTime);
        view->setField(pdal::Dimension::Id::Classification, pid, static_cast<std::uint8_t>(1));
    }

    std::cout << "[self-test] generated synthetic view with " << view->size() << " points.\n";
    OverlapParams params;
    params.cellSize = 1.0;
    params.minPointsPerSwath = 1;
    params.swathKey = SwathKeyMode::PointSourceIdChannel;

    OverlapStats stats;
    auto mask = buildOverlapMask(*view, params, &stats);
    if (stats.overlapCells != 1 || countMarked(mask) != 3)
    {
        std::cerr << "[self-test] Unexpected overlap mask stats: cells="
                  << stats.overlapCells << " marked=" << countMarked(mask) << '\n';
        passed = false;
    }

    EqualizeOverlapStats eqStats;
    auto thinned = equalizeOverlapOnly(view, mask, 1.0, 42, &eqStats);
    if (eqStats.overlapPoints != 3 || eqStats.keptOverlapPoints >= eqStats.overlapPoints)
    {
        std::cerr << "[self-test] Overlap equalization stats off: overlap="
                  << eqStats.overlapPoints << " kept=" << eqStats.keptOverlapPoints << '\n';
        passed = false;
    }
    if (thinned && thinned->size() >= view->size())
    {
        std::cerr << "[self-test] Expected equalizeOverlapOnly to drop some points.\n";
        passed = false;
    }

    JoinOverlapStats joinStats;
    auto joined = joinOverlapByScanAngle(view, mask, &joinStats);
    if (joinStats.overlapPoints != 3 || joinStats.keptOverlapPoints != 2)
    {
        std::cerr << "[self-test] Overlap join stats off: overlap="
                  << joinStats.overlapPoints << " kept=" << joinStats.keptOverlapPoints << '\n';
        passed = false;
    }
    if (!joined)
    {
        std::cerr << "[self-test] joinOverlapByScanAngle returned null view.\n";
        passed = false;
    }
    else
    {
        const auto dimScanAngle = pdal::Dimension::Id::ScanAngleRank;
        bool foundPreferred = false;
        std::size_t totalKept = 0;
        for (pdal::PointId pid = 0; pid < joined->size(); ++pid)
        {
            const double angle = joined->getFieldAs<double>(dimScanAngle, pid);
            if (std::abs(angle - 2.0) < 1e-6)
            {
                foundPreferred = true;
            }
            ++totalKept;
        }
        if (!foundPreferred)
        {
            std::cerr << "[self-test] Expected joiner to keep point with scan angle 2.0.\n";
            passed = false;
        }
        if (totalKept != 4)
        {
            std::cerr << "[self-test] Expected joined view to contain 4 points, got "
                      << totalKept << '\n';
            passed = false;
        }
    }

    if (passed)
        std::cout << "[self-test] Overlap features OK.\n";

    return passed;
}

} // namespace overlap
