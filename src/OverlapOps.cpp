#include "OverlapOps.hpp"

#include <pdal/Dimension.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

namespace overlap
{
namespace
{
struct SwathBucket
{
    SwathKeyMode mode;
    std::uint16_t psid = 0;
    std::uint8_t channel = 0;
    std::vector<pdal::PointId> points;
    std::vector<double> weights;
};

SwathBucket* findBucket(std::vector<SwathBucket>& buckets,
                        SwathKeyMode mode,
                        std::uint16_t psid,
                        std::uint8_t channel)
{
    for (auto& b : buckets)
    {
        if (b.mode != mode)
            continue;
        if (b.psid == psid && (mode != SwathKeyMode::PointSourceIdChannel || b.channel == channel))
            return &b;
    }
    return nullptr;
}

struct Candidate
{
    pdal::PointId id;
    double absAngle;
    double angle;
    int returnNumber;
    double intensity;
    double gpsTime;
    std::uint16_t psid;
    std::uint8_t channel;
};

bool isBetter(const Candidate& a, const Candidate& b)
{
    if (a.absAngle != b.absAngle)
        return a.absAngle < b.absAngle;
    if (a.returnNumber != b.returnNumber)
        return a.returnNumber < b.returnNumber;
    if (a.intensity != b.intensity)
        return a.intensity > b.intensity;
    if (a.gpsTime != b.gpsTime)
        return a.gpsTime < b.gpsTime;
    return a.id < b.id;
}

} // namespace

pdal::PointViewPtr equalizeOverlapOnly(pdal::PointViewPtr view,
                                       const OverlapMask& mask,
                                       double targetDensity,
                                       std::uint64_t seed,
                                       EqualizeOverlapStats* stats)
{
    if (!view)
        throw std::invalid_argument("equalizeOverlapOnly received null PointViewPtr");
    if (mask.pointMask.size() != view->size())
        throw std::runtime_error("Overlap mask does not match point count");
    if (targetDensity < 0.0)
        throw std::invalid_argument("Target density must be non-negative");

    auto output = view->makeNew();
    if (!output)
        throw std::runtime_error("Failed to allocate output PointView");

    if (view->size() == 0)
        return output;

    const auto dimPsid = pdal::Dimension::Id::PointSourceId;
    const auto dimScanner = pdal::Dimension::Id::ScanChannel;

    auto layout = view->layout();
    const bool hasScanner = layout && layout->hasDim(dimScanner);

    const std::size_t cellCount = mask.nx * mask.ny;
    std::vector<std::vector<pdal::PointId>> flaggedPoints(cellCount);
    std::vector<char> keep(view->size(), 0);

    std::size_t overlapPoints = 0;

    for (pdal::PointId pid = 0; pid < view->size(); ++pid)
    {
        if (mask.pointMask[pid])
        {
            ++overlapPoints;
            const std::int32_t storedId = mask.cellIndexPerPoint[pid];
            if (storedId < 0)
                throw std::runtime_error("Overlap mask contains negative cell index");
            const auto cellId = static_cast<std::size_t>(storedId);
            if (cellId >= flaggedPoints.size())
                throw std::runtime_error("Overlap mask cell index out of range");
            flaggedPoints[cellId].push_back(pid);
        }
        else
        {
            keep[pid] = 1;
        }
    }

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    const double cellArea = mask.cellSize * mask.cellSize;

    std::size_t keptOverlap = 0;

    for (std::size_t cellId = 0; cellId < flaggedPoints.size(); ++cellId)
    {
        const auto& ids = flaggedPoints[cellId];
        if (ids.empty())
            continue;

        std::vector<SwathBucket> buckets;
        buckets.reserve(ids.size());
        std::size_t total = ids.size();

        for (pdal::PointId pid : ids)
        {
            const auto psid = view->getFieldAs<std::uint16_t>(dimPsid, pid);
            const auto channel = hasScanner ? view->getFieldAs<std::uint8_t>(dimScanner, pid)
                                            : static_cast<std::uint8_t>(0);
            SwathBucket* bucket = findBucket(buckets, mask.swathKeyUsed, psid, channel);
            if (!bucket)
            {
                SwathBucket newBucket;
                newBucket.mode = mask.swathKeyUsed;
                newBucket.psid = psid;
                newBucket.channel = channel;
                buckets.push_back(newBucket);
                bucket = &buckets.back();
            }
            bucket->points.push_back(pid);
            bucket->weights.push_back(dist(rng));
        }

        std::size_t cellKept = 0;
        const double desiredTotal = targetDensity * cellArea;

        for (auto& bucket : buckets)
        {
            const std::size_t count = bucket.points.size();
            if (count == 0)
                continue;

            double expected = desiredTotal * (static_cast<double>(count) / static_cast<double>(total));
            if (!std::isfinite(expected))
                expected = 0.0;

            std::size_t quota = static_cast<std::size_t>(std::llround(expected));
            quota = std::min<std::size_t>(quota, count);

            if (quota >= count)
            {
                for (pdal::PointId pid : bucket.points)
                {
                    keep[pid] = 1;
                    ++cellKept;
                }
                continue;
            }

            if (quota == 0)
                continue;

            std::vector<std::size_t> indices(count);
            std::iota(indices.begin(), indices.end(), 0);
            std::nth_element(indices.begin(), indices.begin() + static_cast<std::ptrdiff_t>(quota), indices.end(),
                             [&](std::size_t a, std::size_t b) {
                                 return bucket.weights[a] < bucket.weights[b];
                             });
            for (std::size_t i = 0; i < quota; ++i)
            {
                keep[bucket.points[indices[i]]] = 1;
                ++cellKept;
            }
        }
        keptOverlap += cellKept;
    }

    for (pdal::PointId pid = 0; pid < view->size(); ++pid)
    {
        if (keep[pid])
            output->appendPoint(*view, pid);
    }

    if (stats)
    {
        stats->overlapPoints = overlapPoints;
        stats->keptOverlapPoints = keptOverlap;
        stats->droppedOverlapPoints = overlapPoints >= keptOverlap ? overlapPoints - keptOverlap : 0;
    }

    return output;
}

pdal::PointViewPtr joinOverlapByScanAngle(pdal::PointViewPtr view,
                                          const OverlapMask& mask,
                                          JoinOverlapStats* stats)
{
    if (!view)
        throw std::invalid_argument("joinOverlapByScanAngle received null PointViewPtr");
    if (mask.pointMask.size() != view->size())
        throw std::runtime_error("Overlap mask does not match point count");

    auto output = view->makeNew();
    if (!output)
        throw std::runtime_error("Failed to allocate output PointView");

    auto layout = view->layout();
    const bool hasScanAngleRank = layout && layout->hasDim(pdal::Dimension::Id::ScanAngleRank);
    const bool hasReturnNumber = layout && layout->hasDim(pdal::Dimension::Id::ReturnNumber);
    const bool hasIntensity = layout && layout->hasDim(pdal::Dimension::Id::Intensity);
    const bool hasGps = layout && layout->hasDim(pdal::Dimension::Id::GpsTime);
    const bool hasScanner = layout && layout->hasDim(pdal::Dimension::Id::ScanChannel);

    const auto dimPsid = pdal::Dimension::Id::PointSourceId;
    const auto dimScanner = pdal::Dimension::Id::ScanChannel;
    std::vector<std::vector<pdal::PointId>> flaggedPoints(mask.nx * mask.ny);
    std::vector<char> keep(view->size(), 0);

    std::size_t overlapPoints = 0;

    for (pdal::PointId pid = 0; pid < view->size(); ++pid)
    {
        if (mask.pointMask[pid])
        {
            ++overlapPoints;
            const std::int32_t storedId = mask.cellIndexPerPoint[pid];
            if (storedId < 0)
                throw std::runtime_error("Overlap mask contains negative cell index");
            const auto cellId = static_cast<std::size_t>(storedId);
            if (cellId >= flaggedPoints.size())
                throw std::runtime_error("Overlap mask cell index out of range");
            flaggedPoints[cellId].push_back(pid);
        }
        else
        {
            keep[pid] = 1;
        }
    }

    std::size_t keptOverlap = 0;

    for (std::size_t cellId = 0; cellId < flaggedPoints.size(); ++cellId)
    {
        auto& ids = flaggedPoints[cellId];
        if (ids.empty())
            continue;

        Candidate best{};
        bool hasBest = false;

        for (pdal::PointId pid : ids)
        {
            double angle = 0.0;
            if (hasScanAngleRank)
                angle = static_cast<double>(view->getFieldAs<std::int16_t>(pdal::Dimension::Id::ScanAngleRank, pid));

            Candidate cand;
            cand.id = pid;
            cand.angle = angle;
            cand.absAngle = std::abs(angle);
            cand.returnNumber = hasReturnNumber ? view->getFieldAs<std::uint8_t>(pdal::Dimension::Id::ReturnNumber, pid)
                                                : std::numeric_limits<int>::max();
            cand.intensity = hasIntensity ? static_cast<double>(view->getFieldAs<std::uint16_t>(pdal::Dimension::Id::Intensity, pid))
                                          : 0.0;
            cand.gpsTime = hasGps ? view->getFieldAs<double>(pdal::Dimension::Id::GpsTime, pid) : static_cast<double>(pid);
            cand.psid = view->getFieldAs<std::uint16_t>(dimPsid, pid);
            cand.channel = hasScanner ? view->getFieldAs<std::uint8_t>(dimScanner, pid)
                                      : static_cast<std::uint8_t>(0);

            if (!hasBest || isBetter(cand, best))
            {
                best = cand;
                hasBest = true;
            }
        }

        if (hasBest)
        {
            auto matchesBest = [&](pdal::PointId pid) {
                if (view->getFieldAs<std::uint16_t>(dimPsid, pid) != best.psid)
                    return false;
                if (mask.swathKeyUsed == SwathKeyMode::PointSourceIdChannel)
                {
                    const std::uint8_t ch = hasScanner
                                                ? view->getFieldAs<std::uint8_t>(dimScanner, pid)
                                                : static_cast<std::uint8_t>(0);
                    return ch == best.channel;
                }
                return true;
            };

            for (pdal::PointId pid : ids)
            {
                if (matchesBest(pid))
                {
                    if (!keep[pid])
                    {
                        keep[pid] = 1;
                        ++keptOverlap;
                    }
                }
            }
        }
    }

    for (pdal::PointId pid = 0; pid < view->size(); ++pid)
    {
        if (keep[pid])
            output->appendPoint(*view, pid);
    }

    if (stats)
    {
        stats->overlapPoints = overlapPoints;
        stats->keptOverlapPoints = keptOverlap;
        stats->droppedOverlapPoints = overlapPoints >= keptOverlap ? overlapPoints - keptOverlap : 0;
    }

    return output;
}

} // namespace overlap
