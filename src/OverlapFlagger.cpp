#include "OverlapFlagger.hpp"

#include <pdal/Dimension.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

namespace overlap
{
namespace
{
struct SwathKey
{
    std::uint16_t psid = 0;
    std::uint8_t channel = 0;
    std::uint16_t file = 0;

    bool operator==(const SwathKey& other) const noexcept
    {
        return psid == other.psid && channel == other.channel && file == other.file;
    }
};

struct SwathKeyHash
{
    std::size_t operator()(const SwathKey& key) const noexcept
    {
        std::size_t seed = static_cast<std::size_t>(key.psid);
        seed ^= static_cast<std::size_t>(key.channel) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        seed ^= static_cast<std::size_t>(key.file) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        return seed;
    }
};

struct SwathCount
{
    SwathKey key;
    std::size_t count = 0;
};

inline SwathKey makeSwathKey(SwathKeyMode mode,
                             std::uint16_t psid,
                             std::uint16_t fileId,
                             std::uint8_t channel)
{
    SwathKey key;
    key.psid = psid;
    if (mode == SwathKeyMode::PointSourceIdChannel)
        key.channel = channel;
    else if (mode == SwathKeyMode::PointSourceIdFile)
        key.file = fileId;
    return key;
}

inline std::size_t flattenIndex(std::size_t ix, std::size_t iy, std::size_t nx)
{
    return iy * nx + ix;
}

inline std::pair<std::size_t, std::size_t> cellCoords(std::size_t idx, std::size_t nx)
{
    return {idx % nx, idx / nx};
}

} // namespace

OverlapMask buildOverlapMask(const pdal::PointView& view,
                             const OverlapParams& params,
                             OverlapStats* stats)
{
    if (params.cellSize <= 0.0)
        throw std::invalid_argument("Overlap cell size must be positive");

    OverlapMask mask;
    mask.cellSize = params.cellSize;
    mask.swathKeyUsed = params.swathKey;

    if (stats)
    {
        stats->overlapCells = 0;
        stats->overlapPointsFlagged = 0;
        stats->overlapPointCandidates = 0;
        stats->usedOverlapDimension = false;
        stats->downgradedToClassification = false;
        stats->swathKeyUsed = params.swathKey;
    }

    if (view.size() == 0)
        return mask;

    auto layout = view.layout();
    const auto dimX = pdal::Dimension::Id::X;
    const auto dimY = pdal::Dimension::Id::Y;
    const auto dimPsid = pdal::Dimension::Id::PointSourceId;
    const auto dimScanner = pdal::Dimension::Id::ScanChannel;

    const bool hasScanner = layout && layout->hasDim(dimScanner);
    SwathKeyMode mode = params.swathKey;
    if (mode == SwathKeyMode::PointSourceIdChannel && !hasScanner)
    {
        std::cerr << "Warning: ScanChannel dimension not present; using PointSourceId for swath key.\n";
        mode = SwathKeyMode::PointSourceId;
    }
    if (mode == SwathKeyMode::PointSourceIdFile)
    {
        std::cerr << "Warning: psid+file swath key not supported for this view; falling back to PointSourceId.\n";
        mode = SwathKeyMode::PointSourceId;
    }
    mask.swathKeyUsed = mode;
    if (stats)
        stats->swathKeyUsed = mode;

    const auto srs = view.spatialReference();
    const bool geographic = !srs.empty() && srs.isGeographic();
    mask.geographic = geographic;

    constexpr double earthRadius = 6378137.0;
    constexpr double pi = 3.14159265358979323846;
    constexpr double degToRad = pi / 180.0;

    double originLon = 0.0;
    double originLat = 0.0;
    double cosLat0 = 1.0;

    if (geographic)
    {
        originLon = view.getFieldAs<double>(dimX, 0);
        originLat = view.getFieldAs<double>(dimY, 0);
        cosLat0 = std::cos(originLat * degToRad);
        if (std::abs(cosLat0) < 1e-8)
            cosLat0 = (cosLat0 < 0 ? -1.0 : 1.0) * 1e-8;
    }

    mask.originLon = originLon;
    mask.originLat = originLat;
    mask.cosLat0 = cosLat0;

    std::vector<double> localX(view.size());
    std::vector<double> localY(view.size());

    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();

    for (pdal::PointId pid = 0; pid < view.size(); ++pid)
    {
        double x = view.getFieldAs<double>(dimX, pid);
        double y = view.getFieldAs<double>(dimY, pid);

        double lx = x;
        double ly = y;
        if (geographic)
        {
            const double dLon = (x - originLon) * degToRad;
            const double dLat = (y - originLat) * degToRad;
            lx = earthRadius * dLon * cosLat0;
            ly = earthRadius * dLat;
        }

        localX[pid] = lx;
        localY[pid] = ly;

        minX = std::min(minX, lx);
        minY = std::min(minY, ly);
        maxX = std::max(maxX, lx);
        maxY = std::max(maxY, ly);
    }

    const double cellSize = params.cellSize;
    const double spanX = std::max(cellSize, maxX - minX + cellSize * 0.5);
    const double spanY = std::max(cellSize, maxY - minY + cellSize * 0.5);

    const std::size_t nx = static_cast<std::size_t>(std::ceil(spanX / cellSize));
    const std::size_t ny = static_cast<std::size_t>(static_cast<std::size_t>(std::ceil(spanY / cellSize)));

    mask.originX = minX;
    mask.originY = minY;
    mask.nx = nx;
    mask.ny = ny;

    std::vector<std::vector<SwathCount>> swathCounts(nx * ny);
    std::vector<std::vector<pdal::PointId>> cellPoints(nx * ny);
    mask.cellIndexPerPoint.assign(view.size(), -1);

    auto pushSwath = [](std::vector<SwathCount>& counts, const SwathKey& key) {
        for (auto& entry : counts)
        {
            if (entry.key == key)
            {
                ++entry.count;
                return;
            }
        }
        counts.push_back(SwathCount{key, 1});
    };

    for (pdal::PointId pid = 0; pid < view.size(); ++pid)
    {
        const double lx = localX[pid] - minX;
        const double ly = localY[pid] - minY;
        std::size_t ix = static_cast<std::size_t>(std::floor(lx / cellSize));
        std::size_t iy = static_cast<std::size_t>(std::floor(ly / cellSize));
        if (ix >= nx)
            ix = nx - 1;
        if (iy >= ny)
            iy = ny - 1;
        const std::size_t cellId = flattenIndex(ix, iy, nx);

        mask.cellIndexPerPoint[pid] = static_cast<std::int32_t>(cellId);
        cellPoints[cellId].push_back(pid);

        const auto psid = view.getFieldAs<std::uint16_t>(dimPsid, pid);
        const auto channel = hasScanner ? view.getFieldAs<std::uint8_t>(dimScanner, pid)
                                        : static_cast<std::uint8_t>(0);
        const SwathKey key = makeSwathKey(mode, psid, 0, channel);
        pushSwath(swathCounts[cellId], key);
    }

    std::vector<std::uint8_t> cellMask(nx * ny, 0);
    for (std::size_t idx = 0; idx < swathCounts.size(); ++idx)
    {
        std::size_t qualifying = 0;
        for (const auto& entry : swathCounts[idx])
        {
            if (entry.count >= params.minPointsPerSwath)
                ++qualifying;
            if (qualifying >= 2)
            {
                cellMask[idx] = 1;
                break;
            }
        }
    }

    if (params.dilate > 0)
    {
        std::vector<std::uint8_t> dilated = cellMask;
        for (std::size_t idx = 0; idx < cellMask.size(); ++idx)
        {
            if (!cellMask[idx])
                continue;
            const auto [ix, iy] = cellCoords(idx, nx);
            for (int dx = -params.dilate; dx <= params.dilate; ++dx)
            {
                for (int dy = -params.dilate; dy <= params.dilate; ++dy)
                {
                    const std::int64_t nxIdx = static_cast<std::int64_t>(ix) + dx;
                    const std::int64_t nyIdx = static_cast<std::int64_t>(iy) + dy;
                    if (nxIdx < 0 || nyIdx < 0 || nxIdx >= static_cast<std::int64_t>(nx) ||
                        nyIdx >= static_cast<std::int64_t>(ny))
                        continue;
                    dilated[flattenIndex(static_cast<std::size_t>(nxIdx),
                                          static_cast<std::size_t>(nyIdx), nx)] = 1;
                }
            }
        }
        cellMask.swap(dilated);
    }

    std::vector<std::uint8_t> pointMask(view.size(), 0);
    std::size_t candidatePoints = 0;
    std::size_t overlapCellCount = 0;

    for (std::size_t idx = 0; idx < cellMask.size(); ++idx)
    {
        if (!cellMask[idx])
            continue;
        ++overlapCellCount;
        for (pdal::PointId pid : cellPoints[idx])
        {
            if (!pointMask[pid])
            {
                pointMask[pid] = 1;
                ++candidatePoints;
            }
        }
    }

    mask.cellMask = std::move(cellMask);
    mask.pointMask = std::move(pointMask);

    if (stats)
    {
        stats->overlapCells = overlapCellCount;
        stats->overlapPointCandidates = candidatePoints;
    }

    return mask;
}

void applyOverlapMask(pdal::PointView& view,
                      const OverlapMask& mask,
                      bool overwrite,
                      OverlapStats* stats)
{
    if (view.size() != mask.pointMask.size())
        throw std::runtime_error("Overlap mask does not match point view size");

    auto layout = view.layout();
    const auto dimOverlap = pdal::Dimension::Id::Overlap;
    const auto dimClassFlags = pdal::Dimension::Id::ClassFlags;
    const auto dimClassification = pdal::Dimension::Id::Classification;

    const bool hasOverlapDim = layout && layout->hasDim(dimOverlap);
    const bool hasClassFlags = layout && layout->hasDim(dimClassFlags);
    const bool hasClassification = layout && layout->hasDim(dimClassification);

    std::size_t flaggedCount = 0;
    bool usedOverlapDim = false;
    bool downgradedToClass = false;

    if (hasOverlapDim)
    {
        for (pdal::PointId pid = 0; pid < view.size(); ++pid)
        {
            const bool mark = mask.pointMask[pid] != 0;
            auto current = view.getFieldAs<std::uint8_t>(dimOverlap, pid);

            if (!overwrite && !mark)
            {
                flaggedCount += current ? 1 : 0;
                continue;
            }

            std::uint8_t finalValue = overwrite ? static_cast<std::uint8_t>(mark)
                                                : static_cast<std::uint8_t>(current | (mark ? 1 : 0));
            view.setField(dimOverlap, pid, finalValue);
            flaggedCount += finalValue ? 1 : 0;
        }
        usedOverlapDim = true;
    }
    else if (hasClassFlags)
    {
        const std::uint8_t overlapBit = static_cast<std::uint8_t>(1u << 3);
        for (pdal::PointId pid = 0; pid < view.size(); ++pid)
        {
            const bool mark = mask.pointMask[pid] != 0;
            auto current = view.getFieldAs<std::uint8_t>(dimClassFlags, pid);

            if (!overwrite && !mark)
            {
                flaggedCount += (current & overlapBit) ? 1 : 0;
                continue;
            }

            if (overwrite)
                current = mark ? static_cast<std::uint8_t>(current | overlapBit)
                                : static_cast<std::uint8_t>(current & ~overlapBit);
            else if (mark)
                current |= overlapBit;

            view.setField(dimClassFlags, pid, current);
            flaggedCount += (current & overlapBit) ? 1 : 0;
        }
    }
    else if (hasClassification)
    {
        for (pdal::PointId pid = 0; pid < view.size(); ++pid)
        {
            const bool mark = mask.pointMask[pid] != 0;
            if (!mark)
                continue;
            view.setField(dimClassification, pid, static_cast<std::uint8_t>(12));
            ++flaggedCount;
        }
        downgradedToClass = true;
    }
    else
    {
        if (stats)
            stats->downgradedToClassification = true;
        std::cerr << "Warning: Unable to flag overlap (no Overlap, ClassFlags, or Classification dimensions).\n";
    }

    if (stats)
    {
        stats->overlapCells = 0;
        for (auto value : mask.cellMask)
            stats->overlapCells += value ? 1 : 0;
        stats->overlapPointCandidates = 0;
        for (auto value : mask.pointMask)
            stats->overlapPointCandidates += value ? 1 : 0;
        stats->overlapPointsFlagged = flaggedCount;
        stats->usedOverlapDimension = usedOverlapDim;
        stats->downgradedToClassification = downgradedToClass;
        stats->swathKeyUsed = mask.swathKeyUsed;
    }
}

} // namespace overlap
