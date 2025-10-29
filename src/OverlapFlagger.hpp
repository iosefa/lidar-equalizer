#pragma once

#include <pdal/PointView.hpp>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace overlap
{

enum class SwathKeyMode
{
    PointSourceId,
    PointSourceIdChannel,
    PointSourceIdFile
};

struct OverlapParams
{
    double cellSize = 1.0;
    std::size_t minPointsPerSwath = 1;
    int dilate = 0;
    SwathKeyMode swathKey = SwathKeyMode::PointSourceId;
};

struct OverlapMask
{
    double cellSize = 1.0;
    bool geographic = false;
    double originX = 0.0;
    double originY = 0.0;
    double originLon = 0.0;
    double originLat = 0.0;
    double cosLat0 = 1.0;
    std::size_t nx = 0;
    std::size_t ny = 0;
    SwathKeyMode swathKeyUsed = SwathKeyMode::PointSourceId;
    std::vector<std::int32_t> cellIndexPerPoint;
    std::vector<std::uint8_t> pointMask;
    std::vector<std::uint8_t> cellMask;
};

struct OverlapStats
{
    std::size_t overlapCells = 0;
    std::size_t overlapPointsFlagged = 0;
    std::size_t overlapPointCandidates = 0;
    bool usedOverlapDimension = false;
    bool downgradedToClassification = false;
    SwathKeyMode swathKeyUsed = SwathKeyMode::PointSourceId;
};

OverlapMask buildOverlapMask(const pdal::PointView& view,
                             const OverlapParams& params,
                             OverlapStats* stats = nullptr);

void applyOverlapMask(pdal::PointView& view,
                      const OverlapMask& mask,
                      bool overwrite,
                      OverlapStats* stats = nullptr);

} // namespace overlap

