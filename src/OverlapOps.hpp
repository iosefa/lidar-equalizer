#pragma once

#include "OverlapFlagger.hpp"

#include <pdal/PointView.hpp>

#include <cstdint>

namespace overlap
{

struct EqualizeOverlapStats
{
    std::size_t overlapPoints = 0;
    std::size_t keptOverlapPoints = 0;
    std::size_t droppedOverlapPoints = 0;
};

struct JoinOverlapStats
{
    std::size_t overlapPoints = 0;
    std::size_t keptOverlapPoints = 0;
    std::size_t droppedOverlapPoints = 0;
    std::size_t groundPreserved = 0;
};

pdal::PointViewPtr equalizeOverlapOnly(pdal::PointViewPtr view,
                                       const OverlapMask& mask,
                                       double targetDensity,
                                       std::uint64_t seed,
                                       EqualizeOverlapStats* stats = nullptr);

pdal::PointViewPtr joinOverlapByScanAngle(pdal::PointViewPtr view,
                                          const OverlapMask& mask,
                                          bool keepGround,
                                          JoinOverlapStats* stats = nullptr);

} // namespace overlap
