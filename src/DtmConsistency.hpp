#pragma once

#include <pdal/PointView.hpp>

#include <cstddef>
#include <cstdint>
#include <string>

struct DtmConsistencyOptions
{
    std::string dtmPath;
    double threshold = 1.0;            // meters
    std::uint8_t reclassValue = 8;     // ASPRS Model Key default
    enum class Action
    {
        Reclassify,
        Delete
    } action = Action::Reclassify;
};

struct DtmConsistencyStats
{
    std::size_t groundTested = 0;
    std::size_t reclassified = 0;
    std::size_t outsideRaster = 0;
    std::size_t nodataSkipped = 0;
    std::size_t removed = 0;
};

pdal::PointViewPtr applyDtmConsistency(pdal::PointViewPtr view,
                                       const DtmConsistencyOptions& options,
                                       DtmConsistencyStats* stats = nullptr);
