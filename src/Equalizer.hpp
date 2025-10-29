#pragma once

#include <pdal/PointView.hpp>

#include <cstddef>
#include <cstdint>

enum class ClassScope {
    All,
    NonGround,
    Ground
};

struct EqualizeOptions {
    double cell = 1.0;
    double target = 6.0;
    std::uint64_t seed = 42;
    ClassScope scope = ClassScope::All;
};

pdal::PointViewPtr proportionalEqualize(pdal::PointViewPtr view, const EqualizeOptions& opt);
