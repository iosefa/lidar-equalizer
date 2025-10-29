#include "Equalizer.hpp"

#include <pdal/Dimension.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace
{
struct CellKey
{
    std::int64_t ix = 0;
    std::int64_t iy = 0;

    bool operator==(const CellKey& other) const noexcept
    {
        return ix == other.ix && iy == other.iy;
    }
};

struct CellKeyHash
{
    std::size_t operator()(const CellKey& key) const noexcept
    {
        std::size_t seed = static_cast<std::size_t>(key.ix);
        seed ^= static_cast<std::size_t>(key.iy) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        return seed;
    }
};

struct PointRef
{
    pdal::PointId id;
    double weight;
};

struct CellInfo
{
    std::size_t total = 0;
    std::unordered_map<std::uint16_t, std::vector<PointRef>> byPsid;
};

inline CellKey makeKey(double x, double y, double cellSize)
{
    const double inv = 1.0 / cellSize;
    const auto ix = static_cast<std::int64_t>(std::floor(x * inv));
    const auto iy = static_cast<std::int64_t>(std::floor(y * inv));
    return {ix, iy};
}
} // namespace

pdal::PointViewPtr proportionalEqualize(pdal::PointViewPtr view, const EqualizeOptions& opt)
{
    if (!view)
        throw std::invalid_argument("proportionalEqualize received a null PointViewPtr");
    if (opt.cell <= 0.0)
        throw std::invalid_argument("Cell size must be positive");
    if (opt.target < 0.0)
        throw std::invalid_argument("Target density must be non-negative");

    auto output = view->makeNew();
    if (!output)
        throw std::runtime_error("Failed to allocate output PointView");

    if (view->size() == 0)
        return output;

    std::unordered_map<CellKey, CellInfo, CellKeyHash> grid;
    grid.reserve(static_cast<std::size_t>(view->size() / 4 + 1));

    std::mt19937_64 rng(opt.seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    const auto dimX = pdal::Dimension::Id::X;
    const auto dimY = pdal::Dimension::Id::Y;
    const auto dimPsid = pdal::Dimension::Id::PointSourceId;
    const auto dimClass = pdal::Dimension::Id::Classification;

    std::vector<bool> keep(view->size(), false);
    const auto srs = view->spatialReference();
    const bool geographic = !srs.empty() && srs.isGeographic();

    constexpr double earthRadius = 6378137.0;
    constexpr double pi = 3.14159265358979323846;
    constexpr double degToRad = pi / 180.0;
    double originLon = 0.0;
    double originLat = 0.0;
    double cosLat0 = 1.0;

    if (geographic && view->size() > 0)
    {
        originLon = view->getFieldAs<double>(dimX, 0);
        originLat = view->getFieldAs<double>(dimY, 0);
        cosLat0 = std::cos(originLat * degToRad);
        if (std::abs(cosLat0) < 1e-8)
            cosLat0 = (cosLat0 < 0 ? -1.0 : 1.0) * 1e-8;
    }

    for (pdal::PointId pid = 0; pid < view->size(); ++pid)
    {
        const double x = view->getFieldAs<double>(dimX, pid);
        const double y = view->getFieldAs<double>(dimY, pid);
        const auto psid = view->getFieldAs<std::uint16_t>(dimPsid, pid);
        const auto classification = view->getFieldAs<std::uint8_t>(dimClass, pid);

        const bool isGround = classification == 2;
        bool candidate = true;
        switch (opt.scope)
        {
        case ClassScope::All:
            candidate = true;
            break;
        case ClassScope::Ground:
            candidate = isGround;
            break;
        case ClassScope::NonGround:
            candidate = !isGround;
            break;
        }

        if (!candidate)
        {
            keep[pid] = true;
            continue;
        }

        double keyX = x;
        double keyY = y;
        if (geographic)
        {
            const double dLon = (x - originLon) * degToRad;
            const double dLat = (y - originLat) * degToRad;
            keyX = earthRadius * dLon * cosLat0;
            keyY = earthRadius * dLat;
        }

        auto& cell = grid[makeKey(keyX, keyY, opt.cell)];
        cell.total += 1;
        auto& bucket = cell.byPsid[psid];
        bucket.push_back(PointRef{pid, dist(rng)});
    }

    const double cellArea = opt.cell * opt.cell;

    for (auto& [key, cell] : grid)
    {
        (void)key;
        if (cell.total == 0)
            continue;

        const double desiredTotal = opt.target * cellArea;

        for (auto& [psid, refs] : cell.byPsid)
        {
            (void)psid;
            const std::size_t count = refs.size();
            if (count == 0)
                continue;

            double expected = desiredTotal * (static_cast<double>(count) / static_cast<double>(cell.total));
            if (!std::isfinite(expected))
                expected = 0.0;

            std::size_t quota = static_cast<std::size_t>(std::llround(expected));
            quota = std::min<std::size_t>(quota, count);

            if (quota == count)
                continue;

            if (quota == 0)
            {
                refs.clear();
                continue;
            }

            auto pivot = refs.begin() + static_cast<std::ptrdiff_t>(quota);
            std::nth_element(refs.begin(), pivot, refs.end(),
                             [](const PointRef& a, const PointRef& b) { return a.weight < b.weight; });
            refs.resize(quota);
        }

        for (const auto& [psid, refs] : cell.byPsid)
        {
            (void)psid;
            for (const auto& ref : refs)
                keep[ref.id] = true;
        }
    }

    for (pdal::PointId pid = 0; pid < view->size(); ++pid)
    {
        if (keep[pid])
            output->appendPoint(*view, pid);
    }

    return output;
}
