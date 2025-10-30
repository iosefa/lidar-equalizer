#include "DtmConsistency.hpp"

#include <gdal_priv.h>
#include <cpl_conv.h>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{
struct GdalDatasetHandle
{
    GDALDataset* handle = nullptr;
    explicit GdalDatasetHandle(GDALDataset* ds) : handle(ds) {}
    ~GdalDatasetHandle()
    {
        if (handle)
            GDALClose(handle);
    }
};

inline bool nearlyEqual(double a, double b)
{
    return std::abs(a - b) <= 1e-9 * std::max(std::abs(a), std::abs(b));
}
} // namespace

pdal::PointViewPtr applyDtmConsistency(pdal::PointViewPtr view,
                                       const DtmConsistencyOptions& options,
                                       DtmConsistencyStats* stats)
{
    if (!view)
        throw std::invalid_argument("DTM consistency filter requires valid PointView");
    if (options.dtmPath.empty())
        throw std::invalid_argument("DTM consistency filter requires --dtm path");
    if (options.threshold < 0.0)
        throw std::invalid_argument("DTM threshold must be non-negative");

    auto layout = view->layout();
    if (!layout || !layout->hasDim(pdal::Dimension::Id::Classification))
        throw std::runtime_error("PointView missing Classification dimension for DTM consistency filtering");

    GDALAllRegister();
    GDALDataset* dsRaw = static_cast<GDALDataset*>(
        GDALOpen(options.dtmPath.c_str(), GA_ReadOnly));
    if (!dsRaw)
        throw std::runtime_error("Failed to open DTM raster: " + options.dtmPath);

    GdalDatasetHandle dataset(dsRaw);

    double geoTransform[6];
    if (dataset.handle->GetGeoTransform(geoTransform) != CE_None)
        throw std::runtime_error("DTM raster lacks affine geotransform");

    double invTransform[6];
    if (!GDALInvGeoTransform(geoTransform, invTransform))
        throw std::runtime_error("Failed to invert DTM geotransform");

    const int width = dataset.handle->GetRasterXSize();
    const int height = dataset.handle->GetRasterYSize();
    if (width <= 0 || height <= 0)
        throw std::runtime_error("DTM raster has invalid size");

    GDALRasterBand* band = dataset.handle->GetRasterBand(1);
    if (!band)
        throw std::runtime_error("DTM raster missing band 1");

    int success = 0;
    const double nodataValue = band->GetNoDataValue(&success);
    const bool hasNoData = success != 0;

    std::vector<double> raster(width * height);
    const CPLErr err = band->RasterIO(GF_Read,
                                      0, 0,
                                      width, height,
                                      raster.data(),
                                      width, height,
                                      GDT_Float64,
                                      0, 0);
    if (err != CE_None)
        throw std::runtime_error("Failed to read DTM raster data");

    auto sampleDtm = [&](double x, double y, double& value) -> bool {
        const double pixel = invTransform[0] + invTransform[1] * x + invTransform[2] * y;
        const double line = invTransform[3] + invTransform[4] * x + invTransform[5] * y;
        const int col = static_cast<int>(std::floor(pixel));
        const int row = static_cast<int>(std::floor(line));

        if (col < 0 || col >= width || row < 0 || row >= height)
            return false;

        value = raster[static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
                       static_cast<std::size_t>(col)];
        if (hasNoData && (value == nodataValue || nearlyEqual(value, nodataValue)))
            return false;
        return true;
    };

    std::size_t groundTested = 0;
    std::size_t reclassified = 0;
    std::size_t outsideRaster = 0;
    std::size_t nodataSkipped = 0;
    std::size_t removed = 0;

    const auto dimX = pdal::Dimension::Id::X;
    const auto dimY = pdal::Dimension::Id::Y;
    const auto dimZ = pdal::Dimension::Id::Z;
    const auto dimClass = pdal::Dimension::Id::Classification;

    std::vector<char> keepMask;
    if (options.action == DtmConsistencyOptions::Action::Delete)
        keepMask.assign(view->size(), 1);

    for (pdal::PointId pid = 0; pid < view->size(); ++pid)
    {
        auto cls = view->getFieldAs<std::uint8_t>(dimClass, pid);
        if (cls != 2)
            continue;

        ++groundTested;

        const double x = view->getFieldAs<double>(dimX, pid);
        const double y = view->getFieldAs<double>(dimY, pid);
        const double z = view->getFieldAs<double>(dimZ, pid);

        double dtmZ = 0.0;
        if (!sampleDtm(x, y, dtmZ))
        {
            if (stats)
            {
                const double pixel = invTransform[0] + invTransform[1] * x + invTransform[2] * y;
                const double line = invTransform[3] + invTransform[4] * x + invTransform[5] * y;
                const int col = static_cast<int>(std::floor(pixel));
                const int row = static_cast<int>(std::floor(line));
                if (col < 0 || col >= width || row < 0 || row >= height)
                    ++outsideRaster;
                else
                    ++nodataSkipped;
            }
            continue;
        }

        if (std::abs(z - dtmZ) > options.threshold)
        {
            if (options.action == DtmConsistencyOptions::Action::Delete)
            {
                keepMask[pid] = 0;
                ++removed;
            }
            else
            {
                view->setField(dimClass, pid, options.reclassValue);
                ++reclassified;
            }
        }
    }

    pdal::PointViewPtr result = view;
    if (options.action == DtmConsistencyOptions::Action::Delete)
    {
        auto filtered = view->makeNew();
        for (pdal::PointId pid = 0; pid < view->size(); ++pid)
        {
            if (!keepMask[pid])
                continue;
            filtered->appendPoint(*view, pid);
        }
        result = filtered;
    }

    if (stats)
    {
        stats->groundTested = groundTested;
        stats->reclassified = reclassified;
        stats->outsideRaster = outsideRaster;
        stats->nodataSkipped = nodataSkipped;
        stats->removed = removed;
    }

    return result;
}
