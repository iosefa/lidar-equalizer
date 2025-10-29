#include "Equalizer.hpp"

#include <pdal/Options.hpp>
#include <pdal/Metadata.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/StageFactory.hpp>
#include <pdal/io/BufferReader.hpp>

#include "OverlapFlagger.hpp"
#include "OverlapOps.hpp"
#include "OverlapSelfTest.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cstring>

namespace
{
using pdal::MetadataNode;

void printUsage(const char* argv0)
{
    std::cerr << "Usage: " << argv0
              << " input.laz output_equalized.laz [cell_size] [target_density] [seed]"
              << " [--class-scope=all|nonground|ground]"
              << " [--flag-overlap-only|--equalize-overlap-only|--join-overlap-by-scan-angle]"
              << "\n       or: " << argv0 << " --self-test-overlap\n";
}

MetadataNode findChild(const MetadataNode& root, const std::string& key)
{
    return root.findChild(key);
}

template <typename T>
void addOptionIfAvailable(pdal::Options& opts, const MetadataNode& root,
                          const std::string& key)
{
    auto node = findChild(root, key);
    if (!node.empty())
        opts.add(key, node.value<T>());
}

inline MetadataNode readerMetadata(pdal::Stage* reader)
{
    if (!reader)
        return MetadataNode{};
    MetadataNode root = reader->getMetadata();
    MetadataNode meta = root.findChild("metadata");
    return meta.empty() ? root : meta;
}
} // namespace

std::string toLowerCopy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

ClassScope parseScope(const std::string& value)
{
    const std::string lower = toLowerCopy(value);
    if (lower == "all")
        return ClassScope::All;
    if (lower == "nonground")
        return ClassScope::NonGround;
    if (lower == "ground")
        return ClassScope::Ground;
    throw std::invalid_argument("Invalid class scope: " + value +
                                " (expected all|nonground|ground)");
}

std::string scopeToString(ClassScope scope)
{
    switch (scope)
    {
    case ClassScope::All:
        return "all";
    case ClassScope::NonGround:
        return "nonground";
    case ClassScope::Ground:
        return "ground";
    }
    return "all";
}

bool parseBoolOption(const std::string& value)
{
    const std::string lower = toLowerCopy(value);
    if (lower == "1" || lower == "true" || lower == "yes" || lower == "on")
        return true;
    if (lower == "0" || lower == "false" || lower == "no" || lower == "off")
        return false;
    throw std::invalid_argument("Invalid boolean value: " + value);
}

overlap::SwathKeyMode parseSwathKey(const std::string& value)
{
    const std::string lower = toLowerCopy(value);
    if (lower == "psid" || lower == "pointsourceid")
        return overlap::SwathKeyMode::PointSourceId;
    if (lower == "psid+channel" || lower == "psid-channel" || lower == "psid_channel")
        return overlap::SwathKeyMode::PointSourceIdChannel;
    if (lower == "psid+file" || lower == "psid-file" || lower == "psid_file")
        return overlap::SwathKeyMode::PointSourceIdFile;
    if (lower == "channel" || lower == "scanner-channel" || lower == "scanner_channel")
        return overlap::SwathKeyMode::PointSourceIdChannel; // will be handled later if unavailable
    throw std::invalid_argument("Invalid swath key: " + value + " (expected psid|psid-channel|psid-file|channel)");
}

std::string swathKeyToString(overlap::SwathKeyMode mode)
{
    switch (mode)
    {
    case overlap::SwathKeyMode::PointSourceId:
        return "psid";
    case overlap::SwathKeyMode::PointSourceIdChannel:
        return "psid-channel";
    case overlap::SwathKeyMode::PointSourceIdFile:
        return "psid-file";
    }
    return "psid";
}

enum class Mode
{
    EqualizeAll,
    OverlapFlagOnly,
    EqualizeOverlapOnly,
    JoinOverlapByScanAngle
};

int main(int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "--self-test-overlap")
        return overlap::runOverlapSelfTest() ? 0 : 1;

    if (argc < 3)
    {
        printUsage(argv[0]);
        return 1;
    }

    const std::string input = argv[1];
    const std::string output = argv[2];

    double cell = 1.0;
    double target = 6.0;
    std::uint64_t seed = 42ULL;
    ClassScope scope = ClassScope::All;
    overlap::OverlapParams overlapParams;
    bool overwriteOverlap = false;
    Mode mode = Mode::EqualizeAll;

    auto setMode = [&](Mode desired, bool enable) {
        if (!enable)
            return;
        if (mode != Mode::EqualizeAll && mode != desired)
            throw std::invalid_argument("Processing modes are mutually exclusive");
        mode = desired;
    };

    try
    {
        int argIndex = 3;

        auto isFlag = [](const char* s) { return std::strncmp(s, "--", 2) == 0; };

        if (argIndex < argc && !isFlag(argv[argIndex]))
        {
            cell = std::stod(argv[argIndex]);
            ++argIndex;
        }
        if (argIndex < argc && !isFlag(argv[argIndex]))
        {
            target = std::stod(argv[argIndex]);
            ++argIndex;
        }
        if (argIndex < argc && !isFlag(argv[argIndex]))
        {
            seed = std::stoull(argv[argIndex]);
            ++argIndex;
        }

        overlapParams.cellSize = cell;

        while (argIndex < argc)
        {
            std::string arg(argv[argIndex]);
            if (arg.rfind("--class-scope", 0) == 0)
            {
                std::string value;
                auto eq = arg.find('=');
                if (eq != std::string::npos)
                    value = arg.substr(eq + 1);
                else
                {
                    if (argIndex + 1 >= argc)
                        throw std::invalid_argument("--class-scope requires a value");
                    value = argv[++argIndex];
                }
                scope = parseScope(value);
            }
            else if (arg.rfind("--cell-size", 0) == 0)
            {
                std::string value;
                auto eq = arg.find('=');
                if (eq != std::string::npos)
                    value = arg.substr(eq + 1);
                else
                {
                    if (argIndex + 1 >= argc)
                        throw std::invalid_argument("--cell-size requires a value");
                    value = argv[++argIndex];
                }
                cell = std::stod(value);
                overlapParams.cellSize = cell;
            }
            else if (arg.rfind("--flag-overlap-only", 0) == 0)
            {
                bool val = true;
                auto eq = arg.find('=');
                if (eq != std::string::npos)
                    val = parseBoolOption(arg.substr(eq + 1));
                else if (argIndex + 1 < argc && !isFlag(argv[argIndex + 1]))
                    val = parseBoolOption(argv[++argIndex]);
                setMode(Mode::OverlapFlagOnly, val);
            }
            else if (arg.rfind("--equalize-overlap-only", 0) == 0)
            {
                bool val = true;
                auto eq = arg.find('=');
                if (eq != std::string::npos)
                    val = parseBoolOption(arg.substr(eq + 1));
                else if (argIndex + 1 < argc && !isFlag(argv[argIndex + 1]))
                    val = parseBoolOption(argv[++argIndex]);
                setMode(Mode::EqualizeOverlapOnly, val);
            }
            else if (arg.rfind("--join-overlap-by-scan-angle", 0) == 0)
            {
                bool val = true;
                auto eq = arg.find('=');
                if (eq != std::string::npos)
                    val = parseBoolOption(arg.substr(eq + 1));
                else if (argIndex + 1 < argc && !isFlag(argv[argIndex + 1]))
                    val = parseBoolOption(argv[++argIndex]);
                setMode(Mode::JoinOverlapByScanAngle, val);
            }
            else if (arg.rfind("--overwrite-overlap", 0) == 0)
            {
                bool val = true;
                auto eq = arg.find('=');
                if (eq != std::string::npos)
                    val = parseBoolOption(arg.substr(eq + 1));
                else if (argIndex + 1 < argc && !isFlag(argv[argIndex + 1]))
                    val = parseBoolOption(argv[++argIndex]);
                overwriteOverlap = val;
            }
            else if (arg.rfind("--min-points-per-psid", 0) == 0)
            {
                std::string value;
                auto eq = arg.find('=');
                if (eq != std::string::npos)
                    value = arg.substr(eq + 1);
                else
                {
                    if (argIndex + 1 >= argc)
                        throw std::invalid_argument("--min-points-per-psid requires a value");
                    value = argv[++argIndex];
                }
                overlapParams.minPointsPerSwath = static_cast<std::size_t>(std::stoul(value));
            }
            else if (arg.rfind("--overlap-dilate", 0) == 0)
            {
                std::string value;
                auto eq = arg.find('=');
                if (eq != std::string::npos)
                    value = arg.substr(eq + 1);
                else
                {
                    if (argIndex + 1 >= argc)
                        throw std::invalid_argument("--overlap-dilate requires a value");
                    value = argv[++argIndex];
                }
                overlapParams.dilate = std::stoi(value);
                if (overlapParams.dilate < 0)
                    overlapParams.dilate = 0;
            }
            else if (arg.rfind("--swath-key", 0) == 0)
            {
                std::string value;
                auto eq = arg.find('=');
                if (eq != std::string::npos)
                    value = arg.substr(eq + 1);
                else
                {
                    if (argIndex + 1 >= argc)
                        throw std::invalid_argument("--swath-key requires a value");
                    value = argv[++argIndex];
                }
                overlapParams.swathKey = parseSwathKey(value);
            }
            else
            {
                throw std::invalid_argument("Unknown argument: " + arg);
            }
            ++argIndex;
        }

        auto tStart = std::chrono::steady_clock::now();

        pdal::StageFactory factory;
        std::string readerDriver = pdal::StageFactory::inferReaderDriver(input);
        if (readerDriver == "readers.copc" || readerDriver.empty())
            readerDriver = "readers.las";

        auto reader = factory.createStage(readerDriver);
        if (!reader)
            throw std::runtime_error("Failed to create PDAL reader stage");

        pdal::Options readerOpts;
        readerOpts.add("filename", input);
        reader->setOptions(readerOpts);

        pdal::PointTable table;
        reader->prepare(table);
        pdal::PointViewSet viewSet = reader->execute(table);
        if (viewSet.empty())
            throw std::runtime_error("Input dataset produced no point views");

        pdal::PointViewPtr view;
        for (const auto& candidate : viewSet)
        {
            if (!candidate)
                continue;
            if (!view)
                view = candidate;
            else
                view->append(*candidate);
        }
        if (!view)
            throw std::runtime_error("Failed to obtain point view from reader");

        const MetadataNode meta = readerMetadata(reader);

        auto inferWriter = [](const std::string& path) -> std::string {
            auto endsWith = [](const std::string& s, const std::string& suffix) {
                return s.size() >= suffix.size() &&
                       std::equal(suffix.rbegin(), suffix.rend(), s.rbegin(),
                                  [](char a, char b) {
                                      return std::tolower(static_cast<unsigned char>(a)) ==
                                             std::tolower(static_cast<unsigned char>(b));
                                  });
            };
            if (endsWith(path, ".copc.laz") || endsWith(path, ".copc.las"))
                return "writers.copc";
            if (endsWith(path, ".laz") || endsWith(path, ".las"))
                return "writers.las";
            return pdal::StageFactory::inferWriterDriver(path);
        };

        const std::string writerDriver = inferWriter(output);

        auto prepareWriterOptions = [&](pdal::Options& writerOpts) {
            writerOpts.add("filename", output);
            if (writerDriver != "writers.copc")
                writerOpts.add("forward", "all");

            if (writerDriver == "writers.las")
            {
                addOptionIfAvailable<int>(writerOpts, meta, "minor_version");
                addOptionIfAvailable<int>(writerOpts, meta, "major_version");
                addOptionIfAvailable<int>(writerOpts, meta, "dataformat_id");
                addOptionIfAvailable<double>(writerOpts, meta, "scale_x");
                addOptionIfAvailable<double>(writerOpts, meta, "scale_y");
                addOptionIfAvailable<double>(writerOpts, meta, "scale_z");
                addOptionIfAvailable<double>(writerOpts, meta, "offset_x");
                addOptionIfAvailable<double>(writerOpts, meta, "offset_y");
                addOptionIfAvailable<double>(writerOpts, meta, "offset_z");
                addOptionIfAvailable<int>(writerOpts, meta, "global_encoding");
                addOptionIfAvailable<int>(writerOpts, meta, "filesource_id");
                addOptionIfAvailable<int>(writerOpts, meta, "creation_year");
                addOptionIfAvailable<int>(writerOpts, meta, "creation_doy");
                addOptionIfAvailable<std::string>(writerOpts, meta, "project_id");
                addOptionIfAvailable<std::string>(writerOpts, meta, "software_id");
                writerOpts.add("extra_dims", "all");
            }
            else if (writerDriver == "writers.copc")
            {
                addOptionIfAvailable<double>(writerOpts, meta, "scale_x");
                addOptionIfAvailable<double>(writerOpts, meta, "scale_y");
                addOptionIfAvailable<double>(writerOpts, meta, "scale_z");
                addOptionIfAvailable<double>(writerOpts, meta, "offset_x");
                addOptionIfAvailable<double>(writerOpts, meta, "offset_y");
                addOptionIfAvailable<double>(writerOpts, meta, "offset_z");
                addOptionIfAvailable<int>(writerOpts, meta, "global_encoding");
                addOptionIfAvailable<int>(writerOpts, meta, "filesource_id");
                addOptionIfAvailable<int>(writerOpts, meta, "creation_year");
                addOptionIfAvailable<int>(writerOpts, meta, "creation_doy");
                addOptionIfAvailable<std::string>(writerOpts, meta, "project_id");
                addOptionIfAvailable<std::string>(writerOpts, meta, "software_id");
            }
        };

        auto writeViewToFile = [&](pdal::PointViewPtr pv) {
            pdal::BufferReader bufferReader;
            bufferReader.addView(pv);
            if (!pv->spatialReference().empty())
                bufferReader.setSpatialReference(pv->spatialReference());

            pdal::StageFactory writerFactory;
            auto writerStage = writerFactory.createStage(writerDriver);
            if (!writerStage)
                throw std::runtime_error("Failed to create PDAL writer stage");

            writerStage->setInput(bufferReader);

            pdal::Options writerOpts;
            prepareWriterOptions(writerOpts);
            writerStage->setOptions(writerOpts);
            if (!pv->spatialReference().empty())
                writerStage->setSpatialReference(pv->spatialReference());

            pdal::PointTable outTable;
            writerStage->prepare(outTable);
            writerStage->execute(outTable);
        };

        auto maybeBuildMask = [&]() {
            overlap::OverlapStats maskStats;
            auto mask = overlap::buildOverlapMask(*view, overlapParams, &maskStats);
            return std::make_pair(mask, maskStats);
        };

        const auto modeStart = std::chrono::steady_clock::now();

        if (mode == Mode::OverlapFlagOnly)
        {
            auto [mask, maskStats] = maybeBuildMask();
            overlap::OverlapStats applyStats;
            overlap::applyOverlapMask(*view, mask, overwriteOverlap, &applyStats);
            writeViewToFile(view);

            const auto elapsed =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - modeStart).count();

            std::cout << "Mode         : overlap-flag" << '\n';
            std::cout << "Input points : " << view->size() << '\n';
            std::cout << "Overlap cells: " << maskStats.overlapCells << '\n';
            std::cout << "Points flagged: " << applyStats.overlapPointsFlagged
                      << " (candidates " << maskStats.overlapPointCandidates << ")" << '\n';
            std::cout << "Swath key    : " << swathKeyToString(maskStats.swathKeyUsed) << '\n';
            if (applyStats.downgradedToClassification)
                std::cout << "Warning: tagged overlap via Classification=12 (legacy format)." << '\n';
            std::cout << "Done in " << std::fixed << std::setprecision(1) << elapsed << " s." << '\n';
            return 0;
        }

        if (mode == Mode::EqualizeOverlapOnly)
        {
            auto [mask, maskStats] = maybeBuildMask();
            overlap::EqualizeOverlapStats eqStats;
            auto equalized = overlap::equalizeOverlapOnly(view, mask, target, seed, &eqStats);
            writeViewToFile(equalized);

            const auto elapsed =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - modeStart).count();

            std::cout << "Mode         : equalize-overlap-only" << '\n';
            std::cout << "Input points : " << view->size() << '\n';
            std::cout << "Overlap cells: " << maskStats.overlapCells << '\n';
            std::cout << "Overlap kept : " << eqStats.keptOverlapPoints
                      << " / " << eqStats.overlapPoints
                      << " (dropped " << eqStats.droppedOverlapPoints << ")" << '\n';
            std::cout << "Output points: " << equalized->size() << '\n';
            std::cout << "Done in " << std::fixed << std::setprecision(1) << elapsed << " s." << '\n';
            return 0;
        }

        if (mode == Mode::JoinOverlapByScanAngle)
        {
            auto [mask, maskStats] = maybeBuildMask();
            overlap::JoinOverlapStats joinStats;
            auto joined = overlap::joinOverlapByScanAngle(view, mask, &joinStats);
            writeViewToFile(joined);

            const auto elapsed =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - modeStart).count();

            std::cout << "Mode         : join-overlap-by-scan-angle" << '\n';
            std::cout << "Input points : " << view->size() << '\n';
            std::cout << "Overlap cells: " << maskStats.overlapCells << '\n';
            std::cout << "Overlap kept : " << joinStats.keptOverlapPoints
                      << " / " << joinStats.overlapPoints
                      << " (dropped " << joinStats.droppedOverlapPoints << ")" << '\n';
            std::cout << "Output points: " << joined->size() << '\n';
            std::cout << "Done in " << std::fixed << std::setprecision(1) << elapsed << " s." << '\n';
            return 0;
        }

        EqualizeOptions options;
        options.cell = cell;
        options.target = target;
        options.seed = seed;
        options.scope = scope;

        auto equalized = proportionalEqualize(view, options);
        writeViewToFile(equalized);

        const auto inputPoints = static_cast<std::uint64_t>(view->size());
        const auto keptPoints = static_cast<std::uint64_t>(equalized->size());
        const double pct = inputPoints == 0 ? 0.0
                                            : (100.0 * static_cast<double>(keptPoints) /
                                               static_cast<double>(inputPoints));
        const auto elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count();

        std::cout << "Scope        : " << scopeToString(options.scope) << '\n';
        std::cout << "Input points : " << inputPoints << '\n';
        std::cout << "Kept points  : " << keptPoints << " ("
                  << std::fixed << std::setprecision(0) << pct << "%)\n";
        std::cout << "Target density: " << std::fixed << std::setprecision(1) << options.target
                  << " pts/m^2\n";
        std::cout << "Cell size    : " << std::fixed << std::setprecision(2) << options.cell
                  << " m\n";
        std::cout << "Done in " << std::fixed << std::setprecision(1) << elapsed << " s.\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
