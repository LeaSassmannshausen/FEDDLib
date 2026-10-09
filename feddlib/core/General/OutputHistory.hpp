#ifndef FEDD_OUTPUT_HISTORY_HPP
#define FEDD_OUTPUT_HISTORY_HPP

#include "feddlib/core/FEDDCore.hpp"
#include <Teuchos_CommHelpers.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace FEDD {
namespace output {

/// Execute file operations on one rank and propagate failures before collective I/O.
inline void onRank(const Teuchos::Comm<int>& comm, int rank, const std::function<void()>& operation)
{
    std::string error;
    if (comm.getRank() == rank) {
        try { operation(); }
        catch (const std::exception& e) { error = e.what(); }
    }
    int size = static_cast<int>(error.size());
    Teuchos::broadcast(comm, rank, 1, &size);
    error.resize(size);
    if (size) {
        Teuchos::broadcast(comm, rank, size, &error[0]);
        throw std::runtime_error("Output history: " + error);
    }
}

inline bool resume(const ParameterListPtr_Type& parameters)
{
    return !parameters.is_null() && parameters->sublist("Exporter").get("Resume output", false) &&
           parameters->sublist("Timestepping Parameter").get("Restart", false);
}

inline double restartTime(const ParameterListPtr_Type& parameters)
{
    const double time = parameters->sublist("Timestepping Parameter").get("Time step", 0.);
    if (!std::isfinite(time) || time < 0.) throw std::runtime_error("Invalid output restart time");
    return time;
}

inline std::string read(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot read " + filename);
    std::ostringstream contents;
    contents << input.rdbuf();
    if (input.bad()) throw std::runtime_error("Cannot read " + filename);
    return contents.str();
}

/// Replace text through a temporary file, leaving the old file intact on write failure.
inline void replace(const std::string& filename, const std::string& contents)
{
    const auto temporary = filename + ".output-tmp";
    std::ofstream stream(temporary, std::ios::trunc);
    stream << contents;
    stream.close();
    if (!stream) throw std::runtime_error("Cannot write " + temporary);
    std::filesystem::rename(temporary, filename);
}

/** @brief Archive only files owned by the exporters, in one folder per run.
 * Both options default to false. Resume output takes precedence over archiving.
 * The shared parameter list caches the numbered archive directory so text and
 * ParaView files from the same run stay together. Checkpoints are never moved.
 */
inline void archive(const ParameterListPtr_Type& parameters, const Teuchos::Comm<int>& comm,
                    const std::vector<std::string>& files)
{
    if (parameters.is_null() || parameters->sublist("Exporter").get("Resume output", false) ||
        !parameters->sublist("Exporter").get("Keep old output", false)) return;
    std::string directory = parameters->sublist("Exporter").get("_Old output archive", std::string(""));
    onRank(comm, 0, [&] {
        bool exists = false;
        for (const auto& file : files) exists = exists || std::filesystem::exists(file);
        if (!exists) return;
        if (directory.empty()) {
            const auto parent = std::filesystem::current_path() / "old_output";
            std::filesystem::create_directories(parent);
            for (unsigned i = 1; ; ++i) {
                const auto candidate = parent / ("run_" + std::to_string(i));
                if (std::filesystem::create_directory(candidate)) {
                    directory = candidate.string();
                    break;
                }
            }
        }
        // Preflight all destinations before moving any file.
        for (const auto& file : files) {
            if (!std::filesystem::exists(file)) continue;
            auto relative = std::filesystem::path(file).lexically_normal();
            if (relative.is_absolute()) relative = relative.lexically_relative(std::filesystem::current_path());
            if (relative.empty() || *relative.begin() == "..")
                throw std::runtime_error("Archive output paths must be inside the working directory: " + file);
            const auto destination = std::filesystem::path(directory) / relative;
            if (std::filesystem::exists(destination))
                throw std::runtime_error("Archive destination already exists: " + destination.string());
        }
        for (const auto& file : files) {
            if (!std::filesystem::exists(file)) continue;
            auto relative = std::filesystem::path(file).lexically_normal();
            if (relative.is_absolute()) relative = relative.lexically_relative(std::filesystem::current_path());
            const auto destination = std::filesystem::path(directory) / relative;
            std::filesystem::create_directories(destination.parent_path());
            std::filesystem::rename(file, destination);
        }
    });
    int size = static_cast<int>(directory.size());
    Teuchos::broadcast(comm, 0, 1, &size);
    directory.resize(size);
    if (size) Teuchos::broadcast(comm, 0, size, &directory[0]);
    parameters->sublist("Exporter").set("_Old output archive", directory);
}

struct Frame {
    double time;
    int index;
    std::string xml;
    std::set<std::string> groups;
};

/** @brief Parse the temporal collection generated by FEDDLib's XMF exporter.
 * Existing XMF files supply timestamps and dataset names, so legacy output does
 * not require new metadata. Malformed/incomplete collections are rejected.
 */
struct Collection {
    std::string header;
    std::vector<Frame> frames;
    explicit Collection(const std::string& xml) {
        const auto first = xml.find("<!-- Time ");
        if (first == std::string::npos || xml.find("</Xdmf>") == std::string::npos)
            throw std::runtime_error("Missing or incomplete XMF temporal collection");
        header = xml.substr(0, first);
        const std::regex timePattern("<Time[^>]*Value=\"([^\"]+)\"");
        const std::regex indexPattern("Iteration ([0-9]+)");
        const std::regex groupPattern(":/([^<\\s]+)/Values");
        std::size_t position = first;
        double previous = -std::numeric_limits<double>::infinity();
        int previousIndex = -1;
        while (position != std::string::npos) {
            auto end = xml.find("</Grid>", position);
            if (end == std::string::npos) throw std::runtime_error("Incomplete XMF frame");
            end += std::string("</Grid>").size();
            Frame frame;
            frame.xml = xml.substr(position, end - position) + "\n\n";
            std::smatch match;
            if (!std::regex_search(frame.xml, match, timePattern)) throw std::runtime_error("Missing XMF time");
            frame.time = std::stod(match[1]);
            if (!std::regex_search(frame.xml, match, indexPattern)) throw std::runtime_error("Missing XMF frame index");
            frame.index = std::stoi(match[1]);
            if (!std::isfinite(frame.time) || frame.time <= previous || frame.index <= previousIndex)
                throw std::runtime_error("XMF frames must have increasing times and indices");
            for (std::sregex_iterator i(frame.xml.begin(), frame.xml.end(), groupPattern), last; i != last; ++i)
                frame.groups.insert((*i)[1]);
            previous = frame.time;
            previousIndex = frame.index;
            frames.push_back(std::move(frame));
            position = xml.find("<!-- Time ", end);
        }
    }
    std::string through(double time) const {
        std::string result = header;
        for (const auto& frame : frames)
            if (frame.time <= time + 1.e-12) result += frame.xml;
        return result;
    }
};

} // namespace output
} // namespace FEDD
#endif
