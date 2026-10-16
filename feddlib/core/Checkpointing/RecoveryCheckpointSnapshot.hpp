#ifndef FEDD_RECOVERY_CHECKPOINT_SNAPSHOT_HPP
#define FEDD_RECOVERY_CHECKPOINT_SNAPSHOT_HPP

#include "CheckpointMetadata.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <sys/stat.h>
#include <unistd.h>

namespace FEDD {
namespace checkpoint {

// Create the directory for a recovery checkpoint, including any parent directories.
inline void makeRecoveryDirectory(const std::string& path) {
    if (path.empty() || path == ".") return;
    const auto slash = path.find_last_of('/');
    if (slash != std::string::npos && slash > 0) makeRecoveryDirectory(path.substr(0, slash));
    if (::mkdir(path.c_str(), 0777) != 0 && errno != EEXIST)
        throw std::runtime_error("Cannot create recovery directory " + path);
    struct stat info;
    if (::stat(path.c_str(), &info) != 0 || !S_ISDIR(info.st_mode))
        throw std::runtime_error("Recovery path is not a directory: " + path);
}

// Write a Teuchos::ParameterList to an XML file, throwing an exception if the write fails.
inline void writeRecoveryXml(const std::string& path, const Teuchos::ParameterList& data) {
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("Cannot write recovery information " + path);
    Teuchos::writeParameterListToXmlOStream(data, stream);
    stream.close();
    if (!stream) throw std::runtime_error("Cannot finish recovery information " + path);
}

/** @brief Independent copies of the fields and histories consumed by existing restart readers.
 * Capture components at their established checkpoint boundaries. No simulation
 * vector is accessed during writing; a failed solve may already have modified it.
 */
template<class SC, class LO, class GO, class NO>
class RecoveryCheckpointSnapshot {
public:
    using Vector = MultiVector<SC,LO,GO,NO>;
    using ConstVectorPtr = Teuchos::RCP<const Vector>;
    explicit RecoveryCheckpointSnapshot(double time) : time_(time) {}
    double time() const { return time_; }

    void addVector(const std::string& file, double time, ConstVectorPtr vector) {
        vectors_[file][std::to_string(std::max(0., time))] = Teuchos::rcp(new Vector(vector));
    }
    void addManifest(const Teuchos::ParameterList& schema, const Teuchos::ParameterList* clock = nullptr) {
        xml_[manifestName(schema, time_)] = atTime(schema, time_, clock);
    }
    long long stepNumber() const { return xml_.begin()->second.get<long long>("Step number"); }
    void addScalarWriter(const std::string& filename, std::function<void(const std::string&)> writer) {
        scalars_[filename] = std::move(writer);
    }

    /** @brief Reject incomplete captures before starting any collective HDF5 write. */
    void validate() const {
        if (xml_.empty()) throw std::logic_error("Recovery snapshot has no compatibility metadata");
        for (const auto& manifest : xml_) {
            const auto& required = manifest.second.sublist("Required datasets");
            for (auto it = required.begin(); it != required.end(); ++it) {
                const auto& item = required.sublist(required.name(it));
                const auto filename = item.get<std::string>("File");
                const auto file = vectors_.find(filename.substr(0, filename.size() - 3));
                if (file == vectors_.end() || !file->second.count(item.get<std::string>("Key")))
                    throw std::logic_error("Incomplete recovery snapshot: " + filename + "/" + item.get<std::string>("Key"));
            }
            if (manifest.second.isSublist("Required scalar state")) {
                const auto& requiredScalars = manifest.second.sublist("Required scalar state");
                for (auto it = requiredScalars.begin(); it != requiredScalars.end(); ++it)
                    if (!scalars_.count(requiredScalars.get<std::string>(requiredScalars.name(it))))
                        throw std::logic_error("Incomplete recovery scalar history");
            }
        }
    }

    /** @brief Write and close all component files, then validate their persisted datasets. */
    void write(const std::string& directory, const Teuchos::Comm<int>& comm) const {
        validate();
        onRoot(comm, [&] { makeRecoveryDirectory(directory); });
        for (const auto& file : vectors_) {
            HDF5Export<SC,LO,GO,NO> exporter(file.second.begin()->second->getMap(), joinPath(directory, file.first));
            for (const auto& dataset : file.second) exporter.writeVariablesHDF5(dataset.first, dataset.second);
            exporter.closeExporter();
        }
        onRoot(comm, [&] {
            for (const auto& entry : scalars_) entry.second(joinPath(directory, entry.first));
            for (const auto& entry : xml_) {
                writeRecoveryXml(joinPath(directory, entry.first), entry.second);
                const auto& required = entry.second.sublist("Required datasets");
                for (auto it = required.begin(); it != required.end(); ++it) {
                    const auto& item = required.sublist(required.name(it));
                    H5Handle h5(H5Fopen(joinPath(directory, item.get<std::string>("File")).c_str(),
                                       H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
                    if (h5 < 0) throw std::runtime_error("Cannot inspect recovery checkpoint");
                    validateVectorDataset(h5, item.get<std::string>("Key"),
                                          std::stoull(item.get<std::string>("Global DOFs")));
                }
            }
            Teuchos::ParameterList complete;
            complete.set("Complete", true).set("Physical time", time_);
            writeRecoveryXml(joinPath(directory, "Complete.xml"), complete);
        });
    }

    /** @brief Remove only files created by this snapshot, after a replacement is published. */
    void remove(const std::string& directory) const {
        for (const auto& entry : vectors_) std::remove(joinPath(directory, entry.first + ".h5").c_str());
        for (const auto& entry : scalars_) std::remove(joinPath(directory, entry.first).c_str());
        for (const auto& entry : xml_) std::remove(joinPath(directory, entry.first).c_str());
        std::remove(joinPath(directory, "Complete.xml").c_str());
        ::rmdir(directory.c_str());
    }

// Private members 
private:
    double time_;
    std::map<std::string, std::map<std::string, ConstVectorPtr>> vectors_;
    std::map<std::string, Teuchos::ParameterList> xml_;
    std::map<std::string, std::function<void(const std::string&)>> scalars_;
};

} // namespace checkpoint
} // namespace FEDD

#endif
