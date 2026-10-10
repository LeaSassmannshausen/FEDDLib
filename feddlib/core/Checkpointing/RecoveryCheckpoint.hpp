#ifndef FEDD_RECOVERY_CHECKPOINT_HPP
#define FEDD_RECOVERY_CHECKPOINT_HPP

#include "CheckpointMetadata.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include <cerrno>
#include <cstdio>
#include <map>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>

namespace FEDD { namespace checkpoint {

/** @brief Recovery criteria, independent of solver and model implementation.
 * Nonconvergence activates recovery even when cancellation is disabled. After
 * such a step the continued trajectory cannot replace the last reliable state.
 */
class RecoveryCheckpointPolicy {
public:
    explicit RecoveryCheckpointPolicy(double fraction = 0.5) : fraction_(fraction) {
        if (!std::isfinite(fraction) || fraction < 0. || fraction > 1.)
            throw std::logic_error("Recovery linear iteration fraction must be between zero and one");
    }
    template<class Report> bool triggered(const Report& report) const {
        return !report.reliable() || (report.maximumLinearIterations > 0 &&
            report.averageLinearIterations() > fraction_ * report.maximumLinearIterations);
    }
    template<class Report> void observe(const Report& report) {
        active_ = triggered(report);
        trusted_ = trusted_ && report.reliable();
    }
    bool active() const { return active_; }
    bool trusted() const { return trusted_; }
private:
    double fraction_;
    bool active_ = false;
    bool trusted_ = true;
};

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
private:
    double time_;
    std::map<std::string, std::map<std::string, ConstVectorPtr>> vectors_;
    std::map<std::string, Teuchos::ParameterList> xml_;
    std::map<std::string, std::function<void(const std::string&)>> scalars_;
};

/** @brief Own recovery generations and publish pointers only to complete snapshots.
 * The protected generation precedes the first difficult solve. The latest
 * generation follows reliable steps while the warning is active. Ignoring a
 * failed solve leaves both recovery generations protected from that trajectory.
 * Call collectively at known timestep boundaries, never from a signal handler.
 */
template<class SC, class LO, class GO, class NO>
class RecoveryCheckpointManager {
public:
    using Snapshot = RecoveryCheckpointSnapshot<SC,LO,GO,NO>;
    using SnapshotPtr = Teuchos::RCP<Snapshot>;
    RecoveryCheckpointManager(const ParameterListPtr_Type& parameters, const Teuchos::Comm<int>& comm)
        : parameters_(parameters), comm_(comm), policy_(parameters->sublist("Timestepping Parameter")
            .get("Recovery linear iteration fraction", 0.5)),
          enabled_(parameters->sublist("Timestepping Parameter").get("Failure recovery", false)),
          directory_(parameters->sublist("Timestepping Parameter").get("Recovery directory", std::string("recoveryCheckpoints"))) {}
    bool enabled() const { return enabled_; }

    /** @brief Register a complete snapshot of the state preceding the next attempted solve. */
    void beforeSolve(SnapshotPtr snapshot) {
        if (!enabled_ || !policy_.trusted()) return;
        snapshot->validate();
        lastValid_ = snapshot;
        const int interval = parameters_->sublist("Timestepping Parameter").get("Recovery interval steps", 1);
        if (interval < 0) throw std::logic_error("Recovery interval steps must be nonnegative");
        const auto step = snapshot->stepNumber();
        if (policy_.active() || (interval > 0 && step % interval == 0)) publishLatest(snapshot);
    }

    /** @brief Preserve recovery data before applying Cancel MaxNonLinIts.
     * @return True when the caller must stop. False permits the legacy continuation.
     */
    template<class Report> bool afterSolve(const Report& report, double attemptedTime) {
        if (!enabled_) return false;
        const bool trigger = policy_.triggered(report);
        if (trigger && protected_.empty()) {
            if (lastValid_.is_null()) throw std::logic_error("No valid state captured for recovery");
            protected_ = publish(lastValid_);
        }
        if (trigger) publishLatest(lastValid_);
        policy_.observe(report);
        const bool stop = !report.reliable() && parameters_->sublist("Parameter").get("Cancel MaxNonLinIts", false);
        if (trigger || !policy_.trusted()) {
            Teuchos::ParameterList info;
            info.set("Attempted time", attemptedTime).set("Converged", report.converged);
            info.set("Linear solves converged", report.linearSolvesConverged).set("Finite", report.finite);
            info.set("Nonlinear iterations", report.nonlinearIterations).set("Linear solve count", report.linearSolveCount);
            info.set("Maximum nonlinear iterations", parameters_->sublist("Parameter").get("MaxNonLinIts", 10));
            info.set("Trigger reason", std::string(!trigger ? "Unreliable continuation after earlier failure" :
                !report.finite ? "Nonfinite state" :
                !report.linearSolvesConverged ? "Linear nonconvergence" :
                !report.converged ? "Nonlinear nonconvergence" : "Average linear iteration count"));
            info.set("Total linear iterations", report.totalLinearIterations);
            info.set("Average linear iterations", report.averageLinearIterations());
            info.set("Maximum linear iterations", report.maximumLinearIterations).set("Final criterion", report.finalCriterion);
            info.set("Stopped", stop).set("Trajectory reliable", policy_.trusted());
            info.set("Protected directory", protected_).set("Latest directory", latest_);
            info.set("Recovery time", lastValid_->time());
            onRoot(comm_, [&] {
                makeRecoveryDirectory(directory_);
                writeRecoveryXml(joinPath(directory_, "RecoveryStatus.xml.tmp"), info);
                if (std::rename(joinPath(directory_, "RecoveryStatus.xml.tmp").c_str(),
                                joinPath(directory_, "RecoveryStatus.xml").c_str()) != 0)
                    throw std::runtime_error("Cannot publish recovery status");
            });
            if (comm_.getRank() == 0)
                std::cout << "Recovery: timestep " << attemptedTime << ", nonlinear iterations " << report.nonlinearIterations
                    << ", average linear iterations " << report.averageLinearIterations()
                    << (stop ? ": stopping" : ": continuing") << ". Reliable restart state: " << latest_ << std::endl;
        }
        return stop;
    }

    /** @brief Publish an accepted-step BDF state immediately after completing the step. */
    void acceptedState(SnapshotPtr snapshot) {
        if (enabled_ && policy_.trusted()) beforeSolve(snapshot);
    }
private:
    std::string publish(SnapshotPtr snapshot) {
        const std::string name = "generation_" + std::to_string(generation_++);
        const std::string target = joinPath(directory_, name);
        const std::string temporary = target + ".tmp";
        onRoot(comm_, [&] {
            makeRecoveryDirectory(directory_);
            struct stat info;
            if (::stat(target.c_str(), &info) == 0 || ::stat(temporary.c_str(), &info) == 0)
                throw std::runtime_error("Recovery generation already exists: use a new Recovery directory");
        });
        snapshot->write(temporary, comm_);
        onRoot(comm_, [&] {
            if (std::rename(temporary.c_str(), target.c_str()) != 0)
                throw std::runtime_error("Cannot publish recovery generation " + target);
        });
        return target;
    }
    void publishLatest(SnapshotPtr snapshot) {
        if (!latestSnapshot_.is_null() && latestSnapshot_->time() == snapshot->time()) {
            writeLatestIndex(latest_, snapshot->time());
            return;
        }
        const auto next = publish(snapshot);
        writeLatestIndex(next, snapshot->time());
        onRoot(comm_, [&] {
            if (!latest_.empty() && latest_ != protected_) latestSnapshot_->remove(latest_);
        });
        latest_ = next;
        latestSnapshot_ = snapshot;
    }
    void writeLatestIndex(const std::string& directory, double time) {
        onRoot(comm_, [&] {
            Teuchos::ParameterList index;
            index.set("Directory", directory).set("Physical time", time).set("Protected directory", protected_);
            writeRecoveryXml(joinPath(directory_, "Latest.xml.tmp"), index);
            if (std::rename(joinPath(directory_, "Latest.xml.tmp").c_str(), joinPath(directory_, "Latest.xml").c_str()) != 0)
                throw std::runtime_error("Cannot publish latest recovery checkpoint");
        });
    }
    ParameterListPtr_Type parameters_;
    const Teuchos::Comm<int>& comm_;
    RecoveryCheckpointPolicy policy_;
    bool enabled_;
    std::string directory_, protected_, latest_;
    unsigned generation_ = 0;
    SnapshotPtr lastValid_, latestSnapshot_;
};
} } // namespace FEDD::checkpoint
#endif
