#ifndef FEDD_RECOVERY_CHECKPOINT_MANAGER_HPP
#define FEDD_RECOVERY_CHECKPOINT_MANAGER_HPP

#include "RecoveryCheckpointPolicy.hpp"
#include "RecoveryCheckpointSnapshot.hpp"
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>

namespace FEDD {
namespace checkpoint {

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

} // namespace checkpoint
} // namespace FEDD

#endif
