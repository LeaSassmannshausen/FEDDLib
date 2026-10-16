#ifndef FEDD_RECOVERY_CHECKPOINT_POLICY_HPP
#define FEDD_RECOVERY_CHECKPOINT_POLICY_HPP

#include <cmath>
#include <stdexcept>

namespace FEDD {
namespace checkpoint {

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

} // namespace checkpoint
} // namespace FEDD

#endif
