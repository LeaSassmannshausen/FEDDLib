#ifndef FEDD_CHECKPOINT_TIMING_HPP
#define FEDD_CHECKPOINT_TIMING_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <memory>
#include <string>

namespace FEDD {
namespace checkpoint {

/** @brief Optional wall-clock measurements of the production checkpoint paths.
 * Enable General/Checkpoint timings to register Teuchos counters named
 * "FEDD checkpoint - <operation>". Timers are local and introduce no MPI
 * barriers; a benchmark can reduce their elapsed times across ranks afterwards.
 * Nested counters are inclusive and must not be added to their parent counters.
 * No counters are registered when the option is disabled (the default).
 */
class ScopedTimer {
public:
    ScopedTimer(const Teuchos::RCP<Teuchos::ParameterList>& parameters, const std::string& operation,
                bool active = true)
    {
        if (!active || parameters.is_null() ||
            !parameters->sublist("General").get("Checkpoint timings", false))
            return;
        const std::string name = "FEDD checkpoint - " + operation;
        auto counter = Teuchos::TimeMonitor::lookupCounter(name);
        if (counter.is_null()) counter = Teuchos::TimeMonitor::getNewCounter(name);
        monitor_.reset(new Teuchos::TimeMonitor(*counter));
    }

private:
    std::unique_ptr<Teuchos::TimeMonitor> monitor_;
};

} // namespace checkpoint
} // namespace FEDD
#endif
