#ifndef FEDD_CHECKPOINT_TIME_STATE_HPP
#define FEDD_CHECKPOINT_TIME_STATE_HPP

#include "feddlib/core/FEDDCore.hpp"
#include <cmath>
#include <stdexcept>

namespace FEDD { namespace checkpoint {

/** @brief Actual clock/history information, independent of a future timestep schedule.
 * dt is the next planned increment; dt_prev is the last completed increment.
 * Explicit history times avoid reconstructing old states from a new input dt.
 * The underscore-prefixed lists are internal runtime data, not user options.
 */
inline const Teuchos::ParameterList* clockState(const ParameterListPtr_Type& parameters, double time)
{
    const auto& settings = parameters->sublist("Timestepping Parameter");
    for (const char* name : {"_Checkpoint time state", "_Restart time state"}) {
        if (settings.isSublist(name)) {
            const auto& state = settings.sublist(name);
            if (std::abs(state.get<double>("Physical time") - time) <= 1.e-12)
                return &state;
        }
    }
    return nullptr;
}

/// Return a saved history timestamp, falling back only for uniform legacy data.
inline double historyTime(const ParameterListPtr_Type& parameters, double time, int index, double fallbackDt)
{
    if (const auto* state = clockState(parameters, time)) {
        if (state->isSublist("History times"))
            return state->sublist("History times").get<double>(std::to_string(index));
        return time - index * state->get<double>("dt_prev");
    }
    return time - index * fallbackDt;
}

/// Increment of the last completed step, used to finalize structural derivatives.
inline double previousIncrement(const ParameterListPtr_Type& parameters, double time, double fallbackDt)
{
    if (const auto* state = clockState(parameters, time)) return state->get<double>("dt_prev");
    return fallbackDt;
}

} }
#endif
