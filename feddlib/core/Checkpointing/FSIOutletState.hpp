#ifndef FEDD_FSI_OUTLET_STATE_HPP
#define FEDD_FSI_OUTLET_STATE_HPP

#include "CheckpointFiles.hpp"
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace FEDD {

/** @brief Persistent pressure-boundary history at the start of an FSI timestep.
 * The previous flow rate is the value retained by the last pressure evaluation,
 * not the flow computed from the checkpoint's current velocity. Transition areas
 * retain their original capture time across restarts. Zero areas are allowed
 * only while the corresponding initialization/capture flag is false.
 */
struct FSIOutletState
{
    double initialInletArea = 0.;
    double initialOutletArea = 0.;
    double transitionOutletArea = 0.;
    double currentFlowRate = 0.;
    double previousFlowRate = 0.;
    double pressure = 0.;
    bool initialized = false;
    bool transitionCaptured = false;
};

namespace checkpoint {

inline std::string outletStateName(double time)
{
    return "FSIOutletState_" + std::to_string(time) + ".xml";
}

/** @brief Describe the pressure model settings that affect FSI continuation.
 * Defaults match the FE pressure-boundary routines. Solver and output options
 * deliberately do not participate in compatibility checks.
 */
inline Teuchos::ParameterList outletConfiguration(Teuchos::ParameterList& parameters)
{
    auto& fluid = parameters.sublist("Parameter Fluid");
    auto& general = parameters.sublist("General");
    const std::string model = fluid.get("Pressure Boundary Condition", "None");
    if (model != "None" && model != "Resistance" && model != "Absorbing" && model != "Absorbing Paper")
        throw std::logic_error("Unsupported FSI pressure boundary condition: " + model);
    Teuchos::ParameterList result;
    result.set("Model", model);
    if (model == "None") return result;
    result.set("State format version", 1);
    result.set("Average Flowrate", fluid.get("Average Flowrate", false));
    result.set("Flag Inlet Fluid", general.get("Flag Inlet Fluid", 4));
    result.set("Flag Outlet Fluid", general.get("Flag Outlet Fluid", 5));
    result.set("Normal Scale", fluid.get("Normal Scale", 1.));
    result.set("BC Ramp", fluid.get("BC Ramp", 0.1));
    if (model == "Resistance") {
        result.set("Resistance", fluid.get("Resistance", 1.));
        result.set("Viscosity", fluid.get("Viscosity", parameters.sublist("Parameter").get("Viscosity", 0.49)));
        result.set("Density", fluid.get("Density", 1.));
        result.set("Reference fluid pressure", fluid.get("Reference fluid pressure", 119.9));
    }
    else {
        const bool paper = model == "Absorbing Paper";
        result.set("Poisson Ratio", fluid.get("Poisson Ratio", 0.49));
        result.set("E", fluid.get("E", 12.));
        result.set("Wall thickness", fluid.get("Wall thickness", 0.0006));
        result.set("Density", fluid.get("Density", paper ? 1000. : 1.));
        result.set("Reference fluid pressure", fluid.get("Reference fluid pressure", 10666.));
        result.set("Max Ramp Time", fluid.get("Max Ramp Time", 0.1));
        result.set("Flowrate", fluid.get("Flowrate", paper ? 3. : 3.e-6));
        result.set(paper ? "Unsteady Start" : "Heart Beat Start",
                   fluid.get(paper ? "Unsteady Start" : "Heart Beat Start", 0.2));
    }
    return result;
}

/** @brief Validate scalar history without modifying the simulation.
 * An enabled model must have initialized history at nonzero checkpoint times.
 * Uncaptured transition areas are valid before the transition is evaluated.
 */
inline void validateOutletState(const FSIOutletState& state, const std::string& model, double time)
{
    for (double value : {state.initialInletArea, state.initialOutletArea, state.transitionOutletArea,
                         state.currentFlowRate, state.previousFlowRate, state.pressure})
        if (!std::isfinite(value)) throw std::runtime_error("FSI outlet state contains a nonfinite value");
    if (time > 0. && !state.initialized)
        throw std::runtime_error("FSI outlet history is not initialized at a nonzero checkpoint time");
    if (model != "Resistance" && state.initialized &&
        (state.initialInletArea <= 0. || state.initialOutletArea <= 0.))
        throw std::runtime_error("FSI absorbing boundary requires positive initial areas");
    if (state.transitionCaptured && (!state.initialized || state.transitionOutletArea <= 0.))
        throw std::runtime_error("FSI outlet transition area is invalid");
}

/** @brief Read and validate an outlet snapshot on the calling rank only.
 * The caller must propagate root I/O failures before broadcasting the result.
 * Missing state is an error even when legacy field-only restart is permitted.
 */
inline FSIOutletState readOutletState(const std::string& filename, const std::string& model, double time)
{
    if (!std::ifstream(filename)) throw std::runtime_error("Missing FSI outlet checkpoint " + filename);
    const auto data = Teuchos::getParametersFromXmlFile(filename);
    if (data->get<int>("Format version") != 1 || data->get<std::string>("Model") != model ||
        data->get<std::string>("Layout") != "legacy-start-of-step-v1" ||
        !std::isfinite(data->get<double>("Physical time")) ||
        std::abs(data->get<double>("Physical time") - time) >
            100. * std::numeric_limits<double>::epsilon() * std::max(1., std::abs(time)))
        throw std::runtime_error("Incompatible FSI outlet checkpoint " + filename);
    FSIOutletState state;
    state.initialInletArea = data->get<double>("Initial inlet area");
    state.initialOutletArea = data->get<double>("Initial outlet area");
    state.transitionOutletArea = data->get<double>("Transition outlet area");
    state.currentFlowRate = data->get<double>("Current flow rate");
    state.previousFlowRate = data->get<double>("Previous flow rate");
    state.pressure = data->get<double>("Outlet pressure");
    state.initialized = data->get<bool>("Initialized");
    state.transitionCaptured = data->get<bool>("Transition captured");
    validateOutletState(state, model, time);
    return state;
}

/** @brief Write a complete scalar snapshot on the calling rank only.
 * Teuchos serializes double parameters with 17 significant digits, preserving
 * their binary values on read-back. This file is separate from compatibility
 * metadata because it contains evolving simulation state.
 */
inline void writeOutletState(const FSIOutletState& state, const std::string& filename,
                             const std::string& model, double time)
{
    validateOutletState(state, model, time);
    Teuchos::ParameterList data;
    data.set("Format version", 1).set("Model", model).set("Physical time", time);
    data.set("Layout", "legacy-start-of-step-v1");
    data.set("Initial inlet area", state.initialInletArea).set("Initial outlet area", state.initialOutletArea);
    data.set("Transition outlet area", state.transitionOutletArea);
    data.set("Current flow rate", state.currentFlowRate).set("Previous flow rate", state.previousFlowRate);
    data.set("Outlet pressure", state.pressure).set("Initialized", state.initialized);
    data.set("Transition captured", state.transitionCaptured);
    std::ofstream output(filename);
    if (!output) throw std::runtime_error("Cannot write FSI outlet checkpoint " + filename);
    Teuchos::writeParameterListToXmlOStream(data, output);
    output.close();
    if (!output) throw std::runtime_error("Cannot finish FSI outlet checkpoint " + filename);
}

} // namespace checkpoint
} // namespace FEDD
#endif
