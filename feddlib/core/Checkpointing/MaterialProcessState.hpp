#ifndef FEDD_MATERIAL_PROCESS_STATE_HPP
#define FEDD_MATERIAL_PROCESS_STATE_HPP

#include <Teuchos_ParameterList.hpp>
#include <stdexcept>
#include <string>

namespace FEDD { namespace checkpoint {

// Monotone markers for committed material history. This is deliberately not
// the Gauss-point history needed to restart after one of these processes.
struct MaterialProcessState {
    bool active = false;
    bool growth = false;
    bool reorientation = false;

    void record(bool a, bool g, bool r) {
        active = active || a;
        growth = growth || g;
        reorientation = reorientation || r;
    }
    Teuchos::ParameterList metadata() const {
        Teuchos::ParameterList result;
        result.set("Active response took place", active)
              .set("Growth took place", growth)
              .set("Reorientation took place", reorientation);
        return result;
    }
};

inline bool tracksMaterialProcesses(const Teuchos::ParameterList& parameters)
{
    const auto& controls = parameters.sublist("Parameter");
    const auto& solid = parameters.isSublist("Parameter Solid") ? parameters.sublist("Parameter Solid") : controls;
    const std::string fallback = solid.isParameter("Material model") ? solid.get<std::string>("Material model") :
        (controls.isParameter("Material model") ? controls.get<std::string>("Material model") : "");
    const std::string model = controls.isParameter("Structure Model") ? controls.get<std::string>("Structure Model") :
        (solid.isParameter("Structure Model") ? solid.get<std::string>("Structure Model") : fallback);
    return model == "SCI_SMC_Active_Growth_Reorientation" || model == "SCI_SMC_MLCK";
}

inline void requirePassiveMaterialRestart(const Teuchos::ParameterList& saved)
{
    if (!saved.isSublist("Material process state"))
        throw std::runtime_error("Missing material process state: this checkpoint cannot establish that material history is unchanged. Recreate a checkpoint with process tracking enabled.");
    const auto& state = saved.sublist("Material process state");
    for (const char* name : {"Active response took place", "Growth took place", "Reorientation took place"}) {
        if (!state.isType<bool>(name))
            throw std::runtime_error(std::string("Missing or invalid material process marker: ") + name);
        if (state.get<bool>(name))
            throw std::runtime_error(std::string("Cannot restart: ") + name +
                ". Restoring AceGen material history is not implemented; restart from a checkpoint before these processes.");
    }
}

} }
#endif
