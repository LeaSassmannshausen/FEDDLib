#ifndef FEDD_CHECKPOINT_METADATA_HPP
#define FEDD_CHECKPOINT_METADATA_HPP

#include "CheckpointFiles.hpp"
#include "feddlib/core/General/HDF5VectorInfo.hpp"
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_TestForException.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace FEDD {
namespace checkpoint {

/// Broadcast a root-only I/O failure before any rank proceeds to collective reads.
inline void onRoot(const Teuchos::Comm<int>& comm, const std::function<void()>& operation)
{
    std::string error;
    if (comm.getRank() == 0) {
        try { operation(); }
        catch (const std::exception& exception) { error = exception.what(); }
        catch (...) { error = "Unknown checkpoint I/O error"; }
    }
    int size = static_cast<int>(error.size());
    Teuchos::broadcast(comm, 0, 1, &size);
    error.resize(size);
    if (size) {
        Teuchos::broadcast(comm, 0, size, &error[0]);
        throw std::runtime_error("Restart compatibility error: " + error);
    }
}

/** @brief Problem-independent description of one registered checkpoint field.
 * Global counts and index bases use strings to preserve the map ordinal range
 * in XML. The caller computes the fingerprint from the reference mesh.
 */
struct FieldDescription
{
    std::string name;
    std::string feType;
    int dimension;
    int components;
    std::string globalDofs;
    std::string indexBase;
    std::string referenceFingerprint;
};

/** @brief Build the versioned compatibility schema from settings and fields.
 * Owns the metadata keys and supported integration/history rules, independently
 * of Problem. Field order is retained because it identifies component manifests.
 *
 * @param[in,out] parameters Simulation settings; missing defaults are inserted
 *                          using the usual Teuchos ParameterList convention.
 * @param[in] fields Registered fields, described before any ALE mesh motion.
 * @param[in] role Empty for the main problem, or an FSI component role.
 * @throws std::logic_error For unsupported integration settings.
 */
inline Teuchos::ParameterList makeSchema(Teuchos::ParameterList& parameters,
                                        const std::vector<FieldDescription>& fields,
                                        const std::string& role = "")
{
    auto& time = parameters.sublist("Timestepping Parameter");
    TEUCHOS_TEST_FOR_EXCEPTION(time.sublist("Timestepping Intervalls").get("Number of Segments", 0) != 0,
        std::logic_error, "Checkpoint metadata version 1 supports fixed timesteps only.");
    const double dt = time.get("dt", 0.01);
    TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(dt) || dt <= 0., std::logic_error,
                               "Checkpoint metadata requires a finite positive dt.");
    const bool fsi = parameters.sublist("Parameter").get("FSI", false);
    const std::string method = role == "FSI structure" ? "Newmark" : time.get("Class", "Multistep");
    TEUCHOS_TEST_FOR_EXCEPTION(method != "Multistep" && method != "Newmark", std::logic_error,
                               "Checkpoint metadata does not support integration class " << method);
    Teuchos::ParameterList schema;
    schema.set("Format version", 1).set("FSI", fsi).set("Role", role);
    auto& integration = schema.sublist("Integration");
    integration.set("Class", method).set("dt", dt);
    integration.set("Layout", fsi || method == "Newmark" || role == "FSI fluid"
                     ? "legacy-start-of-step-v1" : "accepted-step-bdf-v1");
    int history = 2;
    if (method == "Multistep") {
        const int order = time.get("BDF", 1);
        TEUCHOS_TEST_FOR_EXCEPTION(order != 1 && order != 2, std::logic_error,
                                   "Checkpoint metadata supports BDF1 and BDF2 only.");
        const bool extrapolation = parameters.sublist("General").get("Linearization", "FixedPoint") == "Extrapolation";
        history = extrapolation ? std::max(order, 2) : order;
        integration.set("BDF", order).set("Extrapolation", extrapolation);
    }
    if (method == "Newmark" || fsi) {
        integration.set("beta", time.get("beta", 0.25));
        integration.set("gamma", time.get("gamma", 0.5));
    }
    integration.set("Solution history", history);
    if (fsi) integration.set("Geometry Explicit", parameters.sublist("Parameter").get("Geometry Explicit", true));
    auto& descriptions = schema.sublist("Fields");
    for (const auto& field : fields) {
        descriptions.sublist(field.name)
            .set("FE type", field.feType).set("Dimension", field.dimension)
            .set("Components", field.components).set("Global DOFs", field.globalDofs)
            .set("Index base", field.indexBase)
            .set("Reference mesh and DOF fingerprint", field.referenceFingerprint);
    }
    return schema;
}

/// Component names prevent FSI fluid/structure manifests from overwriting each other.
inline std::string manifestName(const Teuchos::ParameterList& schema, double time)
{
    std::string name = "Checkpoint";
    const auto& fields = schema.sublist("Fields");
    for (auto it = fields.begin(); it != fields.end(); ++it)
        name += "_" + fields.name(it);
    return name + "_" + std::to_string(time) + ".xml";
}

/// Compare only compatibility information, deliberately excluding solver/output settings.
inline void compare(const Teuchos::ParameterList& saved, const Teuchos::ParameterList& expected,
                    const std::string& path = "")
{
    for (auto it = expected.begin(); it != expected.end(); ++it) {
        const std::string name = expected.name(it);
        const std::string key = path.empty() ? name : path + "/" + name;
        if (!saved.isParameter(name) && !saved.isSublist(name))
            throw std::runtime_error(key + ": missing metadata");
        const auto& entry = expected.entry(it);
        if (entry.isList()) {
            if (!saved.isSublist(name)) throw std::runtime_error(key + ": expected a sublist");
            compare(saved.sublist(name), expected.sublist(name), key);
        }
        else if (entry.isType<double>()) {
            if (!saved.isType<double>(name)) throw std::runtime_error(key + ": invalid type");
            const double a = saved.get<double>(name), b = expected.get<double>(name);
            if (!std::isfinite(a) || !std::isfinite(b) ||
                std::abs(a - b) > 100. * std::numeric_limits<double>::epsilon() * std::max(std::abs(a), std::abs(b)))
                throw std::runtime_error(key + ": saved " + std::to_string(a) + ", expected " + std::to_string(b));
        }
        else if (!(saved.getEntry(name) == entry)) {
            std::ostringstream message;
            message << key << ": saved " << saved.getEntry(name) << ", expected " << entry;
            throw std::runtime_error(message.str());
        }
    }
    for (auto it = saved.begin(); it != saved.end(); ++it) {
        const std::string name = saved.name(it);
        if (!expected.isParameter(name) && !expected.isSublist(name))
            throw std::runtime_error(path + "/" + name + ": unexpected metadata entry");
    }
}

/** @brief Describe the exact files and time keys consumed by the current restore routines.
 * Version 1 retains decimal HDF5 keys and the existing Newmark/FSI clock conventions.
 * BDF startup needs only the history available since time zero.
 */
inline Teuchos::ParameterList atTime(const Teuchos::ParameterList& schema, double time)
{
    Teuchos::ParameterList result(schema);
    const auto& integration = schema.sublist("Integration");
    const double dt = integration.get<double>("dt");
    if (!std::isfinite(dt) || dt <= 0. || !std::isfinite(time) || time < 0. ||
        !std::isfinite(time / dt) || time / dt >= static_cast<double>(std::numeric_limits<long long>::max()))
        throw std::runtime_error("Checkpoint requires a finite nonnegative time and positive fixed dt");
    const auto step = std::llround(time / dt);
    if (std::abs(time / dt - step) > 100. * std::numeric_limits<double>::epsilon() * std::max(1., time / dt))
        throw std::runtime_error("Checkpoint metadata version 1 requires a time on the fixed-dt grid");
    result.set("Physical time", time);
    result.set("Step number", static_cast<long long>(step));
    auto& required = result.sublist("Required datasets");
    const auto add = [&](const std::string& file, const std::string& field, double stateTime) {
        if (stateTime < -1.e-12) return;
        stateTime = std::max(0., stateTime);
        const std::string key = std::to_string(stateTime);
        auto& item = required.sublist(file + "/" + key);
        item.set("File", file + ".h5").set("Key", key).set("Field", field);
        item.set("Global DOFs", schema.sublist("Fields").sublist(field).get<std::string>("Global DOFs"));
    };
    const bool fsi = schema.get<bool>("FSI");
    const std::string role = schema.get<std::string>("Role");
    const bool newmark = integration.get<std::string>("Class") == "Newmark";
    const int history = integration.get<int>("Solution history");
    const auto& fields = schema.sublist("Fields");
    for (auto it = fields.begin(); it != fields.end(); ++it) {
        const std::string field = fields.name(it);
        if (field == "d_f" && fsi) {
            add("Solutiond_f", field, time);
            continue;
        }
        // A structure subproblem imports primary d_s from the coupled Solution file.
        add("Solution" + field, field, time);
        if (!newmark) {
            for (int j = 1; j < history; ++j) add("Solution" + field, field, time - j * dt);
        }
        if (role == "FSI fluid")
            for (int j = 0; j < integration.get<int>("BDF"); ++j)
                add("Rhs" + field, field, time - j * dt);
        if (newmark || (fsi && field == "d_s")) {
            add("SolutionNewmark" + field, field, time);
            add("SolutionNewmark" + field, field, time - dt);
            add("ds_Velocity", field, time);
            add("ds_Acceleration", field, time);
        }
    }
    if (fsi) {
        // Fluid subproblem names differ from the monolithic velocity name.
        for (const std::string field : {std::string("u_f"), std::string("p")}) {
            const std::string fluidName = field == "u_f" ? "u" : "p";
            for (int j = 0; j < history; ++j) add("Solution" + fluidName, field, time - j * dt);
            for (int j = 0; j < integration.get<int>("BDF"); ++j)
                add("Rhs" + fluidName, field, time - j * dt);
        }
    }
    return result;
}

/** @brief Validate a manifest and every required dataset before simulation values change.
 * Missing manifests are accepted only with explicit "Allow legacy restart".
 * Legacy mode still inspects all required HDF5 fields and history, but cannot prove mesh identity.
 */
inline void validate(const Teuchos::ParameterList& schema, const ParameterListPtr_Type& parameters,
                     double time, const Teuchos::Comm<int>& comm)
{
    const auto expected = atTime(schema, time);
    onRoot(comm, [&] {
        const std::string file = restartFile(parameters, manifestName(schema, time));
        std::ifstream input(file);
        if (!input) {
            if (!parameters->sublist("Timestepping Parameter").get("Allow legacy restart", false))
                throw std::runtime_error("Missing manifest " + file +
                    ". For an old checkpoint explicitly set 'Allow legacy restart' to true.");
            std::cerr << "Legacy restart: no manifest; mesh/integration compatibility cannot be verified.\n";
        }
        else {
            std::ostringstream xml;
            xml << input.rdbuf();
            compare(*Teuchos::getParametersFromXmlString(xml.str()), expected);
        }
        const auto& datasets = expected.sublist("Required datasets");
        for (auto it = datasets.begin(); it != datasets.end(); ++it) {
            const auto& item = datasets.sublist(datasets.name(it));
            const std::string filename = restartFile(parameters, item.get<std::string>("File"));
            if (!std::ifstream(filename)) throw std::runtime_error("Missing checkpoint file " + filename);
            H5Handle h5(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
            if (h5 < 0) throw std::runtime_error("Cannot open checkpoint file " + filename);
            try {
                validateVectorDataset(h5, item.get<std::string>("Key"),
                    std::stoull(item.get<std::string>("Global DOFs")));
            }
            catch (const std::exception& exception) {
                throw std::runtime_error(filename + ": " + exception.what());
            }
        }
    });
}

/// Write descriptive metadata; this is not an atomic checkpoint completion marker.
inline void write(const Teuchos::ParameterList& schema, const ParameterListPtr_Type& parameters,
                  double time, const Teuchos::Comm<int>& comm)
{
    const auto description = atTime(schema, time);
    onRoot(comm, [&] {
        const std::string filename = checkpointFile(parameters, manifestName(schema, time));
        std::ofstream output(filename);
        if (!output) throw std::runtime_error("Cannot write checkpoint manifest " + filename);
        Teuchos::writeParameterListToXmlOStream(description, output);
        output.close();
        if (!output) throw std::runtime_error("Cannot finish checkpoint manifest " + filename);
    });
}

} // namespace checkpoint
} // namespace FEDD
#endif
