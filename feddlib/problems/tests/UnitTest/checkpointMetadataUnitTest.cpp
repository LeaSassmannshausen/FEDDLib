#include "feddlib/core/Checkpointing/CheckpointMetadata.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <filesystem>
#include <iostream>
#include <map>

using namespace FEDD;

// Small real HDF5 fixtures exercise validation, independently of a PDE solve.
int main(int argc, char** argv)
{
    Teuchos::GlobalMPISession session(&argc, &argv);
    auto comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    auto parameters = Teuchos::rcp(new Teuchos::ParameterList);
    const std::string directory = "checkpointMetadataUnitTest-data";
    parameters->sublist("Timestepping Parameter")
        .set("Checkpoint directory", directory).set("Restart directory", directory)
        .set("Initial solution directory", directory);
    checkpoint::onRoot(*comm, [&] {
        std::filesystem::remove_all(directory);
        std::filesystem::create_directory(directory);
    });
    Teuchos::ParameterList schema;
    schema.set("Format version", 1).set("FSI", false).set("Role", std::string(""));
    schema.sublist("Integration").set("Class", std::string("Multistep"))
        .set("dt", 0.0025).set("Layout", std::string("accepted-step-bdf-v1"))
        .set("BDF", 2).set("Extrapolation", false).set("Solution history", 2);
    for (const std::string field : {std::string("u"), std::string("p")})
        schema.sublist("Fields").sublist(field)
            .set("FE type", std::string("P1")).set("Dimension", 2)
            .set("Components", field == "u" ? 2 : 1)
            .set("Global DOFs", std::string(field == "u" ? "8" : "4"))
            .set("Index base", std::string("0"))
            .set("Reference mesh and DOF fingerprint", std::string("mesh-A"));
    const double time = 0.01;
    const auto makeMap = [&](int count) -> Teuchos::RCP<const Map<>> {
        const int local = count / comm->getSize() + (comm->getRank() < count % comm->getSize() ? 1 : 0);
        return Teuchos::rcp(new Map<default_lo,default_go,default_no>(count, local, 0, comm));
    };
    const auto fixture = [&] {
        for (const std::string field : {std::string("u"), std::string("p")}) {
            auto map = makeMap(field == "u" ? 8 : 4);
            auto values = Teuchos::rcp(new MultiVector<>(map));
            values->putScalar(1.);
            HDF5Export<> output(map, checkpointFile(parameters, "Solution" + field));
            output.writeVariablesHDF5(std::to_string(time), values);
            output.writeVariablesHDF5(std::to_string(time - 0.0025), values);
            output.closeExporter();
        }
        checkpoint::write(schema, parameters, time, *comm);
    };
    int failures = 0;
    const auto rejects = [&](const std::string& label, const std::string& diagnostic,
                             const std::function<void()>& operation) {
        int matched = 0;
        try { operation(); }
        catch (const std::exception& exception) {
            matched = std::string(exception.what()).find(diagnostic) != std::string::npos;
            if (!matched && comm->getRank() == 0) std::cerr << exception.what() << '\n';
        }
        int allMatched = 0;
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &matched, &allMatched);
        if (!allMatched) ++failures;
        if (comm->getRank() == 0)
            std::cout << (allMatched ? "PASS " : "FAIL ") << label << '\n';
    };
    {
        auto resumed = Teuchos::rcp(new Teuchos::ParameterList(*parameters));
        auto& timeSettings = resumed->sublist("Timestepping Parameter");
        timeSettings.set("Restart", true).set("Checkpointing", true);
        rejects("checkpoint output preserves restart input", "must be different", [&] {
            checkpoint::validateRestartOutputDirectory(resumed, *comm);
        });
        timeSettings.set("Checkpoint directory", directory + "/../" + directory);
        rejects("relative alias of restart directory", "must be different", [&] {
            checkpoint::validateRestartOutputDirectory(resumed, *comm);
        });
        checkpoint::onRoot(*comm, [&] {
            std::filesystem::create_directory_symlink(".", directory + "/alias");
        });
        timeSettings.set("Checkpoint directory", directory + "/alias");
        rejects("symlink alias of restart directory", "must be different", [&] {
            checkpoint::validateRestartOutputDirectory(resumed, *comm);
        });
        timeSettings.set("Checkpointing", false);
        resumed->sublist("General").set("Safe all solution", true);
        rejects("save-all output preserves restart input", "must be different", [&] {
            checkpoint::validateRestartOutputDirectory(resumed, *comm);
        });
        timeSettings.set("Checkpoint directory", directory + "/resumed");
        checkpoint::validateRestartOutputDirectory(resumed, *comm);
        timeSettings.set("Checkpoint directory", directory);
        resumed->sublist("General").set("Safe all solution", false);
        checkpoint::validateRestartOutputDirectory(resumed, *comm);
        timeSettings.set("Restart", false).set("Checkpointing", true);
        checkpoint::validateRestartOutputDirectory(resumed, *comm);
        if (comm->getRank() == 0)
            std::cout << "PASS separate output, read-only restart and fresh checkpoint run\n";
    }
    const auto validate = [&](const Teuchos::ParameterList& required) {
        checkpoint::validate(required, parameters, time, *comm);
    };
    const auto validateInitial = [&](const Teuchos::ParameterList& required) {
        checkpoint::validateInitialSolution(required, parameters, time, *comm);
    };
    // Independently specified version-1 schemas must remain compatible after
    // moving schema construction out of Problem, including component roles.
    const auto checksSchema = [&](const std::string& label, const Teuchos::ParameterList& expected,
                                  Teuchos::ParameterList settings,
                                  const std::vector<checkpoint::FieldDescription>& fields) {
        int matched = 1;
        try {
            const auto actual = checkpoint::makeSchema(settings, fields, expected.get<std::string>("Role"));
            auto current = expected;
            current.set("Format version", 2);
            checkpoint::compare(current, actual);
            checkpoint::compare(checkpoint::atTime(current, time), checkpoint::atTime(actual, time));
            if (checkpoint::manifestName(expected, time) != checkpoint::manifestName(actual, time))
                throw std::runtime_error("Changed component manifest name/field order");
        }
        catch (const std::exception& exception) {
            matched = 0;
            if (comm->getRank() == 0) std::cerr << exception.what() << '\n';
        }
        int allMatched = 0;
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &matched, &allMatched);
        if (!allMatched) ++failures;
        if (comm->getRank() == 0)
            std::cout << (allMatched ? "PASS " : "FAIL ") << label << '\n';
    };
    Teuchos::ParameterList settings;
    settings.sublist("Timestepping Parameter").set("Class", std::string("Multistep"))
        .set("dt", 0.0025).set("BDF", 2);
    const std::vector<checkpoint::FieldDescription> fields = {
        {"u", "P1", 2, 2, "8", "0", "mesh-A"},
        {"p", "P1", 2, 1, "4", "0", "mesh-A"}
    };
    checksSchema("version-1 BDF schema", schema, settings, fields);
    {
        auto bdf1 = schema;
        bdf1.sublist("Integration").set("BDF", 1).set("Solution history", 1);
        auto bdf1Settings = settings;
        bdf1Settings.sublist("Timestepping Parameter").set("BDF", 1);
        checksSchema("version-1 BDF1 schema", bdf1, bdf1Settings, fields);
        bdf1.sublist("Integration").set("Extrapolation", true).set("Solution history", 2);
        bdf1Settings.sublist("General").set("Linearization", std::string("Extrapolation"));
        checksSchema("version-1 BDF1 extrapolation history", bdf1, bdf1Settings, fields);
        auto fluid = schema;
        fluid.set("Role", std::string("FSI fluid"));
        fluid.sublist("Integration").set("Layout", std::string("legacy-start-of-step-v1"));
        checksSchema("version-1 FSI fluid schema", fluid, settings, fields);
    }
    fixture();
    validate(schema);
    if (comm->getRank() == 0) std::cout << "PASS compatible checkpoint\n";
    {
        auto changed = schema;
        // Source time need not be on the new dt grid, nor match its BDF order.
        changed.sublist("Integration").set("dt", 0.003).set("BDF", 1).set("Solution history", 1);
        validateInitial(changed);
        if (comm->getRank() == 0) std::cout << "PASS initial fields with new integration settings\n";
        changed.set("Format version", 999);
        rejects("initial solution unknown format", "Format version", [&] { validateInitial(changed); });
        rejects("nonfinite initial source time", "finite and nonnegative", [&] {
            checkpoint::validateInitialSolution(schema, parameters, std::numeric_limits<double>::infinity(), *comm);
        });
    }
    for (const std::string name : {std::string("Reference mesh and DOF fingerprint"), std::string("FE type"), std::string("Global DOFs")}) {
        auto changed = schema;
        changed.sublist("Fields").sublist("u").set(name, std::string("incompatible"));
        rejects("incompatible " + name, name, [&] { validate(changed); });
        rejects("initial solution incompatible " + name, name, [&] { validateInitial(changed); });
    }
    {
        auto changed = schema;
        changed.sublist("Integration").set("BDF", 1);
        rejects("changed BDF order", "BDF", [&] { validate(changed); });
        changed = schema;
        changed.sublist("Integration").set("dt", 0.005);
        validate(changed);
        if (parameters->sublist("Timestepping Parameter").sublist("_Restart time state").get<double>("dt_prev") != 0.0025) ++failures;
        if (comm->getRank() == 0) std::cout << "PASS changed timestep restores version-1 source increment\n";
        changed = schema;
        changed.set("Format version", 999);
        rejects("unknown version", "Format version", [&] { validate(changed); });
    }
    const auto editVelocity = [&](const std::function<void(hid_t)>& edit) {
        checkpoint::onRoot(*comm, [&] {
            checkpoint::H5Handle file(H5Fopen((directory + "/Solutionu.h5").c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
            if (file < 0) throw std::runtime_error("Cannot edit fixture");
            edit(file);
        });
    };
    editVelocity([&](hid_t file) { H5Ldelete(file, std::to_string(time - 0.0025).c_str(), H5P_DEFAULT); });
    rejects("missing BDF history", "missing required history/field group", [&] { validate(schema); });
    validateInitial(schema);
    if (comm->getRank() == 0) std::cout << "PASS initial solution needs no BDF history\n";
    editVelocity([&](hid_t file) { H5Ldelete(file, std::to_string(time).c_str(), H5P_DEFAULT); });
    rejects("initial solution missing primary field", "missing required history/field group", [&] { validateInitial(schema); });
    fixture();
    for (const std::string property : {std::string("GlobalLength"), std::string("NumVectors")}) {
        editVelocity([&](hid_t file) {
            checkpoint::H5Handle dataset(H5Dopen(file, (std::to_string(time) + "/" + property).c_str(), H5P_DEFAULT), H5Dclose);
            int wrong = 99;
            H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &wrong);
        });
        rejects("inconsistent " + property, property, [&] { validate(schema); });
        rejects("initial solution inconsistent " + property, property, [&] { validateInitial(schema); });
        fixture();
    }
    editVelocity([&](hid_t file) {
        const std::string path = std::to_string(time) + "/Values";
        H5Ldelete(file, path.c_str(), H5P_DEFAULT);
        hsize_t shape[2] = {2, 4}; // Same total entries, incompatible vector layout.
        checkpoint::H5Handle space(H5Screate_simple(2, shape, nullptr), H5Sclose);
        checkpoint::H5Handle dataset(H5Dcreate(file, path.c_str(), H5T_NATIVE_DOUBLE, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
    });
    rejects("inconsistent Values shape", "Values dimensions", [&] { validate(schema); });
    fixture();
    editVelocity([&](hid_t file) {
        const std::string path = std::to_string(time) + "/Values";
        H5Ldelete(file, path.c_str(), H5P_DEFAULT);
        hsize_t shape[2] = {1, 8};
        checkpoint::H5Handle space(H5Screate_simple(2, shape, nullptr), H5Sclose);
        checkpoint::H5Handle dataset(H5Dcreate(file, path.c_str(), H5T_NATIVE_INT, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
    });
    rejects("unsupported datatype", "64-bit floating-point", [&] { validate(schema); });
    fixture();
    rejects("reader destination map mismatch", "GlobalLength", [&] {
        auto wrongMap = makeMap(7);
        HDF5Import<> input(wrongMap, directory + "/Solutionu");
        input.readVariablesHDF5(std::to_string(time));
    });
    checkpoint::onRoot(*comm, [&] {
        std::filesystem::remove(directory + "/" + checkpoint::manifestName(schema, time));
    });
    rejects("missing manifest", "Missing manifest", [&] { validate(schema); });
    rejects("initial solution missing manifest", "Missing initial solution manifest", [&] { validateInitial(schema); });
    parameters->sublist("Timestepping Parameter").set("Allow legacy restart", true);
    validate(schema);
    rejects("initial solution cannot bypass mesh checks through legacy mode", "Missing initial solution manifest", [&] { validateInitial(schema); });
    if (comm->getRank() == 0) std::cout << "PASS explicit legacy restart\n";
    // Legacy mode must not bypass dataset checks.
    editVelocity([&](hid_t file) { H5Ldelete(file, std::to_string(time - 0.0025).c_str(), H5P_DEFAULT); });
    rejects("legacy restart missing history", "missing required history/field group", [&] { validate(schema); });
    // Newmark and coupled FSI add distinct auxiliary data to the same validator.
    parameters->sublist("Timestepping Parameter").set("Allow legacy restart", false);
    const auto fixtureFor = [&](const Teuchos::ParameterList& selected) {
        const auto description = checkpoint::atTime(selected, time, checkpoint::clockState(parameters, time));
        const auto& required = description.sublist("Required datasets");
        std::map<std::string, std::vector<Teuchos::ParameterList>> files;
        for (auto it = required.begin(); it != required.end(); ++it) {
            const auto& item = required.sublist(required.name(it));
            files[item.get<std::string>("File")].push_back(item);
        }
        for (const auto& file : files) {
            auto map = makeMap(std::stoi(file.second[0].get<std::string>("Global DOFs")));
            auto values = Teuchos::rcp(new MultiVector<>(map));
            values->putScalar(1.);
            HDF5Export<> output(map, directory + "/" + file.first.substr(0, file.first.size() - 3));
            for (const auto& item : file.second)
                output.writeVariablesHDF5(item.get<std::string>("Key"), values);
            output.closeExporter();
        }
        checkpoint::write(selected, parameters, time, *comm);
    };
    auto newmark = schema;
    newmark.sublist("Fields").remove("p");
    auto& newmarkIntegration = newmark.sublist("Integration");
    newmarkIntegration.remove("BDF");
    newmarkIntegration.remove("Extrapolation");
    newmarkIntegration.set("Class", std::string("Newmark"))
        .set("Layout", std::string("legacy-start-of-step-v1"))
        .set("beta", 0.25).set("gamma", 0.5);
    auto newmarkSettings = settings;
    newmarkSettings.sublist("Timestepping Parameter").set("Class", std::string("Newmark"));
    checksSchema("version-1 Newmark schema", newmark, newmarkSettings, {fields[0]});
    auto structure = newmark;
    structure.set("Role", std::string("FSI structure"));
    checksSchema("version-1 FSI structure overrides integration class", structure, settings, {fields[0]});
    fixtureFor(newmark);
    rejects("initial solution excludes Newmark", "standalone multistep Navier-Stokes", [&] { validateInitial(newmark); });
    validate(newmark);
    for (const std::string setting : {std::string("beta"), std::string("gamma")}) {
        auto changed = newmark;
        changed.sublist("Integration").set(setting, 0.1);
        rejects("changed Newmark " + setting, setting, [&] { validate(changed); });
    }
    checkpoint::onRoot(*comm, [&] {
        std::filesystem::remove(directory + "/ds_Acceleration.h5");
    });
    rejects("missing Newmark acceleration", "ds_Acceleration", [&] { validate(newmark); });
    auto fsi = schema;
    fsi.set("FSI", true);
    fsi.sublist("Integration").set("Layout", std::string("legacy-start-of-step-v1"))
        .set("beta", 0.25).set("gamma", 0.5).set("Geometry Explicit", true);
    const auto field = fsi.sublist("Fields").sublist("u");
    fsi.sublist("Fields").remove("u");
    for (const std::string name : {std::string("u_f"), std::string("d_s"), std::string("lambda"), std::string("d_f")})
        fsi.sublist("Fields").sublist(name).setParameters(field);
    auto fsiSettings = settings;
    fsiSettings.sublist("Parameter").set("FSI", true);
    std::vector<checkpoint::FieldDescription> fsiFields = {fields[1]};
    for (const std::string name : {std::string("u_f"), std::string("d_s"), std::string("lambda"), std::string("d_f")}) {
        auto description = fields[0];
        description.name = name;
        fsiFields.push_back(description);
    }
    checksSchema("version-1 coupled FSI schema", fsi, fsiSettings, fsiFields);
    fixtureFor(fsi);
    rejects("initial solution excludes FSI", "standalone multistep Navier-Stokes", [&] { validateInitial(fsi); });
    validate(fsi);
    checkpoint::onRoot(*comm, [&] {
        std::filesystem::remove(directory + "/Solutiond_f.h5");
    });
    rejects("missing FSI ALE displacement", "Solutiond_f", [&] { validate(fsi); });
    fixtureFor(fsi);
    checkpoint::onRoot(*comm, [&] {
        checkpoint::H5Handle file(H5Fopen((directory + "/Rhsu.h5").c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        H5Ldelete(file, std::to_string(time - 0.0025).c_str(), H5P_DEFAULT);
    });
    rejects("missing FSI moving-mesh mass history", "Rhsu", [&] { validate(fsi); });
    fixtureFor(fsi);
    auto changedFsi = fsi;
    changedFsi.sublist("Integration").set("Geometry Explicit", false);
    rejects("changed FSI geometry integration", "Geometry Explicit", [&] { validate(changedFsi); });

    auto outletSettings = fsiSettings;
    outletSettings.sublist("Parameter Fluid")
        .set("Pressure Boundary Condition", std::string("Absorbing Paper"))
        .set("Average Flowrate", true).set("Unsteady Start", 0.006);
    const auto outletFsi = checkpoint::makeSchema(outletSettings, fsiFields);
    fixtureFor(outletFsi);
    FSIOutletState outlet;
    outlet.initialInletArea = 0.41234567890123456;
    outlet.initialOutletArea = 0.39876543210987654;
    outlet.transitionOutletArea = 0.40000000000000013;
    outlet.currentFlowRate = 1.2345678901234567e-6;
    outlet.previousFlowRate = 9.876543210987654e-7;
    outlet.pressure = 123.45678901234567;
    outlet.initialized = true;
    outlet.transitionCaptured = true;
    const std::string outletFile = directory + "/" + checkpoint::outletStateName(time);
    const auto writeOutlet = [&] {
        checkpoint::onRoot(*comm, [&] {
            checkpoint::writeOutletState(outlet, outletFile, "Absorbing Paper", time);
        });
    };
    writeOutlet();
    validate(outletFsi);
    checkpoint::onRoot(*comm, [&] {
        const auto read = checkpoint::readOutletState(outletFile, "Absorbing Paper", time);
        if (read.initialInletArea != outlet.initialInletArea || read.initialOutletArea != outlet.initialOutletArea ||
            read.transitionOutletArea != outlet.transitionOutletArea || read.currentFlowRate != outlet.currentFlowRate ||
            read.previousFlowRate != outlet.previousFlowRate || read.pressure != outlet.pressure ||
            read.initialized != outlet.initialized || read.transitionCaptured != outlet.transitionCaptured)
            throw std::runtime_error("FSI outlet scalar round-trip lost precision");
    });
    if (comm->getRank() == 0) std::cout << "PASS full-precision FSI outlet state round-trip\n";
    const auto editOutlet = [&](const std::function<void(Teuchos::ParameterList&)>& edit) {
        checkpoint::onRoot(*comm, [&] {
            auto data = Teuchos::getParametersFromXmlFile(outletFile);
            edit(*data);
            std::ofstream output(outletFile);
            Teuchos::writeParameterListToXmlOStream(*data, output);
        });
    };
    editOutlet([](Teuchos::ParameterList& data) { data.set("Initial outlet area", 0.); });
    rejects("invalid absorbing reference area", "positive initial areas", [&] { validate(outletFsi); });
    writeOutlet();
    editOutlet([](Teuchos::ParameterList& data) { data.set("Initialized", false); });
    rejects("uninitialized resumed outlet history", "not initialized", [&] { validate(outletFsi); });
    writeOutlet();
    editOutlet([](Teuchos::ParameterList& data) { data.set("Transition captured", false); });
    rejects("missing captured outlet transition", "captured transition area", [&] { validate(outletFsi); });
    writeOutlet();
    editOutlet([](Teuchos::ParameterList& data) { data.remove("Previous flow rate"); });
    rejects("missing previous outlet flow rate", "Previous flow rate", [&] { validate(outletFsi); });
    writeOutlet();
    auto changedOutlet = outletFsi;
    changedOutlet.sublist("FSI outlet").set("Unsteady Start", 0.005);
    rejects("changed outlet transition time", "Unsteady Start", [&] { validate(changedOutlet); });
    auto nonfiniteOutlet = outlet;
    nonfiniteOutlet.pressure = std::numeric_limits<double>::infinity();
    rejects("nonfinite outlet state", "nonfinite", [&] {
        checkpoint::validateOutletState(nonfiniteOutlet, "Absorbing Paper", time);
    });

    // Density and kinematic viscosity both affect the resistance traction.
    // Rebuilding the schema from changed input must reject the old checkpoint.
    auto resistanceSettings = fsiSettings;
    resistanceSettings.sublist("Parameter Fluid")
        .set("Pressure Boundary Condition", std::string("Resistance"))
        .set("Density", 0.5).set("Viscosity", 0.25);
    const auto resistanceFsi = checkpoint::makeSchema(resistanceSettings, fsiFields);
    fixtureFor(resistanceFsi);
    checkpoint::onRoot(*comm, [&] {
        checkpoint::writeOutletState(outlet, outletFile, "Resistance", time);
    });
    validate(resistanceFsi);
    for (const std::string setting : {std::string("Density"), std::string("Viscosity")}) {
        auto changedSettings = resistanceSettings;
        changedSettings.sublist("Parameter Fluid").set(setting, 0.125);
        const auto changedResistance = checkpoint::makeSchema(changedSettings, fsiFields);
        rejects("changed resistance " + setting, setting, [&] { validate(changedResistance); });
    }
    // Real version-2 fixtures retain actual history independently of the next dt.
    auto variableSchema = checkpoint::makeSchema(settings, fields);
    auto& runtime = parameters->sublist("Timestepping Parameter").sublist("_Checkpoint time state");
    runtime.set("Physical time", time).set("dt", 0.003).set("dt_prev", 0.0025).set("Step number", 4LL);
    runtime.sublist("History times").set("0", time).set("1", 0.0075).set("2", 0.006);
    fixtureFor(variableSchema);
    auto changedNextDt = variableSchema;
    changedNextDt.sublist("Integration").set("dt", 0.0017);
    validate(changedNextDt);
    validateInitial(changedNextDt);
    if (comm->getRank() == 0) std::cout << "PASS version-2 history with changed next dt and off-grid time\n";
    const auto metadataFile = directory + "/" + checkpoint::manifestName(variableSchema, time);
    const auto editClock = [&](const std::function<void(Teuchos::ParameterList&)>& edit) {
        checkpoint::onRoot(*comm, [&] {
            auto description = Teuchos::getParametersFromXmlFile(metadataFile);
            edit(description->sublist("Time state"));
            std::ofstream output(metadataFile);
            Teuchos::writeParameterListToXmlOStream(*description, output);
        });
    };
    editClock([](Teuchos::ParameterList& state) { state.set("dt_prev", 0.); });
    rejects("invalid saved incoming increment", "dt_prev", [&] { validate(changedNextDt); });
    fixtureFor(variableSchema);
    editClock([](Teuchos::ParameterList& state) { state.sublist("History times").set("1", 0.009); });
    rejects("inconsistent saved history time", "history times", [&] { validate(changedNextDt); });
    fixtureFor(variableSchema);
    editClock([](Teuchos::ParameterList& state) { state.set("Step number", -1LL); });
    rejects("negative saved step number", "Step number", [&] { validate(changedNextDt); });
    fixtureFor(variableSchema);
    editVelocity([&](hid_t file) { H5Ldelete(file, "0.007500", H5P_DEFAULT); });
    rejects("missing version-2 actual history", "missing required history/field group", [&] { validate(changedNextDt); });
    checkpoint::onRoot(*comm, [&] { std::filesystem::remove_all(directory); });
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
