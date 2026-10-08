#include "feddlib/core/General/CheckpointMetadata.hpp"
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
        .set("Checkpoint directory", directory).set("Restart directory", directory);
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
    const auto validate = [&](const Teuchos::ParameterList& required) {
        checkpoint::validate(required, parameters, time, *comm);
    };
    fixture();
    validate(schema);
    if (comm->getRank() == 0) std::cout << "PASS compatible checkpoint\n";
    for (const std::string name : {std::string("Reference mesh and DOF fingerprint"), std::string("FE type"), std::string("Global DOFs")}) {
        auto changed = schema;
        changed.sublist("Fields").sublist("u").set(name, std::string("incompatible"));
        rejects("incompatible " + name, name, [&] { validate(changed); });
    }
    {
        auto changed = schema;
        changed.sublist("Integration").set("BDF", 1);
        rejects("changed BDF order", "BDF", [&] { validate(changed); });
        changed = schema;
        changed.sublist("Integration").set("dt", 0.005);
        rejects("changed timestep", "dt", [&] { validate(changed); });
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
    fixture();
    for (const std::string property : {std::string("GlobalLength"), std::string("NumVectors")}) {
        editVelocity([&](hid_t file) {
            checkpoint::H5Handle dataset(H5Dopen(file, (std::to_string(time) + "/" + property).c_str(), H5P_DEFAULT), H5Dclose);
            int wrong = 99;
            H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &wrong);
        });
        rejects("inconsistent " + property, property, [&] { validate(schema); });
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
    parameters->sublist("Timestepping Parameter").set("Allow legacy restart", true);
    validate(schema);
    if (comm->getRank() == 0) std::cout << "PASS explicit legacy restart\n";
    // Legacy mode must not bypass dataset checks.
    editVelocity([&](hid_t file) { H5Ldelete(file, std::to_string(time - 0.0025).c_str(), H5P_DEFAULT); });
    rejects("legacy restart missing history", "missing required history/field group", [&] { validate(schema); });
    // Newmark and coupled FSI add distinct auxiliary data to the same validator.
    parameters->sublist("Timestepping Parameter").set("Allow legacy restart", false);
    const auto fixtureFor = [&](const Teuchos::ParameterList& selected) {
        const auto description = checkpoint::atTime(selected, time);
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
    fixtureFor(newmark);
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
    fixtureFor(fsi);
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
    checkpoint::onRoot(*comm, [&] { std::filesystem::remove_all(directory); });
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
