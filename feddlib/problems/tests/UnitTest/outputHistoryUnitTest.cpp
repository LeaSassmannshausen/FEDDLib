#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/General/ExporterTxt.hpp"
#include "feddlib/core/General/OutputHistory.hpp"
#include "feddlib/core/General/HDF5VectorInfo.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <iostream>

using namespace FEDD;

int main(int argc, char** argv)
{
    Teuchos::GlobalMPISession session(&argc, &argv);
    auto comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    auto meshParameters = Teuchos::rcp(new Teuchos::ParameterList);
    meshParameters->set("Mesh 1 Name", std::string("outputHistory-square.mesh"));
    auto domain = Teuchos::rcp(new Domain<>(comm, 2));
    MeshPartitioner<>::DomainPtrArray_Type domains(1, domain);
    MeshPartitioner<> partitioner(domains, meshParameters, "P1", 2);
    partitioner.readAndPartition();
    const auto originalDirectory = std::filesystem::current_path();
    const std::string directory = "outputHistoryUnitTest-data";
    output::onRank(*comm, 0, [&] {
        std::filesystem::remove_all(directory);
        std::filesystem::create_directory(directory);
    });
    std::filesystem::current_path(directory);
    auto settings = [] {
        auto p = Teuchos::rcp(new Teuchos::ParameterList);
        p->sublist("Timestepping Parameter").set("dt", 0.1).set("Restart", false).set("Time step", 0.1);
        p->sublist("Exporter").set("Write new mesh", true);
        return p;
    };
    auto points = domain->getMesh()->getPointsUnique();
    const auto originalPoints = *points;
    auto pressure = Teuchos::rcp(new MultiVector<>(domain->getMapUnique()));
    auto velocity = Teuchos::rcp(new MultiVector<>(domain->getMapVecFieldUnique()));
    const auto paraview = [&](const std::string& name, ParameterListPtr_Type parameters,
                              const std::vector<double>& times, double offset, bool withDt = true) {
        ExporterParaView<> exporter;
        exporter.setup(name, domain->getMesh(), "P1", parameters);
        Teuchos::RCP<const MultiVector<>> p = pressure, u = velocity;
        exporter.addVariable(p, "pressure", "Scalar", 1, domain->getMapUnique());
        exporter.addVariable(u, "velocity", "Vector", 2, domain->getMapUnique());
        for (double time : times) {
            pressure->putScalar(time + offset);
            velocity->putScalar(time + offset);
            for (std::size_t i = 0; i < points->size(); ++i)
                (*points)[i][0] = originalPoints[i][0] + time * (offset == 0. ? 1. : 2.);
            exporter.updatePoints();
            if (withDt) exporter.save(time, 0.1);
            else exporter.save(time);
        }
        exporter.closeExporter();
    };
    const auto text = [&](ParameterListPtr_Type parameters, const std::vector<double>& times, double offset) {
        ExporterTxt outlet, iterations;
        outlet.setup("pressureOutlet", comm, 0, parameters, true);
        iterations.setup("iterations", comm, 1, parameters);
        for (double time : times) {
            outlet.exportData(time, time + offset);
            iterations.exportDataAtTime(time, time + offset);
        }
        outlet.closeExporter();
        iterations.closeExporter();
    };
    const auto require = [](bool condition, const std::string& message) {
        if (!condition) throw std::runtime_error(message);
    };
    const auto firstValue = [&](hid_t file, const std::string& group) {
        checkpoint::H5Handle dataset(H5Dopen(file, (group + "/Values").c_str(), H5P_DEFAULT), H5Dclose);
        require(dataset >= 0, "Missing " + group);
        checkpoint::H5Handle space(H5Dget_space(dataset), H5Sclose);
        std::vector<double> values(H5Sget_simple_extent_npoints(space));
        require(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) >= 0, "Cannot read values");
        return values.front();
    };
    const auto checkSeries = [&](const std::string& name, std::size_t count) {
        output::Collection series(output::read(name + ".xmf"));
        require(series.frames.size() == count, "Wrong XMF frame count");
        checkpoint::H5Handle file(H5Fopen((name + ".h5").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
        require(file >= 0, "Cannot open output HDF5");
        for (const auto& frame : series.frames)
            for (const auto& group : frame.groups)
                require(H5Lexists(file, (group + "/Values").c_str(), H5P_DEFAULT) > 0, "Dangling XMF reference");
    };

    auto original = settings();
    paraview("fields", original, {0., 0.1, 0.2, 0.3}, 0.);
    text(original, {0., 0.1, 0.2, 0.3}, 0.);
    auto resumed = settings();
    resumed->sublist("Timestepping Parameter").set("Restart", true);
    resumed->sublist("Exporter").set("Resume output", true).set("Keep old output", true);
    paraview("fields", resumed, {0.1, 0.2, 0.3, 0.4}, 100.);
    text(resumed, {0.1, 0.2, 0.3, 0.4}, 100.);
    output::onRank(*comm, 0, [&] {
        checkSeries("fields", 5);
        require(output::Collection(output::read("fields_times.xmf")).frames.size() == 5, "Wrong dt collection");
        require(!std::filesystem::exists("old_output"), "Resume must take precedence over archive");
        checkpoint::H5Handle file(H5Fopen("fields.h5", H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
        require(std::abs(firstValue(file, "pressure.00001") - 0.1) < 1.e-14, "Earlier pressure overwritten");
        require(std::abs(firstValue(file, "pressure.00002") - 100.2) < 1.e-14, "Future pressure was not replaced");
        require(std::abs(firstValue(file, "velocity.00002") - 100.2) < 1.e-14, "Vector output differs");
        require(std::abs(firstValue(file, "PointsX2") - firstValue(file, "PointsX1") - 0.3) < 1.e-13, "Moving mesh was not resumed");
        std::istringstream values(output::read("pressureOutlet.txt"));
        double t, v; unsigned count = 0;
        while (values >> t >> v) {
            require(std::abs(v - (count < 2 ? t : t + 100.)) < 1.e-12, "Text values differ");
            ++count;
        }
        require(count == 5, "Duplicate/missing text rows");
        std::istringstream times(output::read("iterations.txt.times"));
        count = 0; while (times >> t) ++count;
        require(count == 5, "Value-only timestamps differ");
        std::cout << "PASS rewind, duplicate suppression, moving mesh, scalar/vector and text history\n";
    });
    // Appending at the final frame must preserve all old frames.
    resumed->sublist("Timestepping Parameter").set("Time step", 0.4);
    paraview("fields", resumed, {0.4, 0.5}, 100., false);
    text(resumed, {0.4, 0.5}, 100.);
    output::onRank(*comm, 0, [&] {
        checkSeries("fields", 6);
        require(output::Collection(output::read("fields_times.xmf")).frames.size() == 6, "Existing dt collection not extended");
    });

    // One archive contains the old HDF5, both XMF collections and all text files.
    auto archived = settings();
    archived->sublist("Timestepping Parameter").set("Restart", true);
    archived->sublist("Exporter").set("Keep old output", true);
    paraview("fields", archived, {0.1}, 200.);
    text(archived, {0.1}, 200.);
    output::onRank(*comm, 0, [&] {
        const auto folder = archived->sublist("Exporter").get<std::string>("_Old output archive");
        for (const auto& file : {"fields.h5", "fields.xmf", "fields_times.xmf", "pressureOutlet.txt", "iterations.txt", "iterations.txt.times"})
            require(std::filesystem::exists(std::filesystem::path(folder) / file), "Missing archived output");
        checkSeries("fields", 1);
        require(output::Collection(output::read(folder + "/fields.xmf")).frames.size() == 6, "Archived history changed");
        std::cout << "PASS shared output archive\n";
    });
    // Incompatible existing series must fail collectively without changing it.
    auto incompatible = settings();
    incompatible->sublist("Timestepping Parameter").set("Restart", true);
    incompatible->sublist("Exporter").set("Resume output", true);
    std::string before;
    output::onRank(*comm, 0, [&] { before = output::read("fields.xmf"); });
    int rejected = 0;
    try {
        ExporterParaView<> wrong;
        wrong.setup("fields", domain->getMesh(), "P1", incompatible);
        Teuchos::RCP<const MultiVector<>> p = pressure;
        wrong.addVariable(p, "different-pressure", "Scalar", 1, domain->getMapUnique());
        wrong.save(0.1);
    } catch (const std::exception& e) {
        rejected = std::string(e.what()).find("fields do not match") != std::string::npos;
    }
    int allRejected = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &rejected, &allRejected);
    require(allRejected, "Incompatible fields did not fail on every rank");
    output::onRank(*comm, 0, [&] {
        require(output::read("fields.xmf") == before, "Rejected resume changed the series");
        checkSeries("fields", 1);
        std::cout << "PASS collective compatibility rejection preserves output\n";
    });
    auto secondArchive = settings();
    secondArchive->sublist("Exporter").set("Keep old output", true);
    paraview("fields", secondArchive, {0.2}, 300.);
    output::onRank(*comm, 0, [&] {
        const auto firstFolder = archived->sublist("Exporter").get<std::string>("_Old output archive");
        const auto secondFolder = secondArchive->sublist("Exporter").get<std::string>("_Old output archive");
        require(firstFolder != secondFolder, "Prior archive reused");
        require(output::Collection(output::read(firstFolder + "/fields.xmf")).frames.size() == 6, "Prior archive overwritten");
        std::cout << "PASS repeated archiving preserves earlier archives\n";
    });
    paraview("defaults", settings(), {0., 0.1}, 0.);
    paraview("defaults", settings(), {0.2}, 0.);
    output::onRank(*comm, 0, [&] { checkSeries("defaults", 1); });

    // Every-other-step output must keep the original export cadence after resume.
    auto sparse = settings();
    sparse->sublist("Exporter").set("Write new mesh", false).set("Export every X timesteps", 2);
    paraview("sparse", sparse, {0., 0.1, 0.2, 0.3, 0.4}, 0.);
    sparse->sublist("Timestepping Parameter").set("Restart", true).set("Time step", 0.3);
    sparse->sublist("Exporter").set("Resume output", true);
    paraview("sparse", sparse, {0.3, 0.4, 0.5, 0.6}, 100.);
    output::onRank(*comm, 0, [&] { checkSeries("sparse", 4); });
    sparse->sublist("Exporter").set("Export every X timesteps", 1);
    rejected = 0;
    try { paraview("sparse", sparse, {0.3}, 0.); }
    catch (const std::exception& e) {
        rejected = std::string(e.what()).find("Export interval differs") != std::string::npos;
    }
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &rejected, &allRejected);
    require(allRejected, "Changed export cadence did not fail collectively");

    // Legacy pressure logs already contain time; value-only logs use time.txt.
    output::onRank(*comm, 0, [&] {
        output::replace("legacy.txt", "0.1 1\n0.2 2\n0.3 3\n");
        output::replace("time.txt", "0\n0.1\n0.2\n0.3\n");
        output::replace("legacyIts.txt", "1\n2\n3\n");
    });
    auto legacy = settings();
    legacy->sublist("Timestepping Parameter").set("Restart", true);
    legacy->sublist("Exporter").set("Resume output", true);
    ExporterTxt legacyOutlet, legacyIts;
    legacyOutlet.setup("legacy", comm, 0, legacy, true);
    legacyIts.setup("legacyIts", comm, 1, legacy);
    legacyOutlet.exportData(0.2, 20.);
    legacyIts.exportDataAtTime(0.2, 20.);
    legacyOutlet.closeExporter(); legacyIts.closeExporter();
    output::onRank(*comm, 0, [&] {
        std::istringstream rows(output::read("legacy.txt"));
        double t, v;
        require(bool(rows >> t >> v) && std::abs(t - 0.1) < 1.e-14 && v == 1., "Legacy first row differs");
        require(bool(rows >> t >> v) && std::abs(t - 0.2) < 1.e-14 && v == 20., "Legacy resumed row differs");
        require(!(rows >> t), "Legacy future rows retained");
        require(output::read("legacyIts.txt") == "1\n20\n", "Legacy value-only rows differ");
        std::cout << "PASS legacy text series and sparse export cadence\n";
    });
    std::filesystem::current_path(originalDirectory);
    if (comm->getRank() == 0) std::cout << "PASS output history\n";
}
