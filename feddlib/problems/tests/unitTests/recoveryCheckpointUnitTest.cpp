#include "feddlib/core/Checkpointing/RecoveryCheckpoint.hpp"
#include "feddlib/problems/Solver/NonlinearSolveReport.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <filesystem>

using namespace FEDD;

// Boundary criteria and real distributed checkpoint publication without a PDE solve.
int main(int argc, char** argv)
{
    Teuchos::GlobalMPISession session(&argc, &argv);
    auto comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    int failures = 0;
    const auto check = [&](bool passed, const char* label) {
        int local = passed ? 1 : 0, global = 0;
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &local, &global);
        failures += !global;
        if (comm->getRank() == 0) std::cout << (global ? "PASS " : "FAIL ") << label << '\n';
    };
    NonlinearSolveReport report;
    report.converged = true;
    report.maximumLinearIterations = 100;
    checkpoint::RecoveryCheckpointPolicy policy;
    check(!policy.triggered(report) && report.averageLinearIterations() == 0., "no linear solves");
    report.linearSolveCount = 2;
    report.totalLinearIterations = 100;
    check(!policy.triggered(report), "average exactly maxIter/2 does not trigger");
    report.totalLinearIterations = 101;
    check(policy.triggered(report), "average above maxIter/2 triggers");
    report.totalLinearIterations = 0;
    report.nonlinearIterations = 10;
    check(!policy.triggered(report), "convergence on final Newton iteration is allowed");
    report.converged = false;
    check(policy.triggered(report), "nonconvergence triggers independently of linear count");
    policy.observe(report);
    report.converged = true;
    policy.observe(report);
    check(!policy.trusted(), "later convergence cannot replace the pre-failure state");
    report.linearSolvesConverged = false;
    check(policy.triggered(report), "one failed linear solve triggers");
    report.linearSolvesConverged = true;
    report.finite = false;
    check(policy.triggered(report), "nonfinite solution triggers");

    Teuchos::ParameterList solverParameters;
    auto& thyra = solverParameters.sublist("ThyraSolver");
    thyra.set("Linear Solver Type", std::string("Belos"));
    auto& belos = thyra.sublist("Linear Solver Types").sublist("Belos");
    belos.set("Solver Type", std::string("Pseudo Block GMRES"));
    belos.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Maximum Iterations", 77);
    belos.sublist("Solver Types").sublist("Block GMRES").set("Maximum Iterations", 222);
    check(configuredLinearIterationLimit(solverParameters) == 77, "iteration limit uses selected solver");
    Thyra::BelosLinearOpWithSolveFactory<default_sc> factory;
    factory.setParameterList(Teuchos::rcpFromRef(belos));
    check(configuredLinearIterationLimit(solverParameters) == 77,
          "iteration limit survives Trilinos string-to-enum conversion");

    auto parameters = Teuchos::rcp(new Teuchos::ParameterList);
    const std::string directory = "recoveryCheckpointUnitTest-data";
    parameters->sublist("Timestepping Parameter").set("Failure recovery", true)
        .set("Recovery directory", directory).set("Class", std::string("Multistep"))
        .set("BDF", 2).set("dt", 0.0025);
    checkpoint::onRoot(*comm, [&] { std::filesystem::remove_all(directory); });
    const auto schema = checkpoint::makeSchema(*parameters, {{"u", "P1", 2, 2, "8", "0", "test-mesh"}});
    Teuchos::RCP<const Map<>> map = Teuchos::rcp(new Map<>(8, 8 / comm->getSize(), 0, comm));
    auto values = Teuchos::rcp(new MultiVector<>(map));
    values->putScalar(1.);
    // A duplicate dataset exercises the actual HDF5 failure path collectively.
    checkpoint::onRoot(*comm, [&] { checkpoint::makeRecoveryDirectory(directory); });
    {
        HDF5Export<> duplicate(map, joinPath(directory, "duplicate"));
        duplicate.writeVariablesHDF5("state", values);
        bool failedTogether = false;
        try { duplicate.writeVariablesHDF5("state", values); }
        catch (const std::runtime_error& error) {
            failedTogether = std::string(error.what()).find("HDF5 vector write failed") != std::string::npos;
        }
        check(failedTogether, "HDF5 vector write errors reach every MPI rank");
        duplicate.closeExporter();
    }
    using Snapshot = checkpoint::RecoveryCheckpointSnapshot<default_sc,default_lo,default_go,default_no>;
    const auto capture = [&](double time) {
        auto snapshot = Teuchos::rcp(new Snapshot(time));
        snapshot->addVector("Solutionu", time, values);
        snapshot->addVector("Solutionu", time - 0.0025, values);
        snapshot->addManifest(schema);
        return snapshot;
    };
    checkpoint::RecoveryCheckpointManager<default_sc,default_lo,default_go,default_no> manager(parameters, *comm);
    auto snapshot = capture(0.01);
    values->putScalar(9.); // Captured state must remain independent of live vectors.
    manager.beforeSolve(snapshot);
    report.finite = true;
    report.linearSolveCount = 2;
    report.totalLinearIterations = 101;
    manager.afterSolve(report, 0.0125);
    const auto latest = Teuchos::getParametersFromXmlFile(joinPath(directory, "Latest.xml"));
    const auto latestDirectory = latest->get<std::string>("Directory");
    HDF5Import<> input(map, joinPath(latestDirectory, "Solutionu"));
    const auto saved = input.readVariablesHDF5("0.010000");
    const auto data = saved->getData(0);
    bool independent = true;
    for (int i = 0; i < data.size(); ++i) independent &= data[i] == 1.;
    check(independent, "checkpoint writes immutable captured values");

    auto interrupted = capture(0.0125);
    interrupted->addScalarWriter("injected.xml", [](const std::string&) {
        throw std::runtime_error("injected checkpoint write failure");
    });
    bool rejected = false;
    try { manager.beforeSolve(interrupted); }
    catch (const std::runtime_error& error) {
        rejected = std::string(error.what()).find("injected checkpoint write failure") != std::string::npos;
    }
    check(rejected, "partial checkpoint failure reaches every MPI rank");
    const auto retained = Teuchos::getParametersFromXmlFile(joinPath(directory, "Latest.xml"));
    check(retained->get<std::string>("Directory") == latestDirectory &&
          std::ifstream(joinPath(latestDirectory, "Complete.xml")).good(),
          "partial write preserves previous complete generation");
    check(!std::ifstream(joinPath(directory, "generation_2/Complete.xml")).good(),
          "partial generation is never published");

    auto continuationParameters = Teuchos::rcp(new Teuchos::ParameterList(*parameters));
    const auto continuationDirectory = joinPath(directory, "continued");
    continuationParameters->sublist("Timestepping Parameter").set("Recovery directory", continuationDirectory);
    continuationParameters->sublist("Parameter").set("Cancel MaxNonLinIts", false);
    checkpoint::RecoveryCheckpointManager<default_sc,default_lo,default_go,default_no>
        continuing(continuationParameters, *comm);
    continuing.beforeSolve(capture(0.01));
    report.totalLinearIterations = 2;
    report.nonlinearIterations = 1;
    report.converged = false;
    check(!continuing.afterSolve(report, 0.0125), "ignored failure saves recovery without stopping");
    report.converged = true;
    check(!continuing.afterSolve(report, 0.015), "later converged step can continue");
    const auto status = Teuchos::getParametersFromXmlFile(joinPath(continuationDirectory, "RecoveryStatus.xml"));
    check(!status->get<bool>("Trajectory reliable") && status->get<double>("Recovery time") == 0.01 &&
          status->get<std::string>("Trigger reason") == "Unreliable continuation after earlier failure",
          "continued trajectory retains the earlier failure and reliable restart time");
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
