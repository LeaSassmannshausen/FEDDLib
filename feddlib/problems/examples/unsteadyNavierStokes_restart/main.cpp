#include <Tpetra_Core.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_StackedTimer.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace FEDD;
using Teuchos::RCP;
using Teuchos::rcp;

typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;
typedef Domain<SC, LO, GO, NO> Domain_Type;
typedef MeshPartitioner<SC, LO, GO, NO> MeshPartitioner_Type;
typedef NavierStokes<SC, LO, GO, NO> NavierStokes_Type;

namespace {

/// Impose no slip on the three velocity components.
void zeroDirichlet3D(double*, double* result, double, const double*)
{
    result[0] = result[1] = result[2] = 0.;
}

/** @brief Biquadratic inlet profile for the supplied 3D BFS mesh.
 * @param[in] x Boundary coordinates; the inlet spans y,z in [0, height].
 * @param[out] result Velocity, directed along the positive x axis.
 * @param[in] parameters Peak velocity followed by inlet height.
 */
void inflowParabolic3D(double* x, double* result, double, const double* parameters)
{
    const double height = parameters[1];
    result[0] = 16. * parameters[0] * x[1] * (height - x[1]) * x[2] * (height - x[2])
                / (height * height * height * height);
    result[1] = result[2] = 0.;
}

} // namespace

/** @brief Run or restart unsteady Navier-Stokes on a 3D unstructured mesh.
 *
 * Derived from the unsteadyNavierStokes example. The optional restart XML
 * overrides the base problem parameters. Problem initialization restores primary
 * fields; DAESolverInTime restores the integration clock and BDF history.
 * Checkpoints are written after accepted timesteps, including the final step.
 */
int main(int argc, char* argv[])
{
    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

    std::string problemFile = "parametersProblem.xml";
    std::string restartFile;
    std::string precFile = "parametersPrec.xml";
    std::string solverFile = "parametersSolver.xml";
    Teuchos::CommandLineProcessor commandLine;
    commandLine.setOption("problemfile", &problemFile, "Base problem XML.");
    commandLine.setOption("restartfile", &restartFile, "Optional restart XML overriding the base problem.");
    commandLine.setOption("precfile", &precFile, "Monolithic preconditioner XML.");
    commandLine.setOption("solverfile", &solverFile, "Linear solver XML.");
    commandLine.recogniseAllOptions(true);
    commandLine.throwExceptions(false);
    const Teuchos::CommandLineProcessor::EParseCommandLineReturn parseResult = commandLine.parse(argc, argv);
    if (parseResult == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
        return EXIT_SUCCESS;
    if (parseResult != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
        return EXIT_FAILURE;

    ParameterListPtr_Type parameters = Teuchos::getParametersFromXmlFile(problemFile);
    if (!restartFile.empty())
        parameters->setParameters(*Teuchos::getParametersFromXmlFile(restartFile));
    parameters->setParameters(*Teuchos::getParametersFromXmlFile(precFile));
    parameters->setParameters(*Teuchos::getParametersFromXmlFile(solverFile));

    Teuchos::ParameterList& meshParameters = parameters->sublist("Parameter");
    const int dim = meshParameters.get("Dimension", 3);
    const std::string meshType = meshParameters.get("Mesh Type", "unstructured");
    const std::string feTypeVelocity = meshParameters.get("Discretization Velocity", "P2");
    const std::string feTypePressure = meshParameters.get("Discretization Pressure", "P1");
    TEUCHOS_TEST_FOR_EXCEPTION(dim != 3 || meshType != "unstructured", std::logic_error,
                               "unsteadyNavierStokes_restart requires Dimension=3 and Mesh Type=unstructured.");
    TEUCHOS_TEST_FOR_EXCEPTION(feTypePressure != "P1" || (feTypeVelocity != "P1" && feTypeVelocity != "P2"),
                               std::logic_error, "Choose P1/P1 or P2/P1 velocity/pressure discretization.");
    TEUCHOS_TEST_FOR_EXCEPTION(parameters->sublist("Timestepping Parameter").get("Class", "Multistep") != "Multistep",
                               std::logic_error, "This restart example uses Multistep (BDF) time integration.");

    RCP<Teuchos::StackedTimer> stackedTimer = rcp(new Teuchos::StackedTimer("Unsteady Navier-Stokes restart", true));
    Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
    {
        const bool restart = parameters->sublist("Timestepping Parameter").get("Restart", false);
        if (comm->getRank() == 0)
            std::cout << "3D unstructured Navier-Stokes: " << (restart ? "restart" : "fresh run")
                      << " on " << comm->getSize() << " MPI ranks." << std::endl;

        RCP<Domain_Type> pressureDomain = rcp(new Domain_Type(comm, dim));
        RCP<Domain_Type> velocityDomain;
        {
            Teuchos::TimeMonitor meshTimer(*Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
            MeshPartitioner_Type::DomainPtrArray_Type domains(1);
            domains[0] = pressureDomain;
            ParameterListPtr_Type partitionerParameters = Teuchos::sublist(parameters, "Mesh Partitioner");
            MeshPartitioner_Type partitioner(domains, partitionerParameters, "P1", dim);
            partitioner.readAndPartition();
            if (feTypeVelocity == "P2") {
                velocityDomain = rcp(new Domain_Type(comm, dim));
                velocityDomain->buildP2ofP1Domain(pressureDomain);
            }
            else
                velocityDomain = pressureDomain;
        }

        // BFS3dCC.mesh: 1 = walls, 2 = inlet, 3 = natural outflow.
        // Marker 4 is an optional obstacle boundary, as in the original example.
        RCP<BCBuilder<SC, LO, GO, NO>> boundaries = rcp(new BCBuilder<SC, LO, GO, NO>());
        std::vector<double> inflowParameters = {
            meshParameters.get("MaxVelocity", 1.5), meshParameters.get("Inlet height", 1.)};
        boundaries->addBC(zeroDirichlet3D, 1, 0, velocityDomain, "Dirichlet", dim);
        boundaries->addBC(zeroDirichlet3D, 4, 0, velocityDomain, "Dirichlet", dim);
        boundaries->addBC(inflowParabolic3D, 2, 0, velocityDomain, "Dirichlet", dim, inflowParameters);

        NavierStokes_Type problem(velocityDomain, feTypeVelocity, pressureDomain, feTypePressure, parameters);
        problem.addBoundaries(boundaries);
        problem.initializeProblem();
        problem.assemble();
        problem.setBoundariesRHS();

        SmallMatrix<int> timeBlocks(2);
        timeBlocks[0][0] = 1;
        timeBlocks[0][1] = 1;
        timeBlocks[1][0] = 0;
        timeBlocks[1][1] = 0;
        DAESolverInTime<SC, LO, GO, NO> timeSolver(parameters, comm);
        timeSolver.defineTimeStepping(timeBlocks);
        timeSolver.setProblem(problem);
        timeSolver.setupTimeStepping();
        {
            Teuchos::TimeMonitor solveTimer(*Teuchos::TimeMonitor::getNewCounter("main: Solve problem"));
            timeSolver.advanceInTime();
        }
        problem.infoParameter();
    }

    Teuchos::TimeMonitor::report(std::cout);
    stackedTimer->stop("Unsteady Navier-Stokes restart");
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report(std::cout, comm, options);
    return EXIT_SUCCESS;
}
