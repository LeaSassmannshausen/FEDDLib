#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/General/CheckpointFiles.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_DefaultPlatform.hpp>
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
typedef Domain<SC,LO,GO,NO> Domain_Type;
typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
typedef NavierStokes<SC,LO,GO,NO> NavierStokes_Type;
typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;

namespace {

void zeroDirichlet2D(double*, double* result, double, const double*)
{
    result[0] = 0.;
    result[1] = 0.;
}

void zeroDirichlet3D(double*, double* result, double, const double*)
{
    result[0] = 0.;
    result[1] = 0.;
    result[2] = 0.;
}

void inflow2D(double* x, double* result, double, const double* parameters)
{
    const double height = parameters[1];
    result[0] = 4. * parameters[0] * x[1] * (height - x[1]) / (height * height);
    result[1] = 0.;
}

void inflow3D(double* x, double* result, double, const double* parameters)
{
    const double height = parameters[1];
    result[0] = 16. * parameters[0] * x[1] * (height - x[1]) * x[2] * (height - x[2])
                / (height * height * height * height);
    result[1] = 0.;
    result[2] = 0.;
}


bool compareRestart(NavierStokes_Type& problem, ParameterListPtr_Type parameters,
                    RCP<const Teuchos::Comm<int>> comm)
{
    const auto& timeParameters = parameters->sublist("Timestepping Parameter");
    const double finalTime = timeParameters.get<double>("Final time");
    const double tolerance = timeParameters.get<double>("Restart tolerance");
    const char* checkpointNames[] = {"Solutionu", "Solutionp"};
    const char* fieldNames[] = {"fluid velocity", "fluid pressure"};
    bool passed = true;

    for (int block = 0; block < 2; ++block) {
        auto solution = problem.getSolution()->getBlock(block);
        const std::string referenceFile = timeParameters.isParameter("Reference directory")
            ? joinPath(timeParameters.get<std::string>("Reference directory"), checkpointNames[block])
            : restartFile(parameters, checkpointNames[block]);
        HDF5Import<SC,LO,GO,NO> importer(solution->getMap(), referenceFile);
        auto reference = importer.readVariablesHDF5(std::to_string(finalTime));
        MultiVector_Type error(solution->getMap());
        error.update(1., *solution, -1., *reference, 0.);

        Teuchos::Array<SC> errorNorm(1), solutionNorm(1);
        error.norm2(errorNorm);
        solution->norm2(solutionNorm);
        const double relativeError = errorNorm[0] / solutionNorm[0];
        if (comm->getRank() == 0) {
            std::cout << "Restart absolute error (" << fieldNames[block] << ", l2): "
                      << errorNorm[0] << std::endl;
            std::cout << "Restart relative error (" << fieldNames[block] << "): "
                      << relativeError << " (tolerance " << tolerance << ")" << std::endl;
        }
        if (!(relativeError <= tolerance))
            passed = false;
    }
    return passed;
}

} // namespace

// Test for unsteady Navier Stokes restart. 
// The test runs a short simulation, writes checkpoints, and then restarts from the first checkpoint to compare against the original run. 
// The test passes if the restarted solution matches the original within a specified tolerance.
// Test Case: 2D and 3D BFS flow with P1-P1 discretization, using BDF2 time stepping.
int main(int argc, char* argv[])
{
    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);
    auto comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    std::string problemFile = "parametersProblem2D.xml";
    std::string overrideFile;
    bool validateOnly = false;
    bool perturbMesh = false;
    Teuchos::CommandLineProcessor commandLine;
    commandLine.setOption("problemfile", &problemFile, "2D or 3D BFS case parameters.");
    commandLine.setOption("overridefile", &overrideFile, "Optional parameter overrides for this test phase.");
    commandLine.setOption("restartfile", &overrideFile, "Alias for --overridefile for restart runs.");
    commandLine.setOption("validate-checkpoint", "solve", &validateOnly, "Validate/restore without advancing time.");
    commandLine.setOption("perturb-mesh", "original-mesh", &perturbMesh, "Change a coordinate for compatibility rejection testing.");
    commandLine.throwExceptions(false);
    const auto parseResult = commandLine.parse(argc, argv);
    if (parseResult == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
        return EXIT_SUCCESS;
    if (parseResult != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
        return EXIT_FAILURE;

    // All phases use the same case and solver settings, with separate run/restart overrides.
    auto parameters = Teuchos::getParametersFromXmlFile("parametersProblem.xml");
    parameters->setParameters(*Teuchos::getParametersFromXmlFile(problemFile));
    if (!overrideFile.empty())
        parameters->setParameters(*Teuchos::getParametersFromXmlFile(overrideFile));
    parameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersPrec.xml"));
    parameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersSolver.xml"));

    const int dim = parameters->sublist("Parameter").get<int>("Dimension");
    TEUCHOS_TEST_FOR_EXCEPTION(dim != 2 && dim != 3, std::logic_error, "BFS restart tests require dimension 2 or 3.");
    parameters->sublist("ThyraPreconditioner").sublist("Preconditioner Types")
        .sublist("FROSch").set("DofsPerNode1", dim);
    const bool restart = parameters->sublist("Timestepping Parameter").get<bool>("Restart");
    if (comm->getRank() == 0)
        std::cout << dim << "D BFS restart test: " << (restart ? "restart" : "uninterrupted reference")
                  << " on " << comm->getSize() << " ranks, from "
                  << parameters->sublist("Timestepping Parameter").get("Time step", 0.)
                  << " to " << parameters->sublist("Timestepping Parameter").get<double>("Final time")
                  << std::endl;

    // Definition of differen domains for velocity and pressure, with a single partitioner for both. 
    // The velocity domain is used to define the boundary conditions.                  
    RCP<Domain_Type> pressureDomain = rcp(new Domain_Type(comm, dim));
    RCP<Domain_Type> velocityDomain = rcp(new Domain_Type(comm, dim));
    MeshPartitioner_Type::DomainPtrArray_Type domains(1);
    domains[0] = pressureDomain;
    auto partitionerParameters = Teuchos::sublist(parameters, "Mesh Partitioner");
    partitionerParameters->set("Build Edge List", true);
    partitionerParameters->set("Build Surface List", true);
    MeshPartitioner_Type partitioner(domains, partitionerParameters, "P1", dim);
    partitioner.readAndPartition();
    if (perturbMesh && pressureDomain->getMapUnique()->getNodeNumElements() > 0 &&
        pressureDomain->getMapUnique()->getGlobalElement(0) == 0)
        pressureDomain->getPointsUnique()->at(0).at(0) += 0.001;
    velocityDomain = pressureDomain; // We use only P1-P1 discretization for this test, so the velocity and pressure domains are the same.

    // Define boundary conditions for the velocity field. 
    RCP<BCBuilder<SC,LO,GO,NO>> boundaries = rcp(new BCBuilder<SC,LO,GO,NO>());
    std::vector<double> inflowParameters = {parameters->sublist("Parameter").get<double>("MaxVelocity"), 1.};
    auto zero = dim == 2 ? zeroDirichlet2D : zeroDirichlet3D;
    auto inflow = dim == 2 ? inflow2D : inflow3D;
    // 1 = walls, 2 = inlet, 3 = natural outflow; 4 is an optional obstacle boundary.
    boundaries->addBC(zero, 1, 0, velocityDomain, "Dirichlet", dim);
    boundaries->addBC(zero, 4, 0, velocityDomain, "Dirichlet", dim);
    boundaries->addBC(inflow, 2, 0, velocityDomain, "Dirichlet", dim, inflowParameters);

    NavierStokes_Type problem(velocityDomain, "P1", pressureDomain, "P1", parameters);
    problem.addBoundaries(boundaries);
    problem.initializeProblem();
    if (validateOnly) {
        if (comm->getRank() == 0) std::cout << "Checkpoint compatibility validation passed." << std::endl;
        return EXIT_SUCCESS;
    }
    problem.assemble();
    problem.setBoundariesRHS();

    SmallMatrix<int> timeBlocks(2);
    timeBlocks[0][0] = 1;
    timeBlocks[0][1] = 1;
    DAESolverInTime<SC,LO,GO,NO> timeSolver(parameters, comm);
    timeSolver.defineTimeStepping(timeBlocks);
    timeSolver.setProblem(problem);
    timeSolver.setupTimeStepping();
    timeSolver.advanceInTime(); // Time stepping is handled by the solver, which calls the problem's assemble and solve methods as needed.

    if (restart)
        return compareRestart(problem, parameters, comm) ? EXIT_SUCCESS : EXIT_FAILURE;
    return EXIT_SUCCESS;
}
