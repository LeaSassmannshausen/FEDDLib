#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/Checkpointing/CheckpointMetadata.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace FEDD;
using Teuchos::rcp;
using NavierStokesType = NavierStokes<>;
using MultiVectorType = MultiVector<>;

namespace {

void zero2D(double*, double* value, double, const double*)
{
    value[0] = value[1] = 0.;
}

void zero3D(double*, double* value, double, const double*)
{
    value[0] = value[1] = value[2] = 0.;
}

void inlet2D(double* x, double* value, double, const double* parameters)
{
    value[0] = 4. * parameters[0] * x[1] * (1. - x[1]);
    value[1] = 0.;
}

void inlet3D(double* x, double* value, double, const double* parameters)
{
    value[0] = 16. * parameters[0] * x[1] * (1. - x[1]) * x[2] * (1. - x[2]);
    value[1] = value[2] = 0.;
}

/** Compare against an independently read field, with exact checks at initialization. */
void compareFields(NavierStokesType& problem, const std::string& directory,
                   double time, double tolerance, const std::string& label,
                   const std::string& actualDirectory = "")
{
    for (int block = 0; block < 2; ++block) {
        const std::string field = block == 0 ? "u" : "p";
        auto actual = problem.getSolution()->getBlock(block);
        if (!actualDirectory.empty()) {
            HDF5Import<> output(actual->getMap(), joinPath(actualDirectory, "Solution" + field));
            actual = output.readVariablesHDF5(std::to_string(time));
        }
        HDF5Import<> input(actual->getMap(), joinPath(directory, "Solution" + field));
        auto expected = input.readVariablesHDF5(std::to_string(time));
        MultiVectorType difference(actual->getMap());
        difference.update(1., *actual, -1., *expected, 0.);
        Teuchos::Array<double> error(1), referenceNorm(1);
        difference.norm2(error);
        expected->norm2(referenceNorm);
        const double relative = referenceNorm[0] > 0. ? error[0] / referenceNorm[0] : error[0];
        const bool passed = std::isfinite(error[0]) && std::isfinite(referenceNorm[0]) &&
                            error[0] <= tolerance * referenceNorm[0];
        if (problem.comm_->getRank() == 0)
            std::cout << std::setprecision(16) << label << " (" << (block == 0 ? "velocity" : "pressure")
                      << "): absolute l2 error " << error[0] << ", relative error " << relative
                      << " (tolerance " << tolerance << "): " << (passed ? "PASS" : "FAIL") << '\n';
        TEUCHOS_TEST_FOR_EXCEPTION(!passed, std::runtime_error, label << ": field " << field << " differs");
    }
}

} // namespace

int main(int argc, char** argv)
{
    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession session(&argc, &argv, &blackhole);
    auto comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    std::string problemFile = "parametersProblem2D.xml", phase = "load";
    bool perturbMesh = false;
    Teuchos::CommandLineProcessor commandLine;
    commandLine.setOption("problemfile", &problemFile, "2D or 3D mesh parameters.");
    commandLine.setOption("phase", &phase, "source, manual, validate, or load.");
    commandLine.setOption("perturb-mesh", "original-mesh", &perturbMesh, "Modify a coordinate before compatibility validation.");
    commandLine.throwExceptions(false);
    const auto result = commandLine.parse(argc, argv);
    if (result == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return EXIT_SUCCESS;
    if (result != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) return EXIT_FAILURE;
    try {
        TEUCHOS_TEST_FOR_EXCEPTION(phase != "source" && phase != "manual" && phase != "validate" && phase != "load",
                                   std::logic_error, "Unknown initial solution test phase");
        auto parameters = Teuchos::getParametersFromXmlFile("parametersProblem.xml");
        parameters->setParameters(*Teuchos::getParametersFromXmlFile(problemFile));
        parameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersPrec.xml"));
        parameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersSolver.xml"));
        auto& time = parameters->sublist("Timestepping Parameter");
        const double sourceTime = time.get<double>("Initial solution time");
        const std::string sourceDirectory = time.get<std::string>("Initial solution directory");
        time.set("Initial solution", phase == "load" || phase == "validate");
        time.set("Checkpoint directory", phase == "source" ? sourceDirectory : phase == "manual" ? "manualReference" : "initialRun");
        if (phase == "source") {
            // Write a steady state with only primary fields. Its integration
            // metadata differs from the new run, and it contains no BDF history.
            time.set("dt", 0.0025).set("BDF", 1);
            parameters->sublist("Parameter").set("Viscosity", 0.01);
        }
        if (phase == "validate") time.set("Final time", 0.);
        checkpoint::onRoot(*comm, [&] {
            std::filesystem::create_directories(time.get<std::string>("Checkpoint directory"));
        });

        const int dim = parameters->sublist("Parameter").get<int>("Dimension");
        TEUCHOS_TEST_FOR_EXCEPTION(dim != 2 && dim != 3, std::logic_error, "Expected 2D or 3D BFS geometry");
        parameters->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").set("DofsPerNode1", dim);
        auto domain = rcp(new Domain<>(comm, dim));
        MeshPartitioner<>::DomainPtrArray_Type domains(1);
        domains[0] = domain;
        auto partitionerSettings = Teuchos::sublist(parameters, "Mesh Partitioner");
        partitionerSettings->set("Build Edge List", true).set("Build Surface List", true);
        MeshPartitioner<> partitioner(domains, partitionerSettings, "P1", dim);
        partitioner.readAndPartition();
        if (perturbMesh && domain->getMapUnique()->getNodeNumElements() > 0 && domain->getMapUnique()->getGlobalElement(0) == 0)
            domain->getPointsUnique()->at(0).at(0) += 0.001;
        auto boundaries = rcp(new BCBuilder<>());
        auto zero = dim == 2 ? zero2D : zero3D;
        boundaries->addBC(zero, 1, 0, domain, "Dirichlet", dim);
        boundaries->addBC(zero, 4, 0, domain, "Dirichlet", dim);
        std::vector<double> inletParameters{parameters->sublist("Parameter").get<double>("MaxVelocity")};
        boundaries->addBC(dim == 2 ? inlet2D : inlet3D, 2, 0, domain, "Dirichlet", dim,
                          inletParameters);
        NavierStokesType problem(domain, "P1", domain, "P1", parameters);
        problem.addBoundaries(boundaries);
        problem.initializeProblem();
        if (phase == "manual") {
            // Independent baseline: assign only field values through the vector
            // interface, without invoking any checkpoint initialization logic.
            for (int block = 0; block < 2; ++block) {
                auto destination = problem.getSolution()->getBlockNonConst(block);
                HDF5Import<> input(destination->getMap(), joinPath(sourceDirectory, block == 0 ? "Solutionu" : "Solutionp"));
                auto initial = input.readVariablesHDF5(std::to_string(sourceTime));
                destination->update(1., *initial, 0.);
            }
        }
        if (phase != "source") compareFields(problem, sourceDirectory, sourceTime, 0., "Initial field loading");
        problem.assemble();
        problem.setBoundariesRHS();
        if (phase == "source") {
            problem.calculateNonLinResidualVec("reverse");
            const double initialResidual = problem.calculateResidualNorm();
            NonLinearSolver<> nonlinearSolver("Newton");
            nonlinearSolver.solve(problem);
            problem.calculateNonLinResidualVec("reverse");
            const double residual = problem.calculateResidualNorm();
            TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(residual) ||
                                       !(residual <= parameters->sublist("Parameter").get<double>("relNonLinTol") * initialResidual),
                                       std::runtime_error, "Initial steady solution did not converge");
            for (int block = 0; block < 2; ++block) {
                auto values = problem.getSolution()->getBlock(block);
                HDF5Export<> output(values->getMap(), joinPath(sourceDirectory, block == 0 ? "Solutionu" : "Solutionp"));
                output.writeVariablesHDF5(std::to_string(sourceTime), values);
                output.closeExporter();
            }
            problem.writeCheckpointMetadata(sourceTime);
            if (comm->getRank() == 0) std::cout << "PASS converged steady source with primary fields only\n";
            return EXIT_SUCCESS;
        }

        SmallMatrix<int> timeBlocks(2);
        timeBlocks[0][0] = timeBlocks[0][1] = 1;
        DAESolverInTime<> solver(parameters, comm);
        solver.defineTimeStepping(timeBlocks);
        solver.setProblem(problem);
        solver.setupTimeStepping();
        TEUCHOS_TEST_FOR_EXCEPTION(solver.timeSteppingTool_->currentTime() != 0. ||
                                   !solver.problemTime_->solutionPreviousTimesteps_.empty(),
                                   std::runtime_error, "Initial solution restored a clock or timestep history");
        compareFields(problem, sourceDirectory, sourceTime, 0., "Fields after timestep setup");
        solver.advanceInTime();
        TEUCHOS_TEST_FOR_EXCEPTION(std::abs(solver.timeSteppingTool_->currentTime() - time.get<double>("Final time")) > 1.e-14,
                                   std::runtime_error, "Wrong final time for new simulation");
        if (phase == "validate") {
            TEUCHOS_TEST_FOR_EXCEPTION(!solver.problemTime_->solutionPreviousTimesteps_.empty(),
                                       std::runtime_error, "Zero-duration initialization advanced a timestep");
            compareFields(problem, sourceDirectory, sourceTime, 0., "Zero-duration initialization");
            if (comm->getRank() == 0) std::cout << "PASS initialization at time zero without advancing or solving\n";
        }
        if (phase == "load") {
            const double dt = time.get<double>("dt");
            for (double stateTime : {0., dt, 2. * dt}) {
                // Check the written t=0 state and both startup steps, including
                // the first BDF1 step and the subsequent BDF2 step.
                compareFields(problem, "manualReference", stateTime,
                              stateTime == 0. ? 0. : time.get<double>("Initial solution tolerance"),
                              "Direct initialization comparison at t=" + std::to_string(stateTime), "initialRun");
            }
        }
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error) {
        if (comm->getRank() == 0) std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
