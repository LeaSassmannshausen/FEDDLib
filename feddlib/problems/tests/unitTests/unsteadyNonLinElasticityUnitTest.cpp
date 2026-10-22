#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/Checkpointing/CheckpointFiles.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "TransientReference.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Tpetra_Core.hpp>
#include <cstdlib>
#include <stdexcept>
#include <string>

using namespace FEDD;
using Teuchos::RCP;
using Teuchos::rcp;

typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;
typedef Domain<SC,LO,GO,NO> Domain_Type;
typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;

namespace {

void zeroDirichlet(double*, double* result, double, const double*)
{
    result[0] = result[1] = 0.;
}

// The surface assembler supplies time at index 0 and the boundary flag at index 5.
void surfaceLoad(double*, double* result, double* parameters)
{
    result[0] = parameters[5] == 2. ? parameters[1] : 0.;
    result[1] = 0.;
}

} // namespace

/**
 * @brief Compare the square_h02 Newmark case at t = 0.02 with frozen references.
 *
 * The extra step to 0.0225 finalizes the displacement/derivative history at 0.02
 * under the current start-of-step Newmark update convention, as in the restart test.
 */
int main(int argc, char* argv[])
{
    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    auto comm = Tpetra::getDefaultComm();
    std::string referenceDirectory = "ReferenceSolutions";
    bool writeReference = false;
    Teuchos::CommandLineProcessor commandLine;
    commandLine.setOption("reference-directory", &referenceDirectory, "Directory of stored unit-test references.");
    commandLine.setOption("write-reference", "compare-reference", &writeReference,
                          "Explicitly generate reference files instead of comparing them.");
    commandLine.throwExceptions(false);
    const auto parseResult = commandLine.parse(argc, argv);
    if (parseResult == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
        return EXIT_SUCCESS;
    if (parseResult != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
        return EXIT_FAILURE;
    auto parameters = Teuchos::getParametersFromXmlFile("parametersProblemUnsteadyNonLinElasticity.xml");
    parameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersPrecUnsteadyNonLinElasticity.xml"));
    parameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersSolverUnsteadyNonLinElasticity.xml"));
    parameters->sublist("Timestepping Parameter").set("Checkpointing", false);
    parameters->sublist("Timestepping Parameter").set("Restart", false);
    // References represent t=0.02. Keep that duration independent of the shorter
    // restart fixture, with one extra increment to finalize the Newmark buffers.
    parameters->sublist("Timestepping Parameter").set("Final time", 0.02 +
        parameters->sublist("Timestepping Parameter").get<double>("dt"));
    parameters->sublist("General").set("Safe all solution", false);

    RCP<Domain_Type> linearDomain = rcp(new Domain_Type(comm, 2));
    MeshPartitioner_Type::DomainPtrArray_Type domains(1);
    domains[0] = linearDomain;
    MeshPartitioner_Type partitioner(domains, Teuchos::sublist(parameters, "Mesh Partitioner"), "P1", 2);
    partitioner.readAndPartition();
    const std::string feType = parameters->sublist("Parameter").get<std::string>("Discretization");
    TEUCHOS_TEST_FOR_EXCEPTION(feType != "P1" && feType != "P2", std::logic_error,
                               "The elasticity reference test requires P1 or P2 elements.");
    RCP<Domain_Type> domain = linearDomain;
    if (feType == "P2") {
        domain = rcp(new Domain_Type(comm, 2));
        domain->buildP2ofP1Domain(linearDomain);
    }
    domain->preProcessMesh(true, true);

    RCP<BCBuilder<SC,LO,GO,NO>> boundaries = rcp(new BCBuilder<SC,LO,GO,NO>());
    boundaries->addBC(zeroDirichlet, 4, 0, domain, "Dirichlet", 2);
    NonLinElasticity<SC,LO,GO,NO> problem(domain, feType, parameters);
    problem.addRhsFunction(surfaceLoad);
    problem.addParemeterRhs(parameters->sublist("Parameter").get<double>("Volume force"));
    problem.addParemeterRhs(0.); // Constant load, no ramp.
    problem.addParemeterRhs(parameters->sublist("Timestepping Parameter").get<double>("dt"));
    problem.addParemeterRhs(0.); // Constant surface quadrature degree.
    problem.addBoundaries(boundaries);
    problem.initializeProblem();
    problem.assemble();

    SmallMatrix<int> timeBlocks(1);
    timeBlocks[0][0] = 1;
    DAESolverInTime<SC,LO,GO,NO> timeSolver(parameters, comm);
    timeSolver.defineTimeStepping(timeBlocks);
    timeSolver.setProblem(problem);
    timeSolver.setupTimeStepping();
    timeSolver.advanceInTime();
    // The Newmark update still runs before the next solve. At the end of the
    // extra step to 0.0225, these three buffers consistently represent t = 0.02.
    const auto states = {timeSolver.problemTime_->getSolutionPreviousTimestep()->getBlock(0),
                        timeSolver.problemTime_->velocityPreviousTimesteps_.at(0)->getBlock(0),
                        timeSolver.problemTime_->accelerationPreviousTimesteps_.at(0)->getBlock(0)};
    const std::string suffix = "_2d_" + feType + "_4cores";
    const std::string files[] = {"solution_unsteadyNonLinElasticity_displacement" + suffix,
                                "solution_unsteadyNonLinElasticity_velocity" + suffix,
                                "solution_unsteadyNonLinElasticity_acceleration" + suffix};
    const char* fields[] = {"solid displacement", "solid velocity", "solid acceleration"};
    bool passed = true;
    int field = 0;
    for (const auto& state : states) {
        passed = TransientReference::check(state, files[field], fields[field],
                                          referenceDirectory, writeReference) && passed;
        ++field;
    }
    return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
