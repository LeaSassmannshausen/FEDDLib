#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/Checkpointing/CheckpointFiles.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "TransientReference.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/Preconditioner.hpp"
#include "feddlib/problems/specific/FSI.hpp"
#include "feddlib/problems/specific/Geometry.hpp"
#include "feddlib/problems/specific/LinElas.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Tpetra_Core.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace FEDD;
using Teuchos::RCP;
using Teuchos::rcp;

typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;
typedef Domain<SC,LO,GO,NO> Domain_Type;
typedef RCP<Domain_Type> DomainPtr_Type;
typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
typedef BCBuilder<SC,LO,GO,NO> BCBuilder_Type;
typedef FSI<SC,LO,GO,NO> FSI_Type;
typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;

namespace {

void zeroDirichlet(double*, double* result, double, const double*)
{
    result[0] = 0.;
    result[1] = 0.;
}

void inflow(double* x, double* result, double time, const double* parameters)
{
    const double height = parameters[1];
    const double ramp = time < 0.5 ? 0.5 * (1. - std::cos(2. * M_PI * time)) : 1.;
    result[0] = 4. * 1.5 * parameters[0] * x[1] * (height - x[1]) / (height * height) * ramp;
    result[1] = 0.;
}

} // namespace

/** @brief Run the 2D Turek FSI case and compare its state at t = 0.02 with frozen references. */
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
    auto parameters = Teuchos::getParametersFromXmlFile("parametersProblemFSI.xml");
    // Match the restart case at its comparison time, without checkpoint I/O.
    parameters->sublist("Timestepping Parameter").set("Checkpointing", false);
    parameters->sublist("Timestepping Parameter").set("Restart", false);
    parameters->sublist("Timestepping Parameter").set("Final time", 0.02);
    parameters->sublist("General").set("Safe all solution", false);
    parameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersSolverFSI.xml"));

    auto fluidParameters = Teuchos::getParametersFromXmlFile("parametersPrecFluidMono.xml");
    fluidParameters->sublist("Parameter").setParameters(parameters->sublist("Parameter Fluid"));
    fluidParameters->sublist("Timestepping Parameter").setParameters(parameters->sublist("Timestepping Parameter"));

    auto solidParameters = Teuchos::getParametersFromXmlFile("parametersPrecStructure.xml");
    solidParameters->sublist("Parameter").setParameters(parameters->sublist("Parameter Solid"));
    solidParameters->sublist("Timestepping Parameter").setParameters(parameters->sublist("Timestepping Parameter"));

    auto geometryParameters = Teuchos::getParametersFromXmlFile("parametersPrecGeometry.xml");
    geometryParameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersSolverGeometry.xml"));
    geometryParameters->sublist("General").set("Preconditioner Method", "MonolithicConstPrec");
    geometryParameters->sublist("Parameter").set("Model", "Laplace");

    // This test uses the 2D Turek case: P2 velocity/displacement and P1 pressure.
    const int dim = 2;
    const std::string discType = "P2";
    DomainPtr_Type fluidP1 = rcp(new Domain_Type(comm, dim));
    DomainPtr_Type solidP1 = rcp(new Domain_Type(comm, dim));
    MeshPartitioner_Type::DomainPtrArray_Type domains(2);
    domains[0] = fluidP1;
    domains[1] = solidP1;
    auto partitionerParameters = Teuchos::sublist(parameters, "Mesh Partitioner");
    partitionerParameters->set("Build Edge List", true);
    partitionerParameters->set("Build Surface List", true);
    MeshPartitioner_Type partitioner(domains, partitionerParameters, "P1", dim);
    partitioner.readAndPartition();

    DomainPtr_Type fluidP2 = rcp(new Domain_Type(comm, dim));
    DomainPtr_Type solidP2 = rcp(new Domain_Type(comm, dim));
    fluidP2->buildP2ofP1Domain(fluidP1);
    solidP2->buildP2ofP1Domain(solidP1);
    vec_int_Type interfaceIds(1, 5);
    fluidP1->identifyInterfaceParallelAndDistance(solidP1, interfaceIds);
    fluidP2->identifyInterfaceParallelAndDistance(solidP2, interfaceIds);
    fluidP2->buildInterfaceMaps();
    solidP2->buildInterfaceMaps();

    DomainPtr_Type interfaceDomain = rcp(new Domain_Type(comm));
    interfaceDomain->setDummyInterfaceDomain(fluidP2);
    fluidP2->setReferenceConfiguration();
    fluidP1->setReferenceConfiguration();

    RCP<SmallMatrix<int>> timeBlocks = rcp(new SmallMatrix<int>(4));
    (*timeBlocks)[0][0] = 1;
    (*timeBlocks)[0][1] = 1;
    (*timeBlocks)[2][2] = 1;
    FSI_Type fsi(fluidP2, discType, fluidP1, "P1", solidP2, discType,
                 interfaceDomain, discType, fluidP2, discType,
                 fluidParameters, solidParameters, parameters, geometryParameters, timeBlocks);

    RCP<BCBuilder_Type> boundaries = rcp(new BCBuilder_Type());
    RCP<BCBuilder_Type> fluidBoundaries = rcp(new BCBuilder_Type());
    RCP<BCBuilder_Type> solidBoundaries = rcp(new BCBuilder_Type());
    RCP<BCBuilder_Type> geometryBoundaries = rcp(new BCBuilder_Type());
    RCP<BCBuilder_Type> interfaceBoundaries = rcp(new BCBuilder_Type());
    std::vector<double> inflowParameters = {
        parameters->sublist("Parameter").get<double>("MeanVelocity"), 0.41
    };

    // Fluid flags: 1 = wall, 2 = inflow, 3 = outflow, 4 = obstacle, 5 = interface.
    for (int flag : {1, 4}) {
        boundaries->addBC(zeroDirichlet, flag, 0, fluidP2, "Dirichlet", dim);
        fluidBoundaries->addBC(zeroDirichlet, flag, 0, fluidP2, "Dirichlet", dim);
    }
    boundaries->addBC(inflow, 2, 0, fluidP2, "Dirichlet", dim, inflowParameters);
    fluidBoundaries->addBC(inflow, 2, 0, fluidP2, "Dirichlet", dim, inflowParameters);
    boundaries->addBC(zeroDirichlet, 1, 2, solidP2, "Dirichlet", dim);
    solidBoundaries->addBC(zeroDirichlet, 1, 0, solidP2, "Dirichlet", dim);
    for (int flag = 1; flag <= 5; ++flag)
        geometryBoundaries->addBC(zeroDirichlet, flag, 0, fluidP2, "Dirichlet", dim);
    interfaceBoundaries->addBC(zeroDirichlet, 5, 0, fluidP2, "Dirichlet", dim);

    // The subproblems also need boundary conditions when forming the time integration RHS.
    fsi.problemFluid_->addBoundaries(fluidBoundaries);
    fsi.problemStructure_->addBoundaries(solidBoundaries);
    fsi.problemGeometry_->addBoundaries(geometryBoundaries);
    fsi.getPreconditioner()->setFaCSIBCFactory(interfaceBoundaries);
    fsi.addBoundaries(boundaries);
    fsi.initializeProblem();
    fsi.initializeGE();
    fsi.assemble();

    DAESolverInTime<SC,LO,GO,NO> timeSolver(parameters, comm);
    timeSolver.defineTimeStepping(*timeBlocks);
    timeSolver.setProblem(fsi);
    timeSolver.setupTimeStepping();
    timeSolver.advanceInTime();

    bool passed = true;
    const char* files[] = {"solution_fsi_velocity_2d_P2_4cores",
                           "solution_fsi_pressure_2d_P1_4cores",
                           "solution_fsi_displacement_2d_P2_4cores"};
    const char* fields[] = {"fluid velocity", "fluid pressure", "solid displacement"};
    for (int block = 0; block < 3; ++block)
        passed = TransientReference::check(fsi.getSolution()->getBlock(block), files[block],
                    fields[block], referenceDirectory, writeReference) && passed;
    return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
