#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/General/CheckpointFiles.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Tpetra_Core.hpp>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
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
typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
typedef TimeProblem<SC,LO,GO,NO> TimeProblem_Type;

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

/**
 * @brief Compare the continued Newmark state with the uninterrupted reference.
 *
 * The current Newmark loop updates derivatives and writes checkpoints at the
 * beginning of the following step. Both phases therefore advance one extra
 * step. The latest displacement history and derivative buffers then all belong
 * to Comparison time; the primary displacement already belongs to the next step.
 * Compare these runtime buffers, rather than rereading the restarted input.
 */
bool compareRestart(TimeProblem_Type& problem, ParameterListPtr_Type parameters,
                    RCP<const Teuchos::Comm<int>> comm)
{
    const auto& timeParameters = parameters->sublist("Timestepping Parameter");
    const double time = timeParameters.get<double>("Comparison time");
    const double tolerance = timeParameters.get<double>("Restart tolerance");
    const char* checkpointNames[] = {"Solutiond_s", "ds_Velocity", "ds_Acceleration"};
    const char* fieldNames[] = {"solid displacement", "solid velocity", "solid acceleration"};
    const auto states = {problem.getSolutionPreviousTimestep()->getBlock(0),
                         problem.velocityPreviousTimesteps_.at(0)->getBlock(0),
                         problem.accelerationPreviousTimesteps_.at(0)->getBlock(0)};
    bool passed = true;
    int field = 0;
    for (const auto& solution : states) {
        HDF5Import<SC,LO,GO,NO> importer(solution->getMap(), restartFile(parameters, checkpointNames[field]));
        auto reference = importer.readVariablesHDF5(std::to_string(time));
        MultiVector_Type error(solution->getMap());
        error.update(1., *solution, -1., *reference, 0.);
        Teuchos::Array<SC> errorNorm(1), referenceNorm(1);
        error.norm2(errorNorm);
        reference->norm2(referenceNorm);
        const double relativeError = referenceNorm[0] > 0. ? errorNorm[0] / referenceNorm[0] : errorNorm[0];
        // Every loaded field must be nonzero: otherwise a missing surface load
        // could make a trivial all-zero simulation look like a successful restart.
        if (!(referenceNorm[0] > 0.) || !std::isfinite(relativeError) || relativeError > tolerance)
            passed = false;
        if (comm->getRank() == 0) {
            std::cout << std::setprecision(16)
                      << "Reference norm (" << fieldNames[field] << ", l2): " << referenceNorm[0] << '\n'
                      << "Restart absolute error (" << fieldNames[field] << ", l2): " << errorNorm[0] << '\n'
                      << "Restart relative error (" << fieldNames[field] << "): " << relativeError
                      << " (tolerance " << tolerance << ")" << std::endl;
        }
        ++field;
    }
    return passed;
}

} // namespace

int main(int argc, char* argv[])
{
    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    auto comm = Tpetra::getDefaultComm();
    std::string restartFileName;
    Teuchos::CommandLineProcessor commandLine;
    commandLine.setOption("restartfile", &restartFileName, "Restart parameter overrides.");
    commandLine.throwExceptions(false);
    const auto parseResult = commandLine.parse(argc, argv);
    if (parseResult == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
        return EXIT_SUCCESS;
    if (parseResult != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
        return EXIT_FAILURE;

    auto parameters = Teuchos::getParametersFromXmlFile("parametersProblem.xml");
    if (!restartFileName.empty())
        parameters->setParameters(*Teuchos::getParametersFromXmlFile(restartFileName));
    parameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersPrec.xml"));
    parameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersSolver.xml"));
    const bool restart = parameters->sublist("Timestepping Parameter").get<bool>("Restart");
    if (comm->getRank() == 0)
        std::cout << "2D nonlinear elasticity Newmark restart test: "
                  << (restart ? "restart" : "uninterrupted reference")
                  << " on " << comm->getSize() << " ranks" << std::endl;

    RCP<Domain_Type> linearDomain = rcp(new Domain_Type(comm, 2));
    MeshPartitioner_Type::DomainPtrArray_Type domains(1);
    domains[0] = linearDomain;
    MeshPartitioner_Type partitioner(domains, Teuchos::sublist(parameters, "Mesh Partitioner"), "P1", 2);
    partitioner.readAndPartition();
    RCP<Domain_Type> domain = rcp(new Domain_Type(comm, 2));
    domain->buildP2ofP1Domain(linearDomain);
    domain->preProcessMesh(true, true);

    RCP<BCBuilder<SC,LO,GO,NO>> boundaries = rcp(new BCBuilder<SC,LO,GO,NO>());
    boundaries->addBC(zeroDirichlet, 4, 0, domain, "Dirichlet", 2);
    NonLinElasticity<SC,LO,GO,NO> problem(domain, "P2", parameters);
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
    if (restart)
        return compareRestart(*timeSolver.problemTime_, parameters, comm) ? EXIT_SUCCESS : EXIT_FAILURE;
    return EXIT_SUCCESS;
}
