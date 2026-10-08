#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Tpetra_Core.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace FEDD;
using Teuchos::RCP;
using Teuchos::rcp;

using SC = default_sc;
using LO = default_lo;
using GO = default_go;
using NO = default_no;
using Domain_Type = Domain<SC, LO, GO, NO>;
using MeshPartitioner_Type = MeshPartitioner<SC, LO, GO, NO>;
using MultiVector_Type = MultiVector<SC, LO, GO, NO>;

namespace {

void load2D(double*, double* result, double* parameters)
{
    result[0] = 0.;
    result[1] = parameters[1];
}

void load3D(double*, double* result, double* parameters)
{
    result[0] = parameters[1];
    result[1] = result[2] = 0.;
}

void zeroDirichlet2D(double*, double* result, double, const double*)
{
    result[0] = result[1] = 0.;
}

void zeroDirichlet3D(double*, double* result, double, const double*)
{
    result[0] = result[1] = result[2] = 0.;
}

} // namespace

/**
 * @brief Compare steady nonlinear elasticity against the stored P1 references.
 *
 * The registered cases use structured unit squares/cubes with H/h = 6 on four
 * ranks (2D) or eight ranks (3D). Loading acts in y in 2D and x in 3D.
 * The comparison retains the absolute infinity-norm tolerance of 1e-11.
 * @author Christian Hochmuth
 */
int main(int argc, char* argv[])
{
    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    auto comm = Tpetra::getDefaultComm();

    std::string problemFile = "parametersProblemNonLinElasticity.xml";
    std::string precFile = "parametersPrecNonLinElasticity.xml";
    std::string solverFile = "parametersSolver.xml";
    int dim = 2;
    Teuchos::CommandLineProcessor commandLine;
    commandLine.setOption("problemfile", &problemFile, "Problem parameters.");
    commandLine.setOption("precfile", &precFile, "Preconditioner parameters.");
    commandLine.setOption("solverfile", &solverFile, "Linear solver parameters.");
    commandLine.setOption("dim", &dim, "Dimension: 2 or 3.");
    commandLine.throwExceptions(false);
    const auto parseResult = commandLine.parse(argc, argv);
    if (parseResult == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
        return EXIT_SUCCESS;
    if (parseResult != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
        return EXIT_FAILURE;
    TEUCHOS_TEST_FOR_EXCEPTION(dim != 2 && dim != 3, std::logic_error,
                               "Nonlinear elasticity unit tests require dimension 2 or 3.");

    auto parameters = Teuchos::getParametersFromXmlFile(problemFile);
    parameters->setParameters(*Teuchos::getParametersFromXmlFile(precFile));
    parameters->setParameters(*Teuchos::getParametersFromXmlFile(solverFile));
    // Assembly and the preconditioner must use the selected test dimension.
    parameters->sublist("Parameter").set("Dimension", dim);
    parameters->sublist("ThyraPreconditioner").sublist("Preconditioner Types")
        .sublist("FROSch").set("DofsPerNode1", dim);

    const std::string feType = parameters->sublist("Parameter").get("Discretization", "P1");
    const std::string meshType = parameters->sublist("Parameter").get("Mesh Type", "structured");
    const int coarseRanks = parameters->sublist("General").get("Mpi Ranks Coarse", 0);
    const int size = comm->getSize() - coarseRanks;
    RCP<Domain_Type> domain;
    if (meshType == "structured") {
        const int subdomainsPerAxis = static_cast<int>(std::pow(size, 1. / dim)
            + 100. * Teuchos::ScalarTraits<double>::eps());
        const int resolution = parameters->sublist("Parameter").get("H/h", 5);
        std::vector<double> origin(dim, 0.);
        domain = dim == 2 ? rcp(new Domain_Type(origin, 1., 1., comm))
                          : rcp(new Domain_Type(origin, 1., 1., 1., comm));
        domain->buildMesh(1, "Square", dim, feType, subdomainsPerAxis, resolution, coarseRanks);
    } else if (meshType == "unstructured") {
        auto linearDomain = rcp(new Domain_Type(comm, dim));
        MeshPartitioner_Type::DomainPtrArray_Type domains(1);
        domains[0] = linearDomain;
        MeshPartitioner_Type partitioner(domains, Teuchos::sublist(parameters, "Mesh Partitioner"), "P1", dim);
        partitioner.readAndPartition();
        domain = linearDomain;
        if (feType == "P2") {
            domain = rcp(new Domain_Type(comm, dim));
            domain->buildP2ofP1Domain(linearDomain);
        }
    } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unsupported mesh type: " << meshType);
    }

    auto boundaries = rcp(new BCBuilder<SC, LO, GO, NO>());
    boundaries->addBC(dim == 2 ? zeroDirichlet2D : zeroDirichlet3D,
                      meshType == "structured" ? 2 : 1, 0, domain, "Dirichlet", dim);
    NonLinElasticity<SC, LO, GO, NO> elasticity(domain, feType, parameters);
    elasticity.addBoundaries(boundaries);
    elasticity.addRhsFunction(dim == 2 ? load2D : load3D);
    elasticity.addParemeterRhs(parameters->sublist("Parameter").get("Volume force", 0.));
    elasticity.addParemeterRhs(0.); // Constant load quadrature degree.
    elasticity.initializeProblem();
    elasticity.assemble();
    elasticity.setBoundaries();
    NonLinearSolver<SC, LO, GO, NO> solver(parameters->sublist("General").get("Linearization", "Newton"));
    solver.solve(elasticity);

    auto solution = elasticity.getSolution()->getBlock(0);
    const std::string referenceFile = "ReferenceSolutions/solution_nonLinElasticity_"
        + std::to_string(dim) + "d_" + feType + "_" + std::to_string(size) + "cores";
    HDF5Import<SC, LO, GO, NO> importer(solution->getMap(), referenceFile);
    auto reference = importer.readVariablesHDF5("solution");
    MultiVector_Type error(solution->getMap());
    error.update(1., *solution, -1., *reference, 0.);
    Teuchos::Array<SC> norm(1);
    error.normInf(norm);
    constexpr double tolerance = 1.e-11;
    if (comm->getRank() == 0)
        std::cout << "Nonlinear elasticity reference error (infinity): " << norm[0]
                  << " (tolerance " << tolerance << ")" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(!(norm[0] <= tolerance), std::logic_error,
                               "Nonlinear elasticity reference comparison failed: " << norm[0]);
    return EXIT_SUCCESS;
}
