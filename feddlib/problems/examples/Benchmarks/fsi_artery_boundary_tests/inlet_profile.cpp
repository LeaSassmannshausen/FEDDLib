#include <Tpetra_Core.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/Helper.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include <cmath>
#include <iostream>

using namespace FEDD;

// Solve -Delta_Gamma profile = 1 on the inlet, profile = 0 on its rim.
// Keep the full fluid map so the HDF5 solution can be imported directly by FSI.
int main(int argc, char** argv)
{
    Tpetra::ScopeGuard scope(&argc, &argv);
    auto comm = Tpetra::getDefaultComm();
    std::string mesh, feType = "P2", output;
    Teuchos::CommandLineProcessor options;
    options.setOption("mesh", &mesh, "Fluid .mesh file.");
    options.setOption("discretization", &feType, "P1 or P2.");
    options.setOption("output", &output, "Output basename, without .h5.");
    options.throwExceptions(false);
    const auto parsed = options.parse(argc, argv);
    if (parsed == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
    TEUCHOS_TEST_FOR_EXCEPTION(parsed != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ||
        mesh.empty() || output.empty() || (feType != "P1" && feType != "P2"), std::logic_error,
        "Supply --mesh, --discretization=P1|P2 and --output.");

    auto parameters = Teuchos::rcp(new Teuchos::ParameterList);
    parameters->set("Mesh 1 Name", mesh).set("Build Edge List", true).set("Build Surface List", true);
    auto linear = Teuchos::rcp(new Domain<>(comm, 3));
    MeshPartitioner<>::DomainPtrArray_Type domains(1, linear);
    MeshPartitioner<> partitioner(domains, parameters, "P1", 3);
    partitioner.readAndPartition(15, "mm", true);
    auto domain = linear;
    if (feType == "P2") {
        domain = Teuchos::rcp(new Domain<>(comm, 3));
        domain->buildP2ofP1Domain(linear);
    }
    const auto unique = domain->getMapUnique(), repeated = domain->getMapRepeated();
    const auto points = domain->getPointsRepeated();
    const auto flags = domain->getBCFlagRepeated();
    const auto elements = domain->getElementsC();
    auto matrix = Teuchos::rcp(new Matrix<>(unique, domain->getApproxEntriesPerRow()));
    auto rhsRepeated = Teuchos::rcp(new MultiVector<>(repeated));
    rhsRepeated->putScalar(0.);
    auto rhsData = rhsRepeated->getDataNonConst(0);
    vec2D_dbl_ptr_Type phi;
    vec3D_dbl_ptr_Type gradient;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type);
    Helper::getPhi(phi, weights, 2, feType, 2);
    Helper::getDPhi(gradient, weights, 2, feType, 2);
    long long localFaces = 0;
    for (UN element = 0; element < elements->numberElements(); ++element) {
        auto& volume = elements->getElement(element);
        for (int face = 0; face < volume.numSubElements(); ++face) {
            auto& surface = volume.getSubElements()->getElement(face);
            if (surface.getFlag() != 4 || volume.getSubElements()->getDimension() != 2) continue;
            ++localFaces;
            const auto nodes = surface.getVectorNodeList();
            double metric00 = 0., metric01 = 0., metric11 = 0.;
            for (int d = 0; d < 3; ++d) {
                const double a = (*points)[nodes[1]][d] - (*points)[nodes[0]][d];
                const double b = (*points)[nodes[2]][d] - (*points)[nodes[0]][d];
                metric00 += a*a; metric01 += a*b; metric11 += b*b;
            }
            const double determinant = metric00*metric11 - metric01*metric01;
            TEUCHOS_TEST_FOR_EXCEPTION(determinant <= 0., std::runtime_error, "Degenerate inlet triangle.");
            const double scaling = std::sqrt(determinant);
            for (UN i = 0; i < nodes.size(); ++i) {
                if ((*flags)[nodes[i]] != 4) continue; // Homogeneous Dirichlet rim.
                Teuchos::Array<default_go> columns;
                Teuchos::Array<double> entries;
                for (UN j = 0; j < nodes.size(); ++j) {
                    if ((*flags)[nodes[j]] != 4) continue;
                    double value = 0.;
                    for (UN q = 0; q < weights->size(); ++q) {
                        const auto& a = (*gradient)[q][i];
                        const auto& b = (*gradient)[q][j];
                        value += (*weights)[q] * (metric11*a[0]*b[0] - metric01*(a[0]*b[1]+a[1]*b[0])
                            + metric00*a[1]*b[1]) / determinant;
                    }
                    columns.push_back(repeated->getGlobalElement(nodes[j]));
                    entries.push_back(value*scaling);
                }
                matrix->insertGlobalValues(repeated->getGlobalElement(nodes[i]), columns(), entries());
                for (UN q = 0; q < weights->size(); ++q)
                    rhsData[nodes[i]] += scaling*(*weights)[q]*(*phi)[q][i];
            }
        }
    }
    long long faces = 0;
    Teuchos::reduceAll<int, long long>(*comm, Teuchos::REDUCE_SUM, 1, &localFaces, &faces);
    TEUCHOS_TEST_FOR_EXCEPTION(faces == 0, std::runtime_error, "No inlet faces with flag 4.");
    const auto uniqueFlags = domain->getBCFlagUnique();
    for (default_lo node = 0; node < unique->getNodeNumElements(); ++node) {
        if ((*uniqueFlags)[node] == 4) continue;
        Teuchos::Array<default_go> column(1, unique->getGlobalElement(node));
        Teuchos::Array<double> value(1, 1.);
        matrix->insertGlobalValues(column[0], column(), value());
    }
    matrix->fillComplete();
    auto rhs = Teuchos::rcp(new MultiVector<>(unique));
    rhs->putScalar(0.);
    rhs->exportFromVector(rhsRepeated, false, "Add");
    auto solution = Teuchos::rcp(new MultiVector<>(unique));
    solution->putScalar(0.);
    using MV = MultiVector<>::TpetraMultiVector_Type;
    using OP = Tpetra::Operator<default_sc, default_lo, default_go, default_no>;
    auto problem = Teuchos::rcp(new Belos::LinearProblem<double, MV, OP>(matrix->getTpetraMatrix(),
        solution->getTpetraMultiVectorNonConst(), rhs->getTpetraMultiVector()));
    TEUCHOS_TEST_FOR_EXCEPTION(!problem->setProblem(), std::runtime_error, "Cannot initialize inlet solve.");
    auto solverParameters = Teuchos::rcp(new Teuchos::ParameterList);
    solverParameters->set("Convergence Tolerance", 1.e-12).set("Maximum Iterations", 5000);
    Belos::PseudoBlockCGSolMgr<double, MV, OP> solver(problem, solverParameters);
    TEUCHOS_TEST_FOR_EXCEPTION(solver.solve() != Belos::Converged, std::runtime_error, "Inlet Laplace solve failed.");
    auto residual = Teuchos::rcp(new MultiVector<>(unique));
    matrix->apply(*solution, *residual);
    residual->update(-1., *rhs, 1.);
    Teuchos::Array<double> error(1), norm(1);
    residual->norm2(error()); rhs->norm2(norm());
    TEUCHOS_TEST_FOR_EXCEPTION(norm[0] == 0. || error[0] > 1.e-10*norm[0], std::runtime_error,
        "Inlet Laplace residual exceeds tolerance.");
    const double relativeResidual = error[0]/norm[0];
    const double peak = solution->getMax();
    TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(peak) || peak <= 0., std::runtime_error, "Invalid inlet profile.");
    solution->scale(1./peak);
    for (const auto value : solution->getData(0))
        TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(value) || value < -1.e-10 || value > 1.+1.e-10,
            std::runtime_error, "Inlet profile is outside [0,1].");
    HDF5Export<> exporter(unique, output);
    exporter.writeVariablesHDF5("solution", solution);
    exporter.closeExporter();
    // Check the actual import format and global node numbering before finishing.
    HDF5Import<> importer(unique, output);
    auto imported = importer.readVariablesHDF5("solution");
    auto difference = Teuchos::rcp(new MultiVector<>(unique));
    difference->update(1., *imported, -1., *solution, 0.);
    difference->norm2(error());
    TEUCHOS_TEST_FOR_EXCEPTION(error[0] != 0., std::runtime_error, "Inlet HDF5 round trip differs.");
    importer.closeImporter();
    if (comm->getRank() == 0)
        std::cout << "Inlet profile: " << feType << ", " << faces << " faces, " << unique->getGlobalNumElements()
                  << " global nodes, " << solver.getNumIters() << " iterations, relative residual "
                  << relativeResidual << ", exact HDF5 round trip, peak 1\n";
}
