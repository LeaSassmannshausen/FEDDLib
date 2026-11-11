#ifndef FEDD_STRESS_RECOVERY_HPP
#define FEDD_STRESS_RECOVERY_HPP

#include "feddlib/core/General/GaussPointStress.hpp"
#include "feddlib/core/LinearAlgebra/BlockMultiVector.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/Mesh/Mesh.hpp"
#include <BelosTpetraAdapter.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <cmath>

namespace FEDD {

inline double gaussStressValue(const GaussPointStress& q, int field)
{
    return field == 0 ? q.vonMises : q.cauchy[field-1];
}

// Consistent, current-volume L2 projection: M_ij = integral(lambda_i lambda_j),
// b_i = integral(lambda_i q). Only corner nodes are unknowns. The remaining
// P2 nodes receive the linear interpolation afterwards, retaining P2 geometry.
template<class SC, class LO, class GO, class NO>
Teuchos::RCP<BlockMultiVector<SC,LO,GO,NO>> projectGaussPointStresses(
    Teuchos::RCP<Mesh<SC,LO,GO,NO>> mesh, const ElementGaussPointStress& data)
{
    using Vector = MultiVector<SC,LO,GO,NO>;
    using MatrixType = Matrix<SC,LO,GO,NO>;
    using MV = typename Vector::TpetraMultiVector_Type;
    using OP = Tpetra::Operator<SC,LO,GO,NO>;
    const auto unique = mesh->getMapUnique(), repeated = mesh->getMapRepeated();
    const auto elements = mesh->getElementsC();
    TEUCHOS_TEST_FOR_EXCEPTION(data.size() != elements->numberElements(), std::logic_error,
        "Gauss-point stresses do not match the solid elements.");
    auto cornerRep = Teuchos::rcp(new Vector(repeated, 1));
    cornerRep->putScalar(0.);
    for (UN T = 0; T < data.size(); ++T) {
        const auto nodes = elements->getElement(T).getVectorNodeList();
        TEUCHOS_TEST_FOR_EXCEPTION(nodes.size() != 10 || data[T].empty(), std::logic_error,
            "P1 stress recovery requires P2 tetrahedra with Gauss-point data.");
        for (int i = 0; i < 4; ++i) cornerRep->getDataNonConst(0)[nodes[i]] = 1.;
    }
    auto corners = Teuchos::rcp(new Vector(unique, 1));
    corners->putScalar(0.); corners->exportFromVector(cornerRep, false, "Add");
    auto mass = Teuchos::rcp(new MatrixType(unique, 64));
    auto rhsRep = Teuchos::rcp(new Vector(repeated, 10));
    rhsRep->putScalar(0.);
    for (UN T = 0; T < data.size(); ++T) {
        const auto nodes = elements->getElement(T).getVectorNodeList();
        double localMass[4][4]{};
        for (const auto& q : data[T]) {
            TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(q.weight), std::runtime_error,
                "Non-finite weight in P1 stress recovery.");
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j)
                    localMass[i][j] += q.weight*q.shapeP1[i]*q.shapeP1[j];
                for (int field = 0; field < 10; ++field)
                    rhsRep->getDataNonConst(field)[nodes[i]] += q.weight*q.shapeP1[i]*gaussStressValue(q, field);
            }
        }
        Teuchos::Array<GO> columns(4);
        for (int i = 0; i < 4; ++i) columns[i] = repeated->getGlobalElement(nodes[i]);
        for (int i = 0; i < 4; ++i) {
            Teuchos::Array<SC> values(4);
            for (int j = 0; j < 4; ++j) values[j] = localMass[i][j];
            mass->insertGlobalValues(columns[i], columns(), values());
        }
    }
    for (LO i = 0; i < unique->getNodeNumElements(); ++i)
        if (corners->getData(0)[i] == 0.) {
            Teuchos::Array<GO> column(1, unique->getGlobalElement(i));
            Teuchos::Array<SC> one(1, 1.);
            mass->insertGlobalValues(column[0], column(), one());
        }
    mass->fillComplete();
    auto rhs = Teuchos::rcp(new Vector(unique, 10));
    rhs->putScalar(0.); rhs->exportFromVector(rhsRep, false, "Add");
    auto solution = Teuchos::rcp(new Vector(unique, 10));
    solution->putScalar(0.);
    auto problem = Teuchos::rcp(new Belos::LinearProblem<SC,MV,OP>(mass->getTpetraMatrix(),
        solution->getTpetraMultiVectorNonConst(), rhs->getTpetraMultiVector()));
    auto diagonal = mass->buildDiagonalInverse("Diagonal");
    problem->setLeftPrec(diagonal->getTpetraMatrix());
    TEUCHOS_TEST_FOR_EXCEPTION(!problem->setProblem(), std::runtime_error,
        "Cannot initialize the P1 stress projection.");
    auto controls = Teuchos::rcp(new Teuchos::ParameterList);
    controls->set("Convergence Tolerance", 1.e-11).set("Maximum Iterations", 1000);
    Belos::PseudoBlockCGSolMgr<SC,MV,OP> solver(problem, controls);
    TEUCHOS_TEST_FOR_EXCEPTION(solver.solve() != Belos::Converged, std::runtime_error,
        "The P1 stress projection did not converge.");
    auto residual = Teuchos::rcp(new Vector(unique, 10));
    mass->apply(*solution, *residual); residual->update(-1., *rhs, 1.);
    Teuchos::Array<typename Teuchos::ScalarTraits<SC>::magnitudeType> errors(10), norms(10);
    residual->norm2(errors()); rhs->norm2(norms());
    for (int field = 0; field < 10; ++field)
        TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(errors[field]) || errors[field] > 1.e-9*norms[field]+1.e-25,
            std::runtime_error, "P1 stress projection residual exceeds tolerance.");

    // All incident elements give the same midside value. Sum and divide on
    // the unique map so this also works when an edge owner has no local cell.
    auto valuesRep = Teuchos::rcp(new Vector(repeated, 10));
    valuesRep->importFromVector(solution, false, "Insert");
    auto interpolatedRep = Teuchos::rcp(new Vector(repeated, 10));
    auto countRep = Teuchos::rcp(new Vector(repeated, 1));
    interpolatedRep->putScalar(0.); countRep->putScalar(0.);
    const int edges[6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
    for (UN T = 0; T < data.size(); ++T) {
        const auto nodes = elements->getElement(T).getVectorNodeList();
        for (int e = 0; e < 6; ++e) {
            countRep->getDataNonConst(0)[nodes[4+e]] += 1.;
            for (int field = 0; field < 10; ++field)
                interpolatedRep->getDataNonConst(field)[nodes[4+e]] += .5*(
                    valuesRep->getData(field)[nodes[edges[e][0]]]+valuesRep->getData(field)[nodes[edges[e][1]]]);
        }
    }
    auto interpolated = Teuchos::rcp(new Vector(unique, 10));
    auto count = Teuchos::rcp(new Vector(unique, 1));
    interpolated->putScalar(0.); count->putScalar(0.);
    interpolated->exportFromVector(interpolatedRep, false, "Add");
    count->exportFromVector(countRep, false, "Add");
    auto result = Teuchos::rcp(new BlockMultiVector<SC,LO,GO,NO>(10));
    for (int field = 0; field < 10; ++field) {
        auto value = Teuchos::rcp(new Vector(unique, 1));
        for (LO i = 0; i < unique->getNodeNumElements(); ++i)
            value->getDataNonConst(0)[i] = count->getData(0)[i] != 0.
                ? interpolated->getData(field)[i]/count->getData(0)[i] : solution->getData(field)[i];
        result->addBlock(value, field);
    }
    return result;
}

}
#endif
