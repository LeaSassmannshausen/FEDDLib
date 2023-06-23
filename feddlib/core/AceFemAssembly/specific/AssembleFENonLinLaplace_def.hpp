#ifndef ASSEMBLEFENONLINLAPLACE_DEF_hpp
#define ASSEMBLEFENONLINLAPLACE_DEF_hpp

#include "AssembleFENonLinLaplace_decl.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include <Teuchos_Array.hpp>

namespace FEDD {
// TODO extract construction of (1+u^2)\nabla u to a separate function

/*!

 \brief Constructor for AssembleFENonLinLaplace

@param[in] flag Flag of element
@param[in] nodesRefConfig Nodes of element in reference configuration
@param[in] params Parameterlist for current problem

*/
template <class SC, class LO, class GO, class NO>
AssembleFENonLinLaplace<SC, LO, GO, NO>::AssembleFENonLinLaplace(
    int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,
    tuple_disk_vec_ptr_Type tuple)
    : AssembleFE<SC, LO, GO, NO>(flag, nodesRefConfig, params, tuple) {}

/*!

 \brief Assemble Jacobian as per Gateaux derivative of the laplacian
@param[in] &elementMatrix

*/

template <class SC, class LO, class GO, class NO>
void AssembleFENonLinLaplace<SC, LO, GO, NO>::assembleJacobian() {

  int nodesElement = this->nodesRefConfig_.size();
  int dofs = std::get<2>(this->diskTuple_->at(0));
  int dofsElement = nodesElement * dofs;
  SmallMatrixPtr_Type elementMatrix =
      Teuchos::rcp(new SmallMatrix_Type(dofsElement));
  assemblyNonLinLaplacian(elementMatrix);

  this->jacobian_ = elementMatrix;
}

/*!

 \brief Assembly function for \f$ \int_T \nabla v \cdot \nabla u ~dx\f$

@param[in] &elementMatrix

*/
template <class SC, class LO, class GO, class NO>
void AssembleFENonLinLaplace<SC, LO, GO, NO>::assemblyNonLinLaplacian(
    SmallMatrixPtr_Type &elementMatrix) {

  int dim = this->getDim();
  int Grad = 2; // Needs to be fixed
  string FEType = std::get<1>(this->diskTuple_->at(0));
  int dofs = std::get<2>(this->diskTuple_->at(0));
  // Same as this->getNodesRefConfig().size();
  int numNodes = std::get<3>(this->diskTuple_->at(0));
  int dofsElement = numNodes * dofs;
  UN deg = Helper::determineDegree(dim, FEType, Grad);

  vec3D_dbl_ptr_Type dPhi;
  vec2D_dbl_ptr_Type phi;
  vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
  vec_dbl_Type uLoc(weights->size(), -1.);
  // At each quad node p_i gradient vector is duLoc[w] = [\partial_x phi,
  // \partial_y phi]
  vec2D_dbl_Type duLoc(weights->size(), vec_dbl_Type(dim, 0));

  Helper::getDPhi(dPhi, weights, dim, FEType, deg);
  Helper::getPhi(phi, weights, dim, FEType, deg);

  SC detB;
  SC absDetB;
  SmallMatrix<SC> B(dim);
  SmallMatrix<SC> Binv(dim);
  buildTransformation(B);
  detB = B.computeInverse(Binv);
  absDetB = std::fabs(detB);
  vec3D_dbl_Type dPhiTrans(
      dPhi->size(), vec2D_dbl_Type(dPhi->at(0).size(), vec_dbl_Type(dim, 0.)));
  applyBTinv(dPhi, dPhiTrans, Binv);

  // Build vector of current solution at quadrature nodes
  for (int w = 0; w < phi->size(); w++) { // quadrature nodes
    uLoc[w] = 0.;
    for (int i = 0; i < phi->at(0).size();
         i++) { // each basis function at the current quadrature node
      uLoc[w] += (*this->solution_)[i] * phi->at(w).at(i);
      for (int d = 0; d < dim; d++) {
        duLoc[w][d] += (*this->solution_)[i] * dPhiTrans[w][i][d];
      }
    }
  }

  // Build local stiffness matrix
  Teuchos::Array<SC> value(1, 0.);
  for (UN i = 0; i < numNodes; i++) {
    for (UN j = 0; j < numNodes; j++) {
      value[0] = 0.;
      for (UN w = 0; w < weights->size(); w++) {
        for (UN d = 0; d < dim; d++) {
          value[0] += weights->at(w) *
                      ((1 + uLoc[w] * uLoc[w]) * dPhiTrans[w][i][d] +
                       2 * uLoc[w] * phi->at(w).at(i) * duLoc[w][d]) *
                      dPhiTrans[w][j][d];
        }
      }
      value[0] *= absDetB;
      (*elementMatrix)[i][j] = value[0];
    }
  }
}

/*!

 \brief Assembly function for \f$ \int_T f ~ v ~dx \f$, we need to

@param[in] &elementVector

*/
template <class SC, class LO, class GO, class NO>
void AssembleFENonLinLaplace<SC, LO, GO, NO>::assembleRHS() {

  int dim = this->getDim();
  int Grad = 2; // Needs to be fixed
  string FEType = std::get<1>(this->diskTuple_->at(0));
  int dofs = std::get<2>(this->diskTuple_->at(0));
  // Same as this->getNodesRefConfig().size();
  int numNodes = std::get<3>(this->diskTuple_->at(0));
  int dofsElement = numNodes * dofs;
  UN deg = Helper::determineDegree(dim, FEType, Grad);

  vec3D_dbl_ptr_Type dPhi;
  vec2D_dbl_ptr_Type phi;
  vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
  vec_dbl_Type uLoc(weights->size(), -1.);

  Helper::getDPhi(dPhi, weights, dim, FEType, deg);
  Helper::getPhi(phi, weights, dim, FEType, deg);

  SC detB;
  SC absDetB;
  SmallMatrix<SC> B(dim);
  SmallMatrix<SC> Binv(dim);
  buildTransformation(B);
  detB = B.computeInverse(Binv);
  absDetB = std::fabs(detB);
  vec3D_dbl_Type dPhiTrans(
      dPhi->size(), vec2D_dbl_Type(dPhi->at(0).size(), vec_dbl_Type(dim, 0.)));
  applyBTinv(dPhi, dPhiTrans, Binv);

  // for now just const!
  double x;
  std::vector<double> paras0(1);
  std::vector<double> valueFunc(dim);
  SC *paras = &(paras0[0]);
  this->rhsFunc_(&x, &valueFunc[0], paras);

  // Build vector of current solution at quadrature nodes
  for (int w = 0; w < phi->size(); w++) { // quadrature nodes
    uLoc[w] = 0.;
    for (int i = 0; i < phi->at(0).size();
         i++) { // each basis function at the current quadrature node
      uLoc[w] += (*this->solution_)[i] * phi->at(w).at(i);
    }
  }

  // Build local stiffness matrx
  for (UN i = 0; i < numNodes; i++) { // Iterate basis functions (columns)
    // $fv$ term
    Teuchos::Array<SC> valueOne(1, 0.);
    // $(1 + u_0^2)\nabla u_0 \nabla v$ term
    Teuchos::Array<SC> valueTwo(1, 0.);
    // for storing matrix_column \cdot current_solution reduction operation
    Teuchos::Array<SC> reductionValue(1, 0.);
    for (UN j = 0; j < numNodes; j++) { // Iterate test functions (rows)
      valueOne[0] = 0.;
      valueTwo[0] = 0.;
      for (UN w = 0; w < dPhiTrans.size(); w++) { // Iterate quadrature nodes
        for (UN d = 0; d < dim; d++) {
          valueOne[0] += weights->at(w) *
                         phi->at(w).at(i); // Note that phi does not require
                                           // transformation. absDetB suffices.
          valueTwo[0] += weights->at(w) * (1 + uLoc[w] * uLoc[w]) *
                         dPhiTrans[w][i][d] * dPhiTrans[w][j][d];
        }
      }
      valueOne[0] *= absDetB * valueFunc[0];
      valueTwo[0] *= absDetB;
      // Perform reduction of current valueTwo and solution
      reductionValue[0] = valueTwo[0] * (*this->solution_)[j];
    }
    (*this->rhsVec_)[i] = valueOne[0] - reductionValue[0];
  }
}

/*!

 \brief Building Transformation

@param[in] &B

*/

template <class SC, class LO, class GO, class NO>
void AssembleFENonLinLaplace<SC, LO, GO, NO>::buildTransformation(
    SmallMatrix<SC> &B) {

  TEUCHOS_TEST_FOR_EXCEPTION((B.size() < 2 || B.size() > 3), std::logic_error,
                             "Initialize SmallMatrix for transformation.");
  UN index;
  UN index0 = 0;
  for (UN j = 0; j < B.size(); j++) {
    index = j + 1;
    for (UN i = 0; i < B.size(); i++) {
      B[i][j] = this->nodesRefConfig_.at(index).at(i) -
                this->nodesRefConfig_.at(index0).at(i);
    }
  }
}

/*!

 \brief Assembly function for \f$ \int_T f ~ v ~dx \f$, we need to

@param[in] &elementVector

*/

template <class SC, class LO, class GO, class NO>
void AssembleFENonLinLaplace<SC, LO, GO, NO>::applyBTinv(
    vec3D_dbl_ptr_Type &dPhiIn, vec3D_dbl_Type &dPhiOut,
    SmallMatrix<SC> &Binv) {
  UN dim = Binv.size();
  for (UN w = 0; w < dPhiIn->size(); w++) {
    for (UN i = 0; i < dPhiIn->at(w).size(); i++) {
      for (UN d1 = 0; d1 < dim; d1++) {
        for (UN d2 = 0; d2 < dim; d2++) {
          dPhiOut[w][i][d1] += dPhiIn->at(w).at(i).at(d2) * Binv[d2][d1];
        }
      }
    }
  }
}

} // namespace FEDD
#endif
