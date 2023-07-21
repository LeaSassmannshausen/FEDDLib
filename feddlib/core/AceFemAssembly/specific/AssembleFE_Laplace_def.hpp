#ifndef ASSEMBLEFE_LAPLACE_DEF_hpp
#define ASSEMBLEFE_LAPLACE_DEF_hpp

#include "AssembleFE_Laplace_decl.hpp"

namespace FEDD {


/*!

 \brief Constructor for AssembleFE_Laplace

@param[in] flag Flag of element
@param[in] nodesRefConfig Nodes of element in reference configuration
@param[in] params Parameterlist for current problem
*/
template <class SC, class LO, class GO, class NO>
AssembleFE_Laplace<SC,LO,GO,NO>::AssembleFE_Laplace(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params,tuple)
{

}

/*!

 \brief Assembly Jacobian is simply assemblyLaplacian for Laplace Problem
 for scalar-valued function (dim=1) and also for vector-valued function (dim > 1)
@param[in] &elementMatrix
*/ 

template <class SC, class LO, class GO, class NO>
void AssembleFE_Laplace<SC,LO,GO,NO>::assembleJacobian() {

	int nodesElement = this->nodesRefConfig_.size();
	int dofs = std::get<2>(this->diskTuple_->at(0));
	int dofsElement = nodesElement*dofs;
	SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( dofsElement));

	assemblyLaplacian(elementMatrix);
	this->jacobian_ = elementMatrix ;
}

/*!
 \brief Assembly function for \f$ \int_T \nabla v \cdot \nabla u ~dx\f$ 
@param[in] &elementMatrix
*/
template <class SC, class LO, class GO, class NO>
void AssembleFE_Laplace<SC,LO,GO,NO>::assemblyLaplacian(SmallMatrixPtr_Type &elementMatrix) {

	int dim = this->getDim();
	int numNodes= std::get<3>(this->diskTuple_->at(0));//this->getNodesRefConfig().size();
	int Grad =2; // Needs to be fixed	
	string FEType = std::get<1>(this->diskTuple_->at(0));
	int dofs = std::get<2>(this->diskTuple_->at(0));

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    UN deg = Helper::determineDegree(dim,FEType,Grad);
    Helper::getDPhi(dPhi, weights, dim, FEType, deg);
    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    Helper::buildTransformation(B, this->nodesRefConfig_);
    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);

    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    Helper::applyBTinv( dPhi, dPhiTrans, Binv );
    for (UN i=0; i < numNodes; i++) {
        Teuchos::Array<SC> value( dPhiTrans[0].size(), 0. );
        for (UN j=0; j < numNodes; j++) {
            for (UN w=0; w<dPhiTrans.size(); w++) {
                for (UN d=0; d<dim; d++){
                    value[j] += weights->at(w) * dPhiTrans[w][i][d] * dPhiTrans[w][j][d];
                }
            }
            value[j] *= absDetB;
			
			for (UN d=0; d<dofs; d++) {
              (*elementMatrix)[i*dofs +d][j*dofs+d] = value[j];
            }
        }

    }

}

/*!

 \brief  

*/
template <class SC, class LO, class GO, class NO>
void AssembleFE_Laplace<SC,LO,GO,NO>::assembleRHS() {
	int dim = this->getDim();
	int Grad =1; // Needs to be fixed	
	int numNodes= std::get<3>(this->diskTuple_->at(0));//this->getNodesRefConfig().size();
	string FEType = std::get<1>(this->diskTuple_->at(0));
	vec_dbl_Type elementVector(numNodes);

    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
   
    UN deg = Helper::determineDegree(dim,FEType,Grad);
    Helper::getPhi(phi, weights, dim, FEType, deg);
    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    Helper::buildTransformation(B, this->nodesRefConfig_);
    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);

    std::vector<double> paras0(1);

    double x;
    SC value;

    //for now just const!
    std::vector<double> valueFunc(dim);
    SC* paras = &(paras0[0]);
    
    this->rhsFunc_( &x, &valueFunc[0], paras );

    for (UN i=0; i < phi->at(0).size(); i++) {
  	    value = Teuchos::ScalarTraits<SC>::zero();
		for (UN w=0; w<weights->size(); w++){
       		value += weights->at(w) * phi->at(w).at(i);
		}
        value *= absDetB *valueFunc[0];
        elementVector[i] += value;
    }

	(*this->rhsVec_) = elementVector;
}



}
#endif

