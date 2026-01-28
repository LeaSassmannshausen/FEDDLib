#ifndef ASSEMBLEFE_NONLINELAS_DEF_hpp
#define ASSEMBLEFE_NONLINELAS_DEF_hpp

#include "AssembleFE_NonLinElas_decl.hpp"
#include "feddlib/core/AceFemAssembly/AceInterface/NeoHookQuadraticTets2.hpp"
#include <vector>
#include <iostream>

namespace FEDD {


/*!

 \brief Constructor for AssembleFE_NonLinElas

@param[in] flag Flag of element
@param[in] nodesRefConfig Nodes of element in reference configuration
@param[in] params Parameterlist for current problem
@param[in} tuple Vector of tuples with Discretization information

*/
template <class SC, class LO, class GO, class NO>
AssembleFE_NonLinElas<SC,LO,GO,NO>::AssembleFE_NonLinElas(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params,tuple),
v_(1060),
d_(2),
ul_(30),
ul0_(30),
xl_(30),
s_(900),
p_(30),
ht_(10),
hp_(10)
{
	/// Extracting values from ParameterList params:
	E_ = this->params_->sublist("Parameter").get("E",1000.0); // the last value is the dafault value, in case no parameter is set
    //lambda_ = this->params_->sublist("Parameter").get("lambda",1.);
    poissonRatio_ = this->params_->sublist("Parameter").get("Poisson Ratio",0.4e-0);

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	FEType_ = std::get<1>(this->diskTuple_->at(0)); // FEType of Disk
	dofs_ = std::get<2>(this->diskTuple_->at(0)); // Degrees of freedom per node
	numNodes_ = std::get<3>(this->diskTuple_->at(0)); // Number of nodes of element

	dofsElement_ = dofs_*numNodes_; // "Dimension of return matrix"

	this->rhsVec_.reset( new vec_dbl_Type ( dofsElement_,0.) );

	this->elementMatrix_.reset(new SmallMatrix_Type( dofsElement_)); // Matrix we fill with entries.

	this->solution_.reset( new vec_dbl_Type (dofsElement_,0.) );
}

/*!

 \brief Assembly Jacobian

@param[in] &elementMatrix

*/ 

template <class SC, class LO, class GO, class NO>
void AssembleFE_NonLinElas<SC,LO,GO,NO>::assembleJacobian() {

	assemblyNonLinElas(elementMatrix_); // Function that fills the matrix. We pass though a pointer that will be filled.

	this->jacobian_ = elementMatrix_ ; // We init the jacobian matrix with the matrix we just build.
}

/*!

 \brief Assembly function 

@param[in] &elementMatrix

*/
template <class SC, class LO, class GO, class NO>
void AssembleFE_NonLinElas<SC,LO,GO,NO>::assemblyNonLinElas(SmallMatrixPtr_Type &elementMatrix) {

	/// We can access the following values we initialized/extracted in the constructor:
	// dofs_
	// FEType_
	// numNodes_
	// dofsElement_
	// mu_
	// poissonRatio_

	/// Writing entries in the element matrix for nodes 1,2..n , n=numNodes_
	/// 1_x 1_y 1_z 2_x 2_y 2_z .... n_x n_y n_z 
	
	// std::vector<double> v(1060); //Working vector, size defined by AceGen-FEAP
	// std::vector<double> d(2); // Material parameters
	// std::vector<double> ul(30); // The solution vector(or displacement in this case)
	// std::vector<double> ul0(30); // Currently unused but must be passed to match FEAP template
	// std::vector<double> xl(30); // Nodal Positions in reference coordinates
	// std::vector<double> s(900); // Element Stiffness Matrix [Output from skr]
	// std::vector<double> p(30); // Residual vector [Output from skr]
	// std::vector<double> ht(10); // History parameters currently unused
	// std::vector<double> hp(10); // History parameters currently unused

	std::fill(v_.begin(), v_.end(), 0.0);
	std::fill(s_.begin(), s_.end(), 0.0);
	std::fill(p_.begin(), p_.end(), 0.0);
	std::fill(ht_.begin(), ht_.end(), 0.0);
	std::fill(hp_.begin(), hp_.end(), 0.0);

	d_[0] = this->E_; // TODO: Check order if there is a problem
	d_[1] = this->poissonRatio_;

	// for(int i=0;i<30;i++)
	// 	ul_[i] = (*this->solution_)[i]; // What is the order? I need it in the form (u1,v1,w1,u2,v2,w2,...)

    std::copy_n(this->solution_->begin(), 30, ul_.begin());

	// int count = 0;
	// for(int i=0;i<this->numNodes_;i++)
	// 	for(int j=0;j<this->dofs_;j++){
	// 		xl_[count] = this->getNodesRefConfig()[i][j];
	// 		count++;}	
	const auto& nodesRef = this->getNodesRefConfig();
	auto it = xl_.begin();
	for(int i = 0; i < this->numNodes_; i++) {
		it = std::copy_n(nodesRef[i].begin(), this->dofs_, it);
	}

	// std::cout << "[DEBUG] SKR-Jacobian Calls after this line!" << std::endl;
	if(!this->isComputed_){
		skr2(v_.data(), d_.data(), ul_.data(), ul0_.data(), xl_.data(), s_.data(), p_.data(), ht_.data(), hp_.data());
		this->isComputed_ = true;
	}
	// std::cout << "[DEBUG] SKR-Jacobian Call successful!" << std::endl;
	// Note: FEAP/Fortran returns matrices unrolled in column major form. This must be converted for use here.

	/* std::cout << "[DEBUG] Printing FEAP output Stiffness vector" << std::endl;
	for (int i=0;i<900;i++)
		std::cout << s[i] << " ";*/

    for (UN i=0; i < this->dofsElement_; i++) {
        for (UN j=0; j < this->dofsElement_; j++) {
            (*elementMatrix)[i][j] = -s_[this->dofsElement_*j+i]; // Rolling into a matrix using column major (m*j+i)
        }
    }

}

/*!

 \brief Assembly RHS


*/
template <class SC, class LO, class GO, class NO>
void AssembleFE_NonLinElas<SC,LO,GO,NO>::assembleRHS() {

	// [Efficiency] Need to know which is called first: assembleRHS() or assembleJacobian(), so that multiple calls to skr() may be avoided.
	//this->rhsVec_.reset( new vec_dbl_Type ( dofsElement_,0.) );

	// Note skr() computes both elementMatrix_ and rhsVec_
	// std::vector<double> v(1060); //Working vector, size defined by AceGen-FEAP
	// std::vector<double> d(2); // Material parameters
	// std::vector<double> ul(30); // The solution vector(or displacement in this case)
	// std::vector<double> ul0(30); // Currently unused but must be passed to match FEAP template
	// std::vector<double> xl(30); // Nodal Positions in reference coordinates
	// std::vector<double> s(900); // Element Stiffness Matrix [Output from skr]
	// std::vector<double> p(30); // Residual vector [Output from skr]
	// std::vector<double> ht(10); // History parameters currently unused
	// std::vector<double> hp(10); // History parameters currently unused

	std::fill(v_.begin(), v_.end(), 0.0);
	std::fill(s_.begin(), s_.end(), 0.0);
	std::fill(p_.begin(), p_.end(), 0.0);
	std::fill(ht_.begin(), ht_.end(), 0.0);
	std::fill(hp_.begin(), hp_.end(), 0.0);

	d_[0] = this->E_; // TODO: Check order if there is a problem
	d_[1] = this->poissonRatio_;

	// for(int i=0;i<30;i++){
	// 	ul_[i] = (*this->solution_)[i];
	// }
    std::copy_n(this->solution_->begin(), 30, ul_.begin());


	// int count = 0;
	// for(int i=0;i<this->numNodes_;i++)
	// 	for(int j=0;j<this->dofs_;j++){
	// 		xl_[count] = this->getNodesRefConfig()[i][j];
	// 		count++;}

	const auto& nodesRef = this->getNodesRefConfig();
	auto it = xl_.begin();
	for(int i = 0; i < this->numNodes_; i++) {
		it = std::copy_n(nodesRef[i].begin(), this->dofs_, it);
	}

	if(!this->isComputed_){
		skr2(v_.data(), d_.data(), ul_.data(), ul0_.data(), xl_.data(), s_.data(), p_.data(), ht_.data(), hp_.data());
		this->isComputed_ = true;
	}

	for(int i=0; i< p_.size(); i++)
		(*this->rhsVec_)[i] = -p_[i];
}
template <class SC, class LO, class GO, class NO>
void AssembleFE_NonLinElas<SC,LO,GO,NO>:: updateParameter(std::string type, double value){
	if(type == "E")
		this->E_ = value;


}


}
#endif

