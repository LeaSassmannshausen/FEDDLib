#ifndef ASSEMBLEFESTOKES_DEF_hpp
#define ASSEMBLEFESTOKES_DEF_hpp

#include "AssembleFEStokes_decl.hpp"

namespace FEDD {

/*
Stokes is a version of Navier-Stokes neglecting the convective terms:

− \nabla \cdot ( \nabla \sigma) + \nabla p = f, 
\nabla \cdot u = 0, 

important is that the shear-stress related term with the Cauchy stress tensor \sigma =2 \mu  0.5(\nabla u + (\nabla u)^T )
can be written in two forms as viscosity \mu is for Newtonian fluids constants 

1. Laplace formulation: \nabla \cdot ( \nabla \sigma) =  \mu \Delta u

2. Stress-Divergence formulation with symmetrized velocity gradient (strain-rate tensor): \nabla \cdot ( \nabla \sigma) = \nabla \cdot (2 \mu  0.5(\nabla u + (\nabla u)^T ))

Information for the two ansatzes (see Pacheco et. al "A mixed finite element formulation for generalised Newtonian fluid flows with appropriate natural outflow boundary conditions"):
"The most popular treatments of the viscous stress are the Laplace formulation and the stress-divergence
form. The former is computationally cheaper and simpler [2], whereas the latter can be applied straightforwardly to more general (e.g., non-Newtonian) fluids [3]. Although completely equivalent in the Dirichlet case, they can yield different solutions when natural
boundary conditions are considered [4]. In the stress-divergence formulation, assuming zero natural outflow data is known for generating unphysical flow behaviour that can have
a considerable upstream effect, and for that reason the Laplace formulation is preferred when possible" 
*/


template <class SC, class LO, class GO, class NO>
AssembleFEStokes<SC,LO,GO,NO>::AssembleFEStokes(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params,tuple)
{
	int locVelocity=0;
	int locPressure=0;		
	if(std::get<0>(this->diskTuple_->at(0))=="Velocity"){
		locVelocity=0;
		locPressure=1;
	}
	else if(std::get<0>(this->diskTuple_->at(1))=="Velocity"){
		locVelocity=1;
		locPressure=0;
	}
	else
    	TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "No discretisation Information for Velocity in Navier Stokes Element." );
		

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	FETypeVelocity_ = std::get<1>(this->diskTuple_->at(locVelocity));
	FETypePressure_ =std::get<1>(this->diskTuple_->at(locPressure));

	dofsVelocity_ = std::get<2>(this->diskTuple_->at(locVelocity));
	dofsPressure_ =std::get<2>(this->diskTuple_->at(locPressure));

	numNodesVelocity_ = std::get<3>(this->diskTuple_->at(locVelocity));
	numNodesPressure_=std::get<3>(this->diskTuple_->at(locPressure));

	dofsElementVelocity_ = dofsVelocity_*numNodesVelocity_;
	dofsElementPressure_  = dofsPressure_*numNodesPressure_;	

	//this->solution_ = vec_dbl_Type(dofsElementVelocity_); //dofsElementPressure_+
	this->solutionVelocity_ = vec_dbl_Type(dofsElementVelocity_);
	this->solutionPressure_ = vec_dbl_Type(dofsElementPressure_);

 	viscosity_ = this->params_->sublist("Parameter").get("Viscosity",1.);
    density_ = this->params_->sublist("Parameter").get("Density",1.);

	dofsElement_ = dofsElementVelocity_+ dofsElementPressure_;

	SmallMatrix_Type coeff(2); // 2 \mu  0.5(\nabla u + (\nabla u)^T )
	coeff[0][0]=1.; coeff[0][1] = 1.; coeff[1][0] = 1.; coeff[1][1] = 1.; // we keep it constant like this for now. For BDF time disc. okay.
	coeff_ = coeff;

    linearization_ = this->params_->sublist("General").get("Linearization","FixedPoint"); // Information to assemble Jacobian accordingly
}

template <class SC, class LO, class GO, class NO>
void AssembleFEStokes<SC,LO,GO,NO>::setCoeff(SmallMatrix_Type coeff) {
	// We only substitute the coefficients if the matrix has the same 
	// size. In some non timedepenent cases the coeff matrix can be empty. 
	// We prevent that case.
	if(coeff.size() == 2)
		coeff_ = coeff;		

}

template <class SC, class LO, class GO, class NO>
void AssembleFEStokes<SC,LO,GO,NO>::assembleJacobian() {

	SmallMatrixPtr_Type elementMatrixN =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));
	SmallMatrixPtr_Type elementMatrixW =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));

    // As Stokes is linear we do not have to divide between first and following Newton-Steps because the matrices stay constant
	//if(this->newtonStep_ ==0){
		SmallMatrixPtr_Type elementMatrixA =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));
		SmallMatrixPtr_Type elementMatrixB =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));

		constantMatrix_.reset(new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));



        if ( this->params_->sublist("Parameter").get("Symmetric gradient",false) )
        {
        assemblyStress_Divergence(elementMatrixA); // 2 \mu  \nabla \cdot ( 0.5(\nabla u + (\nabla u)^T ))
        }
        else
		{
        assemblyLaplacian(elementMatrixA); // \mu \Delta u
        }

		elementMatrixA->scale(viscosity_); // KINEMATIC VISCOSITY THEREFORE WE HAVE TO MULTIPLY WITH DENSITY
		elementMatrixA->scale(density_);

		constantMatrix_->add( (*elementMatrixA),(*constantMatrix_));

		assemblyDivAndDivT(elementMatrixB); // For Matrix B
		elementMatrixB->scale(-1.);

		constantMatrix_->add( (*elementMatrixB),(*constantMatrix_));
    //}

	ANB_.reset(new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_)); // A + B + N
	ANB_->add( (*constantMatrix_),(*ANB_));

	this->jacobian_.reset(new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));
	this->jacobian_->add((*ANB_),(*this->jacobian_));

}



template <class SC, class LO, class GO, class NO>
void AssembleFEStokes<SC,LO,GO,NO>::assemblyLaplacian(SmallMatrixPtr_Type &elementMatrix) {

	int dim = this->getDim();
	int numNodes= numNodesVelocity_;
	int Grad =2; // Needs to be fixed	
	string FEType = FETypeVelocity_;
	int dofs = dofsVelocity_;

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
    Helper::applyBTinv( dPhi, dPhiTrans, Binv );( dPhi, dPhiTrans, Binv );
  	
    for (UN i=0; i < numNodes; i++) {
        Teuchos::Array<SC> value( dPhiTrans[0].size(), 0. );
        for (UN j=0; j < numNodes; j++) {
            for (UN w=0; w<dPhiTrans.size(); w++) {
                for (UN d=0; d<dim; d++){
                    value[j] += weights->at(w) * dPhiTrans[w][i][d] * dPhiTrans[w][j][d];
                }
            }
            value[j] *= absDetB;
			 /*if (std::fabs(value[j]) < pow(10,-14)) {
		            value[j] = 0.;
		        }*/
			for (UN d=0; d<dofs; d++) {
              (*elementMatrix)[i*dofs +d][j*dofs+d] = value[j];
            }
        }

    }
}



// Implementation from weak formulation of 2 \mu (1/2*(\nabla \phi_j + \nabla \phi_j^T)  :  \nabla \phi_i  )
// BUT viscosity is multiplied above because it is a constant  
template <class SC, class LO, class GO, class NO>
void AssembleFEStokes<SC,LO,GO,NO>::assemblyStress_Divergence(SmallMatrixPtr_Type &elementMatrix) {

	int dim = this->getDim();
	int numNodes= this->numNodesVelocity_;
	int Grad =2; // Needs to be fixed	
	string FEType = this->FETypeVelocity_;
	int dofs = this->dofsVelocity_; // For pressure it would be 1 

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    UN deg = Helper::determineDegree(dim,FEType,Grad); //  for P1 3
    Helper::getDPhi(dPhi, weights, dim, FEType, deg);  //  for deg 5 we get weight vector with 7 entries weights->at(7)
    // Example Values: dPhi->size() = 7 so number of quadrature points, dPhi->at(0).size() = 3 number of local element points, dPhi->at(0).at(0).size() = 2 as we have dim 2 therefore we have 2 derivatives (xi/eta in natural coordinates)
    // Phi is defined on reference element

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    Helper::buildTransformation(B, this->nodesRefConfig_);
    detB = B.computeInverse(Binv); // The function computeInverse returns a double value corrsponding to determinant of B     //B.print();
    absDetB = std::fabs(detB); // absolute value of B


    // dPhiTrans are the transorfmed basisfunctions, so B^(-T) * \grad_phi bzw. \grad_phi^T * B^(-1)
    // Corresponds to \hat{grad_phi}.
    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    Helper::applyBTinv( dPhi, dPhiTrans, Binv ); // so dPhiTrans corresponds now to our basisfunction in natural coordinates
    //dPhiTrans.size() = 7 so number of quadrature points, dPhiTrans[0].size() = 3 number local element points, dPhiTrans[0][0].size() = 2 as we have in dim=2 case two derivatives

    TEUCHOS_TEST_FOR_EXCEPTION(dim == 1,std::logic_error, "AssemblyStress Not implemented for dim=1");
    //***************************************************************************
    //***************************************************************************
    if (dim == 2)
    {
    
    //************************************
    // Compute entries    
    // Initialize some helper vectors/matrices
    double v11, v12, v21, v22, value1_j, value2_j , value1_i, value2_i, viscosity_atw;
       // Construct element matrices 
    for (UN i=0; i < numNodes; i++) 
    {
       // Teuchos::Array<SC> value(dPhiTrans[0].size(), 0. ); // dPhiTrans[0].size() is 3        
        for (UN j=0; j < numNodes; j++) 
        {
        // Reset values
        v11 = 0.0;v12 = 0.0;v21 = 0.0;v22 = 0.0;

            // So in general compute the components of eta*[ dPhiTrans_i : ( dPhiTrans_j + (dPhiTrans_j)^T )]
            for (UN w=0; w<dPhiTrans.size(); w++) 
            {

                 value1_j = dPhiTrans[w][j][0]; // so this corresponds to d\phi_j/dx
                 value2_j = dPhiTrans[w][j][1]; // so this corresponds to d\phi_j/dy

                 value1_i = dPhiTrans[w][i][0]; // so this corresponds to d\phi_i/dx
                 value2_i = dPhiTrans[w][i][1]; // so this corresponds to d\phi_i/dy

                 v11 = v11 +  weights->at(w) *(2.0*value1_i*value1_j+value2_i*value2_j);
                 v12 = v12 +  weights->at(w) *(value2_i*value1_j);
                 v21 = v21 +  weights->at(w) *(value1_i*value2_j);
                 v22 = v22 +  weights->at(w) *(2.0*value2_i*value2_j+value1_i*value1_j);
            
            } // loop end quadrature points
                //multiply determinant from transformation
            v11 *= absDetB; 
            v12 *= absDetB;
            v21 *= absDetB;
            v22 *= absDetB;
            
            // Put values on the right position in element matrix - d=2 because we are in two dimensional case
            // [v11  v12  ]
            // [v21  v22  ]
            (*elementMatrix)[i*dofs][j*dofs]   = v11; // d=0, first dimension
            (*elementMatrix)[i*dofs][j*dofs+1] = v12;
            (*elementMatrix)[i*dofs+1][j*dofs] = v21;
            (*elementMatrix)[i*dofs +1][j*dofs+1] =v22; //d=1, second dimension

        } // loop end over j node 

    } // loop end over i node


    } // end if dim 2

        else if (dim == 3)
    {
          //************************************#
         // Initialize some helper vectors/matrices
    double v11, v12, v13, v21, v22, v23, v31, v32, v33, value1_j, value2_j, value3_j , value1_i, value2_i, value3_i, viscosity_atw;


    // Construct element matrices 
     for (UN i=0; i < numNodes; i++) {
       // Teuchos::Array<SC> value(dPhiTrans[0].size(), 0. ); // dPhiTrans[0].size() is 3        
      
        for (UN j=0; j < numNodes; j++) {
        // Reset values
        v11 = 0.0;v12 = 0.0;v13=0.0; v21 = 0.0;v22 = 0.0;v23=0.0;v31=0.0;v32=0.0;v33=0.0;


            // So in general compute the components of eta*[ dPhiTrans_i : ( dPhiTrans_j + (dPhiTrans_j)^T )]
            for (UN w=0; w<dPhiTrans.size(); w++) {

                        value1_j = dPhiTrans.at(w).at(j).at(0); // so this corresponds to d\phi_j/dx
                        value2_j = dPhiTrans.at(w).at(j).at(1); // so this corresponds to d\phi_j/dy
                        value3_j = dPhiTrans.at(w).at(j).at(2); // so this corresponds to d\phi_j/dz


                        value1_i = dPhiTrans.at(w).at(i).at(0); // so this corresponds to d\phi_i/dx
                        value2_i = dPhiTrans.at(w).at(i).at(1); // so this corresponds to d\phi_i/dy
                        value3_i = dPhiTrans.at(w).at(i).at(2); // so this corresponds to d\phi_i/dz

                    
                       // Construct entries - we go over all quadrature points and if j is updated we set v11 etc. again to zero
                        v11 = v11 + weights->at(w)*(2.0*value1_j*value1_i+value2_j*value2_i+value3_j*value3_i);
                        v12 = v12 + weights->at(w)*(value2_i*value1_j);
                        v13 = v13 + weights->at(w)*(value3_i*value1_j);

                        v21 = v21 + weights->at(w)*(value1_i*value2_j);
                        v22=  v22 + weights->at(w)*(value1_i*value1_j+2.0*value2_j*value2_i+value3_j*value3_i);
                        v23 = v23 + weights->at(w)*(value3_i*value2_j);

                        v31 = v31 + weights->at(w)*(value1_i*value3_j);
                        v32 = v32 + weights->at(w)*(value2_i*value3_j);
                        v33 = v33 + weights->at(w)*(value1_i*value1_j+value2_i*value2_j+2.0*value3_i*value3_j);

                    }// loop end quadrature points

                         //multiply determinant from transformation
                    v11 *= absDetB ;
                    v12 *= absDetB ;
                    v13 *= absDetB ;
                    v21 *= absDetB ;
                    v22 *= absDetB ;
                    v23 *= absDetB ;
                    v31 *= absDetB ;
                    v32 *= absDetB ;
                    v33 *= absDetB ;

                   // Put values on the right position in element matrix 
                   // [v11  v12  v13]
                   // [v21  v22  v23]
                   // [v31  v32  v33]
            (*elementMatrix)[i*dofs][j*dofs]   = v11; // d=0, first dimension
            (*elementMatrix)[i*dofs][j*dofs+1] = v12;
            (*elementMatrix)[i*dofs][j*dofs+2] = v13;
            (*elementMatrix)[i*dofs+1][j*dofs] = v21;
            (*elementMatrix)[i*dofs +1][j*dofs+1] =v22; //d=1, second dimension
            (*elementMatrix)[i*dofs +1][j*dofs+2] =v23; //d=1, second dimension
            (*elementMatrix)[i*dofs+2][j*dofs] = v31;
            (*elementMatrix)[i*dofs +2][j*dofs+1] =v32; //d=2, third dimension
            (*elementMatrix)[i*dofs +2][j*dofs+2] =v33; //d=2, third dimension

                }// loop end over j node 
            }// loop end over i node 
        }// end if dim==3

}

template <class SC, class LO, class GO, class NO>
void AssembleFEStokes<SC,LO,GO,NO>::assemblyDivAndDivT(SmallMatrixPtr_Type &elementMatrix) {

    vec3D_dbl_ptr_Type 	dPhi;
    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
	int dim = this->dim_;

    UN deg =2; // Helper::determineDegree2( dim, FETypeVelocity_, FETypePressure_, Grad, Std);

    Helper::getDPhi(dPhi, weights, dim, FETypeVelocity_, deg);

    //if (FETypePressure_=="P1-disc-global")
      //  Helper::getPhiGlobal(phi, weights, dim, FETypePressure_, deg);
    if (FETypePressure_=="P1-disc" && FETypeVelocity_=="Q2" )
        Helper::getPhi(phi, weights, dim, FETypePressure_, deg, FETypeVelocity_);
    else
        Helper::getPhi(phi, weights, dim, FETypePressure_, deg);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);


	Helper::buildTransformation(B, this->nodesRefConfig_);
    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);

    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    Helper::applyBTinv( dPhi, dPhiTrans, Binv );

    Teuchos::Array<GO> rowIndex( 1, 0 );
    Teuchos::Array<SC> value(1, 0.);


    for (UN i=0; i < phi->at(0).size(); i++) {
        if (FETypePressure_=="P0")
            rowIndex[0] = GO ( 0 );
        else
            rowIndex[0] = GO ( i );

        for (UN j=0; j < dPhiTrans[0].size(); j++) {
            for (UN d=0; d<dim; d++){
                value[0] = 0.;
                for (UN w=0; w<dPhiTrans.size(); w++)
                    value[0] += weights->at(w) * phi->at(w)[i] * dPhiTrans[w][j][d];
                value[0] *= absDetB;


				(*elementMatrix)[rowIndex[0]+dofsVelocity_*numNodesVelocity_][dofsVelocity_ * j + d] +=value[0];
				(*elementMatrix)[dofsVelocity_ * j + d][dofsVelocity_*numNodesVelocity_+rowIndex[0]] +=value[0];
            }
        }
    }

}



// Assemble RHS with updated solution coming from Fixed Point Iter or der Newton.
template <class SC, class LO, class GO, class NO>
void AssembleFEStokes<SC,LO,GO,NO>::assembleRHS(){

	SmallMatrixPtr_Type elementMatrixN =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));

	ANB_.reset(new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_)); // A + B + N
	ANB_->add( (*constantMatrix_),(*ANB_));

	this->rhsVec_.reset( new vec_dbl_Type ( dofsElement_,0.) );
	// Multiplying ANB_ * solution // ANB Matrix without nonlinear part.
	int s=0,t=0;
	for(int i=0 ; i< ANB_->size();i++){
		if (i >= dofsElementVelocity_)
			s=1;
		for(int j=0; j < ANB_->size(); j++){
			if(j >= dofsElementVelocity_)
				t=1;
			(*this->rhsVec_)[i] += (*ANB_)[i][j]*(*this->solution_)[j]*coeff_[s][t];
		}
		t=0;
	}
}
}
#endif


