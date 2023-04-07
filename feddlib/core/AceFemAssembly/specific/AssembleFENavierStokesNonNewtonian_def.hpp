#ifndef ASSEMBLEFENAVIERSTOKESNONNEWTONIAN_DEF_hpp
#define ASSEMBLEFENAVIERSTOKESNONNEWTONIAN_DEF_hpp

#include "AssembleFENavierStokes_decl.hpp"

namespace FEDD {
// All important things are so far defined in AssembleFENavierStokes. Please check there.
/*
This class is specifically for Non-Newtonian fluids where we consider generalized Newtonian models for the viscosity
- Because the viscosity is no longer constant the conventional formulation with the laplacian term can not be considered
(but there is a generalized laplacian version of the equation see "On outflow boundary conditions in finite element simulations of non- Newtonian internal flow" 2021)
Instead we use the stress-divergence formulation of the momentum equation and derive from that the element-wise entries
*/
template <class SC, class LO, class GO, class NO>
AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::AssembleFENavierStokesNonNewtonian(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFENavierStokes<SC,LO,GO,NO>(flag, nodesRefConfig, params,tuple)
{

    ////******************* If we want to save viscosity - it is also possible to compute in paraview**********************************
	dofsElementViscosity_ = this->dofsPressure_*this->numNodesVelocity_; // So it is a scalar quantity but as it depend on the velocity it is defined at the nodes of the velocity
	this->solutionViscosity_ = vec_dbl_Type(dofsElementViscosity_ );    ////**********************************************************************************
    
    // Reading through parameterlist
    shearThinningModel= params->sublist("Material").get("ShearThinningModel","");
    // New: We have to check which material model we use 
	if(shearThinningModel == "Carreau-Yasuda"){
		Teuchos::RCP<CarreauYasuda<SC,LO,GO,NO>> materialModelSpecific(new CarreauYasuda<SC,LO,GO,NO>(params) );
		materialModel = materialModelSpecific;
	}
	else if(shearThinningModel == "Power-Law"){
		Teuchos::RCP<PowerLaw<SC,LO,GO,NO>> materialModelSpecific(new PowerLaw<SC,LO,GO,NO>(params) );
		materialModel = materialModelSpecific;
	}
    else if(shearThinningModel == "Dimless-Carreau"){
		Teuchos::RCP<Dimless_Carreau<SC,LO,GO,NO>> materialModelSpecific(new Dimless_Carreau<SC,LO,GO,NO>(params) );
		materialModel = materialModelSpecific;
	}
	else
    		TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No specific implementation for your material model request. Valid are:Carreau-Yasuda, Power-Law, Dimless-Carreau");
}






template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::assembleJacobian() 
{

    // For nonlinear generalized newtonian stress tensor part
	SmallMatrixPtr_Type elementMatrixN =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
	SmallMatrixPtr_Type elementMatrixW =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));

    // For nonlinear convection
    SmallMatrixPtr_Type elementMatrixNC =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
	SmallMatrixPtr_Type elementMatrixWC =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
    

    // In the first iteration step we initialize the constant matrices
    // So in the case of a newtonian fluid we would have the matrix A with the contributions of the Laplacian term
    // and the matrix B with the mixed-pressure terms. Latter one exists also in the non-newtonian case
	if(this->newtonStep_ ==0){
		SmallMatrixPtr_Type elementMatrixA =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
		SmallMatrixPtr_Type elementMatrixB =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));

		this->constantMatrix_.reset(new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
        // Construct the matrix B from FE formulation - as it is equivalent to Newtonian case we call the same function
	    this->assemblyDivAndDivT(elementMatrixB); // For Matrix B
		elementMatrixB->scale(-1.);
		this->constantMatrix_->add( (*elementMatrixB),(*this->constantMatrix_));

    }
    
    // The other element matrices are not constant so we have to update them in each step
    // As the stress tensor term, considering a non-newtonian constitutive equation, is nonlinear we add its contribution here
  
    // ANB is the FixedPoint Formulation which was named for newtonian fluids. 
    // Matrix A (Laplacian term (here not occuring)), Matrix B for div-Pressure Part, Matrix Nx for nonlinear parts -

	this->ANB_.reset(new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_)); // A + B + N
	this->ANB_->add( ((*this->constantMatrix_)),((*this->ANB_)));

    // Nonlinear advection term \rho (u \cdot \nabla) u 
    // As this class is derived from NavierStokes class we can call already implemented function
    //*************** ADVECTION************************
	/*this->assemblyAdvection(elementMatrixNC);
	elementMatrixNC->scale(this->density_);
	this->ANB_->add( (*elementMatrixNC),((*this->ANB_)));
    */
    // For a non-newtonian fluid we add additional element matrix and fill it with specific contribution
    // Remember that this term is based on the stress-divergence formulation of the momentum equation
    // \nabla \dot \tau with \tau=\eta(\gammaDot)(\nabla u + (\nabla u)^T)
    //*************** STRESS TENSOR************************
    this->assemblyStress(elementMatrixN);
	this->ANB_->add( (*elementMatrixN),((*this->ANB_)));


    this->jacobian_.reset(new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
	this->jacobian_->add((*this->ANB_),(*this->jacobian_));

    // If linearization is not FixdPoint (so NOX or Newton) we add the derivative to the Jacobian matrix. Otherwise the FixedPoint formulation becomes the jacobian.
    if(this->linearization_ != "FixedPoint"){

	    this->assemblyStressDev(elementMatrixW);     // stress tensor
      //  this->assemblyAdvectionInU(elementMatrixWC); // convection
	  //  elementMatrixWC->scale(this->density_);
        
        this->jacobian_->add((*elementMatrixW),(*this->jacobian_)); 
      //  this->jacobian_->add((*elementMatrixWC),(*this->jacobian_));  // int add(SmallMatrix<T> &bMat, SmallMatrix<T> &cMat); //this+B=C elementMatrix + constantMatrix_;
    }
	
    //**************** BOUNDARY TERM *******************************+
    // Because we have stress-divergence form of Navier-Stokes equations in the non-newtonian case
    // we have to add a extra boundary term to get the same outflow boundary condition as in the conventional formulation with
    // the laplacian operator in the equations due to the fact that in the stress-divergence formulation the
    // natural boundary condition is different 
    // We have to check whether it is an element which has edges (2D) / surfaces (3D) corresponding to an Outflow Neumann boundary
    // Then we have to compute contribution 
    if (this->surfaceElement == true) // Only if we consider element with neumann edge
    {
      SmallMatrixPtr_Type elementMatrixNB =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
      this->assemblyNeumannBoundaryTerm(elementMatrixNB);
      this->ANB_->add( (*elementMatrixNB),((*this->ANB_)));

      // If linearization is not FixdPoint (so NOX or Newton) we add the derivative to the Jacobian matrix. Otherwise the FixedPoint formulation becomes the jacobian.
      /*if(this->linearization_ != "FixedPoint")
      {
        SmallMatrixPtr_Type elementMatrixNBW =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
	    this->assemblyNeumannBoundaryTermDev(elementMatrixNBW); //
        this->jacobian_->add((*elementMatrixNBW),(*this->jacobian_));  
      }*/

    }

}

// Extra stress term resulting from chosen non-newtonian constitutive model  - Compute element matrix entries
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::assemblyStress(SmallMatrixPtr_Type &elementMatrix) {

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
  
    buildTransformation(B);
    detB = B.computeInverse(Binv); // The function computeInverse returns a double value corrsponding to determinant of B     //B.print();
    absDetB = std::fabs(detB); // absolute value of B


    // dPhiTrans are the transorfmed basisfunctions, so B^(-T) * \grad_phi bzw. \grad_phi^T * B^(-1)
    // Corresponds to \hat{grad_phi}.
    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    applyBTinv( dPhi, dPhiTrans, Binv ); // so dPhiTrans corresponds now to our basisfunction in natural coordinates
    //dPhiTrans.size() = 7 so number of quadrature points, dPhiTrans[0].size() = 3 number local element points, dPhiTrans[0][0].size() = 2 as we have in dim=2 case two derivatives

    TEUCHOS_TEST_FOR_EXCEPTION(dim == 1,std::logic_error, "AssemblyStress Not implemented for dim=1");
    //***************************************************************************
    //***************************************************************************
    if (dim == 2)
    {
    //************************************
    // Compute shear rate gammaDot, which is a vector because it is evaluated at a gaussian quadrature point 
    // for that compute velocity gradient
    vec_dbl_ptr_Type gammaDot(new vec_dbl_Type(weights->size(),0.0)); //gammaDot->at(j) j=0...weights 
    computeShearRate( dPhiTrans, gammaDot, dim); // updates gammaDot using velocity solution 
    //************************************
    // Compute entries    
    // Initialize some helper vectors/matrices
    double v11, v12, v21, v22, value1_j, value2_j , value1_i, value2_i, viscosity_atw;
    viscosity_atw = 0.;
    
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


                // viscosity function evaluated where we consider the dynamic viscosity!!
                this->materialModel->evaluateFunction(this->params_,  gammaDot->at(w), viscosity_atw);
                //if (i==0 && j==0) viscosity_averageOverT += viscosity_atw; // see below ### 

                 v11 = v11 + viscosity_atw * weights->at(w) *(2.0*value1_i*value1_j+value2_i*value2_j);
                 v12 = v12 + viscosity_atw * weights->at(w) *(value2_i*value1_j);
                 v21 = v21 + viscosity_atw * weights->at(w) *(value1_i*value2_j);
                 v22 = v22 + viscosity_atw * weights->at(w) *(2.0*value2_i*value2_j+value1_i*value1_j);
            
            } // loop end quadrature points
            
            //multiply determinant from transformation
           // if (i==0 && j==0)  viscosity_averageOverT *= absDetB; 
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
    //***************************************************************************
    else if (dim == 3)
    {
          //************************************#

    // Compute shear rate gammaDot, which is a vector because it is evaluated at a gaussian quadrature point 
    // for that compute velocity gradient
    vec_dbl_ptr_Type gammaDot(new vec_dbl_Type(weights->size(),0.0)); //gammaDot->at(j) j=0...weights 
    computeShearRate( dPhiTrans, gammaDot, dim); // updates gammaDot using velcoity solution 

         // Initialize some helper vectors/matrices
    double v11, v12, v13, v21, v22, v23, v31, v32, v33, value1_j, value2_j, value3_j , value1_i, value2_i, value3_i, viscosity_atw;

    viscosity_atw = 0.;

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

                       this->materialModel->evaluateFunction(this->params_,  gammaDot->at(w), viscosity_atw);

                    
                       // Construct entries - we go over all quadrature points and if j is updated we set v11 etc. again to zero
                        v11 = v11 + viscosity_atw *  weights->at(w)*(2.0*value1_j*value1_i+value2_j*value2_i+value3_j*value3_i);
                        v12 = v12 + viscosity_atw *  weights->at(w)*(value2_i*value1_j);
                        v13 = v13 + viscosity_atw *  weights->at(w)*(value3_i*value1_j);

                        v21 = v21 + viscosity_atw *  weights->at(w)*(value1_i*value2_j);
                        v22=  v22 + viscosity_atw *  weights->at(w)*(value1_i*value1_j+2.0*value2_j*value2_i+value3_j*value3_i);
                        v23 = v23 + viscosity_atw *  weights->at(w)*(value3_i*value2_j);

                        v31 = v31 + viscosity_atw *  weights->at(w)*(value1_i*value3_j);
                        v32 = v32 + viscosity_atw *  weights->at(w)*(value2_i*value3_j);
                        v33 = v33 + viscosity_atw *  weights->at(w)*(value1_i*value1_j+value2_i*value2_j+2.0*value3_i*value3_j);

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





//Extra Derivative of Extra stress term resulting from chosen nonlinear non-newtonian model  -----
//Same structure and functions as in assemblyStress 
// (dPhiTrans_i : [-1/4 * deta/dgammaDot * dgammaDot/dTau * (dv^k + (dvh^k)^T)^2 * ( dPhiTrans_j + (dPhiTrans_j)^T)    ]
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::assemblyStressDev(SmallMatrixPtr_Type &elementMatrix) {

	int dim = this->getDim();
	int numNodes= this->numNodesVelocity_;
	int Grad =2; // Needs to be fixed	
	string FEType = this->FETypeVelocity_;
	int dofs = this->dofsVelocity_; // for pressure it would be 1 

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    UN deg = Helper::determineDegree(dim,FEType,Grad); 
    Helper::getDPhi(dPhi, weights, dim, FEType, deg); 
    

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    buildTransformation(B);
    detB = B.computeInverse(Binv); 
    absDetB = std::fabs(detB); 

    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    applyBTinv( dPhi, dPhiTrans, Binv ); 

    TEUCHOS_TEST_FOR_EXCEPTION(dim == 1,std::logic_error, "AssemblyStress Not implemented for dim=1");
    TEUCHOS_TEST_FOR_EXCEPTION(dim == 3,std::logic_error, "AssemblyStress Not implemented for dim=3 not corrected");
    if (dim == 2)
    {
    //************************************
    //************************************
    // Due to the extra term related to the Gaetaeux-derivative there arise prefactors (a, b, fab) which depend on the velocity gradients solutions
    // which therefore also have to be computed here therefore we compute it directly here
    vec_dbl_Type u11(dPhiTrans.size(), -1.); // should correspond to du/dx at each quadrature point
    vec_dbl_Type u12(dPhiTrans.size(), -1.); // should correspond to du/dy at each quadrature point
    vec_dbl_Type u21(dPhiTrans.size(), -1.); // should correspond to dv/dx at each quadrature point
    vec_dbl_Type u22(dPhiTrans.size(), -1.); // should correspond to dv/dy at each quadrature point
    vec_dbl_ptr_Type 	gammaDot(new vec_dbl_Type(weights->size(),0.0)); //gammaDot->at(j) j=0...weights 

    vec_dbl_ptr_Type 	cross_terms(new vec_dbl_Type(weights->size(),0.0));   // prefactor a= 4*(du/dx)^2 + (du/dy+dv/dx)^2
    vec_dbl_ptr_Type 	du_dx(new vec_dbl_Type(weights->size(),0.0));   // prefactor b= 4*(dv/dy)^2 + (du/dy+dv/dx)^2
    vec_dbl_ptr_Type 	dv_dy(new vec_dbl_Type(weights->size(),0.0));   // prefactor b= 4*(dv/dy)^2 + (du/dy+dv/dx)^2

    for (UN w=0; w<dPhiTrans.size(); w++)
     { //quads points
      // set again to zero 
      u11[w] = 0.0;
      u12[w] = 0.0;
      u21[w] = 0.0;
      u22[w] = 0.0;
            for (UN i=0; i < dPhiTrans[0].size(); i++) 
            { // loop unrolling
                LO index1 = dim * i + 0; // x
                LO index2 = dim * i + 1; // y 
                // uLoc[d][w] += this->solution_[index] * phi->at(w).at(i);
                u11[w] += this->solution_[index1] * dPhiTrans[w][i][0]; // u*dphi_dx
                u12[w] += this->solution_[index1] * dPhiTrans[w][i][1]; // because we are in 2D , 0 and 1 
                u21[w] += this->solution_[index2] * dPhiTrans[w][i][0];
                u22[w] += this->solution_[index2] * dPhiTrans[w][i][1];
                
            }
            gammaDot->at(w) = sqrt(2.0*u11[w]*u11[w]+ 2.0*u22[w]*u22[w] + (u12[w]+u21[w])*(u12[w]+u21[w]));

            cross_terms->at(w)= 0.5*(u12[w]+u21[w]);
            du_dx->at(w)= (u11[w]);
            dv_dy->at(w)= (u22[w]);

     }
     //*******************************
        
        // Initialize some helper vectors/matrices
        double v11, v12, v21, v22, value1_j, value2_j , value1_i, value2_i, deta_dgamma_dgamma_dtau ;
    
        deta_dgamma_dgamma_dtau =0.;

    // Construct element matrices 
    for (UN i=0; i < numNodes; i++) {
       // Teuchos::Array<SC> value(dPhiTrans[0].size(), 0. ); // dPhiTrans[0].size() is 3        
      
        for (UN j=0; j < numNodes; j++) {
        // Reset values
        v11 = 0.0;v12 = 0.0;v21 = 0.0;v22 = 0.0;

            // Only the part  deta/dgammaDot is different for all shear thinning models (because we make the assumption of incompressibility)? NO ALSO DIFFERENT FOR EXAMPLE FOR CASSON SO CONSIDERING YIELD STRESS
            // but we put the two terms together because then we can multiply them together and get e.g. for carreau yasuda  : gammaDot^{a-2.0} which is for a=2.0 equals 0 and we do not have to worry about the problem what if gammaDot = 0.0
            for (UN w=0; w<dPhiTrans.size(); w++) {

                 value1_j = dPhiTrans[w][j][0]; // so this corresponds to d\phi_j/dx
                 value2_j = dPhiTrans[w][j][1]; // so this corresponds to d\phi_j/dy

                 value1_i = dPhiTrans[w][i][0]; // so this corresponds to d\phi_i/dx
                 value2_i = dPhiTrans[w][i][1]; // so this corresponds to d\phi_i/dy

                // viscosity function evaluated where we consider the dynamic viscosity!!
                this->materialModel->evaluateDerivative(this->params_,  gammaDot->at(w), deta_dgamma_dgamma_dtau);
	
                //if (i==0 && j==0) viscosity_averageOverT += viscosity_atw; // see below ### 

                 v11 = v11 + -2.0*deta_dgamma_dgamma_dtau * weights->at(w) *(value1_i*value1_j*du_dx->at(w)*du_dx->at(w) + value1_i*value2_j*du_dx->at(w)*cross_terms->at(w) + value2_i*value1_j*du_dx->at(w)*cross_terms->at(w) + value2_i*value2_j*cross_terms->at(w)*cross_terms->at(w));
                 v12 = v12 + -2.0*deta_dgamma_dgamma_dtau * weights->at(w) *(value1_i*value1_j*du_dx->at(w)*cross_terms->at(w) + value1_i*value2_j*du_dx->at(w)*dv_dy->at(w) + value2_i*value1_j*cross_terms->at(w)*cross_terms->at(w) + value2_i*value2_j*dv_dy->at(w)*cross_terms->at(w));
                 v21 = v21 + -2.0*deta_dgamma_dgamma_dtau * weights->at(w) *(value1_i*value1_j*du_dx->at(w)*cross_terms->at(w) + value1_i*value2_j*cross_terms->at(w)*cross_terms->at(w) + value2_i*value1_j*du_dx->at(w)*dv_dy->at(w) + value2_i*value2_j*dv_dy->at(w)*cross_terms->at(w));
                 v22 = v22 + -2.0*deta_dgamma_dgamma_dtau * weights->at(w) *(value1_i*value1_j*cross_terms->at(w)*cross_terms->at(w) + value1_i*value2_j*dv_dy->at(w)*cross_terms->at(w) + value2_i*value1_j*dv_dy->at(w)*cross_terms->at(w) + value2_i*value2_j*dv_dy->at(w)*dv_dy->at(w));

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
    //***************************************************************************
    //***************************************************************************
    else if (dim == 3)
    {
          //************************************
    // Compute shear rate gammaDot, which is a vector because it is evaluated at a gaussian quadrature point 
    // for that compute velocity gradient
    vec_dbl_Type u11(dPhiTrans.size(), -1.); // should correspond to du/dx at each quadrature point
    vec_dbl_Type u12(dPhiTrans.size(), -1.); // should correspond to du/dy at each quadrature point
    vec_dbl_Type u13(dPhiTrans.size(), -1.); // should correspond to du/dz at each quadrature point
    
    vec_dbl_Type u21(dPhiTrans.size(), -1.); // should correspond to dv/dx at each quadrature point
    vec_dbl_Type u22(dPhiTrans.size(), -1.); // should correspond to dv/dy at each quadrature point
    vec_dbl_Type u23(dPhiTrans.size(), -1.); // should correspond to dv/dz at each quadrature point

    vec_dbl_Type u31(dPhiTrans.size(), -1.); // should correspond to dw/dx at each quadrature point
    vec_dbl_Type u32(dPhiTrans.size(), -1.); // should correspond to dw/dy at each quadrature point
    vec_dbl_Type u33(dPhiTrans.size(), -1.); // should correspond to dw/dz at each quadrature point

    vec_dbl_ptr_Type gammaDot(new vec_dbl_Type(weights->size(),0.0)); //gammaDot->at(j) j=0...weights 

    vec_dbl_ptr_Type 	a1(new vec_dbl_Type(weights->size(),0.0)); // prefactor a= 4*(du/dx)^2 + (du/dy+dv/dx)^2 +  (dw/dx+du/dz)^2
    vec_dbl_ptr_Type 	b1(new vec_dbl_Type(weights->size(),0.0)); // prefactor b= 4*(dv/dy)^2 + (du/dy+dv/dx)^2 +  (dw/dy+dv/dz)^2
    vec_dbl_ptr_Type 	c1(new vec_dbl_Type(weights->size(),0.0)); // prefactor c= 4*(dw/dz)^2 + (dw/dx+du/dz)^2 +  (dw/dy+dv/dz)^2
    vec_dbl_ptr_Type 	d1(new vec_dbl_Type(weights->size(),0.0)); // prefactor d= 2*(du/dy+dv/dx)*(du/dx+dv/dy) + (dw/dx+du/dz)*(dw/dy+dv/dz)
    vec_dbl_ptr_Type 	e1(new vec_dbl_Type(weights->size(),0.0)); // prefactor e= 2*(dw/dx+du/dz)*(du/dx+dw/dz) + (du/dy+dv/dx)*(dw/dy+dv/dz)
    vec_dbl_ptr_Type 	g1(new vec_dbl_Type(weights->size(),0.0)); // prefactor g= 2*(dw/dy+dv/dz)*(dv/dy+dw/dz) + (du/dy+dv/dx)*(dw/dx+du/dz)
   
    for (UN w=0; w<dPhiTrans.size(); w++)
     { //quads points
   
      u11[w] = 0.0;
      u12[w] = 0.0;
      u13[w] = 0.0;
      u21[w] = 0.0;
      u22[w] = 0.0;
      u23[w] = 0.0;
      u31[w] = 0.0;
      u32[w] = 0.0;
      u33[w] = 0.0;
      
            for (UN i=0; i < dPhiTrans[0].size(); i++) 
            {
                LO index1 = dim * i + 0; //x
                LO index2 = dim * i + 1; //y 
                LO index3 = dim * i + 2; //z
                // uLoc[d][w] += this->solution_[index] * phi->at(w).at(i);
                u11[w] += this->solution_[index1] * dPhiTrans[w][i][0]; // u*dphi_dx
                u12[w] += this->solution_[index1] * dPhiTrans[w][i][1]; // because we are in 3D , 0 and 1, 2 
                u13[w] += this->solution_[index1] * dPhiTrans[w][i][2]; 
                u21[w] += this->solution_[index2] * dPhiTrans[w][i][0]; // v*dphi_dx
                u22[w] += this->solution_[index2] * dPhiTrans[w][i][1];
                u23[w] += this->solution_[index2] * dPhiTrans[w][i][2];
                u31[w] += this->solution_[index3] * dPhiTrans[w][i][0]; // w*dphi_dx
                u32[w] += this->solution_[index3] * dPhiTrans[w][i][1];
                u33[w] += this->solution_[index3] * dPhiTrans[w][i][2];
                
            }
            gammaDot->at(w) = sqrt(2.0*u11[w]*u11[w]+ 2.0*u22[w]*u22[w] + 2.0*u33[w]*u33[w] +  (u12[w]+u21[w])*(u12[w]+u21[w])   + (u13[w]+u31[w])*(u13[w]+u31[w]) + (u23[w]+u32[w])*(u23[w]+u32[w]) );
            a1->at(w)= 4.0*(u11[w]*u11[w]) + (u12[w]+u21[w])*(u12[w]+u21[w])+(u31[w]+u13[w])*(u31[w]+u13[w]);
            b1->at(w)= 4.0*(u22[w]*u22[w]) + (u12[w]+u21[w])*(u12[w]+u21[w])+(u32[w]+u23[w])*(u32[w]+u23[w]);
            c1->at(w)= 4.0*(u33[w]*u33[w]) + (u31[w]+u13[w])*(u31[w]+u13[w])+(u32[w]+u23[w])*(u32[w]+u23[w]);
            d1->at(w)= 2.0*(u12[w]+u21[w])*(u11[w]+u22[w])+(u31[w]+u13[w])*(u32[w]+u23[w]);
            e1->at(w)= 2.0*(u31[w]+u13[w])*(u11[w]+u33[w])+(u12[w]+u21[w])*(u32[w]+u23[w]);
            g1->at(w)= 2.0*(u32[w]+u23[w])*(u22[w]+u33[w])+(u12[w]+u21[w])*(u31[w]+u13[w]);

      }


      // Initialize some helper vectors/matrices
        double v11, v12, v13, v21, v22, v23, v31, v32, v33, value1_j, value2_j, value3_j , value1_i, value2_i, value3_i,deta_dgamma_dgamma_dtau;
        SmallMatrix<double> e1i(dim);
        SmallMatrix<double> e2i(dim);
        SmallMatrix<double> e3i(dim);
        SmallMatrix<double> e1j(dim);
        SmallMatrix<double> e2j(dim);
        SmallMatrix<double> e3j(dim);

         deta_dgamma_dgamma_dtau =0.;

     // Construct element matrices 
     for (UN i=0; i < numNodes; i++) 
     {
       // Teuchos::Array<SC> value(dPhiTrans[0].size(), 0. ); // dPhiTrans[0].size() is 3        
      
        for (UN j=0; j < numNodes; j++) 
        {
         // Reset values
         v11 = 0.0;v12 = 0.0;v13=0.0; v21 = 0.0;v22 = 0.0;v23=0.0;v31=0.0;v32=0.0;v33=0.0;

            // So in general compute the components of eta*[ dPhiTrans_i : ( dPhiTrans_j + (dPhiTrans_j)^T )]
            for (UN w=0; w<dPhiTrans.size(); w++) 
            {

                        value1_j = dPhiTrans.at(w).at(j).at(0); // so this corresponds to d\phi_j/dx
                        value2_j = dPhiTrans.at(w).at(j).at(1); // so this corresponds to d\phi_j/dy
                        value3_j = dPhiTrans.at(w).at(j).at(2); // so this corresponds to d\phi_j/dz

                        value1_i = dPhiTrans.at(w).at(i).at(0); // so this corresponds to d\phi_i/dx
                        value2_i = dPhiTrans.at(w).at(i).at(1); // so this corresponds to d\phi_i/dy
                        value3_i = dPhiTrans.at(w).at(i).at(2); // so this corresponds to d\phi_i/dz

                        e1j[0][0] = 2.*value1_j;
                        e1j[0][1] = value2_j;
                        e1j[0][2] = value3_j;
                        e1j[1][0] = value2_j;
                        e1j[2][0] = value3_j;

                        e1i[0][0] = value1_i*a1->at(w);
                        e1i[0][1] = value2_i*a1->at(w);
                        e1i[0][2] = value3_i*a1->at(w);
                        e1i[1][0] = value1_i*d1->at(w);
                        e1i[1][1] = value2_i*d1->at(w);
                        e1i[1][2] = value3_i*d1->at(w);
                        e1i[2][0] = value1_i*e1->at(w);
                        e1i[2][1] = value2_i*e1->at(w);
                        e1i[2][2] = value3_i*e1->at(w);
                        
                        e2j[1][0] = value1_j;
                        e2j[1][1] = 2.*value2_j;
                        e2j[1][2] = value3_j;
                        e2j[0][1] = value1_j;
                        e2j[2][1] = value3_j;

                        e2i[0][0] = value1_i*d1->at(w);
                        e2i[0][1] = value2_i*d1->at(w);
                        e2i[0][2] = value3_i*d1->at(w);
                        e2i[1][0] = value1_i*b1->at(w);
                        e2i[1][1] = value2_i*b1->at(w);
                        e2i[1][2] = value3_i*b1->at(w);
                        e2i[2][0] = value1_i*g1->at(w);
                        e2i[2][1] = value2_i*g1->at(w);
                        e2i[2][2] = value3_i*g1->at(w);

                        e3j[2][0] = value1_j;
                        e3j[2][1] = value2_j;
                        e3j[2][2] = 2.*value3_j;
                        e3j[0][2] = value1_j;
                        e3j[1][2] = value2_j;

                        e3i[0][0] = value1_i*e1->at(w);
                        e3i[0][1] = value2_i*e1->at(w);
                        e3i[0][2] = value3_i*e1->at(w);
                        e3i[1][0] = value1_i*g1->at(w);
                        e3i[1][1] = value2_i*g1->at(w);
                        e3i[1][2] = value3_i*g1->at(w);
                        e3i[2][0] = value1_i*c1->at(w);
                        e3i[2][1] = value2_i*c1->at(w);
                        e3i[2][2] = value3_i*c1->at(w);

                       this->materialModel->evaluateDerivative(this->params_,  gammaDot->at(w), deta_dgamma_dgamma_dtau);
	
                       // Construct entries - we go over all quadrature points and if j is updated we set v11 etc. again to zero
                        v11 = v11 + (-0.25)*deta_dgamma_dgamma_dtau  *  weights->at(w) * e1i.innerProduct(e1j);
                        v12 = v12 + (-0.25)*deta_dgamma_dgamma_dtau  *  weights->at(w) * e1i.innerProduct(e2j);
                        v13 = v13 + (-0.25)*deta_dgamma_dgamma_dtau  *  weights->at(w) * e1i.innerProduct(e3j);

                        v21 = v21 + (-0.25)*deta_dgamma_dgamma_dtau  *  weights->at(w) * e2i.innerProduct(e1j);
                        v22 = v22 + (-0.25)*deta_dgamma_dgamma_dtau  *  weights->at(w) * e2i.innerProduct(e2j);
                        v23 = v23 + (-0.25)*deta_dgamma_dgamma_dtau  *  weights->at(w) * e2i.innerProduct(e3j);

                        v31 = v31 + (-0.25)*deta_dgamma_dgamma_dtau  *  weights->at(w) * e3i.innerProduct(e1j);
                        v32 = v32 + (-0.25)*deta_dgamma_dgamma_dtau  *  weights->at(w) * e3i.innerProduct(e2j);
                        v33 = v33 + (-0.25)*deta_dgamma_dgamma_dtau  *  weights->at(w) * e3i.innerProduct(e3j);

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

                   // Put values on the right position in element matrix - d=2 because we are in two dimensional case
                   // [v11  v12  v13]
                   // [v21  v22  v23]
                   // [v31  v32  v33]
                  (*elementMatrix)[i*dofs][j*dofs]      = v11; // d=0, first dimension
                  (*elementMatrix)[i*dofs][j*dofs+1]    = v12;
                  (*elementMatrix)[i*dofs][j*dofs+2]    = v13;
                  (*elementMatrix)[i*dofs+1][j*dofs]    = v21;
                  (*elementMatrix)[i*dofs +1][j*dofs+1] = v22; //d=1, second dimension
                  (*elementMatrix)[i*dofs +1][j*dofs+2] = v23; //d=1, second dimension
                  (*elementMatrix)[i*dofs+2][j*dofs]    = v31;
                  (*elementMatrix)[i*dofs +2][j*dofs+1] = v32; //d=2, third dimension
                  (*elementMatrix)[i*dofs +2][j*dofs+2] = v33; //d=2, third dimension

        }// loop end over j node 
    
    }// loop end over i node 

    }// end if dim = 3

}





// Boundary integral over Neumann boundary resulting from the fact that we want our outflow boundary conditions 
// as in the case for reduced stress tensor - therefore we have to subtract the boundary integral
// int_NeumannBoundary ( \nabla u)^T n \cdot w dNeumannBoundary
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::assemblyNeumannBoundaryTerm(SmallMatrixPtr_Type &elementMatrix) //, dblVecPtr normalVector) 
{


	int dim = this->getDim();
	int numNodes= this->numNodesVelocity_;
	string FEType = this->FETypeVelocity_;
	int dofs = this->dofsVelocity_; // 

    vec3D_dbl_ptr_Type 	dPhi; // derivative of basisfunction
    vec2D_dbl_ptr_Type phi; // basisfunction 

     // Compute phi and derivative of phi at quadrature points
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    buildTransformation(B);
    detB = B.computeInverse(Binv); // The function computeInverse returns a double value corrsponding to determinant of B
    absDetB = std::fabs(detB); // absolute value of B    

    // Now we have to compute the weights and the quadrature points for our line integral (2D), surface integral (3D)
    // where we want to evaluate our dPhi, phi 
    vec2D_dbl_ptr_Type QuadPts;
    vec_dbl_ptr_Type QuadW(new vec_dbl_Type(this->surfaceElement_QuadratureWeightsPhysicalSpace->size(),0.0));
    // So we are now mapping back the quadrature points from the physical edge element to the reference element 
    // via the inverse mapping xi = Binv(x-tau) where tau is p0
    // tau = [this->nodesRefConfig_.at(0).at(0) ;  this->nodesRefConfig_.at(0).at(1)]
    QuadPts.reset(new vec2D_dbl_Type(this->surfaceElement_QuadraturePointsPhysicalSpace->size(),vec_dbl_Type(dim,0.0)));
  
    for(int l=0; l< this->surfaceElement_QuadraturePointsPhysicalSpace->size(); l++)
    {
        QuadW->at(l)=this->surfaceElement_QuadratureWeightsPhysicalSpace->at(l);
        for (int p = 0; p < dim; p++)
        {
            for (int q=0; q < dim; q++)
            {
            QuadPts->at(l).at(p)+= Binv[p][q]*(this->surfaceElement_QuadraturePointsPhysicalSpace->at(l).at(q)-this->nodesRefConfig_.at(0).at(q) );
            }
        }    
    }

    // SO NOW WE HAVE TO CATCH THE CASE THAT IF WE PROJECTED ONTO THE DIAGONAL LINE WE HAVE TO MULTIPLY BY THE LENGTH SO SQRT(2)
    // HOW TO CHECK IF WE ARE ON THE DIAGONAL LINE? We check if the x and y component of a quafrature point a nonzero
    if (  QuadPts->at(0).at(0) != 0. && QuadPts->at(0).at(1) != 0.) // if the x and y component of quadrature point are non-zero we are on the diagonal but also only if the quadrature point was not defined in corners of element
    {
        for(int l=0; l< this->surfaceElement_QuadraturePointsPhysicalSpace->size(); l++)
        {
        QuadW->at(l)=std::sqrt(2.0)*QuadW->at(l);
         for (int q=0; q < dim; q++)
            {
                QuadPts->at(l).at(q)=std::sqrt(2.0)*QuadPts->at(l).at(q);
            }
        }
    }

    Helper::getPhi(phi, QuadW,QuadPts, dim, FEType);
    Helper::getDPhi(dPhi, QuadW,QuadPts, dim, FEType); 

    // In 2D for our line integral we get change in length of edges
    // detB_line = std::sqrt( std::pow( B[0][1] , 2.0) + std::pow( B[1][1], 2.0) );
    SC detB_line;
    detB_line=this->surfaceElement_MappingChangeInArea;
    // dPhiTrans are the transorfmed basisfunctions, so B^(-T) * \grad_phi bzw. \grad_phi^T * B^(-1) Corresponds to \hat{grad_phi}.
    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    applyBTinv( dPhi, dPhiTrans, Binv ); // so dPhiTrans corresponds now to our basisfunction in natural coordinates


    // Compute shear rate gammaDot, which is a vector because it is evaluated at a gaussian quadrature point for that compute velocity gradient
    vec_dbl_ptr_Type gammaDot(new vec_dbl_Type(QuadW->size(),0.0)); //gammaDot->at(j) j=0...weights 
    computeShearRate( dPhiTrans, gammaDot, dim); // updates gammaDot using velocity solution 
    double viscosity_atw = 0.;

    TEUCHOS_TEST_FOR_EXCEPTION(dim == 1,std::logic_error, "AssemblyNeumannBoundaryTerm Not implemented for dim=1");
    // 2D 
    if (dim == 2)
    {
            double v11, v12, v21, v22 ; // helper values for entries 

                // loop over basis functions
                for (UN i=0; i < phi->at(0).size(); i++) 
                {
                    for (UN j=0; j < numNodes; j++)
                    {
                        // Reset values
                        v11 = 0.0;v12 = 0.0;v21 = 0.0;v22 = 0.0;

                        // loop over basis functions quadrature points
                        for (UN w=0; w<phi->size(); w++) 
                        {
                            this->materialModel->evaluateFunction(this->params_,  gammaDot->at(w), viscosity_atw);

                            v11 = v11 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][0] * this->surfaceElement_OutwardNormal[0] * (*phi)[w][i]); // xx contribution: 
                            v12 = v12 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][0] * this->surfaceElement_OutwardNormal[1] * (*phi)[w][i]); // xy contribution:  
                            v21 = v21 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][1] * this->surfaceElement_OutwardNormal[0] * (*phi)[w][i]); // yx contribution:  
                            v22 = v22 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][1] * this->surfaceElement_OutwardNormal[1] * (*phi)[w][i]); // yy contribution:                     
                  
                        } // End loop over quadrature points

                        //multiply determinant from transformation
                        v11 *= detB_line; 
                        v12 *= detB_line;
                        v21 *= detB_line; 
                        v22 *= detB_line;
            
                        // Put values on the right position in element matrix - d=2 because we are in two dimensional case
                        // [v11  v12 ]
                        // [v21  v22 ]
                        (*elementMatrix)[i*dofs][j*dofs]   = v11; // d=0, first dimension
                        (*elementMatrix)[i*dofs][j*dofs+1] = v12;  //
                        (*elementMatrix)[i*dofs+1][j*dofs] = v21;
                        (*elementMatrix)[i*dofs +1][j*dofs+1] =v22; //d=1, second dimension
                                                

                    } // End loop over j nodes

                } // End loop over i nodes
    } // End dim==2
    else if (dim==3)
    {
    double v11, v12, v13, v21, v22, v23, v31, v32, v33; // helper values for entries
                    // loop over basis functions
                for (UN i=0; i < phi->at(0).size(); i++) 
                {
                    for (UN j=0; j < numNodes; j++)
                    {
                        // Reset values
                         v11 = 0.0;v12 = 0.0;v13=0.0; v21 = 0.0;v22 = 0.0;v23=0.0;v31=0.0;v32=0.0;v33=0.0;

                        // loop over basis functions quadrature points
                        for (UN w=0; w<phi->size(); w++) 
                        {
                            this->materialModel->evaluateFunction(this->params_,  gammaDot->at(w), viscosity_atw);

                            v11 = v11 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][0] * this->surfaceElement_OutwardNormal[0] * (*phi)[w][i]); // xx contribution: 
                            v12 = v12 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][0] * this->surfaceElement_OutwardNormal[1] * (*phi)[w][i]); // xy contribution: 
                            v13 = v13 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][0] * this->surfaceElement_OutwardNormal[2] * (*phi)[w][i]); // xz contribution:  
                             
                            
                            v21 = v21 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][1] * this->surfaceElement_OutwardNormal[0] * (*phi)[w][i]); // yx contribution:  
                            v22 = v22 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][1] * this->surfaceElement_OutwardNormal[1] * (*phi)[w][i]); // yy contribution:  
                            v23 = v23 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][1] * this->surfaceElement_OutwardNormal[2] * (*phi)[w][i]); // yz contribution:                     
                                     
                            v31 = v31 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][2] * this->surfaceElement_OutwardNormal[0] * (*phi)[w][i]); // zx contribution:  
                            v32 = v32 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][2] * this->surfaceElement_OutwardNormal[1] * (*phi)[w][i]); // zy contribution:  
                            v33 = v33 + -1.0*(viscosity_atw *  QuadW->at(w)* dPhiTrans[w][j][2] * this->surfaceElement_OutwardNormal[2] * (*phi)[w][i]); // zz contribution:  

                  

                        } // End loop over quadrature points

                        //multiply determinant from transformation
                        v11 *= detB_line; 
                        v12 *= detB_line;
                        v13 *= detB_line;
                        v21 *= detB_line; 
                        v22 *= detB_line;
                        v23 *= detB_line;
                        v31 *= detB_line; 
                        v32 *= detB_line;
                        v33 *= detB_line;
            
                        // Put values on the right position in element matrix - d=2 because we are in two dimensional case
                        // [v11  v12  v13]
                        // [v21  v22  v23]
                        // [v31  v32  v33]
                       (*elementMatrix)[i*dofs][j*dofs]      = v11; // d=0, first dimension
                       (*elementMatrix)[i*dofs][j*dofs+1]    = v12;
                       (*elementMatrix)[i*dofs][j*dofs+2]    = v13;
                       (*elementMatrix)[i*dofs+1][j*dofs]    = v21;
                       (*elementMatrix)[i*dofs +1][j*dofs+1] = v22; //d=1, second dimension
                       (*elementMatrix)[i*dofs +1][j*dofs+2] = v23; //d=1, second dimension
                       (*elementMatrix)[i*dofs+2][j*dofs]    = v31;
                       (*elementMatrix)[i*dofs +2][j*dofs+1] = v32; //d=2, third dimension
                       (*elementMatrix)[i*dofs +2][j*dofs+2] = v33; //d=2, third dimension
                                                

                    } // End loop over j nodes

                } // End loop over i nodes    
    }//end if 3d
    

} // Function End loop 





// Boundary integral over Neumann boundary resulting from the fact that we want our outflow boundary conditions 
// as in the case for reduced stress tensor - therefore we have to subtract the boundary integral
// int_NeumannBoundary ( \nabla u)^T n \cdot w dNeumannBoundary

template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::assemblyNeumannBoundaryTermDev(SmallMatrixPtr_Type &elementMatrix) //, dblVecPtr normalVector) 
{


	int dim = this->getDim();
	int numNodes= this->numNodesVelocity_;
	string FEType = this->FETypeVelocity_;
	int dofs = this->dofsVelocity_; // 

    vec3D_dbl_ptr_Type 	dPhi; // derivative of basisfunction
    vec2D_dbl_ptr_Type phi; // basisfunction 

     // Compute phi and derivative of phi at quadrature points
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    buildTransformation(B);
    detB = B.computeInverse(Binv); // The function computeInverse returns a double value corrsponding to determinant of B
    absDetB = std::fabs(detB); // absolute value of B    


    // Now we have to compute the weights and the quadrature points for our line integral (2D), surface integral (3D)
    // where we want to evaluate our dPhi, phi 
    vec2D_dbl_ptr_Type QuadPts;
    vec_dbl_ptr_Type QuadW(new vec_dbl_Type(this->surfaceElement_QuadratureWeightsPhysicalSpace->size(),0.0));
    // So we are now mapping back the quadrature points from the physical edge element to the reference element 
    // via the inverse mapping xi = Binv(x-tau) where tau is p0
    // tau = [this->nodesRefConfig_.at(0).at(0) ;  this->nodesRefConfig_.at(0).at(1)]
    QuadPts.reset(new vec2D_dbl_Type(this->surfaceElement_QuadraturePointsPhysicalSpace->size(),vec_dbl_Type(dim,0.0)));
  
    for(int l=0; l< this->surfaceElement_QuadraturePointsPhysicalSpace->size(); l++)
    {
        QuadW->at(l)=this->surfaceElement_QuadratureWeightsPhysicalSpace->at(l);
        for (int p = 0; p < dim; p++)
        {
            for (int q=0; q < dim; q++)
            {
            QuadPts->at(l).at(p)+= Binv[p][q]*(this->surfaceElement_QuadraturePointsPhysicalSpace->at(l).at(q)-this->nodesRefConfig_.at(0).at(q) );
            }
        }    
    }
      // SO NOW WE HAVE TO CATCH THE CASE THAT IF WE PROJECTED ONTO THE DIAGONAL LINE WE HAVE TO MULTIPLY BY THE LENGTH SO SQRT(2)
    // HOW TO CHECK IF WE ARE ON THE DIAGONAL LINE? We check if the x and y component of a quafrature point a nonzero
    if (  QuadPts->at(0).at(0) != 0. && QuadPts->at(0).at(1) != 0.) // if the x and y component of quadrature point are non-zero we are on the diagonal but also only if the quadrature point was not defined in corners of element
    {
        for(int l=0; l< this->surfaceElement_QuadraturePointsPhysicalSpace->size(); l++)
        {
        QuadW->at(l)=std::sqrt(2.0)*QuadW->at(l);
         for (int q=0; q < dim; q++)
            {
                QuadPts->at(l).at(q)=std::sqrt(2.0)*QuadPts->at(l).at(q);
            }
        }
    }

    Helper::getPhi(phi, QuadW,QuadPts, dim, FEType);
    Helper::getDPhi(dPhi, QuadW,QuadPts, dim, FEType); 

    // In 2D for our line integral we get change in length of edges
    // detB_line = std::sqrt( std::pow( B[0][1] , 2.0) + std::pow( B[1][1], 2.0) );
    SC detB_line;
    detB_line=this->surfaceElement_MappingChangeInArea;
    
    // dPhiTrans are the transorfmed basisfunctions, so B^(-T) * \grad_phi bzw. \grad_phi^T * B^(-1) Corresponds to \hat{grad_phi}.
    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    applyBTinv( dPhi, dPhiTrans, Binv ); // so dPhiTrans corresponds now to our basisfunction in natural coordinates

    // Compute shear rate gammaDot, which is a vector because it is evaluated at a gaussian quadrature point for that compute velocity gradient
    vec_dbl_ptr_Type gammaDot(new vec_dbl_Type(QuadW->size(),0.0)); //gammaDot->at(j) j=0...weights 
    computeShearRate( dPhiTrans, gammaDot, dim); // updates gammaDot using velocity solution 

    double viscosity_atw = 0.;

    TEUCHOS_TEST_FOR_EXCEPTION(dim == 1,std::logic_error, "AssemblyNeumannBoundaryTerm Not implemented for dim=1");
    
    if (dim == 2)
    {
     double v11, v12, v21, v22,deta_dgamma_dgamma_dtau;
     deta_dgamma_dgamma_dtau =0.;
     //************************************
     //************************************
     // Due to the extra term related to the Gaetaeux-derivative there arise prefactors which depend on the velocity gradients 
     // which therefore also have to be computed here
     vec_dbl_Type u11(dPhiTrans.size(), -1.); // should correspond to du/dx at each quadrature point
     vec_dbl_Type u12(dPhiTrans.size(), -1.); // should correspond to du/dy at each quadrature point
     vec_dbl_Type u21(dPhiTrans.size(), -1.); // should correspond to dv/dx at each quadrature point
     vec_dbl_Type u22(dPhiTrans.size(), -1.); // should correspond to dv/dy at each quadrature point
     vec_dbl_ptr_Type 	gammaDot(new vec_dbl_Type(QuadW->size(),0.0)); //gammaDot->at(j) j=0...weights 

     vec_dbl_ptr_Type 	cross_terms(new vec_dbl_Type(QuadW->size(),0.0));   // prefactor a= 2*(du/dx)^2 + (du/dy+dv/dx)*dv/dx
     vec_dbl_ptr_Type 	du_dx(new vec_dbl_Type(QuadW->size(),0.0));   // prefactor b= (du/dy+dv/dx)*du/dx + 2*(dv/dx)*dv/dy
     vec_dbl_ptr_Type 	du_dy(new vec_dbl_Type(QuadW->size(),0.0));   // prefactor b= (du/dy+dv/dx)*du/dx + 2*(dv/dx)*dv/dy
     vec_dbl_ptr_Type 	dv_dx(new vec_dbl_Type(QuadW->size(),0.0));   // prefactor b= (du/dy+dv/dx)*du/dx + 2*(dv/dx)*dv/dy
     vec_dbl_ptr_Type 	dv_dy(new vec_dbl_Type(QuadW->size(),0.0));   // prefactor b= (du/dy+dv/dx)*du/dx + 2*(dv/dx)*dv/dy

     for (UN w=0; w<dPhiTrans.size(); w++)
     { //quads points
      // set again to zero 
      u11[w] = 0.0;
      u12[w] = 0.0;
      u21[w] = 0.0;
      u22[w] = 0.0;
            for (UN i=0; i < dPhiTrans[0].size(); i++) 
            { // loop unrolling
                LO index1 = dim * i + 0; // x
                LO index2 = dim * i + 1; // y 
                u11[w] += this->solution_[index1] * dPhiTrans[w][i][0]; // u*dphi_dx
                u12[w] += this->solution_[index1] * dPhiTrans[w][i][1]; // because we are in 2D , 0 and 1 
                u21[w] += this->solution_[index2] * dPhiTrans[w][i][0];
                u22[w] += this->solution_[index2] * dPhiTrans[w][i][1];
                
            }
            cross_terms->at(w)= 0.5*(u12[w]+u21[w]);
            du_dx->at(w)= (u11[w]);
            du_dy->at(w)= (u12[w]);
            dv_dy->at(w)= (u22[w]);
            dv_dx->at(w)= (u21[w]);

     }
                // loop over basis functions
     for (UN i=0; i < phi->at(0).size(); i++) 
        {
            for (UN j=0; j < numNodes; j++)
                {
                    // Reset values
                    v11 = 0.0;v12 = 0.0;v21 = 0.0;v22 = 0.0;

                        // loop over basis functions quadrature points
                        for (UN w=0; w<phi->size(); w++) 
                         {
                            this->materialModel->evaluateDerivative(this->params_,  gammaDot->at(w), deta_dgamma_dgamma_dtau);
                            v11 = v11 + deta_dgamma_dgamma_dtau*QuadW->at(w)*(  ( du_dx->at(w))*this->surfaceElement_OutwardNormal[0] + ( dv_dx->at(w)  )*this->surfaceElement_OutwardNormal[1] )*(du_dx->at(w)*dPhiTrans[w][j][0] + cross_terms->at(w)*dPhiTrans[w][j][1]) *((*phi)[w][i]); // xx contribution: 
                            v12 = v12 + deta_dgamma_dgamma_dtau*QuadW->at(w)*(  ( du_dx->at(w))*this->surfaceElement_OutwardNormal[0] + ( dv_dx->at(w)  )*this->surfaceElement_OutwardNormal[1] )*(cross_terms->at(w)*dPhiTrans[w][j][0] + dv_dy->at(w)*dPhiTrans[w][j][1]) *((*phi)[w][i]); // xy contribution:  
                            v21 = v21 + deta_dgamma_dgamma_dtau*QuadW->at(w)*(  ( du_dy->at(w))*this->surfaceElement_OutwardNormal[0] + ( dv_dy->at(w)  )*this->surfaceElement_OutwardNormal[1] )*(du_dx->at(w)*dPhiTrans[w][j][0] + cross_terms->at(w)*dPhiTrans[w][j][1]) *((*phi)[w][i]); // yx contribution:  
                            v22 = v22 + deta_dgamma_dgamma_dtau*QuadW->at(w)*(  ( du_dy->at(w))*this->surfaceElement_OutwardNormal[0] + ( dv_dy->at(w)  )*this->surfaceElement_OutwardNormal[1] )*(cross_terms->at(w)*dPhiTrans[w][j][0] + dv_dy->at(w)*dPhiTrans[w][j][1]) *((*phi)[w][i]);// yy contribution:                     
                  
                         } // End loop over quadrature points

                    //multiply determinant from transformation
                    v11 *= detB_line; 
                    v12 *= detB_line;
                    v21 *= detB_line; 
                    v22 *= detB_line;
            
                    // Put values on the right position in element matrix - d=2 because we are in two dimensional case
                    // [v11  v12 ]cross_terms->at(w)*dPhiTrans[w][j][0] + dv_dy->at(w)*dPhiTrans[w][j][1]
                    // [v21  v22 ]
                    (*elementMatrix)[i*dofs][j*dofs]   = v11; // d=0, first dimension
                    (*elementMatrix)[i*dofs][j*dofs+1] = v12;  //
                    (*elementMatrix)[i*dofs+1][j*dofs] = v21;
                    (*elementMatrix)[i*dofs +1][j*dofs+1] =v22; //d=1, second dimension
                                                

                } // End loop over j nodes

        } // End loop over i nodes
    } // End dim==2
    else if (dim==3)
    {
     vec_dbl_Type u11(dPhiTrans.size(), -1.); // should correspond to du/dx at each quadrature point
     vec_dbl_Type u12(dPhiTrans.size(), -1.); // should correspond to du/dy at each quadrature point
     vec_dbl_Type u13(dPhiTrans.size(), -1.); // should correspond to du/dz at each quadrature point
    
     vec_dbl_Type u21(dPhiTrans.size(), -1.); // should correspond to dv/dx at each quadrature point
     vec_dbl_Type u22(dPhiTrans.size(), -1.); // should correspond to dv/dy at each quadrature point
     vec_dbl_Type u23(dPhiTrans.size(), -1.); // should correspond to dv/dz at each quadrature point

     vec_dbl_Type u31(dPhiTrans.size(), -1.); // should correspond to dw/dx at each quadrature point
     vec_dbl_Type u32(dPhiTrans.size(), -1.); // should correspond to dw/dy at each quadrature point
     vec_dbl_Type u33(dPhiTrans.size(), -1.); // should correspond to dw/dz at each quadrature point

     vec_dbl_ptr_Type gammaDot(new vec_dbl_Type(QuadW->size(),0.0)); //gammaDot->at(j) j=0...weights 

     vec_dbl_ptr_Type 	a1(new vec_dbl_Type(QuadW->size(),0.0)); // prefactor a= 2*(du/dx)^2 + (du/dy+dv/dx)*xv/dy +  (dw/dx+du/dz)*dw/dx
     vec_dbl_ptr_Type 	b1(new vec_dbl_Type(QuadW->size(),0.0)); // prefactor b= 
     vec_dbl_ptr_Type 	c1(new vec_dbl_Type(QuadW->size(),0.0)); // prefactor c= 
     vec_dbl_ptr_Type 	d1(new vec_dbl_Type(QuadW->size(),0.0)); // prefactor d= 
     vec_dbl_ptr_Type 	e1(new vec_dbl_Type(QuadW->size(),0.0)); // prefactor e= 
     vec_dbl_ptr_Type 	f1(new vec_dbl_Type(QuadW->size(),0.0)); // prefactor e= 
     vec_dbl_ptr_Type 	g1(new vec_dbl_Type(QuadW->size(),0.0)); // prefactor g= 
     vec_dbl_ptr_Type 	h1(new vec_dbl_Type(QuadW->size(),0.0)); // prefactor g= 
    vec_dbl_ptr_Type 	i1(new vec_dbl_Type(QuadW->size(),0.0)); // prefactor g= 
   
   
     for (UN w=0; w<dPhiTrans.size(); w++)
     { //quads points
   
      u11[w] = 0.0;
      u12[w] = 0.0;
      u13[w] = 0.0;
      u21[w] = 0.0;
      u22[w] = 0.0;
      u23[w] = 0.0;
      u31[w] = 0.0;
      u32[w] = 0.0;
      u33[w] = 0.0;
      
            for (UN i=0; i < dPhiTrans[0].size(); i++) 
            {
                LO index1 = dim * i + 0; //x
                LO index2 = dim * i + 1; //y 
                LO index3 = dim * i + 2; //z
                // uLoc[d][w] += this->solution_[index] * phi->at(w).at(i);
                u11[w] += this->solution_[index1] * dPhiTrans[w][i][0]; // u*dphi_dx
                u12[w] += this->solution_[index1] * dPhiTrans[w][i][1]; // because we are in 3D , 0 and 1, 2 
                u13[w] += this->solution_[index1] * dPhiTrans[w][i][2]; 
                u21[w] += this->solution_[index2] * dPhiTrans[w][i][0]; // v*dphi_dx
                u22[w] += this->solution_[index2] * dPhiTrans[w][i][1];
                u23[w] += this->solution_[index2] * dPhiTrans[w][i][2];
                u31[w] += this->solution_[index3] * dPhiTrans[w][i][0]; // w*dphi_dx
                u32[w] += this->solution_[index3] * dPhiTrans[w][i][1];
                u33[w] += this->solution_[index3] * dPhiTrans[w][i][2];
                
            }
            gammaDot->at(w) = sqrt(2.0*u11[w]*u11[w]+ 2.0*u22[w]*u22[w] + 2.0*u33[w]*u33[w] +  (u12[w]+u21[w])*(u12[w]+u21[w])   + (u13[w]+u31[w])*(u13[w]+u31[w]) + (u23[w]+u32[w])*(u23[w]+u32[w]) );
        
            a1->at(w)= 2.0*(u11[w]*u11[w]) + (u12[w]+u21[w])*(u21[w])+(u31[w]+u13[w])*(u31[w]);
            b1->at(w)= 2.0*(u21[w]*u22[w]) + (u12[w]+u21[w])*(u11[w])+(u32[w]+u23[w])*(u31[w]);
            c1->at(w)= 2.0*(u33[w]*u31[w]) + (u13[w]+u31[w])*(u11[w])+(u32[w]+u23[w])*(u21[w]);
            d1->at(w)= 2.0*(u11[w]*u12[w]) + (u12[w]+u21[w])*(u22[w])+(u31[w]+u13[w])*(u32[w]);
            e1->at(w)= 2.0*(u22[w]*u22[w]) + (u12[w]+u21[w])*(u12[w])+(u32[w]+u23[w])*(u32[w]);
            f1->at(w)= 2.0*(u33[w]*u32[w]) + (u13[w]+u31[w])*(u12[w])+(u32[w]+u23[w])*(u22[w]);
            g1->at(w)= 2.0*(u11[w]*u13[w]) + (u12[w]+u21[w])*(u23[w])+(u31[w]+u13[w])*(u33[w]);
            h1->at(w)= 2.0*(u22[w]*u23[w]) + (u12[w]+u21[w])*(u13[w])+(u32[w]+u23[w])*(u33[w]);
            i1->at(w)= 2.0*(u33[w]*u33[w]) + (u13[w]+u31[w])*(u13[w])+(u32[w]+u23[w])*(u23[w]);
     }
    double v11, v12, v13, v21, v22, v23, v31, v32, v33,deta_dgamma_dgamma_dtau; // helper values for entries
    deta_dgamma_dgamma_dtau =0.;

                    // loop over basis functions
                for (UN i=0; i < phi->at(0).size(); i++) 
                {
                    for (UN j=0; j < numNodes; j++)
                    {
                        // Reset values
                         v11 = 0.0;v12 = 0.0;v13=0.0; v21 = 0.0;v22 = 0.0;v23=0.0;v31=0.0;v32=0.0;v33=0.0;

                        // loop over basis functions quadrature points
                        for (UN w=0; w<phi->size(); w++) 
                        {
                            this->materialModel->evaluateFunction(this->params_,  gammaDot->at(w), viscosity_atw);

                            v11 = v11 + (0.25*deta_dgamma_dgamma_dtau)*QuadW->at(w)*(  ( 2.0*dPhiTrans[w][j][0]*a1->at(w)+dPhiTrans[w][j][1]*b1->at(w)+dPhiTrans[w][j][2]*c1->at(w))*this->surfaceElement_OutwardNormal[0] + ( dPhiTrans[w][j][1]*a1->at(w)                                         )*this->surfaceElement_OutwardNormal[1] + ( dPhiTrans[w][j][2]*a1->at(w)  )*this->surfaceElement_OutwardNormal[2] ) * ((*phi)[w][i]); // xx contribution: 
                            v12 = v12 + (0.25*deta_dgamma_dgamma_dtau)*QuadW->at(w)*(  ( dPhiTrans[w][j][0]*b1->at(w)                                  )*this->surfaceElement_OutwardNormal[0] + ( dPhiTrans[w][j][0]*a1->at(w) + 2.0*dPhiTrans[w][j][1]*b1->at(w) + dPhiTrans[w][j][2]*c1->at(w)   )*this->surfaceElement_OutwardNormal[1] + ( dPhiTrans[w][j][2]*b1->at(w)  )*this->surfaceElement_OutwardNormal[2] ) * ((*phi)[w][i]); // xy contribution:  
                            v13 = v13 + (0.25*deta_dgamma_dgamma_dtau)*QuadW->at(w)*(  ( dPhiTrans[w][j][0]*c1->at(w)                                  )*this->surfaceElement_OutwardNormal[0] + ( dPhiTrans[w][j][1]*c1->at(w)  )*this->surfaceElement_OutwardNormal[1] + ( dPhiTrans[w][j][0]*a1->at(w) + dPhiTrans[w][j][1]*b1->at(w) + 2.0*dPhiTrans[w][j][2]*c1->at(w)   )*this->surfaceElement_OutwardNormal[2] ) * ((*phi)[w][i]); // xz contribution:  

                            v21 = v21 + (0.25*deta_dgamma_dgamma_dtau)*QuadW->at(w)*(  ( 2.0*dPhiTrans[w][j][0]*d1->at(w)+dPhiTrans[w][j][1]*e1->at(w)+dPhiTrans[w][j][2]*f1->at(w))*this->surfaceElement_OutwardNormal[0] + ( dPhiTrans[w][j][1]*d1->at(w)                                         )*this->surfaceElement_OutwardNormal[1] + ( dPhiTrans[w][j][2]*d1->at(w)  )*this->surfaceElement_OutwardNormal[2] ) * ((*phi)[w][i]); // xx contribution: 
                            v22 = v22 + (0.25*deta_dgamma_dgamma_dtau)*QuadW->at(w)*(  ( dPhiTrans[w][j][0]*e1->at(w)                                  )*this->surfaceElement_OutwardNormal[0] + ( dPhiTrans[w][j][0]*d1->at(w) + 2.0*dPhiTrans[w][j][1]*e1->at(w) + dPhiTrans[w][j][2]*f1->at(w)   )*this->surfaceElement_OutwardNormal[1] + ( dPhiTrans[w][j][2]*e1->at(w)  )*this->surfaceElement_OutwardNormal[2] ) * ((*phi)[w][i]); // xy contribution:  
                            v23 = v23 + (0.25*deta_dgamma_dgamma_dtau)*QuadW->at(w)*(  ( dPhiTrans[w][j][0]*f1->at(w)                                  )*this->surfaceElement_OutwardNormal[0] + ( dPhiTrans[w][j][1]*f1->at(w)  )*this->surfaceElement_OutwardNormal[1] + ( dPhiTrans[w][j][0]*f1->at(w) + dPhiTrans[w][j][1]*e1->at(w) + 2.0*dPhiTrans[w][j][2]*f1->at(w)   )*this->surfaceElement_OutwardNormal[2] ) * ((*phi)[w][i]); // xz contribution:  

                            v31 = v31 + (0.25*deta_dgamma_dgamma_dtau)*QuadW->at(w)*(  ( 2.0*dPhiTrans[w][j][0]*g1->at(w)+dPhiTrans[w][j][1]*h1->at(w)+dPhiTrans[w][j][2]*i1->at(w))*this->surfaceElement_OutwardNormal[0] + ( dPhiTrans[w][j][1]*g1->at(w)                                         )*this->surfaceElement_OutwardNormal[1] + ( dPhiTrans[w][j][2]*g1->at(w)  )*this->surfaceElement_OutwardNormal[2] ) * ((*phi)[w][i]); // xx contribution: 
                            v32 = v32 + (0.25*deta_dgamma_dgamma_dtau)*QuadW->at(w)*(  ( dPhiTrans[w][j][0]*h1->at(w)                                  )*this->surfaceElement_OutwardNormal[0] + ( dPhiTrans[w][j][0]*g1->at(w) + 2.0*dPhiTrans[w][j][1]*h1->at(w) + dPhiTrans[w][j][2]*i1->at(w)   )*this->surfaceElement_OutwardNormal[1] + ( dPhiTrans[w][j][2]*h1->at(w)  )*this->surfaceElement_OutwardNormal[2] ) * ((*phi)[w][i]); // xy contribution:  
                            v33 = v33 + (0.25*deta_dgamma_dgamma_dtau)*QuadW->at(w)*(  ( dPhiTrans[w][j][0]*i1->at(w)                                  )*this->surfaceElement_OutwardNormal[0] + ( dPhiTrans[w][j][1]*i1->at(w)  )*this->surfaceElement_OutwardNormal[1] + ( dPhiTrans[w][j][0]*g1->at(w) + dPhiTrans[w][j][1]*h1->at(w) + 2.0*dPhiTrans[w][j][2]*i1->at(w)   )*this->surfaceElement_OutwardNormal[2] ) * ((*phi)[w][i]); // xz contribution:  
        

                        } // End loop over quadrature points

                        //multiply determinant from transformation
                        v11 *= detB_line; 
                        v12 *= detB_line;
                        v13 *= detB_line;
                        v21 *= detB_line; 
                        v22 *= detB_line;
                        v23 *= detB_line;
                        v31 *= detB_line; 
                        v32 *= detB_line;
                        v33 *= detB_line;
            
                        // Put values on the right position in element matrix - d=2 because we are in two dimensional case
                        // [v11  v12  v13]
                        // [v21  v22  v23]
                        // [v31  v32  v33]
                       (*elementMatrix)[i*dofs][j*dofs]      = v11; // d=0, first dimension
                       (*elementMatrix)[i*dofs][j*dofs+1]    = v12;
                       (*elementMatrix)[i*dofs][j*dofs+2]    = v13;
                       (*elementMatrix)[i*dofs+1][j*dofs]    = v21;
                       (*elementMatrix)[i*dofs +1][j*dofs+1] = v22; //d=1, second dimension
                       (*elementMatrix)[i*dofs +1][j*dofs+2] = v23; //d=1, second dimension
                       (*elementMatrix)[i*dofs+2][j*dofs]    = v31;
                       (*elementMatrix)[i*dofs +2][j*dofs+1] = v32; //d=2, third dimension
                       (*elementMatrix)[i*dofs +2][j*dofs+2] = v33; //d=2, third dimension
                                                

                    } // End loop over j nodes

                } // End loop over i nodes    
    }//end if 3d
    //TEUCHOS_TEST_FOR_EXCEPTION(dim == 3,std::logic_error, "AssemblyNeumannBoundaryTerm Not implemented for dim=3");
    

} // Function End loop 





// "Fixpunkt"- Matrix without jacobian for calculating Ax 
// Here update please to unlinearized System Matrix accordingly.
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::assembleRHS(){

	SmallMatrixPtr_Type elementMatrixN =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
    SmallMatrixPtr_Type elementMatrixNC =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));

	this->ANB_.reset(new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_)); // A + B + N
	this->ANB_->add( (*this->constantMatrix_),(*this->ANB_));

    // Nonlinear stress tensor *******************************
	this->assemblyStress(elementMatrixN);
	this->ANB_->add( (*elementMatrixN),(*this->ANB_));

    // Nonlinear convection term *******************************
    //this->assemblyAdvection(elementMatrixNC);
	//elementMatrixNC->scale(this->density_);
	//this->ANB_->add( (*elementMatrixNC),(*this->ANB_));

    // If boundary element - nonlinear boundar term *******************************
    if (this->surfaceElement == true)
    {
      SmallMatrixPtr_Type elementMatrixNB =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
      this->assemblyNeumannBoundaryTerm(elementMatrixNB);
      this->ANB_->add( (*elementMatrixNB),((*this->ANB_)));
    }


	this->rhsVec_ = vec_dbl_Type(this->dofsElement_,0);
	// Multiplying ANB_ * solution // System Matrix times solution
	int s=0,t=0;
	for(int i=0 ; i< this->ANB_->size();i++){
		if (i >= this->dofsElementVelocity_)
			s=1;
		for(int j=0; j < this->ANB_->size(); j++){
			if(j >= this->dofsElementVelocity_)
				t=1;
			this->rhsVec_[i] += (*this->ANB_)[i][j]*this->solution_[j]*this->coeff_[s][t];
			//cout <<"Solution["<<j <<"]" << this->solution_[i] << endl;
		}
		t=0;
		//cout <<"RHS["<<i <<"]" << this->rhsVec_[i] << endl;
	}
}




/*!

 \brief Building Transformation
@param[in] &B
*/
/* In 2D B=[ nodesRefConfig(1).at(0)-nodesRefConfig(0).at(0)    nodesRefConfig(2).at(0)-nodesRefConfig(0).at(0)    ]
           [ nodesRefConfig(1).at(1)-nodesRefConfig(0).at(1)    nodesRefConfig(2).at(1)-nodesRefConfig(0).at(1)   ]
/*!
    - Triangle numbering 

                    2 
                    * *
                    *   *
              	    5     4
                    *       *
                    *         *
                    0 * * 3 * * 1

    In 2D B=[ 1(x)-0(x)      2(x)-0(x)  ]
            [ 1(y)-0(y)    2(y)-0(y)]
*/
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::buildTransformation(SmallMatrix<SC>& B){

    TEUCHOS_TEST_FOR_EXCEPTION( (B.size()<2 || B.size()>3), std::logic_error, "Initialize SmallMatrix for transformation.");
    UN index;
    UN index0 = 0;
    for (UN j=0; j<B.size(); j++) {
        index = j+1;
        for (UN i=0; i<B.size(); i++) {
            B[i][j] = this->nodesRefConfig_.at(index).at(i) - this->nodesRefConfig_.at(index0).at(i);
        }
    }

}


template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::applyBTinv( vec3D_dbl_ptr_Type& dPhiIn,
                                    vec3D_dbl_Type& dPhiOut,
                                    SmallMatrix<SC>& Binv){
    UN dim = Binv.size();
    for (UN w=0; w<dPhiIn->size(); w++){  // quadrature point
        for (UN i=0; i < dPhiIn->at(w).size(); i++) { // basisfunction iterator 
            for (UN d1=0; d1<dim; d1++) {   // dimension of problem
                for (UN d2=0; d2<dim; d2++) { // dimension of problem
                    dPhiOut[w][i][d1] += dPhiIn->at(w).at(i).at(d2) * Binv[d2][d1];
                }
            }
        }
    }
}


// Compute Shear Rate on quadrature points depending on gradient of velocity solution at nodes
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::computeShearRate(  vec3D_dbl_Type dPhiTrans,
                                    vec_dbl_ptr_Type& gammaDot, int dim){

//****************** TWO DIMENSIONAL *********************************
    if (dim == 2)    
    {            

    vec_dbl_Type u11(dPhiTrans.size(), -1.); // should correspond to du/dx at each quadrature point
    vec_dbl_Type u12(dPhiTrans.size(), -1.); // should correspond to du/dy at each quadrature point
    vec_dbl_Type u21(dPhiTrans.size(), -1.); // should correspond to dv/dx at each quadrature point
    vec_dbl_Type u22(dPhiTrans.size(), -1.); // should correspond to dv/dy at each quadrature point

     for (UN w=0; w<dPhiTrans.size(); w++)
      { //quads points
      // set again to zero 
      u11[w] = 0.0;
      u12[w] = 0.0;
      u21[w] = 0.0;
      u22[w] = 0.0;
            for (UN i=0; i < dPhiTrans[0].size(); i++) 
            { // loop unrolling
                LO index1 = dim * i + 0; // x
                LO index2 = dim * i + 1; // y 
                // uLoc[d][w] += this->solution_[index] * phi->at(w).at(i);
                u11[w] += this->solution_[index1] * dPhiTrans[w][i][0]; // u*dphi_dx
                u12[w] += this->solution_[index1] * dPhiTrans[w][i][1]; // because we are in 2D , 0 and 1 
                u21[w] += this->solution_[index2] * dPhiTrans[w][i][0];
                u22[w] += this->solution_[index2] * dPhiTrans[w][i][1];
                
            }
            gammaDot->at(w) = sqrt(2.0*u11[w]*u11[w]+ 2.0*u22[w]*u22[w] + (u12[w]+u21[w])*(u12[w]+u21[w]));
      }
    }// end if dim == 2
    //****************** THREE DIMENSIONAL *********************************
    else if (dim == 3)
    {
    
    vec_dbl_Type u11(dPhiTrans.size(), -1.); // should correspond to du/dx at each quadrature point
    vec_dbl_Type u12(dPhiTrans.size(), -1.); // should correspond to du/dy at each quadrature point
    vec_dbl_Type u13(dPhiTrans.size(), -1.); // should correspond to du/dz at each quadrature point
    
    vec_dbl_Type u21(dPhiTrans.size(), -1.); // should correspond to dv/dx at each quadrature point
    vec_dbl_Type u22(dPhiTrans.size(), -1.); // should correspond to dv/dy at each quadrature point
    vec_dbl_Type u23(dPhiTrans.size(), -1.); // should correspond to dv/dz at each quadrature point

    vec_dbl_Type u31(dPhiTrans.size(), -1.); // should correspond to dw/dx at each quadrature point
    vec_dbl_Type u32(dPhiTrans.size(), -1.); // should correspond to dw/dy at each quadrature point
    vec_dbl_Type u33(dPhiTrans.size(), -1.); // should correspond to dw/dz at each quadrature point

    for (UN w=0; w<dPhiTrans.size(); w++)
     { //quads points
      // set again to zero 
      u11[w] = 0.0;
      u12[w] = 0.0;
      u13[w] = 0.0;
      u21[w] = 0.0;
      u22[w] = 0.0;
      u23[w] = 0.0;
      u31[w] = 0.0;
      u32[w] = 0.0;
      u33[w] = 0.0;
      
            for (UN i=0; i < dPhiTrans[0].size(); i++) 
            {
                LO index1 = dim * i + 0; //x
                LO index2 = dim * i + 1; //y 
                LO index3 = dim * i + 2; //z
               // uLoc[d][w] += this->solution_[index] * phi->at(w).at(i);
                u11[w] += this->solution_[index1] * dPhiTrans[w][i][0]; // u*dphi_dx
                u12[w] += this->solution_[index1] * dPhiTrans[w][i][1]; // because we are in 3D , 0 and 1, 2 
                u13[w] += this->solution_[index1] * dPhiTrans[w][i][2]; 
                u21[w] += this->solution_[index2] * dPhiTrans[w][i][0]; // v*dphi_dx
                u22[w] += this->solution_[index2] * dPhiTrans[w][i][1];
                u23[w] += this->solution_[index2] * dPhiTrans[w][i][2];
                u31[w] += this->solution_[index3] * dPhiTrans[w][i][0]; // w*dphi_dx
                u32[w] += this->solution_[index3] * dPhiTrans[w][i][1];
                u33[w] += this->solution_[index3] * dPhiTrans[w][i][2];
                
            }
            gammaDot->at(w) = sqrt(2.0*u11[w]*u11[w]+ 2.0*u22[w]*u22[w] + 2.0*u33[w]*u33[w] +  (u12[w]+u21[w])*(u12[w]+u21[w])   + (u13[w]+u31[w])*(u13[w]+u31[w]) + (u23[w]+u32[w])*(u23[w]+u32[w]) );
     }
    } // end if dim == 3
    
}



// So based on the previous solution we can compute viscosity in the center of mass 
// therefore we have to compute center of mass of triangle 
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::computeLocalViscosity()
{
	int dim = this->getDim();
	string FEType = this->FETypeVelocity_;

    SC detB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    buildTransformation(B);
    detB = B.computeInverse(Binv); 

    vec3D_dbl_ptr_Type 	dPhiAtCM;
    vec_dbl_Type 	    CM(dim,0.0); // center of mass - so we want to compute the viscosity value in the middle of the element
   
    // Compute center of mass **********************************************************************************
    TEUCHOS_TEST_FOR_EXCEPTION(dim == 1,std::logic_error, "computeLocalViscosity Not implemented for dim=1");
    if (dim == 2) // center of mass of reference triangle (xi-eta coordinates)
    {
     CM[0]=1.0/3.0;
     CM[1]=1.0/3.0;
    }
    else if(dim == 3) // center of mass of reference tetrahedor (xi - eta - omega coordinates)
    {
     CM[0]=1.0/4.0;
     CM[1]=1.0/4.0;
     CM[2]=1.0/4.0;
    }

    Helper::getDPhiAtCM(dPhiAtCM, dim, FEType, CM); // these are the original coordinates of the reference element
    vec3D_dbl_Type dPhiTransAtCM( dPhiAtCM->size(), vec2D_dbl_Type( dPhiAtCM->at(0).size(), vec_dbl_Type(dim,0.) ) );
    applyBTinv( dPhiAtCM, dPhiTransAtCM, Binv );    // we need transformation because of velocity gradient in shear rate equation

    vec_dbl_ptr_Type gammaDoti(new vec_dbl_Type( dPhiAtCM->size(),0.0)); // Only one value because size is one
    computeShearRate(dPhiTransAtCM, gammaDoti, dim); // updates gammaDot using velcoity solution 
    this->materialModel->evaluateFunction(this->params_,  gammaDoti->at(0), this->solutionViscosity_.at(0));
    

}

}
#endif

