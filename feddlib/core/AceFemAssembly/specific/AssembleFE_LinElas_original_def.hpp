#ifndef AssembleFE_LinElas_ORIGINAL_DEF_hpp
#define AssembleFE_LinElas_ORIGINAL_DEF_hpp

#include "AssembleFE_LinElas_original_decl.hpp"
#include "feddlib/core/AceFemAssembly/AceInterface/NeoHookQuadraticTets.hpp"
#include <vector>
#include <iostream>

namespace FEDD {


/*!

 \brief Constructor for AssembleFE_Laplace

@param[in] flag Flag of element
@param[in] nodesRefConfig Nodes of element in reference configuration
@param[in] params Parameterlist for current problem
@param[in} tuple Vector of tuples with Discretization information

*/
template <class SC, class LO, class GO, class NO>
AssembleFE_LinElas_original<SC,LO,GO,NO>::AssembleFE_LinElas_original(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params,tuple)
{
	/// Extracting values from ParameterList params:
	E_ = this->params_->sublist("Parameter").get("E",3500.0); // the last value is the dafault value, in case no parameter is set
    poissonRatio_ = this->params_->sublist("Parameter").get("Poisson Ratio",0.4e-0);
    mu_ = this->params_->sublist("Parameter").get("mu",384.6153846153146);
    lambda_ = (poissonRatio_*E_)/((1 + poissonRatio_)*(1 - 2*poissonRatio_));

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element, int:Dofs)
	FEType_ = std::get<1>(this->diskTuple_->at(0)); // FEType of Disk
	dofs_ = std::get<2>(this->diskTuple_->at(0)); // Degrees of freedom per node
	numNodes_ = std::get<3>(this->diskTuple_->at(0)); // Number of nodes of element

	dofsElement_ = dofs_*numNodes_; // "Dimension of return matrix"

}

/*!

 \brief Assembly Jacobian is simply assemblyLaplacian for Laplace Problem

@param[in] &elementMatrix

*/ 

template <class SC, class LO, class GO, class NO>
void AssembleFE_LinElas_original<SC,LO,GO,NO>::assembleJacobian() {


	SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( dofsElement_)); // Matrix we fill with entries.

	assemblyLinElas_original(elementMatrix); // Function that fills the matrix. We pass though a pointer that will be filled.

	this->jacobian_ = elementMatrix ; // We init the jacobian matrix with the matrix we just build.
}

/*!

 \brief Assembly function 

@param[in] &elementMatrix

*/
template <class SC, class LO, class GO, class NO>
void AssembleFE_LinElas_original<SC,LO,GO,NO>::assemblyLinElas_original(SmallMatrixPtr_Type &elementMatrix) {

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
  
	// hat man nodesRefConfig_?
    Helper::buildTransformation(B, this->nodesRefConfig_);
    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);
   
	// dPhiTrans sind die transformierten Basifunktionen, also B^(-T) * \grad_phi bzw. \grad_phi^T * B^(-1).
    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    Helper::applyBTinv( dPhi, dPhiTrans, Binv );
    // Fuer Zwischenergebnisse
    SC res;

    // Fuer die Berechnung der Spur
    double res_trace_i, res_trace_j;    
    
    if (dim == 2)
    {

        double v11, v12, v21, v22;
  
        // Matrizen der Groesse (2x2) in denen die einzelnen Epsilon-Tensoren berechnet werden.
        // Siehe unten fuer mehr.
        SmallMatrix<double> epsilonValuesMat1_i(dim), epsilonValuesMat2_i(dim),
        epsilonValuesMat1_j(dim), epsilonValuesMat2_j(dim);

        for (int i = 0; i < dPhi->at(0).size(); i++)
        {
                Teuchos::Array<SC> value11( 1, 0. );
                Teuchos::Array<SC> value12( 1, 0. );
                Teuchos::Array<SC> value21( 1, 0. );
                Teuchos::Array<SC> value22( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhi->at(0).size(); j++)
                {
                    v11 = 0.0; v12 = 0.0; v21 = 0.0; v22 = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {
                        // In epsilonValuesMat1_i (2x2 Matrix) steht fuer die Ansatzfunktion i bzw. \phi_i
                        // der epsilonTensor fuer eine skalare Ansatzfunktion fuer die Richtung 1 (vgl. Mat1).
                        // Also in Mat1_i wird dann also phi_i = (phi_scalar_i, 0) gesetzt und davon \eps berechnet.

                        // Stelle \hat{grad_phi_i} = basisValues_i auf, also B^(-T)*grad_phi_i
                        // GradPhiOnRef( dPhi->at(k).at(i), b_T_inv, basisValues_i );

                        // \eps(v) = \eps(phi_i)
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat1_i, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat2_i, 1); // y-Richtung

                        // Siehe oben, nur fuer j
                        // GradPhiOnRef( DPhi->at(k).at(j), b_T_inv, basisValues_j );

                        // \eps(u) = \eps(phi_j)
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat1_j, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat2_j, 1); // y-Richtung

                        // Nun berechnen wir \eps(u):\eps(v) = \eps(phi_j):\eps(phi_i).
                        // Das Ergebniss steht in res.
                        // Berechne zudem noch die Spur der Epsilon-Tensoren tr(\eps(u)) (j) und tr(\eps(v)) (i)
                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat1_j, res); // x-x
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v11 = v11 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat2_j, res); // x-y
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v12 = v12 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat1_j, res); // y-x
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v21 = v21 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat2_j, res); // y-y
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v22 = v22 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);


                    }
                    // Noch mit der abs(det(B)) skalieren
                    v11 = absDetB * v11;
                    v12 = absDetB * v12;
                    v21 = absDetB * v21;
                    v22 = absDetB * v22;

          			  // Put values on the right position in element matrix - d=2 because we are in two dimensional case
          			  // [v11  v12  ]
          			  // [v21  v22  ]
          		   (*elementMatrix)[i*dofs][j*dofs]   = v11; // d=0, first dimension
           		   (*elementMatrix)[i*dofs][j*dofs+1] = v12;
           		   (*elementMatrix)[i*dofs+1][j*dofs] = v21;
          		   (*elementMatrix)[i*dofs +1][j*dofs+1] =v22; //d=1, second dimension
                }
        }
    }
    else if(dim == 3)
    {

        double v11, v12, v13, v21, v22, v23, v31, v32, v33;

        SmallMatrix<double> epsilonValuesMat1_i(dim), epsilonValuesMat2_i(dim), epsilonValuesMat3_i(dim),
        epsilonValuesMat1_j(dim), epsilonValuesMat2_j(dim), epsilonValuesMat3_j(dim);


            for (int i = 0; i < dPhi->at(0).size(); i++)
            {
                Teuchos::Array<SC> value11( 1, 0. );
                Teuchos::Array<SC> value12( 1, 0. );
                Teuchos::Array<SC> value13( 1, 0. );
                Teuchos::Array<SC> value21( 1, 0. );
                Teuchos::Array<SC> value22( 1, 0. );
                Teuchos::Array<SC> value23( 1, 0. );
                Teuchos::Array<SC> value31( 1, 0. );
                Teuchos::Array<SC> value32( 1, 0. );
                Teuchos::Array<SC> value33( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhi->at(0).size(); j++)
                {
                    v11 = 0.0; v12 = 0.0; v13 = 0.0; v21 = 0.0; v22 = 0.0; v23 = 0.0; v31 = 0.0; v32 = 0.0; v33 = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {

                        // GradPhiOnRef( DPhi->at(k).at(i), b_T_inv, basisValues_i );

                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat1_i, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat2_i, 1); // y-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat3_i, 2); // z-Richtung


                        // GradPhiOnRef( DPhi->at(k).at(j), b_T_inv, basisValues_j );

                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat1_j, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat2_j, 1); // y-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat3_j, 2); // z-Richtung

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat1_j, res); // x-x
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v11 = v11 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat2_j, res); // x-y
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v12 = v12 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat3_j, res); // x-z
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat3_j.trace(res_trace_j);
                        v13 = v13 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat1_j, res); // y-x
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v21 = v21 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat2_j, res); // y-y
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v22 = v22 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat3_j, res); // y-z
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat3_j.trace(res_trace_j);
                        v23 = v23 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                        epsilonValuesMat3_i.innerProduct(epsilonValuesMat1_j, res); // z-x
                        epsilonValuesMat3_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v31 = v31 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                        epsilonValuesMat3_i.innerProduct(epsilonValuesMat2_j, res); // z-y
                        epsilonValuesMat3_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v32 = v32 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                        epsilonValuesMat3_i.innerProduct(epsilonValuesMat3_j, res); // z-z
                        epsilonValuesMat3_i.trace(res_trace_i);
                        epsilonValuesMat3_j.trace(res_trace_j);
                        v33 = v33 + weights->at(k)*(2*mu_*res + lambda_*res_trace_j*res_trace_i);

                    }
                    v11 = absDetB * v11;
                    v12 = absDetB * v12;
                    v13 = absDetB * v13;
                    v21 = absDetB * v21;
                    v22 = absDetB * v22;
                    v23 = absDetB * v23;
                    v31 = absDetB * v31;
                    v32 = absDetB * v32;
                    v33 = absDetB * v33;

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
                }
            }
        
    }
}


template <class SC, class LO, class GO, class NO>
void AssembleFE_LinElas_original<SC,LO,GO,NO>::epsilonTensor(vec_dbl_Type &basisValues, SmallMatrix<SC> &epsilonValues, int activeDof){

    for (int i=0; i<epsilonValues.size(); i++) {
        for (int j=0; j<epsilonValues.size(); j++) {
            epsilonValues[i][j] = 0.;
            if (i==activeDof) {
                epsilonValues[i][j] += 0.5*basisValues.at(j);
            }
            if (j==activeDof) {
                epsilonValues[i][j] += 0.5*basisValues.at(i);
            }
        }
    }
}

}
#endif

