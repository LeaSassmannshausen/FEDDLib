#ifndef Helper_hpp
#define Helper_hpp

//#include "AssembleFE_decl.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/Map.hpp"

namespace FEDD {
class Helper {
  
public:
    enum VarType {Std=0,Grad=1};

    /* Everything related to basisfunctions: */
    typedef default_sc SC;
    typedef default_lo LO;
    typedef default_go GO;
    typedef default_no NO;

    typedef SmallMatrix<SC> SmallMatrix_Type;
    typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;
    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;


    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;

	static void gradPhi(	int Dimension,
                    int intFE,
                    int i,
                    vec_dbl_Type &QuadPts,
                    vec_dbl_ptr_Type &value);
    
    /*! Most of the quadrature formulas can be found in http://code-aster.org/doc/v11/en/man_r/r3/r3.01.01.pdf 01/2021  */
    static void getQuadratureValues(int Dimension,
                            int Degree,
                            vec2D_dbl_ptr_Type &QuadPts,
                            vec_dbl_ptr_Type &QuadW,
                            std::string FEType);
                            
    static vec2D_dbl_Type getQuadratureValuesOnSurface(int dim, 	
    										std::string FEType, 
    										vec_dbl_Type &QuadW, 
    										vec_LO_Type surfaceIDs, 
    										vec2D_dbl_ptr_Type points);
    

    static int getDPhi(	vec3D_dbl_ptr_Type &DPhi,
                	vec_dbl_ptr_Type &weightsDPhi,
                    int Dimension,
                    std::string FEType,
                    int Degree);

    static UN determineDegree(UN dim,
                       std::string FEType,
                       UN degFunc);
                       
    static UN determineDegree(UN dim, std::string FEType, VarType type);


	static UN determineDegree(UN dim, 
								std::string FEType1, 		
								std::string FEType2, 
								VarType type1,
								VarType type2, 
								UN extraDeg = 0);

    static int getPhi(vec2D_dbl_ptr_Type &Phi,
                            vec_dbl_ptr_Type &weightsPhi,
                            int dim,
                            std::string FEType,
                            int Degree,
               			    std::string FETypeQuadPoints="");

	static void phi(int dim,
			  int intFE,
			  int i,
			  vec_dbl_Type &p,
			  double* value);



    static void applyBTinv(vec3D_dbl_ptr_Type& dPhiIn,
		            vec3D_dbl_Type& dPhiOut,
		            SmallMatrix<SC>& Binv);          

    static void buildTransformation(SmallMatrix<SC>& B, vec2D_dbl_Type& nodesRefConfig);
    // 
    static void buildTransformation(const vec_int_Type& element, vec2D_dbl_ptr_Type pointsRep,SmallMatrix<SC>& B,std::string FEType);
    //
    static void buildTransformation(const vec_int_Type& element,
                                          vec2D_dbl_ptr_Type pointsRep,
                                          SmallMatrix<SC>& B,
                                          vec_dbl_Type& b,
                                          std::string FEType);

    static void buildTransformationSurface(const vec_int_Type& element,
                                                 vec2D_dbl_ptr_Type pointsRep,
                                                 SmallMatrix<SC>& B,
                                                 vec_dbl_Type& b,
                                                 std::string FEType);

    static void fillMatrixArray(SmallMatrix<double> &matIn, double* matArrayOut, std::string order,int offset);

    static void buildFullDPhi(vec3D_dbl_ptr_Type dPhi, Teuchos::Array<SmallMatrix<double> >& dPhiMat);

    static void assemblyEmptyMatrix(MatrixPtr_Type &A);
    static void assemblyIdentity(MatrixPtr_Type &A);


    //static int checkFE(int dim,std::string FEType);

    static vec2D_dbl_Type getCoordinates(vec_LO_Type localIDs, vec2D_dbl_ptr_Type points);


private:
	
	Helper(){};

};
}
#endif
