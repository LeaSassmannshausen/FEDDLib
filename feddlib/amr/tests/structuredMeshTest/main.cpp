#ifndef MAIN_TIMER_START
#define MAIN_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Main") + std::string(S))));
#endif

#ifndef MAIN_TIMER_STOP
#define MAIN_TIMER_STOP(A) A.reset();
#endif

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FiniteElement.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/Laplace.hpp"
#include "feddlib/problems/abstract/Problem.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include "feddlib/core/LinearAlgebra/BlockMultiVector.hpp"
#include "feddlib/core/Mesh/Mesh.hpp"
#include "feddlib/core/Mesh/MeshInterface.hpp"
#include "feddlib/core/Mesh/MeshFileReader.hpp"
#include "feddlib/amr/AdaptiveMeshRefinement.hpp"

#include <boost/function.hpp>
#include <chrono> 


/*!
 main of Laplace problem

 @brief Laplace Adaptive
 @author Lea Saßmannshausen
 @version 1.0
 @copyright CH
 */

void oneBC(double* x, double* res, double t, const double* parameters){
    res[0] = 1.;
}
void zeroBC(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
}
void oneFunc(double* x, double* res, double* parameters){
    res[0] = 1.;
}
void dummyFunc(double* x, double* res, double t, const double* parameters){

    return;
}

void dummyFuncExactSol(double* x, double* res){

    return;
}
// #################################################################################

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::REDUCE_MAX;
using Teuchos::outArg;

using namespace FEDD;

int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	typedef Teuchos::RCP<const MultiVector_Type> MultiVectorPtrConst_Type;
	typedef Matrix<SC,LO,GO,NO> Matrix_Type;
	typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;
	typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;
	typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;
	typedef EdgeElements EdgeElements_Type;
	typedef Teuchos::RCP<EdgeElements_Type> EdgeElementsPtr_Type;
	typedef Mesh<SC,LO,GO,NO> MAIN_Type;
	typedef Teuchos::RCP<MAIN_Type > MeshPtr_Type;
	typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
	typedef Teuchos::RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
	typedef typename MAIN_Type::Elements_Type Elements_Type;
	typedef typename MAIN_Type::ElementsPtr_Type ElementsPtr_Type;
    typedef Map<LO,GO,NO> Map_Type;
    typedef typename Map_Type::MapPtr_Type MapPtr_Type;
	typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef Teuchos::RCP<Problem_Type> ProblemPtr_Type;
	typedef boost::function<void(double* x, double* res, double t, const double* parameters)>   BCFunc_Type;  

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");

    std::string vectorLaplace = "false";
    myCLP.setOption("vectorLaplace",&vectorLaplace,"vectorLaplace");
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");
    double length = 1.0;
    myCLP.setOption("length",&length,"length of domain.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    bool vL ( !vectorLaplace.compare("true") );
    
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListSolver);

        // Mesh
        int dim = parameterListProblem->sublist("Parameter").get("Dimension",2);
        int m = parameterListProblem->sublist("Parameter").get("H/h",5);
        std::string FEType = parameterListProblem->sublist("Parameter").get("Discretization","P1");
        std::string meshType = parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        std::string meshName = parameterListProblem->sublist("Parameter").get("Mesh Name","");
        std::string meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");

        int n;
        int size = comm->getSize();
        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        size -= numProcsCoarseSolve;

        int minNumberSubdomains;
        if (!meshType.compare("structured") || !meshType.compare("unstructured_struct")) {
            minNumberSubdomains = 1;
        }
        else if(!meshType.compare("structured_bfs") || !meshType.compare("unstructured_bfs")){
            minNumberSubdomains = (int) 2*length+1;
        }

		Teuchos::RCP<Domain<SC,LO,GO,NO> > domainP1;
        Teuchos::RCP<Domain<SC,LO,GO,NO> > domainP2;
        domainP1.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );

	
    	Teuchos::RCP<Domain<SC,LO,GO,NO> > domain;
        if (!meshType.compare("structured")) {
            TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured mesh.");
            if (dim == 2) {
                n = (int) (std::pow(size,1/2.) + 100.*Teuchos::ScalarTraits<double>::eps()); // 1/H
                std::vector<double> x(2);
                x[0]=0.0;    x[1]=0.0;
                domainP1 = Teuchos::rcp( new Domain<SC,LO,GO,NO>(x, 1., 1., comm) ) ;
                domainP1->buildMesh(1, "Square", dim, "P1", n, m, numProcsCoarseSolve,true);
            }
            else if (dim == 3){
                n = (int) (std::pow(size,1/3.) + 100.*Teuchos::ScalarTraits< SC >::eps()); // 1/H
                std::vector<double> x(3);
                x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
                domainP1 = Teuchos::rcp( new Domain<SC,LO,GO,NO>(x, 1., 1., 1., comm) ) ;
                domainP1->buildMesh(1, "Square", dim, "P1", n, m, numProcsCoarseSolve,true);
            }
        }
		else
		    TEUCHOS_TEST_FOR_EXCEPTION( true , std::logic_error, "Only Structured Meshes!!");

		domainP1->exportNodeFlags();
	
		// MeshRefinement Parameters

	  	MeshPartitioner_Type::DomainPtrArray_Type domainP1RefinedArray(1);
		domainP1RefinedArray[0] = domainP1;
		Teuchos::RCP<Domain<SC,LO,GO,NO> > domainRefined;
				
		Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );


		int maxIter = parameterListProblem->sublist("Mesh Refinement").get("MaxIter",3);

		string modellProblem = parameterListProblem->sublist("Mesh Refinement").get("Modell Problem","Seminar1");	

		AdaptiveMeshRefinement<SC,LO,GO,NO> meshRefiner("Laplace",parameterListProblem,dummyFuncExactSol); // exactSolLShape
		std::vector<std::chrono::duration<double>> meshTiming(maxIter);

		Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Refine Mesh"));

		auto startTotal = std::chrono::high_resolution_clock::now();

		int maxRank = std::get<1>(domainP1->getMesh()->rankRange_);
		int j=0;
		MAIN_TIMER_START(Total," Step 4:	 Total RefinementAlgorithm");
		while(j<maxIter+1 ){

			MAIN_TIMER_START(buildP2," Step 0:	 buildP2Mesh");
			if (FEType=="P2" ) {
				domainP1->setUnstructuredMesh(domainP1->getMesh());
				domainP2.reset( new Domain<SC,LO,GO,NO>( comm, dim ));
				domainP2->buildP2ofP1Domain( domainP1 );
				domain = domainP2;
				domainP2->exportNodeFlags("domainP2");
		    }
			else 
				domain = domainP1; 	
			MAIN_TIMER_STOP(buildP2);		

			MAIN_TIMER_START(Bounds," Step 1:	 bcFactory");


			bcFactory.reset( new BCBuilder<SC,LO,GO,NO>( ) );
			bcFactory->addBC(zeroBC, 1, 0, domain, "Dirichlet", 1);
            bcFactory->addBC(zeroBC, 2, 0, domain, "Dirichlet", 1);
            bcFactory->addBC(zeroBC, 3, 0, domain, "Dirichlet", 1);
			MAIN_TIMER_STOP(Bounds);	
			MAIN_TIMER_START(Solver," Step 2:	 solving PDE");

		   
			Teuchos::RCP<Laplace<SC,LO,GO,NO> > laplace(new Laplace<SC,LO,GO,NO>( domain,FEType,parameterListAll,vL));
			{
				laplace->addBoundaries(bcFactory);
            	laplace->addRhsFunction(oneFunc);
				laplace->initializeProblem();
				laplace->assemble();
				laplace->setBoundaries();
				laplace->solve();
			}
			MAIN_TIMER_STOP(Solver);	
	
			MAIN_TIMER_START(Refinement," Step 3:	 meshRefinement");

			// Refinement
			domainRefined.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
			{
				ProblemPtr_Type problem = Teuchos::rcp_dynamic_cast<Problem_Type>( laplace , true);
				domainRefined = meshRefiner.globalAlgorithm( domainP1,  domain, laplace->getSolution(), problem, oneFunc );
			}

			domainP1 = domainRefined;
			domain = domainP1;
			
			j++;
			MAIN_TIMER_STOP(Refinement);				
			
		}	

		MAIN_TIMER_STOP(Total);	
  		Teuchos::TimeMonitor::report(cout,"Main");


	 
	}
 return(EXIT_SUCCESS);

}






