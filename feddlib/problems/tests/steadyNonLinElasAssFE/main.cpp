#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/LinElas.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"

#include <Teuchos_StackedTimer.hpp>


void zeroDirichlet(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void zeroDirichlet2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;

    return;
}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void zeroDirichletX(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = x[1];
    res[2] = x[2];

    return;
}

void zeroDirichletY(double* x, double* res, double t, const double* parameters)
{
    res[0] = x[0];
    res[1] = 0.;
    res[2] = x[2];


    return;
}

void zeroDirichletZ(double* x, double* res, double t, const double* parameters)
{
    res[0] = x[0];
    res[1] = x[1];
    res[2] = 0.;

    return;
}

void dummyFunc(double* x, double* res, double t, const double* parameters)
{
    return;
}

void rhs2D(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = parameters[1];
    
    return;
}

void rhsY(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = parameters[1];
    res[2] = 0.;
    return;
}

void rhsX(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = parameters[1];
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void rhsYZ(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] =0.;
    res[1] =0.;
    res[2] =0.;   
    double force = parameters[1];

    if(parameters[3] == 5 || parameters[3] == 4){
        res[0] = force;
        res[1] = force;
        res[2] = force;
    }
    return;
}

void rhsInterface(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    double force = parameters[1];
    
    res[0] =0.;
    res[1] =0.;
    res[2] =0.;


    if(parameters[3] == 6){
      	res[0] = force;
        res[1] = force;
        res[2] = force;
    }
 
    return;
}


typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;
using namespace Teuchos;
using namespace std;
int main(int argc, char *argv[])
{

    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef ExporterParaView<SC,LO,GO,NO> ExporterPV_Type;
    typedef RCP<ExporterPV_Type> ExporterPVPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;

    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    // int dim = 2;
    // myCLP.setOption("dim",&dim,"dim");
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        return EXIT_SUCCESS;
    }

    Teuchos::RCP<StackedTimer> stackedTimer =  rcp(new StackedTimer("Steady Nonlinear Elasticity",true));
    TimeMonitor::setStackedTimer(stackedTimer);


    bool verbose (comm->getRank() == 0); // Print-Ausgaben nur auf rank = 0
    if (verbose) {
        cout << "###############################################################" <<endl;
        cout << "############ Starting Steady Nonlinear Elasticity ... ############" <<endl;
        cout << "###############################################################" <<endl;
    }

    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListSolver);

        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",3);
        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        string		meshName    	= parameterListProblem->sublist("Parameter").get("Mesh Name","cube_0_1.mesh");
        string		meshDelimiter   = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        int         n;
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);
        string      FEType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        int size = comm->getSize() - numProcsCoarseSolve;

        Teuchos::RCP<Teuchos::Time> totalTimeAssFE(Teuchos::TimeMonitor::getNewCounter("main: Total Time Solve AssFE"));
        Teuchos::RCP<Teuchos::Time> totalTimeFEDD(Teuchos::TimeMonitor::getNewCounter("main: Total Time Solve FEDD"));

        DomainPtr_Type domain;

        // ########################
        // P1 und P2 Gitter bauen
        // ########################
        int minNumberSubdomains=1;

        if (!meshType.compare("structured")) {
		    TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured mesh.");
            n = (int)(std::pow( size/minNumberSubdomains, 1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
            std::vector<double> x(3);
            x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
            domain.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
        
		    domain->buildMesh( 3,"Square5Element", dim, FEType, n, m, numProcsCoarseSolve);

            domain->preProcessMesh(true,true);
		}
        else if (!meshType.compare("unstructured")) {
            domain.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
            MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
            domainP1Array[0] = domain;
            
            ParameterListPtr_Type pListPartitioner = sublist( parameterListProblem, "Mesh Partitioner" );
            MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
            
            partitionerP1.readAndPartition(15);
            if (FEType=="P2") {
                Teuchos::RCP<Domain<SC,LO,GO,NO> > domainP2;
                domainP2.reset( new Domain_Type( comm, dim ) );
                domainP2->buildP2ofP1Domain( domain );
                domain = domainP2;
            }
        }
        
        // ########################
        // domain->exportNodeFlags();
        // ########################

        TEUCHOS_TEST_FOR_EXCEPTION( dim==2, std::logic_error, "Only 3D tests allowed"); 

        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );
      
        if (!meshType.compare("structured")) { // Case of Cube
            bcFactory->addBC(zeroDirichlet3D, 1, 0, domain, "Dirichlet_X", dim); // x=0
            bcFactory->addBC(zeroDirichlet3D, 2, 0, domain, "Dirichlet_Y", dim); // y=0
            bcFactory->addBC(zeroDirichlet3D, 3, 0, domain, "Dirichlet_Z", dim); // z=0
        
            bcFactory->addBC(zeroDirichlet3D, 0, 0, domain, "Dirichlet", dim);
            bcFactory->addBC(zeroDirichlet3D, 7, 0, domain, "Dirichlet_X_Y", dim); //x,y = 0
            bcFactory->addBC(zeroDirichlet3D, 8, 0, domain, "Dirichlet_Y_Z", dim); // y,z= 0
            bcFactory->addBC(zeroDirichlet3D, 9, 0, domain, "Dirichlet_X_Z", dim); // x,z = 0
        
        }
        else if (!meshType.compare("unstructured")) { // Case of Artery Mesh
            bcFactory->addBC(zeroDirichlet3D, 14, 0, domain, "Dirichlet_Y_Z", dim); // inflow/outflow strip fixed in y direction
            bcFactory->addBC(zeroDirichlet3D, 13, 0, domain, "Dirichlet_X_Z", dim); // inflow/outflow strip fixed in y direction
            bcFactory->addBC(zeroDirichlet3D, 7, 0, domain, "Dirichlet_Z", dim); // inlet fixed in Z direction
            bcFactory->addBC(zeroDirichlet3D, 8, 0, domain, "Dirichlet_Z", dim); // outlet fixed in Z direction
            bcFactory->addBC(zeroDirichlet3D, 9, 0, domain, "Dirichlet_Z", dim); // inlet ring in Z direction
            bcFactory->addBC(zeroDirichlet3D, 1, 0, domain, "Dirichlet_Z", dim); // outer ring of inlet area
            bcFactory->addBC(zeroDirichlet3D, 2, 0, domain, "Dirichlet_Z", dim); // outer ring of outlet area
            bcFactory->addBC(zeroDirichlet3D, 10, 0, domain, "Dirichlet_Z", dim); // outlet ring in Z direction
        }
        
        // LinElas Objekt erstellen
        NonLinElasticity<SC,LO,GO,NO> NonLinElasAssFE( domain, FEType, parameterListAll );

        NonLinElasAssFE.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen
    
        double force = parameterListAll->sublist("Parameter").get("Volume force",0.);
        double degree = 0;

        if (!meshType.compare("structured")) { // Case of Cube
            NonLinElasAssFE.addRhsFunction( rhsYZ );// rhsYZ
        }
        else if (!meshType.compare("unstructured")) { // Case of Artery Mesh
            NonLinElasAssFE.addRhsFunction( rhsInterface );// rhsYZ
        }

        NonLinElasAssFE.addParemeterRhs( force );
        NonLinElasAssFE.addParemeterRhs( degree );
        
        // ######################
        // Matrix assemblieren, RW setzen und System loesen
        // ######################
        NonLinElasAssFE.initializeProblem();
        NonLinElasAssFE.assemble();                
        NonLinElasAssFE.setBoundaries(); // In der Klasse Problem
        NonLinElasAssFE.setBoundariesRHS();

		std::string nlSolverType = parameterListProblem->sublist("General").get("Linearization","FixedPoint");
        NonLinearSolver<SC,LO,GO,NO> nlSolverAssFE( nlSolverType );
        
        {
            Teuchos::TimeMonitor totalTimeMonitorAssFE(*totalTimeAssFE);
            nlSolverAssFE.solve( NonLinElasAssFE );
            comm->barrier();
        }
        
        // Nonlinear Elasticity Objekt erstellen not from AceGen Interface
        parameterListAll->sublist("Parameter").set("Use AceGen Interface",false); 
        NonLinElasticity<SC,LO,GO,NO> NonLinElas( domain, FEType, parameterListAll );

        NonLinElas.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

        if (!meshType.compare("structured")) { // Case of Cube
            NonLinElas.addRhsFunction( rhsYZ );// rhsYZ
        }
        else if (!meshType.compare("unstructured")) { // Case of Artery Mesh
            NonLinElas.addRhsFunction( rhsInterface );// rhsYZ
        }
        
        
        NonLinElas.addParemeterRhs( force );
        NonLinElas.addParemeterRhs( degree );
        // ######################
        // Matrix assemblieren, RW setzen und System loesen
        // ######################
        NonLinElas.initializeProblem();
        NonLinElas.assemble();                
        NonLinElas.setBoundaries(); // In der Klasse Problem
		NonLinElas.setBoundariesRHS();

        NonLinearSolver<SC,LO,GO,NO> nlSolver( nlSolverType );
        
        {
            Teuchos::TimeMonitor totalTimeMonitorFEDD(*totalTimeFEDD);
            nlSolver.solve( NonLinElas );
            comm->barrier();	

        }

		if(comm->getRank() ==0){
			cout << " ############################################### " << endl;
			cout << " Nonlinear Iterations FEDDLib Assembly : " << nlSolver.getNonLinIts() << endl; 
			cout << " Nonlinear Iterations AceGEN Assembly  : " << nlSolverAssFE.getNonLinIts() << endl;
			cout << " ############################################### " << endl;

		}    
       
        MultiVectorConstPtr_Type valuesSolidConst1 = NonLinElas.getSolution()->getBlock(0);
        MultiVectorConstPtr_Type valuesSolidConst2 = NonLinElasAssFE.getSolution()->getBlock(0);
        // Calculating the error per node
        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > errorValues = Teuchos::rcp(new MultiVector<SC,LO,GO,NO>( valuesSolidConst1->getMap() ) ); 
        //this = alpha*A + beta*B + gamma*this
        errorValues->update( 1., valuesSolidConst2, -1. ,valuesSolidConst1, 0.);
        // Taking abs norm
        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > errorValuesAbs = errorValues;

        errorValues->abs(errorValuesAbs);

        if( parameterListProblem->sublist("General").get("ParaViewExport",false) ) {

            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

            exPara->setup( "displacements", domain->getMesh(), FEType );

            exPara->addVariable( valuesSolidConst1, "valuesNonLinElas", "Vector", dim, domain->getMapUnique());

            exPara->addVariable( valuesSolidConst2, "valuesNonLinElasAssFE", "Vector", dim, domain->getMapUnique());
        
            exPara->addVariable( errorValuesAbs, "erroeValues", "Vector", dim, domain->getMapUnique());
            exPara->save(0.0);
            
        }

        Teuchos::Array<SC> norm(1); 
        errorValues->normInf(norm);//const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms);
        double res = norm[0];
        if(comm->getRank() ==0)
            cout << " Inf Norm of Error of Solutions " << res << endl;
        double infNormError = res;
    
        NonLinElas.getSolution()->getBlock(0)->normInf(norm);
        res = norm[0];
        if(comm->getRank() ==0)
            cout << " Relative error Inf-Norm of solution nonlinear elasticity " << infNormError/res << endl;

        NonLinElasAssFE.getSolution()->getBlock(0)->normInf(norm);
        res = norm[0];
        if(comm->getRank() ==0)
            cout << " Relative error Inf-Norm of solutions nonlinear elasticity assemFE " << infNormError/res << endl;
    

        TEUCHOS_TEST_FOR_EXCEPTION( infNormError > 1e-11 , std::logic_error, "Inf Norm of Error between calculated solutions is too great. Exceeded 1e-11. ");

    }
    Teuchos::TimeMonitor::report(cout);
    stackedTimer->stop("Steady Nonlinear Elasticity");
	StackedTimer::OutputOptions options;
	options.output_fraction = options.output_histogram = options.output_minmax = true;
	stackedTimer->report((std::cout),comm,options);

    return(EXIT_SUCCESS);
}
