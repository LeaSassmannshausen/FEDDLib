#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/SCI.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

void rhsDummy2D(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    return;
}

void rhsDummy(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void zeroBC(double* x, double* res, double t, const double* parameters)
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


void zeroDirichlet(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void reactionFunc(double* x, double* res, double* parameters){
	
    double m = 0.0;	
    res[0] = m * x[0];

}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void inflowChem(double* x, double* res, double t, const double* parameters)
{
    res[0] = 1.;
    
    return;
}


void rhsX(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    
    res[0] = parameters[1];
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void rhsY(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] =  parameters[1];
    res[2] = 0.;
    return;
}

void rhsZ(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = parameters[1];
    return;
}

void rhsYZ(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    double force = parameters[1];
    double TRamp = 2.0;
    if(parameters[0] <= TRamp+1e-06)
        force = parameters[0] * force * 1./(TRamp);
    else
        force = parameters[1];

    if(parameters[2] == 5)
        res[1] = force;
    else
        res[1] =0.;
        
    if (parameters[2] == 4)
        res[2] = force;
    else
        res[2] = 0.;
    
    return;
}

// parameter[0] should indicate time
void rhsImp(double* x, double* res, double* parameters){
    
	double r = sqrt(pow(x[0],2)+pow(x[1],2));

	if(parameters[0] <= 2.0 ){
   		 res[0] = (x[0]/r)*parameters[1]*sin(M_PI *1./4*(parameters[0]));
  	 	 res[1] = (x[1]/r)*parameters[1]*sin(M_PI *1./4*(parameters[0]));
	}
	else{ 
	  	res[0] = (x[0]/r)*parameters[1];
		res[1] = (x[1]/r)*parameters[1];
	}

  
    res[2] = 0.;
    return;
}

// parameter[0] should indicate time
void rhsImpTime(double* x, double* res, double* parameters){
    
	double r = sqrt(pow(x[0],2)+pow(x[1],2));
	double T_Ramp = 0.1;
	double a =2.;
	double t= parameters[0];
	double forceS = parameters[1];
	
	if(parameters[2] == 2){
		res[0] =  1./2*cos(1.0*M_PI*x[2]-4.0*t)*(x[0]/r)*forceS +0.5*(x[0]/r)*forceS; //sin(2*M_PI*(x[2]))*(x[0]/r)*forceS*sin(M_PI *1./(2*T_Ramp)*(t));
		res[1] =  1./2*cos(1.0*M_PI*x[2]-4.0*t)*(x[1]/r)*forceS +0.5*(x[1]/r)*forceS;; //sin(2*M_PI*(x[2]))*(x[1]/r)*forceS*sin(M_PI *1./(2*T_Ramp)*(t));	
	}
	else{
		res[0] = 0.0;
		res[1] = 0.0;
		
	}


	/*	if((x[2] < 0.5 + a*(parameters[0]-T_Ramp)/2.) && (x[2] > 0.0+a*(parameters[0]-T_Ramp)/2.)){
			res[0] = sin(a*M_PI*(x[2]-2*(parameters[0]-T_Ramp)/2.))*(x[0]/r)*parameters[1];
			res[1] = sin(a*M_PI*(x[2]-2*(parameters[0]-T_Ramp)/2.))*(x[1]/r)*parameters[1];
		}	*/
		
		
	
	res[2] = 0.0;

    return;
}

void dummyFunc(double* x, double* res, double t, const double* parameters)
{
    return;
}


typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;

using namespace FEDD;
using namespace Teuchos;
using namespace std;

int main(int argc, char *argv[])
{


    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef ExporterParaView<SC,LO,GO,NO> ExporterPV_Type;
    typedef RCP<ExporterPV_Type> ExporterPVPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    string xmlProblemFile = "parametersProblemSCI.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");       
    string xmlSolverFileSCI = "parametersSolverSCI.xml"; 
    myCLP.setOption("solverfileSCI",&xmlSolverFileSCI,".xml file with Inputparameters.");
    
    string xmlPrecFileStructure = "parametersPrecStructure.xml";
    myCLP.setOption("precfileStructure",&xmlPrecFileStructure,".xml file with Inputparameters.");
    string xmlPrecFileChem = "parametersPrecChem.xml";
    myCLP.setOption("precfileChem",&xmlPrecFileChem,".xml file with Inputparameters.");

    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");


    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    bool verbose (comm->getRank() == 0);

    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
       
        ParameterListPtr_Type parameterListSolverSCI = Teuchos::getParametersFromXmlFile(xmlSolverFileSCI);

        ParameterListPtr_Type parameterListPrecStructure = Teuchos::getParametersFromXmlFile(xmlPrecFileStructure);
        ParameterListPtr_Type parameterListPrecChem = Teuchos::getParametersFromXmlFile(xmlPrecFileChem);
  
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);


        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;     
        
        parameterListAll->setParameters(*parameterListSolverSCI);
        parameterListAll->setParameters(*parameterListPrec);

        
        ParameterListPtr_Type parameterListChemAll(new Teuchos::ParameterList(*parameterListPrecChem)) ;
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Chem") );
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter") );
        parameterListChemAll->setParameters(*parameterListPrecChem);

        
        ParameterListPtr_Type parameterListStructureAll(new Teuchos::ParameterList(*parameterListPrecStructure));
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Solid") );
        parameterListStructureAll->setParameters(*parameterListPrecStructure);
        parameterListStructureAll->setParameters(*parameterListProblem);

                 
        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",2);
        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","unstructured");
        
        string      discType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");
        string preconditionerMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        int         n;

        TimePtr_Type totalTime(TimeMonitor_Type::getNewCounter("FEDD - main - Total Time"));
        TimePtr_Type buildMesh(TimeMonitor_Type::getNewCounter("FEDD - main - Build Mesh"));

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);

        int size = comm->getSize() - numProcsCoarseSolve;

        // #####################
        // Mesh bauen und wahlen
        // #####################
    
        if (verbose)
        {
            cout << "###############################################" <<endl;
            cout << "############ Starting SCI  ... ################" <<endl;
            cout << "###############################################" <<endl;
        }

        DomainPtr_Type domainP1chem;
        DomainPtr_Type domainP1struct;
        DomainPtr_Type domainP2chem;
        DomainPtr_Type domainP2struct;
        
        
        DomainPtr_Type domainChem;
        DomainPtr_Type domainStructure;
        
        std::string bcType = parameterListAll->sublist("Parameter").get("BC Type","parabolic");
        
    
        domainP1chem.reset( new Domain_Type( comm, dim ) );
        domainP1struct.reset( new Domain_Type( comm, dim ) );
        domainP2chem.reset( new Domain_Type( comm, dim ) );
        domainP2struct.reset( new Domain_Type( comm, dim ) );
                                
        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(2);
        domainP1Array[0] = domainP1chem;
        domainP1Array[1] = domainP1struct;
    
        ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );                    

        pListPartitioner->set("Build Edge List",true);
        pListPartitioner->set("Build Surface List",true);
                        
        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
        
        partitionerP1.readAndPartition();
                    
        if (!discType.compare("P2")){
            domainP2chem->buildP2ofP1Domain( domainP1chem );
            domainP2struct->buildP2ofP1Domain( domainP1struct );
            
            domainChem = domainP2chem;
            domainStructure = domainP2struct;   
        }        
        else{
            domainStructure = domainP1struct;
            domainChem = domainP1chem;
        }

       // ########################
        // Flags check
        // ########################

		Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaF(new ExporterParaView<SC,LO,GO,NO>());

		Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(domainStructure->getMapUnique()));
		vec_int_ptr_Type BCFlags = domainStructure->getBCFlagUnique();

		Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);
		for(int i=0; i< entries.size(); i++){
			entries[i] = BCFlags->at(i);
		}

		Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

		exParaF->setup("Flags", domainStructure->getMesh(), discType);

		exParaF->addVariable(exportSolutionConst, "Flags", "Scalar", 1,domainStructure->getMapUnique(), domainStructure->getMapUniqueP2());

		exParaF->save(0.0);
		
		exParaF->closeExporter();
		// #####################
		// BC Check
		// #####################
		// Checking BCs
		Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());


		exportSolution.reset(new MultiVector<SC,LO,GO,NO>(domainStructure->getMapVecFieldUnique()));
		exportSolution->putScalar(0.0);
		
	    exportSolutionConst = exportSolution;
		exPara->setup("BC COND", domainStructure->getMesh(), discType);
	    
	    exPara->addVariable(exportSolutionConst, "BC Cond", "Vector", dim, domainStructure->getMapUnique());

		vec_int_ptr_Type flags = domainStructure->getBCFlagUnique();
		vec2D_dbl_ptr_Type nodes = domainStructure->getPointsUnique();

		entries  = exportSolution->getDataNonConst(0);

		double T_Ramp = 2.;
        double dt = parameterListAll->sublist("Timestepping Parameter").get("dt",1.0);
	    double tMax = parameterListAll->sublist("Timestepping Parameter").get("Final time",1.0);
    	double forceS = parameterListProblem->sublist("Parameter").get("Volume force",10.);
		double r=0.;
		vec_dbl_Type res(3);
		double a = 2.;
		for(double t=0.; t < tMax ; t= t+dt){

			for(int i=0; i< nodes->size(); i++){

				if(flags->at(i) == 2){
					vec_dbl_Type x = nodes->at(i);
					r = sqrt(pow(x[0],2)+pow(x[1],2));


					res[0] =  1./2*cos(1.0*M_PI*x[2]-4.0*t)*(x[0]/r)*forceS +0.5*(x[0]/r)*forceS; //sin(2*M_PI*(x[2]))*(x[0]/r)*forceS*sin(M_PI *1./(2*T_Ramp)*(t));
					res[1] =  1./2*cos(1.0*M_PI*x[2]-4.0*t)*(x[1]/r)*forceS +0.5*(x[1]/r)*forceS;; //sin(2*M_PI*(x[2]))*(x[1]/r)*forceS*sin(M_PI *1./(2*T_Ramp)*(t));
						
					res[2] = 0.0;

				
					for(int d=0; d<dim ; d++)
						entries[i*dim+d] = res[d];
				}
			}
	    	exPara->save(t);

		}


        if (parameterListAll->sublist("General").get("ParaView export subdomains",false) ){
            
            if (verbose)
                std::cout << "\t### Exporting subdomains ###\n";

            typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
            typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
            typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
            typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
            typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
            // Same subdomain for solid and chemistry, as they have same domain
            {
                MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type( domainStructure->getElementMap() ) );
                MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
                vecDecomposition->putScalar(comm->getRank()+1.);
                
                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
                
                exPara->setup( "subdomains_solid", domainStructure->getMesh(), "P0" );
                
                exPara->addVariable( vecDecompositionConst, "subdomains", "Scalar", 1, domainStructure->getElementMap());
                exPara->save(0.0);
                exPara->closeExporter();
            }

        }
    
        domainChem->setReferenceConfiguration();

        vec2D_dbl_Type diffusionTensor(dim,vec_dbl_Type(3));
        double D0 = parameterListAll->sublist("Parameter Diffusion").get("D0",1.);
        for(int i=0; i<dim; i++){
            diffusionTensor[0][0] =D0;
            diffusionTensor[1][1] =D0;
            diffusionTensor[2][2] =D0;

            if(i>0){
                diffusionTensor[i][i-1] = 0;
                diffusionTensor[i-1][i] = 0;
            }
            else
                diffusionTensor[i][i+1] = 0;				
        }
 
        
        Teuchos::RCP<SmallMatrix<int>> defTS;

        defTS.reset( new SmallMatrix<int> (2) );

        // Stucture
        (*defTS)[0][0] = 1;
        // Chem
        (*defTS)[1][1] = 1;
			

        SCI<SC,LO,GO,NO> sci(domainStructure, discType,
                                domainChem, discType, diffusionTensor, reactionFunc,
                                parameterListStructureAll,
                                parameterListChemAll,
                                parameterListAll,
                                defTS);
        
        sci.info();
        
            
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) ); 
            
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryChem( new BCBuilder<SC,LO,GO,NO>( ) ); 
        
    
        // Struktur-RW
        
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryStructure( new BCBuilder<SC,LO,GO,NO>( ) );

        if(dim == 2)
        {
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only 3D Test available");                               
        }
        else if(dim == 3)
        {
            bcFactory->addBC(zeroDirichlet3D, 5, 0, domainStructure, "Dirichlet_Y_Z", dim); // donut outflow
            bcFactory->addBC(zeroDirichlet3D, 6, 0, domainStructure, "Dirichlet_X_Z", dim); // donut outflow 
			bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim); // donut outflow
	        bcFactory->addBC(zeroDirichlet3D, 4, 0, domainStructure, "Dirichlet_Z", dim); // donut outflow

             
            bcFactoryStructure->addBC(zeroDirichlet3D, 5, 0, domainStructure, "Dirichlet_Y_Z", dim); // donut outflow
            bcFactoryStructure->addBC(zeroDirichlet3D, 6, 0, domainStructure, "Dirichlet_X_Z", dim); // donut outflow
            bcFactoryStructure->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim); // donut outflow
            bcFactoryStructure->addBC(zeroDirichlet3D, 4, 0, domainStructure, "Dirichlet_Z", dim); // donut outflow           
           
        }
        
        // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
        // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
        if (!sci.problemStructure_.is_null())
            sci.problemStructure_->addBoundaries(bcFactoryStructure);
        else
            sci.problemStructureNonLin_->addBoundaries(bcFactoryStructure);
        
        // RHS dummy for structure
        if (dim==2) {
            if (!sci.problemStructure_.is_null())
                sci.problemStructure_->addRhsFunction( rhsX,0 );
            else
                sci.problemStructureNonLin_->addRhsFunction( rhsX,0 );
            
        }
        else if (dim==3) {
            
            if (!sci.problemStructure_.is_null()){
                sci.problemStructure_->addRhsFunction( rhsImpTime,0 );
                double force = parameterListAll->sublist("Parameter").get("Volume force",1.);
                sci.problemStructure_->addParemeterRhs( force );
                double degree = 0.;
                sci.problemStructure_->addParemeterRhs( degree );

            }
            else{             
                sci.problemStructureNonLin_->addRhsFunction( rhsImpTime,0 );
                double force = parameterListAll->sublist("Parameter").get("Volume force",1.);
                sci.problemStructureNonLin_->addParemeterRhs( force );
                double degree = 0.;
                sci.problemStructureNonLin_->addParemeterRhs( degree );

            }
            

        }
        if (dim==2)
        {
                TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only 3D Test available");                               
                            
        }
        else if(dim==3)
        {

            bcFactory->addBC(inflowChem, 2, 1, domainChem, "Dirichlet", 1); // inflow of Chem
            bcFactoryChem->addBC(inflowChem, 2, 0, domainChem, "Dirichlet", 1); // inflow of Chem
           // bcFactoryChem->addBC(inflowChem, 4, 1, domainChem, "Dirichlet", 1); // inflow of Chem
        }

        // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
        // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
        sci.problemChem_->addBoundaries(bcFactoryChem);
        
          
        // #####################
        // Zeitintegration
        // #####################
        sci.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

        sci.initializeProblem();
        // Matrizen assemblieren
        sci.assemble();
                    
        DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);

        // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
        // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
        daeTimeSolver.defineTimeStepping(*defTS);

        // Uebergebe das (nicht) lineare Problem
        daeTimeSolver.setProblem(sci);

        // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massematrizen auf den Zeilen, welche in
        // defTS definiert worden sind.
        daeTimeSolver.setupTimeStepping();

        daeTimeSolver.advanceInTime();
    }
    TimeMonitor_Type::report(std::cout);

    return(EXIT_SUCCESS);
}