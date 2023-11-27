#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/FSI.hpp"
#include "feddlib/problems/specific/Laplace.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*! Test case for specific artery geometrie or straight tube geometry. Inflow depends on inflow region
	-> artery: Inflow scaled with normal vector on inflow (x,y,z) * laplaceInflow	

*/



void zeroBC(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
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
	if(t>=parameters[0])
    	res[0] = 1.;
    else	
    	res[0] = 0.;
    return;
}

void parabolicInflow3D(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0];

    }

    return;
}

void parabolicInflowDirection3D(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing

    res[0] = 0.;
    res[1] = 0.;
    res[2] = parameters[0] * x[0];

    return;
}

void parabolicInflow3DLin(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * t / parameters[1];
    }
    else
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0];
    }

    return;
}

void parabolicInflow3DArtery(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else
    {
        res[1] = 0.;
        res[0] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0];
    }

    return;
}

// Flowrate Ramp
void parabolicInflow3DFlowRateRamp(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[5] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else
    {
        res[1] = 0.;
        res[0] = 0.;
        res[2] = parameters[5];
    }

    return;
}

void parabolicInflow3DArteryHeartBeat(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    double heartBeatStart = parameters[3];

    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else if(t > heartBeatStart)
    {
    
        double a0    = 11.693284502463376;
        double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
            0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
            0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
        double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
            -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
            0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
                    
        double Q = 0.5*a0;
        

        double t_min = t - fmod(t,1.0)+heartBeatStart-std::floor(t); ; //FlowConditions::t_start_unsteady;
        double t_max = t_min + 1.0; // One heartbeat lasts 1.0 second    
        double y = M_PI * ( 2.0*( t-t_min ) / ( t_max - t_min ) -1.0)  ;
        
        for(int i=0; i< 20; i++)
            Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
        
        
        // Remove initial offset due to FFT
        Q -= 0.026039341343493;
        Q = (Q - 2.85489)/(7.96908-2.85489);

        res[0] = 0.;
        res[1] = 0.;
        res[2] = (parameters[0] / parameters[2] * (x[0] + x[0] * Q)) ;
        
    }
    else
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0];

    }

    return;
}

void flowRate3DArteryHeartBeat(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    double heartBeatStart = parameters[3];

    if(t < parameters[1])
    {
        res[0] = parameters[5] / parameters[2] * x[0] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else if(t > heartBeatStart)
    {
    
        double a0    = 11.693284502463376;
        double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
            0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
            0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
        double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
            -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
            0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
                    
        double Q = 0.5*a0;
        

        double t_min = t - fmod(t,1.0)+heartBeatStart-std::floor(t); ; //FlowConditions::t_start_unsteady;
        double t_max = t_min + 1.0; // One heartbeat lasts 1.0 second    
        double y = M_PI * ( 2.0*( t-t_min ) / ( t_max - t_min ) -1.0)  ;
        
        for(int i=0; i< 20; i++)
            Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
        
        
        // Remove initial offset due to FFT
        Q -= 0.026039341343493;
        Q = (Q - 2.85489)/(7.96908-2.85489);

        res[0] = (parameters[5] / parameters[2] * (x[0] + 1.6*x[0] * Q)) ;
        
    }
    else
    {
        res[0] = parameters[5] / parameters[2] * x[0];

    }

    return;
}
void parabolicInflowSteady(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
   
    res[0] = 0.;
    res[1] = 0.;
    res[2] = parameters[0] / parameters[2] * x[0] * 0.5 * ( ( 1 - cos( 1./2* M_PI*parameters[4]/parameters[1]) ));
    
    

    return;
}



void rhsResistance(double* x, double* res, double* parameters){

    double pressureValue = parameters[1];
    double flag = parameters[2];
    double ramp = parameters[3];


  	res[0] =0.;
    
    if(parameters[0]+1.e-12 < ramp)
        pressureValue = parameters[0]*pressureValue/ramp;
    else
        pressureValue = parameters[1];

    if(flag == 5){
      	res[0] = pressureValue;
        
    }
    
    return;
}

void parabolicInflow3DLinArtery(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * t / parameters[1];
    }
    else
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0];
    }

    return;
}

void rhsDummy(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void oneFunc(double* x, double* res, double* parameters){
    res[0] = 1.;
}

void dummyFunc(double* x, double* res, double t, const double* parameters)
{
    return;
}

void reactionFunc(double* x, double* res, double* parameters){
	
    double m = 0.0;	
    res[0] = m * x[0];

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
    typedef Mesh<SC,LO,GO,NO> Mesh_Type;
    typedef RCP<Mesh_Type> MeshPtr_Type;
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
    string xmlProblemFile = "parametersProblemFSI.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFileGE = "parametersPrecGE.xml"; // GE
    string xmlPrecFileGI = "parametersPrecGI.xml"; // GI
    myCLP.setOption("precfileGE",&xmlPrecFileGE,".xml file with Inputparameters.");
    myCLP.setOption("precfileGI",&xmlPrecFileGI,".xml file with Inputparameters.");
    string xmlSolverFileFSI = "parametersSolverFSI.xml"; // GI
    myCLP.setOption("solverfileFSI",&xmlSolverFileFSI,".xml file with Inputparameters.");
    string xmlSolverFileGeometry = "parametersSolverGeometry.xml"; // GE
    myCLP.setOption("solverfileGeometry",&xmlSolverFileGeometry,".xml file with Inputparameters.");

    string xmlPrecFileFluidMono = "parametersPrecFluidMono.xml";
    string xmlPrecFileFluidTeko = "parametersPrecFluidTeko.xml";
    myCLP.setOption("precfileFluidMono",&xmlPrecFileFluidMono,".xml file with Inputparameters.");
    myCLP.setOption("precfileFluidTeko",&xmlPrecFileFluidTeko,".xml file with Inputparameters.");
    string xmlProblemFileFluid = "parametersProblemFluid.xml";
    myCLP.setOption("problemFileFluid",&xmlProblemFileFluid,".xml file with Inputparameters.");
    string xmlProblemFileStructure = "parametersProblemStructure.xml";
    myCLP.setOption("problemFileStructure",&xmlProblemFileStructure,".xml file with Inputparameters.");
    string xmlPrecFileStructure = "parametersPrecStructure.xml";
    myCLP.setOption("precfileStructure",&xmlPrecFileStructure,".xml file with Inputparameters.");
    string xmlPrecFileGeometry = "parametersPrecGeometry.xml";
    myCLP.setOption("precfileGeometry",&xmlPrecFileGeometry,".xml file with Inputparameters.");
    
    string xmlProbL = "plistProblemLaplace.xml";
    myCLP.setOption("probLaplace",&xmlProbL,".xml file with Inputparameters.");
    string xmlPrecL = "plistPrecLaplace.xml";
    myCLP.setOption("precLaplace",&xmlPrecL,".xml file with Inputparameters.");
    string xmlSolverL = "plistSolverLaplace.xml";
    myCLP.setOption("solverLaplace",&xmlSolverL,".xml file with Inputparameters.");
    
    string xmlPrecFileChem = "parametersPrecChem.xml";
    myCLP.setOption("precfileChem",&xmlPrecFileChem,".xml file with Inputparameters.");
    
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
        ParameterListPtr_Type parameterListStructure = Teuchos::getParametersFromXmlFile(xmlProblemFileStructure);

        ParameterListPtr_Type parameterListSolverFSI = Teuchos::getParametersFromXmlFile(xmlSolverFileFSI);
        ParameterListPtr_Type parameterListSolverGeometry = Teuchos::getParametersFromXmlFile(xmlSolverFileGeometry);
        ParameterListPtr_Type parameterListPrecGeometry = Teuchos::getParametersFromXmlFile(xmlPrecFileGeometry);

        ParameterListPtr_Type parameterListPrecGE = Teuchos::getParametersFromXmlFile(xmlPrecFileGE);
        ParameterListPtr_Type parameterListPrecGI = Teuchos::getParametersFromXmlFile(xmlPrecFileGI);
        ParameterListPtr_Type parameterListPrecFluidMono = Teuchos::getParametersFromXmlFile(xmlPrecFileFluidMono);
        ParameterListPtr_Type parameterListPrecFluidTeko = Teuchos::getParametersFromXmlFile(xmlPrecFileFluidTeko);

        ParameterListPtr_Type parameterListPrecStructure = Teuchos::getParametersFromXmlFile(xmlPrecFileStructure);
        
        ParameterListPtr_Type parameterListPrecChem = Teuchos::getParametersFromXmlFile(xmlPrecFileChem);
        
        bool geometryExplicit = parameterListProblem->sublist("Parameter").get("Geometry Explicit",true);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        if(geometryExplicit)
            parameterListAll->setParameters(*parameterListPrecGE);
        else
            parameterListAll->setParameters(*parameterListPrecGI);
        
        parameterListAll->setParameters(*parameterListSolverFSI);

        
        ParameterListPtr_Type parameterListFluidAll(new Teuchos::ParameterList(*parameterListPrecFluidMono)) ;
        sublist(parameterListFluidAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Fluid") );
        parameterListFluidAll->setParameters(*parameterListPrecFluidTeko);
        parameterListFluidAll->setParameters(*parameterListPrecFluidMono);
        parameterListFluidAll->setParameters(*parameterListSolverFSI);



        
        ParameterListPtr_Type parameterListStructureAll(new Teuchos::ParameterList(*parameterListPrecStructure));
        sublist(parameterListStructureAll, "Parameter Solid")->setParameters( parameterListStructure->sublist("Parameter Solid") );
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Solid") );
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter") );

        //parameterListStructureAll->setParameters(*parameterListStructure);
        parameterListStructureAll->setParameters(*parameterListPrecStructure);

        ParameterListPtr_Type parameterListChemAll(new Teuchos::ParameterList(*parameterListPrecChem));
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Chem") );
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter") );

        ParameterListPtr_Type parameterListSCIAll(new Teuchos::ParameterList(*parameterListPrecStructure));
        parameterListSCIAll->setParameters(*parameterListProblem);
        //sublist(parameterListSCIAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Solid") );
        //sublist(parameterListSCIAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter") );

        
        
        // Fuer das Geometrieproblem, falls GE
        // CH: We might want to add a paramterlist, which defines the Geometry problem
        ParameterListPtr_Type parameterListGeometry(new Teuchos::ParameterList(*parameterListPrecGeometry));
        parameterListGeometry->setParameters(*parameterListSolverGeometry);
        // we only compute the preconditioner for the geometry problem once
        sublist( parameterListGeometry, "General" )->set( "Preconditioner Method", "MonolithicConstPrec" );
        sublist( parameterListGeometry, "Parameter" )->set( "Model", parameterListProblem->sublist("Parameter").get("Model Geometry","Laplace") );
        
        double poissonRatio = parameterListProblem->sublist("Parameter Geometry").get("Poisson Ratio",0.3);
        double mu = parameterListProblem->sublist("Parameter Geometry").get("Mu",2.0e+6);
        double distanceLaplace = parameterListProblem->sublist("Parameter Geometry").get("Distance Laplace",0.1);
        double coefficientLaplace = parameterListProblem->sublist("Parameter Geometry").get("Coefficient Laplace",1000.);
        
        sublist( parameterListGeometry, "Parameter" )->set( "Poisson Ratio", poissonRatio );
        sublist( parameterListGeometry, "Parameter" )->set( "Mu", mu );
        sublist( parameterListGeometry, "Parameter" )->set( "Distance Laplace", distanceLaplace );
        sublist( parameterListGeometry, "Parameter" )->set( "Coefficient Laplace", coefficientLaplace );
            
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
            cout << "############ Starting fsi  ... ################" <<endl;
            cout << "###############################################" <<endl;
        }

        DomainPtr_Type domainP1fluid;
        DomainPtr_Type domainP1struct;
        DomainPtr_Type domainP1chem;
        DomainPtr_Type domainP2fluid;
        DomainPtr_Type domainP2struct;
        DomainPtr_Type domainP2chem;
        
        DomainPtr_Type domainFluidVelocity;
        DomainPtr_Type domainFluidPressure;
        DomainPtr_Type domainChem;
        DomainPtr_Type domainStructure;
        DomainPtr_Type domainGeometry;
        
        std::string bcType = parameterListAll->sublist("Parameter").get("BC Type","Compute Inflow");
        std::string geometryType = parameterListAll->sublist("Parameter").get("Geometry Type","Artery");
        
        
        TimeMonitor_Type totalTimeMonitor(*totalTime);
    
        TimeMonitor_Type buildMeshMonitor(*buildMesh);
        if (verbose)
        {
            cout << " -- Building Mesh ... " << flush;
        }

        domainP1fluid.reset( new Domain_Type( comm, dim ) );
        domainP1struct.reset( new Domain_Type( comm, dim ) );
        domainP1chem.reset(new Domain_Type(comm,dim));
        
        domainP2fluid.reset( new Domain_Type( comm, dim ) );
        domainP2struct.reset( new Domain_Type( comm, dim ) );
        domainP2chem.reset( new Domain_Type(comm,dim));
        
        //                    

        vec_int_Type idsInterface(1,6);
                                
                        
        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(2);
        domainP1Array[0] = domainP1fluid;
        domainP1Array[1] = domainP1struct;
        
        ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
        if (!discType.compare("P2")){
            pListPartitioner->set("Build Edge List",true);
            pListPartitioner->set("Build Surface List",true);
        }
        else{
            pListPartitioner->set("Build Edge List",false);
            pListPartitioner->set("Build Surface List",false);
        }
        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim);
        
        partitionerP1.readAndPartition(15,"cm" , true );
        
        if (!discType.compare("P2")){
            domainP2fluid->buildP2ofP1Domain( domainP1fluid );
            domainP2struct->buildP2ofP1Domain( domainP1struct );
        }
        
        
        

        // Calculate distances is done in: identifyInterfaceParallelAndDistance


        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(domainP2fluid->getMapUnique()));
        vec_int_ptr_Type BCFlags = domainP2fluid->getBCFlagUnique();

        Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);
        for(int i=0; i< entries.size(); i++){
            entries[i] = BCFlags->at(i);
        }

        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

        exPara->setup("FlagsFluid",domainP2fluid->getMesh(), discType);

        exPara->addVariable(exportSolutionConst, "Flags", "Scalar", 1,domainP2fluid->getMapUnique(), domainP2fluid->getMapUniqueP2());

        exPara->save(0.0);

        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara2(new ExporterParaView<SC,LO,GO,NO>());

        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution2(new MultiVector<SC,LO,GO,NO>(domainP2struct->getMapUnique()));
        vec_int_ptr_Type BCFlags2 = domainP2struct->getBCFlagUnique();

        Teuchos::ArrayRCP< SC > entries2  = exportSolution2->getDataNonConst(0);
        for(int i=0; i< entries2.size(); i++){
            entries2[i] = BCFlags2->at(i);
        }

        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst2 = exportSolution2;

        exPara2->setup("FlagsStructure", domainP2struct->getMesh(), discType);

        exPara2->addVariable(exportSolutionConst2, "Flags", "Scalar", 1,domainP2struct->getMapUnique(), domainP2struct->getMapUniqueP2());

        exPara2->save(0.0);

  		 domainP1fluid->identifyInterfaceParallelAndDistance(domainP1struct, idsInterface);
        if (!discType.compare("P2"))
            domainP2fluid->identifyInterfaceParallelAndDistance(domainP2struct, idsInterface);
        
        
        if (!discType.compare("P2"))
        {
            domainFluidVelocity = domainP2fluid;
            domainFluidPressure = domainP1fluid;
            domainStructure = domainP2struct;
            domainGeometry = domainP2fluid;
        }
        else
        {
            domainFluidVelocity = domainP1fluid;
            domainFluidPressure = domainP1fluid;
            domainStructure = domainP1struct;
            domainGeometry = domainP1fluid;
            //                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"P1/P1 for FSI not implemented!");
        }

        if (verbose){
            cout << "done! -- " << endl;
        }
            
                        
        // Baue die Interface-Maps in der Interface-Nummerierung
        domainFluidVelocity->buildInterfaceMaps();
        
        domainStructure->buildInterfaceMaps();

        // domainInterface als dummyDomain mit mapVecFieldRepeated_ als interfaceMapVecFieldUnique_.
        // Wird fuer den Vorkonditionierer und Export gebraucht.
        // mesh is needed for rankRanges
        DomainPtr_Type domainInterface;
        domainInterface.reset( new Domain_Type( comm ) );
        domainInterface->setDummyInterfaceDomain(domainFluidVelocity);

        domainFluidVelocity->setReferenceConfiguration();
        domainFluidPressure->setReferenceConfiguration();
        

        
        /*if (parameterListAll->sublist("General").get("ParaView export subdomains",false) ){
            
            if (verbose)
                std::cout << "\t### Exporting fluid and solid subdomains ###\n";

            typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
            typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
            typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
            typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
            typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

            {
                MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type( domainFluidVelocity->getElementMap() ) );
                MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
                vecDecomposition->putScalar(comm->getRank()+1.);
                
                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
                
                exPara->setup( "subdomains_fluid", domainFluidVelocity->getMesh(), "P0" );
                
                exPara->addVariable( vecDecompositionConst, "subdomains", "Scalar", 1, domainFluidVelocity->getElementMap());
                exPara->save(0.0);
                exPara->closeExporter();
            }
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
        }*/



            // #####################
        // Problem definieren
        // #####################
        Teuchos::RCP<SmallMatrix<int>> defTS;
        if(geometryExplicit)
        {
            // SmallMatrix<int> defTS(4);
            defTS.reset( new SmallMatrix<int> (4) );

            // Fluid
            (*defTS)[0][0] = 1;
            (*defTS)[0][1] = 1;

            // Struktur
            (*defTS)[2][2] = 1;
            
 
        }
        else
        {
            // SmallMatrix<int> defTS(5);
            defTS.reset( new SmallMatrix<int> (5) );

            // Fluid
            (*defTS)[0][0] = 1;
            (*defTS)[0][1] = 1;
            // TODO: [0][4] und [1][4] bei GI + Newton noetig?
            if (verbose)
                std::cout << "### Double check temporal discretization of Shape Derivatives! ###" << std::endl;
            
            (*defTS)[0][4] = 1;
            (*defTS)[1][4] = 1;
            
        
            
        }

       
        FSI<SC,LO,GO,NO> fsi(domainFluidVelocity, discType,
                        domainFluidPressure, "P1",
                        domainStructure, discType,
                        domainInterface, discType,
                        domainGeometry, discType,
                        parameterListFluidAll,
                        parameterListStructureAll,
                        parameterListAll,
                        parameterListGeometry,
                        defTS);


        domainFluidVelocity->info();
        domainFluidPressure->info();
        domainStructure->info();
        domainGeometry->info();

        fsi.info();
                    
        std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter").get("Max Velocity",1.));
        parameter_vec.push_back( parameterListProblem->sublist("Parameter").get("Max Ramp Time",2.) );
        
        std::vector<double> parameter_vec_pressure(1, parameterListProblem->sublist("Parameter").get("Max Ramp Time",1.));
        parameter_vec_pressure.push_back( parameterListProblem->sublist("Parameter").get("Scale Pressure",1.) );

       
        TEUCHOS_TEST_FOR_EXCEPTION(bcType != "Compute Inflow", std::logic_error, "Select a valid boundary condition. Only Compute Inflow available.");

        //#############################################
        //#############################################
        //#### Compute parabolic inflow with laplacian
        //#############################################
        //#############################################
        MultiVectorPtr_Type solutionLaplace;
        MultiVectorConstPtr_Type solutionLaplaceConst;
        {
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryLaplace(new BCBuilder<SC,LO,GO,NO>( ));
            
            bcFactoryLaplace->addBC(zeroBC, 2, 0, domainFluidVelocity, "Dirichlet", 1); //inflow ring
            bcFactoryLaplace->addBC(zeroBC, 3, 0, domainFluidVelocity, "Dirichlet", 1); //outflow ring
            bcFactoryLaplace->addBC(zeroBC, 6, 0, domainFluidVelocity, "Dirichlet", 1); //surface
            
            ParameterListPtr_Type parameterListProblemL = Teuchos::getParametersFromXmlFile(xmlProbL);
            ParameterListPtr_Type parameterListPrecL = Teuchos::getParametersFromXmlFile(xmlPrecL);
            ParameterListPtr_Type parameterListSolverL = Teuchos::getParametersFromXmlFile(xmlSolverL);

            ParameterListPtr_Type parameterListLaplace(new Teuchos::ParameterList(*parameterListProblemL)) ;
            parameterListLaplace->setParameters(*parameterListPrecL);
            parameterListLaplace->setParameters(*parameterListSolverL);
            
            Laplace<SC,LO,GO,NO> laplace( domainFluidVelocity, discType, parameterListLaplace, false );
            {
                laplace.addRhsFunction(oneFunc);
                laplace.addBoundaries(bcFactoryLaplace);
                
                laplace.initializeProblem();
                laplace.assemble();
                laplace.setBoundaries();
                laplace.solve();
            }
            
            //We need the values in the inflow area. Therefore, we use the above bcFactory and the volume flag 10 and the outlet flag 5 and set zero Dirichlet boundary values
            bcFactoryLaplace->addBC(zeroBC, 5, 0, domainFluidVelocity, "Dirichlet", 1);
            bcFactoryLaplace->addBC(zeroBC, 15, 0, domainFluidVelocity, "Dirichlet", 1);
            bcFactoryLaplace->setRHS( laplace.getSolution(), 0.);
            solutionLaplace = laplace.getSolution()->getBlockNonConst(0);
        
            SC maxValue = solutionLaplace->getMax();
            solutionLaplace->scale(1./maxValue); // normalizing solution
            
            solutionLaplaceConst = solutionLaplace; 
 
            parameter_vec.push_back(1.0); // We scaled the solution beforehand, so we dont need the actual maxValue any more and replace it with 1.
            parameter_vec.push_back( parameterListProblem->sublist("Parameter").get("Heart Beat Start",2.) );
            parameter_vec.push_back(parameterListProblem->sublist("Timestepping Parameter").get("dt",0.001));
            parameter_vec.push_back(parameterListProblem->sublist("Parameter").get("Flowrate",3.0));

            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
            
            exPara->setup("parabolicInflow", domainFluidVelocity->getMesh(), discType);
            
//                exPara->setup(domainFluidVelocity->getDimension(), domainFluidVelocity->getNumElementsGlobal(), domainFluidVelocity->getElements(), domainFluidVelocity->getPointsUnique(), domainFluidVelocity->getMapUnique(), domainFluidVelocity->getMapRepeated(), discType, "parabolicInflow", 1, comm);

            MultiVectorConstPtr_Type valuesConst = laplace.getSolution()->getBlock(0);
            exPara->addVariable( valuesConst, "values", "Scalar", 1, domainFluidVelocity->getMapUnique() );

            exPara->save(0.0);
            exPara->closeExporter();

        }
        
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

        // TODO: Vermutlich braucht man keine bcFactoryFluid und bcFactoryStructure,
        // da die RW sowieso auf dem FSI-Problem gesetzt werden.

        // Fluid-RW
        {
            bool zeroPressure = parameterListProblem->sublist("Parameter Fluid").get("Set Outflow Pressure to Zero",false);
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluid( new BCBuilder<SC,LO,GO,NO>( ) );

                            
            bcFactory->addBC(flowRate3DArteryHeartBeat, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplaceConst, true , parabolicInflowDirection3D); // inflow 
            bcFactoryFluid->addBC(flowRate3DArteryHeartBeat, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplaceConst, true,parabolicInflowDirection3D); // inflow 

            bcFactory->addBC(zeroDirichlet3D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec);// solutionLaplaceConst, true , parabolicInflowDirection3D); // inflow 
            bcFactoryFluid->addBC(zeroDirichlet3D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec);// solutionLaplaceConst, true , parabolicInflowDirection3D); // inflow 

            //bcFactory->addBC(zeroDirichlet3D, 2, 0, domainFluidVelocity, "Dirichlet", dim); // inflow ring                
            //bcFactoryFluid->addBC(zeroDirichlet3D, 2, 0, domainFluidVelocity, "Dirichlet", dim); // inflow ring
            //bcFactorySteadyFluid->addBC(zeroDirichlet3D, 2, 0, domainFluidVelocity, "Dirichlet", dim); // inflow ring

            //bcFactory->addBC(zeroDirichlet3D, 3, 0, domainFluidVelocity, "Dirichlet", dim); // outflow ring                
            //bcFactoryFluid->addBC(zeroDirichlet3D, 3, 0, domainFluidVelocity, "Dirichlet", dim); // outflow ring
            //bcFactorySteadyFluid->addBC(zeroDirichlet3D, 6, 0, domainFluidVelocity, "Dirichlet", dim); // Interface
            //bcFactorySteadyFluid->addBC(zeroDirichlet3D, 3, 0, domainFluidVelocity, "Dirichlet", dim); // Outflow ring

            //bcFactory->addBC(pressureBC, 5, 1, domainFluidPressure, "Neumann", 1,parameter_vec_pressure); // outflow                
            //bcFactoryFluid->addBC(pressureBC, 5, 1, domainFluidPressure, "Neumann", 1,parameter_vec_pressure); // outflow
            
            if (zeroPressure) {
                bcFactory->addBC(zeroBC, 3, 1, domainFluidPressure, "Dirichlet", 1); // outflow ring
                bcFactory->addBC(zeroBC, 5, 1, domainFluidPressure, "Dirichlet", 1); // outflow
                
                bcFactoryFluid->addBC(zeroBC, 3, 1, domainFluidPressure, "Dirichlet", 1); // outflow ring
                bcFactoryFluid->addBC(zeroBC, 5, 1, domainFluidPressure, "Dirichlet", 1); // outflow
            }
            
            // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
            // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
            fsi.problemFluid_->addBoundaries(bcFactoryFluid);

        }

        // Struktur-RW
        {
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryStructure( new BCBuilder<SC,LO,GO,NO>( ) );
            bcFactory->addBC(zeroDirichlet3D, 14, 2, domainStructure, "Dirichlet_Y_Z", dim); // inflow/outflow strip fixed in y direction
            bcFactory->addBC(zeroDirichlet3D, 13, 2, domainStructure, "Dirichlet_X_Z", dim); // inflow/outflow strip fixed in y direction
            bcFactory->addBC(zeroDirichlet3D, 7, 2, domainStructure, "Dirichlet_Z", dim); // inlet fixed in Z direction
            bcFactory->addBC(zeroDirichlet3D, 8, 2, domainStructure, "Dirichlet_Z", dim); // outlet fixed in Z direction
            bcFactory->addBC(zeroDirichlet3D, 9, 2, domainStructure, "Dirichlet", dim); // inlet ring in Z direction
            bcFactory->addBC(zeroDirichlet3D, 10, 2, domainStructure, "Dirichlet_Z", dim); // outlet ring in Z direction

            bcFactoryStructure->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Y_Z", dim); 
            bcFactoryStructure->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_X_Z", dim); 
            bcFactoryStructure->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_Z", dim);           
            bcFactoryStructure->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Z", dim); 
            bcFactoryStructure->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Z", dim);           
            bcFactoryStructure->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Z", dim); 
        // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
            // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
        		if (!fsi.problemStructure_.is_null())
                    fsi.problemStructure_->addBoundaries(bcFactoryStructure);
                else
                    fsi.problemStructureNonLin_->addBoundaries(bcFactoryStructure);
            }
            // RHS dummy for structure
  
        if (!fsi.problemStructure_.is_null())
            fsi.problemStructure_->addRhsFunction( rhsDummy );
        else
        fsi.problemStructureNonLin_->addRhsFunction( rhsDummy );
    

        fsi.problemFluid_->addRhsFunction(rhsResistance,0);
        double resistance= parameterListAll->sublist("Parameter Fluid").get("Resistance",0.5);

        fsi.problemFluid_->addParemeterRhs( resistance);
        fsi.problemFluid_->addParemeterRhs( parameterListProblem->sublist("Parameter Fluid").get("Resistance Ramp",0.1));

        // Geometrie-RW separat, falls geometrisch explizit.
        // Bei Geometrisch implizit: Keine RW in die factoryFSI fuer das
        // Geometrie-Teilproblem, da sonst (wg. dem ZeroDirichlet auf dem Interface,
        // was wir brauchen wegen Kopplung der Struktur) der Kopplungsblock C4
        // in derselben Zeile, der nur Werte auf dem Interface haelt, mit eliminiert.
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryGeometry( new BCBuilder<SC,LO,GO,NO>( ) );
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluidInterface;
        if (preconditionerMethod == "FaCSI") // || preconditionerMethod == "FaCSI-Teko")
            bcFactoryFluidInterface = Teuchos::rcp( new BCBuilder<SC,LO,GO,NO>( ) );


        //bcFactoryGeometry->addBC(zeroDirichlet3D, 2, 0, domainGeometry, "Dirichlet", dim); // inlet fixed in Z direction
        //bcFactoryGeometry->addBC(zeroDirichlet3D, 3, 0, domainGeometry, "Dirichlet", dim); // inlet fixed in X direction
        bcFactoryGeometry->addBC(zeroDirichlet3D, 4, 0, domainGeometry, "Dirichlet", dim); // Inlet
        bcFactoryGeometry->addBC(zeroDirichlet3D, 5, 0, domainGeometry, "Dirichlet", dim); // Outlet
        bcFactoryGeometry->addBC(zeroDirichlet3D, 6, 0, domainGeometry, "Dirichlet", dim); // Interface
		/* bcFactoryGeometry->addBC(zeroDirichlet3D, 7, 0, domainGeometry, "Dirichlet", dim); // ?
		bcFactoryGeometry->addBC(zeroDirichlet3D, 8, 0, domainGeometry, "Dirichlet", dim); // ?
		bcFactoryGeometry->addBC(zeroDirichlet3D, 9, 0, domainGeometry, "Dirichlet", dim); // ?
		bcFactoryGeometry->addBC(zeroDirichlet3D, 10, 0, domainGeometry, "Dirichlet", dim); // ? */
        
        // Die RW, welche nicht Null sind in der rechten Seite (nur Interface) setzen wir spaeter per Hand.
        // Hier erstmal Dirichlet Nullrand, wird spaeter von der Sturkturloesung vorgegeben
            bcFactoryGeometry->addBC(zeroDirichlet3D, 6, 0, domainGeometry, "Dirichlet", dim); // interface
        if (preconditionerMethod == "FaCSI" ) //|| preconditionerMethod == "FaCSI-Teko")
                bcFactoryFluidInterface->addBC(zeroDirichlet3D, 6, 0, domainFluidVelocity, "Dirichlet", dim);

        fsi.problemGeometry_->addBoundaries(bcFactoryGeometry);
        if ( preconditionerMethod == "FaCSI")// || preconditionerMethod == "FaCSI-Teko")
            fsi.getPreconditioner()->setFaCSIBCFactory( bcFactoryFluidInterface );

    // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
    // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()


     
        // Matrizen assemblieren
        /*if(parameterListAll->sublist("General").get("Use steady fluid solution",false) == true){
            cout << " Solve Steady State Navier-Stokes " << endl;
            // Defining steady state Navier Stokes problem.
            //this->problemSteadyFluid_->addBoundaries(this->bcFactory_);
            fsi.problemSteadyFluid_->initializeProblem();

            fsi.problemSteadyFluid_->assemble();
            fsi.problemSteadyFluid_->setBoundariesRHS();

            // Solving the problem
            std::string nlSolverType = "NOX";
            NonLinearSolver<SC,LO,GO,NO> nlSolver( nlSolverType );
            nlSolver.solve( *(fsi.problemSteadyFluid_) );

            // Using velocity and pressure as start solution
            fsi.problemFluid_->getSolution()->getBlockNonConst(0)->update(1.0, *fsi.problemSteadyFluid_->getSolution()->getBlockNonConst(0), 1.);
            fsi.problemFluid_->getSolution()->getBlockNonConst(1)->update(1.0, *fsi.problemSteadyFluid_->getSolution()->getBlockNonConst(1), 1.);

            ExporterPVPtr_Type exporterSteadyFluid = Teuchos::rcp(new ExporterPV_Type());
            ExporterPVPtr_Type exporterSteadyPressure = Teuchos::rcp(new ExporterPV_Type());
            
            MeshPtr_Type meshNonConstF = Teuchos::rcp_const_cast<Mesh_Type>( domainFluidVelocity->getMesh() );
            MeshPtr_Type meshNonConstP = Teuchos::rcp_const_cast<Mesh_Type>( domainFluidPressure->getMesh() );

            exporterSteadyFluid->setup("u_f_steady", meshNonConstF, domainFluidVelocity->getFEType(), parameterListAll);
            exporterSteadyPressure->setup("p_steady", meshNonConstP, domainFluidPressure->getFEType(), parameterListAll);

            MultiVectorConstPtr_Type u_f_steady = fsi.problemSteadyFluid_->getSolution()->getBlock(0);            
            MultiVectorConstPtr_Type p_steady = fsi.problemSteadyFluid_->getSolution()->getBlock(1);

            exporterSteadyFluid->addVariable( u_f_steady, "u", "Vector", 3, domainFluidVelocity->getMapUnique() );
            exporterSteadyPressure->addVariable( p_steady, "p", "Scalar", 1, domainFluidPressure->getMapUnique() );

            exporterSteadyFluid->save( 0. );
            exporterSteadyPressure->save(0. );

            exporterSteadyFluid->closeExporter();
            exporterSteadyPressure->closeExporter();

        }*/
           // #####################
        // Zeitintegration
        // #####################
        fsi.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

        fsi.initializeProblem();
        
        fsi.initializeGE();

        fsi.assemble();
    
        DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);

        // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
        // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
        daeTimeSolver.defineTimeStepping(*defTS);

        // Uebergebe das (nicht) lineare Problem
        daeTimeSolver.setProblem(fsi);

        // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massematrizen auf den Zeilen, welche in
        // defTS definiert worden sind.
        daeTimeSolver.setupTimeStepping();

        daeTimeSolver.advanceInTime();
    }
    

    TimeMonitor_Type::report(std::cout);

    return(EXIT_SUCCESS);
}
