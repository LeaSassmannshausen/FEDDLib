#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include <algorithm>
#include <cmath>
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/General/HDF5Import.hpp"

#include "feddlib/problems/specific/FSI.hpp"
#include "feddlib/problems/specific/Laplace.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"
#include "feddlib/problems/specific/LinElas.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include "feddlib/problems/specific/Geometry.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/Solver/Preconditioner.hpp"


/*! Boundary-condition benchmark based on fsi_artery.
    Geometry and boundary labels can be adapted independently in this directory.
    See README.md for the inherited mesh, unit, and inlet-profile assumptions.
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

// BCBuilder passes the imported scalar profile value in x[0].
// The flag-5 inlets lie at the high-z end and point into the fluid along -z.
void inflowParabolicArtery(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = -parameters[0] * x[0];
}

void flowrateParabolicArtery(double* x, double* res, double t, const double* parameters)
{
    double ramp = 1.;
    if (parameters[1] > 0.) {
        const double phase = std::max(0., std::min(1., t/parameters[1]));
        ramp = 0.5 * (1. - std::cos(M_PI * phase));
    }
    res[0] = parameters[2] * ramp;
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

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string xmlProblemFile = "parametersProblemFSI.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlSolverFileFSI = "parametersSolverFSI.xml"; // GI
    myCLP.setOption("solverfileFSI",&xmlSolverFileFSI,".xml file with Inputparameters.");
    string xmlSolverFileGeometry = "parametersSolverGeometry.xml"; // GE
    myCLP.setOption("solverfileGeometry",&xmlSolverFileGeometry,".xml file with Inputparameters.");

    string xmlPrecFileFluidMono = "parametersPrecFluidMono.xml";
    string xmlPrecFileFluidTeko = "parametersPrecFluidTeko.xml";
    string xmlPrecFileFluidBlock = "parametersPrecFluidBlock.xml";

    myCLP.setOption("precfileFluidMono",&xmlPrecFileFluidMono,".xml file with Inputparameters.");
    myCLP.setOption("precfileFluidTeko",&xmlPrecFileFluidTeko,".xml file with Inputparameters.");
    // string xmlProblemFileFluid = "parametersProblemFluid.xml";
    // myCLP.setOption("problemFileFluid",&xmlProblemFileFluid,".xml file with Inputparameters.");
    string xmlPrecFileStructure = "parametersPrecStructure.xml";
    myCLP.setOption("precfileStructure",&xmlPrecFileStructure,".xml file with Inputparameters.");
    string xmlPrecFileGeometry = "parametersPrecGeometry.xml";
    myCLP.setOption("precfileGeometry",&xmlPrecFileGeometry,".xml file with Inputparameters.");
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        return EXIT_SUCCESS;
    }

    bool verbose (comm->getRank() == 0);

    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListSolverFSI = Teuchos::getParametersFromXmlFile(xmlSolverFileFSI);
        ParameterListPtr_Type parameterListSolverGeometry = Teuchos::getParametersFromXmlFile(xmlSolverFileGeometry);
        ParameterListPtr_Type parameterListPrecGeometry = Teuchos::getParametersFromXmlFile(xmlPrecFileGeometry);

        ParameterListPtr_Type parameterListPrecFluidMono = Teuchos::getParametersFromXmlFile(xmlPrecFileFluidMono);
        ParameterListPtr_Type parameterListPrecFluidTeko = Teuchos::getParametersFromXmlFile(xmlPrecFileFluidTeko);
        ParameterListPtr_Type parameterListPrecFluidBlock = Teuchos::getParametersFromXmlFile(xmlPrecFileFluidBlock);



        ParameterListPtr_Type parameterListPrecStructure = Teuchos::getParametersFromXmlFile(xmlPrecFileStructure);
        
        bool geometryExplicit = parameterListProblem->sublist("Parameter").get("Geometry Explicit",true);

        // Insert benchmark defaults before copying settings to the FSI problem.
        parameterListProblem->sublist("General").get("Flag Inlet Fluid", 5);
        parameterListProblem->sublist("General").get("Flag Outlet Fluid", 4);
        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        parameterListAll->setParameters(*parameterListSolverFSI);

        std::string preconditionerMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        ParameterListPtr_Type parametersListPrecFluid;
        
         // We also have the option to use FaCSI with one of the FEDDLib implemented block preconditioners
        if(preconditionerMethod == "FaCSI")
            parametersListPrecFluid = parameterListPrecFluidMono;
        else if(preconditionerMethod == "FaCSI-Teko")
            parametersListPrecFluid=parameterListPrecFluidTeko;
        else if(preconditionerMethod == "FaCSI-Block")
            parametersListPrecFluid=parameterListPrecFluidBlock;

        ParameterListPtr_Type parameterListFluidAll(new Teuchos::ParameterList(*parametersListPrecFluid)) ;
        sublist(parameterListFluidAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Fluid") );
        std::string precTypeFluid = parameterListProblem->sublist("Parameter Fluid").get("Preconditioner Type","Monolithic");

        sublist( parameterListFluidAll, "General" )->set( "Preconditioner Method",precTypeFluid  );
        // Information used for PCD
        sublist( parameterListFluidAll, "General" )->set( "Flag Inlet Fluid",parameterListProblem->sublist("General").get("Flag Inlet Fluid",5) );
        sublist( parameterListFluidAll, "General" )->set( "Flag Outlet Fluid",parameterListProblem->sublist("General").get("Flag Outlet Fluid",4)  );
        sublist( parameterListFluidAll, "General" )->set( "Flag Interface",parameterListProblem->sublist("General").get("Flag Interface",6)  );
        sublist( parameterListFluidAll, "Timestepping Parameter" )->setParameters( parameterListProblem->sublist("Timestepping Parameter") );

        // parameterListFluidAll->setParameters(*parameterListPrecFluidTeko);

        
        ParameterListPtr_Type parameterListStructureAll(new Teuchos::ParameterList(*parameterListPrecStructure));
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Solid") );
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter") );

        parameterListStructureAll->setParameters(*parameterListPrecStructure);


        
        // Fuer das Geometrieproblem, falls GE
        // CH: We might want to add a paramterlist, which defines the Geometry problem
        ParameterListPtr_Type parameterListGeometry(new Teuchos::ParameterList(*parameterListPrecGeometry));
        parameterListGeometry->setParameters(*parameterListSolverGeometry);
        sublist(parameterListGeometry, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Geometry") );

        // we only compute the preconditioner for the geometry problem once
        sublist( parameterListGeometry, "General" )->set( "Preconditioner Method", "MonolithicConstPrec" );
             
        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",3);        
        std::string      discType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");
        int         n;

        TimePtr_Type totalTime(TimeMonitor_Type::getNewCounter("FEDD - main - Total Time"));
        TimePtr_Type buildMesh(TimeMonitor_Type::getNewCounter("FEDD - main - Build Mesh"));

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);

        int size = comm->getSize() - numProcsCoarseSolve;

        // #####################
        // Mesh bauen und wahlen
        // #####################
        {
            if (verbose)
            {
                cout << "###############################################" <<endl;
                cout << "### Starting fsi_artery_boundary_tests ... #####" <<endl;
                cout << "###############################################" <<endl;
            }

            DomainPtr_Type domainP1fluid;
            DomainPtr_Type domainP1struct;
            DomainPtr_Type domainP2fluid;
            DomainPtr_Type domainP2struct;
                    
            DomainPtr_Type domainFluidVelocity;
            DomainPtr_Type domainFluidPressure;
            DomainPtr_Type domainStructure;
            DomainPtr_Type domainGeometry;
            
               
            {
                TimeMonitor_Type totalTimeMonitor(*totalTime);
                {
                    TimeMonitor_Type buildMeshMonitor(*buildMesh);
                    if (verbose)
                    {
                        cout << " -- Building Mesh ... " << flush;
                    }

                    domainP1fluid.reset( new Domain_Type( comm, dim ) );
                    domainP1struct.reset( new Domain_Type( comm, dim ) );
                    domainP2fluid.reset( new Domain_Type( comm, dim ) );
                    domainP2struct.reset( new Domain_Type( comm, dim ) );
                    //                    

                    vec_int_Type idsInterface(1,6);
                    idsInterface.push_back(9);
                    idsInterface.push_back(10);

                    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(2);
                    domainP1Array[0] = domainP1fluid;
                    domainP1Array[1] = domainP1struct;
                    
                    ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
                    if (!discType.compare("P2")){
                        pListPartitioner->set("Build Edge List",true);
                        pListPartitioner->set("Build Surface List",true);
                    }
                    else{
                        pListPartitioner->set("Build Edge List",true);
                        pListPartitioner->set("Build Surface List",true);
                    }
                    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
                    
                    partitionerP1.readAndPartition(15,"mm",true); // We convert the mesh of this test from mm to cm!! 
                    
                    if (!discType.compare("P2")){
                        domainP2fluid->buildP2ofP1Domain( domainP1fluid );
                        domainP2struct->buildP2ofP1Domain( domainP1struct );
                    }
                    
                    
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

                    domainFluidVelocity->exportNodeFlags("Fluid");
                    domainStructure->exportNodeFlags("Solid");

                    // Calculate distances is done in: identifyInterfaceParallelAndDistance
                    domainP1fluid->identifyInterfaceParallelAndDistance(domainP1struct, idsInterface);
                    if (!discType.compare("P2"))
                        domainP2fluid->identifyInterfaceParallelAndDistance(domainP2struct, idsInterface);
                    
                    if (verbose){
                        cout << "done! -- " << endl;
                    }
                }
            }

           
             if (parameterListAll->sublist("General").get("ParaView export subdomains",false) ){
                
                if (verbose)
                    std::cout << "\t### Exporting fluid and solid subdomains ###\n";

               domainFluidVelocity->exportDistribution("Fluid");
               domainStructure->exportDistribution("Solid");

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
                
                // Struktur
                (*defTS)[2][2] = 1;
            }

            FSI<SC,LO,GO,NO> fsi(domainFluidVelocity, discType,
                                 domainFluidPressure, "P1",
                                 domainStructure, discType,
                                 domainInterface, discType,
                                 domainGeometry, discType,
                                 parameterListFluidAll, parameterListStructureAll, parameterListAll,
                                 parameterListGeometry, defTS);


            domainFluidVelocity->info();
            domainFluidPressure->info();
            domainStructure->info();
            domainGeometry->info();
            fsi.info();
                     
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

            // TODO: Vermutlich braucht man keine bcFactoryFluid und bcFactoryStructure,
            // da die RW sowieso auf dem FSI-Problem gesetzt werden.

            // #####################
            // Randwerte
            // Fluid: 5 = inflow, 4 = outflow, 10 = inflow ring, 9 = outflow ring
            // Struktur: 7 = linke (z=0) Seite, 8 = rechte (z=L) seite. 13,14 einzelne Freiheitsgrade festgehalten in x,y Richtung
            // Interface: 6 , 9 , 10  
            // #####################

            // Fluid-RW
            {
                auto& inflow = parameterListProblem->sublist("Parameter Fluid");
                const std::string bcType = inflow.get("BC Type", "parabolic_artery");
                TEUCHOS_TEST_FOR_EXCEPTION(bcType != "parabolic_artery" || dim != 3,
                    std::logic_error, "This benchmark requires the 3D parabolic_artery inlet.");
                std::vector<double> parameter_vec;
                parameter_vec = {1., inflow.get("Inflow Ramp Time", 0.), inflow.get("Flowrate", 3.)};
                for (const double value : parameter_vec)
                    TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(value), std::logic_error,
                        "parabolic_artery parameters must be finite.");
                TEUCHOS_TEST_FOR_EXCEPTION(parameter_vec[2] < 0. || parameter_vec[1] < 0., std::logic_error,
                    "parabolic_artery requires nonnegative flow rate/ramp time.");
                // The wall is the moving FSI interface; its velocity is set by
                // coupling to the solid, rather than a zero-velocity wall BC.
                const int inletFlag = inflow.get("Inflow Flag", 5);
                TEUCHOS_TEST_FOR_EXCEPTION(inletFlag != parameterListProblem->sublist("General").get("Flag Inlet Fluid", 5),
                    std::logic_error, "Inflow Flag must match General/Flag Inlet Fluid.");
                const bool zeroPressure = inflow.get("Set Outflow Pressure to Zero", false);
                Teuchos::RCP<BCBuilder<SC,LO,GO,NO>> bcFactoryFluid(new BCBuilder<SC,LO,GO,NO>());
                std::string meshStem = parameterListProblem->sublist("Mesh Partitioner").get<std::string>("Mesh 1 Name");
                meshStem = meshStem.substr(meshStem.find_last_of("/\\") + 1);
                const auto extension = meshStem.find_last_of('.');
                if (extension != std::string::npos) meshStem.erase(extension);
                const std::string profileFile = inflow.get("Inflow Profile File",
                    "laplace_parabolic_" + meshStem + "_" + discType);
                HDF5Import<SC,LO,GO,NO> importer(domainFluidVelocity->getMapUnique(), profileFile);
                MultiVectorConstPtr_Type solutionLaplace = importer.readVariablesHDF5("solution");
                importer.closeImporter();
                const double profileMaximum = solutionLaplace->getMax();
                TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(profileMaximum) || profileMaximum <= 0.,
                    std::logic_error, "Imported inlet profile must have a positive finite maximum.");
                bcFactory->addBC(inflowParabolicArtery, inletFlag, 0, domainFluidVelocity, "Dirichlet", dim,
                    parameter_vec, solutionLaplace, true, flowrateParabolicArtery);
                bcFactoryFluid->addBC(inflowParabolicArtery, inletFlag, 0, domainFluidVelocity, "Dirichlet", dim,
                    parameter_vec, solutionLaplace, true, flowrateParabolicArtery);
                if (verbose)
                    std::cout << "parabolic_artery: Laplace profile " << profileFile << ".h5, inlet flag "
                              << inletFlag << ", inward direction -z, target flow rate " << parameter_vec[2]
                              << " cm^3/s" << std::endl;

                             
                // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
                // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
                fsi.problemFluid_->addBoundaries(bcFactoryFluid);
            }

            // Struktur-RW
            {
                Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryStructure( new BCBuilder<SC,LO,GO,NO>( ) );
                bcFactory->addBC(zeroDirichlet3D, 14, 2, domainStructure, "Dirichlet_Y_Z", dim); // point in y=0 plane fixed in y-z direction
                bcFactory->addBC(zeroDirichlet3D, 13, 2, domainStructure, "Dirichlet_X_Z", dim); // point in x=0 plane fixed in y direction
                bcFactory->addBC(zeroDirichlet3D, 7, 2, domainStructure, "Dirichlet_Z", dim); // outlet fixed in Z direction
                bcFactory->addBC(zeroDirichlet3D, 8, 2, domainStructure, "Dirichlet_Z", dim); // inlet fixed in Z direction
                bcFactory->addBC(zeroDirichlet3D, 9, 2, domainStructure, "Dirichlet_Z", dim); // outlet ring in Z direction
                bcFactory->addBC(zeroDirichlet3D, 10, 2, domainStructure, "Dirichlet_Z", dim); // inlet ring in Z direction
                bcFactory->addBC(zeroDirichlet3D, 11, 2, domainStructure, "Dirichlet_Z", dim); // outlet outer ring in Z direction
                bcFactory->addBC(zeroDirichlet3D, 12, 2, domainStructure, "Dirichlet_Z", dim); // inlet outer ring in Z direction

                bcFactoryStructure->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Y_Z", dim); 
                bcFactoryStructure->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_X_Z", dim); 
                bcFactoryStructure->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_Z", dim);           
                bcFactoryStructure->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Z", dim); 
                bcFactoryStructure->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Z", dim);  
                bcFactoryStructure->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Z", dim); 
                bcFactoryStructure->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet_Z", dim); 
                bcFactoryStructure->addBC(zeroDirichlet3D, 12, 0, domainStructure, "Dirichlet_Z", dim); 
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
        
            // Geometrie-RW separat, falls geometrisch explizit.
            // Bei Geometrisch implizit: Keine RW in die factoryFSI fuer das
            // Geometrie-Teilproblem, da sonst (wg. dem ZeroDirichlet auf dem Interface,
            // was wir brauchen wegen Kopplung der Struktur) der Kopplungsblock C4
            // in derselben Zeile, der nur Werte auf dem Interface haelt, mit eliminiert.
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryGeometry( new BCBuilder<SC,LO,GO,NO>( ) );
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluidInterface;
            if (preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko" || preconditionerMethod == "FaCSI-Block" )
                bcFactoryFluidInterface = Teuchos::rcp( new BCBuilder<SC,LO,GO,NO>( ) );

            if (preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko" || preconditionerMethod == "FaCSI-Block"){
                bcFactoryFluidInterface->addBC(zeroDirichlet3D, 9, 0, domainFluidVelocity, "Dirichlet", dim);
                bcFactoryFluidInterface->addBC(zeroDirichlet3D, 10, 0, domainFluidVelocity, "Dirichlet", dim);
                bcFactoryFluidInterface->addBC(zeroDirichlet3D, 6, 0, domainFluidVelocity, "Dirichlet", dim);
            }

            if ( preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko" || preconditionerMethod == "FaCSI-Block")
               fsi.getPreconditioner()->setFaCSIBCFactory( bcFactoryFluidInterface );

            bcFactoryGeometry->addBC(zeroDirichlet3D, 9, 0, domainGeometry, "Dirichlet", dim); // outflow ring
            bcFactoryGeometry->addBC(zeroDirichlet3D, 10, 0, domainGeometry, "Dirichlet", dim); // inflow ring
            // Die RW, welche nicht Null sind in der rechten Seite (nur Interface) setzen wir spaeter per Hand.
            // Hier erstmal Dirichlet Nullrand, wird spaeter von der Sturkturloesung vorgegeben
            bcFactoryGeometry->addBC(zeroDirichlet3D, 6, 0, domainGeometry, "Dirichlet", dim); // interface

            fsi.problemGeometry_->addBoundaries(bcFactoryGeometry);
           

            // #####################
            // Zeitintegration
            // #####################
            fsi.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

            fsi.initializeProblem();
            
            fsi.initializeGE();
            // Matrizen assemblieren
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
    }

    TimeMonitor_Type::report(std::cout);
    return EXIT_SUCCESS;
}
