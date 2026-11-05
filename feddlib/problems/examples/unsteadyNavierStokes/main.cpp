#include <Tpetra_Core.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_StackedTimer.hpp>
#include <algorithm>
#include <cmath>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"

/*!
 main of time-dependent Navier-Stokes problem
 
 @brief time-dependent Navier-Stokes main
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;

void zeroDirichlet(double* x, double* res, double t, const double* parameters){

    res[0] = 0.;

    return;
}

void zeroDirichlet2D(double* x, double* res, double t, const double* parameters){

    res[0] = 0.;
    res[1] = 0.;

    return;
}

void couette2D(double* x, double* res, double t, const double* parameters){

    res[0] = 1.*parameters[0];
    res[1] = 0.;

    return;
}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters){

    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void inflowPartialCFD(double* x, double* res, double t, const double* parameters){

    double H = parameters[1];

    if(t < 0.5)
    {
        res[0] = (4.0*1.5*parameters[0]*x[1]*(H-x[1])/(H*H))*((1 - cos(2.0*M_PI*t))/2.0);
        res[1] = 0.;
    }
    else
    {
        res[0] = 4.0*1.5*parameters[0]*x[1]*(H-x[1])/(H*H);
        res[1] = 0.;
    }

    return;
}

void inflowParabolic2D(double* x, double* res, double t, const double* parameters){

    double H = parameters[1];
    res[0] = 4.*parameters[0]*x[1]*(H-x[1])/(H*H);
    res[1] = 0.;
    
    return;
}

void inflowParabolic2DSin(double* x, double* res, double t, const double* parameters){

    double H = parameters[1];
    res[0] = sin(M_PI*t*0.125)*( 6*x[1]*(H-x[1]) ) / (H*H);
    res[1] = 0.;

    return;
}

void inflowParabolic3D(double* x, double* res, double t, const double* parameters){

    double H = parameters[1];
    res[0] = 16*parameters[0]*x[1]*(H-x[1])*x[2]*(H-x[2])/(H*H*H*H);
    res[1] = 0.;
    res[2] = 0.;

    return;
}

// Elliptical artery inlet. Direction and principal axis are orthonormalized
// when registering the boundary condition below.
void inflowParabolicArtery(double* x, double* res, double t, const double* parameters)
{
    const double* direction = parameters + 5;
    const double* axis = parameters + 8;
    const double minorAxis[3] = {
        direction[1]*axis[2] - direction[2]*axis[1],
        direction[2]*axis[0] - direction[0]*axis[2],
        direction[0]*axis[1] - direction[1]*axis[0]};
    double major = 0., minor = 0.;
    for (int i = 0; i < 3; ++i) {
        const double offset = x[i] - parameters[2+i];
        major += offset * axis[i];
        minor += offset * minorAxis[i];
    }
    const double radiusSquared = major*major/(parameters[11]*parameters[11])
                               + minor*minor/(parameters[12]*parameters[12]);
    double ramp = 1.;
    if (parameters[1] > 0.) {
        const double phase = std::max(0., std::min(1., t/parameters[1]));
        ramp = 0.5 * (1. - std::cos(M_PI * phase));
    }
    const double speed = parameters[0] * ramp * std::max(0., 1. - radiusSquared);
    for (int i = 0; i < 3; ++i)
        res[i] = speed * direction[i];
}

void inflow3DRichter(double* x, double* res, double t, const double* parameters)
{
    double H = parameters[1];
    
    if(t < 1.)
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) ) * ((1 - cos(2.0*M_PI*t))/2.0);
        res[1] = 0.;
        res[2] = 0.;
    }
    else
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) );
        res[1] = 0.;
        res[2] = 0.;
    }
    
    return;
}


void dummyFunc(double* x, double* res, double t, const double* parameters){

    return;
}

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;
using namespace Teuchos;

using namespace FEDD;
int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    typedef Teuchos::RCP<Domain<SC,LO,GO,NO> > DomainPtr_Type;

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    bool verbose (comm->getRank() == 0);
    if (verbose) {
        std::cout << "###############################################################" << std::endl;
        std::cout << "################### Unsteady Navier-Stokes ####################" << std::endl;
        std::cout << "###############################################################" << std::endl;
    }

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");

    string xmlTekoPrecFile = "parametersTeko.xml";
    myCLP.setOption("tekoprecfile",&xmlTekoPrecFile,".xml file with Inputparameters.");
    string xmlBlockPrecFile = "parametersPrecBlock.xml";
    myCLP.setOption("blockprecfile",&xmlBlockPrecFile,".xml file with Inputparameters.");
   
    double length = 4.;
    myCLP.setOption("length",&length,"length of domain.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        return EXIT_SUCCESS;
    }

    Teuchos::RCP<StackedTimer> stackedTimer =  rcp(new StackedTimer("Unsteady Navier-Stokes",true));
    TimeMonitor::setStackedTimer(stackedTimer);

    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);

        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);

        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListPrecTeko = Teuchos::getParametersFromXmlFile(xmlTekoPrecFile);

        ParameterListPtr_Type parameterListPrecBlock = Teuchos::getParametersFromXmlFile(xmlBlockPrecFile);

        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",3);
        std::string feTypeV = parameterListProblem->sublist("Parameter").get("Discretization Velocity","P2");
        std::string feTypeP = parameterListProblem->sublist("Parameter").get("Discretization Pressure","P1");
        string		meshType = parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        string		meshName = parameterListProblem->sublist("Parameter").get("Mesh Name","structured");
        string		meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        int 		m = parameterListProblem->sublist("Parameter").get("H/h",5);
        string		linearization = parameterListProblem->sublist("General").get("Linearization","FixedPoint");
        string		precMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        bool computeInflow = parameterListProblem->sublist("Parameter").get("Compute Inflow",false);
        int         n;


        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        if (!precMethod.compare("Monolithic"))
            parameterListAll->setParameters(*parameterListPrec);
        else if(precMethod == "Teko")
            parameterListAll->setParameters(*parameterListPrecTeko);
        else if(precMethod == "Diagonal" || precMethod == "Triangular" || precMethod == "PCD" || precMethod == "LSC")
            parameterListAll->setParameters(*parameterListPrecBlock);

        parameterListAll->setParameters(*parameterListSolver);

        std::string bcType = parameterListProblem->sublist("Parameter").get("BC Type","parabolic");

        int minNumberSubdomains;
        if (!meshType.compare("structured")) {
            minNumberSubdomains = 1;
        }
        else if(!meshType.compare("structured_rec")){
            minNumberSubdomains = length;
        }
        else if(!meshType.compare("structured_bfs")){
            minNumberSubdomains = (int) 2*length+1;
        }

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        int size = comm->getSize() - numProcsCoarseSolve;

        double viscosity = parameterListProblem->sublist("Parameter").get("Viscosity",1.e-3);

        Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
        Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
        Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));
        
        {
            DomainPtr_Type domainPressure;
            DomainPtr_Type domainVelocity;

            Teuchos::TimeMonitor totalTimeMonitor(*totalTime);
            {
                Teuchos::TimeMonitor buildMeshMonitor(*buildMesh);
                if (verbose)
                    cout << "-- Building Mesh ..." << flush;
                
                if (!meshType.compare("structured")) {
                    if (dim == 2) {
                        n = (int) (std::pow(size,1/2.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(2);
                        x[0]=0.0;    x[1]=0.0;
                        domainPressure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
                        domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
                    }
                    else if (dim == 3){
                        n = (int) (std::pow(size,1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(3);
                        x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
                        domainPressure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
                        domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
                    }
                    domainPressure->buildMesh( 1,"Square", dim, feTypeP, n, m, numProcsCoarseSolve);
                    domainVelocity->buildMesh( 1,"Square", dim, feTypeV, n, m, numProcsCoarseSolve);
                }
                if (!meshType.compare("structured_bfs")) {
                            TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured BFS mesh.");
                            if (dim == 2) {
                                n = (int) (std::pow( size/minNumberSubdomains ,1/2.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                                std::vector<double> x(2);
                                x[0]=-1.0;    x[1]=-1.0;
                                domainPressure.reset(new Domain<SC,LO,GO,NO>( x, length+1., 2., comm ) );
                                domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, length+1., 2., comm ) );
                            }
                            else if (dim == 3){
                                n = (int) (std::pow( size/minNumberSubdomains ,1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                                std::vector<double> x(3);
                                x[0]=-1.0;    x[1]=0.0;    x[2]=-1.0;
                                domainPressure.reset(new Domain<SC,LO,GO,NO>( x, length+1., 1., 2., comm));
                                domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, length+1., 1., 2., comm));
                            }
                            domainPressure->buildMesh( 2,"BFS", dim, feTypeP, n, m, numProcsCoarseSolve);
                            domainVelocity->buildMesh( 2,"BFS", dim, feTypeV, n, m, numProcsCoarseSolve);
                }
                else if (!meshType.compare("unstructured")) {
                    domainPressure.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
                    domainVelocity.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
                    
                    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
                    domainP1Array[0] = domainPressure;
                    
                    ParameterListPtr_Type pListPartitioner = sublist( parameterListProblem, "Mesh Partitioner" );
                    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
                    
                    partitionerP1.readAndPartition();
                    
                    domainVelocity->buildP2ofP1Domain( domainPressure );
                    
                    if (feTypeV=="P2")
                        domainVelocity->buildP2ofP1Domain( domainPressure );
                    else
                        domainVelocity = domainPressure;
                    
                }
            }
            domainVelocity->exportNodeFlags();

            std::vector<double> parameter_vec(1);
            if ( !bcType.compare("parabolic") || !bcType.compare("parabolic_benchmark") || !bcType.compare("parabolic_benchmark_sin") )
                parameter_vec[0] = parameterListProblem->sublist("Parameter").get("MaxVelocity",1.5);
            else if ( !bcType.compare("partialCFD") ) //  Fuer CFD3
                parameter_vec[0] = parameterListProblem->sublist("Parameter").get("MeanVelocity",2.);
                        
            // ####################
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

            if (!bcType.compare("parabolic") || !bcType.compare("Couette"))
                parameter_vec.push_back(1.);//height of inflow region
            else if(!bcType.compare("parabolic_benchmark_sin") || !bcType.compare("parabolic_benchmark") || !bcType.compare("partialCFD"))
                parameter_vec.push_back(.41);//height of inflow region
            else if(!bcType.compare("Richter3D"))
                parameter_vec.push_back(.4);
            else if (bcType != "parabolic_artery")
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Select a valid boundary condition.");

            if (bcType == "parabolic_artery") {
                TEUCHOS_TEST_FOR_EXCEPTION(dim != 3, std::logic_error,
                    "parabolic_artery requires a 3D mesh.");
                auto& inflow = parameterListProblem->sublist("Parameter");
                parameter_vec = {
                    inflow.get("MaxVelocity", 1.5), inflow.get("Inflow Ramp Time", 0.),
                    inflow.get("Inflow Center X", 0.), inflow.get("Inflow Center Y", 0.), inflow.get("Inflow Center Z", 0.),
                    inflow.get("Inflow Direction X", 0.), inflow.get("Inflow Direction Y", 0.), inflow.get("Inflow Direction Z", 1.),
                    inflow.get("Inflow Axis X", 1.), inflow.get("Inflow Axis Y", 0.), inflow.get("Inflow Axis Z", 0.),
                    inflow.get("Inflow Radius 1", 1.), inflow.get("Inflow Radius 2", 1.)};
                for (const double value : parameter_vec)
                    TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(value), std::logic_error,
                        "parabolic_artery parameters must be finite.");
                TEUCHOS_TEST_FOR_EXCEPTION(parameter_vec[0] < 0. || parameter_vec[1] < 0. ||
                    parameter_vec[11] <= 0. || parameter_vec[12] <= 0., std::logic_error,
                    "parabolic_artery requires positive radii and nonnegative velocity/ramp time.");
                double directionNorm = 0.;
                for (int i = 0; i < 3; ++i)
                    directionNorm += parameter_vec[5+i]*parameter_vec[5+i];
                directionNorm = std::sqrt(directionNorm);
                TEUCHOS_TEST_FOR_EXCEPTION(directionNorm == 0., std::logic_error,
                    "Inflow Direction must be nonzero and point into the fluid.");
                for (int i = 0; i < 3; ++i)
                    parameter_vec[5+i] /= directionNorm;
                double projection = 0., axisNorm = 0.;
                for (int i = 0; i < 3; ++i)
                    projection += parameter_vec[8+i]*parameter_vec[5+i];
                for (int i = 0; i < 3; ++i) {
                    parameter_vec[8+i] -= projection*parameter_vec[5+i];
                    axisNorm += parameter_vec[8+i]*parameter_vec[8+i];
                }
                axisNorm = std::sqrt(axisNorm);
                TEUCHOS_TEST_FOR_EXCEPTION(axisNorm < 1.e-12, std::logic_error,
                    "Inflow Axis must have a component perpendicular to Inflow Direction.");
                for (int i = 0; i < 3; ++i)
                    parameter_vec[8+i] /= axisNorm;
                const int wallFlag = inflow.get("Wall Flag", 1);
                const int inletFlag = inflow.get("Inflow Flag", 2);
                TEUCHOS_TEST_FOR_EXCEPTION(wallFlag == inletFlag, std::logic_error,
                    "Wall Flag and Inflow Flag must differ.");
                bcFactory->addBC(zeroDirichlet3D, wallFlag, 0, domainVelocity, "Dirichlet", dim);
                bcFactory->addBC(inflowParabolicArtery, inletFlag, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
                if (verbose)
                    std::cout << "parabolic_artery: inlet flag " << inletFlag << ", wall flag " << wallFlag
                              << ", inward direction (" << parameter_vec[5] << ", " << parameter_vec[6]
                              << ", " << parameter_vec[7] << ")" << std::endl;
            }
            else if (!bcType.compare("parabolic") || !bcType.compare("parabolic_benchmark")) {//flag of obstacle
                if (dim==2){
                    bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(inflowParabolic2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
//                    bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                    bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                    bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
                }
                else if (dim==3){
                    bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(inflowParabolic3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
//                    bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                    bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                    bcFactory->addBC(zeroDirichlet3D, 4, 0, domainVelocity, "Dirichlet", dim);
                    
                }
            }
            else if (!bcType.compare("parabolic_benchmark_sin")) {
                if (dim==2){
                    bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(inflowParabolic2DSin, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
                    bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
                }
                else if (dim==3){
                    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Test for 2D only: parabolic_benchmark_sin");
                }
            }
            else if (!bcType.compare("partialCFD")) { // Fuer CFD3 Test
                if (dim==2)
                {
                    bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(inflowPartialCFD, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
//                    bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                    bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                    bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(zeroDirichlet2D, 5, 0, domainVelocity, "Dirichlet", dim);
                }
                else if (dim==3){
                    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No partial CFD test in 3D.");
                }
            }
            else if (!bcType.compare("Richter3D")) {
                bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim); // wall
                bcFactory->addBC(inflow3DRichter, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // inflow
                bcFactory->addBC(zeroDirichlet3D, 3, 0, domainVelocity, "Dirichlet_Z", dim);
                bcFactory->addBC(zeroDirichlet3D, 5, 0, domainVelocity, "Dirichlet", dim);


            }
            else if (!bcType.compare("Couette")){
              bcFactory->addBC(couette2D, 1, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // wall
              bcFactory->addBC(zeroDirichlet2D, 2, 0, domainVelocity, "Dirichlet", dim); // wall
              //bcFactory->addBC(inflow3DRichter, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // inflow
              //bcFactory->addBC(zeroDirichlet3D, 3, 0, domainVelocity, "Dirichlet_Z", dim);
              //bcFactory->addBC(zeroDirichlet3D, 5, 0, domainVelocity, "Dirichlet", dim);
            
            }
            
            int timeDisc = parameterListProblem->sublist("Timestepping Parameter").get("Butcher table",0);

            NavierStokes<SC,LO,GO,NO> navierStokes( domainVelocity, feTypeV, domainPressure, feTypeP, parameterListAll );

            navierStokes.addBoundaries(bcFactory);

            navierStokes.initializeProblem();
            
            navierStokes.assemble();

            navierStokes.setBoundariesRHS();

            DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);
            SmallMatrix<int> defTS(2);
            defTS[0][0] = 1;
            defTS[0][1] = 1;
            defTS[1][0] = 0;
            defTS[1][1] = 0;

            daeTimeSolver.defineTimeStepping(defTS);

            daeTimeSolver.setProblem(navierStokes);

            daeTimeSolver.setupTimeStepping();

            daeTimeSolver.advanceInTime();

            navierStokes.infoParameter();

        }
    }

    Teuchos::TimeMonitor::report(cout);
    stackedTimer->stop("Unsteady Navier-Stokes");
	StackedTimer::OutputOptions options;
	options.output_fraction = options.output_histogram = options.output_minmax = true;
	stackedTimer->report((std::cout),comm,options);

    return EXIT_SUCCESS;
}
