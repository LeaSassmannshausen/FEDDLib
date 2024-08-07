#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/specific/Stokes.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 main of Stokes unit test

 @brief Stokes main

 This laplace unit test compares the current solution of the laplace problem in
 2D and 3D to older already stored solutions. The folling test constructs a
 structured mesh with H/h = 6. The stored solutions make sense for the following
 configurations:
 - 2D on 4 procs with P1 or P2 elements
 - 3D on 8 procs with P1 or P2 elements

 */

void zeroBC(double *x, double *res, double t, const double *parameters) { res[0] = 0.; }
void oneFunc(double *x, double *res, double *parameters) { res[0] = 1.; }

void zeroDirichlet2D(double *x, double *res, double t, const double *parameters) {

    res[0] = 0.;
    res[1] = 0.;

    return;
}

void zeroDirichlet3D(double *x, double *res, double t, const double *parameters) {

    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void inflowParabolic2D(double *x, double *res, double t, const double *parameters) {

    double H = parameters[1];
    res[0] = 4 * parameters[0] * x[1] * (H - x[1]) / (H * H);
    res[1] = 0.;

    return;
}

void inflowParabolic3D(double *x, double *res, double t, const double *parameters) {

    double H = parameters[1];
    res[0] = 16 * parameters[0] * x[1] * (H - x[1]) * x[2] * (H - x[2]) / (H * H * H * H);
    res[1] = 0.;
    res[2] = 0.;

    return;
}

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;

int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC, LO, GO, NO> MeshPartitioner_Type;

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

    Teuchos::RCP<const Teuchos::Comm<int>> comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib", &ulib_str, "Underlying lib");

    int dim = 2;
    myCLP.setOption("dim", &dim, "Dimension");

    std::string FETypeV = "P2";
    myCLP.setOption("FEType", &FETypeV, "Discretization");

    std::string FEType = "P1";

    string xmlPrecFile;
    if(dim==2)
        xmlPrecFile = "parametersPrec_Stokes_2D.xml";
    else if(dim==3)
        xmlPrecFile = "parametersPrec_Stokes_3D.xml";
    
    myCLP.setOption("precfile", &xmlPrecFile, ".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile", &xmlSolverFile, ".xml file with Inputparameters.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc, argv);
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    {
        // ParameterListPtr_Type parameterListProblem =
        // Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListPrec));
        parameterListAll->setParameters(*parameterListSolver);

        // Mesh
        int m = 4;
        std::string meshType = "structured";
        std::string meshDelimiter = " ";

        int n;
        int numProcsCoarseSolve =0;
        int size = comm->getSize();

        int minNumberSubdomains = 1;

        Teuchos::RCP<Domain<SC, LO, GO, NO>> domainPressure;
        Teuchos::RCP<Domain<SC, LO, GO, NO>> domainVelocity;

        if (dim == 2) {
            TEUCHOS_TEST_FOR_EXCEPTION(std::floor(std::sqrt(size)) != std::sqrt(size), std::logic_error, "Wrong number of processors for structured squared mesh.");
            n = (int)(std::pow(size / minNumberSubdomains, 1 / 2.) + 100 * Teuchos::ScalarTraits<double>::eps()); // 1/H
            std::vector<double> x(2);
            x[0] = 0.0;
            x[1] = 0.0;
            domainPressure.reset(new Domain<SC, LO, GO, NO>(x, 1., 1., comm));
            domainVelocity.reset(new Domain<SC, LO, GO, NO>(x, 1., 1., comm));
        } else if (dim == 3) {
            TEUCHOS_TEST_FOR_EXCEPTION(std::floor(std::pow(size, 1 / 3)) != std::pow(size, 1 / 3), std::logic_error, "Wrong number of processors for structured squared mesh.");
            n = (int)(std::pow(size / minNumberSubdomains, 1 / 3.) + 100 * Teuchos::ScalarTraits<double>::eps()); // 1/H
            std::vector<double> x(3);
            x[0] = 0.0;
            x[1] = 0.0;
            x[2] = 0.0;
            domainPressure.reset(new Domain<SC, LO, GO, NO>(x, 1., 1., 1., comm));
            domainVelocity.reset(new Domain<SC, LO, GO, NO>(x, 1., 1., 1., comm));
        }
        domainPressure->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
        domainVelocity->buildMesh(1, "Square", dim, FETypeV, n, m, numProcsCoarseSolve);

        // ####################
        std::vector<double> parameter_vec(1, 1.0); // Velocity
        parameter_vec.push_back(1.);               // height of inflow region

        Teuchos::RCP<BCBuilder<SC, LO, GO, NO>> bcFactory(new BCBuilder<SC, LO, GO, NO>());
        if (dim == 2) { // When it is not LDC ..
            bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
            bcFactory->addBC(inflowParabolic2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
        } else if (dim == 3) {
            bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim);
            bcFactory->addBC(inflowParabolic3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
        }

        Stokes<SC, LO, GO, NO> stokes(domainVelocity, FETypeV, domainPressure, FEType, parameterListAll);

        {
            stokes.addBoundaries(bcFactory);
            stokes.initializeProblem();
            stokes.assemble();
            stokes.setBoundaries();
            stokes.solve();
        }

        bool boolExportSolution = true;
        if (boolExportSolution) {

            // HDF5Export<SC, LO, GO, NO> exporter(stokes.getSolution()->getBlock(0)->getMap(),
            //                                     "ReferenceSolutions/solution_stokes_" + std::to_string(dim) + "d_" + FETypeV + "_" + std::to_string(size) + "cores"); //  Map and file name
            // exporter.writeVariablesHDF5("solution",
            //                             stokes.getSolution()->getBlock(0)); // VariableName and Variable

            // We exclude any other tests, than the one prescribed
            if (dim == 2) {
                TEUCHOS_TEST_FOR_EXCEPTION(!(size == 4 && m == 4), std::logic_error, "The 2D test solutions are only sensible for 4 processors.");
            } else if (dim == 3)
                TEUCHOS_TEST_FOR_EXCEPTION(!(size == 8 && m == 4), std::logic_error, "The 3D test solutions are only sensible for 8 processors.");

            HDF5Import<SC, LO, GO, NO> importer(stokes.getSolution()->getBlock(0)->getMap(),
                                                "ReferenceSolutions/solution_stokes_" + std::to_string(dim) + "d_" + FETypeV + "_" + std::to_string(size) + "cores");
            Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionImported = importer.readVariablesHDF5("solution");

            // We compare the imported solution to the current one
            Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionStokes = stokes.getSolution()->getBlock(0);
            // Calculating the error per node
            Teuchos::RCP<MultiVector<SC, LO, GO, NO>> errorValues = Teuchos::rcp(new MultiVector<SC, LO, GO, NO>(stokes.getSolution()->getBlock(0)->getMap()));
            // this = alpha*A + beta*B + gamma*this
            errorValues->update(1., solutionStokes, -1., solutionImported, 0.);
            // Computing norm
            Teuchos::Array<SC> norm(1);
            errorValues->norm2(norm);
            double normError = norm[0];

            // Output of error
            if (comm->getRank() == 0) {
                cout << " --------------------------------------------------" << endl;
                cout << "  Error Report " << endl;
                cout << "   || solution_current - solution_stored||_2 = " << normError << endl;
                cout << " --------------------------------------------------" << endl;
            }
            // Throwing exception, if error is too great.

            TEUCHOS_TEST_FOR_EXCEPTION(normError > 1.e-12, std::logic_error,
                                       "Difference between current solution and "
                                       "stored solution greater than 1e-12.");

            Teuchos::RCP<ExporterParaView<SC, LO, GO, NO>> exPara(new ExporterParaView<SC, LO, GO, NO>());

            Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> exportSolution = stokes.getSolution()->getBlock(0);

            exPara->setup("solutionStokes", domainVelocity->getMesh(), FEType);

            exPara->addVariable(exportSolution, "u_current", "Vector", dim, domainVelocity->getMapUnique());
            exPara->addVariable(solutionImported, "u_imported", "Vector", dim, domainVelocity->getMapUnique());

            exPara->save(0.0);
        }
    }
    return (EXIT_SUCCESS);
}
