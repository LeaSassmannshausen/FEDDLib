#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NonLinLaplace.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

void zeroDirichlet(double *x, double *res, double t, const double *parameters) {
  res[0] = 0.;

  return;
}

void zeroDirichlet2D(double *x, double *res, double t,
                     const double *parameters) {
  res[0] = 0.;
  res[1] = 0.;

  return;
}

void zeroDirichlet3D(double *x, double *res, double t,
                     const double *parameters) {
  res[0] = 0.;
  res[1] = 0.;
  res[2] = 0.;

  return;
}

void zeroDirichletX(double *x, double *res, double t,
                    const double *parameters) {
  res[0] = 0.;
  res[1] = x[1];
  res[2] = x[2];

  return;
}

void zeroDirichletY(double *x, double *res, double t,
                    const double *parameters) {
  res[0] = x[0];
  res[1] = 0.;
  res[2] = x[2];

  return;
}

void zeroDirichletZ(double *x, double *res, double t,
                    const double *parameters) {
  res[0] = x[0];
  res[1] = x[1];
  res[2] = 0.;

  return;
}

void dummyFunc(double *x, double *res, double t, const double *parameters) {
  return;
}

void rhs2D(double *x, double *res, double *parameters) {
  // parameters[0] is the time, not needed here
  res[0] = 0.;
  res[1] = parameters[1];

  return;
}

void rhsY(double *x, double *res, double *parameters) {
  // parameters[0] is the time, not needed here
  res[0] = 0.;
  res[1] = parameters[1];
  res[2] = 0.;
  return;
}

void rhsX(double *x, double *res, double *parameters) {
  // parameters[0] is the time, not needed here
  res[0] = parameters[1];
  res[1] = 0.;
  res[2] = 0.;
  return;
}

void rhsYZ(double *x, double *res, double *parameters) {
  // parameters[0] is the time, not needed here
  res[0] = 0.;
  double force = parameters[1];

  if (parameters[2] == 5)
    res[1] = force;
  else
    res[1] = 0.;

  if (parameters[2] == 4)
    res[2] = force;
  else
    res[2] = 0.;

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
int main(int argc, char *argv[]) {

  typedef MeshUnstructured<SC, LO, GO, NO> MeshUnstr_Type;
  typedef RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
  typedef Domain<SC, LO, GO, NO> Domain_Type;
  typedef RCP<Domain_Type> DomainPtr_Type;
  typedef ExporterParaView<SC, LO, GO, NO> ExporterPV_Type;
  typedef RCP<ExporterPV_Type> ExporterPVPtr_Type;
  typedef MeshPartitioner<SC, LO, GO, NO> MeshPartitioner_Type;

  typedef Map<LO, GO, NO> Map_Type;
  typedef RCP<Map_Type> MapPtr_Type;
  typedef MultiVector<SC, LO, GO, NO> MultiVector_Type;
  typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
  typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
  typedef BlockMultiVector<SC, LO, GO, NO> BlockMultiVector_Type;
  typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  Teuchos::RCP<const Teuchos::Comm<int>> comm =
      Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

  // Command Line Parameters
  Teuchos::CommandLineProcessor myCLP;
  string ulib_str = "Tpetra";
  myCLP.setOption("ulib", &ulib_str, "Underlying lib");
  // int dim = 2;
  // myCLP.setOption("dim",&dim,"dim");
  string xmlProblemFile = "parametersProblem.xml";
  myCLP.setOption("problemfile", &xmlProblemFile,
                  ".xml file with Inputparameters.");
  string xmlPrecFile = "parametersPrec.xml";
  myCLP.setOption("precfile", &xmlPrecFile, ".xml file with Inputparameters.");
  string xmlSolverFile = "parametersSolver.xml";
  myCLP.setOption("solverfile", &xmlSolverFile,
                  ".xml file with Inputparameters.");

  myCLP.recogniseAllOptions(true);
  myCLP.throwExceptions(false);
  Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn =
      myCLP.parse(argc, argv);
  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
    mpiSession.~GlobalMPISession();
    return 0;
  }

  bool verbose(comm->getRank() == 0); // Only first rank prints
  if (verbose) {
    cout << "###############################################################"
         << endl;
    cout << "############ Starting nonlinear Laplace ... ############" << endl;
    cout << "###############################################################"
         << endl;
  }

  {
    ParameterListPtr_Type parameterListProblem =
        Teuchos::getParametersFromXmlFile(xmlProblemFile);
    ParameterListPtr_Type parameterListPrec =
        Teuchos::getParametersFromXmlFile(xmlPrecFile);
    ParameterListPtr_Type parameterListSolver =
        Teuchos::getParametersFromXmlFile(xmlSolverFile);

    ParameterListPtr_Type parameterListAll(
        new Teuchos::ParameterList(*parameterListProblem));
    parameterListAll->setParameters(*parameterListPrec);
    parameterListAll->setParameters(*parameterListSolver);

    int dim = parameterListProblem->sublist("Parameter").get("Dimension", 2);
    string meshType = parameterListProblem->sublist("Parameter")
                          .get("Mesh Type", "structured");
    string meshName = parameterListProblem->sublist("Parameter")
                          .get("Mesh Name", "square.mesh");
    string meshDelimiter =
        parameterListProblem->sublist("Parameter").get("Mesh Delimiter", " ");
    int m = parameterListProblem->sublist("Parameter").get("H/h", 5);
    string FEType =
        parameterListProblem->sublist("Parameter").get("Discretization", "P1");

    int numProcsCoarseSolve =
        parameterListProblem->sublist("General").get("Mpi Ranks Coarse", 0);
    int size = comm->getSize() - numProcsCoarseSolve;

    Teuchos::RCP<Teuchos::Time> totalTime(
        Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
    Teuchos::RCP<Teuchos::Time> buildMesh(
        Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
    Teuchos::RCP<Teuchos::Time> solveTime(
        Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));

    DomainPtr_Type domain;

    // ########################
    // P1 und P2 Gitter bauen
    // ########################

    domain.reset(new Domain<SC, LO, GO, NO>(comm, dim));
    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
    domainP1Array[0] = domain;

    ParameterListPtr_Type pListPartitioner =
        sublist(parameterListProblem, "Mesh Partitioner");
    MeshPartitioner<SC, LO, GO, NO> partitionerP1(domainP1Array,
                                                  pListPartitioner, "P1", dim);

    partitionerP1.readAndPartition();
    if (FEType == "P2") {
      Teuchos::RCP<Domain<SC, LO, GO, NO>> domainP2;
      domainP2.reset(new Domain_Type(comm, dim));
      domainP2->buildP2ofP1Domain(domain);
      domain = domainP2;
    }

    // ########################
    // Flags setzen
    // ########################

    Teuchos::RCP<BCBuilder<SC, LO, GO, NO>> bcFactory(
        new BCBuilder<SC, LO, GO, NO>());
    if (dim == 2)
      bcFactory->addBC(zeroDirichlet2D, 1, 0, domain, "Dirichlet", dim);
    else if (dim == 3) {

      bcFactory->addBC(zeroDirichlet, 1, 0, domain, "Dirichlet_X", dim);
      bcFactory->addBC(zeroDirichlet, 2, 0, domain, "Dirichlet_Y", dim);
      bcFactory->addBC(zeroDirichlet, 3, 0, domain, "Dirichlet_Z", dim);
      bcFactory->addBC(zeroDirichlet3D, 0, 0, domain, "Dirichlet", dim);
      bcFactory->addBC(zeroDirichlet2D, 7, 0, domain, "Dirichlet_X_Y", dim);
      bcFactory->addBC(zeroDirichlet2D, 8, 0, domain, "Dirichlet_Y_Z", dim);
      bcFactory->addBC(zeroDirichlet2D, 9, 0, domain, "Dirichlet_X_Z", dim);
    }

    // Declare nonlinlaplace object
    NonLinLaplace<SC, LO, GO, NO> NonLinLaplace(domain, FEType,
                                                parameterListAll);

    NonLinLaplace.addBoundaries(bcFactory);

    if (dim == 2)
      NonLinLaplace.addRhsFunction(rhs2D);
    else if (dim == 3)
      NonLinLaplace.addRhsFunction(rhsYZ);

    // ######################
    // Assemble matrix, set boundary conditions and solve
    // ######################
    NonLinLaplace.initializeProblem();
    NonLinLaplace.assemble();
    NonLinLaplace.setBoundaries();
    NonLinLaplace.setBoundariesRHS();

    std::string nlSolverType =
        parameterListProblem->sublist("General").get("Linearization", "NOX");
    NonLinearSolver<SC, LO, GO, NO> nlSolverAssFE(nlSolverType);
    nlSolverAssFE.solve(NonLinLaplace);
    comm->barrier();
  }
  return (EXIT_SUCCESS);
}
