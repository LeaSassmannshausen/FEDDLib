#include <Tpetra_Core.hpp>
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/problems/specific/SCI.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include <algorithm>
#include <utility>
#include <vector>

using namespace FEDD;

void zeroDisplacement(double*, double* result, double, const double*)
{
    result[0] = result[1] = result[2] = 0.;
}

void concentrationInflow(double*, double* result, double time, const double* parameters)
{
    result[0] = time >= parameters[0] ? parameters[1] : 0.;
}

void noReaction(double*, double* result, double*) { result[0] = 0.; }

// Same surface-load convention as the supplied SCI artery example:
// parameters = time, pressure, load increment, ramp end, selected flag, face flag.
void arteryPressure(double*, double* result, double* parameters)
{
    double pressure = 0.;
    if (parameters[5] == parameters[4]) {
        const double ramp = parameters[3] > 0.
            ? std::min(1., (parameters[0] + parameters[2]) / parameters[3]) : 1.;
        pressure = parameters[1] * ramp;
    }
    result[0] = result[1] = result[2] = pressure;
}

int main(int argc, char** argv)
{
    Tpetra::ScopeGuard scope(&argc, &argv);
    auto comm = Tpetra::getDefaultComm();
    std::string problemFile = "parametersProblemSCI.xml", solverFile = "parametersSolverSCI.xml";
    std::string structurePrecFile = "parametersPrecStructure.xml", chemPrecFile = "parametersPrecChem.xml";
    std::string explicitPrecFile = "parametersPrecCE.xml";
    Teuchos::CommandLineProcessor options;
    options.setOption("problemfile", &problemFile, "SCI artery problem settings.");
    options.setOption("solverfileSCI", &solverFile, "SCI linear/nonlinear solver settings.");
    options.setOption("precfileStructure", &structurePrecFile, "Coupled SCI preconditioner settings.");
    options.setOption("precfileChem", &chemPrecFile, "Explicit chemistry preconditioner settings.");
    options.setOption("precfileCE", &explicitPrecFile, "Explicit structure preconditioner settings.");
    options.throwExceptions(false);
    const auto parsed = options.parse(argc, argv);
    if (parsed == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
    if (parsed != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) return 1;

    auto problem = Teuchos::getParametersFromXmlFile(problemFile);
    auto solver = Teuchos::getParametersFromXmlFile(solverFile);
    auto& controls = problem->sublist("Parameter");
    const bool explicitChemistry = controls.get("Chemistry Explicit", false);
    TEUCHOS_TEST_FOR_EXCEPTION(controls.get("Dimension", 3) != 3 ||
        controls.get("Discretization", std::string("P2")) != "P2", std::logic_error,
        "sci_artery requires 3D P2 displacement and concentration fields.");

    auto parameters = Teuchos::rcp(new Teuchos::ParameterList(*problem));
    parameters->setParameters(*solver);
    parameters->setParameters(*Teuchos::getParametersFromXmlFile(explicitChemistry ? explicitPrecFile : structurePrecFile));
    auto structure = Teuchos::rcp(new Teuchos::ParameterList(*parameters));
    structure->sublist("Parameter").setParameters(problem->sublist("Parameter Solid"));
    auto chemistry = Teuchos::rcp(new Teuchos::ParameterList(*problem));
    chemistry->setParameters(*solver);
    chemistry->setParameters(*Teuchos::getParametersFromXmlFile(chemPrecFile));
    chemistry->sublist("Parameter").setParameters(problem->sublist("Parameter Diffusion"));

    auto linear = Teuchos::rcp(new Domain<>(comm, 3));
    MeshPartitioner<>::DomainPtrArray_Type domains(1, linear);
    auto partitioning = Teuchos::sublist(parameters, "Mesh Partitioner");
    partitioning->set("Build Edge List", true).set("Build Surface List", true);
    MeshPartitioner<> partitioner(domains, partitioning, "P1", 3);
    partitioner.readAndPartition(15, "mm", true);
    auto solidDomain = Teuchos::rcp(new Domain<>(comm, 3));
    solidDomain->buildP2ofP1Domain(linear);
    // Both fields share the solid mesh and its global element/node ordering.
    auto chemistryDomain = solidDomain;
    if (explicitChemistry) {
        // A separately moved transport mesh must not move the solid reference mesh.
        chemistryDomain = Teuchos::rcp(new Domain<>(comm, 3));
        chemistryDomain->buildP2ofP1Domain(linear);
        chemistryDomain->setReferenceConfiguration();
    }
    solidDomain->setReferenceConfiguration();
    if (parameters->sublist("General").get("ParaView export subdomains", false))
        solidDomain->exportDistribution("SCI");

    vec2D_dbl_Type diffusionTensor(3, vec_dbl_Type(3, 0.));
    for (int d = 0; d < 3; ++d)
        diffusionTensor[d][d] = problem->sublist("Parameter Diffusion").get("D0", 0.001);
    auto timeBlocks = Teuchos::rcp(new SmallMatrix<int>(explicitChemistry ? 1 : 2));
    (*timeBlocks)[0][0] = 1;
    if (!explicitChemistry) (*timeBlocks)[1][1] = 1;
    SCI<> sci(solidDomain, "P2", chemistryDomain, "P2", diffusionTensor, noReaction,
              structure, chemistry, parameters, timeBlocks);

    auto boundaries = Teuchos::rcp(new BCBuilder<>());
    auto solidBoundaries = Teuchos::rcp(new BCBuilder<>());
    auto chemBoundaries = Teuchos::rcp(new BCBuilder<>());
    // Identical solid constraints to fsi_artery_boundary_tests, on block zero.
    for (const auto& boundary : std::vector<std::pair<int, std::string>>{
             {14, "Dirichlet_Y_Z"}, {13, "Dirichlet_X_Z"}, {7, "Dirichlet_Z"}, {8, "Dirichlet_Z"},
             {9, "Dirichlet_Z"}, {10, "Dirichlet_Z"}, {11, "Dirichlet_Z"}, {12, "Dirichlet_Z"}}) {
        boundaries->addBC(zeroDisplacement, boundary.first, 0, solidDomain, boundary.second, 3);
        solidBoundaries->addBC(zeroDisplacement, boundary.first, 0, solidDomain, boundary.second, 3);
    }
    std::vector<double> inflow = {controls.get("Inflow Start Time", 0.), controls.get("Inflow Concentration", 1.)};
    const int chemFlag = controls.get("Chemistry Flag", 6);
    if (!explicitChemistry)
        boundaries->addBC(concentrationInflow, chemFlag, 1, chemistryDomain, "Dirichlet", 1, inflow);
    chemBoundaries->addBC(concentrationInflow, chemFlag, 0, chemistryDomain, "Dirichlet", 1, inflow);
    sci.problemChem_->addBoundaries(chemBoundaries);
    sci.problemStructureNonLin_->addBoundaries(solidBoundaries);
    sci.problemStructureNonLin_->addRhsFunction(arteryPressure, 0);
    sci.problemStructureNonLin_->addParemeterRhs(controls.get("Volume force", 106.6));
    sci.problemStructureNonLin_->addParemeterRhs(controls.get("Load Step Size", 0.005));
    sci.problemStructureNonLin_->addParemeterRhs(controls.get("Load Ramp End", 0.1));
    sci.problemStructureNonLin_->addParemeterRhs(controls.get("Pressure Flag", 6));
    sci.addBoundaries(boundaries);
    sci.initializeProblem();
    sci.initializeCE();
    sci.assemble();

    DAESolverInTime<> timeSolver(parameters, comm);
    timeSolver.defineTimeStepping(*timeBlocks);
    timeSolver.setProblem(sci);
    timeSolver.setupTimeStepping();
    timeSolver.advanceInTime();
    return 0;
}
