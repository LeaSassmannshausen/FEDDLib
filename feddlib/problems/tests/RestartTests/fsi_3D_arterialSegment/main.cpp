#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/Checkpointing/CheckpointFiles.hpp"
#include "feddlib/core/Checkpointing/CheckpointMetadata.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/Preconditioner.hpp"
#include "feddlib/problems/specific/FSI.hpp"
#include "feddlib/problems/specific/Geometry.hpp"
#include "feddlib/problems/specific/LinElas.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace FEDD;
using Teuchos::RCP;
using Teuchos::rcp;

typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;
typedef Domain<SC,LO,GO,NO> Domain_Type;
typedef RCP<Domain_Type> DomainPtr_Type;
typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
typedef BCBuilder<SC,LO,GO,NO> BCBuilder_Type;
typedef FSI<SC,LO,GO,NO> FSI_Type;
typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;

namespace {

void zeroDirichlet(double*, double* result, double, const double*)
{
    result[0] = result[1] = result[2] = 0.;
}

// BCBuilder supplies the scalar profile as x[0] and normalizes its integral
// to the desired flow rate on the current inlet mesh.
void parabolicInflow(double* x, double* result, double, const double* parameters)
{
    result[0] = result[1] = 0.;
    result[2] = parameters[0] * x[0];
}

void flowRate(double*, double* result, double time, const double* parameters)
{
    const double rampTime = parameters[1];
    const double ramp = time < rampTime ? 0.5 * (1. - std::cos(M_PI * time / rampTime)) : 1.;
    result[0] = parameters[2] * ramp;
}

void constantResistance(double*, double* result, double*) { result[0] = 2.; }

// An affine velocity has constant grad(u)*n. Summed nodal loads on the planar
// outlet must equal its area times the analytic pressure and viscous traction.
bool checkResistanceAssembly(DomainPtr_Type domain, RCP<const Teuchos::Comm<int>> comm)
{
    FE<SC,LO,GO,NO> fe;
    fe.addFE(domain);
    fe.addFE(domain);
    auto velocity = rcp(new MultiVector_Type(domain->getMapVecFieldRepeated()));
    auto values = velocity->getDataNonConst(0);
    const auto points = domain->getPointsRepeated();
    for (unsigned i = 0; i < points->size(); ++i) {
        const auto& x = points->at(i);
        values[3 * i] = 1. + x[0] + 2. * x[1] + 3. * x[2];
        values[3 * i + 1] = 2. + 4. * x[0] + 5. * x[1] + 6. * x[2];
        values[3 * i + 2] = 3. + 7. * x[0] + 8. * x[1] + 9. * x[2];
    }
    auto settings = rcp(new Teuchos::ParameterList);
    settings->sublist("Parameter Fluid").set("Viscosity", 0.25).set("Density", 0.5);
    double flow = 0., area = 0.;
    const int backflow = fe.assemblyFlowRate(3, flow, "P1", 3, 5, velocity);
    fe.assemblyArea(3, area, 5);
    vec_dbl_Type flowRates(2, flow), time(1, 0.);
    auto repeatedLoad = rcp(new MultiVector_Type(domain->getMapVecFieldRepeated()));
    const double pressure = fe.assemblyResistanceBoundary(3, "P1", repeatedLoad, velocity,
        flowRates, time, constantResistance, settings, 0);
    MultiVector_Type load(domain->getMapVecFieldUnique());
    load.exportFromVector(repeatedLoad, false, "Add");
    double local[3] = {0., 0., 0.}, total[3] = {0., 0., 0.};
    const auto data = load.getData(0);
    for (int i = 0; i < data.size(); ++i) local[i % 3] += data[i];
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 3, local, total);
    const double expected[] = {-area * 0.375, -area * 0.75, area * (pressure - 1.125)};
    double errorPositive = 0., errorNegative = 0.;
    for (int i = 0; i < 3; ++i) {
        errorPositive += std::pow(total[i] - expected[i], 2);
        errorNegative += std::pow(total[i] + expected[i], 2);
    }
    const double error = std::sqrt(std::min(errorPositive, errorNegative));
    const bool passed = std::isfinite(error) && error < 1.e-11 &&
        std::abs(pressure - (backflow ? 0. : 2. * flow)) < 1.e-12;
    if (comm->getRank() == 0)
        std::cout << "Resistance affine-velocity traction error: " << error
                  << ": " << (passed ? "PASS" : "FAIL") << '\n';
    return passed;
}

bool compareRestart(FSI_Type& fsi, ParameterListPtr_Type parameters, RCP<const Teuchos::Comm<int>> comm)
{
    const auto& timeParameters = parameters->sublist("Timestepping Parameter");
    const double finalTime = timeParameters.get<double>("Final time");
    const double tolerance = timeParameters.get<double>("Restart tolerance");
    const double absoluteTolerance = timeParameters.get<double>("Restart absolute tolerance");
    const char* checkpointNames[] = {"Solutionu_f", "Solutionp", "Solutiond_s"};
    const char* fieldNames[] = {"fluid velocity", "fluid pressure", "solid displacement"};
    bool passed = true;

    for (int block = 0; block < 3; ++block) {
        auto solution = fsi.getSolution()->getBlock(block);
        HDF5Import<SC,LO,GO,NO> importer(solution->getMap(), restartFile(parameters, checkpointNames[block]));
        auto reference = importer.readVariablesHDF5(std::to_string(finalTime));
        MultiVector_Type error(solution->getMap());
        error.update(1., *solution, -1., *reference, 0.);

        Teuchos::Array<SC> errorNorm(1), referenceNorm(1);
        error.norm2(errorNorm);
        reference->norm2(referenceNorm);
        // A small reference field needs an absolute round-off allowance as
        // well as a relative bound. Always scale by the reference solution.
        const double allowed = absoluteTolerance + tolerance * referenceNorm[0];
        const bool matches = std::isfinite(errorNorm[0]) && std::isfinite(referenceNorm[0]) &&
                             errorNorm[0] <= allowed;
        if (comm->getRank() == 0) {
            std::cout << "Restart absolute error (" << fieldNames[block] << ", l2): "
                      << errorNorm[0] << " (allowed " << allowed << " = " << absoluteTolerance
                      << " + " << tolerance << " * reference l2 norm): "
                      << (matches ? "PASS" : "FAIL") << std::endl;
            std::cout << "Restart relative error (" << fieldNames[block] << "): ";
            if (referenceNorm[0] > 0.) std::cout << errorNorm[0] / referenceNorm[0];
            else std::cout << "undefined (zero reference norm; using absolute bound)";
            std::cout << " (relative tolerance " << tolerance << ")" << std::endl;
        }
        passed = passed && matches;
    }
    const std::string model = parameters->sublist("Parameter Fluid").get("Pressure Boundary Condition", "None");
    if (model != "None") {
        checkpoint::onRoot(*comm, [&] {
            const auto reference = checkpoint::readOutletState(
                restartFile(parameters, checkpoint::outletStateName(finalTime)), model, finalTime);
            const auto& state = fsi.getOutletState();
            if (state.initialized != reference.initialized || state.transitionCaptured != reference.transitionCaptured)
                passed = false;
            const double values[] = {state.initialInletArea, state.initialOutletArea, state.transitionOutletArea,
                                     state.currentFlowRate, state.previousFlowRate, state.pressure};
            const double references[] = {reference.initialInletArea, reference.initialOutletArea, reference.transitionOutletArea,
                                         reference.currentFlowRate, reference.previousFlowRate, reference.pressure};
            const char* names[] = {"initial inlet area", "initial outlet area", "transition outlet area",
                                   "current flow rate", "previous flow rate", "outlet pressure"};
            for (int i = 0; i < 6; ++i) {
                const double error = std::abs(values[i] - references[i]);
                const double allowed = 1.e-14 + tolerance * std::abs(references[i]);
                const bool matches = std::isfinite(error) && error <= allowed;
                std::cout << "Restart absolute error (" << names[i] << "): " << error
                          << " (allowed " << allowed << "): " << (matches ? "PASS" : "FAIL") << '\n';
                passed = passed && matches;
            }
        });
    }
    int localPassed = passed ? 1 : 0, allPassed = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &localPassed, &allPassed);
    return allPassed != 0;
}


} // namespace

// Three-dimensional arterial segment in centimetres. The supplied geometry
// has radius 0.09 cm and length 0.05 cm (0.5 mm); no coordinate scaling is used.
int main(int argc, char* argv[])
{
    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);
    auto comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    const std::string baseProblemFile = "parametersProblemFSI.xml";
    std::string problemFile = baseProblemFile;
    std::string pressureModel = "Absorbing Paper";
    double restartTime = -1.;
    bool saveAll = false, validateOnly = false, incompatibleOutlet = false;
    Teuchos::CommandLineProcessor commandLine;
    commandLine.setOption("problemfile", &problemFile, "Case parameters or restart overrides.");
    commandLine.setOption("pressure-model", &pressureModel, "Resistance, Absorbing, or Absorbing Paper.");
    commandLine.setOption("restart-time", &restartTime, "Override the restart time.");
    commandLine.setOption("save-all", "selected-checkpoints", &saveAll, "Exercise saving every solution.");
    commandLine.setOption("validate-checkpoint", "solve", &validateOnly, "Validate checkpoint before simulation.");
    commandLine.setOption("incompatible-outlet", "compatible-outlet", &incompatibleOutlet, "Change averaging for a rejection test.");
    commandLine.throwExceptions(false);
    const auto parsed = commandLine.parse(argc, argv);
    if (parsed == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return EXIT_SUCCESS;
    if (parsed != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) return EXIT_FAILURE;

    auto parameters = Teuchos::getParametersFromXmlFile(baseProblemFile);
    if (problemFile != baseProblemFile)
        parameters->setParameters(*Teuchos::getParametersFromXmlFile(problemFile));
    parameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersPrecGE.xml"));
    parameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersSolverFSI.xml"));
    parameters->sublist("Parameter Fluid").set("Pressure Boundary Condition", pressureModel);
    if (restartTime >= 0.) parameters->sublist("Timestepping Parameter").set("Time step", restartTime);
    if (incompatibleOutlet) parameters->sublist("Parameter Fluid").set("Average Flowrate", false);
    if (saveAll) {
        parameters->sublist("General").set("Safe all solution", true);
        parameters->sublist("Timestepping Parameter").set("Checkpointing", false);
    }

    auto fluidParameters = Teuchos::getParametersFromXmlFile("parametersPrecFluidMono.xml");
    fluidParameters->sublist("Parameter").setParameters(parameters->sublist("Parameter Fluid"));
    fluidParameters->sublist("Timestepping Parameter").setParameters(parameters->sublist("Timestepping Parameter"));
    fluidParameters->sublist("General").set("Preconditioner Method", "Monolithic");
    fluidParameters->sublist("General").set("Flag Inlet Fluid", 4).set("Flag Outlet Fluid", 5).set("Flag Interface", 6);

    auto solidParameters = Teuchos::getParametersFromXmlFile("parametersPrecStructure.xml");
    solidParameters->sublist("Parameter").setParameters(parameters->sublist("Parameter Solid"));
    solidParameters->sublist("Parameter").set("Use AceGen Interface", true);
    solidParameters->sublist("Parameter Solid").setParameters(parameters->sublist("Parameter Solid"));
    solidParameters->sublist("Timestepping Parameter").setParameters(parameters->sublist("Timestepping Parameter"));

    auto geometryParameters = Teuchos::getParametersFromXmlFile("parametersPrecGeometry.xml");
    geometryParameters->setParameters(*Teuchos::getParametersFromXmlFile("parametersSolverGeometry.xml"));
    geometryParameters->sublist("General").set("Preconditioner Method", "MonolithicConstPrec");
    geometryParameters->sublist("Parameter").setParameters(parameters->sublist("Parameter Geometry"));

    const int dim = 3;
    const std::string discType = "P1";
    DomainPtr_Type fluid = rcp(new Domain_Type(comm, dim));
    DomainPtr_Type solid = rcp(new Domain_Type(comm, dim));
    MeshPartitioner_Type::DomainPtrArray_Type domains(2);
    domains[0] = fluid;
    domains[1] = solid;
    auto partitionerParameters = Teuchos::sublist(parameters, "Mesh Partitioner");
    partitionerParameters->set("Build Edge List", true).set("Build Surface List", true);
    MeshPartitioner_Type partitioner(domains, partitionerParameters, discType, dim);
    partitioner.readAndPartition(15);
    vec_int_Type interfaceIds = {6, 9, 10};
    fluid->identifyInterfaceParallelAndDistance(solid, interfaceIds);
    fluid->buildInterfaceMaps();
    solid->buildInterfaceMaps();
    DomainPtr_Type interfaceDomain = rcp(new Domain_Type(comm));
    interfaceDomain->setDummyInterfaceDomain(fluid);
    fluid->setReferenceConfiguration();
    solid->setReferenceConfiguration();

    RCP<SmallMatrix<int>> timeBlocks = rcp(new SmallMatrix<int>(4));
    (*timeBlocks)[0][0] = (*timeBlocks)[0][1] = (*timeBlocks)[2][2] = 1;
    // Build the reference inlet profile before FSI restores and moves the mesh.
    // Its nodal values must match those used by the uninterrupted run.
    auto profile = rcp(new MultiVector_Type(fluid->getMapUnique()));
    auto profileValues = profile->getDataNonConst(0);
    const auto points = fluid->getPointsUnique();
    const double radius = 0.09;
    for (unsigned i = 0; i < points->size(); ++i) {
        const auto& x = points->at(i);
        profileValues[i] = std::max(0., 1. - (x[0] * x[0] + x[1] * x[1]) / (radius * radius));
    }

    if (pressureModel == "Resistance" && !validateOnly && !checkResistanceAssembly(fluid, comm))
        return EXIT_FAILURE;

    FSI_Type fsi(fluid, discType, fluid, discType, solid, discType,
                 interfaceDomain, discType, fluid, discType,
                 fluidParameters, solidParameters, parameters, geometryParameters, timeBlocks);
    if (validateOnly) return EXIT_SUCCESS;

    auto boundaries = rcp(new BCBuilder_Type());
    auto fluidBoundaries = rcp(new BCBuilder_Type());
    auto solidBoundaries = rcp(new BCBuilder_Type());
    auto geometryBoundaries = rcp(new BCBuilder_Type());
    auto interfaceBoundaries = rcp(new BCBuilder_Type());
    const auto& fluidSettings = parameters->sublist("Parameter Fluid");
    std::vector<double> inflowParameters = {1., fluidSettings.get<double>("Max Ramp Time"),
                                           fluidSettings.get<double>("Flowrate")};
    RCP<const MultiVector_Type> inflowProfile = profile;
    for (auto bc : {boundaries, fluidBoundaries}) {
        bc->addBC(parabolicInflow, 4, 0, fluid, "Dirichlet", dim, inflowParameters, inflowProfile, true, flowRate);
        bc->addBC(zeroDirichlet, 9, 0, fluid, "Dirichlet_Z", dim);
    }

    // The inlet/outlet are fixed longitudinally. Two strips remove transverse
    // rigid-body motion while allowing radial expansion of the arterial wall.
    for (const auto& boundary : std::vector<std::pair<int, std::string>>{
             {14, "Dirichlet_Y_Z"}, {13, "Dirichlet_X_Z"},
             {7, "Dirichlet_Z"}, {8, "Dirichlet_Z"}, {9, "Dirichlet_Z"}, {10, "Dirichlet_Z"}}) {
        boundaries->addBC(zeroDirichlet, boundary.first, 2, solid, boundary.second, dim);
        solidBoundaries->addBC(zeroDirichlet, boundary.first, 0, solid, boundary.second, dim);
    }
    for (int flag : {6, 9, 10}) {
        geometryBoundaries->addBC(zeroDirichlet, flag, 0, fluid, "Dirichlet", dim);
        interfaceBoundaries->addBC(zeroDirichlet, flag, 0, fluid, "Dirichlet", dim);
    }
    fsi.problemFluid_->addBoundaries(fluidBoundaries);
    fsi.problemStructureNonLin_->addBoundaries(solidBoundaries);
    fsi.problemGeometry_->addBoundaries(geometryBoundaries);
    fsi.getPreconditioner()->setFaCSIBCFactory(interfaceBoundaries);
    fsi.addBoundaries(boundaries);
    fsi.initializeProblem();
    fsi.initializeGE();
    fsi.assemble();

    DAESolverInTime<SC,LO,GO,NO> timeSolver(parameters, comm);
    timeSolver.defineTimeStepping(*timeBlocks);
    timeSolver.setProblem(fsi);
    timeSolver.setupTimeStepping();
    timeSolver.advanceInTime();
    if (parameters->sublist("Timestepping Parameter").get<bool>("Restart"))
        return compareRestart(fsi, parameters, comm) ? EXIT_SUCCESS : EXIT_FAILURE;
    return EXIT_SUCCESS;
}
