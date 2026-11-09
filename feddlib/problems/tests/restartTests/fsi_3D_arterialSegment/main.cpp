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
#include <filesystem>
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

/** @brief Sample the actual start-of-step pressure evaluation for GI restart tests.
 * The uninterrupted run writes samples; a resumed run compares against the same
 * physical times. No extra pressure evaluation is performed, so outlet history
 * advances exactly as it does in the production time loop.
 */
class BoundaryCheckedFSI : public FSI_Type {
public:
    using FSI_Type::FSI_Type;

    void setBoundarySamples(const std::string& output, const std::string& reference)
    {
        boundaryOutput_ = output;
        boundaryReference_ = reference;
    }

    void assemble(std::string type = "") const override
    {
        const bool sample = type == "ComputePressureRHSInTime" && !boundaryOutput_.empty();
        const double previousFlow = sample ? this->getOutletState().previousFlowRate : 0.;
        FSI_Type::assemble(type);
        if (!sample) return;

        double inletArea = 0., outletArea = 0., inletFlow = 0.;
        this->feFactory_->assemblyArea(3, inletArea, 4);
        this->feFactory_->assemblyArea(3, outletArea, 5);
        this->feFactory_->assemblyFlowRate(3, inletFlow, this->getFEType(0), 3, 4, this->u_rep_);
        const FSIOutletState& state = this->getOutletState();
        const double time = this->timeSteppingTool_->currentTime();
        Teuchos::ParameterList values;
        values.set("Inlet area", inletArea).set("Outlet area", outletArea)
            .set("Initial outlet area", state.initialOutletArea)
            .set("Inlet flow", inletFlow).set("Outlet flow", state.currentFlowRate)
            .set("Previous outlet flow", previousFlow)
            .set("Averaged outlet flow", 0.5 * (previousFlow + state.currentFlowRate))
            .set("Outlet pressure", state.pressure);
        checkpoint::onRoot(*this->comm_, [&] {
            const std::string filename = "Boundary_" + std::to_string(time) + ".xml";
            if (!boundaryReference_.empty()) {
                const RCP<Teuchos::ParameterList> reference = Teuchos::getParametersFromXmlFile(
                    joinPath(boundaryReference_, filename));
                const double initialArea = reference->get<double>("Initial outlet area");
                TEUCHOS_TEST_FOR_EXCEPTION(std::abs(reference->get<double>("Outlet area") - initialArea)
                        <= 1.e-8 * initialArea, std::runtime_error,
                    "GI boundary reference must have a deformed outlet to exercise mesh restoration.");
                for (Teuchos::ParameterList::ConstIterator it = values.begin(); it != values.end(); ++it) {
                    const std::string name = values.name(it);
                    const double expected = reference->get<double>(name);
                    const double error = std::abs(values.get<double>(name) - expected);
                    const double allowed = 1.e-12 + 1.e-8 * std::abs(expected);
                    std::cout << "GI boundary restart error at " << time << " (" << name << "): "
                              << error << " (allowed " << allowed << "): "
                              << (std::isfinite(error) && error <= allowed ? "PASS" : "FAIL") << '\n';
                    TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(error) || error > allowed,
                                               std::runtime_error, "GI boundary restart mismatch: " << name);
                }
            }
            std::filesystem::create_directories(boundaryOutput_);
            Teuchos::writeParameterListToXmlFile(values, joinPath(boundaryOutput_, filename));
        });
    }

private:
    std::string boundaryOutput_, boundaryReference_;
};

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
    settings->sublist("Parameter Fluid").set("Viscosity", 0.25).set("Density", 0.5)
        .set("Resistance", 2.).set("BC Ramp", 10.);
    double flow = 0., area = 0.;
    const int backflow = fe.assemblyFlowRate(3, flow, "P1", 3, 5, velocity);
    fe.assemblyArea(3, area, 5);
    vec_dbl_Type flowRates(2, flow);
    auto repeatedLoad = rcp(new MultiVector_Type(domain->getMapVecFieldRepeated()));
    const double pressure = fe.assemblyResistanceBoundary(3, "P1", repeatedLoad, velocity,
        flowRates, settings, 0);
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
    bool passed = std::isfinite(error) && error < 1.e-11 &&
        std::abs(pressure - (backflow ? 0. : 2. * flow)) < 1.e-12;
    if (comm->getRank() == 0)
        std::cout << "Resistance affine-velocity traction error: " << error
                  << ": " << (passed ? "PASS" : "FAIL") << '\n';

    double inletArea = 0., inletFlow = 0.;
    fe.assemblyArea(3, inletArea, 4);
    fe.assemblyFlowRate(3, inletFlow, "P1", 3, 4, velocity);
    auto& fluid = settings->sublist("Parameter Fluid");
    fluid.set("Distal Pressure", 1.25);
    const double constantPressure = fe.assemblyResistanceBoundary(3, "P1", repeatedLoad, velocity,
        flowRates, settings);
    passed = passed && std::abs(constantPressure - (1.25 + 2. * flow)) < 1.e-12;

    fluid.set("Resistance Type", std::string("Relative")).set("Unsteady Start", 0.2);
    // Unequal ratios (2 and 0.5) distinguish the mean of ratios from a ratio
    // of summed areas, and distinguish alpha^-2 from alpha^-1 or alpha^-4.
    const double initialInletArea = inletArea / 2., initialOutletArea = area * 2.;
    const vec_dbl_Type history = {0.25 * flow, 0.75 * flow};
    for (bool average : {false, true}) {
        fluid.set("Average Flowrate", average);
        for (bool useInflow : {false, true}) {
            fluid.set("Use Inflow", useInflow);
            for (double time : {0.1, 0.2, 0.3}) {
                const double effectiveR = time < 0.2 ? 2. : 2. / (1.25 * 1.25);
                const double selectedFlow = useInflow && time < 0.2 ? inletFlow
                    : (average ? 0.5 * flow : 1.25 * flow);
                const double expectedPressure = 1.25 + effectiveR * selectedFlow;
                const double actual = fe.assemblyResistanceBoundary(3, "P1", repeatedLoad, velocity,
                    history, settings, 0, time, initialInletArea, initialOutletArea);
                load.putScalar(0.);
                load.exportFromVector(repeatedLoad, false, "Add");
                std::fill(local, local + 3, 0.);
                const auto assembled = load.getData(0);
                for (int i = 0; i < assembled.size(); ++i) local[i % 3] += assembled[i];
                Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 3, local, total);
                const double expectedLoad[] = {-area * 0.375, -area * 0.75, area * (expectedPressure - 1.125)};
                double positive = 0., negative = 0.;
                for (int i = 0; i < 3; ++i) {
                    positive += std::pow(total[i] - expectedLoad[i], 2);
                    negative += std::pow(total[i] + expectedLoad[i], 2);
                }
                const bool matches = std::isfinite(actual) && std::abs(actual - expectedPressure) < 1.e-12 &&
                    std::sqrt(std::min(positive, negative)) < 1.e-11;
                passed = passed && matches;
                if (comm->getRank() == 0)
                    std::cout << "Relative resistance pressure/traction at " << time
                              << " (Average Flowrate=" << average << ", Use Inflow=" << useInflow << "): "
                              << (matches ? "PASS" : "FAIL") << '\n';
            }
        }
    }
    bool rejectedMissingArea = false;
    try {
        fe.assemblyResistanceBoundary(3, "P1", repeatedLoad, velocity, history, settings);
    } catch (const std::logic_error&) { rejectedMissingArea = true; }
    passed = passed && rejectedMissingArea;
    return passed;
}

// Distinct inlet/outlet flows and retained history exercise source selection
// independently of the nonlinear solve and of the frozen restart references.
bool checkAbsorbingPaperFlowSelection(DomainPtr_Type domain, RCP<const Teuchos::Comm<int>> comm)
{
    FE<SC,LO,GO,NO> fe;
    fe.addFE(domain);
    fe.addFE(domain);
    RCP<MultiVector_Type> velocity = rcp(new MultiVector_Type(domain->getMapVecFieldRepeated()));
    Teuchos::ArrayRCP<SC> values = velocity->getDataNonConst(0);
    const vec2D_dbl_ptr_Type points = domain->getPointsRepeated();
    for (unsigned i = 0; i < points->size(); ++i) {
        values[3 * i] = values[3 * i + 1] = 0.;
        values[3 * i + 2] = 1. + points->at(i).at(2);
    }
    double inlet = 0., outlet = 0., area = 0.;
    fe.assemblyFlowRate(3, inlet, "P1", 3, 4, velocity);
    const int backflow = fe.assemblyFlowRate(3, outlet, "P1", 3, 5, velocity);
    fe.assemblyArea(3, area, 5);
    TEUCHOS_TEST_FOR_EXCEPTION(backflow || std::abs(inlet - outlet) <= 1.e-6 * outlet,
        std::runtime_error, "Flow-source check requires distinct forward inlet/outlet flows.");

    ParameterListPtr_Type settings = rcp(new Teuchos::ParameterList);
    settings->sublist("Parameter Fluid").set("Unsteady Start", 0.2)
        .set("Poisson Ratio", 0.).set("E", 1.).set("Wall thickness", 1.)
        .set("Density", 1.).set("Reference fluid pressure", 2.).set("Flowrate", 1.);
    RCP<MultiVector_Type> load = rcp(new MultiVector_Type(domain->getMapVecFieldRepeated()));
    const vec_dbl_Type history = {0.5 * outlet, outlet};
    const double beta = M_PI / area;
    const double root = std::sqrt(beta * std::sqrt(area));
    const double steadyArea = 1. / ((std::sqrt(root * root + 2.) - root) * 2. * std::sqrt(2.));
    const double expectedInletPressure = std::pow(inlet / (2. * std::sqrt(2.) * steadyArea) + root, 2)
        - root * root;
    bool passed = true;
    for (bool average : {false, true}) {
        settings->sublist("Parameter Fluid").set("Average Flowrate", average);
        for (double time : {0.1, 0.2, 0.3}) {
            settings->sublist("Parameter Fluid").remove("Use Inflow", false);
            const double baseline = fe.assemblyAbsorbingBoundaryPaper(3, "P1", load, velocity,
                history, time, area, area, settings);
            settings->sublist("Parameter Fluid").set("Use Inflow", false);
            const double disabled = fe.assemblyAbsorbingBoundaryPaper(3, "P1", load, velocity,
                history, time, area, area, settings);
            settings->sublist("Parameter Fluid").set("Use Inflow", true);
            const double enabled = fe.assemblyAbsorbingBoundaryPaper(3, "P1", load, velocity,
                history, time, area, area, settings);
            const double expected = time < 0.2 ? expectedInletPressure : baseline;
            const bool matches = disabled == baseline && std::isfinite(enabled) &&
                std::abs(enabled - expected) <= 1.e-12 * std::max(1., std::abs(expected));
            passed = passed && matches;
            if (comm->getRank() == 0)
                std::cout << "Absorbing Paper inlet-flow selection at " << time
                          << " (Average Flowrate=" << average << "): "
                          << (matches ? "PASS" : "FAIL") << '\n';
        }
    }
    return passed;
}

bool compareRestart(FSI_Type& fsi, ParameterListPtr_Type parameters, RCP<const Teuchos::Comm<int>> comm)
{
    const auto& timeParameters = parameters->sublist("Timestepping Parameter");
    const double finalTime = timeParameters.get<double>("Final time");
    const double tolerance = timeParameters.get<double>("Restart tolerance");
    const double absoluteTolerance = timeParameters.get<double>("Restart absolute tolerance");
    const char* checkpointNames[] = {"Solutionu_f", "Solutionp", "Solutiond_s", "Solutiond_f"};
    const char* fieldNames[] = {"fluid velocity", "fluid pressure", "solid displacement", "mesh displacement"};
    const int solutionBlocks[] = {0, 1, 2, 4};
    const int fieldCount = parameters->sublist("Parameter").get("Geometry Explicit", true) ? 3 : 4;
    bool passed = true;

    for (int block = 0; block < fieldCount; ++block) {
        auto solution = fsi.getSolution()->getBlock(solutionBlocks[block]);
        const auto referenceFile = timeParameters.isParameter("Reference directory")
            ? joinPath(timeParameters.get<std::string>("Reference directory"), checkpointNames[block])
            : restartFile(parameters, checkpointNames[block]);
        HDF5Import<SC,LO,GO,NO> importer(solution->getMap(), referenceFile);
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
                timeParameters.isParameter("Reference directory")
                    ? joinPath(timeParameters.get<std::string>("Reference directory"), checkpoint::outletStateName(finalTime))
                    : restartFile(parameters, checkpoint::outletStateName(finalTime)), model, finalTime);
            const auto& state = fsi.getOutletState();
            const bool paper = model == "Absorbing Paper";
            if (state.initialized != reference.initialized ||
                (paper && state.transitionCaptured != reference.transitionCaptured) ||
                (!paper && state.transitionCaptured))
                passed = false;
            const double values[] = {state.initialInletArea, state.initialOutletArea, state.transitionOutletArea,
                                     state.currentFlowRate, state.previousFlowRate, state.pressure};
            const double references[] = {reference.initialInletArea, reference.initialOutletArea, reference.transitionOutletArea,
                                         reference.currentFlowRate, reference.previousFlowRate, reference.pressure};
            const char* names[] = {"initial inlet area", "initial outlet area", "transition outlet area",
                                   "current flow rate", "previous flow rate", "outlet pressure"};
            for (int i = 0; i < 6; ++i) {
                // Legacy plain-absorbing references contain an unused captured area.
                if (i == 2 && !paper) continue;
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
    std::string solverFile = "parametersSolverFSI.xml";
    std::string boundarySamples, boundaryReference;
    std::string checkpointOutput;
    double restartTime = -1.;
    int restartSteps = 0;
    bool validateOnly = false, incompatibleOutlet = false;
    Teuchos::CommandLineProcessor commandLine;
    commandLine.setOption("problemfile", &problemFile, "Case parameters or restart overrides.");
    commandLine.setOption("pressure-model", &pressureModel, "Resistance, Absorbing, or Absorbing Paper.");
    commandLine.setOption("solverfile", &solverFile, "Coupled linear solver settings.");
    commandLine.setOption("boundary-samples", &boundarySamples, "Write start-of-step boundary samples for a GI regression.");
    commandLine.setOption("boundary-reference", &boundaryReference, "Compare boundary samples with an uninterrupted GI run.");
    commandLine.setOption("restart-time", &restartTime, "Override the restart time.");
    commandLine.setOption("restart-steps", &restartSteps, "Override restart length with a fixed step count; zero keeps the configured Final time.");
    commandLine.setOption("checkpoint-output", &checkpointOutput, "Enable checkpoint output to this directory.");
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
    parameters->setParameters(*Teuchos::getParametersFromXmlFile(solverFile));
    parameters->sublist("Parameter Fluid").set("Pressure Boundary Condition", pressureModel);
    if (restartTime >= 0.) parameters->sublist("Timestepping Parameter").set("Time step", restartTime);
    if (!checkpointOutput.empty())
        parameters->sublist("Timestepping Parameter").set("Checkpointing", true)
            .set("Checkpoint directory", checkpointOutput);
    if (incompatibleOutlet) parameters->sublist("Parameter Fluid").set("Average Flowrate", false);
    auto& timeParameters = parameters->sublist("Timestepping Parameter");
    TEUCHOS_TEST_FOR_EXCEPTION(restartSteps < 0, std::logic_error, "restart-steps must be nonnegative.");
    if (timeParameters.get("Restart", false) && restartSteps > 0) {
        const double dt = timeParameters.get<double>("dt");
        TEUCHOS_TEST_FOR_EXCEPTION(restartSteps < 1 || !std::isfinite(dt) || dt <= 0., std::logic_error,
                                   "Restart requires a positive dt and at least one restart step.");
        const double startTime = timeParameters.get<double>("Time step");
        timeParameters.set("Final time", startTime + restartSteps * dt);
        if (comm->getRank() == 0 && !validateOnly)
            std::cout << "Arterial restart: " << restartSteps << " timesteps from " << startTime
                      << " to " << timeParameters.get<double>("Final time") << std::endl;
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

    const bool geometryExplicit = parameters->sublist("Parameter").get("Geometry Explicit", true);
    RCP<SmallMatrix<int>> timeBlocks = rcp(new SmallMatrix<int>(geometryExplicit ? 4 : 5));
    (*timeBlocks)[0][0] = (*timeBlocks)[0][1] = (*timeBlocks)[2][2] = 1;
    if (!geometryExplicit)
        (*timeBlocks)[0][4] = (*timeBlocks)[1][4] = 1;
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
    if (pressureModel == "Absorbing Paper" && !validateOnly && !checkAbsorbingPaperFlowSelection(fluid, comm))
        return EXIT_FAILURE;

    BoundaryCheckedFSI fsi(fluid, discType, fluid, discType, solid, discType,
                 interfaceDomain, discType, fluid, discType,
                 fluidParameters, solidParameters, parameters, geometryParameters, timeBlocks);
    fsi.setBoundarySamples(boundarySamples, boundaryReference);
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
    const bool restarting = timeParameters.get<bool>("Restart");
    const auto* initialClock = restarting
        ? checkpoint::clockState(parameters, timeParameters.get<double>("Time step")) : nullptr;
    const long long initialStep = initialClock ? initialClock->get<long long>("Step number") : 0;
    timeSolver.advanceInTime();
    if (restarting) {
        if (restartSteps > 0) {
            const auto* finalClock = checkpoint::clockState(parameters, timeParameters.get<double>("Final time"));
            TEUCHOS_TEST_FOR_EXCEPTION(!initialClock || !finalClock ||
                                       finalClock->get<long long>("Step number") - initialStep != restartSteps,
                                       std::runtime_error, "The arterial restart did not execute the requested number of timesteps.");
        }
        return compareRestart(fsi, parameters, comm) ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
