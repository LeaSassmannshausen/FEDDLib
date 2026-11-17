#include <Tpetra_Core.hpp>
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/specific/FSI.hpp"
#include "feddlib/problems/specific/FSCI.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/Preconditioner.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace FEDD;
using Teuchos::rcp;
using SC = default_sc;
using LO = default_lo;
using GO = default_go;
using NO = default_no;
using DomainType = Domain<SC,LO,GO,NO>;
using DomainPtr = Teuchos::RCP<DomainType>;
using Vector = MultiVector<SC,LO,GO,NO>;
using VectorPtr = Teuchos::RCP<Vector>;
using FSIType = FSI<SC,LO,GO,NO>;
using FSCIType = FSCI<SC,LO,GO,NO>;
using Boundaries = BCBuilder<SC,LO,GO,NO>;
using CommConstPtr_Type = Teuchos::RCP<const Teuchos::Comm<int>>;

/*!
 @brief Compare mechanics-only FSI with zero-concentration FSCI.

 As in the AssembleFE comparison tests, main() constructs and solves the two
 problem types explicitly, then compares their results. Each solve owns fresh
 domains, clocks and component parameter lists: the first ALE solve must not
 change the second solve's initial configuration.
*/

namespace {
constexpr int numberOfTimeSteps = 2;

// Command-line filenames follow the AssembleFE tests. Both formulations read
// the same physical/solver inputs; only their formulation and block flags differ.
struct InputFiles {
    std::string problem = "parameters.xml";
    std::string solver = "parametersSolver.xml";
    std::string fluidPreconditioner = "parametersPrecFluidMono.xml";
    std::string solidPreconditioner = "parametersPrecStructure.xml";
    std::string geometryPreconditioner = "parametersPrecGeometry.xml";
    std::string geometrySolver = "parametersSolverGeometry.xml";
};

struct ComponentParameters {
    ParameterListPtr_Type fluid, solid, geometry;
};

struct Domains {
    DomainPtr fluidPressure, fluidVelocity, solid, interfaceDomain;
};

struct ComparisonSettings {
    double absoluteTolerance, relativeTolerance, concentrationTolerance;
};

// Boundary callbacks. The vector-valued zero callback is used only for the
// mechanical fields; chemistry has no imposed concentration or source.
void zero(double*, double* result, double, const double*)
{
    result[0] = result[1] = result[2] = 0.;
}
void noReaction(double*, double* result, double*) { result[0] = 0.; }
void concentrationBoundary(double*, double* result, double, const double*) { result[0] = .2; }
// BCBuilder passes the scalar nodal inlet-profile value in profile[0], and
// replaces parameters[0] with the scale needed to obtain the prescribed flux.
// Inlet flag 4 is at z=0; positive z points into the arterial segment.
void inlet(double* profile, double* result, double, const double* parameters)
{
    result[0] = result[1] = 0.;
    result[2] = parameters[0]*profile[0];
}
// parameters = {profile scale, ramp duration, target volumetric flow rate}.
// This callback returns Q(t), not a maximum velocity. A zero ramp is a step.
void flowRate(double*, double* result, double time, const double* parameters)
{
    const double phase = parameters[1] > 0. ? std::min(1., std::max(0., time/parameters[1])) : 1.;
    result[0] = parameters[2]*.5*(1.-std::cos(std::acos(-1.)*phase));
}

// Deep copies retain both accepted steps after the live solution is overwritten.
struct Frame {
    double time, pressure, concentrationNorm;
    std::vector<VectorPtr> fields;
    VectorPtr concentration;
};
struct Result {
    std::vector<Frame> frames;
    std::vector<double> linear, nonlinear;
    double dt;
};

// Sample accepted solutions before UpdateTime changes material/Newmark history.
// These small subclasses only observe accepted steps; the inherited FSI/FSCI
// assembly and solution methods are unchanged. In particular, do not call the
// pressure assembler here: that would advance the explicit outlet state again.
template<class Base>
class ObservedProblem : public Base {
public:
    using Base::Base;
    std::vector<Frame> frames;
    void computeValuesOfInterestAndExport() override
    {
        Frame frame;
        frame.time = this->timeSteppingTool_->currentTime() + this->timeSteppingTool_->get_dt();
        frame.pressure = this->getPressureOutlet();
        frame.concentrationNorm = 0.;
        // Store the common fields in the fixed comparison order
        // {u_f, p, d_s, lambda, d_f}. For GE the live solution's block 4 is
        // chemistry in FSCI, so always take d_f from the geometry subproblem.
        for (int block = 0; block < 5; ++block) {
            auto values = block == 4 ? this->problemGeometry_->getSolution()->getBlock(0)
                                     : this->getSolution()->getBlock(block);
            auto copy = rcp(new Vector(values->getMap()));
            copy->update(1., *values, 0.);
            frame.fields.push_back(copy);
        }
        if (auto fsci = dynamic_cast<FSCIType*>(this)) {
            Teuchos::Array<double> norm(1);
            fsci->getSCIProblem()->getSolution()->getBlock(1)->normInf(norm());
            frame.concentrationNorm = norm[0];
            auto c = fsci->getSCIProblem()->getSolution()->getBlock(1);
            frame.concentration = rcp(new Vector(c->getMap()));
            frame.concentration->update(1., *c, 0.);
        }
        frames.push_back(frame);
        Base::computeValuesOfInterestAndExport();
    }
};

// Keep these names explicit at the construction sites in main(). The wrapper
// forwards every constructor argument to the corresponding production class.
using FSIValidationProblem = ObservedProblem<FSIType>;
using FSCIValidationProblem = ObservedProblem<FSCIType>;

// The DAE solver writes average linear iterations and Newton counts to separate
// value-only files. Their .times companions associate each row with an accepted
// physical time, rather than relying on a line number alone.
std::vector<double> readLog(const std::string& file, double dt)
{
    std::ifstream values(file), times(file + ".times");
    TEUCHOS_TEST_FOR_EXCEPTION(!values || !times, std::runtime_error, "Missing iteration log " << file);
    std::vector<double> result;
    double value, time;
    while (values >> value) {
        TEUCHOS_TEST_FOR_EXCEPTION(!(times >> time) || !std::isfinite(time) || !std::isfinite(value) || value < 0. ||
            std::abs(time - (result.size()+1)*dt) > 1.e-12, std::runtime_error,
            "Invalid iteration log or timestamp in " << file);
        result.push_back(value);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!values.eof() || (times >> time) || !times.eof() || result.size() != numberOfTimeSteps,
        std::runtime_error, "Expected exactly two accepted steps in " << file);
    return result;
}

// Each formulation receives independent lists. Copying the same named solid
// settings into both locations is required by the FSI component and the native
// AceGen material adapter, respectively.
ComponentParameters buildComponentParameters(const ParameterListPtr_Type& parameters, const InputFiles& files)
{
    auto fluidParameters = Teuchos::getParametersFromXmlFile(files.fluidPreconditioner);
    fluidParameters->sublist("Parameter") = parameters->sublist("Parameter Fluid");
    fluidParameters->sublist("General").set("Preconditioner Method", "Monolithic")
        .set("Flag Inlet Fluid", 4).set("Flag Outlet Fluid", 5).set("Flag Interface", 6);
    auto solidParameters = Teuchos::getParametersFromXmlFile(files.solidPreconditioner);
    solidParameters->sublist("Parameter") = parameters->sublist("Parameter Solid");
    solidParameters->sublist("Parameter Solid") = parameters->sublist("Parameter Solid");
    auto geometryParameters = Teuchos::getParametersFromXmlFile(files.geometryPreconditioner);
    geometryParameters->setParameters(*Teuchos::getParametersFromXmlFile(files.geometrySolver));
    geometryParameters->sublist("Parameter") = parameters->sublist("Parameter Geometry");
    geometryParameters->sublist("General").set("Preconditioner Method", "MonolithicConstPrec");

    return {fluidParameters, solidParameters, geometryParameters};
}

// The restart mesh uses centimetres already: radius 0.09, length 0.05.
// Do not scale these coordinates. Pressure uses the original P1 mesh. Velocity
// uses P2 here, as do the solid and chemistry required by the native material.
Domains buildDomains(const ParameterListPtr_Type& parameters, CommConstPtr_Type comm)
{
    const int dim = 3;
    auto fluidP1 = rcp(new DomainType(comm, dim));
    auto solidP1 = rcp(new DomainType(comm, dim));
    using Partitioner = MeshPartitioner<SC,LO,GO,NO>;
    Partitioner::DomainPtrArray_Type domains(2);
    domains[0] = fluidP1;
    domains[1] = solidP1;
    auto partitionerParameters = Teuchos::sublist(parameters, "Mesh Partitioner");
    partitionerParameters->set("Build Edge List", true).set("Build Surface List", true);
    Partitioner partitioner(domains, partitionerParameters, "P1", dim);
    partitioner.readAndPartition(15);
    auto fluid = rcp(new DomainType(comm, dim));
    auto solid = rcp(new DomainType(comm, dim));
    fluid->buildP2ofP1Domain(fluidP1);
    solid->buildP2ofP1Domain(solidP1);
    vec_int_Type interfaceIds{6, 9, 10};
    fluidP1->identifyInterfaceParallelAndDistance(solidP1, interfaceIds);
    fluid->identifyInterfaceParallelAndDistance(solid, interfaceIds);
    fluid->buildInterfaceMaps();
    solid->buildInterfaceMaps();
    auto interfaceDomain = rcp(new DomainType(comm));
    interfaceDomain->setDummyInterfaceDomain(fluid);
    // Keep the material reference configuration and the ALE initial mesh
    // identical in both runs. These domains are recreated for the next run.
    fluidP1->setReferenceConfiguration();
    fluid->setReferenceConfiguration();
    solid->setReferenceConfiguration();

    return {fluidP1, fluid, solid, interfaceDomain};
}

// Coupled solution ordering:
// FSI  GE: [u_f, p, d_s, lambda]       GI: [..., d_f]
// FSCI GE: [u_f, p, d_s, lambda, c]    GI: [..., d_f, c]
// The entries below identify blocks with temporal terms. The SCI material
// already supplies Mc/dt and the chemical rate residual. The FSCI DAE path
// recognizes the chemistry block and avoids adding a second chemical mass term.
Teuchos::RCP<SmallMatrix<int>> buildTimeSteppingDefinition(bool geometryExplicit, bool includeChemistry)
{
    const int chemistryBlock = geometryExplicit ? 4 : 5;
    auto timeBlocks = rcp(new SmallMatrix<int>(chemistryBlock + (includeChemistry ? 1 : 0)));
    (*timeBlocks)[0][0] = (*timeBlocks)[0][1] = (*timeBlocks)[2][2] = 1;
    if (!geometryExplicit) (*timeBlocks)[0][4] = (*timeBlocks)[1][4] = 1;
    if (includeChemistry) (*timeBlocks)[chemistryBlock][chemistryBlock] = 1;
    return timeBlocks;
}

// Common mechanical boundary setup. Return the solid BCs so main() can also
// attach a copy to SCI's local displacement block for the FSCI run.
Teuchos::RCP<Boundaries> addMechanicalBoundaries(FSIType& problem, const Domains& domains,
                                                const ParameterListPtr_Type& parameters, bool loadChemistry = false)
{
    const int dim = 3;
    auto fluid = domains.fluidVelocity;
    auto solid = domains.solid;
    auto boundaries = rcp(new Boundaries);
    auto fluidBoundaries = rcp(new Boundaries);
    auto solidBoundaries = rcp(new Boundaries);
    auto geometryBoundaries = rcp(new Boundaries);
    auto interfaceBoundaries = rcp(new Boundaries);
    // Fluid inlet: interpolate 1-r^2/R^2 on the P2 nodes, then let BCBuilder
    // normalize its actual mesh flux to Q(t). Flag 9 is the inlet/interface rim.
    // The moving fluid wall obtains its velocity through the FSI coupling.
    auto profile = rcp(new Vector(fluid->getMapUnique()));
    auto values = profile->getDataNonConst(0);
    const auto points = fluid->getPointsUnique();
    for (unsigned i = 0; i < points->size(); ++i) {
        const auto& x = points->at(i);
        values[i] = std::max(0., 1. - (x[0]*x[0]+x[1]*x[1])/(.09*.09));
    }
    const auto& fluidSettings = parameters->sublist("Parameter Fluid");
    vec_dbl_Type inflow{1., fluidSettings.get<double>("Max Ramp Time"), fluidSettings.get<double>("Flowrate")};
    Teuchos::RCP<const Vector> inletProfile = profile;
    for (auto bc : {boundaries, fluidBoundaries}) {
        bc->addBC(inlet, 4, 0, fluid, "Dirichlet", dim, inflow, inletProfile, true, flowRate);
        bc->addBC(zero, 9, 0, fluid, "Dirichlet_Z", dim);
    }
    // Outlet flag 5: the FSI time-step assembly supplies the resistance
    // traction from "Parameter Fluid". Its applied pressure is recorded by
    // ObservedProblem along with the accepted solution at each step.
    // Solid: axial constraints on the ends/rims, plus transverse constraints
    // on flags 13 and 14 to remove rigid motion while allowing radial expansion.
    // The coupled displacement is block 2; the solid subproblem uses block 0.
    for (const auto& bc : std::vector<std::pair<int, std::string>>{
             {14, "Dirichlet_Y_Z"}, {13, "Dirichlet_X_Z"}, {7, "Dirichlet_Z"},
             {8, "Dirichlet_Z"}, {9, "Dirichlet_Z"}, {10, "Dirichlet_Z"}}) {
        boundaries->addBC(zero, bc.first, 2, solid, bc.second, dim);
        solidBoundaries->addBC(zero, bc.first, 0, solid, bc.second, dim);
    }
    // Geometry/interface preconditioner: flags 6, 9, 10 form the interface.
    // Geometry's nonzero displacement is supplied by coupling. These BCs must
    // stay on its subproblem, or global elimination would remove that coupling.
    for (int flag : {6, 9, 10}) {
        geometryBoundaries->addBC(zero, flag, 0, fluid, "Dirichlet", dim);
        interfaceBoundaries->addBC(zero, flag, 0, fluid, "Dirichlet", dim);
    }
    problem.problemFluid_->addBoundaries(fluidBoundaries);
    problem.problemStructureNonLin_->addBoundaries(solidBoundaries);
    problem.problemGeometry_->addBoundaries(geometryBoundaries);
    problem.getPreconditioner()->setFaCSIBCFactory(interfaceBoundaries);
    if (loadChemistry)
        boundaries->addBC(concentrationBoundary, 7, parameters->sublist("Parameter").get("Geometry Explicit", true) ? 4 : 5,
            solid, "Dirichlet", 1);
    problem.addBoundaries(boundaries);
    return solidBoundaries;
}

// Initialization and time integration are shared, just as both AssembleFE
// variants use the same DAE workflow. This helper does not construct a problem.
template<class Problem>
Result advanceTwoSteps(Problem& problem, ParameterListPtr_Type parameters,
                       Teuchos::RCP<SmallMatrix<int>>& timeBlocks, CommConstPtr_Type comm)
{
    // initializeGE prepares the separate ALE solve when geometry is explicit.
    // In GI the geometry block instead belongs to the coupled Newton system.
    problem.initializeProblem();
    problem.initializeGE();
    problem.assemble();
    DAESolverInTime<SC,LO,GO,NO> solver(parameters, comm);
    solver.defineTimeStepping(*timeBlocks);
    solver.setProblem(problem);
    solver.setupTimeStepping();
    solver.advanceInTime();

    Result result;
    result.frames = problem.frames;
    result.dt = parameters->sublist("Timestepping Parameter").get<double>("dt");
    const auto suffix = parameters->sublist("General").get<std::string>("Export Suffix");
    result.linear = readLog("linearIterations" + suffix + ".txt", result.dt);
    result.nonlinear = readLog("newtonIterations" + suffix + ".txt", result.dt);
    return result;
}

bool compare(const Result& fsi, const Result& fsci, const ComparisonSettings& tolerances, CommConstPtr_Type comm)
{
    TEUCHOS_TEST_FOR_EXCEPTION(fsi.frames.size() != numberOfTimeSteps || fsci.frames.size() != numberOfTimeSteps,
        std::runtime_error, "Both problems must complete exactly two timesteps.");
    const char* names[] = {"u_f", "p", "d_s", "lambda", "d_f"};
    bool passed = true;
    auto check = [&](const std::string& name, double error, double allowed) {
        const bool matches = std::isfinite(error) && std::isfinite(allowed) && error <= allowed;
        if (comm->getRank() == 0)
            std::cout << name << ": error=" << error << ", allowed=" << allowed
                      << (matches ? " PASS" : " FAIL") << std::endl;
        passed = passed && matches;
    };
    for (int step = 0; step < numberOfTimeSteps; ++step) {
        const auto& a = fsi.frames[step];
        const auto& b = fsci.frames[step];
        if (comm->getRank() == 0) std::cout << "Comparing accepted step " << step+1 << ":" << std::endl;
        // Check against the requested clock as well as against the other run;
        // two equally stale component clocks must not pass this comparison.
        check("FSI accepted time", std::abs(a.time-(step+1)*fsi.dt), 1.e-12);
        check("FSCI accepted time", std::abs(b.time-(step+1)*fsci.dt), 1.e-12);
        check("matching time", std::abs(a.time-b.time), 1.e-12);
        for (int block = 0; block < 5; ++block) {
            // Import by global ID; do not assume equal local mesh partitioning.
            auto difference = rcp(new Vector(b.fields[block]->getMap()));
            difference->importFromVector(a.fields[block]);
            difference->update(-1., *b.fields[block], 1.);
            // L2 measures the whole distributed field; Linf also catches a
            // localized discrepancy. Print each reference/computed norm to make
            // the scale of the absolute and relative errors easy to inspect.
            for (bool infinityNorm : {false, true}) {
                Teuchos::Array<double> error(1), reference(1), computed(1);
                if (infinityNorm) {
                    difference->normInf(error());
                    a.fields[block]->normInf(reference());
                    b.fields[block]->normInf(computed());
                } else {
                    difference->norm2(error());
                    a.fields[block]->norm2(reference());
                    b.fields[block]->norm2(computed());
                }
                const std::string label = std::string(names[block]) + (infinityNorm ? " Linf" : " L2");
                if (comm->getRank() == 0) {
                    std::cout << label << ": FSI=" << reference[0] << ", FSCI=" << computed[0]
                              << ", relative error=";
                    if (reference[0] > 0.) std::cout << error[0]/reference[0];
                    else std::cout << "n/a (zero reference; absolute bound applies)";
                    std::cout << std::endl;
                }
                check(label, error[0], tolerances.absoluteTolerance + tolerances.relativeTolerance*reference[0]);
            }
        }
        // This is the pressure actually applied during the accepted step.
        // The resistance BC is explicit, so it uses the previous flow state.
        if (comm->getRank() == 0)
            std::cout << "pressureOutlet: FSI=" << a.pressure << ", FSCI=" << b.pressure << std::endl;
        check("pressureOutlet", std::abs(a.pressure-b.pressure),
              tolerances.absoluteTolerance+tolerances.relativeTolerance*std::abs(a.pressure));
        check("zero concentration", b.concentrationNorm, tolerances.concentrationTolerance);
        if (step == numberOfTimeSteps-1) {
            Teuchos::Array<double> velocity(1);
            a.fields[0]->norm2(velocity());
            check("nonzero flow", velocity[0] > 1.e-12 && std::abs(a.pressure) > 1.e-12 ? 0. : 1., 0.);
        }
        // The solver exports the average linear iterations per Newton update.
        // Also compare the reconstructed total number of Krylov iterations.
        const double totalA = fsi.linear[step]*fsi.nonlinear[step];
        const double totalB = fsci.linear[step]*fsci.nonlinear[step];
        if (comm->getRank() == 0)
            std::cout << "Iterations FSI/FSCI: nonlinear " << fsi.nonlinear[step] << "/"
                      << fsci.nonlinear[step] << ", linear total " << totalA << "/" << totalB
                      << ", linear average " << fsi.linear[step] << "/" << fsci.linear[step] << std::endl;
        check("integer nonlinear count FSI", std::abs(fsi.nonlinear[step]-std::round(fsi.nonlinear[step])), 0.);
        check("integer nonlinear count FSCI", std::abs(fsci.nonlinear[step]-std::round(fsci.nonlinear[step])), 0.);
        check("nonlinear iterations", std::abs(fsi.nonlinear[step]-fsci.nonlinear[step]), 0.);
        check("linear iterations (average)", std::abs(fsi.linear[step]-fsci.linear[step]), 1.e-10);
        check("linear iterations (total)", std::abs(totalA-totalB), 1.e-8);
        check("integer linear total FSI", std::abs(totalA-std::round(totalA)), 1.e-8);
        check("integer linear total FSCI", std::abs(totalB-std::round(totalB)), 1.e-8);
        check("nonzero solver work", fsi.nonlinear[step] > 0. && fsci.nonlinear[step] > 0. &&
              totalA > 0. && totalB > 0. ? 0. : 1., 0.);
    }
    int local = passed ? 1 : 0, global = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &local, &global);
    return global != 0;
}
// Restart regression uses the same small artery and assembly path, but with a
// nonzero concentration boundary. This detects lost c_n in the material adapter,
// which a zero-chemistry FSI/FSCI equivalence test cannot expose.
bool checkPassiveRestart(ParameterListPtr_Type shared, const InputFiles& files,
                         bool geometryExplicit, CommConstPtr_Type comm)
{
    const double dt = shared->sublist("Timestepping Parameter").get<double>("dt");
    const std::string directory = std::string("passive_FSCI_") + (geometryExplicit ? "GE_" : "GI_") +
        std::to_string(comm->getSize());
    auto run = [&](bool restart, bool validateOnly = false) {
        auto parameters = rcp(new Teuchos::ParameterList(*shared));
        parameters->sublist("Parameter").set("FSI", false).set("FSCI", true)
            .set("Geometry Explicit", geometryExplicit);
        parameters->sublist("General").set("Preconditioner Method", "FaCSCI")
            .set("Export Data", false).set("Export Suffix", "_" + directory);
        parameters->sublist("Exporter").set("Resume output", false);
        parameters->sublist("Timestepping Parameter").set("Checkpointing", !restart).set("Restart", restart)
            .set("Time step", dt).set("Final time", 2.*dt).set("Number Checkpoints", 1)
            .set("Checkpoint directory", directory).set("Restart directory", directory);
        parameters->sublist("Timestepping Parameter").sublist("Checkpoints").set("1", dt);
        parameters->sublist("Parameter Diffusion").set("D0", .001);
        parameters->sublist("Parameter Solid").sublist("1").set("D0", .001);
        auto components = buildComponentParameters(parameters, files);
        auto domains = buildDomains(parameters, comm);
        auto timeBlocks = buildTimeSteppingDefinition(geometryExplicit, true);
        auto sciParameters = rcp(new Teuchos::ParameterList(*parameters));
        sciParameters->setParameters(*Teuchos::getParametersFromXmlFile(files.solidPreconditioner));
        sciParameters->set("Number of blocks", 2);
        sciParameters->sublist("Parameter").set("SCI", true).set("FSCI", false);
        sciParameters->sublist("General").set("Preconditioner Method", "Monolithic");
        auto chemistry = rcp(new Teuchos::ParameterList(*sciParameters));
        chemistry->set("Number of blocks", 1);
        chemistry->sublist("Parameter") = parameters->sublist("Parameter Diffusion");
        vec2D_dbl_Type diffusion(3, vec_dbl_Type(3, 0.));
        for (int i = 0; i < 3; ++i) diffusion[i][i] = .001;

        FSCIValidationProblem fsci(domains.fluidVelocity, "P2", domains.fluidPressure, "P1",
            domains.solid, "P2", domains.solid, "P2", domains.interfaceDomain, "P2",
            domains.fluidVelocity, "P2", diffusion, noReaction, components.fluid, components.solid,
            chemistry, sciParameters, parameters, components.geometry, timeBlocks);
        if (validateOnly) return std::vector<Frame>{};
        auto boundaries = addMechanicalBoundaries(fsci, domains, parameters, true);
        auto sciBoundaries = rcp(new Boundaries(*boundaries));
        sciBoundaries->addBC(concentrationBoundary, 7, 1, domains.solid, "Dirichlet", 1);
        fsci.getSCIProblem()->addBoundaries(sciBoundaries);
        auto chemBoundaries = rcp(new Boundaries);
        chemBoundaries->addBC(concentrationBoundary, 7, 0, domains.solid, "Dirichlet", 1);
        fsci.getSCIProblem()->problemChem_->addBoundaries(chemBoundaries);
        fsci.initializeProblem();
        fsci.initializeGE();
        fsci.assemble();
        DAESolverInTime<SC,LO,GO,NO> solver(parameters, comm);
        solver.defineTimeStepping(*timeBlocks);
        solver.setProblem(fsci);
        solver.setupTimeStepping();
        solver.advanceInTime();
        return fsci.frames;
    };

    const auto reference = run(false);
    const auto restarted = run(true);
    TEUCHOS_TEST_FOR_EXCEPTION(reference.size() != 2 || restarted.size() != 1,
        std::runtime_error, "Passive FSCI restart must resume at dt and complete the second step.");
    const auto& a = reference.back();
    const auto& b = restarted.back();
    TEUCHOS_TEST_FOR_EXCEPTION(a.concentrationNorm < .01 || b.concentrationNorm < .01 ||
        std::abs(b.time-2.*dt) > 1.e-12, std::runtime_error, "Missing concentration or wrong restart clock.");
    bool passed = std::abs(a.pressure-b.pressure) < 1.e-8;
    for (int i = 0; i < 6; ++i) {
        auto expected = i == 5 ? a.concentration : a.fields[i];
        auto actual = i == 5 ? b.concentration : b.fields[i];
        auto difference = rcp(new Vector(expected->getMap()));
        difference->update(1., *expected, 0.);
        difference->update(-1., *actual, 1.);
        Teuchos::Array<double> error(1), scale(1);
        difference->normInf(error());
        expected->normInf(scale());
        const double tolerance = 1.e-10 + 1.e-7*scale[0];
        passed = passed && std::isfinite(error[0]) && error[0] <= tolerance;
        if (comm->getRank() == 0)
            std::cout << "Passive restart field " << i << ": error=" << error[0]
                      << ", tolerance=" << tolerance << std::endl;
    }

    // Check the coupled manifest before any component restart read: marking any
    // of the three irreversible processes must reject the otherwise valid run.
    const std::string manifest = directory + "/Checkpoint_u_f_p_d_s_lambda_d_f_c_" + std::to_string(dt) + ".xml";
    for (const char* marker : {"Active response took place", "Growth took place", "Reorientation took place"}) {
        checkpoint::onRoot(*comm, [&] {
            auto saved = Teuchos::getParametersFromXmlFile(manifest);
            checkpoint::requirePassiveMaterialRestart(*saved);
            saved->sublist("Material process state").set(marker, true);
            Teuchos::writeParameterListToXmlFile(*saved, manifest);
        });
        bool rejected = false;
        try { run(true, true); }
        catch (const std::runtime_error& error) {
            rejected = std::string(error.what()).find(std::string("Cannot restart: ") + marker) != std::string::npos;
        }
        checkpoint::onRoot(*comm, [&] {
            auto saved = Teuchos::getParametersFromXmlFile(manifest);
            saved->sublist("Material process state").set(marker, false);
            Teuchos::writeParameterListToXmlFile(*saved, manifest);
        });
        passed = passed && rejected;
    }
    return passed;
}
} // namespace

int main(int argc, char* argv[])
{
    // ---------------------------------------------------------------------
    // MPI and command-line parameters (same entry-point pattern as AssembleFE).
    // ---------------------------------------------------------------------
    Tpetra::ScopeGuard scope(&argc, &argv);
    auto comm = Tpetra::getDefaultComm();
    const bool verbose = comm->getRank() == 0;
    InputFiles files;
    // Default to GI: ALE displacement is solved in the coupled Newton system.
    bool geometryExplicit = false;
    bool passiveRestart = false;
    Teuchos::CommandLineProcessor commandLine;
    commandLine.setOption("problemfile", &files.problem, "Shared physical parameters and validation tolerances.");
    commandLine.setOption("solverfile", &files.solver, "Shared coupled linear solver parameters.");
    commandLine.setOption("precfileFluid", &files.fluidPreconditioner, "Fluid preconditioner parameters.");
    commandLine.setOption("precfileStructure", &files.solidPreconditioner, "Solid/SCI preconditioner parameters.");
    commandLine.setOption("precfileGeometry", &files.geometryPreconditioner, "ALE preconditioner parameters.");
    commandLine.setOption("solverfileGeometry", &files.geometrySolver, "ALE linear solver parameters.");
    commandLine.setOption("geometry-explicit", "geometry-implicit", &geometryExplicit,
                          "Compare both formulations with explicit or implicit ALE geometry.");
    commandLine.setOption("passive-restart", "equivalence", &passiveRestart,
                          "Compare passive FSCI restart with uninterrupted nonzero-chemistry evolution.");
    commandLine.throwExceptions(false);
    const auto parsed = commandLine.parse(argc, argv);
    if (parsed == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return EXIT_SUCCESS;
    if (parsed != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) return EXIT_FAILURE;

    try {
        // -----------------------------------------------------------------
        // Shared inputs. These lists remain untouched by either simulation;
        // each setup below starts from its own copy of the same parameters.
        // -----------------------------------------------------------------
        auto parameterListAll = Teuchos::getParametersFromXmlFile(files.problem);
        parameterListAll->setParameters(*Teuchos::getParametersFromXmlFile(files.solver));
        if (passiveRestart) {
            const bool passed = checkPassiveRestart(parameterListAll, files, geometryExplicit, comm);
            if (verbose) std::cout << "FSCI passive restart: " << (passed ? "PASS" : "FAIL") << std::endl;
            return passed ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        const auto& problemSettings = parameterListAll->sublist("Parameter");
        const auto& timeSettings = parameterListAll->sublist("Timestepping Parameter");
        const auto& material = parameterListAll->sublist("Parameter Solid").sublist("1");
        const int dim = problemSettings.get<int>("Dimension");
        const std::string discType = problemSettings.get<std::string>("Discretization");
        const double dt = timeSettings.get<double>("dt");
        TEUCHOS_TEST_FOR_EXCEPTION(dim != 3 || discType != "P2", std::logic_error,
            "This native-material comparison requires 3D P2 displacement/velocity and P1 pressure.");
        TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(dt) || dt <= 0. ||
            std::abs(timeSettings.get<double>("Final time")-numberOfTimeSteps*dt) > 1.e-12 ||
            timeSettings.get<bool>("Restart") || timeSettings.get<bool>("Checkpointing"), std::logic_error,
            "Validation requires exactly two timesteps from zero, without restart/checkpointing.");
        TEUCHOS_TEST_FOR_EXCEPTION(problemSettings.get<bool>("Chemistry Explicit") ||
            parameterListAll->sublist("Parameter Diffusion").get<double>("D0") != 0. ||
            material.get<double>("D0") != 0. || material.get<double>("ActiveBool") != 0. ||
            material.get<double>("growthBool") != 0. || material.get<double>("reorientationBool") != 0.,
            std::logic_error, "Validation requires implicit, zero-diffusion chemistry and a passive material.");
        TEUCHOS_TEST_FOR_EXCEPTION(!parameterListAll->sublist("General").get<bool>("Export Data"), std::logic_error,
            "Validation needs Export Data=true to compare iteration counts from this run.");
        const auto& validation = parameterListAll->sublist("Validation");
        ComparisonSettings tolerances{validation.get<double>("Absolute tolerance"),
            validation.get<double>("Relative tolerance"), validation.get<double>("Concentration tolerance")};
        for (double tolerance : {tolerances.absoluteTolerance, tolerances.relativeTolerance, tolerances.concentrationTolerance})
            TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(tolerance) || tolerance < 0., std::logic_error,
                "Validation tolerances must be finite and nonnegative.");
        const std::string runSuffix = std::string(geometryExplicit ? "GE_" : "GI_") + std::to_string(comm->getSize());
        if (verbose)
            std::cout << std::setprecision(12) << "FSI/FSCI validation: " << runSuffix << ", "
                      << numberOfTimeSteps << " steps, dt=" << dt << std::endl;

        Result fsiResult, fsciResult;
        // -----------------------------------------------------------------
        // FSI PROBLEM: mechanics-only reference, with FaCSI preconditioning.
        // -----------------------------------------------------------------
        {
            if (verbose) std::cout << "\n========== FSI: setup and solve ==========\n";
            auto parameterListFSI = rcp(new Teuchos::ParameterList(*parameterListAll));
            parameterListFSI->sublist("Parameter").set("FSI", true).set("FSCI", false)
                .set("Geometry Explicit", geometryExplicit);
            parameterListFSI->sublist("General").set("Preconditioner Method", "FaCSI")
                .set("Export Suffix", "_FSI_" + runSuffix);
            auto components = buildComponentParameters(parameterListFSI, files);
            auto domains = buildDomains(parameterListFSI, comm);
            auto timeBlocks = buildTimeSteppingDefinition(geometryExplicit, false);

            FSIValidationProblem fsi(domains.fluidVelocity, discType,
                domains.fluidPressure, "P1", domains.solid, discType,
                domains.interfaceDomain, discType, domains.fluidVelocity, discType,
                components.fluid, components.solid, parameterListFSI, components.geometry, timeBlocks);
            addMechanicalBoundaries(fsi, domains, parameterListFSI);
            fsiResult = advanceTwoSteps(fsi, parameterListFSI, timeBlocks, comm);
        }

        // -----------------------------------------------------------------
        // FSCI PROBLEM: the same mechanics plus a zero-concentration equation.
        // FaCSCI preconditions the additional solid/chemistry block. With no
        // diffusion, reaction or concentration input, c must remain zero.
        // -----------------------------------------------------------------
        {
            if (verbose) std::cout << "\n========== FSCI: setup and solve ==========\n";
            auto parameterListFSCI = rcp(new Teuchos::ParameterList(*parameterListAll));
            parameterListFSCI->sublist("Parameter").set("FSI", false).set("FSCI", true)
                .set("Geometry Explicit", geometryExplicit);
            parameterListFSCI->sublist("General").set("Preconditioner Method", "FaCSCI")
                .set("Export Suffix", "_FSCI_" + runSuffix);
            auto components = buildComponentParameters(parameterListFSCI, files);
            auto domains = buildDomains(parameterListFSCI, comm);
            auto timeBlocks = buildTimeSteppingDefinition(geometryExplicit, true);

            // SCI's local blocks are [d_s, c], independently of their positions
            // in the full FSCI system. Reuse the solid preconditioner settings
            // and add its scalar block, with the same overlap/solver options.
            auto parameterListSCI = rcp(new Teuchos::ParameterList(*parameterListFSCI));
            parameterListSCI->setParameters(*Teuchos::getParametersFromXmlFile(files.solidPreconditioner));
            parameterListSCI->set("Number of blocks", 2);
            parameterListSCI->sublist("Parameter").set("SCI", true).set("FSCI", false);
            parameterListSCI->sublist("General").set("Preconditioner Method", "Monolithic");
            auto parameterListChem = rcp(new Teuchos::ParameterList(*parameterListSCI));
            parameterListChem->set("Number of blocks", 1);
            parameterListChem->sublist("Parameter") = parameterListFSCI->sublist("Parameter Diffusion");
            vec2D_dbl_Type diffusionTensor(dim, vec_dbl_Type(dim, 0.));

            // Wall chemistry lives on the same P2 domain as the solid.
            FSCIValidationProblem fsci(domains.fluidVelocity, discType,
                domains.fluidPressure, "P1", domains.solid, discType,
                domains.solid, discType, domains.interfaceDomain, discType,
                domains.fluidVelocity, discType, diffusionTensor, noReaction,
                components.fluid, components.solid, parameterListChem, parameterListSCI,
                parameterListFSCI, components.geometry, timeBlocks);
            auto solidBoundaries = addMechanicalBoundaries(fsci, domains, parameterListFSCI);
            fsci.getSCIProblem()->addBoundaries(rcp(new Boundaries(*solidBoundaries)));
            // No concentration Dirichlet condition, flux, reaction or source.
            // The transient chemical mass term keeps the zero initial state.
            fsci.getSCIProblem()->problemChem_->addBoundaries(rcp(new Boundaries));
            fsciResult = advanceTwoSteps(fsci, parameterListFSCI, timeBlocks, comm);
        }

        // -----------------------------------------------------------------
        // Compare accepted solutions, outlet loads and solver work at each step.
        // Results own deep copies, so destroying the first problem is harmless.
        // -----------------------------------------------------------------
        if (verbose) std::cout << "\n========== FSI/FSCI comparison ==========\n";
        const bool passed = compare(fsiResult, fsciResult, tolerances, comm);
        if (verbose) std::cout << "FSI_FSCI_Validation_Test: " << (passed ? "PASS" : "FAIL") << std::endl;
        return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "FSI_FSCI_Validation_Test, rank " << comm->getRank() << ": " << error.what() << std::endl;
        return EXIT_FAILURE;
    }
}
