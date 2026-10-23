#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/Checkpointing/CheckpointMetadata.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_StackedTimer.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace FEDD;
using Teuchos::RCP;
using Teuchos::rcp;
using DomainType = Domain<default_sc,default_lo,default_go,default_no>;
using NavierStokesType = NavierStokes<default_sc,default_lo,default_go,default_no>;
using TimeSolver = DAESolverInTime<default_sc,default_lo,default_go,default_no>;

namespace {

void zero2D(double*, double* result, double, const double*)
{
    result[0] = result[1] = 0.;
}

void zero3D(double*, double* result, double, const double*)
{
    result[0] = result[1] = result[2] = 0.;
}

void inflow2D(double* x, double* result, double, const double* parameters)
{
    // The structured BFS inlet occupies y in [0,1].
    result[0] = 4. * parameters[0] * x[1] * (1. - x[1]);
    result[1] = 0.;
}

void inflow3D(double* x, double* result, double, const double* parameters)
{
    // In 3D the inlet occupies y,z in [0,1]; the step is below z=0.
    result[0] = 16. * parameters[0] * x[1] * (1. - x[1]) * x[2] * (1. - x[2]);
    result[1] = result[2] = 0.;
}

/** @brief Hide normal library progress output during measurement, restoring it on exit. */
class QuietOutput {
public:
    explicit QuietOutput(bool quiet) : original_(quiet ? std::cout.rdbuf(sink_.rdbuf()) : nullptr) {}
    ~QuietOutput() { if (original_) std::cout.rdbuf(original_); }
private:
    Teuchos::oblackholestream sink_;
    std::streambuf* original_;
};

/** @brief Add a benchmark region above the nested library TimeMonitor counters.
 * These regions introduce no synchronization. If TimeMonitor disables stacked
 * timing after detecting overlapping timers, leave the damaged stack alone;
 * the benchmark subsequently reports the error instead of publishing a report.
 */
class StackedRegion {
public:
    explicit StackedRegion(const std::string& name)
        : timer_(Teuchos::TimeMonitor::getStackedTimer()), name_(name)
    {
        if (!timer_.is_null()) timer_->start(name_);
    }
    ~StackedRegion()
    {
        if (!timer_.is_null() && Teuchos::TimeMonitor::getStackedTimer() == timer_)
            timer_->stop(name_);
    }
private:
    RCP<Teuchos::StackedTimer> timer_;
    std::string name_;
};

/** @brief Attach one stacked timer for this invocation and restore the prior binding.
 * Checkpoint and solver TimeMonitor scopes automatically populate this stack.
 * MPI aggregation/reporting happens collectively after stopping the root timer;
 * writing the aggregated report is a rank-zero operation.
 */
class StackedTiming {
public:
    StackedTiming() : previous_(Teuchos::TimeMonitor::getStackedTimer()),
                     timer_(rcp(new Teuchos::StackedTimer("BFS restart performance")))
    {
        Teuchos::TimeMonitor::setStackedTimer(timer_);
    }
    ~StackedTiming() { Teuchos::TimeMonitor::setStackedTimer(previous_); }

    void check() const
    {
        TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::TimeMonitor::getStackedTimer() != timer_,
                                   std::runtime_error, "Stacked timing was disabled by a library timer overlap; no valid stacked report can be produced.");
    }

    void report(RCP<const Teuchos::Comm<int>> comm, const std::filesystem::path& directory)
    {
        check();
        timer_->stop("BFS restart performance");
        Teuchos::StackedTimer::OutputOptions options;
        options.output_fraction = options.output_minmax = options.output_proc_minmax = true;
        options.align_columns = true;
        std::ostringstream output;
        timer_->report(output, comm, options);
        checkpoint::onRoot(*comm, [&] {
            std::ofstream file(directory / "stacked_timings.txt");
            file << output.str();
            file.close();
            if (!file) throw std::runtime_error("Cannot write stacked_timings.txt.");
            std::cout << '\n' << output.str();
        });
    }

private:
    RCP<Teuchos::StackedTimer> previous_, timer_;
};

/** @brief Execute a stage inside a stacked timer, synchronizing its boundaries. */
template<class Operation>
void measure(const std::string& name, const Teuchos::Comm<int>& comm, Operation operation)
{
    comm.barrier();
    {
        StackedRegion region(name);
        operation();
    }
    comm.barrier();
}

/** @brief Execute a Navier-Stokes trajectory with nested stage and library timers.
 * Each phase constructs its own mesh, problem and time solver, so startup costs
 * are reported separately from I/O and continuation. No reference is loaded and
 * no solution-error comparison is performed.
 */
void runPhase(ParameterListPtr_Type parameters, RCP<const Teuchos::Comm<int>> comm,
              int dim, int subdivisions, int ratio, int length, bool quiet)
{
    QuietOutput silence(quiet);
    RCP<DomainType> domain;
    measure("mesh_build", *comm, [&] {
        const std::vector<double> origin = dim == 2
            ? std::vector<double>{-1., -1.} : std::vector<double>{-1., 0., -1.};
        domain = dim == 2 ? rcp(new DomainType(origin, length + 1., 2., comm))
                         : rcp(new DomainType(origin, length + 1., 1., 2., comm));
        domain->buildMesh(2, "BFS", dim, "P1", subdivisions, ratio);
    });

    auto boundaries = rcp(new BCBuilder<default_sc,default_lo,default_go,default_no>());
    std::vector<double> inlet = {parameters->sublist("Parameter").get<double>("MaxVelocity")};
    auto zero = dim == 2 ? zero2D : zero3D;
    boundaries->addBC(zero, 1, 0, domain, "Dirichlet", dim);
    boundaries->addBC(zero, 4, 0, domain, "Dirichlet", dim);
    boundaries->addBC(dim == 2 ? inflow2D : inflow3D, 2, 0, domain, "Dirichlet", dim, inlet);

    RCP<NavierStokesType> problem = rcp(new NavierStokesType(domain, "P1", domain, "P1", parameters));
    problem->addBoundaries(boundaries);
    measure("problem_initialization", *comm, [&] { problem->initializeProblem(); });
    measure("spatial_assembly", *comm, [&] {
        problem->assemble();
        problem->setBoundariesRHS();
    });
    SmallMatrix<int> timeBlocks(2);
    timeBlocks[0][0] = timeBlocks[0][1] = 1;
    RCP<TimeSolver> solver = rcp(new TimeSolver(parameters, comm));
    measure("time_setup", *comm, [&] {
        solver->defineTimeStepping(timeBlocks);
        solver->setProblem(*problem);
        solver->setupTimeStepping();
    });
    measure("advance", *comm, [&] { solver->advanceInTime(); });
    // Preserve destruction order and close files before the next phase reads them.
    // files_close is nested here; other destruction is the teardown remainder.
    measure("teardown", *comm, [&] {
        solver = Teuchos::null;
        problem = Teuchos::null;
        boundaries = Teuchos::null;
        domain = Teuchos::null;
    });
}

} // namespace

/** @brief Structured BFS weak-scaling benchmark of Navier-Stokes checkpoint operations. */
int main(int argc, char* argv[])
{
    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession session(&argc, &argv, &blackhole);
    auto comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    std::string problemFile = "parametersProblem.xml", precFile = "parametersPrec.xml";
    std::string solverFile = "parametersSolver.xml", outputDirectory = "performance_output";
    std::string benchmarkFile= "parametersBenchmark.xml", aspect = "all";
    int dim = 2, ratio = 3, length = 4, steps = 4, continuationSteps = 2, repetitions = 1;
    bool quiet = true, baseline = true, resumeOutput = true, initialSolution = true;
    Teuchos::CommandLineProcessor commandLine;
    commandLine.setOption("benchmarkfile", &benchmarkFile, "XML benchmark controls; explicit command-line options override its values.");
    commandLine.setOption("aspect", &aspect, "all, baseline, checkpoint_write, restart, initial_solution or resume_output.");
    commandLine.setOption("problemfile", &problemFile, "Physics and time integration parameters.");
    commandLine.setOption("precfile", &precFile, "Preconditioner settings (replace for cluster runs if needed).");
    commandLine.setOption("solverfile", &solverFile, "Linear solver settings.");
    commandLine.setOption("dimension", &dim, "Structured BFS dimension, 2 or 3.");
    commandLine.setOption("H-h", &ratio, "H/h: mesh subdivisions per subdomain edge (fixed for weak scaling).");
    commandLine.setOption("length", &length, "Downstream BFS length; rank count must be (2*length+1)*k^dimension.");
    commandLine.setOption("steps", &steps, "Producer timesteps; its final accepted state is checkpointed.");
    commandLine.setOption("restart-steps", &continuationSteps, "Timesteps in each continuation / initial-solution phase.");
    commandLine.setOption("repetitions", &repetitions, "Independent samples; each keeps its own output directory.");
    commandLine.setOption("output-directory", &outputDirectory, "New directory for checkpoints, exports and stacked_timings.txt.");
    commandLine.setOption("quiet", "verbose", &quiet, "Only print timing output (default), or also library progress.");
    commandLine.setOption("baseline", "no-baseline", &baseline, "Also measure a fresh run without checkpoint/output I/O.");
    commandLine.setOption("resume-output", "no-resume-output", &resumeOutput, "Also measure continuation of ParaView/text output.");
    commandLine.setOption("initial-solution", "no-initial-solution", &initialSolution, "Also measure loading primary fields at time zero.");
    commandLine.throwExceptions(false);
    const auto parsed = commandLine.parse(argc, argv);
    if (parsed == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return EXIT_SUCCESS;
    if (parsed != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) return EXIT_FAILURE;

    try {
        if (!benchmarkFile.empty()) {
            auto controls = Teuchos::getParametersFromXmlFile(benchmarkFile);
            dim = controls->get("Dimension", dim);
            ratio = controls->get("H/h", ratio);
            length = controls->get("BFS length", length);
            steps = controls->get("Steps", steps);
            continuationSteps = controls->get("Restart steps", continuationSteps);
            repetitions = controls->get("Repetitions", repetitions);
            aspect = controls->get("Aspect", aspect);
            quiet = controls->get("Quiet", quiet);
            outputDirectory = controls->get("Output directory", outputDirectory);
            // Input references are relative to the benchmark XML; output remains
            // relative to the launch directory (typically shared cluster scratch).
            const auto inputDirectory = std::filesystem::absolute(benchmarkFile).parent_path();
            const auto inputPath = [&](const char* key, const std::string& fallback) {
                return (inputDirectory / controls->get(key, fallback)).lexically_normal().string();
            };
            problemFile = inputPath("Problem file", problemFile);
            precFile = inputPath("Preconditioner file", precFile);
            solverFile = inputPath("Solver file", solverFile);
            // Teuchos only changes options explicitly present on the command line.
            // Reparse to apply those overrides to the defaults just read from XML.
            if (commandLine.parse(argc, argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
                return EXIT_FAILURE;
        }
        bool checkpointRun = true, restartRun = true;
        if (aspect != "all") {
            TEUCHOS_TEST_FOR_EXCEPTION(aspect != "baseline" && aspect != "checkpoint_write" && aspect != "restart" &&
                                       aspect != "initial_solution" && aspect != "resume_output", std::logic_error,
                                       "Unknown benchmark Aspect '" << aspect << "'. Use all, baseline, checkpoint_write, restart, initial_solution or resume_output.");
            baseline = aspect == "baseline";
            checkpointRun = !baseline;
            restartRun = aspect == "restart";
            resumeOutput = aspect == "resume_output";
            initialSolution = aspect == "initial_solution";
        }
        TEUCHOS_TEST_FOR_EXCEPTION((dim != 2 && dim != 3) || ratio < 2 || length < 1 || steps < 2 ||
                                   continuationSteps < 1 || repetitions < 1, std::logic_error,
                                   "Require dimension 2/3, H/h >= 2, length >= 1, steps >= 2, restart-steps >= 1 and repetitions >= 1.");
        const long long blocks = 2LL * length + 1;
        const int subdivisions = static_cast<int>(std::llround(std::pow(comm->getSize() / double(blocks), 1. / dim)));
        long long required = blocks;
        for (int i = 0; i < dim; ++i) required *= subdivisions;
        TEUCHOS_TEST_FOR_EXCEPTION(subdivisions < 1 || required != comm->getSize(), std::logic_error,
                                   "Structured BFS requires ranks = (2*length+1)*k^dimension for integer k >= 1. "
                                   "For length=4: 2D ranks 9,36,81,144,...; 3D ranks 9,72,243,576,...");
        auto base = Teuchos::getParametersFromXmlFile(problemFile);
        base->setParameters(*Teuchos::getParametersFromXmlFile(precFile));
        base->setParameters(*Teuchos::getParametersFromXmlFile(solverFile));
        auto& time = base->sublist("Timestepping Parameter");
        const double dt = time.get<double>("dt");
        TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(dt) || dt <= 0. || time.get("Class", "Multistep") != "Multistep" ||
                                   time.get("BDF", 2) != 2, std::logic_error,
                                   "This benchmark requires fixed positive dt and BDF2 multistep integration.");
        time.set("Restart", false).set("Initial solution", false).set("Failure recovery", false);
        time.set("Timestepping type", "non-adaptive");
        time.sublist("Timestepping Intervalls").set("Number of Segments", 0);
        base->sublist("Parameter").set("Dimension", dim).set("Cancel MaxNonLinIts", true);
        base->sublist("General").set("Checkpoint timings", true).set("Safe all solution", false);
        base->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").set("DofsPerNode1", dim);
        // Record the resolved run controls with the physical/solver parameters.
        auto& benchmark = base->sublist("Benchmark");
        benchmark.set("Dimension", dim).set("H/h", ratio).set("BFS length", length)
            .set("Steps", steps).set("Restart steps", continuationSteps).set("Repetitions", repetitions)
            .set("Aspect", aspect).set("Quiet", quiet)
            .set("MPI ranks", comm->getSize()).set("Subdomains per unit", subdivisions)
            .set("Output directory", std::filesystem::absolute(outputDirectory).string())
            .set("Problem file", std::filesystem::absolute(problemFile).string())
            .set("Preconditioner file", std::filesystem::absolute(precFile).string())
            .set("Solver file", std::filesystem::absolute(solverFile).string())
            .set("Run baseline", baseline).set("Run checkpoint producer", checkpointRun)
            .set("Run restart", restartRun).set("Run output continuation", resumeOutput)
            .set("Run initial solution", initialSolution);
        const auto originalDirectory = std::filesystem::current_path();
        const auto root = std::filesystem::absolute(outputDirectory);
        checkpoint::onRoot(*comm, [&] {
            if (std::filesystem::exists(root))
                throw std::runtime_error("Output directory already exists; choose a new --output-directory to preserve previous measurements.");
            std::filesystem::create_directories(root);
        });
        comm->barrier();
        checkpoint::onRoot(*comm, [&] {
            Teuchos::writeParameterListToXmlFile(*base, (root / "parametersEffective.xml").string());
        });
        StackedTiming stackedTiming;
        const double checkpointTime = steps * dt;
        for (int repeat = 1; repeat <= repetitions; ++repeat) {
            StackedRegion repetition("repeat_" + std::to_string(repeat));
            const auto repeatDirectory = root / ("repeat_" + std::to_string(repeat));
            const auto producerDirectory = repeatDirectory / "checkpoint_run";
            const auto checkpointDirectory = producerDirectory / "checkpoints";
            const auto execute = [&](const std::string& phase, bool restart, bool initial, bool output, bool write) {
                const auto directory = phase == "resume_output" ? producerDirectory : repeatDirectory / phase;
                checkpoint::onRoot(*comm, [&] { std::filesystem::create_directories(directory); });
                comm->barrier();
                std::filesystem::current_path(directory);
                auto parameters = rcp(new Teuchos::ParameterList(*base));
                auto& phaseTime = parameters->sublist("Timestepping Parameter");
                const double finalTime = restart ? checkpointTime + continuationSteps * dt
                    : initial ? continuationSteps * dt : checkpointTime;
                phaseTime.set("Restart", restart).set("Initial solution", initial)
                    .set("Restart directory", checkpointDirectory.string()).set("Time step", restart ? checkpointTime : 0.)
                    .set("Initial solution directory", checkpointDirectory.string()).set("Initial solution time", checkpointTime)
                    .set("Final time", finalTime).set("Checkpointing", write).set("Number Checkpoints", 1)
                    .set("Checkpoint directory", (directory / (restart ? "continued_checkpoints" : "checkpoints")).string());
                phaseTime.sublist("Checkpoints").set("1", finalTime);
                parameters->sublist("General").set("ParaViewExport", output).set("Export Data", output);
                parameters->sublist("Exporter").set("Resume output", phase == "resume_output").set("Keep old output", false);
                checkpoint::onRoot(*comm, [&] {
                    Teuchos::writeParameterListToXmlFile(*parameters, (directory / (phase + ".xml")).string());
                });
                {
                    StackedRegion phaseRegion(phase);
                    runPhase(parameters, comm, dim, subdivisions, ratio, length, quiet);
                }
                stackedTiming.check();
                std::filesystem::current_path(originalDirectory);
            };
            if (baseline) execute("baseline", false, false, false, false);
            if (checkpointRun) execute("checkpoint_run", false, false, resumeOutput, true);
            if (restartRun) execute("restart", true, false, false, false);
            if (resumeOutput) execute("resume_output", true, false, true, true);
            if (initialSolution) execute("initial_solution", false, true, false, false);
        }
        stackedTiming.report(comm, root);
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error) {
        if (comm->getRank() == 0) std::cerr << "BFS performance benchmark failed: " << error.what() << std::endl;
        return EXIT_FAILURE;
    }
}
