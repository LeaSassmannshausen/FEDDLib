#include "feddlib/problems/Solver/TimeSteppingTools.hpp"
#include "feddlib/core/Checkpointing/CheckpointMetadata.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <iostream>

using namespace FEDD;

// Check clocks and BDF2 against polynomial derivatives, independently of PDE solves.
int main(int argc, char** argv)
{
    Teuchos::GlobalMPISession session(&argc, &argv);
    auto comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    auto settings = Teuchos::rcp(new Teuchos::ParameterList);
    settings->set("Class", "Multistep").set("BDF", 2).set("dt", 0.003).set("Final time", 0.012);
    settings->sublist("Timestepping Intervalls").set("Number of Segments", 2);
    settings->sublist("Timestepping Intervalls").sublist("1").set("Start Time", 0.005).set("dt", 0.0015);
    settings->sublist("Timestepping Intervalls").sublist("2").set("Start Time", 0.0095).set("dt", 0.002);
    TimeSteppingTools clock(settings, comm);
    int failures = 0;
    auto require = [&](bool ok, const char* message) {
        if (!ok) { ++failures; if (comm->getRank() == 0) std::cerr << "FAIL " << message << '\n'; }
    };
    const auto requireHistory = [&](const TimeSteppingTools& tool, const std::vector<double>& expected,
                                    const char* message) {
        bool matches = tool.historyTimes_.size() == 3;
        for (int i = 0; matches && i < expected.size(); ++i)
            matches = std::abs(tool.historyTimes_[i] - expected[i]) < 1.e-14;
        require(matches, message);
    };
    const std::vector<double> ends = {0.003, 0.005, 0.0065, 0.008, 0.0095, 0.0115, 0.012};
    std::vector<double> acceptedTimes = {0.};
    double previousTime = -clock.get_dt_prev();
    for (double end : ends) {
        clock.prepareStep();
        const double t = clock.currentTime(), h = clock.get_dt();
        require(std::abs(t + h - end) < 1.e-14, "land on interval/final boundary before assembly");
        const double derivative = (clock.getInformationBDF(0) * end * end
            - clock.getInformationBDF(2) * t * t
            - clock.getInformationBDF(3) * previousTime * previousTime) / h;
        require(std::abs(derivative - 2. * end) < 1.e-13, "variable BDF2 differentiates a quadratic");
        previousTime = t;
        clock.advanceTime();
        require(std::abs(clock.get_dt_prev() - h) < 1.e-14, "retain actual completed increment");
        acceptedTimes.insert(acceptedTimes.begin(), end);
        acceptedTimes.resize(std::min<std::size_t>(3, acceptedTimes.size()));
        requireHistory(clock, acceptedTimes, "retain every available timestamp across interval and final-step changes");
    }
    require(clock.step_ == 7, "step counter independent of time/dt");
    auto resumed = Teuchos::rcp(new Teuchos::ParameterList);
    resumed->set("Class", "Multistep").set("BDF", 2).set("dt", 0.0015)
        .set("Final time", 0.02).set("Restart", true).set("Time step", 0.012);
    resumed->sublist("_Restart time state") = settings->sublist("_Checkpoint time state");
    TimeSteppingTools restart(resumed, comm);
    require(std::abs(restart.get_dt_prev() - 0.0005) < 1.e-14, "restore shortened incoming increment");
    require(std::abs(restart.getInformationBDF(0) - 1.75) < 1.e-13, "use new dt / saved dt_prev");
    require(restart.step_ == 7, "restore saved timestep number");
    const std::vector<double> savedHistory = {0.012, 0.0115, 0.0095};
    requireHistory(restart, savedHistory, "restore all three nonuniform timestamps independently of new dt");
    auto all = Teuchos::rcp(new Teuchos::ParameterList);
    all->sublist("Timestepping Parameter") = *resumed;
    auto schema = checkpoint::makeSchema(*all, {});
    auto description = checkpoint::atTime(schema, 0.012, &resumed->sublist("_Checkpoint time state"));
    for (int i = 0; i < savedHistory.size(); ++i)
        require(std::abs(description.sublist("Time state").sublist("History times").get<double>(std::to_string(i)) - savedHistory[i]) < 1.e-14,
                "metadata retains every actual history timestamp rather than subtracting new dt");

    // Clocks created before problem initialization must not invent restart history,
    // even when legacy reading is allowed. The validator supplies the legacy clock.
    auto pendingSettings = Teuchos::rcp(new Teuchos::ParameterList(*resumed));
    pendingSettings->remove("_Restart time state");
    pendingSettings->set("Allow legacy restart", true);
    TimeSteppingTools pending(pendingSettings, comm);
    require(pending.historyTimes_.empty() && !pendingSettings->isSublist("_Checkpoint time state"),
            "pending restart publishes no fabricated clock");
    bool rejectedCoefficients = false, rejectedAdvance = false;
    try { pending.getInformationBDF(0); }
    catch (const std::logic_error&) { rejectedCoefficients = true; }
    try { pending.advanceTime(); }
    catch (const std::logic_error&) { rejectedAdvance = true; }
    require(rejectedCoefficients && rejectedAdvance, "restart cannot proceed before its clock is validated");
    pendingSettings->sublist("_Restart time state") = settings->sublist("_Checkpoint time state");
    pending.prepareStep();
    requireHistory(pending, savedHistory, "restore history when validation finishes after clock construction");

    restart.advanceTime();
    restart.prepareStep();
    requireHistory(restart, {0.0135, 0.012, 0.0115}, "continuation shifts actual history without reloading old timestamps");

    // A fresh clock uses its selected first increment, rather than the default dt.
    auto freshSettings = Teuchos::rcp(new Teuchos::ParameterList);
    freshSettings->set("Class", "Multistep").set("BDF", 2).set("dt", 0.003).set("Final time", 0.004);
    freshSettings->sublist("Timestepping Intervalls").set("Number of Segments", 1);
    freshSettings->sublist("Timestepping Intervalls").sublist("1").set("Start Time", 0.).set("dt", 0.0015);
    TimeSteppingTools fresh(freshSettings, comm);
    require(std::abs(fresh.get_dt_prev() - 0.0015) < 1.e-14, "initialize from the interval active at time zero");
    requireHistory(fresh, {0., -0.0015, -0.003}, "only fresh startup creates negative history placeholders");

    resumed->set("Restart", false);
    restart.setParameter();
    require(restart.currentTime() == 0. && restart.step_ == 0, "fresh reinitialization resets the restart clock");
    requireHistory(restart, {0., -0.0015, -0.003}, "fresh reinitialization ignores old saved restart history");
    if (comm->getRank() == 0) std::cout << (failures ? "FAIL" : "PASS") << " variable timestep clocks and BDF coefficients\n";
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
