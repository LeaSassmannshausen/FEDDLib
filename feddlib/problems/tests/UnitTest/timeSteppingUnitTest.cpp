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
    const std::vector<double> ends = {0.003, 0.005, 0.0065, 0.008, 0.0095, 0.0115, 0.012};
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
    auto all = Teuchos::rcp(new Teuchos::ParameterList);
    all->sublist("Timestepping Parameter") = *resumed;
    auto schema = checkpoint::makeSchema(*all, {});
    auto description = checkpoint::atTime(schema, 0.012, &resumed->sublist("_Checkpoint time state"));
    require(std::abs(description.sublist("Time state").sublist("History times").get<double>("1") - 0.0115) < 1.e-14,
            "metadata retains actual history time rather than subtracting new dt");
    if (comm->getRank() == 0) std::cout << (failures ? "FAIL" : "PASS") << " variable timestep clocks and BDF coefficients\n";
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
