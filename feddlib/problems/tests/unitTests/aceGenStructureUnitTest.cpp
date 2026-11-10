#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include <Tpetra_Core.hpp>
#include <cmath>
#include <iostream>
#include <aceinterface.hpp>
#include <filesystem>

using namespace FEDD;
using LO = default_lo;

void require(bool condition, const char* message)
{
    TEUCHOS_TEST_FOR_EXCEPTION(!condition, std::runtime_error, message);
}

void setPassiveMaterial(const ParameterListPtr_Type& parameters, int flag)
{
    // Native model defaults, in getDomainDataNames() order. Keep this fixture
    // independent of benchmark inputs and supply the current numbered material.
    const double defaults[] = {30., 1.2, .9, .0443, -.0443, 50.,
        .5131, .1624, .4, .134, .00166, .000066, 148.262, .001006,
        11.857, 26.68, 151.73775, 2.75662, 11.52507, 1.27631, 3.08798,
        1.5, 1., -.0010694, 100., .0009735, 200., -.0002323, .0000699,
        50., 1000., -.00001, .07797910389823641, .13677878297366072,
        .6673453235597163, .00006, 0., 1000., .6, 5., 1., 1., 3., 3.,
        1.000882, 1.234826, 1.141419, .98, .98, .98, .15, 0., 0., 0.,
        0., .13333333333333333, 2.958083520579672, 1., 1., 2.};
    AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10 model(18);
    require(model.getNumberOfDomainData() == sizeof(defaults)/sizeof(double), "Material parameter fixture is outdated");
    auto& material = parameters->sublist("Parameter Solid").sublist("1");
    material.set("Volume Flag", flag);
    for (int i = 0; i < model.getNumberOfDomainData(); ++i) {
        const std::string name = model.getDomainDataNames()[i];
        const auto first = name.find('-') + 1;
        material.set(name.substr(first, name.find('_') - first), defaults[i]);
    }
}

int main(int argc, char** argv)
{
    Tpetra::ScopeGuard scope(&argc, &argv);
    auto comm = Tpetra::getDefaultComm();
    auto domain = Teuchos::rcp(new Domain<>(std::vector<double>{0., 0., 0.}, 1., 1., 1., comm));
    domain->buildMesh(1, "Square", 3, "P2", 1, 1, 0);
    // AceGen uses signed integration weights; the structured generator has
    // mixed orientations. Keep this arterial-material fixture off its axis.
    domain->preProcessMesh(false, true);
    for (auto points : {domain->getPointsUnique(), domain->getPointsRepeated()})
        for (auto& point : *points) {
            point[0] += 2.;
            point[1] += 1.;
        }
    for (const std::string model : {"Neo-Hooke", "SCI_NH", "SCI_SMC_MLCK", "SCI_SMC_Active_Growth_Reorientation"}) {
        const std::string structureModel = model == "Neo-Hooke" ? "SCI_NH" : model;
        auto parameters = Teuchos::rcp(new Teuchos::ParameterList);
        parameters->sublist("Parameter").set("Dimension", 3).set("Material model", model);
        if (model == "Neo-Hooke")
            parameters->sublist("Parameter").set("FSCI", true).set("Structure Model", structureModel);
        parameters->sublist("Timestepping Parameter").set("Class", std::string("Multistep"))
            .set("dt", 0.01).set("Final time", 0.1);
        auto& intervals = parameters->sublist("Timestepping Parameter").sublist("Timestepping Intervalls");
        intervals.set("Number of Segments", 1);
        intervals.sublist("1").set("Start Time", 0.01).set("dt", 0.02);
        // Explicit parameters for the newer AceGen model; passive response only.
        parameters->sublist("Parameter Solid").set("Gamma1Min", 0.07797910389823641)
            .set("Gamma3Min", 0.13677878297366072).set("ActiveBool", 0.)
            .set("growthBool", 0.).set("reorientationBool", 0.).set("growthDirectionBool", 0.)
            .set("kFiberDistribution", 0.13333333333333333).set("P", 2.958083520579672)
            .set("PostProcessingPower", 2.);
        if (structureModel == "SCI_SMC_Active_Growth_Reorientation") {
            setPassiveMaterial(parameters, domain->getElementsC()->getElement(0).getFlag());
            parameters->sublist("Parameter Solid").sublist("1").set("ActiveStartTime", .015);
        }
        NonLinElasticity<> problem(domain, "P2", parameters);
        problem.addBoundaries(Teuchos::rcp(new BCBuilder<>()));
        problem.initializeProblem();
        problem.assemble(); // No updateConcentration(): mechanics-only zero field.
        auto displacement = problem.getSolution()->getBlockNonConst(0)->getDataNonConst(0);
        auto points = domain->getPointsUnique();
        for (LO i = 0; i < points->size(); ++i) {
            displacement[3*i] = 0.02 * (*points)[i][0] + 0.03 * (*points)[i][2];
            displacement[3*i+1] = 0.01 * (*points)[i][1];
        }
        problem.calculateNonLinResidualVec("standard");
        auto standard = Teuchos::rcp(new MultiVector<>(problem.getResidualVector()->getBlock(0)));
        problem.calculateNonLinResidualVec("reverse");
        standard->update(1., *problem.getResidualVector()->getBlock(0), 1.);
        for (const auto value : standard->getData(0))
            require(std::abs(value) < 1.e-10, "AceGen standard/reverse residual signs differ");

        const auto names = problem.getPostprocessingNames();
        auto first = problem.getPostProcessingData();
        if (structureModel != "SCI_SMC_Active_Growth_Reorientation") {
            require(names.empty() && first->size() == 0, "Unsupported fields must return an empty result");
        } else {
            require(names.size() == 17 && first->size() == 17, "Missing material output fields");
            require(names[11] == "Stretch1" && names[12] == "Stretch2", "Stretch field names");
            auto second = problem.getPostProcessingData();
            double stretchDifference = 0.;
            for (UN j = 0; j < first->size(); ++j) {
                auto a = first->getBlock(j)->getData(0);
                auto b = second->getBlock(j)->getData(0);
                for (LO i = 0; i < a.size(); ++i)
                    TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(a[i]) || !std::isfinite(b[i]) ||
                        std::abs(a[i] - b[i]) >= 1.e-10, std::runtime_error,
                        "Post-processing must be finite and repeatable: " << names[j]
                        << " node " << i << " values " << a[i] << ", " << b[i]);
            }
            auto stretch1 = first->getBlock(11)->getData(0);
            auto stretch2 = first->getBlock(12)->getData(0);
            auto coverage = Teuchos::rcp(new MultiVector<>(domain->getMapUnique(), 1));
            Teuchos::rcp_const_cast<FE<>>(problem.getFEFactory())->postProcessing(0, coverage);
            const auto covered = coverage->getData(0);
            for (LO i = 0; i < stretch1.size(); ++i) {
                if (covered[i] == 0.) {
                    require(stretch1[i] == 0. && stretch2[i] == 0., "Zero-weight nodes must project to zero");
                    continue;
                }
                TEUCHOS_TEST_FOR_EXCEPTION(stretch1[i] <= 0. || stretch2[i] <= 0., std::runtime_error,
                    "Stretches must be positive on covered nodes: node " << i
                    << " values " << stretch1[i] << ", " << stretch2[i]);
                stretchDifference += std::abs(stretch1[i] - stretch2[i]);
            }
            require(stretchDifference > 1.e-8, "Stretch2 must not duplicate Stretch1 under shear");
        }
        problem.assemble("UpdateTime");
        problem.assemble("UpdateTime");
        require(std::abs(problem.timeSteppingTool_->currentTime() - 0.03) < 1.e-14,
            "Solid clock must follow variable timesteps");
        require(problem.timeSteppingTool_->step_ == 2, "Solid clock must maintain step history");
        if (structureModel == "SCI_SMC_Active_Growth_Reorientation") {
            const std::string directory = "aceGenStructureUnitTest-process";
            std::filesystem::create_directories(directory);
            parameters->sublist("Timestepping Parameter").set("Checkpointing", true)
                .set("Checkpoint directory", directory);
            problem.writeCheckpointMetadata(.03);
            const auto saved = Teuchos::getParametersFromXmlFile(directory + "/" +
                checkpoint::manifestName(problem.checkpointMetadataForWrite(), .03));
            require(saved->sublist("Material process state").get<bool>("Active response took place"),
                "Checkpoint writer lost the material's active marker");
        }
        if (model == "SCI_NH") {
            const std::string directory = "aceGenStructureUnitTest-restart";
            std::filesystem::create_directories(directory);
            auto& settings = parameters->sublist("Timestepping Parameter");
            settings.set("Checkpointing", true).set("Checkpoint directory", directory)
                .set("Restart directory", directory);
            HDF5Export<> output(domain->getMapVecFieldUnique(), directory + "/Solutiond_s");
            output.writeVariablesHDF5(std::to_string(.03), problem.getSolution()->getBlock(0));
            output.closeExporter();
            problem.writeCheckpointMetadata(.03);
            auto resumedParameters = Teuchos::rcp(new Teuchos::ParameterList(*parameters));
            resumedParameters->sublist("Timestepping Parameter").set("Restart", true)
                .set("Time step", .03).set("Checkpointing", false);
            NonLinElasticity<> resumed(domain, "P2", resumedParameters);
            resumed.addBoundaries(Teuchos::rcp(new BCBuilder<>()));
            resumed.initializeProblem();
            resumed.assemble();
            require(resumed.timeSteppingTool_->step_ == 2, "Solid did not restore its checkpoint counter before assembly");
            resumed.assemble("UpdateTime");
            require(std::abs(resumed.timeSteppingTool_->currentTime() - .05) < 1.e-14 &&
                resumed.timeSteppingTool_->step_ == 3, "Solid continuation clock is wrong");
        }
    }
    // Exercise a real material continuation: the restarted adapter must use
    // time 0.1 and dt 0.001, and reproduce the passive stress and tangent.
    for (const bool computeFibers : {false, true}) {
        using Factory = AssembleFEFactory<default_sc,default_lo,default_go,default_no>;
        auto parameters = Teuchos::rcp(new Teuchos::ParameterList);
        parameters->sublist("Parameter").set("Dimension", 3);
        parameters->sublist("Timestepping Parameter").set("dt", .05);
        const auto& meshElement = domain->getElementsC()->getElement(0);
        const int flag = meshElement.getFlag();
        setPassiveMaterial(parameters, flag);
        parameters->sublist("Parameter Solid").sublist("1").set("Compute FA", computeFibers);
        vec2D_dbl_Type coordinates;
        for (const auto node : meshElement.getVectorNodeList())
            coordinates.push_back(domain->getPointsRepeated()->at(node));
        auto tuple = Teuchos::rcp(new tuple_disk_vec_Type{{"Solid", "P2", 3, 10}, {"Chemistry", "P2", 1, 10}});
        Factory factory;
        auto original = factory.build("SCI_SMC_Active_Growth_Reorientation", flag, coordinates, parameters, tuple);
        vec_dbl_Type solution(40, 0.);
        for (int i = 0; i < 10; ++i) {
            solution[3*i] = .02 * coordinates[i][0] + .03 * coordinates[i][2];
            solution[3*i+1] = .01 * coordinates[i][1];
        }
        original->updateSolution(solution);
        original->assembleJacobian();
        original->advanceInTime(.05);
        original->assembleJacobian();
        original->advanceInTime(.001);
        require(!original->getMaterialProcessState().active, "Passive material was marked active");
        auto restartedParameters = Teuchos::rcp(new Teuchos::ParameterList(*parameters));
        auto& stepping = restartedParameters->sublist("Timestepping Parameter");
        stepping.set("Restart", true).set("Time step", .1);
        stepping.sublist("_Restart time state").set("Physical time", .1).set("dt", .001);
        auto restarted = factory.build("SCI_SMC_Active_Growth_Reorientation", flag, coordinates, restartedParameters, tuple);
        require(std::abs(restarted->getTimeStep() - .1) < 1.e-14, "Material restart time is wrong");
        require(std::abs(restarted->getTimeIncrement() - .001) < 1.e-14, "Material restart increment is wrong");
        restarted->updateSolution(solution);
        original->assembleJacobian();
        restarted->assembleJacobian();
        original->assembleRHS();
        restarted->assembleRHS();
        for (int i = 0; i < 40; ++i) {
            const double a = original->getRHS()->at(i), b = restarted->getRHS()->at(i);
            require(std::isfinite(a) && std::isfinite(b) && std::abs(a-b) < 1.e-9*(1.+std::abs(a)),
                "Passive material residual changed across restart");
            for (int j = 0; j < 40; ++j) {
                const double x = (*original->getJacobian())[i][j], y = (*restarted->getJacobian())[i][j];
                require(std::isfinite(x) && std::isfinite(y) && std::abs(x-y) < 1.e-9*(1.+std::abs(x)),
                    "Passive material tangent changed across restart");
            }
        }
        // Once a process has changed history its marker must survive switches.
        auto processParameters = Teuchos::rcp(new Teuchos::ParameterList(*parameters));
        processParameters->sublist("Parameter Solid").sublist("1")
            .set("ActiveBool", 1.).set("growthBool", 1.).set("reorientationBool", 1.);
        auto marked = factory.build("SCI_SMC_Active_Growth_Reorientation", flag, coordinates, processParameters, tuple);
        marked->advanceInTime(.001);
        marked->advanceInTime(.001);
        const auto state = marked->getMaterialProcessState();
        require(state.active && state.growth && state.reorientation, "Material process markers were lost after switching off");
        auto boundaryParameters = Teuchos::rcp(new Teuchos::ParameterList(*parameters));
        auto& reorientation = boundaryParameters->sublist("Parameter Solid").sublist("1")
            .sublist("Timestepping Intervalls Reorientation");
        reorientation.set("Number of Segments", 1);
        reorientation.sublist("1").set("Start Time", .05).set("End Time", .1);
        auto boundary = factory.build("SCI_SMC_Active_Growth_Reorientation", flag, coordinates, boundaryParameters, tuple);
        boundary->advanceInTime(.001);
        require(boundary->getMaterialProcessState().reorientation, "Entering a process phase must exclude passive restart");
        boundary->advanceInTime(.001);
        require(boundary->getMaterialProcessState().reorientation, "Completed reorientation was not recorded");
    }
    std::cout << "AceGen structure checks passed\n";
}
