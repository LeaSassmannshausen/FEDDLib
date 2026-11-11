#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/StressPostProcessing.hpp"
#include <Tpetra_Core.hpp>
#include <cmath>
#include <iostream>
#include <aceinterface.hpp>
#include <filesystem>
#include <fstream>
#include <sstream>

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

void checkElementStressExport(NonLinElasticity<>& problem,
                             const Teuchos::RCP<Domain<>>& domain,
                             const ParameterListPtr_Type& parameters)
{
    auto fields = problem.getPostProcessingData(true);
    const auto gauss = problem.getGaussPointStressData();
    require(gauss.size() == domain->getElementsC()->numberElements(), "Missing Gauss-point elements");
    require(fields->size() == 4, "Element output must contain only the four stress fields");
    auto repeated = Teuchos::rcp(new MultiVector<>(domain->getMapVecFieldRepeated(), 1));
    repeated->importFromVector(problem.getSolution()->getBlock(0));
    const auto displacement = repeated->getData(0);
    const int columns[] = {10, 11, 12, 13};
    for (UN j = 0; j < fields->size(); ++j) {
        require(fields->getBlock(j)->getMap()->getTpetraMap()->isSameAs(*domain->getElementMap()->getTpetraMap()),
            "Element output has the wrong map");
        for (double value : fields->getBlock(j)->getData(0))
            require(std::isfinite(value), "Element stress must be finite");
    }

    // Independent reference: for power one, summing all ten nodal integrals
    // eliminates the shape functions by partition of unity. This must agree
    // with the element average even though the problem uses nodal power two.
    auto referenceParameters = Teuchos::rcp(new Teuchos::ParameterList(*parameters));
    referenceParameters->sublist("Parameter Solid").sublist("1").set("PostProcessingPower", 1.);
    auto discretization = Teuchos::rcp(new tuple_disk_vec_Type{
        tuple_ssii_Type("Solid", "P2", 3, 10), tuple_ssii_Type("Chemistry", "P2", 1, 10)});
    AssembleFEFactory<default_sc, default_lo, default_go, default_no> factory;
    for (UN T = 0; T < domain->getElementsC()->numberElements(); ++T) {
        auto element = domain->getElementsC()->getElement(T);
        auto nodes = element.getVectorNodeList();
        vec2D_dbl_Type positions;
        auto solution = Teuchos::rcp(new vec_dbl_Type(40, 0.));
        for (UN i = 0; i < nodes.size(); ++i) {
            positions.push_back(domain->getPointsRepeated()->at(nodes[i]));
            for (int k = 0; k < 3; ++k)
                solution->at(3*i+k) = displacement[3*nodes[i]+k];
        }
        auto material = factory.build("SCI_SMC_Active_Growth_Reorientation", element.getFlag(),
            positions, referenceParameters, discretization);
        material->setGlobalElementID(domain->getElementMap()->getGlobalElement(T));
        material->updateSolution(*solution);
        material->postProcessing();
        auto data = material->getPostProcessingData();
        double volume = 0.;
        for (const auto& row : *data) volume += row[0];
        double gaussVolume = 0.;
        for (const auto& q : gauss[T]) {
            require(q.weight > 0., "The four-point rule must have positive volume weights");
            gaussVolume += q.weight;
            const auto& x = q.reference;
            require(std::abs(q.displacement[0] - (.02*x[0]+.03*x[2]+.002*x[0]*x[2])) < 1.e-12,
                "Gauss-point coordinate/displacement interpolation is inconsistent");
        }
        require(std::abs(gaussVolume-std::abs(volume)) < 1.e-12, "Gauss-point volume is incorrect");
        for (int field = 0; field < 10; ++field) {
            double nativeIntegral = 0., gaussIntegral = 0.;
            const int column = field == 0 ? 10 : field;
            for (const auto& row : *data) nativeIntegral += row[column];
            for (const auto& q : gauss[T]) gaussIntegral += q.weight*gaussStressValue(q, field);
            const double expected = nativeIntegral/volume;
            require(std::abs(gaussIntegral/gaussVolume-expected) < 1.e-9*(1.+std::abs(expected)),
                "Raw Gauss-point stresses disagree with native task-5 output");
        }
        for (UN j = 0; j < fields->size(); ++j) {
            double integral = 0.;
            for (const auto& row : *data) integral += row[columns[j]];
            const double expected = integral / volume;
            require(std::abs(fields->getBlock(j)->getData(0)[T] - expected) < 1.e-9 * (1. + std::abs(expected)),
                "Element stress differs from the partition-of-unity volume average");
        }
    }

    const std::string filename = "aceGenStructureUnitTest-element-stresses";
    ExporterParaView<> exporter;
    exporter.setup(filename, domain->getMesh(), "P0");
    const auto names = problem.getPostprocessingNames();
    for (UN j = 0; j < fields->size(); ++j) {
        Teuchos::RCP<const MultiVector<>> field = fields->getBlock(j);
        exporter.addVariable(field, names[j], "Scalar", 1, domain->getElementMap());
    }
    exporter.save(.01);
    exporter.closeExporter();

    // A consistent projection must reproduce a linear field exactly, including
    // values at P2 midsides, even with nonuniform current-volume weights.
    auto linear = gauss;
    for (auto& element : linear)
        for (auto& q : element) {
            const double value = 2.+.3*q.reference[0]-.2*q.reference[1]+.1*q.reference[2];
            q.vonMises = value;
            for (int field = 0; field < 9; ++field) q.cauchy[field] = (field+2)*value;
        }
    auto recovered = projectGaussPointStresses<default_sc,default_lo,default_go,default_no>(domain->getMesh(), linear);
    for (int field = 0; field < 10; ++field)
        for (LO i = 0; i < domain->getMapUnique()->getNodeNumElements(); ++i) {
            const auto& x = domain->getPointsUnique()->at(i);
            const double expected = (field+1)*(2.+.3*x[0]-.2*x[1]+.1*x[2]);
            require(std::abs(recovered->getBlock(field)->getData(0)[i]-expected) < 1.e-9,
                "P1 stress recovery does not reproduce a linear field");
        }
    auto exportParameters = Teuchos::rcp(new Teuchos::ParameterList(*parameters));
    exportParameters->sublist("Exporter").set("Export Gauss point stresses", true)
        .set("Export P1 stresses", true).set("Export every X timesteps", 1);
    StressPostProcessing<> stressExporter;
    stressExporter.save(domain->getMesh(), gauss, exportParameters, .01);
    stressExporter.save(domain->getMesh(), gauss, exportParameters, .02);
    stressExporter.close();
    if (domain->getComm()->getRank() == 0) {
        std::ifstream file(filename + ".xmf");
        std::ostringstream contents;
        contents << file.rdbuf();
        const auto xml = contents.str();
        require(xml.find("Center=\"Cell\"") != std::string::npos, "Stress output must be cell data");
        require(xml.find("vonMisesStress") != std::string::npos, "Element von Mises output is missing");
        require(xml.find("Stretch1") == std::string::npos, "Element exporter must contain stresses only");
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
            displacement[3*i] = 0.02 * (*points)[i][0] + 0.03 * (*points)[i][2]
                + 0.002 * (*points)[i][0] * (*points)[i][2];
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
            checkElementStressExport(problem, domain, parameters);
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
            // Repeat the integral check after material history has advanced
            // into the native active-response time interval.
            const auto gauss = problem.getGaussPointStressData();
            const auto averaged = problem.getPostProcessingData(true);
            for (UN T = 0; T < gauss.size(); ++T) {
                double volume = 0., mises = 0.;
                for (const auto& q : gauss[T]) { volume += q.weight; mises += q.weight*q.vonMises; }
                const double expected = averaged->getBlock(0)->getData(0)[T];
                require(std::abs(mises/volume-expected) < 1.e-9*(1.+std::abs(expected)),
                    "Gauss stresses disagree after advancing material history");
            }
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
