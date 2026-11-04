#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include <Tpetra_Core.hpp>
#include <cmath>
#include <iostream>

using namespace FEDD;
using LO = default_lo;

void require(bool condition, const char* message)
{
    TEUCHOS_TEST_FOR_EXCEPTION(!condition, std::runtime_error, message);
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
        parameters->sublist("Timestepping Parameter").set("dt", 0.01).set("Final time", 0.1);
        auto& intervals = parameters->sublist("Timestepping Parameter").sublist("Timestepping Intervalls");
        intervals.set("Number of Segments", 1);
        intervals.sublist("1").set("Start Time", 0.01).set("dt", 0.02);
        // Explicit parameters for the newer AceGen model; passive response only.
        parameters->sublist("Parameter Solid").set("Gamma1Min", 0.07797910389823641)
            .set("Gamma3Min", 0.13677878297366072).set("ActiveBool", 0.)
            .set("growthBool", 0.).set("reorientationBool", 0.).set("growthDirectionBool", 0.)
            .set("kFiberDistribution", 0.13333333333333333).set("P", 2.958083520579672)
            .set("PostProcessingPower", 2.);
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
    }
    std::cout << "AceGen structure checks passed\n";
}
