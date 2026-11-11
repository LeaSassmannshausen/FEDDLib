#ifndef FEDD_STRESS_POST_PROCESSING_HPP
#define FEDD_STRESS_POST_PROCESSING_HPP

#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/General/StressRecovery.hpp"

namespace FEDD {

// Keeps the diagnostic point cloud and the recovered continuous field in
// separate ParaView series, using the existing cadence and restart handling.
template<class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class StressPostProcessing {
    using Vector = MultiVector<SC,LO,GO,NO>;
    using BlockVector = BlockMultiVector<SC,LO,GO,NO>;
    using Exporter = ExporterParaView<SC,LO,GO,NO>;
    using MeshType = Mesh<SC,LO,GO,NO>;
    Teuchos::RCP<Exporter> gaussExporter_, p1Exporter_;
    Teuchos::RCP<BlockVector> gaussValues_, p1Values_;
    Teuchos::RCP<Vector> displacement_;
    Teuchos::RCP<const Map<LO,GO,NO>> gaussMap_;
    vec2D_dbl_ptr_Type gaussCoordinates_;
    const vec_string_Type names_{"vonMisesStress", "Sxx", "Sxy", "Sxz", "Syx", "Syy", "Syz", "Szx", "Szy", "Szz"};
public:
    ~StressPostProcessing() { close(); }

    void close() {
        if (!gaussExporter_.is_null()) { gaussExporter_->closeExporter(); gaussExporter_ = Teuchos::null; }
        if (!p1Exporter_.is_null()) { p1Exporter_->closeExporter(); p1Exporter_ = Teuchos::null; }
    }

    void save(Teuchos::RCP<MeshType> mesh, const ElementGaussPointStress& data,
              ParameterListPtr_Type parameters, double time)
    {
        auto& options = parameters->sublist("Exporter");
        const int interval = options.get("Export every X timesteps", 1);
        TEUCHOS_TEST_FOR_EXCEPTION(data.size() != mesh->getElementsC()->numberElements(),
            std::logic_error, "Gauss-point stress output does not match the solid mesh.");
        if (options.get("Export Gauss point stresses", false)) {
            int localPoints = data.empty() ? 0 : data.front().size(), pointsPerElement = 0;
            Teuchos::reduceAll(*mesh->getComm(), Teuchos::REDUCE_MAX, 1, &localPoints, &pointsPerElement);
            TEUCHOS_TEST_FOR_EXCEPTION(pointsPerElement <= 0, std::logic_error, "No Gauss-point stresses to export.");
            for (const auto& element : data)
                TEUCHOS_TEST_FOR_EXCEPTION(element.size() != pointsPerElement, std::logic_error,
                    "Point-cloud output requires the same quadrature rule in every solid element.");
            if (gaussExporter_.is_null()) {
                gaussMap_ = mesh->getElementMap()->buildVecFieldMap(pointsPerElement);
                gaussCoordinates_ = Teuchos::rcp(new vec2D_dbl_Type(data.size()*pointsPerElement, vec_dbl_Type(3)));
                gaussValues_ = Teuchos::rcp(new BlockVector(13));
                for (int field = 0; field < 13; ++field)
                    gaussValues_->addBlock(Teuchos::rcp(new Vector(gaussMap_, 1)), field);
                displacement_ = Teuchos::rcp(new Vector(gaussMap_->buildVecFieldMap(3), 1));
                for (UN T = 0; T < data.size(); ++T)
                    for (int g = 0; g < pointsPerElement; ++g)
                        for (int d = 0; d < 3; ++d)
                            gaussCoordinates_->at(T*pointsPerElement+g)[d] = data[T][g].reference[d];
                gaussExporter_ = Teuchos::rcp(new Exporter);
                gaussExporter_->setupPointCloud("solid_stress_gauss", gaussMap_, gaussCoordinates_, parameters);
                auto names = names_;
                names.insert(names.end(), {"elementID", "gaussPointID", "currentVolumeWeight"});
                for (int field = 0; field < 13; ++field) {
                    Teuchos::RCP<const Vector> value = gaussValues_->getBlock(field);
                    gaussExporter_->addVariable(value, names[field], "Scalar", 1, gaussMap_);
                }
                Teuchos::RCP<const Vector> displacement = displacement_;
                gaussExporter_->addVariable(displacement, "displacement", "Vector", 3, gaussMap_);
            }
            TEUCHOS_TEST_FOR_EXCEPTION(gaussMap_->getNodeNumElements() != data.size()*pointsPerElement,
                std::logic_error, "Gauss-point output changed its point count during the run.");
            for (UN T = 0; T < data.size(); ++T)
                for (int g = 0; g < pointsPerElement; ++g) {
                    const UN index = T*pointsPerElement+g;
                    const auto& q = data[T][g];
                    for (int field = 0; field < 10; ++field)
                        gaussValues_->getBlockNonConst(field)->getDataNonConst(0)[index] = gaussStressValue(q, field);
                    gaussValues_->getBlockNonConst(10)->getDataNonConst(0)[index] = mesh->getElementMap()->getGlobalElement(T);
                    gaussValues_->getBlockNonConst(11)->getDataNonConst(0)[index] = g;
                    gaussValues_->getBlockNonConst(12)->getDataNonConst(0)[index] = q.weight;
                    for (int d = 0; d < 3; ++d) displacement_->getDataNonConst(0)[3*index+d] = q.displacement[d];
                }
            gaussExporter_->save(time);
        }
        if (options.get("Export P1 stresses", false)) {
            auto projected = projectGaussPointStresses<SC,LO,GO,NO>(mesh, data);
            if (p1Exporter_.is_null()) {
                p1Values_ = projected;
                p1Exporter_ = Teuchos::rcp(new Exporter);
                // A P1 field represented on the original P2 geometry: each
                // midside value is the mean of its two corner values.
                p1Exporter_->setup("solid_stress_p1", mesh, "P2", interval, parameters);
                for (int field = 0; field < 10; ++field) {
                    Teuchos::RCP<const Vector> value = p1Values_->getBlock(field);
                    p1Exporter_->addVariable(value, names_[field], "Scalar", 1, mesh->getMapUnique());
                }
            } else p1Values_->update(1., *projected, 0.);
            p1Exporter_->save(time);
        }
    }
};

}
#endif
