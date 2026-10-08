#ifndef FEDD_TEST_TRANSIENT_REFERENCE_HPP
#define FEDD_TEST_TRANSIENT_REFERENCE_HPP

#include "feddlib/core/Checkpointing/CheckpointFiles.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace TransientReference {

using Vector = FEDD::MultiVector<default_sc, default_lo, default_go, default_no>;

/**
 * @brief Compare a simulated field against a frozen unit-test reference.
 *
 * Normal runs read the existing "solution" dataset and never regenerate it.
 * Explicit reference-generation runs write to the selected directory instead.
 * Report absolute l2, relative l2 and absolute infinity errors. Require finite,
 * nonzero reference data and a relative l2 error at most 1e-12. The infinity
 * error must be at most 1e-11 + 1e-12 * ||reference||_infinity, so the absolute
 * bound scales with the field magnitude while retaining a floor for small fields.
 *
 * @param solution Runtime state at the documented comparison time.
 * @param file Stem of the reference file, without the .h5 suffix.
 * @param field Human-readable field name for the error report.
 * @param directory Directory containing the frozen references.
 * @param writeReference Whether explicit reference generation was requested.
 */
inline bool check(Teuchos::RCP<const Vector> solution, const std::string& file,
                  const std::string& field, const std::string& directory,
                  bool writeReference)
{
    const auto comm = solution->getMap()->getComm();
    const std::string path = FEDD::joinPath(directory, file);
    Teuchos::Array<default_sc> solutionNorm(1);
    solution->norm2(solutionNorm);
    if (!std::isfinite(solutionNorm[0]) || !(solutionNorm[0] > 0.)) {
        if (comm->getRank() == 0)
            std::cerr << "Invalid or zero computed field: " << field << std::endl;
        return false;
    }
    if (writeReference) {
        FEDD::HDF5Export<default_sc, default_lo, default_go, default_no> exporter(solution->getMap(), path);
        exporter.writeVariablesHDF5("solution", solution);
        exporter.closeExporter();
        if (comm->getRank() == 0)
            std::cout << "Wrote reference for " << field << ": " << path << ".h5" << std::endl;
        return true;
    }

    FEDD::HDF5Import<default_sc, default_lo, default_go, default_no> importer(solution->getMap(), path);
    auto reference = importer.readVariablesHDF5("solution");
    Vector error(solution->getMap());
    error.update(1., *solution, -1., *reference, 0.);
    Teuchos::Array<default_sc> errorL2(1), errorInf(1), referenceL2(1), referenceInf(1);
    error.norm2(errorL2);
    error.normInf(errorInf);
    reference->norm2(referenceL2);
    reference->normInf(referenceInf);
    constexpr double absoluteTolerance = 1.e-11;
    constexpr double relativeTolerance = 1.e-12;
    const double allowedInfError = absoluteTolerance + relativeTolerance * referenceInf[0];
    const double relative = referenceL2[0] > 0. ? errorL2[0] / referenceL2[0] : errorL2[0];
    const bool passed = std::isfinite(referenceL2[0]) && referenceL2[0] > 0.
        && std::isfinite(referenceInf[0]) && std::isfinite(allowedInfError)
        && std::isfinite(relative) && std::isfinite(errorInf[0])
        && relative <= relativeTolerance && errorInf[0] <= allowedInfError;
    if (comm->getRank() == 0)
        std::cout << std::setprecision(16)
                  << "Reference absolute error (" << field << ", l2): " << errorL2[0] << '\n'
                  << "Reference relative error (" << field << "): " << relative
                  << " (tolerance " << relativeTolerance << ")\n"
                  << "Reference norm (" << field << ", infinity): " << referenceInf[0] << '\n'
                  << "Reference absolute error (" << field << ", infinity): " << errorInf[0]
                  << " (allowed " << allowedInfError << " = " << absoluteTolerance
                  << " + " << relativeTolerance << " * reference infinity norm): "
                  << (passed ? "PASS" : "FAIL") << std::endl;
    return passed;
}

} // namespace TransientReference

#endif
