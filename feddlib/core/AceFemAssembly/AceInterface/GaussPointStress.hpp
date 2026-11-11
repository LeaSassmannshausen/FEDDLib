#ifndef FEDD_ACEGEN_GAUSS_POINT_STRESS_HPP
#define FEDD_ACEGEN_GAUSS_POINT_STRESS_HPP

#include "feddlib/core/General/GaussPointStress.hpp"

namespace FEDD {
namespace acegen {

// Direct native task-6 output. The implementation owns all native scratch
// state and computes coordinates independently of AceGen's coordinate buffer.
std::vector<GaussPointStress> gaussPointStress(
    const std::vector<double>& positions, const double* displacement,
    const double* concentration, const double* rates, const double* acceleration,
    std::vector<double> domainData, std::vector<double> history,
    double tolerance, double dt, double time, int integrationCode, int elementID);

}
}
#endif
