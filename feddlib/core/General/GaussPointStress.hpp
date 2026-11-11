#ifndef FEDD_GAUSS_POINT_STRESS_HPP
#define FEDD_GAUSS_POINT_STRESS_HPP

#include <array>
#include <vector>

namespace FEDD {

struct GaussPointStress {
    std::array<double, 3> reference{}, displacement{};
    std::array<double, 4> shapeP1{};
    std::array<double, 9> cauchy{}; // xx, xy, xz, yx, yy, yz, zx, zy, zz
    double weight = 0.; // Current-volume quadrature weight.
    double vonMises = 0.;
};

using ElementGaussPointStress = std::vector<std::vector<GaussPointStress>>;

}
#endif
