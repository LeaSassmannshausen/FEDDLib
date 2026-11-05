#include "feddlib/core/General/BCBuilder.hpp"
#include <Tpetra_Core.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace FEDD;

namespace {
const double cosine = std::cos(0.6), sine = std::sin(0.6);
double inletHeight = 1.;
bool zeroProfile = false;

void profile(double* x, double* result, double, const double* parameters)
{
    const double y = (-sine*x[0] + cosine*x[1]) / inletHeight;
    const double speed = zeroProfile ? 0. : parameters[0] * 4.*y*(1.-y);
    result[0] = speed * cosine;
    result[1] = speed * sine;
}

void flowRate(double*, double* result, double t, const double* parameters)
{
    result[0] = parameters[1] * t;
}

void importedProfile(double* x, double* result, double, const double* parameters)
{
    result[0] = parameters[0] * x[0];
    result[1] = 0.;
}

void require(bool condition, const char* message)
{
    TEUCHOS_TEST_FOR_EXCEPTION(!condition, std::runtime_error, message);
}
}

int main(int argc, char** argv)
{
    Tpetra::ScopeGuard scope(&argc, &argv);
    auto comm = Tpetra::getDefaultComm();
    const int subdivisions = static_cast<int>(std::sqrt(comm->getSize()));
    require(subdivisions*subdivisions == comm->getSize(), "Use a square number of MPI ranks");
    for (const std::string feType : {"P1", "P2"}) {
        inletHeight = 1.;
        auto domain = Teuchos::rcp(new Domain<>(std::vector<double>{0., 0.}, 1., 1., comm));
        domain->buildMesh(1, "Square", 2, feType, subdivisions, 4);
        // The structured 2D generator does not create boundary subelements.
        auto elements = domain->getElementsC();
        const auto pointsRepeated = domain->getPointsRepeated();
        for (default_lo element = 0; element < elements->numberElements(); ++element) {
            vec_LO_Type nodes;
            for (auto node : elements->getElement(element).getVectorNodeList())
                if (std::abs((*pointsRepeated)[node][0]) < 1.e-12)
                    nodes.push_back(node);
            if (nodes.size() != (feType == "P1" ? 2 : 3)) continue;
            std::sort(nodes.begin(), nodes.end(), [&](default_lo a, default_lo b) {
                return (*pointsRepeated)[a][1] < (*pointsRepeated)[b][1];
            });
            if (feType == "P2") std::swap(nodes[1], nodes[2]);
            auto& triangle = elements->getElement(element);
            if (!triangle.subElementsInitialized()) triangle.initializeSubElements(feType, 1);
            FiniteElement edge(nodes, 2);
            triangle.addSubElement(edge);
        }
        // Exercise the original imported scalar-profile overload as well.
        {
            auto imported = Teuchos::rcp(new MultiVector<>(domain->getMapUnique()));
            const auto points = domain->getPointsUnique();
            auto data = imported->getDataNonConst(0);
            for (default_lo node = 0; node < points->size(); ++node) {
                const double y = (*points)[node][1];
                data[node] = 4.*y*(1.-y);
            }
            Teuchos::RCP<const MultiVector<>> external = imported;
            vec_dbl_Type parameters{1., 3.};
            BCBuilder<> legacy;
            legacy.addBC(importedProfile, 2, 0, domain, "Dirichlet", 2, parameters, external, true, flowRate);
            auto velocity = Teuchos::rcp(new MultiVector<>(domain->getMapVecFieldUnique()));
            velocity->putScalar(0.);
            auto rhs = Teuchos::rcp(new BlockMultiVector<>(1));
            rhs->addBlock(velocity, 0);
            legacy.setRHS(rhs, 0.25);
            auto repeated = Teuchos::rcp(new MultiVector<>(domain->getMapVecFieldRepeated()));
            repeated->importFromVector(velocity, false, "Insert");
            FE<> fe;
            fe.addFE(domain);
            double flux = 0.;
            fe.assemblyFlowRate(2, flux, feType, 2, 2, repeated);
            require(std::abs(flux-0.75) < 1.e-10, "Imported scalar-profile flow rate changed");
        }
        for (auto points : {domain->getPointsUnique(), domain->getPointsRepeated()})
            for (auto& point : *points) {
                const double x = point[0], y = point[1];
                point[0] = cosine*x - sine*y;
                point[1] = sine*x + cosine*y;
            }
        BCBuilder<> boundaries;
        vec_dbl_Type parameters{1., 3.};
        Teuchos::RCP<const MultiVector<>> external;
        boundaries.addBC(profile, 2, 0, domain, "Dirichlet", 2, parameters, external, true, flowRate);
        auto values = Teuchos::rcp(new MultiVector<>(domain->getMapVecFieldUnique()));
        auto rhs = Teuchos::rcp(new BlockMultiVector<>(1));
        rhs->addBlock(values, 0);
        auto difference = Teuchos::rcp(new BlockMultiVector<>(1));
        auto offset = Teuchos::rcp(new MultiVector<>(domain->getMapVecFieldUnique()));
        offset->putScalar(2.);
        difference->addBlock(offset, 0);
        auto repeated = Teuchos::rcp(new MultiVector<>(domain->getMapVecFieldRepeated()));
        FE<> fe;
        fe.addFE(domain);
        for (double height : {1., 1.5}) {
            const double ratio = height/inletHeight;
            for (auto points : {domain->getPointsUnique(), domain->getPointsRepeated()})
                for (auto& point : *points) {
                    const double x = cosine*point[0] + sine*point[1];
                    const double y = (-sine*point[0] + cosine*point[1]) * ratio;
                    point[0] = cosine*x - sine*y;
                    point[1] = sine*x + cosine*y;
                }
            inletHeight = height;
            for (double time : {0., 0.25, 1.}) {
                for (int mode = 0; mode < 3; ++mode) {
                    values->putScalar(0.);
                    if (mode == 0) boundaries.setRHS(rhs, time);
                    if (mode == 1) boundaries.setBCMinusVector(rhs, difference, time);
                    if (mode == 2) boundaries.setVectorMinusBC(rhs, difference, time);
                    // Recover the prescribed velocity from the two residual forms.
                    auto data = values->getDataNonConst(0);
                    const auto flags = domain->getBCFlagUnique();
                    const auto points = domain->getPointsUnique();
                    const double intervals = 4.*subdivisions;
                    const double integral = height * 2./3. * (feType == "P1" ? 1.-1./(intervals*intervals) : 1.);
                    vec_dbl_Type expected(2), unitParameters{3.*time/integral};
                    for (default_lo node = 0; node < points->size(); ++node) {
                        if ((*flags)[node] != 2) continue;
                        if (mode == 1) { data[2*node] += 2.; data[2*node+1] += 2.; }
                        if (mode == 2) { data[2*node] = 2.-data[2*node]; data[2*node+1] = 2.-data[2*node+1]; }
                        auto point = (*points)[node];
                        profile(point.data(), expected.data(), time, unitParameters.data());
                        require(std::abs(data[2*node]-expected[0]) < 1.e-10 &&
                            std::abs(data[2*node+1]-expected[1]) < 1.e-10,
                            "Flow-normalized boundary values differ from the exact discrete profile integral");
                    }
                    repeated->importFromVector(values, false, "Insert");
                    double flux = 0.;
                    fe.assemblyFlowRate(2, flux, feType, 2, 2, repeated);
                    require(std::abs(flux-3.*time) < 1.e-10, "Prescribed flux differs from target");
                }
            }
        }
        zeroProfile = true;
        bool rejected = false;
        try { boundaries.setRHS(rhs, 1.); }
        catch (const std::logic_error&) { rejected = true; }
        zeroProfile = false;
        require(rejected, "A zero-flux profile must be rejected");
    }
    if (comm->getRank() == 0)
        std::cout << "Flow-rate BC checks passed: P1/P2, rotated inlet, ramp, moving mesh and residual forms\n";
}
