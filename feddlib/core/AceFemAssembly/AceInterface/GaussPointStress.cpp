#include "GaussPointStress.hpp"
#include "feddlib/core/core_config.h"

#ifdef FEDD_HAVE_ACEGENINTERFACE
#include <ace2.h>
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_TestForException.hpp>
#include <cmath>
#include <cstring>
#include <string>
#include <unistd.h>

namespace FEDD {
namespace acegen {

// AceGen expects eigenvectors in rows. LAPACK's column-major result has that
// layout when read as the row-major array passed by the generated C routine.
inline void stressEigenSystem(int* n, double* matrix, double* values, int* info)
{
    Teuchos::LAPACK<int, double> lapack;
    std::vector<double> work(8 * (*n));
    lapack.SYEV('V', 'U', *n, matrix, *n, values, work.data(), work.size(), info);
}

inline double stressDeterminant(const double a[3][3])
{
    return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
         - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
         + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
}

// Wrap the installed library's public C entry point; no change to Interface2
// is required. All native state, including the updated history, is a local copy.
std::vector<GaussPointStress> gaussPointStress(
    const std::vector<double>& positions, const double* displacement,
    const double* concentration, const double* rates, const double* acceleration,
    std::vector<double> domainData, std::vector<double> history,
    double tolerance, double dt, double time, int integrationCode, int elementID)
{
    ElementSpec es{};
    std::vector<int> idata(IData_Length, 0);
    idata[ID_MaxMessages] = 1000000;
    idata[ID_CurrentElement] = idata[ID_DebugElement] = elementID;
    TEUCHOS_TEST_FOR_EXCEPTION(SMTSetElSpecSMCActiveGrowthReorientation(
        &es, idata.data(), integrationCode, -1) != 0, std::runtime_error,
        "Cannot initialize AceGen Gauss-point stress output.");
    es.IntPoints = SMTIntPoints(&es.id.IntCode, &es.id.NoIntPoints);
    TEUCHOS_TEST_FOR_EXCEPTION(es.id.NoNodes != 40 || positions.size() != 30 ||
        es.id.NoTimeStorage != history.size() || es.id.NoIntPoints <= 0,
        std::logic_error, "Unsupported AceGen Gauss-point stress layout.");
    SMTStructure smt{};
    smt.NumericalLibrary.SolveSymmetricEigensystem = &stressEigenSystem;
    es.SMT = &smt;
    std::string log = ".ErrorLog.raw." + std::to_string(getpid()) + ".txt";
    es.OutFileName = log.data();
    es.Data = domainData.data();
    std::vector<NodeData> nodes(es.id.NoNodes);
    std::vector<NodeSpec> specs(es.id.NoNodes);
    std::vector<NodeData*> nd(es.id.NoNodes);
    std::vector<NodeSpec*> ns(es.id.NoNodes);
    std::vector<std::array<double, 3>> x(es.id.NoNodes), at(es.id.NoNodes), ap(es.id.NoNodes);
    int c = 0, rate = 0, acc = 0;
    for (int i = 0; i < es.id.NoNodes; ++i) {
        nd[i] = &nodes[i]; ns[i] = &specs[i];
        nodes[i].X = x[i].data(); nodes[i].at = at[i].data(); nodes[i].ap = ap[i].data();
        if (std::strcmp(es.NodeID[i], "D") == 0) {
            for (int d = 0; d < 3; ++d) {
                x[i][d] = positions[3*i+d]; at[i][d] = displacement[3*i+d];
            }
        } else if (std::strcmp(es.NodeID[i], "C") == 0) at[i][0] = concentration[c++];
        else if (std::strcmp(es.NodeID[i], "CD") == 0) at[i][0] = rates[rate++];
        else if (std::strcmp(es.NodeID[i], "DDD") == 0) {
            for (int d = 0; d < 3; ++d) at[i][d] = acceleration[3*acc+d];
            ++acc;
        }
    }
    auto previous = history;
    ElementData ed{};
    ed.hp = history.data(); ed.ht = previous.data();
    std::vector<double> rdata(RData_Length, 0.), work(es.id.WorkingVectorSize, 0.);
    rdata[RD_SubIterationTolerance] = tolerance;
    rdata[RD_TimeIncrement] = dt; rdata[RD_Time] = time;
    std::vector<GaussPointStress> result(es.id.NoIntPoints);
    std::vector<double*> stresses(es.id.NoIntPoints);
    for (int g = 0; g < es.id.NoIntPoints; ++g) stresses[g] = result[g].cauchy.data();
    // Task 6's coordinate indices overlap; give it a safe scratch buffer and
    // compute the actual coordinates below using the native P2 node ordering.
    std::vector<double> coordinateScratch(3*es.id.NoIntPoints+3, 0.);
    int task = 6;
    TasksSMCActiveGrowthReorientation(work.data(), &es, &ed, ns.data(), nd.data(),
        rdata.data(), idata.data(), coordinateScratch.data(), stresses.data(),
        &task, nullptr, nullptr, nullptr, nullptr, nullptr);
    TEUCHOS_TEST_FOR_EXCEPTION(idata[ID_ErrorStatus] >= 2, std::runtime_error,
        "AceGen Gauss-point stress output failed in element " << elementID);
    const int edges[6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
    const double dl[4][3] = {{1,0,0},{-1,-1,-1},{0,0,1},{0,1,0}};
    for (int g = 0; g < es.id.NoIntPoints; ++g) {
        auto& q = result[g];
        const double* gp = es.IntPoints + 4*g;
        q.shapeP1 = {{gp[0], 1.-gp[0]-gp[1]-gp[2], gp[2], gp[1]}};
        double phi[10]{}, derivative[10][3]{}, referenceJacobian[3][3]{}, currentJacobian[3][3]{};
        for (int i = 0; i < 4; ++i) {
            phi[i] = q.shapeP1[i]*(2.*q.shapeP1[i]-1.);
            for (int d = 0; d < 3; ++d) derivative[i][d] = (4.*q.shapeP1[i]-1.)*dl[i][d];
        }
        for (int e = 0; e < 6; ++e) {
            const int i = edges[e][0], j = edges[e][1];
            phi[4+e] = 4.*q.shapeP1[i]*q.shapeP1[j];
            for (int d = 0; d < 3; ++d)
                derivative[4+e][d] = 4.*(dl[i][d]*q.shapeP1[j]+q.shapeP1[i]*dl[j][d]);
        }
        for (int i = 0; i < 10; ++i)
            for (int d = 0; d < 3; ++d) {
                q.reference[d] += phi[i]*positions[3*i+d];
                q.displacement[d] += phi[i]*displacement[3*i+d];
                for (int k = 0; k < 3; ++k) {
                    referenceJacobian[d][k] += derivative[i][k]*positions[3*i+d];
                    currentJacobian[d][k] += derivative[i][k]*(positions[3*i+d]+displacement[3*i+d]);
                }
            }
        const double referenceDet = stressDeterminant(referenceJacobian);
        const double currentDet = stressDeterminant(currentJacobian);
        TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(currentDet) || !std::isfinite(referenceDet) ||
            referenceDet == 0. || currentDet/referenceDet <= 0., std::runtime_error,
            "Invalid element geometry in Gauss-point stress output: " << elementID);
        q.weight = std::abs(currentDet)*gp[3];
        double trace = (q.cauchy[0]+q.cauchy[4]+q.cauchy[8])/3., deviatorNorm = 0.;
        for (int i = 0; i < 9; ++i) {
            TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(q.cauchy[i]), std::runtime_error,
                "Non-finite Gauss-point stress in element " << elementID);
            double value = q.cauchy[i] - (i%4 == 0 ? trace : 0.);
            deviatorNorm += value*value;
        }
        q.vonMises = std::sqrt(1.5*deviatorNorm);
    }
    return result;
}

}
}

#endif
