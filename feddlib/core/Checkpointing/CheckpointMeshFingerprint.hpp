#ifndef FEDD_CHECKPOINT_MESH_FINGERPRINT_HPP
#define FEDD_CHECKPOINT_MESH_FINGERPRINT_HPP

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FiniteElement.hpp"
#include <Teuchos_CommHelpers.hpp>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace FEDD {
namespace checkpoint {

/// Hash one mesh record; aggregate record hashes independently of traversal order.
inline unsigned long long hashRecord(const std::string& record, unsigned long long seed)
{
    for (unsigned char byte : record) {
        seed ^= byte;
        seed *= 1099511628211ULL;
    }
    return seed;
}

/** @brief Compute a partition-independent reference mesh and field DOF identity.
 * Sums two independent 64-bit record hashes over uniquely owned nodes, DOFs and
 * elements, including boundary markers and subelements. MPI ownership may change;
 * global numbering, coordinates and connectivity must match. These checksums
 * describe compatibility, not cryptographic integrity.
 *
 * @param[in] domain Field domain, including its reference mesh and boundary data.
 * @param[in] map Unique scalar or vector field map.
 * @param[in] comm Communicator owning the domain and map.
 * @return Two hexadecimal checksums separated by a colon.
 * @pre Call collectively before ALE moves the reference coordinates.
 */
template<class SC, class LO, class GO, class NO>
std::string meshAndDofFingerprint(const Domain<SC, LO, GO, NO>& domain,
                                  const Map<LO, GO, NO>& map,
                                  const Teuchos::Comm<int>& comm)
{
    unsigned long long local[2] = {0, 0}, global[2] = {0, 0};
    const auto record = [&](const std::string& value) {
        local[0] += hashRecord(value, 14695981039346656037ULL);
        local[1] += hashRecord(value, 7809847782465536322ULL);
    };
    const auto nodes = domain.getMapUnique();
    const auto points = domain.getPointsUnique();
    const auto flags = domain.getBCFlagUnique();
    for (UN j = 0; j < nodes->getNodeNumElements(); ++j) {
        std::ostringstream value;
        value << "node " << nodes->getGlobalElement(j) << std::hexfloat;
        for (double coordinate : points->at(j)) value << ' ' << coordinate;
        if (!flags.is_null()) value << " flag " << flags->at(j);
        record(value.str());
    }
    for (UN j = 0; j < map.getNodeNumElements(); ++j)
        record("dof " + std::to_string(map.getGlobalElement(j)));
    const auto elements = domain.getElementsC();
    // FSI's dummy interface domain has nodes/DOFs but no volume element map.
    if (elements->numberElements() > 0) {
        const auto elementMap = domain.getElementMap();
        const auto repeated = domain.getMapRepeated();
        std::function<std::string(FiniteElement)> elementRecord = [&](FiniteElement element) {
            std::ostringstream value;
            value << "flag " << element.getFlag();
            for (auto node : element.getVectorNodeList()) value << ' ' << repeated->getGlobalElement(node);
            std::vector<std::string> surfaces;
            auto children = element.getSubElements();
            if (!children.is_null())
                for (UN k = 0; k < children->numberElements(); ++k)
                    surfaces.push_back(elementRecord(children->getElement(k)));
            std::sort(surfaces.begin(), surfaces.end());
            for (const auto& surface : surfaces) value << " [" << surface << ']';
            return value.str();
        };
        for (UN j = 0; j < elementMap->getNodeNumElements(); ++j)
            record("element " + std::to_string(elementMap->getGlobalElement(j)) + " " + elementRecord(elements->getElement(j)));
    }
    Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 2, local, global);
    std::ostringstream fingerprint;
    fingerprint << std::hex << global[0] << ':' << global[1];
    return fingerprint.str();
}

} // namespace checkpoint
} // namespace FEDD

#endif
