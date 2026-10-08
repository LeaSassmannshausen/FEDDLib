#ifndef FEDD_HDF5_VECTOR_INFO_HPP
#define FEDD_HDF5_VECTOR_INFO_HPP

#include <hdf5.h>
#include <stdexcept>
#include <string>

namespace FEDD {
namespace checkpoint {

/// Own an HDF5 handle, including during compatibility-check exceptions.
class H5Handle {
    hid_t id_;
    herr_t (*close_)(hid_t);
public:
    H5Handle(hid_t id, herr_t (*close)(hid_t)) : id_(id), close_(close) {}
    ~H5Handle() { if (id_ >= 0) close_(id_); }
    H5Handle(const H5Handle&) = delete;
    H5Handle& operator=(const H5Handle&) = delete;
    operator hid_t() const { return id_; }
};

/** @brief Inspect a vector without reading or modifying its values.
 * Supports the row-vector layout written by FEDDLib and legacy Epetra files.
 * Throws before any hyperslab read if properties, dimensions or type disagree.
 */
inline void validateVectorDataset(hid_t file, const std::string& key,
                                  unsigned long long length, int vectors = 1)
{
    const auto fail = [&](const std::string& reason) {
        throw std::runtime_error("Checkpoint dataset '" + key + "': " + reason);
    };
    if (H5Lexists(file, key.c_str(), H5P_DEFAULT) <= 0)
        fail("missing required history/field group");
    H5Handle group(H5Gopen(file, key.c_str(), H5P_DEFAULT), H5Gclose);
    if (group < 0) fail("cannot open group");
    const auto scalar = [&](const char* name) {
        if (H5Lexists(group, name, H5P_DEFAULT) <= 0)
            fail(std::string("missing ") + name);
        H5Handle dataset(H5Dopen(group, name, H5P_DEFAULT), H5Dclose);
        if (dataset < 0) fail(std::string("cannot open ") + name);
        H5Handle space(H5Dget_space(dataset), H5Sclose);
        H5Handle type(H5Dget_type(dataset), H5Tclose);
        if (space < 0 || type < 0 || H5Sget_simple_extent_npoints(space) != 1 ||
            H5Tget_class(type) != H5T_INTEGER)
            fail(std::string("invalid scalar property ") + name);
        long long value = 0;
        if (H5Dread(dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) < 0)
            fail(std::string("cannot read ") + name);
        return value;
    };
    const auto savedLength = scalar("GlobalLength");
    if (savedLength < 0 || static_cast<unsigned long long>(savedLength) != length)
        fail("GlobalLength mismatch: saved " + std::to_string(savedLength) +
             ", expected " + std::to_string(length));
    if (scalar("NumVectors") != vectors)
        fail("NumVectors mismatch");
    if (H5Lexists(group, "Values", H5P_DEFAULT) <= 0) fail("missing Values dataset");
    H5Handle dataset(H5Dopen(group, "Values", H5P_DEFAULT), H5Dclose);
    if (dataset < 0) fail("cannot open Values");
    H5Handle space(H5Dget_space(dataset), H5Sclose);
    H5Handle type(H5Dget_type(dataset), H5Tclose);
    hsize_t dimensions[2] = {0, 0};
    if (space < 0 || H5Sget_simple_extent_ndims(space) != 2 ||
        H5Sget_simple_extent_dims(space, dimensions, nullptr) < 0 ||
        dimensions[0] != static_cast<hsize_t>(vectors) || dimensions[1] != length)
        fail("Values dimensions disagree with GlobalLength/NumVectors");
    if (type < 0 || H5Tget_class(type) != H5T_FLOAT || H5Tget_size(type) != sizeof(double))
        fail("Values must contain 64-bit floating-point data");
}

} // namespace checkpoint
} // namespace FEDD
#endif
