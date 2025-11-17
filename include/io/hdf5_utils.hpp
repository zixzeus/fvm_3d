#ifndef HDF5_UTILS_HPP
#define HDF5_UTILS_HPP

#include "core/field3d.hpp"
#include <hdf5.h>
#include <string>
#include <vector>
#include <stdexcept>

namespace fvm3d::io::hdf5_utils {

/**
 * @brief Get variable names based on number of variables
 *
 * Returns a const reference to static arrays, avoiding allocations.
 *
 * @param num_vars Number of state variables (5 for Euler, 8/9 for MHD)
 * @return Const reference to vector of variable names
 */
const std::vector<std::string>& get_variable_names(int num_vars);

/**
 * @brief Extract interior cells from state field to flat array
 *
 * Extracts data from a single variable, excluding ghost cells.
 * Inlined for performance in hot paths.
 *
 * @param state State field (with ghost cells)
 * @param var_index Variable index to extract
 * @param nx, ny, nz Interior cell dimensions (without ghosts)
 * @param nghost Number of ghost cells
 * @param output Output vector (will be resized to nx*ny*nz)
 */
inline void extract_interior_cells(
    const core::StateField3D& state,
    int var_index,
    int nx, int ny, int nz,
    int nghost,
    std::vector<double>& output
) {
    output.resize(nx * ny * nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx = i * ny * nz + j * nz + k;
                output[idx] = state(var_index, i + nghost, j + nghost, k + nghost);
            }
        }
    }
}

/**
 * @brief Insert flat array data back into state field interior cells
 *
 * Inlined for performance in hot paths.
 *
 * @param state State field (with ghost cells)
 * @param var_index Variable index to write to
 * @param data Input data array (size must be nx*ny*nz)
 * @param nx, ny, nz Interior cell dimensions (without ghosts)
 * @param nghost Number of ghost cells
 */
inline void insert_interior_cells(
    core::StateField3D& state,
    int var_index,
    const std::vector<double>& data,
    int nx, int ny, int nz,
    int nghost
) {
    if (data.size() != static_cast<size_t>(nx * ny * nz)) {
        throw std::runtime_error("Data size mismatch in insert_interior_cells");
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx = i * ny * nz + j * nz + k;
                state(var_index, i + nghost, j + nghost, k + nghost) = data[idx];
            }
        }
    }
}

/**
 * @brief Write a scalar attribute to HDF5 location
 *
 * @param loc_id HDF5 location (file, group, or dataset)
 * @param name Attribute name
 * @param value Attribute value (int, double, etc.)
 */
template<typename T>
void write_scalar_attribute(hid_t loc_id, const std::string& name, T value);

/**
 * @brief Read a scalar attribute from HDF5 location
 *
 * @param loc_id HDF5 location (file, group, or dataset)
 * @param name Attribute name
 * @return Attribute value
 */
template<typename T>
T read_scalar_attribute(hid_t loc_id, const std::string& name);

/**
 * @brief Write a string attribute to HDF5 location
 *
 * @param loc_id HDF5 location
 * @param name Attribute name
 * @param value String value
 */
void write_string_attribute(hid_t loc_id, const std::string& name, const std::string& value);

/**
 * @brief Read a string attribute from HDF5 location
 *
 * @param loc_id HDF5 location
 * @param name Attribute name
 * @return String value
 */
std::string read_string_attribute(hid_t loc_id, const std::string& name);

// Template implementations (must be in header for templates)

template<typename T>
void write_scalar_attribute(hid_t loc_id, const std::string& name, T value) {
    hid_t attr_space = H5Screate(H5S_SCALAR);

    // Determine HDF5 type based on T
    hid_t hdf5_type;
    hid_t file_type;

    if constexpr (std::is_same_v<T, int>) {
        hdf5_type = H5T_NATIVE_INT;
        file_type = H5T_STD_I32LE;
    } else if constexpr (std::is_same_v<T, double>) {
        hdf5_type = H5T_NATIVE_DOUBLE;
        file_type = H5T_IEEE_F64LE;
    } else if constexpr (std::is_same_v<T, float>) {
        hdf5_type = H5T_NATIVE_FLOAT;
        file_type = H5T_IEEE_F32LE;
    }

    hid_t attr = H5Acreate(loc_id, name.c_str(), file_type, attr_space,
                          H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, hdf5_type, &value);
    H5Aclose(attr);
    H5Sclose(attr_space);
}

template<typename T>
T read_scalar_attribute(hid_t loc_id, const std::string& name) {
    hid_t attr = H5Aopen(loc_id, name.c_str(), H5P_DEFAULT);

    // Determine HDF5 type based on T
    hid_t hdf5_type;
    if constexpr (std::is_same_v<T, int>) {
        hdf5_type = H5T_NATIVE_INT;
    } else if constexpr (std::is_same_v<T, double>) {
        hdf5_type = H5T_NATIVE_DOUBLE;
    } else if constexpr (std::is_same_v<T, float>) {
        hdf5_type = H5T_NATIVE_FLOAT;
    }

    T value;
    H5Aread(attr, hdf5_type, &value);
    H5Aclose(attr);

    return value;
}

} // namespace fvm3d::io::hdf5_utils

#endif // HDF5_UTILS_HPP
