#include "io/hdf5_utils.hpp"
#include <stdexcept>

namespace fvm3d::io::hdf5_utils {

std::vector<std::string> get_variable_names(int num_vars) {
    if (num_vars == 5) {
        return {"rho", "rho_u", "rho_v", "rho_w", "E"};
    } else if (num_vars == 8) {
        return {"rho", "rho_u", "rho_v", "rho_w", "E", "Bx", "By", "Bz"};
    } else if (num_vars == 9) {
        return {"rho", "rho_u", "rho_v", "rho_w", "E", "Bx", "By", "Bz", "psi"};
    } else {
        throw std::runtime_error("Unsupported number of variables: " + std::to_string(num_vars));
    }
}

void extract_interior_cells(
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

void insert_interior_cells(
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

void write_string_attribute(hid_t loc_id, const std::string& name, const std::string& value) {
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t str_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_type, value.length());

    hid_t attr = H5Acreate(loc_id, name.c_str(), str_type, attr_space,
                          H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, str_type, value.c_str());

    H5Aclose(attr);
    H5Tclose(str_type);
    H5Sclose(attr_space);
}

std::string read_string_attribute(hid_t loc_id, const std::string& name) {
    hid_t attr = H5Aopen(loc_id, name.c_str(), H5P_DEFAULT);
    hid_t attr_type = H5Aget_type(attr);
    size_t str_size = H5Tget_size(attr_type);

    std::vector<char> buffer(str_size + 1, '\0');
    H5Aread(attr, attr_type, buffer.data());

    H5Tclose(attr_type);
    H5Aclose(attr);

    return std::string(buffer.data());
}

} // namespace fvm3d::io::hdf5_utils
