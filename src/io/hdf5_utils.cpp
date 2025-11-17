#include "io/hdf5_utils.hpp"
#include <stdexcept>

namespace fvm3d::io::hdf5_utils {

const std::vector<std::string>& get_variable_names(int num_vars) {
    // Use static const to avoid repeated allocations
    static const std::vector<std::string> vars_5 = {"rho", "rho_u", "rho_v", "rho_w", "E"};
    static const std::vector<std::string> vars_8 = {"rho", "rho_u", "rho_v", "rho_w", "E", "Bx", "By", "Bz"};
    static const std::vector<std::string> vars_9 = {"rho", "rho_u", "rho_v", "rho_w", "E", "Bx", "By", "Bz", "psi"};

    if (num_vars == 5) {
        return vars_5;
    } else if (num_vars == 8) {
        return vars_8;
    } else if (num_vars == 9) {
        return vars_9;
    } else {
        throw std::runtime_error("Unsupported number of variables: " + std::to_string(num_vars));
    }
}

// extract_interior_cells and insert_interior_cells are now inline in the header

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
