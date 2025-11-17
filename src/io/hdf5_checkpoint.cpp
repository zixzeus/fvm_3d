#include "io/hdf5_checkpoint.hpp"
#include <hdf5.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <vector>
#include <sys/stat.h>

namespace fvm3d::io {

HDF5Checkpoint::HDF5Checkpoint(const std::string& filename)
    : filename_(filename) {
}

void HDF5Checkpoint::save(
    const StateField3D& state,
    const Grid3D& grid,
    double time,
    int step_count,
    const std::string& description
) {
    // Create HDF5 file
    hid_t file_id = H5Fcreate(
        filename_.c_str(),
        H5F_ACC_TRUNC,
        H5P_DEFAULT,
        H5P_DEFAULT
    );

    if (file_id < 0) {
        throw std::runtime_error("Failed to create HDF5 file: " + filename_);
    }

    try {
        // Save each component
        save_grid(reinterpret_cast<void*>(&file_id), grid);
        save_state(reinterpret_cast<void*>(&file_id), state);
        save_metadata(
            reinterpret_cast<void*>(&file_id),
            time,
            step_count,
            description
        );

        std::cout << "Checkpoint saved to " << filename_
                  << " (time=" << std::scientific << std::setprecision(6) << time
                  << ", step=" << step_count << ")\n";
    } catch (...) {
        H5Fclose(file_id);
        throw;
    }

    // Close file
    H5Fclose(file_id);
}

bool HDF5Checkpoint::load(
    StateField3D& state,
    double& time,
    int& step_count
) {
    // Open HDF5 file for reading
    hid_t file_id = H5Fopen(
        filename_.c_str(),
        H5F_ACC_RDONLY,
        H5P_DEFAULT
    );

    if (file_id < 0) {
        std::cerr << "Failed to open HDF5 file: " << filename_ << std::endl;
        return false;
    }

    try {
        std::string description;
        if (!load_metadata(
            reinterpret_cast<void*>(&file_id),
            time,
            step_count,
            description
        )) {
            H5Fclose(file_id);
            return false;
        }

        if (!load_state(reinterpret_cast<void*>(&file_id), state)) {
            H5Fclose(file_id);
            return false;
        }

        std::cout << "Checkpoint loaded from " << filename_
                  << " (time=" << std::scientific << std::setprecision(6) << time
                  << ", step=" << step_count << ")\n";

        H5Fclose(file_id);
        return true;

    } catch (...) {
        H5Fclose(file_id);
        return false;
    }
}

bool HDF5Checkpoint::file_exists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

bool HDF5Checkpoint::read_metadata(
    const std::string& filename,
    double& time,
    int& step_count,
    std::string& description
) {
    if (!file_exists(filename)) {
        return false;
    }

    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        return false;
    }

    HDF5Checkpoint checkpoint(filename);
    bool success = checkpoint.load_metadata(
        reinterpret_cast<void*>(&file_id),
        time,
        step_count,
        description
    );

    H5Fclose(file_id);
    return success;
}

void HDF5Checkpoint::save_grid(void* file_id_void, const Grid3D& grid) {
    hid_t file_id = *reinterpret_cast<hid_t*>(file_id_void);

    // Create grid group
    hid_t group_id = H5Gcreate(
        file_id,
        "/grid",
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT
    );

    if (group_id < 0) {
        throw std::runtime_error("Failed to create /grid group in HDF5 file");
    }

    try {
        const auto& geom = grid.geometry();

        // Save geometry as dataset (9 doubles: xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz)
        hsize_t dims_geom[1] = {9};
        hid_t space_geom = H5Screate_simple(1, dims_geom, nullptr);

        double geom_data[9] = {
            geom.xmin, geom.ymin, geom.zmin,
            geom.xmin + geom.Lx, geom.ymin + geom.Ly, geom.zmin + geom.Lz,
            geom.dx, geom.dy, geom.dz
        };

        hid_t dset_geom = H5Dcreate(
            group_id,
            "geometry",
            H5T_IEEE_F64LE,
            space_geom,
            H5P_DEFAULT,
            H5P_DEFAULT,
            H5P_DEFAULT
        );

        H5Dwrite(dset_geom, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, geom_data);
        H5Dclose(dset_geom);
        H5Sclose(space_geom);

        // Save grid dimensions as attributes
        hid_t attr_space = H5Screate(H5S_SCALAR);
        int int_attr;

        int_attr = grid.nx_local();
        hid_t attr_nx = H5Acreate(group_id, "nx", H5T_STD_I32LE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_nx, H5T_NATIVE_INT, &int_attr);
        H5Aclose(attr_nx);

        int_attr = grid.ny_local();
        hid_t attr_ny = H5Acreate(group_id, "ny", H5T_STD_I32LE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_ny, H5T_NATIVE_INT, &int_attr);
        H5Aclose(attr_ny);

        int_attr = grid.nz_local();
        hid_t attr_nz = H5Acreate(group_id, "nz", H5T_STD_I32LE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_nz, H5T_NATIVE_INT, &int_attr);
        H5Aclose(attr_nz);

        int_attr = grid.nghost();
        hid_t attr_nghost = H5Acreate(group_id, "nghost", H5T_STD_I32LE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_nghost, H5T_NATIVE_INT, &int_attr);
        H5Aclose(attr_nghost);

        H5Sclose(attr_space);

    } catch (...) {
        H5Gclose(group_id);
        throw;
    }

    H5Gclose(group_id);
}

void HDF5Checkpoint::save_state(void* file_id_void, const StateField3D& state) {
    hid_t file_id = *reinterpret_cast<hid_t*>(file_id_void);

    // Create state group
    hid_t group_id = H5Gcreate(
        file_id,
        "/state",
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT
    );

    if (group_id < 0) {
        throw std::runtime_error("Failed to create /state group in HDF5 file");
    }

    try {
        int nx = state.nx();
        int ny = state.ny();
        int nz = state.nz();

        hsize_t dims[3] = {(hsize_t)nx, (hsize_t)ny, (hsize_t)nz};
        hid_t space_id = H5Screate_simple(3, dims, nullptr);

        const char* var_names[5] = {"rho", "rho_u", "rho_v", "rho_w", "E"};

        // Save each variable
        for (int v = 0; v < 5; v++) {
            hid_t dset_id = H5Dcreate(
                group_id,
                var_names[v],
                H5T_IEEE_F64LE,
                space_id,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT
            );

            // Write data for this variable
            std::vector<double> var_data(nx * ny * nz);
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    for (int k = 0; k < nz; k++) {
                        var_data[i * ny * nz + j * nz + k] = state(v, i, j, k);
                    }
                }
            }

            H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, var_data.data());
            H5Dclose(dset_id);
        }

        H5Sclose(space_id);

    } catch (...) {
        H5Gclose(group_id);
        throw;
    }

    H5Gclose(group_id);
}

void HDF5Checkpoint::save_metadata(
    void* file_id_void,
    double time,
    int step_count,
    const std::string& description
) {
    hid_t file_id = *reinterpret_cast<hid_t*>(file_id_void);

    // Create metadata group
    hid_t group_id = H5Gcreate(
        file_id,
        "/metadata",
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT
    );

    if (group_id < 0) {
        throw std::runtime_error("Failed to create /metadata group in HDF5 file");
    }

    try {
        hid_t attr_space = H5Screate(H5S_SCALAR);

        // Save time
        hid_t attr_time = H5Acreate(group_id, "time", H5T_IEEE_F64LE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_time, H5T_NATIVE_DOUBLE, &time);
        H5Aclose(attr_time);

        // Save step count
        hid_t attr_step = H5Acreate(group_id, "step_count", H5T_STD_I32LE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_step, H5T_NATIVE_INT, &step_count);
        H5Aclose(attr_step);

        // Save description (if not empty)
        if (!description.empty()) {
            hid_t str_type = H5Tcopy(H5T_C_S1);
            H5Tset_size(str_type, description.length());
            hid_t attr_desc = H5Acreate(group_id, "description", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(attr_desc, str_type, description.c_str());
            H5Aclose(attr_desc);
            H5Tclose(str_type);
        }

        H5Sclose(attr_space);

    } catch (...) {
        H5Gclose(group_id);
        throw;
    }

    H5Gclose(group_id);
}

bool HDF5Checkpoint::load_grid(void* file_id_void) {
    hid_t file_id = *reinterpret_cast<hid_t*>(file_id_void);

    hid_t group_id = H5Gopen(file_id, "/grid", H5P_DEFAULT);
    if (group_id < 0) {
        std::cerr << "Failed to open /grid group\n";
        return false;
    }

    try {
        // Grid info is stored for reference but not loaded
        // (Grid is already initialized in FVMSolver3D)
        H5Gclose(group_id);
        return true;
    } catch (...) {
        H5Gclose(group_id);
        return false;
    }
}

bool HDF5Checkpoint::load_state(void* file_id_void, StateField3D& state) {
    hid_t file_id = *reinterpret_cast<hid_t*>(file_id_void);

    hid_t group_id = H5Gopen(file_id, "/state", H5P_DEFAULT);
    if (group_id < 0) {
        std::cerr << "Failed to open /state group\n";
        return false;
    }

    try {
        const char* var_names[5] = {"rho", "rho_u", "rho_v", "rho_w", "E"};
        int nx = state.nx();
        int ny = state.ny();
        int nz = state.nz();

        for (int v = 0; v < 5; v++) {
            hid_t dset_id = H5Dopen(group_id, var_names[v], H5P_DEFAULT);
            if (dset_id < 0) {
                std::cerr << "Failed to open dataset " << var_names[v] << "\n";
                H5Gclose(group_id);
                return false;
            }

            std::vector<double> var_data(nx * ny * nz);
            H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, var_data.data());

            // Copy data back to state field
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    for (int k = 0; k < nz; k++) {
                        state(v, i, j, k) = var_data[i * ny * nz + j * nz + k];
                    }
                }
            }

            H5Dclose(dset_id);
        }

        H5Gclose(group_id);
        return true;

    } catch (...) {
        H5Gclose(group_id);
        return false;
    }
}

bool HDF5Checkpoint::load_metadata(
    void* file_id_void,
    double& time,
    int& step_count,
    std::string& description
) {
    hid_t file_id = *reinterpret_cast<hid_t*>(file_id_void);

    hid_t group_id = H5Gopen(file_id, "/metadata", H5P_DEFAULT);
    if (group_id < 0) {
        std::cerr << "Failed to open /metadata group\n";
        return false;
    }

    try {
        // Read time
        hid_t attr_time = H5Aopen(group_id, "time", H5P_DEFAULT);
        if (attr_time < 0) {
            std::cerr << "Failed to open time attribute\n";
            H5Gclose(group_id);
            return false;
        }
        H5Aread(attr_time, H5T_NATIVE_DOUBLE, &time);
        H5Aclose(attr_time);

        // Read step_count
        hid_t attr_step = H5Aopen(group_id, "step_count", H5P_DEFAULT);
        if (attr_step < 0) {
            std::cerr << "Failed to open step_count attribute\n";
            H5Gclose(group_id);
            return false;
        }
        H5Aread(attr_step, H5T_NATIVE_INT, &step_count);
        H5Aclose(attr_step);

        // Read description (optional)
        hid_t attr_desc = H5Aopen(group_id, "description", H5P_DEFAULT);
        if (attr_desc >= 0) {
            hid_t desc_type = H5Aget_type(attr_desc);
            size_t str_len = H5Tget_size(desc_type);
            char* desc_buf = new char[str_len + 1];
            H5Aread(attr_desc, desc_type, desc_buf);
            desc_buf[str_len] = '\0';
            description = std::string(desc_buf);
            delete[] desc_buf;
            H5Tclose(desc_type);
            H5Aclose(attr_desc);
        }

        H5Gclose(group_id);
        return true;

    } catch (...) {
        H5Gclose(group_id);
        return false;
    }
}

} // namespace fvm3d::io
