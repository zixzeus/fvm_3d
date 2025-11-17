#include "io/mpi_hdf5_checkpoint.hpp"
#include <iostream>
#include <stdexcept>
#include <cstring>

namespace fvm3d::io {

void MPIHDFCheckpoint::save(
    const std::string& filename,
    const core::StateField3D& state,
    const core::Grid3D& grid,
    const parallel::MPIDomainDecomposer& decomposer,
    double time,
    int step_count,
    const std::string& description,
    MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    // Create parallel HDF5 file
    hid_t file_id = create_parallel_file(filename, comm);

    // Get dimensions
    auto global_grid = decomposer.global_grid();
    int global_nx = global_grid[0];
    int global_ny = global_grid[1];
    int global_nz = global_grid[2];

    auto local_cells = decomposer.local_interior_cells();
    int local_nx = local_cells[0];
    int local_ny = local_cells[1];
    int local_nz = local_cells[2];

    auto idx_range = decomposer.global_index_range();
    int offset_x = idx_range.i_min;
    int offset_y = idx_range.j_min;
    int offset_z = idx_range.k_min;

    int num_vars = state.nvars();

    // Write grid metadata (root only)
    if (rank == 0) {
        hid_t grid_group = H5Gcreate(file_id, "/grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        write_attribute(grid_group, "global_nx", global_nx);
        write_attribute(grid_group, "global_ny", global_ny);
        write_attribute(grid_group, "global_nz", global_nz);

        const auto& geom = grid.geometry();
        write_attribute(grid_group, "xmin", geom.xmin);
        write_attribute(grid_group, "ymin", geom.ymin);
        write_attribute(grid_group, "zmin", geom.zmin);
        write_attribute(grid_group, "Lx", geom.Lx);
        write_attribute(grid_group, "Ly", geom.Ly);
        write_attribute(grid_group, "Lz", geom.Lz);

        H5Gclose(grid_group);
    }

    // Write decomposition info (root only)
    if (rank == 0) {
        hid_t decomp_group = H5Gcreate(file_id, "/decomposition", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        auto proc_grid = decomposer.process_grid();
        write_attribute(decomp_group, "px", proc_grid[0]);
        write_attribute(decomp_group, "py", proc_grid[1]);
        write_attribute(decomp_group, "pz", proc_grid[2]);

        int nprocs;
        MPI_Comm_size(comm, &nprocs);
        write_attribute(decomp_group, "nprocs", nprocs);

        H5Gclose(decomp_group);
    }

    // Synchronize before writing state
    MPI_Barrier(comm);

    // Create state group
    hid_t state_group;
    if (rank == 0) {
        state_group = H5Gcreate(file_id, "/state", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(state_group);
    }
    MPI_Barrier(comm);

    // Write each variable as a 3D dataset using collective I/O
    std::vector<std::string> var_names;
    if (num_vars == 5) {
        var_names = {"rho", "rho_u", "rho_v", "rho_w", "E"};
    } else if (num_vars == 8) {
        var_names = {"rho", "rho_u", "rho_v", "rho_w", "E", "Bx", "By", "Bz"};
    } else if (num_vars == 9) {
        var_names = {"rho", "rho_u", "rho_v", "rho_w", "E", "Bx", "By", "Bz", "psi"};
    } else {
        throw std::runtime_error("Unsupported number of variables: " + std::to_string(num_vars));
    }

    int nghost = grid.nghost();

    for (int v = 0; v < num_vars; v++) {
        // Extract interior cells only (no ghost cells)
        std::vector<double> local_data(local_nx * local_ny * local_nz);

        for (int i = 0; i < local_nx; i++) {
            for (int j = 0; j < local_ny; j++) {
                for (int k = 0; k < local_nz; k++) {
                    int idx = i * local_ny * local_nz + j * local_nz + k;
                    local_data[idx] = state(v, i + nghost, j + nghost, k + nghost);
                }
            }
        }

        // Write dataset
        std::string dataset_name = "/state/" + var_names[v];
        write_3d_dataset(file_id, dataset_name, local_data.data(),
                        global_nx, global_ny, global_nz,
                        local_nx, local_ny, local_nz,
                        offset_x, offset_y, offset_z);
    }

    // Write metadata (root only)
    if (rank == 0) {
        hid_t meta_group = H5Gcreate(file_id, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        write_attribute(meta_group, "time", time);
        write_attribute(meta_group, "step_count", step_count);
        write_attribute(meta_group, "num_vars", num_vars);

        if (!description.empty()) {
            write_string_attribute(meta_group, "description", description);
        }

        H5Gclose(meta_group);
    }

    // Close file
    H5Fclose(file_id);

    if (rank == 0) {
        std::cout << "Parallel checkpoint saved: " << filename << "\n";
    }
}

bool MPIHDFCheckpoint::load(
    const std::string& filename,
    core::StateField3D& state,
    const core::Grid3D& grid,
    const parallel::MPIDomainDecomposer& decomposer,
    double& time,
    int& step_count,
    MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    // Open parallel HDF5 file
    hid_t file_id = open_parallel_file(filename, comm);
    if (file_id < 0) {
        if (rank == 0) {
            std::cerr << "Failed to open checkpoint file: " << filename << "\n";
        }
        return false;
    }

    // Read metadata (all ranks)
    hid_t meta_group = H5Gopen(file_id, "/metadata", H5P_DEFAULT);
    time = read_attribute<double>(meta_group, "time");
    step_count = read_attribute<int>(meta_group, "step_count");
    int num_vars = read_attribute<int>(meta_group, "num_vars");
    H5Gclose(meta_group);

    // Verify number of variables matches
    if (num_vars != state.nvars()) {
        if (rank == 0) {
            std::cerr << "Variable count mismatch: file has " << num_vars
                      << ", state has " << state.nvars() << "\n";
        }
        H5Fclose(file_id);
        return false;
    }

    // Get dimensions
    auto global_grid = decomposer.global_grid();
    int global_nx = global_grid[0];
    int global_ny = global_grid[1];
    int global_nz = global_grid[2];

    auto local_cells = decomposer.local_interior_cells();
    int local_nx = local_cells[0];
    int local_ny = local_cells[1];
    int local_nz = local_cells[2];

    auto idx_range = decomposer.global_index_range();
    int offset_x = idx_range.i_min;
    int offset_y = idx_range.j_min;
    int offset_z = idx_range.k_min;

    // Variable names
    std::vector<std::string> var_names;
    if (num_vars == 5) {
        var_names = {"rho", "rho_u", "rho_v", "rho_w", "E"};
    } else if (num_vars == 8) {
        var_names = {"rho", "rho_u", "rho_v", "rho_w", "E", "Bx", "By", "Bz"};
    } else if (num_vars == 9) {
        var_names = {"rho", "rho_u", "rho_v", "rho_w", "E", "Bx", "By", "Bz", "psi"};
    }

    int nghost = grid.nghost();

    // Read each variable
    for (int v = 0; v < num_vars; v++) {
        std::vector<double> local_data(local_nx * local_ny * local_nz);

        std::string dataset_name = "/state/" + var_names[v];
        read_3d_dataset(file_id, dataset_name, local_data.data(),
                       global_nx, global_ny, global_nz,
                       local_nx, local_ny, local_nz,
                       offset_x, offset_y, offset_z);

        // Copy to state (interior cells)
        for (int i = 0; i < local_nx; i++) {
            for (int j = 0; j < local_ny; j++) {
                for (int k = 0; k < local_nz; k++) {
                    int idx = i * local_ny * local_nz + j * local_nz + k;
                    state(v, i + nghost, j + nghost, k + nghost) = local_data[idx];
                }
            }
        }
    }

    H5Fclose(file_id);

    if (rank == 0) {
        std::cout << "Parallel checkpoint loaded: " << filename << "\n";
    }

    return true;
}

bool MPIHDFCheckpoint::read_metadata(
    const std::string& filename,
    double& time,
    int& step_count,
    int& num_vars,
    std::string& description,
    MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
        // Open file (serial access for metadata read)
        hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) {
            std::cerr << "Failed to open file: " << filename << "\n";
            return false;
        }

        // Read metadata
        hid_t meta_group = H5Gopen(file_id, "/metadata", H5P_DEFAULT);
        time = read_attribute<double>(meta_group, "time");
        step_count = read_attribute<int>(meta_group, "step_count");
        num_vars = read_attribute<int>(meta_group, "num_vars");

        // Try to read description (optional)
        htri_t exists = H5Aexists(meta_group, "description");
        if (exists > 0) {
            description = read_string_attribute(meta_group, "description");
        } else {
            description = "";
        }

        H5Gclose(meta_group);
        H5Fclose(file_id);
    }

    // Broadcast metadata to all ranks
    MPI_Bcast(&time, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&step_count, 1, MPI_INT, 0, comm);
    MPI_Bcast(&num_vars, 1, MPI_INT, 0, comm);

    return true;
}

hid_t MPIHDFCheckpoint::create_parallel_file(const std::string& filename, MPI_Comm comm) {
    // Set up parallel file access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);

    // Create file
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    return file_id;
}

hid_t MPIHDFCheckpoint::open_parallel_file(const std::string& filename, MPI_Comm comm) {
    // Set up parallel file access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);

    // Open file
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    return file_id;
}

void MPIHDFCheckpoint::write_3d_dataset(
    hid_t file_id,
    const std::string& dataset_name,
    const double* data,
    int global_nx, int global_ny, int global_nz,
    int local_nx, int local_ny, int local_nz,
    int offset_x, int offset_y, int offset_z)
{
    // Create dataspace for global dataset
    hsize_t global_dims[3] = {
        static_cast<hsize_t>(global_nx),
        static_cast<hsize_t>(global_ny),
        static_cast<hsize_t>(global_nz)
    };
    hid_t global_space = H5Screate_simple(3, global_dims, NULL);

    // Create dataset
    hid_t dataset = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE,
                              global_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create dataspace for local chunk
    hsize_t local_dims[3] = {
        static_cast<hsize_t>(local_nx),
        static_cast<hsize_t>(local_ny),
        static_cast<hsize_t>(local_nz)
    };
    hid_t local_space = H5Screate_simple(3, local_dims, NULL);

    // Select hyperslab in global dataspace
    hsize_t start[3] = {
        static_cast<hsize_t>(offset_x),
        static_cast<hsize_t>(offset_y),
        static_cast<hsize_t>(offset_z)
    };
    hsize_t count[3] = {
        static_cast<hsize_t>(local_nx),
        static_cast<hsize_t>(local_ny),
        static_cast<hsize_t>(local_nz)
    };
    H5Sselect_hyperslab(global_space, H5S_SELECT_SET, start, NULL, count, NULL);

    // Set up collective I/O
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write data
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, local_space, global_space, plist_id, data);

    // Cleanup
    H5Pclose(plist_id);
    H5Sclose(local_space);
    H5Sclose(global_space);
    H5Dclose(dataset);
}

void MPIHDFCheckpoint::read_3d_dataset(
    hid_t file_id,
    const std::string& dataset_name,
    double* data,
    int global_nx, int global_ny, int global_nz,
    int local_nx, int local_ny, int local_nz,
    int offset_x, int offset_y, int offset_z)
{
    // Open dataset
    hid_t dataset = H5Dopen(file_id, dataset_name.c_str(), H5P_DEFAULT);
    hid_t global_space = H5Dget_space(dataset);

    // Create dataspace for local chunk
    hsize_t local_dims[3] = {
        static_cast<hsize_t>(local_nx),
        static_cast<hsize_t>(local_ny),
        static_cast<hsize_t>(local_nz)
    };
    hid_t local_space = H5Screate_simple(3, local_dims, NULL);

    // Select hyperslab in global dataspace
    hsize_t start[3] = {
        static_cast<hsize_t>(offset_x),
        static_cast<hsize_t>(offset_y),
        static_cast<hsize_t>(offset_z)
    };
    hsize_t count[3] = {
        static_cast<hsize_t>(local_nx),
        static_cast<hsize_t>(local_ny),
        static_cast<hsize_t>(local_nz)
    };
    H5Sselect_hyperslab(global_space, H5S_SELECT_SET, start, NULL, count, NULL);

    // Set up collective I/O
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Read data
    H5Dread(dataset, H5T_NATIVE_DOUBLE, local_space, global_space, plist_id, data);

    // Cleanup
    H5Pclose(plist_id);
    H5Sclose(local_space);
    H5Sclose(global_space);
    H5Dclose(dataset);
}

template<typename T>
void MPIHDFCheckpoint::write_attribute(hid_t loc_id, const std::string& name, T value) {
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t type;

    if constexpr (std::is_same_v<T, int>) {
        type = H5T_NATIVE_INT;
    } else if constexpr (std::is_same_v<T, double>) {
        type = H5T_NATIVE_DOUBLE;
    } else {
        throw std::runtime_error("Unsupported attribute type");
    }

    hid_t attr = H5Acreate(loc_id, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, &value);
    H5Aclose(attr);
    H5Sclose(space);
}

template<typename T>
T MPIHDFCheckpoint::read_attribute(hid_t loc_id, const std::string& name) {
    hid_t attr = H5Aopen(loc_id, name.c_str(), H5P_DEFAULT);
    T value;

    hid_t type;
    if constexpr (std::is_same_v<T, int>) {
        type = H5T_NATIVE_INT;
    } else if constexpr (std::is_same_v<T, double>) {
        type = H5T_NATIVE_DOUBLE;
    }

    H5Aread(attr, type, &value);
    H5Aclose(attr);
    return value;
}

void MPIHDFCheckpoint::write_string_attribute(hid_t loc_id, const std::string& name, const std::string& value) {
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, value.length() + 1);

    hid_t attr = H5Acreate(loc_id, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, value.c_str());
    H5Aclose(attr);
    H5Tclose(type);
    H5Sclose(space);
}

std::string MPIHDFCheckpoint::read_string_attribute(hid_t loc_id, const std::string& name) {
    hid_t attr = H5Aopen(loc_id, name.c_str(), H5P_DEFAULT);
    hid_t type = H5Aget_type(attr);
    size_t size = H5Tget_size(type);

    char* buffer = new char[size + 1];
    H5Aread(attr, type, buffer);
    buffer[size] = '\0';

    std::string result(buffer);
    delete[] buffer;

    H5Tclose(type);
    H5Aclose(attr);
    return result;
}

// Explicit template instantiations
template void MPIHDFCheckpoint::write_attribute<int>(hid_t, const std::string&, int);
template void MPIHDFCheckpoint::write_attribute<double>(hid_t, const std::string&, double);
template int MPIHDFCheckpoint::read_attribute<int>(hid_t, const std::string&);
template double MPIHDFCheckpoint::read_attribute<double>(hid_t, const std::string&);

} // namespace fvm3d::io
