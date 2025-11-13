#pragma once

#include "core/grid3d.hpp"
#include "core/field3d.hpp"
#include "parallel/mpi_domain_decomposer.hpp"
#include <string>
#include <hdf5.h>
#include <mpi.h>

namespace fvm3d::io {

/**
 * Parallel HDF5 checkpoint/restart for MPI-distributed data.
 *
 * Features:
 * - Each MPI rank writes its local subdomain
 * - Uses collective I/O for performance
 * - Stores grid geometry and decomposition info
 * - Supports restart from different process counts (TODO)
 *
 * File structure:
 * /grid/
 *   - global_nx, global_ny, global_nz (attributes)
 *   - xmin, ymin, zmin, Lx, Ly, Lz (attributes)
 * /state/
 *   - rho, rho_u, rho_v, rho_w, E, [Bx, By, Bz, psi] (3D datasets)
 * /metadata/
 *   - time, step_count, num_vars (attributes)
 *   - description (optional string)
 * /decomposition/
 *   - px, py, pz (process grid)
 *   - nprocs (total number of processes)
 */
class MPIHDFCheckpoint {
public:
    /**
     * Save checkpoint from MPI-distributed state.
     *
     * @param filename: Output HDF5 file path
     * @param state: Local state field (this rank's subdomain)
     * @param grid: Local grid information
     * @param decomposer: Domain decomposition info
     * @param time: Current simulation time
     * @param step_count: Current step number
     * @param description: Optional description string
     * @param comm: MPI communicator (default: MPI_COMM_WORLD)
     */
    static void save(
        const std::string& filename,
        const core::StateField3D& state,
        const core::Grid3D& grid,
        const parallel::MPIDomainDecomposer& decomposer,
        double time,
        int step_count,
        const std::string& description = "",
        MPI_Comm comm = MPI_COMM_WORLD
    );

    /**
     * Load checkpoint to MPI-distributed state.
     *
     * @param filename: Input HDF5 file path
     * @param state: Local state field to fill (must be pre-allocated)
     * @param grid: Local grid information
     * @param decomposer: Domain decomposition info (must match file)
     * @param time: Output simulation time
     * @param step_count: Output step number
     * @param comm: MPI communicator (default: MPI_COMM_WORLD)
     * @return true if successful, false otherwise
     */
    static bool load(
        const std::string& filename,
        core::StateField3D& state,
        const core::Grid3D& grid,
        const parallel::MPIDomainDecomposer& decomposer,
        double& time,
        int& step_count,
        MPI_Comm comm = MPI_COMM_WORLD
    );

    /**
     * Read metadata from checkpoint file (root process only).
     *
     * @param filename: Input HDF5 file path
     * @param time: Output simulation time
     * @param step_count: Output step number
     * @param num_vars: Output number of variables
     * @param description: Output description string
     * @param comm: MPI communicator
     * @return true if successful, false otherwise
     */
    static bool read_metadata(
        const std::string& filename,
        double& time,
        int& step_count,
        int& num_vars,
        std::string& description,
        MPI_Comm comm = MPI_COMM_WORLD
    );

private:
    /**
     * Create HDF5 file with parallel access.
     */
    static hid_t create_parallel_file(const std::string& filename, MPI_Comm comm);

    /**
     * Open HDF5 file with parallel access.
     */
    static hid_t open_parallel_file(const std::string& filename, MPI_Comm comm);

    /**
     * Write 3D dataset using collective I/O.
     */
    static void write_3d_dataset(
        hid_t file_id,
        const std::string& dataset_name,
        const double* data,
        int global_nx, int global_ny, int global_nz,
        int local_nx, int local_ny, int local_nz,
        int offset_x, int offset_y, int offset_z
    );

    /**
     * Read 3D dataset using collective I/O.
     */
    static void read_3d_dataset(
        hid_t file_id,
        const std::string& dataset_name,
        double* data,
        int global_nx, int global_ny, int global_nz,
        int local_nx, int local_ny, int local_nz,
        int offset_x, int offset_y, int offset_z
    );

    /**
     * Write scalar attribute (root process only).
     */
    template<typename T>
    static void write_attribute(hid_t loc_id, const std::string& name, T value);

    /**
     * Read scalar attribute.
     */
    template<typename T>
    static T read_attribute(hid_t loc_id, const std::string& name);

    /**
     * Write string attribute (root process only).
     */
    static void write_string_attribute(hid_t loc_id, const std::string& name, const std::string& value);

    /**
     * Read string attribute.
     */
    static std::string read_string_attribute(hid_t loc_id, const std::string& name);
};

} // namespace fvm3d::io
