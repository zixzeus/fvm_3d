#pragma once

#include "core/grid3d.hpp"
#include "core/field3d.hpp"
#include "parallel/mpi_domain_decomposer.hpp"
#include "io/vtk_writer.hpp"
#include <string>
#include <mpi.h>

namespace fvm3d::io {

/**
 * Parallel VTK writer for MPI-distributed 3D FVM data.
 *
 * Creates:
 * - Individual .vti files for each MPI process (subdomain)
 * - A master .pvti file (parallel VTK) that references all subdomains
 *
 * The .pvti file can be opened directly in ParaView to visualize
 * the entire parallel dataset as a single unified grid.
 *
 * Usage (in MPI context):
 *   ParallelVTKWriter writer("output", decomposer);
 *   writer.set_local_grid(local_grid);
 *   writer.add_scalar_field("density", rho_data, nx, ny, nz);
 *   writer.write_parallel();  // Each rank writes its subdomain + root writes .pvti
 */
class ParallelVTKWriter {
public:
    /**
     * Constructor.
     * @param base_filename: Base name without extension (e.g., "output")
     * @param decomposer: MPI domain decomposer (for subdomain info)
     * @param comm: MPI communicator (default: MPI_COMM_WORLD)
     */
    ParallelVTKWriter(
        const std::string& base_filename,
        const parallel::MPIDomainDecomposer& decomposer,
        MPI_Comm comm = MPI_COMM_WORLD
    );

    /**
     * Set local grid geometry for this process.
     */
    void set_local_grid(const core::Grid3D& local_grid);

    /**
     * Add scalar field (same interface as VTKWriter).
     */
    void add_scalar_field(
        const std::string& name,
        const double* data,
        int nx, int ny, int nz
    );

    /**
     * Add vector field (same interface as VTKWriter).
     */
    void add_vector_field(
        const std::string& name,
        const double* data_x,
        const double* data_y,
        const double* data_z,
        int nx, int ny, int nz
    );

    /**
     * Add scalar field from StateField3D.
     */
    void add_scalar_from_state(
        const std::string& name,
        const core::StateField3D& field,
        int var_index,
        const core::Grid3D& grid
    );

    /**
     * Add vector field from StateField3D.
     */
    void add_vector_from_state(
        const std::string& name,
        const core::StateField3D& field,
        int var_x, int var_y, int var_z,
        const core::Grid3D& grid
    );

    /**
     * Write parallel VTK files.
     * Each rank writes its .vti file, rank 0 writes the .pvti master file.
     */
    void write_parallel();

    /**
     * Clear all fields.
     */
    void clear_fields();

private:
    std::string base_filename_;
    const parallel::MPIDomainDecomposer& decomposer_;
    MPI_Comm comm_;
    int rank_, size_;

    VTKWriter local_writer_;  // Writer for this process's subdomain
    bool grid_set_;

    // Global grid information (from decomposer)
    int global_nx_, global_ny_, global_nz_;
    double global_xmin_, global_ymin_, global_zmin_;
    double global_dx_, global_dy_, global_dz_;

    // Field names (tracked for .pvti metadata)
    std::vector<std::string> scalar_field_names_;
    std::vector<std::string> vector_field_names_;

    /**
     * Generate filename for this rank's subdomain.
     * Format: base_filename_rank000.vti
     */
    std::string get_piece_filename(int rank) const;

    /**
     * Write the master .pvti file (rank 0 only).
     */
    void write_master_file();

    /**
     * Get extent for a specific rank.
     */
    struct Extent {
        int i_min, i_max;
        int j_min, j_max;
        int k_min, k_max;
    };
    Extent get_rank_extent(int rank) const;
};

/**
 * Utility function: Export MPI-distributed state to parallel VTK.
 *
 * @param base_filename: Base filename (without extension)
 * @param state: Local conservative state field
 * @param local_grid: Local grid geometry
 * @param decomposer: Domain decomposer
 * @param num_vars: Number of variables
 * @param physics_type: "euler" or "mhd"
 * @param time: Current simulation time
 * @param comm: MPI communicator
 */
void export_parallel_state_to_vtk(
    const std::string& base_filename,
    const core::StateField3D& state,
    const core::Grid3D& local_grid,
    const parallel::MPIDomainDecomposer& decomposer,
    int num_vars,
    const std::string& physics_type,
    double time = 0.0,
    MPI_Comm comm = MPI_COMM_WORLD
);

} // namespace fvm3d::io
