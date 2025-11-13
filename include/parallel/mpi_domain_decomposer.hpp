#pragma once

#include <mpi.h>
#include <array>
#include <vector>
#include <cstdint>

namespace fvm3d::parallel {

/**
 * Domain decomposition for MPI parallelization.
 *
 * Divides 3D computational domain into rectangular subdomains using
 * Cartesian process topology (rank → (px, py, pz)).
 *
 * Example: 8 processes decomposed as 2×2×2:
 * ┌─────────────────────────────────┐
 * │ Rank layout (X-Y plane, Z fixed)│
 * ├─────────────────────────────────┤
 * │  0 │ 1                          │
 * ├─────┼──────                     │
 * │  2 │ 3       (Y direction)      │
 * └─────┴──────────────────────────┘
 */
class MPIDomainDecomposer {
public:
    /**
     * Create Cartesian domain decomposition.
     *
     * @param nx_global, ny_global, nz_global: Total grid dimensions
     * @param np_x, np_y, np_z: Number of processes in each direction
     *   If 0, automatic determination based on total process count
     */
    MPIDomainDecomposer(
        int nx_global, int ny_global, int nz_global,
        int np_x = 0, int np_y = 0, int np_z = 0
    );

    /**
     * Get the Cartesian communicator.
     */
    MPI_Comm cartesian_comm() const { return cartesian_comm_; }

    /**
     * Get process grid dimensions.
     */
    std::array<int, 3> process_grid() const {
        return {np_x_, np_y_, np_z_};
    }

    /**
     * Get coordinates of current process in 3D process grid.
     */
    std::array<int, 3> process_coords() const {
        return {my_x_, my_y_, my_z_};
    }

    /**
     * Get local grid dimensions (interior cells only, excluding ghosts).
     */
    std::array<int, 3> local_interior_cells() const {
        return {nx_local_, ny_local_, nz_local_};
    }

    /**
     * Get local grid dimensions (including ghost cells).
     */
    std::array<int, 3> local_total_cells() const {
        int nghost = 2;  // Standard ghost cell width
        return {nx_local_ + 2*nghost, ny_local_ + 2*nghost, nz_local_ + 2*nghost};
    }

    /**
     * Get global grid dimensions.
     */
    std::array<int, 3> global_grid() const {
        return {nx_global_, ny_global_, nz_global_};
    }

    /**
     * Get global index range for this process.
     * @return: {{i_min, i_max}, {j_min, j_max}, {k_min, k_max}}
     */
    struct IndexRange {
        int i_min, i_max;
        int j_min, j_max;
        int k_min, k_max;
    };
    IndexRange global_index_range() const;

    /**
     * Get neighbor process ranks.
     * @return: {{x_minus, x_plus}, {y_minus, y_plus}, {z_minus, z_plus}}
     */
    struct Neighbors {
        int x_minus, x_plus;
        int y_minus, y_plus;
        int z_minus, z_plus;
    };
    Neighbors neighbors() const;

    /**
     * Check if this process has a neighbor in given direction.
     * @param direction: 0=X, 1=Y, 2=Z
     * @param side: -1=minus, +1=plus
     */
    bool has_neighbor(int direction, int side) const;

    /**
     * Get neighbor rank in given direction.
     * Returns -1 if no neighbor (boundary).
     */
    int neighbor_rank(int direction, int side) const;

    /**
     * Check if process is on boundary in given direction.
     */
    bool is_boundary(int direction) const;

    /**
     * Print decomposition info (root process only).
     */
    void print_info() const;

private:
    MPI_Comm cartesian_comm_;

    // Global dimensions
    int nx_global_, ny_global_, nz_global_;

    // Process grid dimensions
    int np_x_, np_y_, np_z_;

    // Local process coordinates
    int my_x_, my_y_, my_z_;

    // Local domain sizes (interior cells)
    int nx_local_, ny_local_, nz_local_;

    /**
     * Compute local cell count for a process coordinate.
     * Balances load by distributing remainder cells.
     */
    static int compute_local_size(
        int global_size, int num_procs, int my_rank
    );

    /**
     * Create Cartesian communicator and determine rank layout.
     */
    void create_cartesian_topology();
};

} // namespace fvm3d::parallel
