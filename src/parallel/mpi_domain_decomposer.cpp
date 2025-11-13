#include "parallel/mpi_domain_decomposer.hpp"
#include "parallel/mpi_utils.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>

namespace fvm3d::parallel {

MPIDomainDecomposer::MPIDomainDecomposer(
    int nx_global, int ny_global, int nz_global,
    int np_x, int np_y, int np_z
)
    : nx_global_(nx_global), ny_global_(ny_global), nz_global_(nz_global) {

    // Auto-determine process grid if not specified
    if (np_x == 0 || np_y == 0 || np_z == 0) {
        int total_procs = MPIUtils::size();
        int cube_root = static_cast<int>(std::cbrt(total_procs) + 0.5);

        // Find factorization close to cubic
        if (np_x == 0) np_x = cube_root;
        if (np_y == 0) np_y = cube_root;
        if (np_z == 0) np_z = (total_procs + np_x * np_y - 1) / (np_x * np_y);

        // Adjust to match total process count
        while (np_x * np_y * np_z > total_procs) {
            if (np_z > 1) np_z--;
            else if (np_y > 1) np_y--;
            else np_x--;
        }
    }

    np_x_ = np_x;
    np_y_ = np_y;
    np_z_ = np_z;

    if (np_x_ * np_y_ * np_z_ != MPIUtils::size()) {
        throw std::invalid_argument(
            "Process grid dimensions don't match total process count"
        );
    }

    // Create Cartesian topology
    create_cartesian_topology();

    // Compute local domain sizes
    nx_local_ = compute_local_size(nx_global, np_x_, my_x_);
    ny_local_ = compute_local_size(ny_global, np_y_, my_y_);
    nz_local_ = compute_local_size(nz_global, np_z_, my_z_);
}

int MPIDomainDecomposer::compute_local_size(
    int global_size, int num_procs, int my_rank
) {
    int base_size = global_size / num_procs;
    int remainder = global_size % num_procs;

    // First 'remainder' processes get one extra cell
    if (my_rank < remainder) {
        return base_size + 1;
    } else {
        return base_size;
    }
}

void MPIDomainDecomposer::create_cartesian_topology() {
    // Define dimensions (periodic boundaries)
    int dims[3] = {np_x_, np_y_, np_z_};
    int periods[3] = {0, 0, 0};  // Non-periodic in all directions
    int reorder = 1;  // Allow MPI to reorder processes for optimization

    // Create Cartesian communicator
    int error = MPI_Cart_create(
        MPI_COMM_WORLD,
        3,           // 3D Cartesian
        dims,
        periods,
        reorder,
        &cartesian_comm_
    );
    MPIUtils::check_mpi_error(error, "MPI_Cart_create");

    // Get my coordinates in the Cartesian grid
    int my_coords[3];
    error = MPI_Cart_coords(cartesian_comm_, MPIUtils::rank(), 3, my_coords);
    MPIUtils::check_mpi_error(error, "MPI_Cart_coords");

    my_x_ = my_coords[0];
    my_y_ = my_coords[1];
    my_z_ = my_coords[2];
}

MPIDomainDecomposer::IndexRange MPIDomainDecomposer::global_index_range() const {
    int i_min = 0, j_min = 0, k_min = 0;

    // Compute start indices
    for (int i = 0; i < my_x_; i++) {
        i_min += compute_local_size(nx_global_, np_x_, i);
    }
    for (int j = 0; j < my_y_; j++) {
        j_min += compute_local_size(ny_global_, np_y_, j);
    }
    for (int k = 0; k < my_z_; k++) {
        k_min += compute_local_size(nz_global_, np_z_, k);
    }

    return {
        i_min, i_min + nx_local_,
        j_min, j_min + ny_local_,
        k_min, k_min + nz_local_
    };
}

MPIDomainDecomposer::Neighbors MPIDomainDecomposer::neighbors() const {
    Neighbors nb;
    nb.x_minus = neighbor_rank(0, -1);
    nb.x_plus = neighbor_rank(0, +1);
    nb.y_minus = neighbor_rank(1, -1);
    nb.y_plus = neighbor_rank(1, +1);
    nb.z_minus = neighbor_rank(2, -1);
    nb.z_plus = neighbor_rank(2, +1);
    return nb;
}

bool MPIDomainDecomposer::has_neighbor(int direction, int side) const {
    return neighbor_rank(direction, side) != MPI_PROC_NULL;
}

int MPIDomainDecomposer::neighbor_rank(int direction, int side) const {
    int coords[3] = {my_x_, my_y_, my_z_};
    int disp = (side > 0) ? 1 : -1;
    coords[direction] += disp;

    // Check boundary
    int np_dims[3] = {np_x_, np_y_, np_z_};
    if (coords[direction] < 0 || coords[direction] >= np_dims[direction]) {
        return MPI_PROC_NULL;
    }

    int rank;
    int error = MPI_Cart_rank(cartesian_comm_, coords, &rank);
    MPIUtils::check_mpi_error(error, "MPI_Cart_rank");

    return rank;
}

bool MPIDomainDecomposer::is_boundary(int direction) const {
    int coords[3] = {my_x_, my_y_, my_z_};
    int np_dims[3] = {np_x_, np_y_, np_z_};
    int coord = coords[direction];
    int np = np_dims[direction];
    return (coord == 0) || (coord == np - 1);
}

void MPIDomainDecomposer::print_info() const {
    if (!MPIUtils::is_root()) {
        return;
    }

    std::cout << "\n=== MPI Domain Decomposition ===" << std::endl;
    std::cout << "Global domain: " << nx_global_ << " x "
             << ny_global_ << " x " << nz_global_ << std::endl;
    std::cout << "Process grid: " << np_x_ << " x "
             << np_y_ << " x " << np_z_ << std::endl;
    std::cout << "Total processes: " << np_x_ * np_y_ * np_z_ << std::endl;

    // Print per-process info (only root, but show all)
    std::cout << "\nLocal domain sizes:" << std::endl;
    std::cout << "  Process (0,0,0): " << compute_local_size(nx_global_, np_x_, 0) << " x "
             << compute_local_size(ny_global_, np_y_, 0) << " x "
             << compute_local_size(nz_global_, np_z_, 0) << std::endl;
}

} // namespace fvm3d::parallel
