#include "spatial/reconstruction/constant_reconstruction.hpp"
#include <stdexcept>

namespace fvm3d::spatial {

ConstantReconstruction::ConstantReconstruction(int num_vars)
    : ReconstructionMethod("constant", 1), num_vars_(num_vars) {
    if (num_vars <= 0) {
        throw std::invalid_argument("ConstantReconstruction: num_vars must be positive");
    }
}

void ConstantReconstruction::reconstruct(
    const core::Field3D<double>& U,
    int i, int j, int k,
    int direction,
    Eigen::VectorXd& U_L,
    Eigen::VectorXd& U_R
) const {
    // For constant reconstruction:
    // Left state at interface (i+1/2, j, k) = cell-centered value at (i, j, k)
    // Right state at interface (i+1/2, j, k) = cell-centered value at (i+1, j, k)

    U_L.resize(num_vars_);
    U_R.resize(num_vars_);

    // Extract left state (current cell) with SIMD vectorization
    #pragma omp simd
    for (int var = 0; var < num_vars_; ++var) {
        U_L(var) = U(var, i, j, k);
    }

    // Extract right state (next cell in the specified direction)
    int i_next = i, j_next = j, k_next = k;
    if (direction == 0) {
        i_next = i + 1;
    } else if (direction == 1) {
        j_next = j + 1;
    } else if (direction == 2) {
        k_next = k + 1;
    }

    // SIMD vectorization for extracting right state
    #pragma omp simd
    for (int var = 0; var < num_vars_; ++var) {
        U_R(var) = U(var, i_next, j_next, k_next);
    }
}

void ConstantReconstruction::reconstruct_all(
    const core::Field3D<double>& U,
    int direction,
    core::Field3D<double>& U_left,
    core::Field3D<double>& U_right
) const {
    int nx = U.nx();
    int ny = U.ny();
    int nz = U.nz();

    // Determine number of interfaces in each direction
    int n_interfaces_x = (direction == 0) ? nx + 1 : nx;
    int n_interfaces_y = (direction == 1) ? ny + 1 : ny;
    int n_interfaces_z = (direction == 2) ? nz + 1 : nz;

    // Allocate output fields
    U_left = core::Field3D<double>(num_vars_, n_interfaces_x, n_interfaces_y, n_interfaces_z);
    U_right = core::Field3D<double>(num_vars_, n_interfaces_x, n_interfaces_y, n_interfaces_z);

    // Loop over all interfaces
    for (int i = 0; i < n_interfaces_x; ++i) {
        for (int j = 0; j < n_interfaces_y; ++j) {
            for (int k = 0; k < n_interfaces_z; ++k) {
                Eigen::VectorXd U_L, U_R;

                // Get cell indices for reconstruction
                int i_cell = (direction == 0 && i > 0) ? i - 1 : i;
                int j_cell = (direction == 1 && j > 0) ? j - 1 : j;
                int k_cell = (direction == 2 && k > 0) ? k - 1 : k;

                reconstruct(U, i_cell, j_cell, k_cell, direction, U_L, U_R);

                // Store in output fields with SIMD vectorization
                #pragma omp simd
                for (int var = 0; var < num_vars_; ++var) {
                    U_left(var, i, j, k) = U_L(var);
                    U_right(var, i, j, k) = U_R(var);
                }
            }
        }
    }
}

} // namespace fvm3d::spatial
