#include "spatial/reconstruction/muscl_reconstruction.hpp"
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace fvm3d::spatial {

MUSCLReconstruction::MUSCLReconstruction(
    int num_vars,
    const std::string& limiter,
    double kappa
) : LimitedReconstructionMethod("muscl", 2, limiter),
    num_vars_(num_vars),
    kappa_(kappa) {

    if (num_vars <= 0) {
        throw std::invalid_argument("MUSCLReconstruction: num_vars must be positive");
    }

    if (kappa < -1.0 || kappa > 1.0) {
        throw std::invalid_argument("MUSCLReconstruction: kappa must be in [-1, 1]");
    }
}

void MUSCLReconstruction::reconstruct(
    const core::Field3D<double>& U,
    int i, int j, int k,
    int direction,
    Eigen::VectorXd& U_L,
    Eigen::VectorXd& U_R
) const {
    U_L.resize(num_vars_);
    U_R.resize(num_vars_);

    // Reconstruct each variable independently
    // Standard MUSCL uses 3-point stencil for each side
    for (int var = 0; var < num_vars_; ++var) {
        // Get cell values based on direction
        double U_minus1, U_center, U_plus1, U_plus2;

        if (direction == 0) {  // x-direction
            U_minus1 = get_cell_value(U, var, i - 1, j, k);
            U_center = get_cell_value(U, var, i, j, k);
            U_plus1 = get_cell_value(U, var, i + 1, j, k);
            U_plus2 = get_cell_value(U, var, i + 2, j, k);
        } else if (direction == 1) {  // y-direction
            U_minus1 = get_cell_value(U, var, i, j - 1, k);
            U_center = get_cell_value(U, var, i, j, k);
            U_plus1 = get_cell_value(U, var, i, j + 1, k);
            U_plus2 = get_cell_value(U, var, i, j + 2, k);
        } else {  // z-direction
            U_minus1 = get_cell_value(U, var, i, j, k - 1);
            U_center = get_cell_value(U, var, i, j, k);
            U_plus1 = get_cell_value(U, var, i, j, k + 1);
            U_plus2 = get_cell_value(U, var, i, j, k + 2);
        }

        // Compute limited slopes for left and right cells
        double slope_L = compute_limited_slope(U_minus1, U_center, U_plus1);
        double slope_R = compute_limited_slope(U_center, U_plus1, U_plus2);

        // MUSCL reconstruction at interface i+1/2
        // U_L: extrapolate from cell i (using slope_L)
        U_L(var) = U_center + 0.5 * slope_L;

        // U_R: extrapolate from cell i+1 (using slope_R)
        U_R(var) = U_plus1 - 0.5 * slope_R;
    }
}

void MUSCLReconstruction::reconstruct_all(
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

    // Loop over all interfaces (excluding ghost cells that don't have enough stencil)
    int i_start = (direction == 0) ? 2 : 0;
    int j_start = (direction == 1) ? 2 : 0;
    int k_start = (direction == 2) ? 2 : 0;

    int i_end = (direction == 0) ? nx - 1 : n_interfaces_x;
    int j_end = (direction == 1) ? ny - 1 : n_interfaces_y;
    int k_end = (direction == 2) ? nz - 1 : n_interfaces_z;

    for (int i = i_start; i < i_end; ++i) {
        for (int j = j_start; j < j_end; ++j) {
            for (int k = k_start; k < k_end; ++k) {
                Eigen::VectorXd U_L, U_R;

                // Get cell index (interface at i+1/2 uses cells i and i+1)
                int i_cell = (direction == 0) ? i : i;
                int j_cell = (direction == 1) ? j : j;
                int k_cell = (direction == 2) ? k : k;

                reconstruct(U, i_cell, j_cell, k_cell, direction, U_L, U_R);

                // Store in output fields
                for (int var = 0; var < num_vars_; ++var) {
                    U_left(var, i, j, k) = U_L(var);
                    U_right(var, i, j, k) = U_R(var);
                }
            }
        }
    }
}

double MUSCLReconstruction::compute_limited_slope(
    double U_minus,
    double U_center,
    double U_plus
) const {
    // Compute consecutive differences
    double delta_minus = U_center - U_minus;
    double delta_plus = U_plus - U_center;

    // Avoid division by zero
    const double epsilon = 1.0e-14;
    if (std::abs(delta_plus) < epsilon) {
        return 0.0;
    }

    // Compute slope ratio
    double r = delta_minus / delta_plus;

    // Apply limiter
    double phi = apply_limiter(r);

    // Return limited slope
    return phi * delta_plus;
}

double MUSCLReconstruction::get_cell_value(
    const core::Field3D<double>& U,
    int var, int i, int j, int k
) const {
    // Bounds checking
    if (i < 0 || i >= U.nx() || j < 0 || j >= U.ny() || k < 0 || k >= U.nz()) {
        throw std::out_of_range("MUSCLReconstruction: cell index out of bounds");
    }
    return U(var, i, j, k);
}

} // namespace fvm3d::spatial
