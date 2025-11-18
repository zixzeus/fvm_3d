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

void MUSCLReconstruction::set_reconstruction_config(
    const ReconstructionConfig& config,
    const std::shared_ptr<physics::PhysicsBase>& physics
) {
    config_ = config;
    physics_ = physics;
    use_mixed_reconstruction_ = true;

    // Validate configuration
    if (config_.var_types.size() != static_cast<size_t>(num_vars_)) {
        throw std::invalid_argument(
            "MUSCLReconstruction: config size mismatch with num_vars"
        );
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

    // If using mixed reconstruction, convert states to primitive first
    Eigen::VectorXd V_minus1, V_center, V_plus1, V_plus2;

    if (use_mixed_reconstruction_) {
        // Extract conservative states
        Eigen::VectorXd U_minus1_cons(num_vars_), U_center_cons(num_vars_);
        Eigen::VectorXd U_plus1_cons(num_vars_), U_plus2_cons(num_vars_);

        for (int var = 0; var < num_vars_; ++var) {
            if (direction == 0) {
                U_minus1_cons(var) = get_cell_value(U, var, i - 1, j, k);
                U_center_cons(var) = get_cell_value(U, var, i, j, k);
                U_plus1_cons(var) = get_cell_value(U, var, i + 1, j, k);
                U_plus2_cons(var) = get_cell_value(U, var, i + 2, j, k);
            } else if (direction == 1) {
                U_minus1_cons(var) = get_cell_value(U, var, i, j - 1, k);
                U_center_cons(var) = get_cell_value(U, var, i, j, k);
                U_plus1_cons(var) = get_cell_value(U, var, i, j + 1, k);
                U_plus2_cons(var) = get_cell_value(U, var, i, j + 2, k);
            } else {
                U_minus1_cons(var) = get_cell_value(U, var, i, j, k - 1);
                U_center_cons(var) = get_cell_value(U, var, i, j, k);
                U_plus1_cons(var) = get_cell_value(U, var, i, j, k + 1);
                U_plus2_cons(var) = get_cell_value(U, var, i, j, k + 2);
            }
        }

        // Convert to primitive variables
        V_minus1 = physics_->conservative_to_primitive(U_minus1_cons);
        V_center = physics_->conservative_to_primitive(U_center_cons);
        V_plus1 = physics_->conservative_to_primitive(U_plus1_cons);
        V_plus2 = physics_->conservative_to_primitive(U_plus2_cons);

        // Reconstruct in primitive form
        Eigen::VectorXd V_L(num_vars_), V_R(num_vars_);

        for (int var = 0; var < num_vars_; ++var) {
            double val_minus1, val_center, val_plus1, val_plus2;

            // Extract value based on reconstruction type
            if (config_.var_types[var] == ReconstructionVariableType::PRIMITIVE) {
                // Use primitive variables
                val_minus1 = V_minus1(var);
                val_center = V_center(var);
                val_plus1 = V_plus1(var);
                val_plus2 = V_plus2(var);
            } else {
                // Use conservative variables
                if (direction == 0) {
                    val_minus1 = get_cell_value(U, var, i - 1, j, k);
                    val_center = get_cell_value(U, var, i, j, k);
                    val_plus1 = get_cell_value(U, var, i + 1, j, k);
                    val_plus2 = get_cell_value(U, var, i + 2, j, k);
                } else if (direction == 1) {
                    val_minus1 = get_cell_value(U, var, i, j - 1, k);
                    val_center = get_cell_value(U, var, i, j, k);
                    val_plus1 = get_cell_value(U, var, i, j + 1, k);
                    val_plus2 = get_cell_value(U, var, i, j + 2, k);
                } else {
                    val_minus1 = get_cell_value(U, var, i, j, k - 1);
                    val_center = get_cell_value(U, var, i, j, k);
                    val_plus1 = get_cell_value(U, var, i, j, k + 1);
                    val_plus2 = get_cell_value(U, var, i, j, k + 2);
                }
            }

            // Compute limited slopes
            double slope_L = compute_limited_slope(val_minus1, val_center, val_plus1);
            double slope_R = compute_limited_slope(val_center, val_plus1, val_plus2);

            // MUSCL reconstruction
            V_L(var) = val_center + 0.5 * slope_L;
            V_R(var) = val_plus1 - 0.5 * slope_R;
        }

        // ========== POSITIVITY-PRESERVING LIMITER (OpenMHD Style) ==========
        // After reconstruction, ensure density and pressure remain positive.
        // If not, reduce slopes to preserve positivity (fallback to first-order).
        // This prevents artificial mass/energy injection in primitive_to_conservative.
        // Reference: Waagan (2009) "A positive MUSCL-Hancock scheme for ideal MHD"

        const double rho_floor = physics_->rho_floor();
        const double p_floor = physics_->p_floor();

        // Check and fix density (index 0)
        if (V_L(0) < rho_floor || V_R(0) < rho_floor) {
            // Fallback to first-order (zero slope) for this cell
            V_L(0) = V_center(0);
            V_R(0) = V_plus1(0);
        }

        // Check and fix pressure (index 4 for MHD/Euler)
        if (num_vars_ >= 5 && (V_L(4) < p_floor || V_R(4) < p_floor)) {
            // Fallback to first-order (zero slope) for this cell
            V_L(4) = V_center(4);
            V_R(4) = V_plus1(4);
        }

        // Convert back to conservative form
        U_L = physics_->primitive_to_conservative(V_L);
        U_R = physics_->primitive_to_conservative(V_R);

    } else {
        // Original behavior: reconstruct conservative variables directly
        #pragma omp simd
        for (int var = 0; var < num_vars_; ++var) {
            double U_minus1, U_center, U_plus1, U_plus2;

            if (direction == 0) {
                U_minus1 = get_cell_value(U, var, i - 1, j, k);
                U_center = get_cell_value(U, var, i, j, k);
                U_plus1 = get_cell_value(U, var, i + 1, j, k);
                U_plus2 = get_cell_value(U, var, i + 2, j, k);
            } else if (direction == 1) {
                U_minus1 = get_cell_value(U, var, i, j - 1, k);
                U_center = get_cell_value(U, var, i, j, k);
                U_plus1 = get_cell_value(U, var, i, j + 1, k);
                U_plus2 = get_cell_value(U, var, i, j + 2, k);
            } else {
                U_minus1 = get_cell_value(U, var, i, j, k - 1);
                U_center = get_cell_value(U, var, i, j, k);
                U_plus1 = get_cell_value(U, var, i, j, k + 1);
                U_plus2 = get_cell_value(U, var, i, j, k + 2);
            }

            double slope_L = compute_limited_slope(U_minus1, U_center, U_plus1);
            double slope_R = compute_limited_slope(U_center, U_plus1, U_plus2);

            U_L(var) = U_center + 0.5 * slope_L;
            U_R(var) = U_plus1 - 0.5 * slope_R;
        }
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
