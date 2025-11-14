#pragma once

#include "reconstruction_base.hpp"

namespace fvm3d::spatial {

/**
 * MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws) reconstruction.
 *
 * Second-order accurate reconstruction using limited linear interpolation.
 * The interface states are reconstructed using:
 *
 * U_L(i+1/2) = U(i) + 0.5 * phi(r_i) * (U(i) - U(i-1))
 * U_R(i+1/2) = U(i+1) - 0.5 * phi(r_{i+1}) * (U(i+2) - U(i+1))
 *
 * where phi(r) is the slope limiter function and r is the ratio of consecutive slopes.
 *
 * Properties:
 * - Order: 2 (in smooth regions)
 * - TVD: Yes (with appropriate limiter)
 * - Limiters: minmod, van Leer, superbee, MC
 * - Ghost cells: 2
 *
 * Advantages:
 * - Second-order accurate in smooth regions
 * - TVD property prevents spurious oscillations
 * - Widely used and well-tested
 *
 * Disadvantages:
 * - Reduces to first-order at extrema
 * - More computational cost than constant reconstruction
 */
class MUSCLReconstruction : public LimitedReconstructionMethod {
public:
    /**
     * Constructor for MUSCL reconstruction.
     * @param num_vars: Number of conserved variables
     * @param limiter: Slope limiter ("minmod", "van_leer", "superbee", "mc")
     * @param kappa: MUSCL parameter (-1 to 1, default 1/3 for 3rd order at smooth extrema)
     */
    explicit MUSCLReconstruction(
        int num_vars,
        const std::string& limiter = "minmod",
        double kappa = 1.0/3.0
    );

    void reconstruct(
        const core::Field3D<double>& U,
        int i, int j, int k,
        int direction,
        Eigen::VectorXd& U_L,
        Eigen::VectorXd& U_R
    ) const override;

    void reconstruct_all(
        const core::Field3D<double>& U,
        int direction,
        core::Field3D<double>& U_left,
        core::Field3D<double>& U_right
    ) const override;

private:
    int num_vars_;
    double kappa_;  // MUSCL parameter for interpolation

    /**
     * Compute limited slope at a cell.
     * @param U_minus: Value at i-1 (or j-1, k-1)
     * @param U_center: Value at i (or j, k)
     * @param U_plus: Value at i+1 (or j+1, k+1)
     * @return Limited slope
     */
    double compute_limited_slope(double U_minus, double U_center, double U_plus) const;

    /**
     * Extract cell value for a specific variable and location.
     */
    double get_cell_value(
        const core::Field3D<double>& U,
        int var, int i, int j, int k
    ) const;
};

} // namespace fvm3d::spatial
