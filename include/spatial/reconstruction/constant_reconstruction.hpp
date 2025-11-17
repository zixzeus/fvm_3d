#pragma once

#include "reconstruction_base.hpp"

namespace fvm3d::spatial {

/**
 * First-order constant (piecewise constant) reconstruction.
 *
 * This is the simplest reconstruction method where the left and right
 * interface states are simply set to the cell-centered value.
 *
 * U_L(i+1/2) = U(i)
 * U_R(i+1/2) = U(i+1)
 *
 * Properties:
 * - Order: 1
 * - TVD: Yes (no oscillations possible)
 * - Limiters: None (not needed)
 * - Ghost cells: 1 (standard for finite volume)
 *
 * Advantages:
 * - Simple and robust
 * - Always TVD
 * - Minimal computational cost
 *
 * Disadvantages:
 * - Only first-order accurate
 * - Excessive numerical diffusion for smooth flows
 */
class ConstantReconstruction : public ReconstructionMethod {
public:
    /**
     * Constructor for constant reconstruction.
     * @param num_vars: Number of conserved variables
     */
    explicit ConstantReconstruction(int num_vars);

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
};

} // namespace fvm3d::spatial
