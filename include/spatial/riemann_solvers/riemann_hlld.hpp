#pragma once

#include "spatial/riemann_solvers/riemann_solver.hpp"
#include "physics/resistive_mhd3d_advanced.hpp"
#include <Eigen/Dense>
#include <memory>

namespace fvm3d::spatial {

/**
 * HLLD (Harten-Lax-van Leer-Discontinuities) Riemann Solver for MHD.
 *
 * Uses AdvancedResistiveMHD3D physics object for flux calculation.
 * Supports 9-variable GLM-MHD systems.
 *
 * Conservative variables:
 *   9-var: [ρ, ρu, ρv, ρw, E, Bx, By, Bz, ψ]  (with GLM divergence cleaning)
 */
class HLLDSolver : public RiemannSolver {
public:
    /**
     * Constructor with MHD physics object.
     * @param mhd: Advanced resistive MHD physics object
     */
    explicit HLLDSolver(const std::shared_ptr<physics::AdvancedResistiveMHD3D>& mhd);

    /**
     * Solve 1D MHD Riemann problem and return flux.
     * @param U_L: Left conservative state
     * @param U_R: Right conservative state
     * @param direction: 0=X, 1=Y, 2=Z
     * @return: Numerical flux at interface
     */
    Eigen::VectorXd solve(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override;

    /**
     * Get maximum wave speed.
     */
    double max_wave_speed(const Eigen::VectorXd& U_L, const Eigen::VectorXd& U_R, int direction) const override;

    /**
     * Get solver name.
     */
    std::string name() const override { return "HLLD"; }

private:
    std::shared_ptr<physics::AdvancedResistiveMHD3D> mhd_;  // MHD physics object

    /**
     * Rotate state to align with normal direction.
     */
    Eigen::VectorXd rotate_to_normal(const Eigen::VectorXd& U, int direction) const;

    /**
     * Rotate state back from normal direction.
     */
    Eigen::VectorXd rotate_from_normal(const Eigen::VectorXd& U, int direction) const;
};

} // namespace fvm3d::spatial
