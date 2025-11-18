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

    /**
     * Estimate contact discontinuity speed using Roe averaging.
     */
    double estimate_s_contact(
        double rho_L, double u_L, double p_L, double S_L,
        double rho_R, double u_R, double p_R, double S_R
    ) const;

    /**
     * Estimate middle pressure using Roe-type averaging.
     */
    double estimate_p_middle(
        const Eigen::VectorXd& V_L, const Eigen::VectorXd& V_R
    ) const;

    /**
     * Compute left intermediate state U_*L using full HLLD formula.
     * Includes proper Alfvén wave handling.
     */
    Eigen::VectorXd compute_state_L(
        const Eigen::VectorXd& U_L, const Eigen::VectorXd& V_L,
        const Eigen::VectorXd& U_hll, const Eigen::VectorXd& V_hll,
        double S_L, double S_M, double p_m, double B_x,
        double rho_Lstar
    ) const;

    /**
     * Compute right intermediate state U_*R using full HLLD formula.
     * Includes proper Alfvén wave handling.
     */
    Eigen::VectorXd compute_state_R(
        const Eigen::VectorXd& U_R, const Eigen::VectorXd& V_R,
        const Eigen::VectorXd& U_hll, const Eigen::VectorXd& V_hll,
        double S_R, double S_M, double p_m, double B_x,
        double rho_Rstar
    ) const;

    /**
     * Compute central region state U_** (between Alfvén waves).
     * Uses Roe averaging for tangential components.
     */
    Eigen::VectorXd compute_central_state(
        const Eigen::VectorXd& U_Lstar, const Eigen::VectorXd& U_Rstar,
        const Eigen::VectorXd& V_Lstar, const Eigen::VectorXd& V_Rstar,
        double S_M, double B_x, int direction
    ) const;

    /**
     * Compute Alfvén wave speeds aL1 (= S_M - c_A) and aR1 (= S_M + c_A).
     */
    void compute_alfven_speeds(
        double B_x, double rho_Lstar, double rho_Rstar,
        double& aL1, double& aR1
    ) const;
};

} // namespace fvm3d::spatial
