#pragma once

#include "riemann_solver.hpp"

namespace fvm3d::spatial {

/**
 * HLLC (Harten-Lax-van Leer-Contact) Riemann solver.
 * Improves upon HLL by properly capturing contact discontinuities.
 * Uses three-wave model: left shock, contact discontinuity, right shock.
 *
 * More accurate than HLL for problems with strong contact discontinuities
 * (e.g., different densities with same velocity and pressure).
 */
class HLLCSolver : public RiemannSolver {
public:
    Eigen::VectorXd solve(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override;

    double max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override;

    std::string name() const override { return "HLLC"; }

private:
    /**
     * Estimate left and right wave speeds.
     */
    double estimate_s_left(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const;

    double estimate_s_right(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const;

    /**
     * Estimate contact discontinuity speed and middle pressure.
     * Uses Roe-type averaging formulas.
     */
    double estimate_s_contact(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction,
        double S_L,
        double S_R
    ) const;

    double estimate_p_middle(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const;

    Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U,
        int direction
    ) const;

    /**
     * Compute HLLC states in left and right regions.
     */
    Eigen::VectorXd compute_u_left_hllc(
        const Eigen::VectorXd& U_L,
        double S_L,
        double S_M,
        double p_m,
        int direction
    ) const;

    Eigen::VectorXd compute_u_right_hllc(
        const Eigen::VectorXd& U_R,
        double S_R,
        double S_M,
        double p_m,
        int direction
    ) const;
};

} // namespace fvm3d::spatial
