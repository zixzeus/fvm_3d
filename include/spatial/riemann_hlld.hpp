#pragma once

#include "spatial/riemann_solver.hpp"
#include "physics/resistive_mhd3d.hpp"
#include <memory>

namespace fvm3d::spatial {

/**
 * HLLD (Harten-Lax-van Leer-Discontinuities) Riemann Solver for MHD.
 *
 * The HLLD solver is specifically designed for magnetohydrodynamics equations.
 * It properly resolves:
 * - Fast magnetosonic shocks
 * - Slow magnetosonic waves
 * - Alfvén waves
 * - Contact discontinuities
 *
 * The solver uses a 7-wave approximate Riemann solution:
 * U_L | S_fL | U_fL* | S_aL | U_aL* | S_c | U_aR* | S_aR | U_fR* | S_fR | U_R
 *
 * Where:
 *   S_f  : Fast magnetosonic wave speeds
 *   S_a  : Alfvén wave speeds
 *   S_c  : Contact discontinuity speed
 *   *    : Intermediate states between waves
 */
class HLLDSolver : public RiemannSolver {
public:
    /**
     * Constructor with MHD physics.
     */
    explicit HLLDSolver(const std::shared_ptr<physics::ResistiveMHD3D>& mhd);

    /**
     * Solve 1D MHD Riemann problem.
     * @param U_L: Left state
     * @param U_R: Right state
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
    std::shared_ptr<physics::ResistiveMHD3D> mhd_;

    /**
     * Estimate fastest wave speed using Roe average.
     */
    struct WaveSpeeds {
        double S_fL;   // Left fast magnetosonic
        double S_aL;   // Left Alfvén
        double S_c;    // Contact discontinuity
        double S_aR;   // Right Alfvén
        double S_fR;   // Right fast magnetosonic
    };

    /**
     * Compute wave speed estimates.
     */
    WaveSpeeds estimate_wave_speeds(
        double rho_L, double u_L, double v_L, double w_L, double p_L,
        double Bx, double By_L, double Bz_L,
        double rho_R, double u_R, double v_R, double w_R, double p_R,
        double By_R, double Bz_R
    ) const;

    /**
     * Compute intermediate states in the HLLD structure.
     * Returns intermediate states for state reconstruction.
     */
    struct HLLDStates {
        Eigen::VectorXd U_L;
        Eigen::VectorXd U_fL_star;
        Eigen::VectorXd U_aL_star;
        Eigen::VectorXd U_aR_star;
        Eigen::VectorXd U_fR_star;
        Eigen::VectorXd U_R;
    };

    /**
     * Compute HLLD intermediate states.
     */
    HLLDStates compute_hlld_states(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const WaveSpeeds& speeds,
        int direction
    ) const;

    /**
     * Rotate state to align with normal direction.
     * This allows solving 1D Riemann problem in normal direction.
     */
    Eigen::VectorXd rotate_to_normal(const Eigen::VectorXd& U, int direction) const;

    /**
     * Rotate state back from normal direction.
     */
    Eigen::VectorXd rotate_from_normal(const Eigen::VectorXd& U, int direction) const;
};

} // namespace fvm3d::spatial
