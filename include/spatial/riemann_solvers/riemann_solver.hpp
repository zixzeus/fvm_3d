#pragma once

#include <Eigen/Dense>
#include <memory>
#include <string>
#include "physics/physics_base.hpp"

namespace fvm3d::spatial {

/**
 * Abstract base class for Riemann solvers.
 */
class RiemannSolver {
public:
    /**
     * Constructor with physics object.
     * @param physics: Physics equation system for getting constants and methods
     */
    explicit RiemannSolver(const std::shared_ptr<physics::PhysicsBase>& physics)
        : physics_(physics) {}

    virtual ~RiemannSolver() = default;

    /**
     * Compute the numerical flux at an interface.
     * @param U_L: left conservative state
     * @param U_R: right conservative state
     * @param direction: 0=X, 1=Y, 2=Z
     * @return: numerical flux
     */
    virtual Eigen::VectorXd solve(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const = 0;

    /**
     * Get the name of the Riemann solver.
     */
    virtual std::string name() const = 0;

    virtual int num_variables() const { return physics_->num_variables(); }

    /**
     * Estimate left and right wave speeds using Davis estimates.
     *
     * CRITICAL FIX: Uses asymmetric wave speed estimates that account for
     * flow velocity asymmetry. Previous implementations assuming symmetric waves
     * (s_left = -s_right) caused excessive dissipation or instability.
     *
     * @param U_L: Left conservative state
     * @param U_R: Right conservative state
     * @param direction: Direction (0=x, 1=y, 2=z)
     * @return Pair of (left_wave_speed, right_wave_speed)
     */
    std::pair<double, double> estimate_wave_speeds(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const {
        // Convert to primitive variables
        Eigen::VectorXd V_L = physics_->conservative_to_primitive(U_L);
        Eigen::VectorXd V_R = physics_->conservative_to_primitive(U_R);

        // Extract velocity in normal direction
        double u_L = V_L(direction + 1);
        double u_R = V_R(direction + 1);

        // Get characteristic wave speeds using physics object
        // max_wave_speed returns |u_n| + c, so we extract c by subtracting |u_n|
        double max_speed_L = physics_->max_wave_speed(U_L, direction);
        double max_speed_R = physics_->max_wave_speed(U_R, direction);

        double c_L = max_speed_L - std::abs(u_L);
        double c_R = max_speed_R - std::abs(u_R);

        // Davis estimates (standard for Riemann solvers)
        double s_left = std::min(u_L - c_L, u_R - c_R);
        double s_right = std::max(u_L + c_L, u_R + c_R);

        return {s_left, s_right};
    }

protected:
    std::shared_ptr<physics::PhysicsBase> physics_;  ///< Physics object for constants and methods

    /**
     * Compute speed of sound.
     * Uses physics object's gamma constant.
     */
    double sound_speed(double rho, double p) const {
        double gamma = physics_->gamma();
        double rho_floor = physics_->rho_floor();
        double p_floor = physics_->p_floor();

        rho = std::max(rho, rho_floor);
        p = std::max(p, p_floor);
        return std::sqrt(gamma * p / rho);
    }
};

} // namespace fvm3d::spatial
