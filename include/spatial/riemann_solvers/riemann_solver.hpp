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
     * Get maximum wave speed for CFL condition.
     */
    virtual double max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const = 0;

    /**
     * Get the name of the Riemann solver.
     */
    virtual std::string name() const = 0;

    virtual int num_variables() const { return physics_->num_variables(); }

protected:
    std::shared_ptr<physics::PhysicsBase> physics_;  ///< Physics object for constants and methods

    /**
     * Extract primitive variables from conservative variables.
     * Uses physics object's conversion method.
     */
    void conservative_to_primitive(
        const Eigen::VectorXd& U,
        double& rho, double& u, double& v, double& w, double& p
    ) const {
        Eigen::VectorXd V = physics_->conservative_to_primitive(U);
        rho = V(0);
        u = V(1);
        v = V(2);
        w = V(3);
        p = V(4);
    }

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

    /**
     * Extract velocity component in a given direction.
     */
    double velocity_in_direction(const Eigen::VectorXd& U, int direction) const {
        double rho = std::max(U(0), physics_->rho_floor());
        return U(direction + 1) / rho;
    }
};

} // namespace fvm3d::spatial
