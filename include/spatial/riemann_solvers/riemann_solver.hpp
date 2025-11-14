#pragma once

#include <Eigen/Dense>
#include <memory>
#include <string>

namespace fvm3d::spatial {

/**
 * Abstract base class for Riemann solvers.
 */
class RiemannSolver {
public:
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

    virtual int num_variables() const { return 5; }

protected:
    static constexpr double GAMMA = 1.4;
    static constexpr double RHO_FLOOR = 1e-10;
    static constexpr double P_FLOOR = 1e-11;

    /**
     * Extract primitive variables from conservative variables.
     */
    void conservative_to_primitive(
        const Eigen::VectorXd& U,
        double& rho, double& u, double& v, double& w, double& p
    ) const {
        rho = std::max(U(0), RHO_FLOOR);
        u = U(1) / rho;
        v = U(2) / rho;
        w = U(3) / rho;

        double kinetic_energy = 0.5 * (U(1)*U(1) + U(2)*U(2) + U(3)*U(3)) / rho;
        double internal_energy = U(4) / rho - kinetic_energy;
        p = std::max((GAMMA - 1.0) * rho * internal_energy, P_FLOOR);
    }

    /**
     * Compute speed of sound.
     */
    double sound_speed(double rho, double p) const {
        rho = std::max(rho, RHO_FLOOR);
        p = std::max(p, P_FLOOR);
        return std::sqrt(GAMMA * p / rho);
    }

    /**
     * Extract velocity component in a given direction.
     */
    double velocity_in_direction(const Eigen::VectorXd& U, int direction) const {
        double rho = std::max(U(0), RHO_FLOOR);
        return U(direction + 1) / rho;
    }
};

} // namespace fvm3d::spatial
