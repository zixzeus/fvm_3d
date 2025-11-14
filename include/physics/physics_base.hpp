#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <cmath>

namespace fvm3d::physics {

/**
 * Abstract base class for physics equations.
 *
 * Defines the standard interface that all physics equation implementations
 * must follow for integration with the FVM framework.
 *
 * All physics implementations should inherit from this class and implement
 * the pure virtual methods for:
 * - Variable conversion (conservative ↔ primitive)
 * - Flux computation (in all spatial directions)
 * - Wave speed calculation (for CFL condition)
 */
class PhysicsBase {
public:
    /**
     * Constructor.
     * @param name: Human-readable name of the physics equations
     * @param num_variables: Number of conservative variables
     */
    PhysicsBase(const std::string& name, int num_variables)
        : name_(name), num_variables_(num_variables) {}

    virtual ~PhysicsBase() = default;

    // Prevent copying (physics objects should be managed by pointers/references)
    PhysicsBase(const PhysicsBase&) = delete;
    PhysicsBase& operator=(const PhysicsBase&) = delete;

    // Allow moving
    PhysicsBase(PhysicsBase&&) = default;
    PhysicsBase& operator=(PhysicsBase&&) = default;

    // ========== Core Physics Methods (Pure Virtual) ==========

    /**
     * Compute flux in a given direction.
     *
     * This is the primary flux computation method. Implementations should
     * compute the physical flux F(U) in the specified direction.
     *
     * @param U: Conservative variable vector
     * @param direction: Spatial direction (0=x, 1=y, 2=z)
     * @return Flux vector F(U)
     */
    virtual Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U,
        int direction
    ) const = 0;

    /**
     * Compute maximum characteristic wave speed.
     *
     * Used for CFL time step calculation. Should return the maximum
     * eigenvalue magnitude in the specified direction.
     *
     * @param U: Conservative variable vector
     * @param direction: Spatial direction (0=x, 1=y, 2=z)
     * @return Maximum wave speed (always positive)
     */
    virtual double max_wave_speed(
        const Eigen::VectorXd& U,
        int direction
    ) const = 0;

    // ========== Optional Methods (With Default Implementations) ==========

    /**
     * Compute source terms S(U).
     *
     * For conservation laws of the form: ∂U/∂t + ∇·F(U) = S(U)
     * Default implementation returns zero (no source terms).
     *
     * @param U: Conservative variable vector
     * @param x, y, z: Spatial coordinates (for position-dependent sources)
     * @param dx, dy, dz: Grid spacing (for gradient-based sources)
     * @return Source term vector S(U)
     */
    virtual Eigen::VectorXd compute_source(
        const Eigen::VectorXd& U,
        double x, double y, double z,
        double dx, double dy, double dz
    ) const {
        return Eigen::VectorXd::Zero(num_variables_);
    }

    /**
     * Validate physical consistency of state.
     *
     * Checks if the state variables are physically valid (e.g., positive
     * density, positive pressure for gases, etc.).
     * Default implementation only checks for NaN/Inf.
     *
     * @param U: Conservative variable vector
     * @return true if state is physically valid
     */
    virtual bool validate_state(const Eigen::VectorXd& U) const {
        if (U.size() != num_variables_) {
            return false;
        }

        for (int i = 0; i < U.size(); ++i) {
            if (std::isnan(U(i)) || std::isinf(U(i))) {
                return false;
            }
        }

        return true;
    }

    /**
     * Compute characteristic eigenvalues (wave speeds).
     *
     * Returns all eigenvalues of the flux Jacobian ∂F/∂U in the given direction.
     * Default implementation returns max_wave_speed repeated num_variables times.
     *
     * @param U: Conservative variable vector
     * @param direction: Spatial direction (0=x, 1=y, 2=z)
     * @return Vector of eigenvalues
     */
    virtual Eigen::VectorXd compute_eigenvalues(
        const Eigen::VectorXd& U,
        int direction
    ) const {
        double max_speed = max_wave_speed(U, direction);
        return Eigen::VectorXd::Constant(num_variables_, max_speed);
    }

    // ========== Getters ==========

    /**
     * Get human-readable name of physics equations.
     */
    std::string name() const { return name_; }

    /**
     * Get number of conservative variables.
     */
    int num_variables() const { return num_variables_; }

    /**
     * Get names of conservative variables.
     * Override in derived classes for physics-specific names.
     *
     * @return Vector of variable names (e.g., ["rho", "rho_u", "rho_v", "E"])
     */
    virtual std::vector<std::string> get_variable_names() const {
        std::vector<std::string> names;
        names.reserve(num_variables_);
        for (int i = 0; i < num_variables_; ++i) {
            names.push_back("var_" + std::to_string(i));
        }
        return names;
    }

    /**
     * Get names of primitive variables.
     * Override in derived classes for physics-specific names.
     *
     * @return Vector of primitive variable names (e.g., ["rho", "u", "v", "p"])
     */
    virtual std::vector<std::string> get_primitive_names() const {
        std::vector<std::string> names;
        names.reserve(num_variables_);
        for (int i = 0; i < num_variables_; ++i) {
            names.push_back("prim_" + std::to_string(i));
        }
        return names;
    }

    /**
     * Check if this physics has temperature as a variable.
     * Override to return true if temperature can be computed.
     */
    virtual bool has_temperature() const {
        return false;
    }

    /**
     * Check if this physics has pressure as a variable.
     * Override to return true if pressure can be computed.
     */
    virtual bool has_pressure() const {
        return false;
    }

protected:
    std::string name_;         ///< Human-readable name
    int num_variables_;        ///< Number of conservative variables
};

/**
 * Base class for hyperbolic conservation laws.
 *
 * Specialized base class for physics that follow the conservation law form:
 *   ∂U/∂t + ∇·F(U) = S(U)
 *
 * Provides additional methods for flux Jacobian computation and
 * hyperbolicity checking.
 */
class ConservationLaw : public PhysicsBase {
public:
    /**
     * Constructor.
     * @param name: Name of the conservation law
     * @param num_variables: Number of conservative variables
     * @param num_dimensions: Spatial dimensions (default: 3)
     */
    ConservationLaw(const std::string& name, int num_variables, int num_dimensions = 3)
        : PhysicsBase(name, num_variables), num_dimensions_(num_dimensions) {}

    virtual ~ConservationLaw() = default;

    /**
     * Get number of spatial dimensions.
     */
    int num_dimensions() const { return num_dimensions_; }

    /**
     * Compute flux Jacobian matrix ∂F/∂U.
     *
     * Uses finite differences for approximation by default.
     * Override for analytical Jacobian computation.
     *
     * @param U: Conservative variable vector
     * @param direction: Spatial direction
     * @return Jacobian matrix (num_vars × num_vars)
     */
    virtual Eigen::MatrixXd compute_flux_jacobian(
        const Eigen::VectorXd& U,
        int direction
    ) const {
        constexpr double eps = 1e-8;
        Eigen::MatrixXd jacobian(num_variables_, num_variables_);

        Eigen::VectorXd f0 = compute_flux(U, direction);

        for (int i = 0; i < num_variables_; ++i) {
            Eigen::VectorXd U_pert = U;
            U_pert(i) += eps;
            Eigen::VectorXd f_pert = compute_flux(U_pert, direction);
            jacobian.col(i) = (f_pert - f0) / eps;
        }

        return jacobian;
    }

    /**
     * Check if system is hyperbolic at given state.
     *
     * A system is hyperbolic if all eigenvalues of the flux Jacobian are real.
     *
     * @param U: Conservative variable vector
     * @param direction: Spatial direction
     * @return true if system has real eigenvalues
     */
    virtual bool is_hyperbolic(
        const Eigen::VectorXd& U,
        int direction
    ) const {
        try {
            Eigen::MatrixXd jacobian = compute_flux_jacobian(U, direction);
            Eigen::VectorXcd eigenvals = jacobian.eigenvalues();

            // Check if all eigenvalues are real (imaginary part near zero)
            constexpr double tol = 1e-10;
            for (int i = 0; i < eigenvals.size(); ++i) {
                if (std::abs(eigenvals(i).imag()) > tol) {
                    return false;
                }
            }
            return true;
        } catch (...) {
            return false;
        }
    }

protected:
    int num_dimensions_;  ///< Number of spatial dimensions
};

} // namespace fvm3d::physics
