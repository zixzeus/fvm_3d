#pragma once

#include "core/field3d.hpp"
#include <Eigen/Dense>
#include <string>

namespace fvm3d::spatial {

/**
 * Abstract base class for spatial reconstruction methods.
 *
 * Reconstruction methods compute left and right interface states from
 * cell-centered conservative variables. These interface states are then
 * used by flux calculators to compute numerical fluxes.
 *
 * Design Philosophy:
 * - Unified interface for all reconstruction methods
 * - Supports 1st order (constant), 2nd order (linear), and higher order schemes
 * - Includes slope limiters for TVD/ENO/WENO properties
 *
 * Common Methods:
 * - Constant (1st order): Piecewise constant reconstruction
 * - MUSCL (2nd order): Monotone Upstream-centered Scheme for Conservation Laws
 * - WENO (high order): Weighted Essentially Non-Oscillatory
 * - ENO (high order): Essentially Non-Oscillatory
 */
class ReconstructionMethod {
public:
    /**
     * Constructor.
     * @param name: Name of the reconstruction method
     * @param order: Order of accuracy (1, 2, 3, 5, etc.)
     */
    explicit ReconstructionMethod(const std::string& name, int order)
        : name_(name), order_(order) {}

    virtual ~ReconstructionMethod() = default;

    /**
     * Reconstruct left and right states at interface i+1/2.
     *
     * Given cell-centered values at i-1, i, i+1, i+2, etc.,
     * reconstruct left state U_L(i+1/2) and right state U_R(i+1/2).
     *
     * @param U: Cell-centered conservative variables
     * @param i: Cell index
     * @param j: Cell index (y-direction)
     * @param k: Cell index (z-direction)
     * @param direction: 0=x, 1=y, 2=z
     * @param U_L: Output left state
     * @param U_R: Output right state
     */
    virtual void reconstruct(
        const core::Field3D<double>& U,
        int i, int j, int k,
        int direction,
        Eigen::VectorXd& U_L,
        Eigen::VectorXd& U_R
    ) const = 0;

    /**
     * Reconstruct all interface states in given direction.
     *
     * More efficient than calling reconstruct() in a loop.
     *
     * @param U: Cell-centered conservative variables (nvars × nx × ny × nz)
     * @param direction: 0=x, 1=y, 2=z
     * @param U_left: Output left states at all interfaces
     * @param U_right: Output right states at all interfaces
     */
    virtual void reconstruct_all(
        const core::Field3D<double>& U,
        int direction,
        core::Field3D<double>& U_left,
        core::Field3D<double>& U_right
    ) const;

    /**
     * Get name of this reconstruction method.
     */
    std::string name() const { return name_; }

    /**
     * Get order of accuracy.
     */
    int order() const { return order_; }

    /**
     * Check if this method is TVD (Total Variation Diminishing).
     *
     * TVD methods prevent spurious oscillations near discontinuities.
     *
     * @return: True if TVD
     */
    virtual bool is_tvd() const { return order_ <= 1; }

    /**
     * Check if this method uses slope limiters.
     *
     * @return: True if limiters are used
     */
    virtual bool uses_limiters() const { return false; }

    /**
     * Get number of ghost cells required.
     *
     * Higher-order methods need more ghost cells for the stencil.
     *
     * @return: Number of ghost cells needed
     */
    virtual int required_ghost_cells() const {
        return (order_ + 1) / 2;  // Rule of thumb
    }

protected:
    std::string name_;  ///< Name of the method
    int order_;         ///< Order of accuracy
};

/**
 * Base class for reconstruction methods with slope limiters.
 *
 * Slope limiters ensure TVD property by limiting the slope of
 * the reconstruction near discontinuities.
 */
class LimitedReconstructionMethod : public ReconstructionMethod {
public:
    /**
     * Constructor.
     * @param name: Name of the reconstruction method
     * @param order: Order of accuracy
     * @param limiter_name: Name of the slope limiter
     */
    explicit LimitedReconstructionMethod(
        const std::string& name,
        int order,
        const std::string& limiter_name
    ) : ReconstructionMethod(name, order), limiter_name_(limiter_name) {}

    /**
     * This method uses limiters.
     */
    bool uses_limiters() const override { return true; }

    /**
     * Get name of the slope limiter.
     */
    std::string limiter_name() const { return limiter_name_; }

    /**
     * TVD methods with limiters.
     */
    bool is_tvd() const override { return true; }

protected:
    /**
     * Apply slope limiter function.
     *
     * Common limiters: minmod, van Leer, superbee, MC
     *
     * @param r: Ratio of consecutive slopes
     * @return: Limited slope ratio
     */
    double apply_limiter(double r) const;

    std::string limiter_name_;  ///< Name of the slope limiter
};

} // namespace fvm3d::spatial
