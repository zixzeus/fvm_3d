#pragma once

#include <Eigen/Dense>
#include <string>
#include <memory>

namespace fvm3d::spatial {

/**
 * Abstract base class for reconstruction schemes.
 * Reconstructs interface values from cell-averaged values.
 *
 * For a 1D problem with cell-averaged values U_i:
 * - Constant: U_L = U_i, U_R = U_i (0th order)
 * - MUSCL: U_L, U_R computed from neighboring cells (2nd order)
 * - WENO: Higher-order using larger stencils (3rd+ order)
 */
class ReconstructionScheme {
public:
    virtual ~ReconstructionScheme() = default;

    /**
     * Reconstruct left and right interface values.
     * @param U_L: value at left cell center
     * @param U_C: value at center cell center (reference cell)
     * @param U_R: value at right cell center
     * @param U_LL: value at left-left cell center (for higher-order methods)
     * @param U_RR: value at right-right cell center
     * @param U_recon_L: reconstructed value at left interface (output)
     * @param U_recon_R: reconstructed value at right interface (output)
     */
    virtual void reconstruct(
        double U_LL, double U_L, double U_C, double U_R, double U_RR,
        double& U_recon_L, double& U_recon_R
    ) const = 0;

    /**
     * Get the name of the reconstruction scheme.
     */
    virtual std::string name() const = 0;

    /**
     * Get the order of accuracy.
     */
    virtual int order() const = 0;
};

/**
 * Constant reconstruction (no interpolation).
 * Uses cell-averaged values directly: U_L = U_R = U_C
 * Order: 0 (piecewise constant)
 */
class ConstantReconstruction : public ReconstructionScheme {
public:
    void reconstruct(
        double U_LL, double U_L, double U_C, double U_R, double U_RR,
        double& U_recon_L, double& U_recon_R
    ) const override {
        U_recon_L = U_C;
        U_recon_R = U_C;
    }

    std::string name() const override { return "Constant"; }
    int order() const override { return 0; }
};

/**
 * MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws).
 * Reconstructs linear interface values with slope limiting.
 * Order: 2 (second-order accurate)
 *
 * Formula:
 *   U_L = U_C - 0.5 * slope
 *   U_R = U_C + 0.5 * slope
 * where slope is limited by minmod or van Leer limiter
 */
class MUSCLReconstruction : public ReconstructionScheme {
public:
    /**
     * Limiter type for slope computation.
     */
    enum class LimiterType {
        MINMOD,      // Most restrictive (most diffusive)
        VAN_LEER,    // Smooth transition
        SUPERBEE     // Least restrictive (least diffusive)
    };

    MUSCLReconstruction(LimiterType limiter = LimiterType::VAN_LEER)
        : limiter_(limiter) {}

    void reconstruct(
        double U_LL, double U_L, double U_C, double U_R, double U_RR,
        double& U_recon_L, double& U_recon_R
    ) const override;

    std::string name() const override { return "MUSCL"; }
    int order() const override { return 2; }

private:
    LimiterType limiter_;

    /**
     * Minmod limiter (most conservative).
     * minmod(a, b) = sign(a) * max(0, min(|a|, sign(a)*b))
     */
    double minmod(double a, double b) const {
        if (a * b <= 0.0) return 0.0;
        if (std::abs(a) < std::abs(b)) return a;
        return b;
    }

    /**
     * van Leer limiter (smooth, less diffusive than minmod).
     * van_leer(a, b) = (a*b) / (a + b + 1e-16) if a*b > 0, else 0
     */
    double van_leer(double a, double b) const {
        if (a * b <= 0.0) return 0.0;
        return 2.0 * a * b / (a + b);
    }

    /**
     * Superbee limiter (least diffusive).
     * More aggressive limiting than van Leer.
     */
    double superbee(double a, double b) const {
        double minmod_ab = minmod(a, b);
        double minmod_2a = minmod(2.0 * a, b);
        double minmod_a2b = minmod(a, 2.0 * b);

        if (a * b <= 0.0) return 0.0;
        return minmod(minmod_2a, minmod_a2b);
    }

    /**
     * Apply slope limiter based on type.
     */
    double apply_limiter(double slope_L, double slope_R) const {
        switch (limiter_) {
            case LimiterType::MINMOD:
                return minmod(slope_L, slope_R);
            case LimiterType::VAN_LEER:
                return van_leer(slope_L, slope_R);
            case LimiterType::SUPERBEE:
                return superbee(slope_L, slope_R);
            default:
                return van_leer(slope_L, slope_R);
        }
    }
};

/**
 * Reconstruction scheme factory.
 */
class ReconstructionFactory {
public:
    static std::unique_ptr<ReconstructionScheme> create(
        const std::string& name,
        const std::string& limiter_type = "van_leer"
    );

    static std::vector<std::string> supported_schemes();
};

} // namespace fvm3d::spatial
