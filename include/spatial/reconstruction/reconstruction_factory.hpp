#pragma once

#include "reconstruction_base.hpp"
#include <memory>
#include <string>
#include <vector>

namespace fvm3d::spatial {

/**
 * Factory for creating reconstruction methods.
 *
 * Provides a unified interface for creating different reconstruction methods
 * with appropriate parameters.
 *
 * Supported methods:
 * - "constant": First-order piecewise constant (no limiters)
 * - "muscl": Second-order MUSCL with slope limiters
 *
 * Future methods:
 * - "weno3", "weno5": Weighted ENO schemes
 * - "eno": Essentially Non-Oscillatory
 * - "ppm": Piecewise Parabolic Method
 *
 * Example usage:
 *
 *     // Create first-order constant reconstruction
 *     auto recon = ReconstructionFactory::create("constant", 5);
 *
 *     // Create second-order MUSCL with minmod limiter
 *     auto recon = ReconstructionFactory::create("muscl", 5, "minmod");
 *
 *     // Create MUSCL with superbee limiter
 *     auto recon = ReconstructionFactory::create("muscl", 5, "superbee");
 */
class ReconstructionFactory {
public:
    /**
     * Create a reconstruction method.
     *
     * @param name: Method name ("constant", "muscl", "weno3", etc.)
     * @param num_vars: Number of conserved variables
     * @param limiter: Slope limiter for limited methods (default "minmod")
     *                 Options: "minmod", "van_leer", "superbee", "mc"
     * @param kappa: MUSCL parameter for MUSCL method (default 1/3)
     * @return Unique pointer to reconstruction method
     * @throws std::invalid_argument if method name is not recognized
     */
    static std::unique_ptr<ReconstructionMethod> create(
        const std::string& name,
        int num_vars,
        const std::string& limiter = "minmod",
        double kappa = 1.0/3.0
    );

    /**
     * Get list of supported reconstruction methods.
     * @return Vector of method names
     */
    static std::vector<std::string> supported_methods();

    /**
     * Check if a reconstruction method is available.
     * @param name: Method name to check
     * @return True if method is supported
     */
    static bool is_available(const std::string& name);

    /**
     * Print available reconstruction methods to stdout.
     */
    static void print_available();

    /**
     * Get list of supported slope limiters.
     * @return Vector of limiter names
     */
    static std::vector<std::string> supported_limiters();
};

} // namespace fvm3d::spatial
