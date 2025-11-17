#pragma once

#include <vector>
#include <string>

namespace fvm3d::spatial {

/**
 * Enumeration for reconstruction variable type.
 *
 * Determines whether a variable should be reconstructed in
 * conservative or primitive form.
 */
enum class ReconstructionVariableType {
    CONSERVATIVE,  ///< Reconstruct in conservative form (e.g., ρ, Bx, By, Bz, ψ)
    PRIMITIVE      ///< Reconstruct in primitive form (e.g., u, v, w, p)
};

/**
 * Configuration for mixed variable reconstruction.
 *
 * Following OpenMHD strategy:
 * - Density, magnetic fields: reconstruct in conservative form
 * - Velocities, pressure: reconstruct in primitive form
 *
 * This prevents negative density and pressure after reconstruction.
 *
 * Reference: OpenMHD 3D_basic/limiter.f90
 */
struct ReconstructionConfig {
    std::vector<ReconstructionVariableType> var_types;

    /**
     * Default configuration for MHD (9 variables).
     *
     * Variable layout:
     *   [0]: ρ    - CONSERVATIVE
     *   [1]: ρu   - PRIMITIVE (as u)
     *   [2]: ρv   - PRIMITIVE (as v)
     *   [3]: ρw   - PRIMITIVE (as w)
     *   [4]: E    - PRIMITIVE (as p)
     *   [5]: Bx   - CONSERVATIVE
     *   [6]: By   - CONSERVATIVE
     *   [7]: Bz   - CONSERVATIVE
     *   [8]: ψ    - CONSERVATIVE
     */
    static ReconstructionConfig default_mhd() {
        ReconstructionConfig config;
        config.var_types = {
            ReconstructionVariableType::CONSERVATIVE,  // 0: ρ
            ReconstructionVariableType::PRIMITIVE,     // 1: u (from ρu)
            ReconstructionVariableType::PRIMITIVE,     // 2: v (from ρv)
            ReconstructionVariableType::PRIMITIVE,     // 3: w (from ρw)
            ReconstructionVariableType::PRIMITIVE,     // 4: p (from E)
            ReconstructionVariableType::CONSERVATIVE,  // 5: Bx
            ReconstructionVariableType::CONSERVATIVE,  // 6: By
            ReconstructionVariableType::CONSERVATIVE,  // 7: Bz
            ReconstructionVariableType::CONSERVATIVE   // 8: ψ
        };
        return config;
    }

    /**
     * All conservative (legacy behavior).
     */
    static ReconstructionConfig all_conservative(int num_vars) {
        ReconstructionConfig config;
        config.var_types.resize(num_vars, ReconstructionVariableType::CONSERVATIVE);
        return config;
    }

    /**
     * Custom configuration.
     */
    static ReconstructionConfig custom(
        const std::vector<ReconstructionVariableType>& types
    ) {
        ReconstructionConfig config;
        config.var_types = types;
        return config;
    }
};

} // namespace fvm3d::spatial
