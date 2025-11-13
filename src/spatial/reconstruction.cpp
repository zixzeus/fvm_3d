#include "spatial/reconstruction.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fvm3d::spatial {

void MUSCLReconstruction::reconstruct(
    double U_LL, double U_L, double U_C, double U_R, double U_RR,
    double& U_recon_L, double& U_recon_R
) const {
    // Compute slopes at left and right sides of center cell
    double slope_L = U_C - U_L;  // Slope to the left
    double slope_R = U_R - U_C;  // Slope to the right

    // Apply slope limiter
    double limited_slope = apply_limiter(slope_L, slope_R);

    // Reconstruct interface values
    U_recon_L = U_C - 0.5 * limited_slope;  // Left interface
    U_recon_R = U_C + 0.5 * limited_slope;  // Right interface
}

std::unique_ptr<ReconstructionScheme> ReconstructionFactory::create(
    const std::string& name,
    const std::string& limiter_type
) {
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (name_lower == "constant") {
        return std::make_unique<ConstantReconstruction>();
    } else if (name_lower == "muscl") {
        std::string limiter_lower = limiter_type;
        std::transform(limiter_lower.begin(), limiter_lower.end(), limiter_lower.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        MUSCLReconstruction::LimiterType limiter = MUSCLReconstruction::LimiterType::VAN_LEER;

        if (limiter_lower == "minmod") {
            limiter = MUSCLReconstruction::LimiterType::MINMOD;
        } else if (limiter_lower == "van_leer") {
            limiter = MUSCLReconstruction::LimiterType::VAN_LEER;
        } else if (limiter_lower == "superbee") {
            limiter = MUSCLReconstruction::LimiterType::SUPERBEE;
        }

        return std::make_unique<MUSCLReconstruction>(limiter);
    } else {
        throw std::invalid_argument("Unknown reconstruction scheme: " + name);
    }
}

std::vector<std::string> ReconstructionFactory::supported_schemes() {
    return {
        "constant",
        "muscl (with minmod, van_leer, superbee limiters)"
    };
}

} // namespace fvm3d::spatial
