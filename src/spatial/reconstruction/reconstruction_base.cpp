#include "spatial/reconstruction/reconstruction_base.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

void ReconstructionMethod::reconstruct_all(
    const core::Field3D<double>& U,
    int direction,
    core::Field3D<double>& U_left,
    core::Field3D<double>& U_right
) const {
    // Default implementation: loop over all interfaces
    // Subclasses should override for better performance

    int nvars = U.nvars();
    int nx = U.nx();
    int ny = U.ny();
    int nz = U.nz();

    Eigen::VectorXd left_state(nvars);
    Eigen::VectorXd right_state(nvars);

    if (direction == 0) {
        // X-direction: reconstruct at i+1/2 interfaces
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    reconstruct(U, i, j, k, direction, left_state, right_state);

                    // Store in output fields
                    for (int v = 0; v < nvars; v++) {
                        U_left(v, i, j, k) = left_state(v);
                        U_right(v, i, j, k) = right_state(v);
                    }
                }
            }
        }
    } else if (direction == 1) {
        // Y-direction: reconstruct at j+1/2 interfaces
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    reconstruct(U, i, j, k, direction, left_state, right_state);

                    for (int v = 0; v < nvars; v++) {
                        U_left(v, i, j, k) = left_state(v);
                        U_right(v, i, j, k) = right_state(v);
                    }
                }
            }
        }
    } else {
        // Z-direction: reconstruct at k+1/2 interfaces
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    reconstruct(U, i, j, k, direction, left_state, right_state);

                    for (int v = 0; v < nvars; v++) {
                        U_left(v, i, j, k) = left_state(v);
                        U_right(v, i, j, k) = right_state(v);
                    }
                }
            }
        }
    }
}

double LimitedReconstructionMethod::apply_limiter(double r) const {
    // Apply slope limiter based on limiter_name_

    if (limiter_name_ == "minmod") {
        // minmod(r) = max(0, min(1, r))
        return std::max(0.0, std::min(1.0, r));
    }
    else if (limiter_name_ == "van_leer" || limiter_name_ == "vanleer") {
        // van Leer: (r + |r|) / (1 + |r|)
        return (r + std::abs(r)) / (1.0 + std::abs(r));
    }
    else if (limiter_name_ == "superbee") {
        // superbee: max(0, min(2r, 1), min(r, 2))
        return std::max({0.0, std::min(2.0 * r, 1.0), std::min(r, 2.0)});
    }
    else if (limiter_name_ == "mc") {
        // MC (Monotonized Central): min(2r, 0.5(1+r), 2)
        return std::max(0.0, std::min({2.0 * r, 0.5 * (1.0 + r), 2.0}));
    }
    else if (limiter_name_ == "none" || limiter_name_ == "unlimited") {
        // No limiting
        return r;
    }
    else {
        // Default: minmod (most conservative)
        return std::max(0.0, std::min(1.0, r));
    }
}

} // namespace fvm3d::spatial
