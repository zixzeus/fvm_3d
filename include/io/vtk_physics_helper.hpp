#ifndef VTK_PHYSICS_HELPER_HPP
#define VTK_PHYSICS_HELPER_HPP

#include "core/field3d.hpp"
#include "core/grid3d.hpp"
#include <vector>
#include <string>

namespace fvm3d::io {

/**
 * @brief Helper class for computing physics-derived quantities for VTK export
 *
 * This class provides static methods to compute primitive variables and derived
 * quantities from conservative state variables, eliminating code duplication
 * between serial and parallel VTK writers.
 *
 * Performance: Functions are marked inline and optimized for hot path execution.
 */
class VTKPhysicsHelper {
public:
    /**
     * @brief Compute Euler equation primitive variables
     *
     * @param state Conservative state (rho, rho_u, rho_v, rho_w, E)
     * @param grid Computational grid
     * @param pressure Output: pressure field
     * @param vel_x Output: x-velocity field
     * @param vel_y Output: y-velocity field
     * @param vel_z Output: z-velocity field
     * @param temperature Output: temperature field (optional, can be nullptr)
     * @param gamma Ratio of specific heats (default: 1.4 for air)
     * @param R Gas constant (default: 1.0)
     */
    static inline void compute_euler_primitives(
        const core::StateField3D& state,
        const core::Grid3D& grid,
        std::vector<double>& pressure,
        std::vector<double>& vel_x,
        std::vector<double>& vel_y,
        std::vector<double>& vel_z,
        std::vector<double>* temperature = nullptr,
        double gamma = 1.4,
        double R = 1.0
    ) __attribute__((hot));

    /**
     * @brief Compute MHD primitive variables
     *
     * @param state Conservative state (rho, rho_u, rho_v, rho_w, E, Bx, By, Bz[, psi])
     * @param grid Computational grid
     * @param pressure Output: pressure field
     * @param vel_x Output: x-velocity field
     * @param vel_y Output: y-velocity field
     * @param vel_z Output: z-velocity field
     * @param B_magnitude Output: magnetic field magnitude
     * @param gamma Ratio of specific heats (default: 5/3 for MHD)
     */
    static inline void compute_mhd_primitives(
        const core::StateField3D& state,
        const core::Grid3D& grid,
        std::vector<double>& pressure,
        std::vector<double>& vel_x,
        std::vector<double>& vel_y,
        std::vector<double>& vel_z,
        std::vector<double>& B_magnitude,
        double gamma = 5.0/3.0
    ) __attribute__((hot));

    /**
     * @brief Extract magnetic field components from state
     *
     * @param state Conservative state containing magnetic field (vars 5, 6, 7)
     * @param grid Computational grid
     * @param Bx Output: x-component of magnetic field
     * @param By Output: y-component of magnetic field
     * @param Bz Output: z-component of magnetic field
     */
    static inline void extract_magnetic_field(
        const core::StateField3D& state,
        const core::Grid3D& grid,
        std::vector<double>& Bx,
        std::vector<double>& By,
        std::vector<double>& Bz
    ) __attribute__((hot));

private:
    // Minimum density threshold to avoid division by zero
    static constexpr double MIN_DENSITY = 1e-10;
};

} // namespace fvm3d::io

#endif // VTK_PHYSICS_HELPER_HPP
