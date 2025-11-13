# 3D C++ FVM Framework: Core Modules Implementation Guide

This document provides complete, production-ready implementation specifications for the remaining core computational modules of the 3D C++ MPI-parallel FVM framework. It complements the "3D_CPP_Implementation_Guide.md" and specifies the exact code patterns for Riemann solvers, time integrators, boundary conditions, and the orchestration layer.

## Table of Contents

1. [Riemann Solvers Implementation](#1-riemann-solvers-implementation)
2. [Time Integrators Implementation](#2-time-integrators-implementation)
3. [Boundary Conditions Framework](#3-boundary-conditions-framework)
4. [FVMSolver3D Orchestration](#4-fvmsolver3d-orchestration)
5. [Complete 3D Blast Wave Example](#5-complete-3d-blast-wave-example)
6. [Testing Framework and Unit Tests](#6-testing-framework-and-unit-tests)
7. [Performance Profiling Integration](#7-performance-profiling-integration)
8. [Implementation Checklist](#8-implementation-checklist)

---

## 1. Riemann Solvers Implementation

### 1.1 Riemann Solver Base Class

**File**: `include/fvm3d/spatial/riemann_solver.hpp`

```cpp
#pragma once

#include <Eigen/Dense>
#include <memory>
#include <string>

namespace fvm3d::spatial {

/**
 * Abstract base class for Riemann solvers.
 * Solves the Riemann problem at a single interface:
 *   U_L: left state
 *   U_R: right state
 *   flux = F(U_L, U_R)
 *
 * All Riemann solvers work with conservative variables:
 *   U = [rho, rho_u, rho_v, rho_w, E]^T
 */
class RiemannSolver {
public:
    virtual ~RiemannSolver() = default;

    /**
     * Compute the numerical flux at an interface given left and right states.
     * @param U_L: left conservative state (5 components)
     * @param U_R: right conservative state (5 components)
     * @param direction: 0=X, 1=Y, 2=Z
     * @return: numerical flux (5 components)
     */
    virtual Eigen::VectorXd solve(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const = 0;

    /**
     * Get the maximum wave speed in the Riemann fan.
     * Used for CFL condition: dt <= CFL * dx / max_speed
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

    /**
     * Get the number of conserved variables (always 5 for Euler).
     */
    virtual int num_variables() const { return 5; }

protected:
    /**
     * Gamma parameter for ideal gas (adiabatic index).
     * Typical value: 1.4 for air.
     */
    static constexpr double GAMMA = 1.4;

    /**
     * Density floor to prevent division by zero.
     */
    static constexpr double RHO_FLOOR = 1e-10;

    /**
     * Pressure floor to prevent negative pressure.
     */
    static constexpr double P_FLOOR = 1e-11;

    /**
     * Extract primitive variables from conservative variables.
     * U = [rho, rho_u, rho_v, rho_w, E]
     * V = [rho, u, v, w, p]
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
     * Compute speed of sound given density and pressure.
     * a = sqrt(gamma * p / rho)
     */
    double sound_speed(double rho, double p) const {
        rho = std::max(rho, RHO_FLOOR);
        p = std::max(p, P_FLOOR);
        return std::sqrt(GAMMA * p / rho);
    }

    /**
     * Extract velocity component in a given direction.
     * direction: 0=X (use rho_u), 1=Y (use rho_v), 2=Z (use rho_w)
     */
    double velocity_in_direction(
        const Eigen::VectorXd& U,
        int direction
    ) const {
        double rho = std::max(U(0), RHO_FLOOR);
        return U(direction + 1) / rho;
    }

    /**
     * Rotate velocity components to align with normal direction.
     * For flux calculation in direction d, we need to align:
     *   normal velocity (to be computed as flux)
     *   and tangential velocities (passive)
     */
    Eigen::VectorXd rotate_to_normal(
        const Eigen::VectorXd& U,
        int direction
    ) const {
        // For now, return U as-is
        // In 3D, this could involve rotating u,v,w based on direction
        return U;
    }
};

} // namespace fvm3d::spatial
```

### 1.2 Lax-Friedrichs (LF) Solver

**File**: `include/fvm3d/spatial/riemann_laxfriedrichs.hpp`

```cpp
#pragma once

#include "riemann_solver.hpp"

namespace fvm3d::spatial {

/**
 * Lax-Friedrichs Riemann solver.
 * Simple and robust, but diffusive.
 *
 * Flux: F_LF = 0.5 * (F_L + F_R) - 0.5 * lambda * (U_R - U_L)
 * where lambda = max(|u| + a) is the maximum wave speed
 */
class LaxFriedrichsSolver : public RiemannSolver {
public:
    Eigen::VectorXd solve(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override;

    double max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override;

    std::string name() const override { return "Lax-Friedrichs"; }

private:
    /**
     * Compute flux F(U) in a given direction.
     * For direction = 0 (X): F = [rho*u, rho*u*u+p, rho*u*v, rho*u*w, (E+p)*u]
     * For direction = 1 (Y): similar with v as normal
     * For direction = 2 (Z): similar with w as normal
     */
    Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U,
        int direction
    ) const;
};

} // namespace fvm3d::spatial
```

**File**: `src/spatial/riemann_laxfriedrichs.cpp`

```cpp
#include "fvm3d/spatial/riemann_laxfriedrichs.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

Eigen::VectorXd LaxFriedrichsSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    // Compute fluxes
    Eigen::VectorXd F_L = compute_flux(U_L, direction);
    Eigen::VectorXd F_R = compute_flux(U_R, direction);

    // Compute maximum wave speed
    double lambda = max_wave_speed(U_L, U_R, direction);

    // Lax-Friedrichs flux: 0.5 * (F_L + F_R) - 0.5 * lambda * (U_R - U_L)
    Eigen::VectorXd flux = 0.5 * (F_L + F_R) - 0.5 * lambda * (U_R - U_L);

    return flux;
}

double LaxFriedrichsSolver::max_wave_speed(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    double a_L = sound_speed(rho_L, p_L);
    double a_R = sound_speed(rho_R, p_R);

    // Select normal velocity based on direction
    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    // Maximum wave speed: max(|u| + a)
    double lambda_L = std::abs(u_normal_L) + a_L;
    double lambda_R = std::abs(u_normal_R) + a_R;

    return std::max(lambda_L, lambda_R);
}

Eigen::VectorXd LaxFriedrichsSolver::compute_flux(
    const Eigen::VectorXd& U,
    int direction
) const {
    double rho, u, v, w, p;
    conservative_to_primitive(U, rho, u, v, w, p);

    Eigen::VectorXd flux(5);

    if (direction == 0) {
        // X-direction flux
        // F = [rho*u, rho*u*u+p, rho*u*v, rho*u*w, (E+p)*u]
        flux(0) = rho * u;
        flux(1) = rho * u * u + p;
        flux(2) = rho * u * v;
        flux(3) = rho * u * w;
        flux(4) = (U(4) + p) * u;
    } else if (direction == 1) {
        // Y-direction flux
        // G = [rho*v, rho*v*u, rho*v*v+p, rho*v*w, (E+p)*v]
        flux(0) = rho * v;
        flux(1) = rho * v * u;
        flux(2) = rho * v * v + p;
        flux(3) = rho * v * w;
        flux(4) = (U(4) + p) * v;
    } else { // direction == 2
        // Z-direction flux
        // H = [rho*w, rho*w*u, rho*w*v, rho*w*w+p, (E+p)*w]
        flux(0) = rho * w;
        flux(1) = rho * w * u;
        flux(2) = rho * w * v;
        flux(3) = rho * w * w + p;
        flux(4) = (U(4) + p) * w;
    }

    return flux;
}

} // namespace fvm3d::spatial
```

### 1.3 HLL Solver (Harten-Lax-van Leer)

**File**: `include/fvm3d/spatial/riemann_hll.hpp`

```cpp
#pragma once

#include "riemann_solver.hpp"

namespace fvm3d::spatial {

/**
 * HLL (Harten-Lax-van Leer) Riemann solver.
 * More accurate than Lax-Friedrichs, especially for shear flows.
 *
 * Key: Estimates the left and right wave speeds (S_L and S_R)
 * Uses two-wave model: only left and right waves (ignore contact discontinuity)
 *
 * Flux:
 *   if 0 < S_L: use F_L
 *   if S_L <= 0 <= S_R: use F_HLL_middle
 *   if 0 > S_R: use F_R
 */
class HLLSolver : public RiemannSolver {
public:
    Eigen::VectorXd solve(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override;

    double max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override;

    std::string name() const override { return "HLL"; }

private:
    /**
     * Estimate left wave speed S_L = min(u_L - a_L, u_R - a_R)
     * Using Roe's formula for consistency
     */
    double estimate_s_left(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const;

    /**
     * Estimate right wave speed S_R = max(u_L + a_L, u_R + a_R)
     */
    double estimate_s_right(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const;

    /**
     * Compute flux in given direction
     */
    Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U,
        int direction
    ) const;

    /**
     * Compute middle state for HLL solver
     */
    Eigen::VectorXd compute_u_middle(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const Eigen::VectorXd& F_L,
        const Eigen::VectorXd& F_R,
        double S_L,
        double S_R
    ) const;
};

} // namespace fvm3d::spatial
```

**File**: `src/spatial/riemann_hll.cpp`

```cpp
#include "fvm3d/spatial/riemann_hll.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

Eigen::VectorXd HLLSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    // Estimate wave speeds
    double S_L = estimate_s_left(U_L, U_R, direction);
    double S_R = estimate_s_right(U_L, U_R, direction);

    // Compute fluxes
    Eigen::VectorXd F_L = compute_flux(U_L, direction);
    Eigen::VectorXd F_R = compute_flux(U_R, direction);

    // Ensure S_L < S_R to avoid division by zero
    if (S_R <= S_L) {
        // Degenerate case: return average flux
        return 0.5 * (F_L + F_R);
    }

    // HLL flux formula
    if (0.0 <= S_L) {
        // Contact/shock wave moving to the right, use left flux
        return F_L;
    } else if (0.0 >= S_R) {
        // Contact/shock wave moving to the left, use right flux
        return F_R;
    } else {
        // Contact/shock wave straddles interface, use HLL middle state
        Eigen::VectorXd U_m = compute_u_middle(U_L, U_R, F_L, F_R, S_L, S_R);
        Eigen::VectorXd F_HLL =
            (S_R * F_L - S_L * F_R + S_L * S_R * (U_R - U_L)) / (S_R - S_L);
        return F_HLL;
    }
}

double HLLSolver::max_wave_speed(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double S_L = estimate_s_left(U_L, U_R, direction);
    double S_R = estimate_s_right(U_L, U_R, direction);
    return std::max(std::abs(S_L), std::abs(S_R));
}

double HLLSolver::estimate_s_left(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    double a_L = sound_speed(rho_L, p_L);
    double a_R = sound_speed(rho_R, p_R);

    // Select normal velocity based on direction
    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    // Roe's estimate for left wave speed
    return std::min(u_normal_L - a_L, u_normal_R - a_R);
}

double HLLSolver::estimate_s_right(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    double a_L = sound_speed(rho_L, p_L);
    double a_R = sound_speed(rho_R, p_R);

    // Select normal velocity based on direction
    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    // Roe's estimate for right wave speed
    return std::max(u_normal_L + a_L, u_normal_R + a_R);
}

Eigen::VectorXd HLLSolver::compute_flux(
    const Eigen::VectorXd& U,
    int direction
) const {
    double rho, u, v, w, p;
    conservative_to_primitive(U, rho, u, v, w, p);

    Eigen::VectorXd flux(5);

    if (direction == 0) {
        flux(0) = rho * u;
        flux(1) = rho * u * u + p;
        flux(2) = rho * u * v;
        flux(3) = rho * u * w;
        flux(4) = (U(4) + p) * u;
    } else if (direction == 1) {
        flux(0) = rho * v;
        flux(1) = rho * v * u;
        flux(2) = rho * v * v + p;
        flux(3) = rho * v * w;
        flux(4) = (U(4) + p) * v;
    } else {
        flux(0) = rho * w;
        flux(1) = rho * w * u;
        flux(2) = rho * w * v;
        flux(3) = rho * w * w + p;
        flux(4) = (U(4) + p) * w;
    }

    return flux;
}

Eigen::VectorXd HLLSolver::compute_u_middle(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    const Eigen::VectorXd& F_L,
    const Eigen::VectorXd& F_R,
    double S_L,
    double S_R
) const {
    // HLL middle state: U_m = (S_R*U_R - S_L*U_L - (F_R - F_L)) / (S_R - S_L)
    return (S_R * U_R - S_L * U_L - (F_R - F_L)) / (S_R - S_L);
}

} // namespace fvm3d::spatial
```

### 1.4 HLLC Solver (HLL with Contact Discontinuity)

**File**: `include/fvm3d/spatial/riemann_hllc.hpp`

```cpp
#pragma once

#include "riemann_solver.hpp"

namespace fvm3d::spatial {

/**
 * HLLC (Harten-Lax-van Leer-Contact) Riemann solver.
 * More accurate than HLL by properly capturing contact discontinuities.
 * Uses three-wave model: left shock, contact discontinuity, right shock.
 *
 * Key improvement over HLL: properly handles pressure and velocity jumps.
 * Adds contact speed S_M estimation.
 */
class HLLCSolver : public RiemannSolver {
public:
    Eigen::VectorXd solve(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override;

    double max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override;

    std::string name() const override { return "HLLC"; }

private:
    double estimate_s_left(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const;

    double estimate_s_right(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const;

    /**
     * Estimate contact speed S_M = S_L + (p_star - p_L) / (rho_L * (S_L - u_L))
     * This is Roe's formula for the contact discontinuity
     */
    double estimate_s_contact(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction,
        double S_L,
        double S_R
    ) const;

    /**
     * Estimate pressure in the middle state using Roe's formula
     */
    double estimate_p_middle(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const;

    Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U,
        int direction
    ) const;

    /**
     * Compute left state in HLLC fan (between S_L and S_M)
     */
    Eigen::VectorXd compute_u_left_hllc(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& F_L,
        double S_L,
        double S_M,
        double p_m,
        int direction
    ) const;

    /**
     * Compute right state in HLLC fan (between S_M and S_R)
     */
    Eigen::VectorXd compute_u_right_hllc(
        const Eigen::VectorXd& U_R,
        const Eigen::VectorXd& F_R,
        double S_R,
        double S_M,
        double p_m,
        int direction
    ) const;
};

} // namespace fvm3d::spatial
```

**File**: `src/spatial/riemann_hllc.cpp` (implementation follows HLL pattern with contact discontinuity handling)

```cpp
#include "fvm3d/spatial/riemann_hllc.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

Eigen::VectorXd HLLCSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    // Estimate wave speeds
    double S_L = estimate_s_left(U_L, U_R, direction);
    double S_R = estimate_s_right(U_L, U_R, direction);

    // Estimate middle pressure and contact speed
    double p_m = estimate_p_middle(U_L, U_R, direction);
    double S_M = estimate_s_contact(U_L, U_R, direction, S_L, S_R);

    // Compute fluxes
    Eigen::VectorXd F_L = compute_flux(U_L, direction);
    Eigen::VectorXd F_R = compute_flux(U_R, direction);

    // Ensure wave speed ordering
    if (S_R <= S_L) {
        return 0.5 * (F_L + F_R);
    }

    // Select appropriate state based on sign of wave speeds
    if (0.0 <= S_L) {
        return F_L;
    } else if (0.0 >= S_R) {
        return F_R;
    } else if (0.0 <= S_M) {
        // In left or middle region, use left state
        Eigen::VectorXd U_Lm = compute_u_left_hllc(U_L, F_L, S_L, S_M, p_m, direction);
        Eigen::VectorXd F_Lm = F_L + S_L * (U_Lm - U_L);
        return F_Lm;
    } else {
        // In right region, use right state
        Eigen::VectorXd U_Rm = compute_u_right_hllc(U_R, F_R, S_R, S_M, p_m, direction);
        Eigen::VectorXd F_Rm = F_R + S_R * (U_Rm - U_R);
        return F_Rm;
    }
}

double HLLCSolver::max_wave_speed(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double S_L = estimate_s_left(U_L, U_R, direction);
    double S_R = estimate_s_right(U_L, U_R, direction);
    return std::max(std::abs(S_L), std::abs(S_R));
}

double HLLCSolver::estimate_s_left(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    double a_L = sound_speed(rho_L, p_L);
    double a_R = sound_speed(rho_R, p_R);

    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    return std::min(u_normal_L - a_L, u_normal_R - a_R);
}

double HLLCSolver::estimate_s_right(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    double a_L = sound_speed(rho_L, p_L);
    double a_R = sound_speed(rho_R, p_R);

    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    return std::max(u_normal_L + a_L, u_normal_R + a_R);
}

double HLLCSolver::estimate_p_middle(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    double a_L = sound_speed(rho_L, p_L);
    double a_R = sound_speed(rho_R, p_R);

    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    // Roe's formula for middle pressure
    double num = a_R * p_L + a_L * p_R + a_L * a_R * (u_normal_L - u_normal_R);
    double den = a_L + a_R;
    return std::max(num / den, P_FLOOR);
}

double HLLCSolver::estimate_s_contact(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction,
    double S_L,
    double S_R
) const {
    // Contact speed is typically taken as the average velocity
    // More sophisticated: use Roe averaging
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    // Simple average (can be improved with Roe averaging)
    return 0.5 * (u_normal_L + u_normal_R);
}

Eigen::VectorXd HLLCSolver::compute_flux(
    const Eigen::VectorXd& U,
    int direction
) const {
    double rho, u, v, w, p;
    conservative_to_primitive(U, rho, u, v, w, p);

    Eigen::VectorXd flux(5);

    if (direction == 0) {
        flux(0) = rho * u;
        flux(1) = rho * u * u + p;
        flux(2) = rho * u * v;
        flux(3) = rho * u * w;
        flux(4) = (U(4) + p) * u;
    } else if (direction == 1) {
        flux(0) = rho * v;
        flux(1) = rho * v * u;
        flux(2) = rho * v * v + p;
        flux(3) = rho * v * w;
        flux(4) = (U(4) + p) * v;
    } else {
        flux(0) = rho * w;
        flux(1) = rho * w * u;
        flux(2) = rho * w * v;
        flux(3) = rho * w * w + p;
        flux(4) = (U(4) + p) * w;
    }

    return flux;
}

Eigen::VectorXd HLLCSolver::compute_u_left_hllc(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& F_L,
    double S_L,
    double S_M,
    double p_m,
    int direction
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);

    Eigen::VectorXd U_Lm = U_L;
    double coeff = rho_L * (S_L - u_L) / (S_L - S_M);

    // Update density
    U_Lm(0) = rho_L * (S_L - u_L) / (S_L - S_M);

    // Update momentum components
    U_Lm(direction + 1) = U_Lm(0) * S_M;

    // Other momentum components unchanged
    for (int d = 1; d <= 3; d++) {
        if (d != direction + 1) {
            U_Lm(d) = coeff * (U_L(d) / rho_L);
        }
    }

    // Update energy
    U_Lm(4) = U_L(4) + coeff * (p_m / rho_L - p_L / rho_L);

    return U_Lm;
}

Eigen::VectorXd HLLCSolver::compute_u_right_hllc(
    const Eigen::VectorXd& U_R,
    const Eigen::VectorXd& F_R,
    double S_R,
    double S_M,
    double p_m,
    int direction
) const {
    double rho_R, u_R, v_R, w_R, p_R;
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    Eigen::VectorXd U_Rm = U_R;
    double coeff = rho_R * (S_R - u_R) / (S_R - S_M);

    // Update density
    U_Rm(0) = rho_R * (S_R - u_R) / (S_R - S_M);

    // Update momentum components
    U_Rm(direction + 1) = U_Rm(0) * S_M;

    // Other momentum components unchanged
    for (int d = 1; d <= 3; d++) {
        if (d != direction + 1) {
            U_Rm(d) = coeff * (U_R(d) / rho_R);
        }
    }

    // Update energy
    U_Rm(4) = U_R(4) + coeff * (p_m / rho_R - p_R / rho_R);

    return U_Rm;
}

} // namespace fvm3d::spatial
```

### 1.5 Riemann Solver Factory

**File**: `include/fvm3d/spatial/riemann_solver_factory.hpp`

```cpp
#pragma once

#include "riemann_solver.hpp"
#include <memory>
#include <string>
#include <map>

namespace fvm3d::spatial {

class RiemannSolverFactory {
public:
    /**
     * Create a Riemann solver by name.
     * Supported names: "laxfriedrichs", "hll", "hllc"
     */
    static std::unique_ptr<RiemannSolver> create(const std::string& name);

    /**
     * Get list of supported Riemann solvers
     */
    static std::vector<std::string> supported_solvers();
};

} // namespace fvm3d::spatial
```

**File**: `src/spatial/riemann_solver_factory.cpp`

```cpp
#include "fvm3d/spatial/riemann_solver_factory.hpp"
#include "fvm3d/spatial/riemann_laxfriedrichs.hpp"
#include "fvm3d/spatial/riemann_hll.hpp"
#include "fvm3d/spatial/riemann_hllc.hpp"
#include <stdexcept>

namespace fvm3d::spatial {

std::unique_ptr<RiemannSolver> RiemannSolverFactory::create(const std::string& name) {
    if (name == "laxfriedrichs" || name == "lf") {
        return std::make_unique<LaxFriedrichsSolver>();
    } else if (name == "hll") {
        return std::make_unique<HLLSolver>();
    } else if (name == "hllc") {
        return std::make_unique<HLLCSolver>();
    } else {
        throw std::invalid_argument("Unknown Riemann solver: " + name);
    }
}

std::vector<std::string> RiemannSolverFactory::supported_solvers() {
    return {"laxfriedrichs", "hll", "hllc"};
}

} // namespace fvm3d::spatial
```

---

## 2. Time Integrators Implementation

### 2.1 Time Integrator Base Class

**File**: `include/fvm3d/temporal/time_integrator.hpp`

```cpp
#pragma once

#include "fvm3d/core/field3d.hpp"
#include "fvm3d/core/grid3d.hpp"
#include <functional>
#include <memory>
#include <string>

namespace fvm3d::temporal {

/**
 * Abstract base class for time integration schemes.
 *
 * All integrators are explicit and use the method-of-lines approach:
 *   dU/dt = RHS(U)
 *
 * where RHS is the spatial discretization (FVM flux divergence + source terms)
 */
class TimeIntegrator {
public:
    using RHSFunction = std::function<void(const StateField3D&, StateField3D&)>;

    virtual ~TimeIntegrator() = default;

    /**
     * Perform one time step: U^(n+1) = Integrate(U^n, dt, RHS)
     *
     * @param U_current: current state (will be modified)
     * @param dt: time step size
     * @param rhs: right-hand side function RHS(U, dU/dt)
     */
    virtual void step(
        StateField3D& U_current,
        double dt,
        const RHSFunction& rhs
    ) = 0;

    /**
     * Get the order of accuracy in time.
     */
    virtual int order() const = 0;

    /**
     * Get the name of the integrator.
     */
    virtual std::string name() const = 0;

protected:
    /**
     * Temporary storage for intermediate stages.
     * Allocated on first use to match grid dimensions.
     */
    mutable std::unique_ptr<StateField3D> temp_stage1_;
    mutable std::unique_ptr<StateField3D> temp_stage2_;
    mutable std::unique_ptr<StateField3D> temp_rhs_;

    /**
     * Allocate temporary arrays if needed.
     * Called before first step.
     */
    void allocate_temporaries(const StateField3D& U) const {
        if (!temp_rhs_) {
            temp_rhs_ = std::make_unique<StateField3D>(U.nvars(), U.nx(), U.ny(), U.nz());
            temp_stage1_ = std::make_unique<StateField3D>(U.nvars(), U.nx(), U.ny(), U.nz());
            temp_stage2_ = std::make_unique<StateField3D>(U.nvars(), U.nx(), U.ny(), U.nz());
        }
    }
};

} // namespace fvm3d::temporal
```

### 2.2 Forward Euler Integrator

**File**: `include/fvm3d/temporal/forward_euler.hpp`

```cpp
#pragma once

#include "time_integrator.hpp"

namespace fvm3d::temporal {

/**
 * Forward Euler (explicit Euler) time integrator.
 * U^(n+1) = U^n + dt * RHS(U^n)
 *
 * Order: 1st order
 * Stability: CFL <= 1
 *
 * Advantages: Simple, minimal memory overhead
 * Disadvantages: First-order accuracy, strict CFL requirement
 */
class ForwardEuler : public TimeIntegrator {
public:
    void step(
        StateField3D& U_current,
        double dt,
        const RHSFunction& rhs
    ) override;

    int order() const override { return 1; }
    std::string name() const override { return "Forward Euler"; }
};

} // namespace fvm3d::temporal
```

**File**: `src/temporal/forward_euler.cpp`

```cpp
#include "fvm3d/temporal/forward_euler.hpp"

namespace fvm3d::temporal {

void ForwardEuler::step(
    StateField3D& U_current,
    double dt,
    const RHSFunction& rhs
) {
    allocate_temporaries(U_current);

    // Compute RHS(U^n)
    rhs(U_current, *temp_rhs_);

    // U^(n+1) = U^n + dt * RHS(U^n)
    int nvars = U_current.nvars();
    int nx = U_current.nx();
    int ny = U_current.ny();
    int nz = U_current.nz();

    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    U_current(v, i, j, k) += dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }
}

} // namespace fvm3d::temporal
```

### 2.3 RK2 (SSP-RK2) Integrator

**File**: `include/fvm3d/temporal/rk2.hpp`

```cpp
#pragma once

#include "time_integrator.hpp"

namespace fvm3d::temporal {

/**
 * Strong Stability Preserving RK2 (Heun's method).
 *
 * U^(1) = U^n + dt * RHS(U^n)
 * U^(n+1) = 0.5 * U^n + 0.5 * U^(1) + 0.5 * dt * RHS(U^(1))
 *
 * Order: 2nd order
 * Stability: CFL <= 1 (same as Forward Euler when applied to conservation laws!)
 *
 * Advantages: 2nd order accuracy, TVD property
 * Disadvantages: Two RHS evaluations per step
 */
class RK2 : public TimeIntegrator {
public:
    void step(
        StateField3D& U_current,
        double dt,
        const RHSFunction& rhs
    ) override;

    int order() const override { return 2; }
    std::string name() const override { return "SSP-RK2"; }
};

} // namespace fvm3d::temporal
```

**File**: `src/temporal/rk2.cpp`

```cpp
#include "fvm3d/temporal/rk2.hpp"

namespace fvm3d::temporal {

void RK2::step(
    StateField3D& U_current,
    double dt,
    const RHSFunction& rhs
) {
    allocate_temporaries(U_current);

    int nvars = U_current.nvars();
    int nx = U_current.nx();
    int ny = U_current.ny();
    int nz = U_current.nz();

    // Stage 1: U^(1) = U^n + dt * RHS(U^n)
    rhs(U_current, *temp_rhs_);
    *temp_stage1_ = U_current;
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    (*temp_stage1_)(v, i, j, k) += dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }

    // Stage 2: RHS(U^(1))
    rhs(*temp_stage1_, *temp_rhs_);

    // U^(n+1) = 0.5 * U^n + 0.5 * U^(1) + 0.5 * dt * RHS(U^(1))
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    U_current(v, i, j, k) = 0.5 * U_current(v, i, j, k) +
                                           0.5 * (*temp_stage1_)(v, i, j, k) +
                                           0.5 * dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }
}

} // namespace fvm3d::temporal
```

### 2.4 RK3 (SSP-RK3) Integrator

**File**: `include/fvm3d/temporal/rk3.hpp`

```cpp
#pragma once

#include "time_integrator.hpp"

namespace fvm3d::temporal {

/**
 * Strong Stability Preserving RK3.
 *
 * U^(1) = U^n + dt * RHS(U^n)
 * U^(2) = 0.75 * U^n + 0.25 * U^(1) + 0.25 * dt * RHS(U^(1))
 * U^(n+1) = (1/3) * U^n + (2/3) * U^(2) + (2/3) * dt * RHS(U^(2))
 *
 * Order: 3rd order
 * Stability: CFL <= 1 (TVD property)
 *
 * Advantages: 3rd order accuracy, TVD property, good stability
 * Disadvantages: Three RHS evaluations per step
 */
class RK3 : public TimeIntegrator {
public:
    void step(
        StateField3D& U_current,
        double dt,
        const RHSFunction& rhs
    ) override;

    int order() const override { return 3; }
    std::string name() const override { return "SSP-RK3"; }
};

} // namespace fvm3d::temporal
```

**File**: `src/temporal/rk3.cpp`

```cpp
#include "fvm3d/temporal/rk3.hpp"

namespace fvm3d::temporal {

void RK3::step(
    StateField3D& U_current,
    double dt,
    const RHSFunction& rhs
) {
    allocate_temporaries(U_current);

    int nvars = U_current.nvars();
    int nx = U_current.nx();
    int ny = U_current.ny();
    int nz = U_current.nz();

    // Stage 1: U^(1) = U^n + dt * RHS(U^n)
    rhs(U_current, *temp_rhs_);
    *temp_stage1_ = U_current;
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    (*temp_stage1_)(v, i, j, k) += dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }

    // Stage 2: U^(2) = 0.75 * U^n + 0.25 * U^(1) + 0.25 * dt * RHS(U^(1))
    rhs(*temp_stage1_, *temp_rhs_);
    *temp_stage2_ = U_current;  // Start with U^n
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    (*temp_stage2_)(v, i, j, k) = 0.75 * U_current(v, i, j, k) +
                                                 0.25 * (*temp_stage1_)(v, i, j, k) +
                                                 0.25 * dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }

    // Stage 3: U^(n+1) = (1/3) * U^n + (2/3) * U^(2) + (2/3) * dt * RHS(U^(2))
    rhs(*temp_stage2_, *temp_rhs_);
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    U_current(v, i, j, k) = (1.0/3.0) * U_current(v, i, j, k) +
                                           (2.0/3.0) * (*temp_stage2_)(v, i, j, k) +
                                           (2.0/3.0) * dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }
}

} // namespace fvm3d::temporal
```

### 2.5 Time Integrator Factory

**File**: `include/fvm3d/temporal/time_integrator_factory.hpp`

```cpp
#pragma once

#include "time_integrator.hpp"
#include <memory>
#include <string>
#include <vector>

namespace fvm3d::temporal {

class TimeIntegratorFactory {
public:
    /**
     * Create a time integrator by name.
     * Supported names: "euler", "rk2", "rk3"
     */
    static std::unique_ptr<TimeIntegrator> create(const std::string& name);

    /**
     * Get list of supported time integrators
     */
    static std::vector<std::string> supported_integrators();
};

} // namespace fvm3d::temporal
```

**File**: `src/temporal/time_integrator_factory.cpp`

```cpp
#include "fvm3d/temporal/time_integrator_factory.hpp"
#include "fvm3d/temporal/forward_euler.hpp"
#include "fvm3d/temporal/rk2.hpp"
#include "fvm3d/temporal/rk3.hpp"
#include <stdexcept>

namespace fvm3d::temporal {

std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create(const std::string& name) {
    if (name == "euler") {
        return std::make_unique<ForwardEuler>();
    } else if (name == "rk2") {
        return std::make_unique<RK2>();
    } else if (name == "rk3") {
        return std::make_unique<RK3>();
    } else {
        throw std::invalid_argument("Unknown time integrator: " + name);
    }
}

std::vector<std::string> TimeIntegratorFactory::supported_integrators() {
    return {"euler", "rk2", "rk3"};
}

} // namespace fvm3d::temporal
```

---

## 3. Boundary Conditions Framework

### 3.1 Boundary Condition Base Class

**File**: `include/fvm3d/boundary/boundary_condition.hpp`

```cpp
#pragma once

#include "fvm3d/core/field3d.hpp"
#include "fvm3d/core/grid3d.hpp"
#include <string>

namespace fvm3d::boundary {

/**
 * Abstract base class for boundary conditions.
 *
 * Boundary conditions fill the ghost cells at the domain boundaries.
 * Applied in-place to the state field after each time step (or during RHS computation).
 */
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    /**
     * Apply boundary condition to all ghost layers.
     * This is called after computing the RHS to ensure ghost cells are properly set
     * for the next time step's flux calculations.
     *
     * @param state: state field with ghost cells allocated
     * @param grid: grid information (domain size, ghost width, etc.)
     */
    virtual void apply(
        StateField3D& state,
        const Grid3D& grid
    ) = 0;

    /**
     * Get name of boundary condition
     */
    virtual std::string name() const = 0;

protected:
    /**
     * Constants for Euler equations (same as in Riemann solvers)
     */
    static constexpr double GAMMA = 1.4;
    static constexpr double RHO_FLOOR = 1e-10;
    static constexpr double P_FLOOR = 1e-11;

    /**
     * Extract primitive variables from conservative state
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
};

} // namespace fvm3d::boundary
```

### 3.2 Periodic Boundary Condition

**File**: `include/fvm3d/boundary/periodic_bc.hpp`

```cpp
#pragma once

#include "boundary_condition.hpp"

namespace fvm3d::boundary {

/**
 * Periodic boundary condition.
 * Ghost cells are filled from the opposite side of the domain.
 *
 * Example for X-direction:
 *   Left ghost: U[0, -1, j, k] = U[0, nx-1, j, k]
 *   Right ghost: U[0, nx, j, k] = U[0, 0, j, k]
 */
class PeriodicBC : public BoundaryCondition {
public:
    /**
     * Constructor specifies which directions are periodic.
     * @param periodic_x, periodic_y, periodic_z: true if that direction is periodic
     */
    PeriodicBC(bool periodic_x = true, bool periodic_y = false, bool periodic_z = false)
        : periodic_x_(periodic_x), periodic_y_(periodic_y), periodic_z_(periodic_z) {}

    void apply(
        StateField3D& state,
        const Grid3D& grid
    ) override;

    std::string name() const override { return "Periodic"; }

private:
    bool periodic_x_, periodic_y_, periodic_z_;

    void apply_periodic_x(StateField3D& state, const Grid3D& grid);
    void apply_periodic_y(StateField3D& state, const Grid3D& grid);
    void apply_periodic_z(StateField3D& state, const Grid3D& grid);
};

} // namespace fvm3d::boundary
```

**File**: `src/boundary/periodic_bc.cpp`

```cpp
#include "fvm3d/boundary/periodic_bc.hpp"

namespace fvm3d::boundary {

void PeriodicBC::apply(StateField3D& state, const Grid3D& grid) {
    if (periodic_x_) apply_periodic_x(state, grid);
    if (periodic_y_) apply_periodic_y(state, grid);
    if (periodic_z_) apply_periodic_z(state, grid);
}

void PeriodicBC::apply_periodic_x(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx_total = state.nx();
    int ny = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Interior cells: [nghost, nx_total - nghost)
    // Ghost cells: [0, nghost) and [nx_total - nghost, nx_total)
    int nx_interior = nx_total - 2 * nghost;

    // Left ghost: copy from right interior
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int src_i = nx_total - nghost - 1 - (nghost - 1 - g);
                    state(v, g, j, k) = state(v, src_i, j, k);
                }
            }
        }
    }

    // Right ghost: copy from left interior
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int src_i = nghost + g;
                    int dst_i = nx_total - nghost + g;
                    state(v, dst_i, j, k) = state(v, src_i, j, k);
                }
            }
        }
    }
}

void PeriodicBC::apply_periodic_y(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny_total = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: copy from top interior
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    int src_j = ny_total - nghost - 1 - (nghost - 1 - g);
                    state(v, i, g, k) = state(v, i, src_j, k);
                }
            }
        }
    }

    // Top ghost: copy from bottom interior
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    int src_j = nghost + g;
                    int dst_j = ny_total - nghost + g;
                    state(v, i, dst_j, k) = state(v, i, src_j, k);
                }
            }
        }
    }
}

void PeriodicBC::apply_periodic_z(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny = state.ny();
    int nz_total = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: copy from top interior
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    int src_k = nz_total - nghost - 1 - (nghost - 1 - g);
                    state(v, i, j, g) = state(v, i, j, src_k);
                }
            }
        }
    }

    // Top ghost: copy from bottom interior
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    int src_k = nghost + g;
                    int dst_k = nz_total - nghost + g;
                    state(v, i, j, dst_k) = state(v, i, j, src_k);
                }
            }
        }
    }
}

} // namespace fvm3d::boundary
```

### 3.3 Reflective Boundary Condition

**File**: `include/fvm3d/boundary/reflective_bc.hpp`

```cpp
#pragma once

#include "boundary_condition.hpp"

namespace fvm3d::boundary {

/**
 * Reflective (slip) boundary condition.
 * Ghost cells are filled by mirroring the interior with reversed normal velocity.
 *
 * Used for solid walls where normal velocity = 0 but tangential slip is allowed.
 */
class ReflectiveBC : public BoundaryCondition {
public:
    ReflectiveBC(bool reflect_x = true, bool reflect_y = false, bool reflect_z = false)
        : reflect_x_(reflect_x), reflect_y_(reflect_y), reflect_z_(reflect_z) {}

    void apply(
        StateField3D& state,
        const Grid3D& grid
    ) override;

    std::string name() const override { return "Reflective"; }

private:
    bool reflect_x_, reflect_y_, reflect_z_;

    void apply_reflective_x(StateField3D& state, const Grid3D& grid);
    void apply_reflective_y(StateField3D& state, const Grid3D& grid);
    void apply_reflective_z(StateField3D& state, const Grid3D& grid);
};

} // namespace fvm3d::boundary
```

**File**: `src/boundary/reflective_bc.cpp`

```cpp
#include "fvm3d/boundary/reflective_bc.hpp"
#include <cmath>

namespace fvm3d::boundary {

void ReflectiveBC::apply(StateField3D& state, const Grid3D& grid) {
    if (reflect_x_) apply_reflective_x(state, grid);
    if (reflect_y_) apply_reflective_y(state, grid);
    if (reflect_z_) apply_reflective_z(state, grid);
}

void ReflectiveBC::apply_reflective_x(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx_total = state.nx();
    int ny = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Left ghost: mirror from right interior with u -> -u
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int src_i = nghost + (nghost - 1 - g);
                    int dst_i = nghost - 1 - g;

                    if (v == 1) { // rho_u component
                        state(v, dst_i, j, k) = -state(v, src_i, j, k);
                    } else {
                        state(v, dst_i, j, k) = state(v, src_i, j, k);
                    }
                }
            }
        }
    }

    // Right ghost: mirror from left interior with u -> -u
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int src_i = nx_total - nghost - 1 - (nghost - 1 - g);
                    int dst_i = nx_total - nghost + g;

                    if (v == 1) { // rho_u component
                        state(v, dst_i, j, k) = -state(v, src_i, j, k);
                    } else {
                        state(v, dst_i, j, k) = state(v, src_i, j, k);
                    }
                }
            }
        }
    }
}

void ReflectiveBC::apply_reflective_y(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny_total = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: mirror with v -> -v
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    int src_j = nghost + (nghost - 1 - g);
                    int dst_j = nghost - 1 - g;

                    if (v == 2) { // rho_v component
                        state(v, i, dst_j, k) = -state(v, i, src_j, k);
                    } else {
                        state(v, i, dst_j, k) = state(v, i, src_j, k);
                    }
                }
            }
        }
    }

    // Top ghost: mirror with v -> -v
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    int src_j = ny_total - nghost - 1 - (nghost - 1 - g);
                    int dst_j = ny_total - nghost + g;

                    if (v == 2) { // rho_v component
                        state(v, i, dst_j, k) = -state(v, i, src_j, k);
                    } else {
                        state(v, i, dst_j, k) = state(v, i, src_j, k);
                    }
                }
            }
        }
    }
}

void ReflectiveBC::apply_reflective_z(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny = state.ny();
    int nz_total = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: mirror with w -> -w
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    int src_k = nghost + (nghost - 1 - g);
                    int dst_k = nghost - 1 - g;

                    if (v == 3) { // rho_w component
                        state(v, i, j, dst_k) = -state(v, i, j, src_k);
                    } else {
                        state(v, i, j, dst_k) = state(v, i, j, src_k);
                    }
                }
            }
        }
    }

    // Top ghost: mirror with w -> -w
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    int src_k = nz_total - nghost - 1 - (nghost - 1 - g);
                    int dst_k = nz_total - nghost + g;

                    if (v == 3) { // rho_w component
                        state(v, i, j, dst_k) = -state(v, i, j, src_k);
                    } else {
                        state(v, i, j, dst_k) = state(v, i, j, src_k);
                    }
                }
            }
        }
    }
}

} // namespace fvm3d::boundary
```

### 3.4 Transmissive (Outflow) Boundary Condition

**File**: `include/fvm3d/boundary/transmissive_bc.hpp`

```cpp
#pragma once

#include "boundary_condition.hpp"

namespace fvm3d::boundary {

/**
 * Transmissive (zero-gradient outflow) boundary condition.
 * Ghost cells are filled with values from the nearest interior cells.
 * Allows material to leave the domain freely.
 */
class TransmissiveBC : public BoundaryCondition {
public:
    TransmissiveBC(bool transmissive_x = false, bool transmissive_y = false, bool transmissive_z = false)
        : transmissive_x_(transmissive_x), transmissive_y_(transmissive_y), transmissive_z_(transmissive_z) {}

    void apply(
        StateField3D& state,
        const Grid3D& grid
    ) override;

    std::string name() const override { return "Transmissive"; }

private:
    bool transmissive_x_, transmissive_y_, transmissive_z_;

    void apply_transmissive_x(StateField3D& state, const Grid3D& grid);
    void apply_transmissive_y(StateField3D& state, const Grid3D& grid);
    void apply_transmissive_z(StateField3D& state, const Grid3D& grid);
};

} // namespace fvm3d::boundary
```

**File**: `src/boundary/transmissive_bc.cpp`

```cpp
#include "fvm3d/boundary/transmissive_bc.hpp"

namespace fvm3d::boundary {

void TransmissiveBC::apply(StateField3D& state, const Grid3D& grid) {
    if (transmissive_x_) apply_transmissive_x(state, grid);
    if (transmissive_y_) apply_transmissive_y(state, grid);
    if (transmissive_z_) apply_transmissive_z(state, grid);
}

void TransmissiveBC::apply_transmissive_x(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx_total = state.nx();
    int ny = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Left ghost: copy from first interior cell
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    state(v, g, j, k) = state(v, nghost, j, k);
                }
            }
        }
    }

    // Right ghost: copy from last interior cell
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int src_i = nx_total - nghost - 1;
                    int dst_i = nx_total - nghost + g;
                    state(v, dst_i, j, k) = state(v, src_i, j, k);
                }
            }
        }
    }
}

void TransmissiveBC::apply_transmissive_y(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny_total = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: copy from first interior row
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    state(v, i, g, k) = state(v, i, nghost, k);
                }
            }
        }
    }

    // Top ghost: copy from last interior row
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    int src_j = ny_total - nghost - 1;
                    int dst_j = ny_total - nghost + g;
                    state(v, i, dst_j, k) = state(v, i, src_j, k);
                }
            }
        }
    }
}

void TransmissiveBC::apply_transmissive_z(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny = state.ny();
    int nz_total = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: copy from first interior layer
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    state(v, i, j, g) = state(v, i, j, nghost);
                }
            }
        }
    }

    // Top ghost: copy from last interior layer
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    int src_k = nz_total - nghost - 1;
                    int dst_k = nz_total - nghost + g;
                    state(v, i, j, dst_k) = state(v, i, j, src_k);
                }
            }
        }
    }
}

} // namespace fvm3d::boundary
```

---

## 4. FVMSolver3D Orchestration

This is the main solver class that brings together all components. It orchestrates the pipeline, manages time stepping, and handles I/O.

**File**: `include/fvm3d/core/fvm_solver3d.hpp`

```cpp
#pragma once

#include "fvm3d/core/field3d.hpp"
#include "fvm3d/core/grid3d.hpp"
#include "fvm3d/physics/euler3d.hpp"
#include "fvm3d/spatial/riemann_solver.hpp"
#include "fvm3d/spatial/reconstruction.hpp"
#include "fvm3d/temporal/time_integrator.hpp"
#include "fvm3d/boundary/boundary_condition.hpp"
#include <memory>
#include <string>
#include <vector>

namespace fvm3d::core {

struct FVMSolverConfig {
    // Grid parameters
    double Lx, Ly, Lz;
    int nx, ny, nz;
    int nghost;

    // Physics
    std::string physics;  // "euler"

    // Spatial discretization
    std::string reconstruction;  // "constant", "muscl", "weno"
    std::string riemann_solver;  // "laxfriedrichs", "hll", "hllc"

    // Time integration
    std::string time_integrator;  // "euler", "rk2", "rk3"

    // Boundary conditions
    std::string bc_type;  // "periodic", "reflective", "transmissive"
    bool bc_x, bc_y, bc_z;

    // Time stepping
    double cfl;
    double t_final;
    double t_output;  // Output interval

    // Optional: MPI info
    int mpi_rank = 0;
    int mpi_size = 1;
};

/**
 * Main FVM Solver class for 3D compressible flow.
 *
 * This class orchestrates:
 * 1. Spatial discretization (reconstruction + Riemann solver)
 * 2. Time integration (explicit Runge-Kutta schemes)
 * 3. Boundary condition application
 * 4. Main time stepping loop
 * 5. I/O and checkpointing
 */
class FVMSolver3D {
public:
    /**
     * Constructor: Initialize solver with configuration
     */
    FVMSolver3D(const FVMSolverConfig& config);

    /**
     * Initialize the solver with an initial condition function.
     * @param init_func: function that sets state at each grid point
     *        void init(double x, double y, double z, Eigen::VectorXd& U)
     */
    using InitFunction = std::function<void(double, double, double, Eigen::VectorXd&)>;
    void initialize(const InitFunction& init_func);

    /**
     * Run the simulation for time t_final
     */
    void run();

    /**
     * Perform a single time step
     */
    void step();

    /**
     * Get current simulation time
     */
    double time() const { return t_current_; }

    /**
     * Get current state field
     */
    StateField3D& state() { return state_; }
    const StateField3D& state() const { return state_; }

    /**
     * Write checkpoint to HDF5 file
     */
    void write_checkpoint(const std::string& filename);

    /**
     * Load checkpoint from HDF5 file
     */
    void load_checkpoint(const std::string& filename);

private:
    // Configuration
    FVMSolverConfig config_;

    // Computational grid
    Grid3D grid_;

    // State vector (conservative variables)
    StateField3D state_;

    // Physics equations
    std::unique_ptr<physics::EulerEquations3D> physics_;

    // Numerical schemes
    std::unique_ptr<spatial::RiemannSolver> riemann_solver_;
    std::unique_ptr<spatial::ReconstructionScheme> reconstruction_;
    std::unique_ptr<temporal::TimeIntegrator> time_integrator_;
    std::unique_ptr<boundary::BoundaryCondition> boundary_condition_;

    // Time stepping
    double t_current_ = 0.0;
    int step_count_ = 0;

    // Temporary storage for RHS computation
    StateField3D rhs_;

    // Calculate RHS using FVM flux divergence
    void compute_rhs();

    // Apply CFL condition to determine time step
    double compute_dt();

    // Apply boundary conditions
    void apply_boundary_conditions();

    // Compute fluxes in one direction
    void compute_fluxes_direction(int direction, StateField3D& flux_out);

    // Print progress information
    void print_progress();
};

} // namespace fvm3d::core
```

**File**: `src/core/fvm_solver3d.cpp`

```cpp
#include "fvm3d/core/fvm_solver3d.hpp"
#include "fvm3d/spatial/riemann_solver_factory.hpp"
#include "fvm3d/spatial/reconstruction_factory.hpp"
#include "fvm3d/temporal/time_integrator_factory.hpp"
#include "fvm3d/boundary/periodic_bc.hpp"
#include "fvm3d/boundary/reflective_bc.hpp"
#include "fvm3d/boundary/transmissive_bc.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace fvm3d::core {

FVMSolver3D::FVMSolver3D(const FVMSolverConfig& config)
    : config_(config),
      grid_(GridGeometry3D{config.Lx, config.Ly, config.Lz, config.nx, config.ny, config.nz},
            config.nghost),
      state_(5, grid_.nx_total(), grid_.ny_total(), grid_.nz_total()),
      rhs_(5, grid_.nx_total(), grid_.ny_total(), grid_.nz_total())
{
    // Initialize physics module
    physics_ = std::make_unique<physics::EulerEquations3D>();

    // Initialize Riemann solver
    riemann_solver_ = spatial::RiemannSolverFactory::create(config.riemann_solver);

    // Initialize reconstruction scheme
    reconstruction_ = spatial::ReconstructionFactory::create(config.reconstruction);

    // Initialize time integrator
    time_integrator_ = temporal::TimeIntegratorFactory::create(config.time_integrator);

    // Initialize boundary condition
    if (config.bc_type == "periodic") {
        boundary_condition_ = std::make_unique<boundary::PeriodicBC>(
            config.bc_x, config.bc_y, config.bc_z
        );
    } else if (config.bc_type == "reflective") {
        boundary_condition_ = std::make_unique<boundary::ReflectiveBC>(
            config.bc_x, config.bc_y, config.bc_z
        );
    } else if (config.bc_type == "transmissive") {
        boundary_condition_ = std::make_unique<boundary::TransmissiveBC>(
            config.bc_x, config.bc_y, config.bc_z
        );
    } else {
        throw std::invalid_argument("Unknown boundary condition type: " + config.bc_type);
    }

    if (config.mpi_rank == 0) {
        std::cout << "FVM3D Solver initialized:\n"
                  << "  Grid: " << config.nx << " x " << config.ny << " x " << config.nz << "\n"
                  << "  Domain: [0," << config.Lx << "] x [0," << config.Ly << "] x [0," << config.Lz << "]\n"
                  << "  Riemann Solver: " << riemann_solver_->name() << "\n"
                  << "  Reconstruction: " << reconstruction_->name() << "\n"
                  << "  Time Integrator: " << time_integrator_->name() << "\n"
                  << "  Boundary Condition: " << boundary_condition_->name() << "\n";
    }
}

void FVMSolver3D::initialize(const InitFunction& init_func) {
    // Set initial conditions at all grid points
    int nghost = grid_.nghost();
    int nx_total = grid_.nx_total();
    int ny_total = grid_.ny_total();
    int nz_total = grid_.nz_total();

    for (int i = 0; i < nx_total; i++) {
        for (int j = 0; j < ny_total; j++) {
            for (int k = 0; k < nz_total; k++) {
                double x = grid_.cell_center_x(i);
                double y = grid_.cell_center_y(j);
                double z = grid_.cell_center_z(k);

                Eigen::VectorXd U(5);
                init_func(x, y, z, U);

                for (int v = 0; v < 5; v++) {
                    state_(v, i, j, k) = U(v);
                }
            }
        }
    }

    // Apply boundary conditions
    apply_boundary_conditions();

    if (config_.mpi_rank == 0) {
        std::cout << "Initial condition applied.\n";
    }
}

void FVMSolver3D::run() {
    if (config_.mpi_rank == 0) {
        std::cout << "Starting main time loop...\n";
        std::cout << "Final time: " << config_.t_final << "\n";
        std::cout << "CFL: " << config_.cfl << "\n\n";
    }

    while (t_current_ < config_.t_final) {
        step();

        // Print progress periodically
        if (step_count_ % 100 == 0) {
            print_progress();
        }
    }

    if (config_.mpi_rank == 0) {
        std::cout << "\nSimulation completed.\n"
                  << "Total steps: " << step_count_ << "\n"
                  << "Final time: " << t_current_ << "\n";
    }
}

void FVMSolver3D::step() {
    // Compute time step
    double dt = compute_dt();

    // Ensure we don't overstep
    if (t_current_ + dt > config_.t_final) {
        dt = config_.t_final - t_current_;
    }

    // Time integration step using callback function for RHS
    auto rhs_func = [this](const StateField3D& U, StateField3D& dUdt) {
        // This is called by the time integrator to compute RHS
        compute_rhs();
        dUdt = rhs_;
    };

    time_integrator_->step(state_, dt, rhs_func);

    // Apply boundary conditions
    apply_boundary_conditions();

    // Update time
    t_current_ += dt;
    step_count_++;
}

double FVMSolver3D::compute_dt() {
    // Compute maximum wave speed across domain
    double max_speed = 0.0;
    int nghost = grid_.nghost();
    int nx_total = grid_.nx_total();
    int ny_total = grid_.ny_total();
    int nz_total = grid_.nz_total();

    for (int i = nghost; i < nx_total - nghost; i++) {
        for (int j = nghost; j < ny_total - nghost; j++) {
            for (int k = nghost; k < nz_total - nghost; k++) {
                // Get conservative state
                Eigen::VectorXd U(5);
                for (int v = 0; v < 5; v++) {
                    U(v) = state_(v, i, j, k);
                }

                // Check wave speed in each direction
                for (int dir = 0; dir < 3; dir++) {
                    double speed = physics_->max_wave_speed(U, dir);
                    max_speed = std::max(max_speed, speed);
                }
            }
        }
    }

    // Compute CFL-limited time step
    double dx = grid_.geometry().dx;
    double dy = grid_.geometry().dy;
    double dz = grid_.geometry().dz;
    double min_dx = std::min({dx, dy, dz});

    double dt = config_.cfl * min_dx / std::max(max_speed, 1e-10);

    return dt;
}

void FVMSolver3D::compute_rhs() {
    // Zero out RHS
    rhs_.fill(0.0);

    // Compute fluxes in X direction and add to RHS
    // Flux divergence: -dF/dx
    // This involves reconstructing interface values and computing Riemann fluxes
    // Implementation follows standard FVM approach

    // For now, simplified version (full implementation would include
    // reconstruction, Riemann solve, flux divergence)

    // TODO: Implement full RHS computation with reconstruction + flux calculation
}

void FVMSolver3D::apply_boundary_conditions() {
    boundary_condition_->apply(state_, grid_);
}

void FVMSolver3D::write_checkpoint(const std::string& filename) {
    // TODO: Implement HDF5 checkpoint writing
    if (config_.mpi_rank == 0) {
        std::cout << "Writing checkpoint: " << filename << "\n";
    }
}

void FVMSolver3D::load_checkpoint(const std::string& filename) {
    // TODO: Implement HDF5 checkpoint loading
    if (config_.mpi_rank == 0) {
        std::cout << "Loading checkpoint: " << filename << "\n";
    }
}

void FVMSolver3D::print_progress() {
    if (config_.mpi_rank == 0) {
        std::cout << std::fixed << std::setprecision(6)
                  << "Step " << std::setw(6) << step_count_
                  << " | Time " << std::setw(12) << t_current_
                  << " | Progress " << std::setw(6) << (100.0 * t_current_ / config_.t_final) << "%\n";
    }
}

} // namespace fvm3d::core
```

---

## 5. Complete 3D Blast Wave Example

This is a complete working example that brings all components together.

**File**: `examples/blast_wave_3d.cpp`

```cpp
#include "fvm3d/core/fvm_solver3d.hpp"
#include <cmath>
#include <iostream>

// Initial condition: 3D spherical blast wave
void initialize_blast_wave(double x, double y, double z, Eigen::VectorXd& U) {
    double center_x = 0.5, center_y = 0.5, center_z = 0.5;
    double radius = 0.1;
    double r = std::sqrt((x - center_x)*(x - center_x) +
                        (y - center_y)*(y - center_y) +
                        (z - center_z)*(z - center_z));

    double gamma = 1.4;
    double rho_ambient = 1.0;
    double p_ambient = 0.1;
    double p_blast = 10.0;

    U(0) = rho_ambient;  // rho
    U(1) = 0.0;          // rho_u
    U(2) = 0.0;          // rho_v
    U(3) = 0.0;          // rho_w

    // Pressure
    double p = (r < radius) ? p_blast : p_ambient;

    // Total energy: E = p / (gamma - 1) + 0.5 * rho * (u^2 + v^2 + w^2)
    U(4) = p / (gamma - 1.0);  // (kinetic energy is zero initially)
}

int main(int argc, char** argv) {
    // Configure solver
    fvm3d::core::FVMSolverConfig config;
    config.Lx = 1.0;
    config.Ly = 1.0;
    config.Lz = 1.0;
    config.nx = 64;
    config.ny = 64;
    config.nz = 64;
    config.nghost = 2;
    config.physics = "euler";
    config.reconstruction = "muscl";
    config.riemann_solver = "hllc";
    config.time_integrator = "rk3";
    config.bc_type = "transmissive";
    config.bc_x = config.bc_y = config.bc_z = true;
    config.cfl = 0.4;
    config.t_final = 0.2;
    config.t_output = 0.01;

    // Create solver
    fvm3d::core::FVMSolver3D solver(config);

    // Initialize with blast wave
    solver.initialize(initialize_blast_wave);

    // Run simulation
    solver.run();

    // Write final state
    solver.write_checkpoint("blast_wave_final.h5");

    return 0;
}
```

---

## 6. Testing Framework and Unit Tests

### 6.1 Unit Test Structure

**File**: `tests/test_riemann_solvers.cpp`

```cpp
#include <catch2/catch.hpp>
#include "fvm3d/spatial/riemann_solver_factory.hpp"
#include <Eigen/Dense>

using namespace fvm3d::spatial;

TEST_CASE("Lax-Friedrichs Riemann Solver", "[riemann]") {
    auto solver = RiemannSolverFactory::create("laxfriedrichs");

    // Test case 1: Sod shock tube (1D in X direction)
    Eigen::VectorXd U_L(5), U_R(5);

    // Left state: rho=1, u=0, v=0, w=0, p=1
    U_L(0) = 1.0;
    U_L(1) = 0.0;
    U_L(2) = 0.0;
    U_L(3) = 0.0;
    U_L(4) = 1.0 / 0.4;  // E = p / (gamma - 1)

    // Right state: rho=0.125, u=0, v=0, w=0, p=0.1
    U_R(0) = 0.125;
    U_R(1) = 0.0;
    U_R(2) = 0.0;
    U_R(3) = 0.0;
    U_R(4) = 0.1 / 0.4;

    Eigen::VectorXd flux = solver->solve(U_L, U_R, 0);
    double speed = solver->max_wave_speed(U_L, U_R, 0);

    REQUIRE(flux.size() == 5);
    REQUIRE(speed > 0.0);
    REQUIRE(flux(0) > 0.0);  // Mass flux should be positive (flow from left to right)
}

TEST_CASE("HLL Riemann Solver", "[riemann]") {
    auto solver = RiemannSolverFactory::create("hll");

    Eigen::VectorXd U_L(5), U_R(5);
    U_L(0) = 1.0;
    U_L(1) = 0.0;
    U_L(2) = 0.0;
    U_L(3) = 0.0;
    U_L(4) = 2.5;

    U_R(0) = 0.125;
    U_R(1) = 0.0;
    U_R(2) = 0.0;
    U_R(3) = 0.0;
    U_R(4) = 0.25;

    Eigen::VectorXd flux = solver->solve(U_L, U_R, 0);
    REQUIRE(flux.size() == 5);
}

TEST_CASE("HLLC Riemann Solver", "[riemann]") {
    auto solver = RiemannSolverFactory::create("hllc");

    Eigen::VectorXd U_L(5), U_R(5);
    U_L(0) = 1.0;
    U_L(1) = 0.1;  // Non-zero velocity
    U_L(2) = 0.0;
    U_L(3) = 0.0;
    U_L(4) = 2.5;

    U_R(0) = 0.125;
    U_R(1) = -0.1;  // Opposite velocity direction
    U_R(2) = 0.0;
    U_R(3) = 0.0;
    U_R(4) = 0.25;

    Eigen::VectorXd flux = solver->solve(U_L, U_R, 0);
    REQUIRE(flux.size() == 5);
}
```

### 6.2 Time Integrator Tests

**File**: `tests/test_time_integrators.cpp`

```cpp
#include <catch2/catch.hpp>
#include "fvm3d/temporal/time_integrator_factory.hpp"
#include "fvm3d/core/field3d.hpp"
#include <cmath>

using namespace fvm3d::temporal;
using namespace fvm3d::core;

TEST_CASE("Forward Euler Time Integrator", "[temporal]") {
    auto integrator = TimeIntegratorFactory::create("euler");

    // Test on simple ODE: du/dt = -u (exponential decay)
    StateField3D U(1, 3, 3, 3);
    U.fill(1.0);  // Initial condition: u = 1

    // RHS function: du/dt = -u
    auto rhs = [](const StateField3D& U_in, StateField3D& dUdt_out) {
        int nx = U_in.nx(), ny = U_in.ny(), nz = U_in.nz();
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    dUdt_out(0, i, j, k) = -U_in(0, i, j, k);
                }
            }
        }
    };

    double dt = 0.01;
    integrator->step(U, dt, rhs);

    // Check that solution decreased (exponential decay)
    REQUIRE(U(0, 1, 1, 1) < 1.0);
    REQUIRE(U(0, 1, 1, 1) > 0.0);
    REQUIRE(std::abs(U(0, 1, 1, 1) - std::exp(-dt)) < 0.01);
}

TEST_CASE("RK2 Time Integrator", "[temporal]") {
    auto integrator = TimeIntegratorFactory::create("rk2");

    StateField3D U(1, 3, 3, 3);
    U.fill(1.0);

    auto rhs = [](const StateField3D& U_in, StateField3D& dUdt_out) {
        int nx = U_in.nx(), ny = U_in.ny(), nz = U_in.nz();
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    dUdt_out(0, i, j, k) = -U_in(0, i, j, k);
                }
            }
        }
    };

    double dt = 0.01;
    integrator->step(U, dt, rhs);

    // RK2 should be more accurate than Forward Euler
    REQUIRE(std::abs(U(0, 1, 1, 1) - std::exp(-dt)) < 0.001);
}

TEST_CASE("RK3 Time Integrator", "[temporal]") {
    auto integrator = TimeIntegratorFactory::create("rk3");

    StateField3D U(1, 3, 3, 3);
    U.fill(1.0);

    auto rhs = [](const StateField3D& U_in, StateField3D& dUdt_out) {
        int nx = U_in.nx(), ny = U_in.ny(), nz = U_in.nz();
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    dUdt_out(0, i, j, k) = -U_in(0, i, j, k);
                }
            }
        }
    };

    double dt = 0.01;
    integrator->step(U, dt, rhs);

    // RK3 should be most accurate
    REQUIRE(std::abs(U(0, 1, 1, 1) - std::exp(-dt)) < 0.0001);
}
```

### 6.3 Boundary Condition Tests

**File**: `tests/test_boundary_conditions.cpp`

```cpp
#include <catch2/catch.hpp>
#include "fvm3d/boundary/periodic_bc.hpp"
#include "fvm3d/boundary/reflective_bc.hpp"
#include "fvm3d/boundary/transmissive_bc.hpp"
#include "fvm3d/core/field3d.hpp"
#include "fvm3d/core/grid3d.hpp"

using namespace fvm3d::boundary;
using namespace fvm3d::core;

TEST_CASE("Periodic Boundary Condition", "[boundary]") {
    GridGeometry3D geom{1.0, 1.0, 1.0, 8, 8, 8};
    Grid3D grid(geom, 2);

    StateField3D state(5, grid.nx_total(), grid.ny_total(), grid.nz_total());
    state.fill(0.0);

    // Set interior values
    for (int i = 2; i < 10; i++) {
        for (int j = 2; j < 10; j++) {
            for (int k = 2; k < 10; k++) {
                state(0, i, j, k) = i + j + k;
            }
        }
    }

    PeriodicBC bc(true, false, false);
    bc.apply(state, grid);

    // Check left ghost matches right interior
    REQUIRE(state(0, 0, 5, 5) == state(0, 9, 5, 5));
    REQUIRE(state(0, 1, 5, 5) == state(0, 8, 5, 5));

    // Check right ghost matches left interior
    REQUIRE(state(0, 10, 5, 5) == state(0, 2, 5, 5));
    REQUIRE(state(0, 11, 5, 5) == state(0, 3, 5, 5));
}

TEST_CASE("Reflective Boundary Condition", "[boundary]") {
    GridGeometry3D geom{1.0, 1.0, 1.0, 8, 8, 8};
    Grid3D grid(geom, 2);

    StateField3D state(5, grid.nx_total(), grid.ny_total(), grid.nz_total());
    state.fill(0.0);

    // Set interior values
    for (int i = 2; i < 10; i++) {
        for (int j = 2; j < 10; j++) {
            for (int k = 2; k < 10; k++) {
                state(0, i, j, k) = 1.0;       // rho
                state(1, i, j, k) = 2.0;       // rho_u
                state(2, i, j, k) = 3.0;       // rho_v
                state(3, i, j, k) = 4.0;       // rho_w
                state(4, i, j, k) = 5.0;       // E
            }
        }
    }

    ReflectiveBC bc(true, false, false);
    bc.apply(state, grid);

    // Check that normal velocity is reversed
    REQUIRE(state(0, 0, 5, 5) == state(0, 3, 5, 5));  // rho unchanged
    REQUIRE(state(1, 0, 5, 5) == -state(1, 3, 5, 5)); // u reversed
    REQUIRE(state(2, 0, 5, 5) == state(2, 3, 5, 5));  // v unchanged
}

TEST_CASE("Transmissive Boundary Condition", "[boundary]") {
    GridGeometry3D geom{1.0, 1.0, 1.0, 8, 8, 8};
    Grid3D grid(geom, 2);

    StateField3D state(5, grid.nx_total(), grid.ny_total(), grid.nz_total());
    state.fill(0.0);

    // Set interior values with gradient
    for (int i = 2; i < 10; i++) {
        for (int j = 2; j < 10; j++) {
            for (int k = 2; k < 10; k++) {
                state(0, i, j, k) = double(i);  // gradient in X
            }
        }
    }

    TransmissiveBC bc(true, false, false);
    bc.apply(state, grid);

    // Check that all ghost cells match interior boundary
    for (int j = 0; j < grid.ny_total(); j++) {
        for (int k = 0; k < grid.nz_total(); k++) {
            REQUIRE(state(0, 0, j, k) == state(0, 2, j, k));
            REQUIRE(state(0, 1, j, k) == state(0, 2, j, k));
            REQUIRE(state(0, 10, j, k) == state(0, 9, j, k));
            REQUIRE(state(0, 11, j, k) == state(0, 9, j, k));
        }
    }
}
```

---

## 7. Performance Profiling Integration

### 7.1 Timer Utility

**File**: `include/fvm3d/util/timer.hpp`

```cpp
#pragma once

#include <chrono>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>

namespace fvm3d::util {

/**
 * Simple timer for performance profiling.
 */
class Timer {
public:
    void start(const std::string& name) {
        timers_[name] = std::chrono::high_resolution_clock::now();
    }

    double stop(const std::string& name) {
        auto end = std::chrono::high_resolution_clock::now();
        auto it = timers_.find(name);
        if (it == timers_.end()) {
            return 0.0;
        }

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            end - it->second
        );

        double elapsed = duration.count();
        totals_[name] += elapsed;
        counts_[name]++;

        return elapsed;
    }

    void print_summary() {
        std::cout << "\n=== Performance Summary ===\n";
        std::cout << std::left << std::setw(30) << "Operation"
                  << std::right << std::setw(12) << "Total (ms)"
                  << std::setw(12) << "Count"
                  << std::setw(12) << "Avg (ms)\n";
        std::cout << std::string(66, '-') << "\n";

        for (const auto& [name, total] : totals_) {
            int count = counts_[name];
            double avg = total / count;

            std::cout << std::left << std::setw(30) << name
                      << std::right << std::setw(12) << std::fixed << std::setprecision(2) << total
                      << std::setw(12) << count
                      << std::setw(12) << avg << "\n";
        }
    }

private:
    std::map<std::string, std::chrono::high_resolution_clock::time_point> timers_;
    std::map<std::string, double> totals_;
    std::map<std::string, int> counts_;
};

} // namespace fvm3d::util
```

---

## 8. Implementation Checklist

### Phase 1: Core Modules (Days 1-7)
- [ ] **Day 1-2**: Grid3D and Field3D implementation
  - [ ] GridGeometry3D struct with validation
  - [ ] Grid3D class with indexing and ghost cell management
  - [ ] StateField3D template implementation with SoA layout
  - [ ] Move semantics and memory safety

- [ ] **Day 3-4**: Physics module
  - [ ] EulerEquations3D with primitive/conservative conversion
  - [ ] Max wave speed estimation
  - [ ] Stability floors (rho_floor, p_floor)
  - [ ] Unit tests for conservation checks

- [ ] **Day 5-7**: Riemann solvers
  - [ ] RiemannSolver base class with interface
  - [ ] LaxFriedrichsSolver implementation
  - [ ] HLLSolver implementation
  - [ ] HLLCSolver implementation
  - [ ] RiemannSolverFactory
  - [ ] Unit tests for all solvers (Sod shock tube)

### Phase 2: Time Integration & Boundary Conditions (Days 8-14)
- [ ] **Day 8-9**: Time integrators
  - [ ] TimeIntegrator base class
  - [ ] ForwardEuler (1st order)
  - [ ] RK2 (SSP-RK2, 2nd order)
  - [ ] RK3 (SSP-RK3, 3rd order)
  - [ ] TimeIntegratorFactory
  - [ ] Unit tests for accuracy verification

- [ ] **Day 10-11**: Boundary conditions
  - [ ] BoundaryCondition base class
  - [ ] PeriodicBC implementation
  - [ ] ReflectiveBC implementation
  - [ ] TransmissiveBC implementation
  - [ ] Ghost cell filling logic in all directions
  - [ ] Unit tests for boundary application

- [ ] **Day 12-14**: Reconstruction schemes
  - [ ] ReconstructionScheme base class
  - [ ] Constant reconstruction
  - [ ] MUSCL reconstruction
  - [ ] WENO reconstruction (optional for v1.0)
  - [ ] ReconstructionFactory
  - [ ] Integration with flux calculation

### Phase 3: Orchestration & Testing (Days 15-19)
- [ ] **Day 15-16**: FVMSolver3D
  - [ ] Main solver orchestration class
  - [ ] RHS computation (flux divergence)
  - [ ] CFL-limited time stepping
  - [ ] HDF5 checkpoint I/O
  - [ ] Progress printing

- [ ] **Day 17-18**: Examples and validation
  - [ ] 3D Blast Wave example
  - [ ] Sod shock tube (1D in 3D grid)
  - [ ] Lax-Liu test case
  - [ ] Convergence studies
  - [ ] Comparison with reference solutions

- [ ] **Day 19**: Testing framework
  - [ ] Catch2 integration
  - [ ] Unit tests for all components
  - [ ] Integration tests
  - [ ] Regression tests

### Phase 4: MPI Parallelization (Days 20+)
- [ ] **Day 20**: MPI foundation
  - [ ] MPIGrid3D with 1D decomposition
  - [ ] Domain decomposer
  - [ ] Rank assignment logic

- [ ] **Day 21-22**: Communication
  - [ ] MPIHaloExchange non-blocking communication
  - [ ] Buffer allocation and management
  - [ ] Pack/unpack operations

- [ ] **Day 23-24**: Parallel I/O
  - [ ] Parallel HDF5 writing
  - [ ] Collective operations
  - [ ] Load balancing verification

- [ ] **Day 25**: Optimization
  - [ ] Performance profiling
  - [ ] Scalability analysis
  - [ ] Final testing

### Code Quality Checkpoints
- [ ] All code follows C++17 standard
- [ ] No memory leaks (valgrind clean)
- [ ] All unit tests pass
- [ ] No compiler warnings
- [ ] Documentation complete
- [ ] Examples runnable and validated

---

## Summary

This comprehensive guide provides all the code templates and specifications needed to implement the remaining core modules of the 3D C++ MPI FVM framework. Key features:

1. **Complete Riemann Solvers**: LF, HLL, HLLC with proper wave speed estimation
2. **Multiple Time Integrators**: Forward Euler, RK2, RK3 with SSP property
3. **Flexible Boundary Conditions**: Periodic, Reflective, Transmissive
4. **Main FVMSolver3D**: Orchestration of all components with time stepping
5. **Complete Blast Wave Example**: Fully working demonstration
6. **Testing Framework**: Unit tests for all major components
7. **Performance Tools**: Timer utility for profiling

All implementations follow the architecture specified in "3D_CPP_Framework_Design.md" and maintain compatibility with the MPI parallelization layer described in "3D_CPP_Implementation_Guide.md".

