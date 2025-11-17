# Flux Calculation and Reconstruction Architecture

## Overview

The FVM 3D solver now features a modular architecture separating flux calculation and spatial reconstruction into independent, interchangeable components.

## Design Philosophy

### Separation of Concerns

1. **Reconstruction Methods**: Compute left/right interface states from cell-centered values
2. **Flux Calculators**: Compute numerical fluxes from interface states
3. **Physics Equations**: Define physical flux functions and variable transformations

```
┌─────────────────┐
│ Cell-Centered   │
│    Values       │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Reconstruction  │  ← Constant, MUSCL, WENO, ENO
│    Method       │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Interface      │
│    States       │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│     Flux        │  ← Riemann solvers, Central schemes
│  Calculator     │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Numerical      │
│    Fluxes       │
└─────────────────┘
```

## Directory Structure

```
include/spatial/
├── flux_calculation/
│   ├── flux_calculator_base.hpp      # Base class for all flux calculators
│   ├── riemann_flux_adapter.hpp      # Adapter for existing Riemann solvers
│   └── flux_calculator_factory.hpp   # Factory to create flux calculators
│
├── reconstruction/
│   ├── reconstruction_base.hpp       # Base class for reconstruction methods
│   └── reconstruction_factory.hpp    # Factory (future)
│
└── [existing files...]
    ├── riemann_solver.hpp            # Existing Riemann solver interface
    ├── riemann_laxfriedrichs.hpp
    ├── riemann_hll.hpp
    ├── riemann_hllc.hpp
    └── riemann_hlld.hpp
```

## Class Hierarchy

### Flux Calculators

```
FluxCalculator (abstract base)
├── RiemannSolverFlux
│   └── RiemannFluxAdapter (wraps existing RiemannSolver)
│       ├── Lax-Friedrichs
│       ├── HLL
│       ├── HLLC
│       └── HLLD
├── DissipativeFluxCalculator
│   └── (Future: Rusanov, local Lax-Friedrichs)
└── CentralFluxCalculator
    └── (Future: central differencing, spectral methods)
```

### Reconstruction Methods

```
ReconstructionMethod (abstract base)
├── ConstantReconstruction (1st order, piecewise constant)
├── LimitedReconstructionMethod
│   ├── MUSCLReconstruction (2nd order with limiters)
│   │   ├── minmod limiter
│   │   ├── van Leer limiter
│   │   ├── superbee limiter
│   │   └── MC limiter
│   └── (Future: ENO, WENO)
└── (Future: High-order methods)
```

## Key Interfaces

### FluxCalculator

```cpp
class FluxCalculator {
public:
    // Core method: compute flux from left/right states
    virtual Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const physics::PhysicsBase& physics,
        int direction
    ) const = 0;

    // Compute max wave speed for CFL condition
    virtual double compute_max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const physics::PhysicsBase& physics,
        int direction
    ) const = 0;

    // Query methods
    virtual bool is_physics_agnostic() const;
    virtual bool supports_physics(const std::string& type) const;
    virtual bool needs_wave_speed() const;
};
```

### ReconstructionMethod

```cpp
class ReconstructionMethod {
public:
    // Core method: reconstruct interface states
    virtual void reconstruct(
        const core::Field3D<double>& U,
        int i, int j, int k,
        int direction,
        Eigen::VectorXd& U_L,
        Eigen::VectorXd& U_R
    ) const = 0;

    // Batch processing (more efficient)
    virtual void reconstruct_all(
        const core::Field3D<double>& U,
        int direction,
        core::Field3D<double>& U_left,
        core::Field3D<double>& U_right
    ) const;

    // Query methods
    int order() const;
    bool is_tvd() const;
    bool uses_limiters() const;
    int required_ghost_cells() const;
};
```

## Usage Examples

### Creating Flux Calculators

```cpp
// Via factory (recommended)
auto flux_calc = FluxCalculatorFactory::create("hllc", "mhd_advanced", 9);

// Check capabilities
if (flux_calc->supports_physics("mhd_advanced")) {
    std::cout << "HLLC supports MHD!\n";
}

// Compute flux at an interface
Eigen::VectorXd flux = flux_calc->compute_flux(U_L, U_R, physics, direction);
double max_speed = flux_calc->compute_max_wave_speed(U_L, U_R, physics, direction);
```

### Creating Reconstruction Methods

```cpp
// Via factory (future)
auto recon = ReconstructionFactory::create("muscl", "minmod");

// Query properties
std::cout << "Order: " << recon->order() << "\n";
std::cout << "TVD: " << (recon->is_tvd() ? "yes" : "no") << "\n";
std::cout << "Ghost cells needed: " << recon->required_ghost_cells() << "\n";

// Reconstruct interface states
Eigen::VectorXd U_L, U_R;
recon->reconstruct(U, i, j, k, direction, U_L, U_R);
```

### Complete Simulation Loop

```cpp
// Setup
auto recon = ReconstructionFactory::create("muscl", "minmod");
auto flux_calc = FluxCalculatorFactory::create("hllc", "euler", 5);
auto physics = std::make_shared<physics::EulerEquations3D>();

// Main loop
for (int step = 0; step < num_steps; ++step) {
    // 1. Reconstruct interface states
    Field3D<double> U_left, U_right;
    recon->reconstruct_all(U_current, direction, U_left, U_right);

    // 2. Compute numerical fluxes
    Field3D<double> fluxes(nvars, nx+1, ny, nz);
    for (int i = 0; i < nx+1; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                Eigen::VectorXd U_L = get_state(U_left, i, j, k);
                Eigen::VectorXd U_R = get_state(U_right, i, j, k);

                auto flux = flux_calc->compute_flux(U_L, U_R, *physics, direction);
                set_state(fluxes, i, j, k, flux);
            }
        }
    }

    // 3. Update solution (time integration)
    time_integrator->step(U_current, dt, fluxes);
}
```

## Benefits

### 1. Modularity
- Mix and match reconstruction methods with flux calculators
- Easy to add new methods without modifying existing code

### 2. Extensibility
- Add WENO/ENO reconstruction without touching flux calculators
- Add AUSM/KFVS flux methods without touching reconstruction
- Add new physics equations - both components work unchanged

### 3. Testing
- Test reconstruction and flux calculation independently
- Unit test each component in isolation
- Easier to verify correctness

### 4. Performance
- Optimize reconstruction and flux calculation separately
- Profile to identify bottlenecks
- Easy to experiment with different combinations

### 5. Backwards Compatibility
- Existing Riemann solvers work via RiemannFluxAdapter
- No need to rewrite existing code
- Gradual migration path

## Migration Guide

### From Old RiemannSolver to New FluxCalculator

**Old code:**
```cpp
auto riemann = RiemannSolverFactory::create("hllc", "euler", 5);
auto flux = riemann->solve(U_L, U_R, direction);
```

**New code (option 1 - via adapter):**
```cpp
auto riemann = RiemannSolverFactory::create("hllc", "euler", 5);
auto flux_calc = std::make_unique<RiemannFluxAdapter>(std::move(riemann));
auto flux = flux_calc->compute_flux(U_L, U_R, physics, direction);
```

**New code (option 2 - via factory):**
```cpp
auto flux_calc = FluxCalculatorFactory::create("hllc", "euler", 5);
auto flux = flux_calc->compute_flux(U_L, U_R, physics, direction);
```

## Future Extensions

### Planned Flux Calculators
- **Central schemes**: No Riemann problem solving, faster for smooth flows
- **AUSM family**: Advection Upstream Splitting Methods
- **KFVS**: Kinetic Flux Vector Splitting
- **Roe-type**: Approximate Riemann solvers

### Planned Reconstruction Methods
- **WENO3/5**: Weighted Essentially Non-Oscillatory
- **ENO**: Essentially Non-Oscillatory
- **PPM**: Piecewise Parabolic Method
- **Compact schemes**: High-order on compact stencils

### Planned Features
- **Adaptive method selection**: Choose method based on local flow features
- **Hybrid schemes**: Mix methods in different regions
- **Characteristic-wise reconstruction**: Reconstruct in characteristic variables
- **Primitive variable reconstruction**: Option to reconstruct p,u,v,w instead of ρ,ρu,E

## References

1. LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*. Cambridge University Press.

2. Toro, E. F. (2009). *Riemann Solvers and Numerical Methods for Fluid Dynamics*. Springer.

3. Shu, C. W. (2009). "High order weighted essentially non-oscillatory schemes for convection dominated problems." *SIAM Review*, 51(1), 82-126.

4. van Leer, B. (1979). "Towards the ultimate conservative difference scheme. V. A second-order sequel to Godunov's method." *Journal of Computational Physics*, 32(1), 101-136.
