# FVM3D MHD Benchmark Tests

Comprehensive benchmark suite for validating the HLLD Riemann solver implementation for magnetohydrodynamics (MHD).

## Overview

These benchmarks provide well-established test cases from the MHD literature that validate:
- Correctness of the HLLD Riemann solver (5-wave structure)
- Proper handling of magnetic field evolution
- Stability and accuracy of numerical schemes
- Alfvén wave capturing capability

## Benchmarks

### 1. Brio-Wu Shock Tube Test

**Type:** 1D MHD Riemann problem
**Executable:** `build/brio_wu_shock_tube`

#### Description

The Brio-Wu test is a classic 1D MHD shock tube problem that validates HLLD solver performance on discontinuous initial conditions. The key feature is a magnetic field discontinuity (flip in tangential component).

#### Initial Conditions (Brio & Wu 1988)

**Left state (x < 0.5):**
- Density: ρ = 1.0
- Velocity: u = v = w = 0
- Pressure: p = 1.0
- Magnetic field: Bx = 0.75, By = 1.0, Bz = 0

**Right state (x > 0.5):**
- Density: ρ = 0.125
- Velocity: u = v = w = 0
- Pressure: p = 0.1
- Magnetic field: Bx = 0.75, By = -1.0, Bz = 0

#### Key Features

1. **5 MHD Waves:** The solution structure contains all characteristic waves:
   - Fast magnetosonic wave (left)
   - Slow magnetosonic wave (left)
   - **Alfvén wave (left)** ← Tests HLLD Alfvén handling
   - Contact discontinuity
   - Alfvén wave (right)
   - Slow magnetosonic wave (right)
   - Fast magnetosonic wave (right)

2. **Magnetic Field Flip:** By changes sign across discontinuity, testing proper magnetic field evolution

3. **Exact Solution:** Published analytical solution available for comparison (Brio & Wu 1988)

#### Running the Test

```bash
cd /tmp
mkdir -p brio_wu_test && cd brio_wu_test
/path/to/fvm_3d/build/brio_wu_shock_tube
```

#### Expected Results

- **Solution structure:** 5 distinct wave regions with proper density/pressure profiles
- **Alfvén waves:** Sharp discontinuities in tangential magnetic field (By)
- **Density range:** [0.125, 1.0]
- **Pressure range:** [0.1, 1.0]
- **Smooth regions:** Between fast and slow waves, density/pressure vary smoothly
- **Divergence cleaning:** GLM variable (psi) remains small (~0)

#### Validation Criteria

✓ Compare density and pressure profiles with published data
✓ Check for Alfvén wave discontinuities in By
✓ Verify energy conservation (within numerical error)
✓ Confirm positivity of density and pressure throughout

### 2. Orszag-Tang Vortex Test

**Type:** 2D smooth MHD problem with periodic boundary conditions
**Executable:** `build/orszag_tang_vortex`

#### Description

The Orszag-Tang vortex is a smooth 2D MHD problem that tests the solver's ability to handle:
- Smooth velocity and magnetic field structures
- Complex magnetic field interactions
- Potential development of turbulence-like structures
- Long-time stability

This problem is excellent for testing overall solution accuracy without shock-induced errors.

#### Initial Conditions (Orszag & Tang 1979)

**Domain:** [0, 2π] × [0, 2π] with **periodic boundary conditions**

**Primitive variables:**
```
Density:     ρ(x,y) = γ²
Velocity:    u(x,y) = -sin(y),  v(x,y) = sin(x),  w = 0
Magnetic:    Bx(x,y) = -sin(y), By(x,y) = sin(2x), Bz = 0
Pressure:    p(x,y) = (γ-1)/(γ·π²) · |B|²/2 + 1.0
```

where γ = 5/3 (adiabatic index)

#### Key Features

1. **Smooth Initial Conditions:** No discontinuities, tests accuracy on smooth problems

2. **Vortex Structure:** Symmetric sine wave patterns create vortex-like flows

3. **Magnetic Reconnection:** Crossing magnetic field lines create reconnection regions

4. **Turbulent Cascade:** Small-scale structures develop naturally from smooth initial conditions

5. **Periodic BC:** Tests periodic boundary condition handling

#### Running the Test

```bash
cd /tmp
mkdir -p orszag_tang_test && cd orszag_tang_test
/path/to/fvm_3d/build/orszag_tang_vortex
```

#### Expected Results

- **Smooth evolution:** Initial smooth structures remain stable
- **Structure development:** Small-scale features develop over time
- **Density positivity:** ρ > 0 throughout simulation
- **Pressure positivity:** p > 0 throughout simulation
- **Energy behavior:** Energy dissipates smoothly due to numerical diffusion
- **Divergence cleaning:** ∇·B remains small (GLM maintains constraint)

#### Validation Criteria

✓ Check density and pressure positivity throughout
✓ Compare with published benchmark results (e.g., Orszag & Tang 1979)
✓ Verify energy evolution (monotonic decrease expected)
✓ Examine divergence of magnetic field (∇·B ≈ 0)
✓ Look for expected small-scale structure development

## Configuration Details

### Grid Resolution

**Brio-Wu Shock Tube:**
- x-direction: 128 cells (1D shock tube)
- y, z-directions: 2 cells each (minimal 3D overhead)
- Domain: [0, 1] in x, [0, 0.01] in y,z

**Orszag-Tang Vortex:**
- x, y-directions: 96 cells each (2D resolution)
- z-direction: 2 cells (minimal 3D overhead)
- Domain: [0, 2π] in x,y; [0, 0.1] in z

### Solver Configuration

**Physics:** mhd_advanced (with GLM divergence cleaning)
- Resistivity model: eta0 = 1e-3, eta1 = 1.67e-2
- GLM parameters: ch = 0.2, cr = 0.2

**Riemann Solver:** HLLD (5-wave for MHD)
- Handles fast MS waves, Alfvén waves, slow MS waves, contact discontinuity

**Reconstruction:** MUSCL with minmod limiter
- Higher-order spatial accuracy
- TVD (Total Variation Diminishing) property

**Time Integration:** RK2 (2nd order Runge-Kutta)
- Explicit time stepping with CFL control

### Parameters

| Parameter | Brio-Wu | Orszag-Tang |
|-----------|---------|-------------|
| CFL | 0.4 | 0.3 |
| Final time | 0.2 | 5.0 |
| Max steps | 5000 | 10000 |
| Output interval | 50 | 100 |

## Output Files

Both benchmarks write output to `./output/` directory:
- VTK files: `solution_*.vtk` (can be visualized with ParaView)
- HDF5 checkpoint files: `checkpoint_*.h5` (intermediate solutions)

## Performance Expectations

### Runtime

- **Brio-Wu (128³ grid):** ~5-10 minutes (Release mode)
- **Orszag-Tang (96×96×2 grid):** ~5-10 minutes (Release mode)

### System Requirements

- Memory: ~500 MB for Brio-Wu, ~800 MB for Orszag-Tang
- CPU: Any modern processor (serial execution)

## Visualization

### Using ParaView

```bash
paraview output/solution_0000.vtk
```

### Key Quantities to Examine

1. **Density (ρ):** Should show shock structures in Brio-Wu, smooth evolution in Orszag-Tang
2. **Pressure (p):** Complementary to density, shows thermodynamic properties
3. **Velocity magnitude (|u|):** Shows flow patterns
4. **Magnetic field magnitude (|B|):** Shows field evolution
5. **Magnetic field components:** Critical for verifying Alfvén wave handling
6. **Divergence cleaning (ψ):** Should remain small (~10^-6 or smaller)

## Troubleshooting

### Benchmark doesn't run

Check that:
- Build directory exists: `build/`
- Executables exist: `build/brio_wu_shock_tube`, `build/orszag_tang_vortex`
- Output directory can be created: Check write permissions
- Physics type is correct: "mhd_advanced"
- Flux calculator is correct: "hlld"

### Divergence in solution

Possible causes:
- CFL number too large (reduce from 0.4 to 0.3)
- Too coarse grid resolution
- Reconstruction limiter too aggressive
- HLLD solver issue (check denominator stability)

### Negative density or pressure

Indicators:
- Numerical issue in solver
- Reconstruction creating unphysical values
- Need stronger limiter or positivity preservation

## References

1. **Brio-Wu Test:**
   - Brio, M., & Wu, C. C. (1988). "An upwind differencing scheme for the equations of ideal magnetohydrodynamics." Journal of Computational Physics, 75(2), 400-422.
   - Provides exact solution for validation

2. **Orszag-Tang Vortex:**
   - Orszag, S. A., & Tang, C. M. (1979). "Small-scale structure of two-dimensional magnetohydrodynamic turbulence." Journal of Fluid Mechanics, 90(1), 129-143.

3. **HLLD Solver:**
   - Miyoshi, T., & Kusano, K. (2005). "A multi-state HLL approximate Riemann solver for ideal magnetohydrodynamics." Journal of Computational Physics, 208(2), 315-344.
   - Standard reference for HLLD implementation

4. **GLM Divergence Cleaning:**
   - Dedner, A., et al. (2002). "Hyperbolic divergence cleaning for the MHD equations." Journal of Computational Physics, 175(2), 645-673.

## Integration with CI/CD

These benchmarks can be integrated into continuous integration pipelines:

```bash
# Run both benchmarks
cmake --build build --target brio_wu_shock_tube orszag_tang_vortex

# Run with timeout
timeout 600s ./build/brio_wu_shock_tube
timeout 600s ./build/orszag_tang_vortex

# Check for successful completion (exit code 0)
```

## Future Enhancements

- [ ] Add exact solution comparison script for Brio-Wu
- [ ] Add turbulence statistics analysis for Orszag-Tang
- [ ] Create automated validation against published results
- [ ] Add performance benchmarking (execution time, memory usage)
- [ ] Implement grid convergence studies
- [ ] Add 3D MHD tests (e.g., 3D Orszag-Tang)
