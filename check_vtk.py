#!/usr/bin/env python3
"""
Quick VTK file checker and visualizer.
Tests if PyVista can read the VTK files and shows basic info.
"""

import sys
import glob

try:
    import pyvista as pv
    import numpy as np
    print("✓ PyVista is installed")
except ImportError:
    print("✗ PyVista is not installed")
    print("\nTo install PyVista, run:")
    print("  pip install pyvista")
    print("  or")
    print("  conda install -c conda-forge pyvista")
    sys.exit(1)

# Find VTK files
files = sorted(glob.glob('output_*.vti'))
if not files:
    print("No output_*.vti files found in current directory")
    sys.exit(1)

print(f"\n✓ Found {len(files)} VTK files")

# Try to read first file
print(f"\nTesting file: {files[0]}")
try:
    mesh = pv.read(files[0])
    print("✓ File loaded successfully!")
    print(f"\nMesh information:")
    print(f"  Points: {mesh.n_points}")
    print(f"  Cells: {mesh.n_cells}")
    print(f"  Bounds: {mesh.bounds}")

    print(f"\nAvailable fields:")
    for name in mesh.array_names:
        arr = mesh[name]
        if len(arr.shape) == 1:
            print(f"  - {name}: min={arr.min():.3e}, max={arr.max():.3e}")
        else:
            print(f"  - {name}: vector field {arr.shape}")

    print("\n✓ VTK files are valid and can be visualized!")
    print("\nTo visualize, run:")
    print("  python visualize_blast_wave.py --mode interactive")
    print("  python visualize_blast_wave.py --mode grid")
    print("  python visualize_blast_wave.py --mode animation")

except Exception as e:
    print(f"✗ Error reading file: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
