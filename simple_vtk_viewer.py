#!/usr/bin/env python
"""
Minimal VTK viewer - loads and displays blast wave data
"""

import pyvista as pv
import glob
import sys
import logging

# Suppress VTK warnings
logging.getLogger().setLevel(logging.ERROR)

def main():
    # Find files
    files = sorted(glob.glob('output_*.vti'))
    if not files:
        print("No output_*.vti files found!")
        sys.exit(1)

    print(f"Found {len(files)} files")

    # Load first file to test
    print(f"\nLoading {files[0]}...")
    try:
        mesh0 = pv.read(files[0])
        print(f"Success! Points: {mesh0.n_points}, Cells: {mesh0.n_cells}")
        print(f"Fields: {mesh0.array_names}")
    except Exception as e:
        print(f"Error loading file: {e}")
        sys.exit(1)

    # Simple visualization of first file
    field = 'density'
    print(f"\nVisualizing {field} from {files[0]}")

    p = pv.Plotter()
    p.add_mesh(mesh0, scalars=field, cmap='jet', show_scalar_bar=True)
    p.add_text(f"{files[0]}\n{field}", font_size=10)
    p.show()

if __name__ == '__main__':
    main()
