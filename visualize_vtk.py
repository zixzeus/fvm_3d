#!/usr/bin/env python
"""
Simple VTK visualization script for blast wave simulation.
"""

import pyvista as pv
import numpy as np
import glob
import sys

# Set PyVista to work without display (optional)
# pv.set_plot_theme("document")

def main():
    # Find all VTK files
    files = sorted(glob.glob('output_*.vti'))
    if not files:
        print("No output_*.vti files found!")
        print("Current directory:", __import__('os').getcwd())
        sys.exit(1)

    print(f"Found {len(files)} VTK files\n")

    # Load all meshes
    meshes = []
    times = []

    for i, fname in enumerate(files):
        try:
            mesh = pv.read(fname)
            meshes.append(mesh)
            # Try to get time from filename if not in data
            time = i * 0.1  # default
            times.append(time)
            print(f"Loaded {fname}: {mesh.n_points} points, {mesh.n_cells} cells")
        except Exception as e:
            print(f"Warning: Could not load {fname}: {e}")

    if not meshes:
        print("No meshes loaded successfully!")
        sys.exit(1)

    print(f"\nAvailable fields: {meshes[0].array_names}")

    # Choose field to visualize
    field = 'density'
    if field not in meshes[0].array_names:
        field = meshes[0].array_names[0]

    print(f"\nVisualizing field: {field}")

    # Get global min/max for consistent colorbar
    vmin = min(m[field].min() for m in meshes)
    vmax = max(m[field].max() for m in meshes)
    print(f"Field range: [{vmin:.3e}, {vmax:.3e}]")

    # Create interactive plotter with slider
    plotter = pv.Plotter()
    plotter.add_text(f"Use slider to change timestep\nField: {field}",
                     position='upper_left', font_size=10)

    # Add initial mesh
    actor = plotter.add_mesh(meshes[0], scalars=field,
                            clim=[vmin, vmax],
                            cmap='jet',
                            show_scalar_bar=True,
                            scalar_bar_args={'title': field})

    time_text = plotter.add_text(f"t = {times[0]:.3f}",
                                 position='lower_right',
                                 font_size=12)

    def update_plot(value):
        """Update callback for slider."""
        idx = int(value)
        plotter.clear()

        plotter.add_mesh(meshes[idx], scalars=field,
                        clim=[vmin, vmax],
                        cmap='jet',
                        show_scalar_bar=True,
                        scalar_bar_args={'title': field})

        plotter.add_text(f"Use slider to change timestep\nField: {field}",
                        position='upper_left', font_size=10)
        plotter.add_text(f"t = {times[idx]:.3f} (step {idx}/{len(meshes)-1})",
                        position='lower_right', font_size=12)

    # Add slider widget
    plotter.add_slider_widget(
        update_plot,
        [0, len(meshes)-1],
        value=0,
        title="Time Step",
        pointa=(0.1, 0.1),
        pointb=(0.4, 0.1),
        style='modern'
    )

    # Set camera
    plotter.camera_position = 'iso'
    plotter.show()


if __name__ == '__main__':
    main()
