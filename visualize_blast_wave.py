#!/usr/bin/env python3
"""
3D Blast Wave Visualization Script

Visualizes VTK output files from the FVM solver using PyVista.

Usage:
    python visualize_blast_wave.py --mode animation  # Create animation
    python visualize_blast_wave.py --mode grid       # Show grid of snapshots
    python visualize_blast_wave.py --mode slice      # Show 2D slices
    python visualize_blast_wave.py --mode interactive # Interactive viewer
"""

import pyvista as pv
import numpy as np
import glob
import argparse
from pathlib import Path


def load_vtk_files(pattern="output_*.vti"):
    """Load all VTK files matching the pattern."""
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files found matching pattern: {pattern}")

    print(f"Found {len(files)} VTK files:")
    for f in files:
        print(f"  - {f}")

    return files


def get_field_info(mesh):
    """Print available fields in the mesh."""
    print("\nAvailable scalar fields:")
    for name in mesh.array_names:
        arr = mesh[name]
        if len(arr.shape) == 1:
            print(f"  - {name}: min={arr.min():.3e}, max={arr.max():.3e}")
        else:
            print(f"  - {name}: vector field, shape={arr.shape}")


def create_animation(files, field='density', output_file='blast_wave.gif',
                     fps=2, colormap='jet'):
    """Create an animated GIF showing the evolution."""
    print(f"\nCreating animation with field: {field}")

    # Create plotter
    plotter = pv.Plotter(off_screen=True)
    plotter.open_gif(output_file, fps=fps)

    for i, file in enumerate(files):
        mesh = pv.read(file)

        # Add time annotation
        time = mesh.field_data.get('TIME', [i * 0.1])[0] if 'TIME' in mesh.field_data else i * 0.1

        plotter.clear()
        plotter.add_mesh(mesh, scalars=field, cmap=colormap,
                        show_edges=False, opacity=1.0,
                        clim=[mesh[field].min(), mesh[field].max()])
        plotter.add_text(f"t = {time:.3f}", position='upper_left', font_size=12)
        plotter.camera_position = 'iso'
        plotter.write_frame()

    plotter.close()
    print(f"Animation saved to: {output_file}")


def show_grid_comparison(files, field='density', rows=2, colormap='jet'):
    """Show multiple snapshots in a grid layout."""
    print(f"\nShowing grid comparison with field: {field}")

    n_files = len(files)
    cols = int(np.ceil(n_files / rows))

    plotter = pv.Plotter(shape=(rows, cols), window_size=(400*cols, 400*rows))

    for i, file in enumerate(files):
        row = i // cols
        col = i % cols
        plotter.subplot(row, col)

        mesh = pv.read(file)
        time = mesh.field_data.get('TIME', [i * 0.1])[0] if 'TIME' in mesh.field_data else i * 0.1

        # Get global min/max for consistent colorbar
        if i == 0:
            global_min = mesh[field].min()
            global_max = mesh[field].max()
        else:
            global_min = min(global_min, mesh[field].min())
            global_max = max(global_max, mesh[field].max())

        plotter.add_mesh(mesh, scalars=field, cmap=colormap,
                        show_edges=False, clim=[global_min, global_max])
        plotter.add_text(f"t = {time:.3f}", font_size=10)
        plotter.camera_position = 'iso'

    plotter.show()


def show_slice_evolution(files, field='density', normal='z', colormap='jet'):
    """Show 2D slices through the center for all timesteps."""
    print(f"\nShowing slice evolution with field: {field}")

    n_files = len(files)
    rows = int(np.ceil(n_files / 3))
    cols = min(3, n_files)

    plotter = pv.Plotter(shape=(rows, cols), window_size=(400*cols, 400*rows))

    for i, file in enumerate(files):
        row = i // cols
        col = i % cols
        plotter.subplot(row, col)

        mesh = pv.read(file)
        time = mesh.field_data.get('TIME', [i * 0.1])[0] if 'TIME' in mesh.field_data else i * 0.1

        # Create slice through center
        if normal == 'x':
            slice_mesh = mesh.slice(normal='x', origin=mesh.center)
        elif normal == 'y':
            slice_mesh = mesh.slice(normal='y', origin=mesh.center)
        else:  # z
            slice_mesh = mesh.slice(normal='z', origin=mesh.center)

        plotter.add_mesh(slice_mesh, scalars=field, cmap=colormap, show_edges=False)
        plotter.add_text(f"t = {time:.3f}", font_size=10)
        plotter.view_xy()

    plotter.show()


def interactive_viewer(files, field='density', colormap='jet'):
    """Interactive viewer with slider to control time."""
    print(f"\nLaunching interactive viewer with field: {field}")

    # Load all meshes
    meshes = [pv.read(f) for f in files]
    times = [m.field_data.get('TIME', [i * 0.1])[0] if 'TIME' in m.field_data
             else i * 0.1 for i, m in enumerate(meshes)]

    # Get global min/max for consistent colorbar
    global_min = min(m[field].min() for m in meshes)
    global_max = max(m[field].max() for m in meshes)

    plotter = pv.Plotter()

    # Add initial mesh
    current_actor = plotter.add_mesh(meshes[0], scalars=field, cmap=colormap,
                                     clim=[global_min, global_max],
                                     show_edges=False, name='mesh')
    text_actor = plotter.add_text(f"t = {times[0]:.3f}", position='upper_left',
                                  font_size=12, name='time_text')

    def update_mesh(value):
        """Callback to update mesh based on slider."""
        idx = int(value)
        plotter.remove_actor('mesh')
        plotter.remove_actor('time_text')

        plotter.add_mesh(meshes[idx], scalars=field, cmap=colormap,
                        clim=[global_min, global_max],
                        show_edges=False, name='mesh')
        plotter.add_text(f"t = {times[idx]:.3f}", position='upper_left',
                        font_size=12, name='time_text')

    plotter.add_slider_widget(update_mesh, [0, len(meshes)-1],
                             value=0, title="Time Step",
                             pointa=(0.1, 0.1), pointb=(0.4, 0.1))

    plotter.camera_position = 'iso'
    plotter.show()


def show_multi_field_comparison(file_idx=0, fields=['density', 'pressure', 'velocity_magnitude']):
    """Show multiple fields for a single timestep."""
    files = load_vtk_files()
    mesh = pv.read(files[file_idx])

    # Calculate velocity magnitude if needed
    if 'velocity_magnitude' in fields and 'velocity_magnitude' not in mesh.array_names:
        if 'velocity' in mesh.array_names:
            vel = mesh['velocity']
            mesh['velocity_magnitude'] = np.linalg.norm(vel, axis=1)

    n_fields = len(fields)
    plotter = pv.Plotter(shape=(1, n_fields), window_size=(400*n_fields, 400))

    for i, field in enumerate(fields):
        plotter.subplot(0, i)
        if field in mesh.array_names:
            plotter.add_mesh(mesh, scalars=field, cmap='jet', show_edges=False)
            plotter.add_text(f"{field}", font_size=10)
            plotter.camera_position = 'iso'
        else:
            plotter.add_text(f"{field}\n(not available)", font_size=10)

    plotter.show()


def main():
    parser = argparse.ArgumentParser(description='Visualize 3D blast wave VTK files')
    parser.add_argument('--mode', type=str, default='interactive',
                       choices=['animation', 'grid', 'slice', 'interactive', 'multi'],
                       help='Visualization mode')
    parser.add_argument('--field', type=str, default='density',
                       help='Field to visualize (default: density)')
    parser.add_argument('--colormap', type=str, default='jet',
                       help='Colormap (default: jet)')
    parser.add_argument('--pattern', type=str, default='output_*.vti',
                       help='File pattern (default: output_*.vti)')
    parser.add_argument('--output', type=str, default='blast_wave.gif',
                       help='Output file for animation (default: blast_wave.gif)')
    parser.add_argument('--fps', type=int, default=2,
                       help='Frames per second for animation (default: 2)')
    parser.add_argument('--normal', type=str, default='z', choices=['x', 'y', 'z'],
                       help='Normal direction for slice mode (default: z)')

    args = parser.parse_args()

    # Load files
    files = load_vtk_files(args.pattern)

    # Show field info for first file
    mesh = pv.read(files[0])
    get_field_info(mesh)

    # Execute visualization based on mode
    if args.mode == 'animation':
        create_animation(files, field=args.field, output_file=args.output,
                        fps=args.fps, colormap=args.colormap)
    elif args.mode == 'grid':
        show_grid_comparison(files, field=args.field, colormap=args.colormap)
    elif args.mode == 'slice':
        show_slice_evolution(files, field=args.field, normal=args.normal,
                            colormap=args.colormap)
    elif args.mode == 'interactive':
        interactive_viewer(files, field=args.field, colormap=args.colormap)
    elif args.mode == 'multi':
        show_multi_field_comparison(file_idx=len(files)//2)


if __name__ == '__main__':
    main()
