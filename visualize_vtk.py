#!/usr/bin/env python3
"""
VTK文件可视化脚本
使用 PyVista 快速查看 .vti 文件

安装：pip install pyvista

用法：
    python visualize_vtk.py output_0000.vti [mode]

模式：
    slice     - 三平面切片（默认）
    volume    - 3D体渲染
    isosurface - 等值面
    contour   - 多个等值面
"""

import pyvista as pv
import sys
import glob
import numpy as np

def visualize_single(filename, mode='slice'):
    """可视化单个VTK文件"""
    print(f"加载: {filename}")
    print(f"模式: {mode}")

    # 读取VTK文件
    mesh = pv.read(filename)

    # 显示数据信息
    print(f"网格尺寸: {mesh.dimensions}")
    print(f"数据字段: {mesh.array_names}")

    # 选择要可视化的字段
    if 'density' in mesh.array_names:
        field_name = 'density'
    else:
        field_name = mesh.array_names[0]

    print(f"显示字段: {field_name}")

    # 获取数据范围
    data_range = mesh.get_data_range(field_name)
    print(f"数据范围: {data_range[0]:.3f} - {data_range[1]:.3f}")

    # 创建绘图器
    plotter = pv.Plotter()

    if mode == 'volume':
        # 3D体渲染 - 真正的3D效果
        print("使用体渲染...")
        plotter.add_volume(
            mesh,
            scalars=field_name,
            cmap='hot',
            opacity='sigmoid',  # 自动透明度映射
            shade=True
        )

    elif mode == 'isosurface':
        # 单个等值面
        print("生成等值面...")
        iso_value = np.mean(data_range)  # 使用中间值
        contour = mesh.contour([iso_value], scalars=field_name)
        plotter.add_mesh(
            contour,
            scalars=field_name,
            cmap='viridis',
            show_edges=False,
            smooth_shading=True
        )

    elif mode == 'contour':
        # 多个等值面
        print("生成多个等值面...")
        n_contours = 5
        contours = mesh.contour(n_contours, scalars=field_name)
        plotter.add_mesh(
            contours,
            scalars=field_name,
            cmap='plasma',
            opacity=0.7,
            show_edges=False
        )

    else:  # slice (默认)
        # 三平面切片
        print("使用三平面切片...")
        slices = mesh.slice_orthogonal()
        plotter.add_mesh(
            slices,
            scalars=field_name,
            cmap='viridis',
            show_edges=False
        )

    plotter.add_scalar_bar(field_name, title=field_name.capitalize())
    plotter.show_axes()
    plotter.show_grid()
    plotter.add_text(f"File: {filename}\nMode: {mode}", position='upper_left', font_size=10)
    plotter.show()

def visualize_animation(filenames):
    """可视化时间序列动画"""
    filenames = sorted(filenames)
    print(f"找到 {len(filenames)} 个时间步")

    # 创建绘图器
    plotter = pv.Plotter()
    plotter.open_gif("animation.gif")  # 可选：保存为GIF

    for filename in filenames:
        print(f"处理: {filename}")
        mesh = pv.read(filename)

        plotter.clear()

        # 显示密度场的切片
        if 'density' in mesh.array_names:
            plotter.add_mesh(mesh.slice_orthogonal(),
                             scalars='density',
                             cmap='viridis',
                             clim=[0.8, 2.0])  # 固定颜色范围
            plotter.add_scalar_bar('Density')

        plotter.show_axes()
        plotter.write_frame()  # 写入GIF帧

    plotter.close()
    print("动画已保存为 animation.gif")

def visualize_3d_volume(filename):
    """体渲染（3D内部结构）"""
    mesh = pv.read(filename)

    plotter = pv.Plotter()
    plotter.add_volume(mesh,
                       scalars='density',
                       cmap='hot',
                       opacity='sigmoid')
    plotter.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python visualize_vtk.py <vtk_file> [mode]")
        print("\n示例:")
        print("  python visualize_vtk.py output_0000.vti           # 三平面切片")
        print("  python visualize_vtk.py output_0000.vti volume    # 3D体渲染")
        print("  python visualize_vtk.py output_0000.vti isosurface # 等值面")
        print("  python visualize_vtk.py output_0000.vti contour   # 多等值面")
        print("  python visualize_vtk.py 'output_*.vti'            # 动画")
        print("\n模式说明:")
        print("  slice      - 三个正交切片平面（默认）")
        print("  volume     - 3D体渲染，半透明显示整个数据体")
        print("  isosurface - 单个等值面（3D表面）")
        print("  contour    - 多个等值面")
        sys.exit(1)

    pattern = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else 'slice'

    # 展开通配符
    files = glob.glob(pattern)

    if not files:
        print(f"错误: 找不到匹配 '{pattern}' 的文件")
        sys.exit(1)

    if len(files) == 1:
        # 单个文件
        visualize_single(files[0], mode)
    else:
        # 多个文件 - 制作动画
        if mode != 'slice':
            print(f"警告: 动画模式下忽略 '{mode}'，使用默认 'slice' 模式")
        visualize_animation(files)
