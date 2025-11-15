# VTK可视化脚本使用说明

本目录包含多个Python脚本用于可视化3D爆炸波模拟的VTK输出文件。

## 前置要求

确保已安装PyVista：
```bash
conda activate base  # 或你的conda环境
pip install pyvista
```

## 运行模拟生成VTK文件

```bash
cd build
./vtk_visualization_example
```

这将在build目录下生成10个VTK文件：`output_0000.vti` 到 `output_0009.vti`

## 可视化脚本

### 1. `simple_vtk_viewer.py` - 最简单的查看器

显示第一个时间步的3D可视化：

```bash
cd build  # 确保在有output_*.vti文件的目录
python ../simple_vtk_viewer.py
```

### 2. `visualize_vtk.py` - 交互式时间序列查看器（推荐）

带滑动条的交互式查看器，可以浏览所有时间步：

```bash
python ../visualize_vtk.py
```

**功能：**
- 使用滑动条切换不同时间步
- 显示密度场演化
- 实时更新colorbar
- 等值面视图

### 3. `visualize_blast_wave.py` - 全功能可视化工具

支持多种可视化模式：

```bash
# 交互式查看器（带时间滑块）
python ../visualize_blast_wave.py --mode interactive --field density

# 网格布局显示所有时间步
python ../visualize_blast_wave.py --mode grid --field density

# 2D切片演化
python ../visualize_blast_wave.py --mode slice --field pressure --normal z

# 创建动画GIF
python ../visualize_blast_wave.py --mode animation --field density --output blast.gif --fps 2
```

**可用字段：**
- `density` - 密度
- `pressure` - 压强
- `temperature` - 温度
- `energy` - 能量
- `velocity` - 速度矢量
- `momentum_x/y/z` - 动量分量

**可用选项：**
```
--mode      : 可视化模式 (interactive, grid, slice, animation)
--field     : 要显示的物理场
--colormap  : 颜色图 (jet, viridis, plasma, etc.)
--fps       : 动画帧率
--normal    : 切片方向 (x, y, z)
```

## 使用ParaView查看（替代方案）

如果Python脚本有问题，也可以直接用ParaView打开：

1. 打开ParaView
2. File → Open → 选择 `output_*.vti` 文件
3. 点击 "Apply"
4. 在 "Coloring" 下拉菜单选择要显示的场（density, pressure等）
5. 使用动画控制按钮浏览时间序列

## 输出文件说明

| 文件 | 时间 (t) | 描述 |
|------|----------|------|
| output_0000.vti | 0.00 | 初始状态（高压区） |
| output_0001.vti | 0.10 | 爆炸波开始传播 |
| output_0002.vti | 0.20 | 激波扩散 |
| output_0003.vti | 0.30 | 激波继续传播 |
| ... | ... | ... |
| output_0009.vti | 0.80 | 最终状态 |

## 故障排除

### 问题：找不到output_*.vti文件
**解决：** 确保在build目录下运行脚本，或使用绝对路径

### 问题：PyVista显示错误
**解决：** 尝试更新PyVista版本
```bash
pip install --upgrade pyvista
```

### 问题：无法显示GUI窗口
**解决：**
- 确保有图形界面环境
- 或者使用 `--mode animation` 生成GIF文件

### 问题：内存不足
**解决：** 只可视化部分时间步，或使用切片模式而非3D体渲染
