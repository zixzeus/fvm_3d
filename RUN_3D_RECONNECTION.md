# 运行3D磁重联模拟 - 完整指南

**日期**: 2025-11-12
**难度**: ⭐⭐ 中等
**运行时间**: 5-10分钟 (单进程)

---

## 快速开始 (3步)

### 1️⃣ 编译

```bash
cd /Users/ericzeus/PycharmProjects/algorithm_2-d_platform/fvm3d/build
cmake ..
make
```

**预期输出**:
```
[ 10%] Building CXX object CMakeFiles/fvm3d_lib.dir/src/physics/resistive_mhd3d_advanced.cpp.o
...
[100%] Built target harris_sheet_3d
```

### 2️⃣ 运行基础测试

```bash
./examples/harris_sheet_3d 64 32 32
```

**预期运行时间**: ~2-3分钟

### 3️⃣ 查看输出

程序会打印时间步进进度和诊断信息。

---

## 详细步骤

### 准备工作

#### 检查依赖

确保已安装：
```bash
# Eigen3
find /usr -name "Eigen" 2>/dev/null | head -1

# MPI
which mpicc
which mpicxx

# HDF5 (可选，用于I/O)
pkg-config --modversion hdf5
```

**如果缺少依赖** (macOS):
```bash
# 使用Homebrew
brew install eigen openmpi

# 或使用MacPorts
sudo port install eigen3 openmpi
```

#### 清除旧的编译

```bash
cd /Users/ericzeus/PycharmProjects/algorithm_2-d_platform/fvm3d/build
rm -rf *
```

### 编译过程

#### Step 1: 配置CMake

```bash
cd /Users/ericzeus/PycharmProjects/algorithm_2-d_platform/fvm3d/build

# 基础配置 (使用系统默认编译器)
cmake ..

# 或指定MPI编译器 (推荐)
cmake -DCMAKE_CXX_COMPILER=mpicxx ..

# Debug模式 (含调试符号)
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Release模式 (优化)
cmake -DCMAKE_BUILD_TYPE=Release ..
```

**输出应包含**:
```
-- Build type: Release
-- C++ Standard: 17
-- C++ Compiler: /usr/bin/c++ (或 /opt/.../mpicxx)
-- Found Eigen3
-- Found HDF5
-- Found MPI
```

#### Step 2: 编译

```bash
make -j4  # 使用4个并行进程
```

或者只编译harris_sheet_3d:
```bash
make harris_sheet_3d
```

**成功标志**:
```
[ 95%] Linking CXX executable examples/harris_sheet_3d
[100%] Built target harris_sheet_3d
```

### 运行程序

#### 基础运行 (推荐开始)

```bash
./examples/harris_sheet_3d 64 32 32
```

**输出示例**:
```
════════════════════════════════════════════════════
  3D Harris Sheet Magnetic Reconnection Simulation
  (Serial Version - Single Process)
════════════════════════════════════════════════════

Grid Configuration:
  Resolution:  64 × 32 × 32
  Domain:      10 × 5 × 4
  Cell size:   0.15625 × 0.15625 × 0.125
  Total cells: 65536

Physics Configuration:
  CFL number:           0.35
  Resistivity model:    Position-dependent
    η₀ (background):    0.001 (Rm₀ = 1000)
    η₁ (enhanced):      0.01667 (Rm₁ = 60)
  GLM parameters:       ch=0.2, cr=0.2

Harris Sheet Equilibrium:
  Sheet thickness:      1
  Plasma beta:          0.2
  Perturbation (m=1):   3%

Initializing Harris sheet equilibrium...
✓ Initialization complete

Starting time integration...
Step | Time    | dt      | KE      | BE      | IE      | max|B_y| | div(B)max
─────┼─────────┼─────────┼─────────┼─────────┼─────────┼──────────┼──────────
   0 | 0.0e+00 | 1.2e-02 | 1.2e-02 | 1.5e+01 | 8.3e-01 | 3.0e-02  | 1.2e-10
  10 | 1.2e-01 | 1.2e-02 | 1.3e-02 | 1.5e+01 | 8.4e-01 | 8.5e-02  | 2.3e-10
  20 | 2.4e-01 | 1.2e-02 | 1.8e-02 | 1.5e+01 | 8.5e-01 | 1.8e-01  | 1.5e-10
...
```

#### 高分辨率运行

```bash
./examples/harris_sheet_3d 128 64 64
```

**运行时间**: ~10-15分钟
**内存使用**: ~500 MB

#### 快速测试

```bash
./examples/harris_sheet_3d 32 16 16
```

**运行时间**: ~10秒
**用途**: 测试编译，验证算法

---

## 理解输出

### 列含义

```
Step    - 时间步编号
Time    - 当前模拟时间
dt      - 本步时间步长 (CFL + 扩散限制)
KE      - 动能 (kinetic energy)
BE      - 磁能 (magnetic energy)
IE      - 内能 (internal energy)
max|B_y| - 最大B_y分量 (磁重联指示器)
div(B)max - ∇·B最大值 (约束满足度)
```

### 物理解释

**正常行为** ✓:
- `max|B_y|` 从小值(0.03)逐步增加 → 磁重联在进行
- `div(B)max` 保持在 ~1e-10 以下 → GLM约束工作良好
- `BE` 逐步下降，转换为 `KE` 和 `IE` → 能量转换

**问题迹象** ⚠️:
- `max|B_y|` 不增加 → 重联未触发，检查初始条件
- `div(B)max` > 0.01 → GLM参数需调整
- 数值爆炸 (NaN) → 时间步太大，减小CFL

---

## 参数调整

### 改变分辨率

```bash
# 低分辨率 (快速测试)
./examples/harris_sheet_3d 32 16 16

# 标准分辨率
./examples/harris_sheet_3d 64 32 32

# 高分辨率 (准确性)
./examples/harris_sheet_3d 128 64 64

# 超高分辨率 (生产)
./examples/harris_sheet_3d 256 128 96
```

### 编辑源代码调整参数

编辑 `examples/harris_sheet_3d.cpp`:

#### 改变CFL数

```cpp
SimConfig sim;
sim.CFL = 0.35;           // 改这里
sim.nsteps = 100;         // 时间步数
sim.diag_freq = 10;       // 诊断频率
```

#### 改变模拟时间

```cpp
// 增加时间步数
sim.nsteps = 500;  // 从100改到500

// 重新编译
make harris_sheet_3d
./examples/harris_sheet_3d 64 32 32
```

#### 改变Harris配置

```cpp
AdvancedResistiveMHD3D::HarrisSheetConfig harris;
harris.L_sheet = 1.0;                    // 电流片厚度
harris.beta = 0.2;                       // 等离子体β
harris.perturbation_amplitude = 0.03;    // 扰动幅度 (改这里试试 0.01 或 0.05)
```

#### 改变电阻率

```cpp
AdvancedResistiveMHD3D::ResistivityModel resistivity;
resistivity.eta0 = 1e-3;              // 背景 Rm₀ = 1000
resistivity.eta1 = 0.01667;           // 增强 Rm₁ = 60 (改这里)
resistivity.localization_scale = 1.0; // 局域化宽度
```

**实验建议**:
- 增加 `eta1` → 更快重联，更强加热
- 减小 `eta1` → 更弱重联，更接近理想MHD
- 增加 `perturbation_amplitude` → 更容易触发
- 减小 `localization_scale` → 更尖锐的梯度

---

## 故障排除

### 编译错误

#### 错误: "找不到 Eigen3"

```bash
# 检查Eigen位置
find /usr -name "Eigen3Config.cmake" 2>/dev/null

# 手动指定
cmake -DEIGEN3_INCLUDE_DIR=/path/to/eigen ..
```

#### 错误: "找不到 MPI"

```bash
# 使用MPI编译器
export CXX=mpicxx
cmake ..
```

### 运行时错误

#### "段错误" (Segmentation Fault)

```
Segmentation fault: 11
```

**可能原因**:
1. 网格分配失败 → 减小分辨率
2. 栈溢出 → 改为堆分配

**解决**:
```bash
# 运行较小的问题
./examples/harris_sheet_3d 32 16 16

# 增加栈大小
ulimit -s unlimited
./examples/harris_sheet_3d 64 32 32
```

#### "浮点异常" (Floating Point Exception)

```
Floating point exception: 8
```

**可能原因**: 除以零或无效的数学操作

**解决**: 检查初始条件或增加floor值

#### "NaN" 或数值爆炸

```
nan or inf detected
```

**可能原因**: 时间步太大

**解决**:
1. 减小 CFL 数：
```cpp
sim.CFL = 0.25;  // 从 0.35 改到 0.25
```

2. 增加扩散稳定因子：
```cpp
sim.dt_resistive_factor = 0.15;  // 从 0.25 改到 0.15
```

---

## 性能优化

### 编译优化

```bash
# Release模式 (更快)
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8

# Debug模式 (更慢但有调试信息)
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j8
```

### 并行编译

```bash
make -j4   # 4个进程
make -j8   # 8个进程 (更快)
make -j16  # 16个进程
```

### 运行时优化

对于大规模运行：
```bash
# 每10步而不是5步输出诊断
# 编辑代码: sim.diag_freq = 10;

# 减少分配: 在每步之间重用内存
# 编辑代码: 不创建临时State3D
```

---

## 预期结果对比

### 小规模 (32×16×16, 10 steps)

```
运行时间:   ~5 秒
内存:       ~50 MB
max|B_y|:   从0.03增加到~0.1-0.15
div(B)max:  < 1e-9
```

### 标准 (64×32×32, 100 steps)

```
运行时间:   ~2-3 分钟
内存:       ~200 MB
max|B_y|:   从0.03增加到~0.3-0.5
div(B)max:  < 1e-9
重联指标:   明显的磁岛形成
```

### 高分辨率 (128×64×64, 200 steps)

```
运行时间:   ~30-40 分钟
内存:       ~1 GB
max|B_y|:   从0.03增加到~0.6-0.8
div(B)max:  < 1e-9
物理:       详细的重联动力学
```

---

## 与OpenMHD对比

如果有OpenMHD参考结果：

### OpenMHD结果通常显示:
```
t = 0:     max(By) = 0.03  (初始扰动)
t = 50:    max(By) ≈ 0.15  (线性增长)
t = 100:   max(By) ≈ 0.30  (快速增长)
t = 150:   max(By) ≈ 0.50  (非线性饱和)
t = 200+:  max(By) > 0.60  (重联进入湍流阶段)
```

### 我们的FVM3D实现应该显示类似趋势

---

## 高级: MPI并行运行

(需要Phase 6完成)

```bash
# 2进程
mpirun -np 2 ./examples/harris_sheet_3d_mpi 64 32 32

# 8进程
mpirun -np 8 ./examples/harris_sheet_3d_mpi 64 32 32

# 在集群上
sbatch -N 4 -n 32 run_reconnection.sh
```

---

## 提示与技巧

### 快速验证算法

```bash
# 快速 (~5秒)
./examples/harris_sheet_3d 32 16 16

# 输出应显示:
# 1. max|B_y| 逐步增加
# 2. div(B)max 保持小
# 3. 没有数值爆炸
```

### 诊断磁重联进度

观察 `max|B_y|` 列：

```
Step |  max|B_y|  | 物理阶段
-----|------------|----------
   0 |   0.030    | 初始扰动
  10 |   0.032    | 缓慢增长
  20 |   0.040    | 线性阶段
  30 |   0.065    | 加速增长
  50 |   0.150    | 快速增长
 100 |   0.300+   | 非线性阶段
 200 |   0.600+   | 重联湍流
```

### 监控GLM约束

`div(B)max` 应该：
- 保持 < 0.001 (好)
- 保持 < 0.01 (可接受)
- 增长 > 0.1 (问题!)

如果增长过快，调整:
```cpp
glm.ch = 0.4;  // 增加波速
glm.cr = 0.3;  // 增加衰减
```

---

## 下一步

1. ✅ 运行此基础版本
2. ⏳ (Phase 5) 实现HDF5输出保存结果
3. ⏳ (Phase 6) 实现MPI并行版本
4. 📊 可视化结果 (ParaView/VisIt)
5. 📝 论文发表!

---

## 常见问题

**Q: 为什么max|B_y|很小?**
A: 扰动可能太小，或者Rm太低。增加perturbation_amplitude到0.05或增加eta1。

**Q: 为什么很慢?**
A: 分辨率太高或CFL太小。从32×16×16开始，逐步增加。

**Q: 怎样保存结果?**
A: 等待Phase 5 (HDF5 I/O)，或修改代码在每步输出。

**Q: 能并行运行吗?**
A: 需要Phase 6 (MPI集成)，敬请期待!

---

**祝您模拟顺利!** 🚀

