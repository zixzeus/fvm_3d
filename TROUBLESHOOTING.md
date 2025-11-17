# FVM3D 故障排除指南

## 常见问题和解决方案

### 1. CFL时间步变为零错误

**错误信息**:
```
compute_dt call #9: min_dt=0.125746, max_wave_speed=0.994071
Error: CFL time step must be positive: 0.000000
```

**原因**:
- 数值格式不稳定，导致状态变量出现非物理值（负压力、NaN、Inf）
- 波速计算失败，进而导致时间步计算失败

**解决方案**:

#### 方案A：使用MUSCL重构（推荐）✅
```cpp
config.reconstruction = "muscl";  // 而不是 "constant"
config.cfl = 0.3;                 // MUSCL的典型值
```

**原理**: MUSCL重构提供二阶精度和更好的激波捕捉能力，对于磁重联等不连续问题更稳定。

#### 方案B：降低CFL数
```cpp
config.cfl = 0.05;  // 非常保守的值
```

#### 方案C：添加状态验证
在`compute_dt()`中添加更严格的检查：
```cpp
// 在 src/core/mpi_fvm_solver3d.cpp:439
double p = (5.0/3.0 - 1.0) * internal_energy;
p = std::max(p, 1e-12);  // 压力地板

// 更严格的检查
if (p > 1e-12 && !std::isnan(p) && !std::isinf(p)) {
    // ... 计算波速
}
```

---

### 2. OpenCL警告（macOS）

**警告信息**:
```
[CL_INVALID_OPERATION] : OpenCL Error : Failed to retrieve device information!
```

**原因**: macOS系统弃用了OpenCL支持，但Eigen库或系统框架仍尝试查询OpenCL设备。

**影响**: **无影响**，这是非致命警告，程序正常运行。

**消除方法** (可选):
```bash
# 设置环境变量禁用OpenCL查询
export CL_DISABLE=1
```

或在代码中添加：
```cpp
// 在main()开始处
putenv("CL_DISABLE=1");
```

---

### 3. MPI进程数与网格不匹配

**警告信息**:
```
Warning: Process 0 has 10000 cells, Process 7 has 8000 cells
```

**原因**: 网格尺寸不能被MPI进程数整除，导致负载不均衡。

**解决方案**:
选择可整除的网格尺寸：
```cpp
// 示例：8个MPI进程 (2×2×2)
config.nx = 64;  // 64/2 = 32 ✅
config.ny = 64;  // 64/2 = 32 ✅
config.nz = 64;  // 64/2 = 32 ✅
```

---

### 4. HDF5并行支持缺失

**错误信息**:
```
Error: H5Pset_fapl_mpio not found
```

**原因**: HDF5编译时未启用并行支持。

**解决方案**:
```bash
# 重新编译HDF5
./configure --enable-parallel --enable-shared \
    CC=mpicc CXX=mpicxx
make && sudo make install

# 验证
h5cc -showconfig | grep "Parallel HDF5"
# 应该显示: Parallel HDF5: yes
```

---

### 5. 磁重联模拟不稳定

**症状**:
- 压力变负
- ∇·B快速增长
- 数值爆炸

**诊断检查清单**:

1. **检查初始条件平衡**:
```cpp
// Harris sheet必须满足压力平衡
harris_config.beta = 1.0;  // 典型值0.5-2.0
harris_config.L_sheet = 0.5;  // 不要太小
```

2. **检查电阻率设置**:
```cpp
resistivity.eta0 = 1e-3;   // 背景电阻率
resistivity.eta1 = 0.01667;  // 增强电阻率
// Rm_min = 1/eta1 ≈ 60 是合理的
```

3. **使用合适的数值格式**:
```cpp
config.reconstruction = "muscl";  // 二阶精度
config.riemann_solver = "hlld";   // MHD专用
config.time_integrator = "rk2";   // 或"rk3"
config.cfl = 0.3;                 // 0.2-0.4之间
```

4. **启用GLM散度清理**:
```cpp
glm.ch = 0.2;  // 散度传播速度
glm.cr = 0.2;  // 阻尼系数
```

---

### 6. 编译错误：找不到Eigen头文件

**错误信息**:
```
fatal error: Eigen/Dense: No such file or directory
```

**解决方案**:
```bash
# Ubuntu/Debian
sudo apt-get install libeigen3-dev

# macOS
brew install eigen

# 或手动安装
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
mkdir build && cd build
cmake ..
sudo make install
```

---

### 7. MPI运行时错误

**错误信息**:
```
ORTE_ERROR_LOG: Not found in file orted/pmix/pmix_server.c
```

**解决方案**:
```bash
# 检查MPI环境
which mpirun
mpirun --version

# 确保编译和运行使用相同的MPI
export PATH=/path/to/mpi/bin:$PATH
export LD_LIBRARY_PATH=/path/to/mpi/lib:$LD_LIBRARY_PATH

# 重新编译
cd build
rm -rf *
cmake .. -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc
make -j4
```

---

## 性能优化建议

### 1. 选择合适的重构方案

| 方案 | 精度 | 稳定性 | 适用场景 |
|------|------|--------|---------|
| constant | 一阶 | 高 | 光滑问题、调试 |
| muscl | 二阶 | 中 | 大部分问题（推荐） |
| weno | 高阶 | 低 | 光滑问题、高精度需求 |

### 2. 选择合适的CFL数

| 时间积分器 | 重构方案 | 推荐CFL | 最大CFL |
|-----------|---------|---------|---------|
| Forward Euler | constant | 0.5 | 0.9 |
| RK2 | constant | 0.8 | 1.0 |
| RK2 | muscl | 0.3 | 0.5 |
| RK3 | muscl | 0.5 | 0.8 |

### 3. MPI并行效率优化

```cpp
// 每进程单元格数建议: 10万-100万
// 示例：256个进程，512³全局网格
config.nx = 512;  // 512³ = 134M cells
config.ny = 512;
config.nz = 512;
// 每进程: 134M/256 = 524K cells ✅

// 避免过度分解
config.px = 8;   // 而不是16
config.py = 8;   // 保持局部域尺寸合理
config.pz = 4;
```

---

## 调试技巧

### 1. 启用详细输出

```cpp
config.verbose = 2;  // 0=静默, 1=基本, 2=详细, 3=调试
```

### 2. 添加状态检查

```cpp
// 在求解器循环中
for (int step = 0; step < max_steps; step++) {
    solver.step();

    // 检查NaN
    const auto& state = solver.state();
    bool has_nan = false;
    for (int v = 0; v < num_vars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    if (std::isnan(state(v,i,j,k))) {
                        std::cout << "NaN at step=" << step
                                  << " v=" << v << " (" << i << "," << j << "," << k << ")\n";
                        has_nan = true;
                    }
                }
            }
        }
    }
    if (has_nan) break;
}
```

### 3. 可视化诊断

```bash
# 使用ParaView查看VTK输出
paraview reconnection_0000.pvti

# 检查：
# - 密度/压力是否为正
# - ∇·B是否接近零
# - 数值是否合理范围
```

---

## 联系和支持

如果问题仍未解决：

1. **检查日志**: 查看完整的错误堆栈
2. **最小复现**: 创建最小的测试案例
3. **提供信息**:
   - 系统环境（OS, 编译器, MPI版本）
   - 完整错误信息
   - 配置参数
   - 初始条件
4. **提交Issue**: 包含上述信息到项目仓库

---

*最后更新: 2025-11-14*
