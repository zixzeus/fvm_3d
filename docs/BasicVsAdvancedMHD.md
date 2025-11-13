# 基础 vs 高级 电阻MHD 实现对比

**日期**: 2025-11-12
**框架**: FVM3D with MPI Parallelization

---

## 快速对比

| 特性 | 基础MHD | 高级MHD (磁重联) |
|------|---------|------------------|
| **主要文件** | `resistive_mhd3d.hpp/.cpp` | `resistive_mhd3d_advanced.hpp/.cpp` |
| **守恒变量** | 8个 | 9个 (+ GLM) |
| **电阻率** | 常数 η | 位置相关 η(x,y,z) |
| **GLM约束** | ✗ (通量形式维持) | ✓ (显式维持) |
| **源项** | 简单 | 完整 |
| **初始条件** | 通用 | Harris平衡 |
| **适用场景** | 通用MHD问题 | 磁重联模拟 |
| **代码行数** | ~390行 | ~1000行 |
| **学习曲线** | 低 | 中等 |
| **计算成本** | 标准 | +~30-50% |
| **物理保真度** | 良好 | 优秀 |

---

## 详细对比

### 1. 守恒变量

**基础MHD (8个)**:
```
U = [ρ, ρu, ρv, ρw, E, Bx, By, Bz]
    [0  1   2   3   4  5   6   7]
```

**高级MHD (9个)**:
```
U = [ρ, ρu, ρv, ρw, E, Bx, By, Bz, ψ]
    [0  1   2   3   4  5   6   7   8]
```

**ψ字段含义**:
- GLM (Generalized Lagrangian Multiplier) 散度清理标量
- 跟踪 ∇·B 的数值误差
- 通过双曲-抛物方程演化
- 自动衰减到零，维持约束

---

### 2. 电阻率模型

#### 基础MHD：常数电阻

```cpp
explicit ResistiveMHD3D(double resistivity = 0.0)
    : eta_(resistivity) {}

// 在源项中使用常数 eta_
S(4) = eta_ * J_sq;  // 所有地方相同
S(5) = eta_ * lapl_Bx;
```

**物理意义**:
- 整个域内统一的电阻性质
- 适合局域被加热情景
- 简单参数化

**应用**:
- MHD冲击波
- 等离子体加热
- 碰撞效应的简单模型

#### 高级MHD：位置相关电阻

```cpp
struct ResistivityModel {
    double eta0;              // 背景电阻: 1e-3
    double eta1;              // 增强电阻: 1.667e-2
    double localization_scale; // 宽度: 1.0

    double operator()(double x, double y, double z) const {
        double r_sq = x*x + y*y;
        double r = std::sqrt(r_sq);
        r = std::min(r, 25.0);
        double sech_sq = 1.0 / std::cosh(r / localization_scale);
        sech_sq *= sech_sq;
        return eta0 + (eta1 - eta0) * sech_sq;
    }
};
```

**数学形式**:
```
η(x,y) = η₀ + (η₁ - η₀)·sech²(r/L_r)

其中:
  r = √(x² + y²)
  η₀ = 1/1000 = 0.001         [背景: 弱]
  η₁ = 1/60 ≈ 0.01667        [增强: 强]
  L_r = 1.0                  [转换宽度]
```

**物理含义**:
- **高η区域** (原点): 强烈欧姆耗散 → 快速重联
- **低η区域** (远处): 接近理想MHD → 磁场线冻结
- **梯度**: 能量集中在原点附近释放

**磁Reynolds数**:
```
背景: Rm₀ = 1/η₀ = 1000  → 磁场冻结，切割动力学占主
中心: Rm₁ = 1/η₁ = 60    → 扩散与对流竞争，重联发生
```

---

### 3. 磁场约束（∇·B = 0）

#### 基础MHD：单一策略

**仅依赖通量形式**:
```
∂B/∂t + ∇·F_B = 0

其中 F_B 设计使得:
  ∇·F_x = 0 (X-分量)

这自动意味着:
  ∇·B = ∂Bx/∂x + ∂By/∂y + ∂Bz/∂z ≈ 0 (机器精度)
```

**优点**:
- ✓ 在无源项时精确保持
- ✓ 无额外计算成本
- ✓ 实现简洁

**缺点**:
- ✗ 源项会产生错误
- ✗ 数值舍入累积
- ✗ 需要谨慎实现源项

#### 高级MHD：双层策略

**策略1：通量形式** (主要)
```
∂B/∂t + ∇·F_B = 0
```

**策略2：GLM显式维持** (辅助)
```
∂ψ/∂t = -ch·(∇·B) - (cr/ch)·ψ

双曲传播:    ψ 波以速度 ch 传播
抛物衰减:    ψ 指数衰减到零
```

**源项附加到ψ方程**:
```
S_ψ = -ch·(∇·B) - (cr/ch)·ψ
```

**磁场交互**:
```
∂Bx/∂t 包含项: ch·∂ψ/∂x
```

**优点**:
- ✓ 显式错误检测和清理
- ✓ 超不稳定数值设置仍可运行
- ✓ 添加物理约束
- ✓ OpenMHD验证方法

**缺点**:
- ✗ 增加计算成本 (~5%)
- ✗ 引入额外参数 (ch, cr)
- ✗ 一个额外变量

---

### 4. 源项处理

#### 基础MHD

**电阻耗散源**:
```cpp
Eigen::VectorXd resistive_source(
    const Eigen::VectorXd& U_center,
    double dx, double dy, double dz,
    const Eigen::VectorXd& neighbors[6]
) const {
    Eigen::VectorXd S = Eigen::VectorXd::Zero(8);

    // 计算电流密度
    double eta = eta_;  // 常数!
    Eigen::Vector3d J = compute_current_density(...);
    double J_sq = J.squaredNorm();

    // 能量耗散
    S(4) = eta * J_sq;

    // 磁场扩散
    Eigen::Vector3d lapl_B = compute_laplacian_B(...);
    S(5) = eta * lapl_B(0);
    S(6) = eta * lapl_B(1);
    S(7) = eta * lapl_B(2);

    return S;
}
```

**特点**:
- 简单直接
- η(x,y,z) 替代为 η (常数)
- 8个源项

#### 高级MHD

**电阻耗散源**:
```cpp
Eigen::VectorXd resistive_source(
    const Eigen::VectorXd& U_center,
    double x, double y, double z,  // 位置用于 η(x,y,z)
    double dx, double dy, double dz,
    const Eigen::VectorXd& neighbors[6]
) const {
    Eigen::VectorXd S = Eigen::VectorXd::Zero(9);

    // 在该位置查询电阻
    double eta = resistivity_model_(x, y, z);  // 位置相关!

    // 计算电流密度
    Eigen::Vector3d J = compute_current_density(...);
    double J_sq = J.squaredNorm();

    // 能量耗散: η·J²
    S(4) = eta * J_sq;

    // 磁场扩散: η·∇²B
    Eigen::Vector3d lapl_B = compute_laplacian_B(...);
    S(5) = eta * lapl_B(0);
    S(6) = eta * lapl_B(1);
    S(7) = eta * lapl_B(2);

    // GLM源: -ch·∇·B - (cr/ch)·ψ
    double div_B = compute_div_B(...);
    double psi = U_center(8);
    double ch = glm_params_.ch;
    double cr = glm_params_.cr;
    S(8) = -ch * div_B - (cr/ch) * psi;

    return S;
}
```

**额外助手函数**:
```cpp
Eigen::Vector3d compute_current_density(...);
double compute_div_B(...);
Eigen::Vector3d compute_laplacian_B(...);
double glm_source(...);
```

**特点**:
- η 随位置变化
- GLM 字段额外源
- 9个源项
- 完整诊断能力

---

### 5. 初始条件

#### 基础MHD

```cpp
// 需要手动指定初始状态
Eigen::VectorXd U0(8);
U0(0) = rho0;        // ρ
U0(1) = rho0 * u0;   // ρu
// ... etc

// 或使用简单方法:
uniform_field_initial(rho0, p0, Bx0, By0, Bz0, perturbation);
```

**可用**:
- ✓ 均匀磁场 + 小扰动
- ✓ 自定义初始化

**不适合**:
- ✗ Harris平衡磁场
- ✗ 磁重联特定设置

#### 高级MHD

```cpp
// 专为磁重联设计的Harris平衡
struct HarrisSheetConfig {
    double L_sheet = 1.0;     // 电流片厚度
    double n0 = 1.0;          // 参考密度
    double p0 = 0.1;          // 参考压力
    double B0 = 1.0;          // 参考磁场
    double beta = 0.2;        // 等离子体β
    double perturbation_amplitude = 0.03;  // 3%触发
};

Eigen::VectorXd U0 = mhd.harris_sheet_initial(x, y, z, harris);
```

**Harris配置包括**:

1. **磁场配置**:
   ```
   Bx(y) = B₀·tanh(y/L)  [剪切场]
   By(x,y) = A·sin(πx)·e^(-y²)  [m=1扰动]
   ```

2. **密度分布** (压力平衡):
   ```
   ρ(y) = ρ₀[1 + (1/β - 1)·sech²(y/L)]
   ```

3. **压力平衡**:
   ```
   p(y) = p₀ - 0.5B₀²tanh²(y/L)/μ₀
   ```

4. **扰动触发**:
   ```
   δBy = 0.03·B₀·sin(πx)·e^(-y²)
   ```

**优点**:
- ✓ 物理上正确的平衡
- ✓ 已知会发生磁重联
- ✓ OpenMHD 参数兼容
- ✓ 可调整参数

---

### 6. 时间步长约束

#### 基础MHD

**仅CFL条件**:
```
dt < CFL·dx / (|u| + c_fast + ch)

其中:
  c_fast = √(a² + B²/(μ₀ρ))  [快速磁声速]
  ch ≈ 0                      [无GLM]
```

**简单**:
- ✓ 单一约束
- ✓ 快速计算

#### 高级MHD

**多重约束**:

1. **双曲CFL**:
   ```
   dt < CFL·dx / (|u| + max(c_fast, ch))
   ```

2. **抛物约束** (来自扩散):
   ```
   dt < 0.25·(dx²)/(η_max·κ_max)

   其中 κ_max 是热扩散系数
   对于 η = 0.01667, dx = 0.01:
     dt < 0.25·0.0001/(0.01667) ≈ 1.5e-4
   ```

3. **GLM约束**:
   ```
   dt·ch/(cr·dx²) < 1  [GLM衰减]
   ```

**结果**:
- ✗ 多个约束需检查
- ✗ 通常由扩散约束主导
- ✗ 时间步长更小 (扩散稳定性)
- ✓ 数值稳定性得到保证

---

### 7. 计算成本对比

#### 每个单元每个时间步

**基础MHD**:
```
理想通量: ~100 FLOPs
电阻源: ~150 FLOPs (求电流、Laplacian)
总计: ~250 FLOPs/cell
```

**高级MHD**:
```
理想通量: ~100 FLOPs
电阻源 (位置相关): ~150 FLOPs
GLM源: ~50 FLOPs (额外div_B计算)
总计: ~300 FLOPs/cell
```

**差异**: +20% 计算成本

#### 内存访问

**基础MHD**:
```
状态数组: 8变量 × 8字节 = 64字节/cell
邻域: 12邻域 × 8变量 × 8字节 = 768字节
临时通量: 3方向 × 8变量 × 8字节 = 192字节
总计: ~1 KB working set
```

**高级MHD**:
```
状态数组: 9变量 × 8字节 = 72字节/cell
邻域: 12邻域 × 9变量 × 8字节 = 864字节
临时通量: 3方向 × 9变量 × 8字节 = 216字节
总计: ~1.1 KB working set
```

**差异**: +11% 内存

---

## 何时使用哪一种？

### 使用**基础MHD**如果：

1. ✅ 模拟**通用MHD问题**
   - MHD激波
   - 等离子体不稳定性
   - 磁场传播

2. ✅ 需要**快速代码执行**
   - 参数扫描
   - 初步研究
   - 方法验证

3. ✅ **学习/教育用途**
   - 理解MHD物理
   - 简单实现
   - 调试容易

4. ✅ **已有通用代码库**
   - 与现有Euler代码集成
   - 最小改动

### 使用**高级MHD**如果：

1. ✅ 模拟**磁重联**
   - 日冕加热
   - 太阳耀斑
   - 核聚变等离子体

2. ✅ 需要**精确物理**
   - 论文质量结果
   - 对比实验数据
   - OpenMHD参数化

3. ✅ 有**足够计算资源**
   - 超级计算集群
   - 可承受 ~30% 额外成本
   - GPU加速可行

4. ✅ 需要**完整诊断**
   - 能量追踪
   - 磁拓扑分析
   - 重联速率测量

---

## 实际迁移

### 从基础→高级

**代码改动最少**:

1. **替换包含**:
   ```cpp
   // #include "physics/resistive_mhd3d.hpp"
   #include "physics/resistive_mhd3d_advanced.hpp"
   using namespace fvm3d::physics;
   ```

2. **替换声明**:
   ```cpp
   // ResistiveMHD3D mhd(eta);
   AdvancedResistiveMHD3D::ResistivityModel resistivity;
   resistivity.eta0 = 1e-3;
   resistivity.eta1 = 0.01667;

   AdvancedResistiveMHD3D mhd(resistivity);
   ```

3. **通量调用** (添加位置参数):
   ```cpp
   // F = mhd.flux_x(U);
   F = mhd.flux_x(U, x, y, z);  // 添加位置
   ```

4. **源项调用** (添加位置参数):
   ```cpp
   // S = mhd.resistive_source(U, dx, dy, dz, ...);
   S = mhd.resistive_source(U, x, y, z, dx, dy, dz, ...);  // 添加位置
   ```

5. **初始条件**:
   ```cpp
   // 如果需要Harris平衡:
   U0 = mhd.harris_sheet_initial(x, y, z);
   ```

---

## 总结对比表

```
┌─────────────────────────────────────────────────────────┐
│          特性             │    基础    │    高级    │
├─────────────────────────────────────────────────────────┤
│ 守恒变量              │ 8         │ 9 (+ψ)    │
│ 电阻率模型            │ 常数 η    │ η(x,y,z) │
│ GLM约束               │ 隐含      │ 显式      │
│ 源项复杂度            │ 低        │ 中等      │
│ 初始条件灵活性        │ 中等      │ 优秀      │
│ 计算成本              │ 标准      │ +20%      │
│ 适用场景              │ 通用      │ 磁重联    │
│ 参数调整              │ 简单      │ 复杂      │
│ 验证成熟度            │ 中等      │ OpenMHD   │
│ 学习成本              │ 低        │ 中等      │
└─────────────────────────────────────────────────────────┘
```

---

## 参考资源

**基础MHD**:
- `docs/Phase4_GlobalReduction_Implementation.md`
- Evan et al. (2009) - MHD Riemann Solvers

**高级MHD**:
- `docs/AdvancedResistiveMHD_Implementation.md`
- Dedner et al. (2002) - GLM divergence cleaning
- Zenitani & Miyoshi (2011) - Magnetic reconnection
- OpenMHD Documentation

---

**推荐**: 开始用**基础MHD**进行快速原型设计和测试，当需要精确磁重联物理时升级到**高级MHD**。

