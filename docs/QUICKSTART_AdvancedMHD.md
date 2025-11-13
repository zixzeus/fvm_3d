# 快速开始：高级磁重联MHD求解器

**OpenMHD风格的FVM3D实现**

---

## 30秒概览

```cpp
#include "physics/resistive_mhd3d_advanced.hpp"
using namespace fvm3d::physics;

// 1. 配置电阻率模型 (位置相关)
AdvancedResistiveMHD3D::ResistivityModel resistivity;
resistivity.eta0 = 1e-3;    // 背景: Rm0=1000
resistivity.eta1 = 0.01667; // 增强: Rm1=60
resistivity.localization_scale = 1.0;

// 2. 配置GLM约束
AdvancedResistiveMHD3D::GLMParameters glm(0.2, 0.2);

// 3. 创建求解器
AdvancedResistiveMHD3D mhd(resistivity, glm);

// 4. 初始化Harris平衡磁场
AdvancedResistiveMHD3D::HarrisSheetConfig harris;
harris.L_sheet = 1.0;
harris.beta = 0.2;
harris.perturbation_amplitude = 0.03;

Eigen::VectorXd U = mhd.harris_sheet_initial(x, y, z, harris);
```

---

## 关键参数速查表

### 电阻率配置

```cpp
ResistivityModel resistivity;
resistivity.eta0 = 1e-3;              // 背景电阻, Rm₀ = 1000
resistivity.eta1 = 0.01667;           // 增强电阻, Rm₁ = 60
resistivity.localization_scale = 1.0; // 转换宽度
```

**含义**:
```
η(x,y) = η₀ + (η₁ - η₀)·sech²(√(x²+y²))
        ↓
原点强烈加热 → 远处接近理想MHD
```

### GLM参数

```cpp
GLMParameters glm;
glm.ch = 0.2;  // 波速 (0.2 × max MHD速度)
glm.cr = 0.2;  // 衰减系数 (e^(-t/1) per timestep)
```

**含义**: ψ 字段错误在0.4-1.0时间步内衰减到零

### Harris片配置

```cpp
HarrisSheetConfig harris;
harris.L_sheet = 1.0;              // 电流片厚度
harris.n0 = 1.0;                   // 参考密度
harris.p0 = 0.1;                   // 参考压力
harris.B0 = 1.0;                   // 参考磁场
harris.beta = 0.2;                 // 等离子体 β (磁场主导)
harris.perturbation_amplitude = 0.03; // 3% 触发
```

---

## 物理方程速查

### 守恒变量 (9个)

| 索引 | 名称 | 含义 | 单位 |
|------|------|------|------|
| 0 | ρ | 质量密度 | - |
| 1-3 | ρu,ρv,ρw | 动量 | - |
| 4 | E | 总能量 | - |
| 5-7 | Bx,By,Bz | 磁场 | - |
| 8 | ψ | GLM约束标量 | - |

### 能量分解

```
E = ρe_int + ½ρ(u²+v²+w²) + ½(Bx²+By²+Bz²)/μ₀
    └─────────┬──────────┘  └──────────────┬──────────────┘
      内能         动能                 磁能
```

### 关键物理方程

**位置相关电阻**:
```
η(x,y,z) = η₀ + (η₁-η₀)·sech²(r/L)
```

**Ohmic加热**:
```
S_E = η·J²  where J = ∇×B
```

**磁场扩散**:
```
S_B = η·∇²B
```

**GLM约束维持**:
```
∂ψ/∂t = -ch·∇·B - (cr/ch)·ψ
```

---

## 通用工作流

### 1. 求解器初始化

```cpp
// 配置
AdvancedResistiveMHD3D::ResistivityModel res;
res.eta0 = 1e-3;
res.eta1 = 0.01667;
AdvancedResistiveMHD3D::GLMParameters glm(0.2, 0.2);

// 创建
AdvancedResistiveMHD3D mhd(res, glm);
```

### 2. 初始化

```cpp
// Harris平衡
AdvancedResistiveMHD3D::HarrisSheetConfig harris;
Eigen::VectorXd U = mhd.harris_sheet_initial(x, y, z, harris);

// 或通用
U = mhd.uniform_field_initial(rho0, p0, Bx0, By0, Bz0, pert);
```

### 3. 时间步长计算

```cpp
// 双曲CFL
double c_fast = mhd.fast_speed(rho, p, Bx, By, Bz);
double u_max = std::sqrt(u*u + v*v + w*w);
double wave_speed = u_max + std::max(c_fast, glm.ch);
double dt_hyperbolic = CFL * dx / wave_speed;

// 抛物扩散约束
double eta_max = res.eta1;  // 最强
double dt_parabolic = 0.25 * (dx*dx) / eta_max;

// 使用最严格
dt = std::min(dt_hyperbolic, dt_parabolic);
```

### 4. 通量计算

```cpp
Eigen::VectorXd F = mhd.flux_x(U, x, y, z);
Eigen::VectorXd G = mhd.flux_y(U, x, y, z);
Eigen::VectorXd H = mhd.flux_z(U, x, y, z);
```

### 5. 源项计算

```cpp
// 包括电阻和GLM
Eigen::VectorXd S = mhd.resistive_source(
    U_center, x, y, z, dx, dy, dz,
    U_xm, U_xp, U_ym, U_yp, U_zm, U_zp
);

// S包含: S_E = η·J²,  S_B = η·∇²B,  S_ψ = GLM
```

### 6. 状态更新 (显式Euler)

```cpp
Eigen::VectorXd dUdt_conv = (F_xp - F_xm)/dx +
                             (G_yp - G_ym)/dy +
                             (H_zp - H_zm)/dz;
U_new = U + dt * (dUdt_conv + S/rho);
```

### 7. 诊断

```cpp
// 检查有效性
if (!mhd.is_valid_state(U)) {
    // 处理故障
}

// 能量诊断
Eigen::Vector4d E = mhd.compute_energies(U);
double kinetic = E(0);
double magnetic = E(1);
double internal = E(2);
double p_total = E(3);
```

---

## 常见设置

### 标准磁重联 (推荐开始)

```cpp
// 电阻率
ResistivityModel res;
res.eta0 = 1e-3;      // 弱背景
res.eta1 = 0.01667;   // Rm₁=60强化

// GLM
GLMParameters glm(0.2, 0.2);

// Harris配置
HarrisSheetConfig harris;
harris.L_sheet = 1.0;
harris.beta = 0.2;
harris.perturbation_amplitude = 0.03;

// 参数
CFL = 0.35;
dt_resistive_scale = 0.25;  // 抛物稳定性
```

### 高Rm重联 (快速重联)

```cpp
// 增加Rm₁
res.eta1 = 0.005;  // Rm₁ ≈ 200 (更快重联)

// 缩小局域化区域
res.localization_scale = 0.5;  // 更尖锐的梯度
```

### 弱重联 (检验稳定性)

```cpp
// 减少扰动
harris.perturbation_amplitude = 0.01;  // 1% instead of 3%

// 增加背景电阻
res.eta0 = 0.005;  // 增加背景耗散
```

---

## 故障排除

### 问题: 数值爆炸

**可能原因**:
- dt太大 → 检查抛物约束: `dt < 0.25·(dx²)/η_max`
- 密度/压力变为负 → 添加floor: `rho = max(rho, 1e-10)`

**解决**:
```cpp
// 验证状态有效性
if (!mhd.is_valid_state(U)) {
    // 应用floor
    U(0) = std::max(U(0), 1e-10);  // ρ
    // 重新计算压力
}
```

### 问题: ∇·B很大

**可能原因**: 源项计算有误或GLM参数不当

**诊断**:
```cpp
double div_B = mhd.compute_div_B(...);
std::cout << "Divergence B: " << div_B << std::endl;

if (div_B > 0.01) {
    // 增加GLM波速
    glm.ch = 0.4;  // 加倍
    glm.cr = 0.4;
}
```

### 问题: 重联不发生

**可能原因**:
- 扰动幅度太小
- Rm太低
- Harris平衡有误

**检查**:
```cpp
// 验证Harris配置
double B_at_y0 = harris.B0 * std::tanh(0.0);  // = 0, 应该!
double By_pert = harris.perturbation_amplitude * std::sin(M_PI*0.5) * 1.0;
// 应该 ≈ 0.03

// 增加扰动
harris.perturbation_amplitude = 0.05;
```

---

## 性能优化

### 内存访问

```cpp
// ✓ 好: 顺序访问
for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
        for(int k=0; k<nz; k++) {
            double rho = U(0, i, j, k);  // 连续
        }
    }
}

// ✗ 差: 跳跃访问
for(int v=0; v<9; v++) {
    for(int i=0; i<nx; i++) {
        // ...
    }
}
```

### 计算成本

```
每单元每时间步:
  理想MHD通量: ~100 FLOPs
  电阻源 (η·J² + η·∇²B): ~150 FLOPs
  GLM源: ~50 FLOPs
  ────────────────────────
  总计: ~300 FLOPs/cell

与基础MHD相比: +20% 成本
```

---

## 有用的宏和常数

```cpp
// 物理常数 (内置)
static constexpr double GAMMA = 5.0/3.0;  // 单原子气体
static constexpr double MU0 = 1.0;        // 归一化单位
static constexpr int nvars = 9;           // 变量数

// 数值stability
static constexpr double RHO_FLOOR = 1e-10;
static constexpr double P_FLOOR = 1e-11;
static constexpr double B_FLOOR = 1e-12;

// 典型值
double CFL = 0.35;            // Courant数
double dt_res_factor = 0.25;  // 抛物稳定性
```

---

## 验证清单

在运行生产模拟前:

- [ ] 验证Harris平衡: `∇·B ≈ 0` 在初始时
- [ ] 检查能量守恒: `dE_total/dt ≈ dissipation`
- [ ] 确认重联增长: `By_max`单调增加
- [ ] 验证CFL条件: `max(|u|+c_f) * dt/dx < 0.35`
- [ ] 检查GLM衰减: `ψ`指数衰减
- [ ] 对比OpenMHD: 能量分布相似

---

## 进一步阅读

1. **详细文档**: `docs/AdvancedResistiveMHD_Implementation.md`
2. **对比指南**: `docs/BasicVsAdvancedMHD.md`
3. **参考论文**:
   - Dedner et al. 2002 (GLM)
   - Zenitani & Miyoshi 2011 (磁重联)
   - Miyoshi & Kusano 2005 (HLLD)

---

**准备开始?** 复制`harris_sheet_initial()`调用并开始时间步进!

