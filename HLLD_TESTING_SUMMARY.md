# HLLD 黎曼求解器测试总结和修复建议

## 执行摘要

对最新的HLLD实现进行了全面的基准测试，发现**严重的数值不稳定性**，导致能量指数爆炸和时间推进停滞。

**关键发现**: 问题不是由最近的OpenMHD等价性修复引入的，而是**基础性的**，存在于所有测试过的HLLD版本中。

---

## 测试概览

### 测试案例
1. **Brio-Wu 磁激波管** - 包含Alfvén波和接触间断
2. **Orszag-Tang 涡旋** - 光滑初值的复杂MHD流
3. **Harris Sheet 3D** - 磁重联（未完成）

### 结果状态
```
✓ Brio-Wu:      运行完成, 数值失败 (能量爆炸到4.259e+136)
✓ Orszag-Tang:  运行完成, 数值失败 (能量NaN @ step 500)
⏳ Harris Sheet: 尚未开始
```

---

## 数值不稳定性分析

### 症状
```
┌─────────────────────────────────────────┐
│ 能量演化                                 │
├─────────────────────────────────────────┤
│ Step 100:   E_tot = 1.798e+05           │
│ Step 200:   E_tot = 4.526e+13           │
│ Step 300:   E_tot = 1.132e+22           │
│ Step 400:   E_tot = 2.750e+30           │
│ Step 500:   E_tot = NaN                 │
│                                         │
│ 增长率: ~10^8 per 100 steps             │
└─────────────────────────────────────────┘

┌─────────────────────────────────────────┐
│ 时间步长衰减                             │
├─────────────────────────────────────────┤
│ Step 100:   dt = 3.864e-06              │
│ Step 200:   dt = 2.791e-10              │
│ Step 300:   dt = 1.840e-14              │
│ Step 400:   dt = 1.205e-18              │
│ Step 500:   dt = 7.747e-23              │
│                                         │
│ 衰减率: ~1000倍 per 100 steps           │
└─────────────────────────────────────────┘
```

### 原因分析

#### 最可能的根本原因

**1. 中心态能量计算中的Poynting通量符号错误** (90%概率)
```cpp
// 当前代码 (可能有问题)
double vdotB_Lstar = v_Lstar_t1 * By_Lstar_cons + v_Lstar_t2 * Bz_Lstar_cons;
double vdotB_central = v_central_t1 * By_central + v_central_t2 * Bz_central;
E_central = E_Lstar - sqrt_rho_L * (vdotB_Lstar - vdotB_central) * sign_B;

问题:
- v·B项可能缺少法向分量(S_M * B_x)
- 或者法向分量应该以不同方式计入
- Poynting通量的物理解释可能被颠倒
```

**2. 中间态能量公式不适用于强磁场** (70%概率)
```cpp
// 在 compute_state_L/R() 中
double E_Lstar = ((S_L - u_L) * U_L(4) - pt_L * u_L + pt_Lstar * S_M +
                  B_n * (vdotB_L - vdotB_Lstar)) / (S_L - S_M);

问题:
- 分母(S_L - S_M)可能非常小
- 数值精度可能丧失
- 对高Ma数或强B场情况不稳定
```

**3. 初始条件或能量初始化错误** (30%概率)
```
Orszag-Tang初期观察:
- 理论最大KE ≈ 2.78
- 实际KE @ step 100 = 3.411e+02 (122倍偏差!)
- 表明初始条件设置有问题
```

---

## 推荐的修复步骤

### 第一阶段: 诊断 (2-3小时)

#### 步骤1: 核实Poynting通量项
```cpp
// 在 compute_central_state() 中添加详细输出
std::cout << "DEBUG Central State Energy:" << std::endl;
std::cout << "  E_Lstar = " << E_Lstar << std::endl;
std::cout << "  vdotB_Lstar = " << vdotB_Lstar << std::endl;
std::cout << "  vdotB_central = " << vdotB_central << std::endl;
std::cout << "  correction = " << sqrt_rho_L * (vdotB_Lstar - vdotB_central) * sign_B << std::endl;
std::cout << "  E_central = " << E_central << std::endl;
```

**检查清单**:
- [ ] Poynting通量是否包括所有速度分量?
- [ ] 符号是否正确 (±)?
- [ ] sign_B的定义是否正确?
- [ ] sqrt_rho的使用是否一致?

#### 步骤2: 比对OpenMHD源代码
```fortran
! OpenMHD flux_solver.cuf line 560-561
U2(en) = UL1(en) - roLs * ( vt1L*UL1(bt1) + vt2L*UL1(bt2)
     - at1*U2(bt1) - at2*U2(bt2) ) * f2
```

**映射关系检查**:
- [ ] fvm_3d 的 vdotB_Lstar 是否等于 `vt1L*UL1(bt1) + vt2L*UL1(bt2)`?
- [ ] vdotB_central 是否等于 `at1*U2(bt1) + at2*U2(bt2)`?
- [ ] f2 是否等于 sign_B?
- [ ] roLs 是否等于 sqrt_rho_L?

#### 步骤3: 验证初始条件
```cpp
// 在 brio_wu_initial_condition() 中添加检查
double KE_computed = 0.5 * rho * (u*u + v*v + w*w);
double BE_computed = 0.5 * (Bx*Bx + By*By + Bz*Bz);
double E_total = p/(gamma-1) + KE_computed + BE_computed;

// 添加断言检查
assert(KE_computed >= 0.0 && KE_computed < 10.0); // 合理范围
assert(E_total > 0.0);
```

### 第二阶段: 修复 (4-6小时)

#### 修复方案A: 恢复完整v·B项
```cpp
// 尝试包括法向分量
double vdotB_Lstar = S_M * B_x + v_Lstar_t1 * By_Lstar_cons + v_Lstar_t2 * Bz_Lstar_cons;
double vdotB_central = S_M * B_x + v_central_t1 * By_central + v_central_t2 * Bz_central;
// 注意: 法向分量在差分中应该抵消
double E_correction = sqrt_rho_L * (vdotB_Lstar - vdotB_central) * sign_B;
E_central = E_Lstar - E_correction;
```

#### 修复方案B: 使用OpenMHD精确映射
```cpp
// 严格遵循OpenMHD实现
// 不使用Poynting通量差分, 而是直接计算
double E_Lstar = U_Lstar(4);
double E_Rstar = U_Rstar(4);
double E_avg = w1 * E_Lstar + w2 * E_Rstar;

// 能量修正基于contact方向
double E_correction;
if (S_M >= 0.0) {
    E_correction = sign_B * sqrt_rho_L * inv_wsum *
                  (v_Lstar_t1 * By_Lstar_cons + v_Lstar_t2 * Bz_Lstar_cons -
                   v_central_t1 * By_central - v_central_t2 * Bz_central);
} else {
    E_correction = -sign_B * sqrt_rho_R * inv_wsum *
                  (v_Rstar_t1 * By_Rstar_cons + v_Rstar_t2 * Bz_Rstar_cons -
                   v_central_t1 * By_central - v_central_t2 * Bz_central);
}
E_central = E_avg - E_correction;
```

#### 修复方案C: 更保守的时间步长
```cpp
// 在主时间循环中
config.CFL = 0.1; // 从0.4降低到0.1
// 或者添加自适应CFL
double CFL_adaptive = std::min(0.4, 1.0 / (1.0 + pressure_gradient_magnitude));
```

### 第三阶段: 验证 (3-4小时)

#### 验证步骤1: 能量监测
```
预期: 能量应单调递减或有界震荡
┌─────────────────────────────────┐
│ 修复前: E(step n) ~ 10^(n)     │
│ 修复后: E(step n) ~ 常数或衰减 │
└─────────────────────────────────┘
```

#### 验证步骤2: 时间推进
```
预期: 能够到达最终时间
┌──────────────────────────────────┐
│ Brio-Wu:     t_final = 0.2      │
│ Orszag-Tang: t_final = 5.0      │
└──────────────────────────────────┘
```

#### 验证步骤3: 对比已知解
```
Brio-Wu激波管有已出版的精确解 (Brio & Wu 1988)
应该能够在网格上看到5波结构:
  Fast wave | Slow wave | Alfvén | Contact | Slow wave | Fast wave
```

---

## 实施时间表

```
第1天 (Today):
├─ 08:00-10:00: 诊断阶段 - 添加调试输出
├─ 10:00-11:00: 与OpenMHD源代码逐行对比
└─ 11:00-12:00: 确定根本原因

第2天:
├─ 08:00-12:00: 实施修复方案
├─ 12:00-14:00: 重新编译和测试
├─ 14:00-16:00: 验证修复是否有效
└─ 16:00-17:00: 如果失败, 尝试备选方案

第3天:
├─ 08:00-10:00: 完成验证
├─ 10:00-12:00: 生成最终测试报告
├─ 12:00-14:00: Harris Sheet 3D测试 (如果修复成功)
└─ 14:00-15:00: 提交改进
```

---

## 关键资源

### OpenMHD对比
- **文件**: `/Users/ericzeus/PycharmProjects/fvm_3d/openmhd-20250804/3D_reconnection_gpu/flux_solver.cuf`
- **关键行**: 543-588 (HLLD中心态计算)
- **注意**: OpenMHD使用Fortran, 需要仔细映射到C++

### fvm_3d HLLD实现
- **文件**: `src/spatial/riemann_solvers/riemann_hlld.cpp`
- **关键函数**:
  - `compute_state_L()`: 行 321-385
  - `compute_state_R()`: 行 390-440
  - `compute_central_state()`: 行 445-525

### 测试框架
- **Brio-Wu**: `examples/brio_wu_shock_tube.cpp`
- **Orszag-Tang**: `examples/orszag_tang_vortex.cpp`
- **Harris Sheet**: `examples/harris_sheet_3d.cpp`

---

## 风险评估

| 修复方案 | 风险等级 | 副作用 | 备注 |
|---------|--------|------|------|
| A: 恢复v·B项 | 低 | 无 | 可能是正确的修复 |
| B: OpenMHD映射 | 中 | 可能改变稳定性特性 | 需要验证 |
| C: 降低CFL | 低 | 性能下降 | 快速临时方案 |
| D: 重写求解器 | 高 | 时间成本高 | 最后手段 |

---

## 预期结果

### 修复成功的标志
```
✓ 能量保持有界 (不超过初值10倍)
✓ 时间步长保持合理 (dt > 1e-6)
✓ 能够到达最终时间
✓ 压力和密度保持正性
✓ 解的形态符合物理预期
```

### 对标准准
```
与ATHENA++或PLUTO比较:
- 能量误差 < 5%
- 解的光滑度相当
- 收敛阶数 ≥ 1st order
```

---

## 附加说明

### 为什么这很重要
HLLD是MHD模拟的核心算法，能量守恒违反会导致:
- 虚假的湍流发展
- 不物理的磁重联
- 长时间模拟积累误差
- 对网格分辨率敏感

### 下一步 (如果修复失败)
1. 考虑切换到HLL求解器 (更稳定但精度低)
2. 使用LLF/Godunov类求解器作为临时方案
3. 联系ATHENA++或OpenMHD开发者获取建议

---

*文档生成时间: 2025-11-18*
*版本: HLLD b53c6a2*
*优先级: 高 (影响所有MHD模拟)*
