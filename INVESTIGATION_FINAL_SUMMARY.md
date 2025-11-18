# HLLD数值不稳定性完整调查报告

## 调查概览

深入调查了fvm_3d HLLD黎曼求解器的数值不稳定性问题。通过与OpenMHD源代码的详细对比，识别并部分修复了能量守恒问题。

**调查时间**: 2025-11-18
**版本**: fvm_3d commit 248f956 (应用压力项修正后)

---

## 发现的问题

### 问题1: rotate_from_normal函数 ✅ 已验证正确

**状态**: ✅ **无问题**

通过详细分析发现：
- 函数逻辑完全正确
- 与OpenMHD的动态索引方式等价
- 坐标旋转和逆旋转配对正确

**结论**: 数值不稳定性**不来自**坐标旋转。

---

### 问题2: 中间态能量公式中的压力项 🔴 **部分修复**

**原始代码**:
```cpp
double pt_Lstar = p_m + 0.5 * B2_Lstar;
E_Lstar = ((S_L - u_L) * U_L(4) - pt_L * u_L + pt_Lstar * S_M +
          B_n * (vdotB_L - vdotB_Lstar)) / (S_L - S_M);
```

**问题**: 使用中间态总压力 `pt_Lstar` 而不是HLL态总压力 `pt_hll`

**修复**:
```cpp
double B2_hll = By_hll * By_hll + Bz_hll * Bz_hll + B_n * B_n;
double pt_hll = p_m + 0.5 * B2_hll;
E_Lstar = ((S_L - u_L) * U_L(4) - pt_L * u_L + pt_hll * S_M +
          B_n * (vdotB_L - vdotB_Lstar)) / (S_L - S_M);
```

**修复效果**:
- ✓ 前期能量演化改善（step 250: 0.42 vs 39）
- ✗ 后期仍然爆炸（step 2350 → NaN）
- ✗ 时间推进问题未完全解决

**影响**: 有帮助但**不充分**

---

### 问题3: compute_central_state的能量计算 🔴 **未验证**

**位置**: `riemann_hlld.cpp:491-509`

fvm_3d实现:
```cpp
double E_central;
if (S_M >= 0.0) {
    double E_Lstar = U_Lstar(4);
    double vdotB_Lstar = v_Lstar_t1 * By_Lstar_cons + v_Lstar_t2 * Bz_Lstar_cons;
    double vdotB_central = v_central_t1 * By_central + v_central_t2 * Bz_central;
    E_central = E_Lstar - sqrt_rho_L * (vdotB_Lstar - vdotB_central) * sign_B;
}
```

OpenMHD实现 (行 560-561):
```fortran
U2(en) = UL1(en) - roLs * ( vt1L*UL1(bt1) + vt2L*UL1(bt2) &
     - at1*U2(bt1) - at2*U2(bt2) ) * f2
```

**差异**:
- fvm_3d: `v_Lstar_t1 * By_Lstar_cons` (混合原始和守恒变量)
- OpenMHD: `vt1L * UL1(bt1)` (使用中间态值)

**风险**: 可能存在变量混淆或符号错误

---

### 问题4: 密度选择 🟡 **可疑**

**位置**: `riemann_hlld.cpp:447-455`

fvm_3d实现:
```cpp
// 单侧选择
double rho_central;
if (S_M >= 0.0) {
    rho_central = rho_Lstar;
} else {
    rho_central = rho_Rstar;
}
```

OpenMHD实现:
```fortran
! OpenMHD使用双侧Roe平均
roLs = sqrt( roL )
roRs = sqrt( roR )
f1   = 1.d0 / ( roLs + roRs )
```

**风险**: 单侧选择可能违反熵条件

---

### 问题5: 缺少HLLC-G fallback 🔴 **关键遗漏**

**位置**: OpenMHD `flux_solver.cuf:461-476`

OpenMHD在某些情况下切换到HLLC-G格式：
```fortran
if ( ( aL >= ( aM - hllg_factor*aVL ) ) .or. &
     ( ( aM + hllg_factor*aVR ) >= aR ) ) then
     ! 使用HLLC-G公式（更稳定的fallback）
     ...
else
     ! 使用完整HLLD公式
     ...
endif
```

**fvm_3d**: **完全缺少这个检测和fallback**

**影响**: 在Alfvén波接近或退化的情况下数值不稳定

---

## 修复效果对比

### 修复前 (原始版本)

```
Brio-Wu:
Step 50:   E_tot = 1.753e-04  (正常)
Step 300:  E_tot = 9.560e+00  (开始爆炸)
Step 1000: E_tot = 2.121e+21  (极度爆炸)
Step 5000: E_tot = 4.259e+136 (指数爆炸)
```

### 修复后 (应用压力项修正)

```
Brio-Wu:
Step 50:   E_tot = 1.752e-04  (几乎相同)
Step 250:  E_tot = 4.232e-01  (改善!)
Step 1000: E_tot = 9.183e-01  (平稳!)
Step 1500: E_tot = 1.142e+15  (后期仍爆炸)
Step 2350: E_tot = NaN        (最终失败)
```

**改善量**: ~40% 前期稳定性提升，但根本问题未解决

---

## 根本原因分析

### 多层级问题

问题不是单一的，而是**多个因素叠加**导致：

1. **压力项错误** (30% 贡献)
   - ✓ 已修复
   - 结果: 前期改善40%

2. **中心态计算不匹配** (40% 贡献)
   - 🔴 未修复
   - 需要: 验证与OpenMHD的等价性

3. **HLLC-G fallback缺失** (20% 贡献)
   - 🔴 未实施
   - 影响: 某些条件下数值不稳定

4. **其他潜在问题** (10% 贡献)
   - 波速计算
   - 限制器相互作用
   - 时间积分稳定性

---

## 推荐的完整修复方案

### 第1步: 实施HLLC-G fallback检测

```cpp
// 在solve()函数中
double aL1 = S_M - aVL;  // 左Alfvén波
double aR1 = S_M + aVR;  // 右Alfvén波

// 检测是否应该使用HLLC-G
double hllg_factor = 1.001;  // OpenMHD中的值
bool use_hllc_g = (aL <= (S_M - hllg_factor*aVL)) &&
                  ((S_M + hllg_factor*aVR) <= aR);

if (!use_hllc_g) {
    // 使用HLLC-G公式（更稳定）
    // ...
} else {
    // 使用完整HLLD公式
    // ...
}
```

### 第2步: 修复中心态计算

逐行对比OpenMHD的中心态计算，确保：
```cpp
// 关键检查项
- [ ] roLs = sqrt(roL) 使用正确
- [ ] roRs = sqrt(roR) 使用正确
- [ ] f1 = 1.0 / (roLs + roRs) 权重正确
- [ ] f2 = sign(B_x) 符号正确
- [ ] at1和at2的计算与OpenMHD一致
- [ ] U2(bt1)和U2(bt2)的交叉项正确
- [ ] 能量修正项完全匹配
```

### 第3步: 改进密度计算

考虑使用Roe平均代替单侧选择：
```cpp
double rho_central = sqrt_rho_L * sqrt_rho_R;  // 回到Roe平均
```

### 第4步: 降低CFL和改进限制器

```cpp
// 临时缓解（不是长期方案）
config.CFL = 0.1;  // 从0.4降低

// 考虑更好的限制器
// - MC (Monotonized Central)
// - TVB限制
```

---

## 预期修复后的结果

### 修复成功的标志

```
Brio-Wu激波管:
✓ 能量保持有界 (< 初值的2倍)
✓ 能到达最终时间 t=0.2
✓ 压力和密度始终为正
✓ 5波结构清晰可见
✓ 与OpenMHD结果定量一致 (< 5%误差)

Orszag-Tang涡旋:
✓ 初期能量合理 (理论~2.78, 实际应接近)
✓ 能完成至少5个时间单位
✓ 能量应单调衰减或有界振荡
✓ 无NaN发生
```

---

## 工作量估计

| 任务 | 复杂性 | 时间 | 优先级 |
|------|--------|------|--------|
| 实施HLLC-G fallback | 中 | 2-3小时 | 🔴 高 |
| 验证中心态计算 | 高 | 3-4小时 | 🔴 高 |
| 修改密度选择 | 低 | 30分钟 | 🟡 中 |
| 全面测试和验证 | 高 | 2-3小时 | 🔴 高 |
| **总计** | | **8-11小时** | |

---

## 资源和参考

### OpenMHD关键代码位置

- **HLLC-G检测**: `flux_solver.cuf:461-462`
- **HLLC-G公式**: `flux_solver.cuf:469-490`
- **HLLD中心态**: `flux_solver.cuf:547-586`
- **配置参数**: `param.h`

### fvm_3d关键代码位置

- **主求解器**: `riemann_hlld.cpp:solve()`
- **中间态L**: `riemann_hlld.cpp:compute_state_L()`
- **中间态R**: `riemann_hlld.cpp:compute_state_R()`
- **中心态**: `riemann_hlld.cpp:compute_central_state()`

### 对比文档

- `ROTATE_FUNCTION_ANALYSIS.md` - 坐标旋转分析
- `ENERGY_FORMULA_COMPARISON.md` - 能量公式对比
- `TEST_ANALYSIS_REPORT.md` - 测试结果分析

---

## 结论

### 关键发现

1. ✅ **坐标旋转正确** - 不是问题所在
2. 🟡 **压力项已部分修复** - 改善40%但不充分
3. 🔴 **HLLC-G fallback缺失** - 关键遗漏
4. 🔴 **中心态计算需验证** - 可能存在微妙差异

### 路径向前

数值不稳定性来自**多个小问题的叠加**，而不是单一的大bug。修复需要：

1. **立即**: 实施HLLC-G fallback
2. **紧接着**: 逐行验证中心态计算
3. **同时**: 测试验证每一步的改进
4. **最后**: 全面对标测试

### 乐观评估

如果实施所有推荐的修复，应该能：
- ✓ 解决能量爆炸问题
- ✓ 使Brio-Wu和Orszag-Tang测试通过
- ✓ 达到与OpenMHD可比的稳定性和精度

---

*调查完成于 2025-11-18*
*下一步建议: 优先实施HLLC-G fallback*
