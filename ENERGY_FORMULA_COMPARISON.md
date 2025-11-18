# 中间态能量公式详细对比分析

## 问题陈述

fvm_3d的能量计算与OpenMHD存在可能的偏差，导致数值不稳定性。本文逐行对比两个实现。

---

## 中间态（U_Lstar）能量公式对比

### 1. OpenMHD实现 (flux_solver.cuf, 行 503-504)

```fortran
UL1(en) = ( ( aL-VL(i,j,k,vn) )*UL(en) - ptL*VL(i,j,k,vn) + pt*aM +
     U_hll(bn)*( vBL - aM*U_hll(bn) - vt1L*UL1(bt1) - vt2L*UL1(bt2)) ) / ( aL - aM )
```

**分解**:
```
分子 = (aL - vL_n) * E_L                           [项1]
       - ptL * vL_n                                [项2]
       + pt * aM                                   [项3]
       + B_n * (vBL - aM*B_n - vt1L*By_L* - vt2L*Bz_L*)  [项4: Poynting通量]

分母 = aL - aM                                     [切线速度]
```

其中：
- `ptL` = `p_L + 0.5*B^2_L` (左侧总压力)
- `pt` = `p_m + 0.5*B^2_hll` (HLL态总压力)
- `aM` = contact速度 (S_M)
- `vBL` = `v_L · B_L` = `u_L*B_n + v_L*By_L + w_L*Bz_L`
- `By_L*` = `UL1(bt1)` = 左侧中间态 **conserved** 磁场
- `Bz_L*` = `UL1(bt2)` = 左侧中间态 **conserved** 磁场
- `vt1L` = 左侧中间态 **primitive** 切向速度1
- `vt2L` = 左侧中间态 **primitive** 切向速度2

**关键点**: Poynting项中使用的是**已计算的中间态值** `UL1(bt1)` 和 `UL1(bt2)`。

---

### 2. fvm_3d实现 (riemann_hlld.cpp, 行 330-331)

```cpp
double E_Lstar = ((S_L - u_L) * U_L(4) - pt_L * u_L + pt_Lstar * S_M +
                  B_n * (vdotB_L - vdotB_Lstar)) / (S_L - S_M);
```

其中：
```cpp
double vdotB_L = u_L * B_n + v_L * By_L + w_L * Bz_L;        // 原始L态v·B
double vdotB_Lstar = S_M * B_n + v_Lstar * By_Lstar + w_Lstar * Bz_Lstar;
                                                               // 中间态L*的v·B
```

**分解**:
```
分子 = (S_L - u_L) * E_L                          [项1]
       - pt_L * u_L                               [项2]
       + pt_Lstar * S_M                           [项3]
       + B_n * (vdotB_L - vdotB_Lstar)           [项4: Poynting通量]

分母 = S_L - S_M                                 [切线速度]

展开项4:
B_n * (u_L*B_n + v_L*By_L + w_L*Bz_L
     - S_M*B_n - v_Lstar*By_Lstar - w_Lstar*Bz_Lstar)

= B_n * [(u_L - S_M)*B_n + (v_L*By_L - v_Lstar*By_Lstar) + (w_L*Bz_L - w_Lstar*Bz_Lstar)]
```

**关键点**: Poynting项包含了 `u_L * B_n` 和 `S_M * B_n`，这在OpenMHD中是分离处理的！

---

## 详细公式对比

### 对应关系映射

| 量 | OpenMHD | fvm_3d | 备注 |
|----|---------|--------|------|
| 左侧快波速 | aL | S_L | ✓ 相同 |
| contact速度 | aM | S_M | ✓ 相同 |
| 左侧密度 | VL(ro) | rho_L | ✓ 相同 |
| 左侧法向速度 | VL(vn) | u_L | ✓ 相同 |
| 左侧总压力 | ptL | pt_L | ✓ 相同 |
| HLL态总压力 | pt | pt_Lstar | ✓ 相同 |
| 左侧v·B | vBL | vdotB_L | ✓ 相同 |
| 法向磁场 | U_hll(bn) | B_n | ✓ 相同 |
| **关键!** 中间态By | UL1(bt1) | By_Lstar | 差异！ |
| **关键!** 中间态Bz | UL1(bt2) | Bz_Lstar | 差异！ |

---

## 公式数值验证

### 场景: 标准HLLD Riemann问题

设定:
```
左侧:   ρ=1.0, u=0, v=0, p=1.0, B=0.75
右侧:   ρ=0.125, u=0, v=0, p=0.1, B=0.75
```

计算过程:
```
1. 计算波速: S_L, S_M, S_R
2. 计算HLL态和中间态磁场: By_Lstar, Bz_Lstar
3. 计算能量修正
```

### OpenMHD Poynting项

```fortran
vBL = u_L*B_n + v_L*By_L + w_L*Bz_L  = 0 + 0 + 0 = 0

! 关键：这里使用的是已计算的中间态值 UL1(bt1)
Poynting = U_hll(bn)*( vBL - aM*U_hll(bn) - vt1L*UL1(bt1) - vt2L*UL1(bt2))
         = B_n * (0 - aM*B_n - 0 - 0)              [因为v=0]
         = B_n * (-aM*B_n)
         = -aM * B_n^2
```

### fvm_3d Poynting项

```cpp
vdotB_L = u_L*B_n + v_L*By_L + w_L*Bz_L = 0
vdotB_Lstar = S_M*B_n + v_Lstar*By_Lstar + w_Lstar*Bz_Lstar
            = S_M*B_n + 0 + 0                       [因为v_Lstar=0]

Poynting = B_n * (vdotB_L - vdotB_Lstar)
         = B_n * (0 - S_M*B_n)
         = -S_M * B_n^2                            [注意: S_M ≠ aM!]
```

**差异分析**:
- OpenMHD: `-aM * B_n^2` (其中aM = S_M = contact速度)
- fvm_3d:  `-S_M * B_n^2`

在这个简单情况下，两者应该相等（因为aM = S_M），但**fvm_3d的公式结构不同！**

---

## 关键差异识别

### 差异1: Poynting通量中的磁场来源

**OpenMHD**:
```fortran
! 显式使用中间态的conserved磁场分量
Poynting = B_n * (vBL - aM*B_n - vt1L*UL1(bt1) - vt2L*UL1(bt2))
                            ↑这里的UL1(bt1)是中间态!
```

**fvm_3d**:
```cpp
Poynting = B_n * (vdotB_L - vdotB_Lstar)
                 ↑原始态      ↑中间态

其中vdotB_Lstar使用:
= S_M*B_n + v_Lstar*By_Lstar + w_Lstar*Bz_Lstar
```

这看起来是等价的，因为：
- `v_Lstar = vt1L` (fvm_3d的标记)
- `By_Lstar = UL1(bt1)` (中间态磁场)

**但是**，让我们检查这些值的计算方式...

### 差异2: 分子中的压力项

**OpenMHD**:
```
项2: -ptL * vL_n                [原始L态]
项3: pt * aM                     [HLL态]
```

**fvm_3d**:
```
项2: -pt_L * u_L                [原始L态] ✓ 相同
项3: pt_Lstar * S_M             [??? 是pt_Lstar还是pt?]
```

**这是一个潜在的bug!** fvm_3d使用的是 `pt_Lstar` (中间态总压力)，而OpenMHD使用的是 `pt` (HLL态总压力)！

---

## 潜在的Bug

### Bug #1: 压力项使用错误的状态

**fvm_3d代码**:
```cpp
double pt_Lstar = p_m + 0.5 * B2_Lstar;  // 中间态磁压力

E_Lstar = ((S_L - u_L) * U_L(4) - pt_L * u_L + pt_Lstar * S_M +  // ← BUG!
          B_n * (vdotB_L - vdotB_Lstar)) / (S_L - S_M);
```

**应该是**:
```cpp
double p_middle = p_m;  // 中间态气压 (从HLL状态)
double B2_hll = By_hll*By_hll + Bz_hll*Bz_hll + B_n*B_n;
double pt = p_middle + 0.5 * B2_hll;  // HLL态总压力

E_Lstar = ((S_L - u_L) * U_L(4) - pt_L * u_L + pt * S_M +  // ← 应该用HLL态
          B_n * (vdotB_L - vdotB_Lstar)) / (S_L - S_M);
```

**影响**: 使用中间态磁场的压力（依赖于磁重构）而不是HLL态的固定磁场，可能导致能量不守恒。

### Bug #2: Poynting项的完整性

OpenMHD的Poynting项：
```fortran
U_hll(bn)*( vBL - aM*U_hll(bn) - vt1L*UL1(bt1) - vt2L*UL1(bt2))
```

展开：
```
= B_n * vBL - B_n^2*aM - B_n*vt1L*UL1(bt1) - B_n*vt2L*UL1(bt2)
```

而fvm_3d的Poynting项：
```cpp
B_n * (vdotB_L - vdotB_Lstar)
    = B_n*vdotB_L - B_n*vdotB_Lstar
    = B_n*(u_L*B_n + v_L*By_L + w_L*Bz_L) - B_n*(S_M*B_n + v_Lstar*By_Lstar + w_Lstar*Bz_Lstar)
    = B_n*u_L*B_n - B_n*S_M*B_n + B_n*v_L*By_L - B_n*v_Lstar*By_Lstar + ...
    = B_n^2*(u_L - S_M) + B_n*v_L*By_L - B_n*v_Lstar*By_Lstar + ...
```

这不等于OpenMHD的展开式！fvm_3d缺少了某些项或项的符号不对。

---

## 推荐修复

### 修复方案：精确映射OpenMHD公式

```cpp
// 中间态能量计算 - 严格遵循OpenMHD
Eigen::VectorXd HLLDSolver::compute_state_L(
    const Eigen::VectorXd& U_L, const Eigen::VectorXd& V_L,
    const Eigen::VectorXd& U_hll, const Eigen::VectorXd& V_hll,
    double S_L, double S_M, double p_m, double B_x,
    double rho_Lstar
) const {
    // ... 前面的代码 ...

    // 关键修复: 使用HLL态的总压力，而不是中间态的
    double B2_hll = V_hll(5)*V_hll(5) + V_hll(6)*V_hll(6) + V_hll(7)*V_hll(7);
    double pt_hll = p_m + 0.5 * B2_hll;  // ← 使用HLL态磁场

    double pt_L = p_L + 0.5 * B2_L;

    // Poynting通量: 遵循OpenMHD公式
    double vBL = u_L * B_n + v_L * By_L + w_L * Bz_L;
    double vBLstar = S_M * B_n + v_Lstar * By_Lstar + w_Lstar * Bz_Lstar;

    // OpenMHD格式的能量公式
    double E_Lstar = ((S_L - u_L) * U_L(4) - pt_L * u_L + pt_hll * S_M +
                      B_n * (vBL - vBLstar)) / (S_L - S_M);

    U_star(4) = E_Lstar;

    return U_star;
}
```

---

## 测试验证

### 验证步骤

1. **添加详细输出**
```cpp
// 在compute_state_L中
std::cout << "DEBUG Energy L*:" << std::endl;
std::cout << "  pt_L = " << pt_L << std::endl;
std::cout << "  pt_hll = " << pt_hll << std::endl;
std::cout << "  pt_Lstar = " << pt_Lstar << std::endl;  // 当前错误的值
std::cout << "  vBL = " << vBL << std::endl;
std::cout << "  vBLstar = " << vBLstar << std::endl;
std::cout << "  E_Lstar (old) = " << E_Lstar_old << std::endl;
std::cout << "  E_Lstar (new) = " << E_Lstar_new << std::endl;
```

2. **对比OpenMHD**
   - 运行相同的Brio-Wu初值条件
   - 比较早期步骤（100, 200）的能量演化
   - 验证能量是否趋势正确

3. **预期结果**
   - ✓ 能量增长曲线应该更平缓
   - ✓ 时间步长应该保持合理
   - ✓ 能够到达最终时间

---

## 总结

### 识别的Bug

| Bug | 严重性 | 位置 | 影响 |
|-----|--------|------|------|
| **压力项** | 🔴 高 | 行 330 | 能量不守恒 |
| **Poynting项** | 🟡 中 | 行 331 | 能量计算不准确 |
| 密度选择 | 🟡 中 | 行 451-455 | 熵条件违反 |

### 下一步

1. 立即修复压力项（使用pt_hll而不是pt_Lstar）
2. 验证Poynting项的完整性
3. 重新测试Brio-Wu和Orszag-Tang
4. 与OpenMHD对标

---

*分析日期: 2025-11-18*
*对标版本: OpenMHD 20250804, fvm_3d commit b53c6a2*
