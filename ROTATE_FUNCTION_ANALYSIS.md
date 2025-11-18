# rotate_from_normal 函数对比分析

## 问题陈述

比较fvm_3d的`rotate_from_normal`函数与OpenMHD的实现，确认是否存在bug。

---

## 变量索引对比

### OpenMHD (CUDA Fortran, 基于param.h)

OpenMHD使用**Fortran 1-based indexing**:
```fortran
integer, parameter :: mx = 1, vx = 1  ! 动量x / 速度x
integer, parameter :: my = 2, vy = 2  ! 动量y / 速度y
integer, parameter :: mz = 3, vz = 3  ! 动量z / 速度z
integer, parameter :: en = 4, pr = 4  ! 总能量 / 压力
integer, parameter :: ro = 5          ! 密度
integer, parameter :: bx = 6          ! 磁场x
integer, parameter :: by = 7          ! 磁场y
integer, parameter :: bz = 8          ! 磁场z
integer, parameter :: ps = 9          ! 散度清理标量

! 变量顺序: [mx, my, mz, en, ro, bx, by, bz, ps]
!        1   2   3   4   5   6   7   8   9  (Fortran索引)
```

### fvm_3d (C++, 基于physics类)

fvm_3d使用**C++ 0-based indexing**:
```cpp
// AdvancedResistiveMHD3D (9个变量)
// 变量顺序: [ρ, ρu, ρv, ρw, E, Bx, By, Bz, ψ]
//          0   1   2   3  4  5   6   7  8  (C++索引)

U(0) = ρ (density)
U(1) = ρu (momentum x)
U(2) = ρv (momentum y)
U(3) = ρw (momentum z)
U(4) = E (total energy)
U(5) = Bx (magnetic field x)
U(6) = By (magnetic field y)
U(7) = Bz (magnetic field z)
U(8) = ψ (GLM psi)
```

### 索引映射

```
fvm_3d C++ Index  ←→  OpenMHD Fortran Index  ←→  Physical Quantity
─────────────────────────────────────────────────────────────────
      0           ←→        5               ←→  Density (ρ)
      1           ←→       mx=1             ←→  Momentum x (ρu)
      2           ←→       my=2             ←→  Momentum y (ρv)
      3           ←→       mz=3             ←→  Momentum z (ρw)
      4           ←→       en=4             ←→  Total Energy (E)
      5           ←→       bx=6             ←→  Magnetic field x
      6           ←→       by=7             ←→  Magnetic field y
      7           ←→       bz=8             ←→  Magnetic field z
      8           ←→       ps=9             ←→  GLM scalar ψ
```

**关键发现**: 顺序完全不同！OpenMHD把密度放在最后(index 5)，而fvm_3d把它放在最前(index 0)。

---

## OpenMHD的方向处理方式

OpenMHD**不使用旋转函数**，而是使用**动态索引化**：

```fortran
! 根据方向动态计算索引
vn  = vx + mod(dir-1,3)    ! 法向速度索引
vt1 = vx + mod(dir  ,3)    ! 切向速度1索引
vt2 = vx + mod(dir+1,3)    ! 切向速度2索引

bn  = bx + mod(dir-1,3)    ! 法向磁场索引
bt1 = bx + mod(dir  ,3)    ! 切向磁场1索引
bt2 = bx + mod(dir+1,3)    ! 切向磁场2索引

mn  = mx + mod(dir-1,3)    ! 法向动量索引
mt1 = mx + mod(dir  ,3)    ! 切向动量1索引
mt2 = mx + mod(dir+1,3)    ! 切向动量2索引
```

### 方向索引计算示例

| 方向 | dir值 | vn计算 | vt1计算 | vt2计算 | 结果 |
|------|------|--------|---------|---------|------|
| **X方向** | 1 | 1+mod(0,3)=1 | 1+mod(1,3)=2 | 1+mod(2,3)=3 | vn=vx, vt1=vy, vt2=vz |
| **Y方向** | 2 | 1+mod(1,3)=2 | 1+mod(2,3)=3 | 1+mod(3,3)=1 | vn=vy, vt1=vz, vt2=vx |
| **Z方向** | 3 | 1+mod(2,3)=3 | 1+mod(3,3)=1 | 1+mod(4,3)=2 | vn=vz, vt1=vx, vt2=vy |

这个巧妙的索引化方式**避免了状态旋转**，通过直接访问正确的分量实现了对所有方向的支持。

---

## fvm_3d的旋转方式

fvm_3d使用**状态交换**来处理方向：

```cpp
Eigen::VectorXd HLLDSolver::rotate_to_normal(const Eigen::VectorXd& U, int direction) const {
    Eigen::VectorXd U_rot = U;

    if (direction == 0) {
        // X是法向: 不旋转
        return U_rot;
    } else if (direction == 1) {
        // Y是法向: 交换x↔y
        std::swap(U_rot(1), U_rot(2));  // ρu ↔ ρv
        std::swap(U_rot(5), U_rot(6));  // Bx ↔ By
    } else if (direction == 2) {
        // Z是法向: 交换x↔z
        std::swap(U_rot(1), U_rot(3));  // ρu ↔ ρw
        std::swap(U_rot(5), U_rot(7));  // Bx ↔ Bz
    }
    return U_rot;
}

Eigen::VectorXd HLLDSolver::rotate_from_normal(const Eigen::VectorXd& U, int direction) const {
    Eigen::VectorXd U_orig = U;

    if (direction == 0) {
        return U_orig;
    } else if (direction == 1) {
        std::swap(U_orig(1), U_orig(2));  // 交换回 ρv ↔ ρu
        std::swap(U_orig(5), U_orig(6));  // 交换回 By ↔ Bx
    } else if (direction == 2) {
        std::swap(U_orig(1), U_orig(3));  // 交换回 ρw ↔ ρu
        std::swap(U_orig(5), U_orig(7));  // 交换回 Bz ↔ Bx
    }
    return U_orig;
}
```

---

## 正确性分析

### ✓ 位置1: 动量分量交换

**fvm_3d**:
```cpp
std::swap(U_rot(1), U_rot(2));  // 交换ρu和ρv
std::swap(U_rot(1), U_rot(3));  // 交换ρu和ρw
```

**验证**: ✅ 正确
- 动量分量在索引1, 2, 3中
- Y方向旋转: (1,2,3) → (2,1,3) ✓
- Z方向旋转: (1,2,3) → (3,2,1) ✓

### ✓ 位置2: 磁场分量交换

**fvm_3d**:
```cpp
std::swap(U_rot(5), U_rot(6));  // 交换Bx和By
std::swap(U_rot(5), U_rot(7));  // 交换Bx和Bz
```

**验证**: ✅ 正确
- 磁场分量在索引5, 6, 7中
- Y方向旋转: (5,6,7) → (6,5,7) ✓
- Z方向旋转: (5,6,7) → (7,6,5) ✓

### ✓ 位置3: 密度和能量分量

**fvm_3d**:
```cpp
// 不交换: U(0) 和 U(4)
```

**验证**: ✅ 正确
- 密度(0)和能量(4)是标量，方向无关
- 不应该交换

### ✓ 位置4: GLM分量ψ

**fvm_3d**:
```cpp
// 不交换: U(8) 保持不变
```

**验证**: ✅ 正确
- ψ是散度清理标量，方向无关
- 不应该交换

---

## 与OpenMHD的等价性验证

### 映射验证（Y方向例子）

**OpenMHD Y方向**:
```fortran
! 计算的索引
vn  = 1 + mod(2-1,3) = 1 + 1 = 2  → vy (访问原始的vy)
vt1 = 1 + mod(2,3)   = 1 + 2 = 3  → vz (访问原始的vz)
vt2 = 1 + mod(3,3)   = 1 + 0 = 1  → vx (访问原始的vx)

bn  = 6 + mod(2-1,3) = 6 + 1 = 7  → by
bt1 = 6 + mod(2,3)   = 6 + 2 = 8  → bz
bt2 = 6 + mod(3,3)   = 6 + 0 = 6  → bx
```

**fvm_3d Y方向旋转**:
```cpp
// 旋转前状态（原始）
U[0]=ρ, U[1]=ρu, U[2]=ρv, U[3]=ρw, U[4]=E, U[5]=Bx, U[6]=By, U[7]=Bz

// rotate_to_normal(dir=1)
swap(U[1], U[2])  // U[1]=ρv, U[2]=ρu
swap(U[5], U[6])  // U[5]=By, U[6]=Bx

// 旋转后状态
U[0]=ρ, U[1]=ρv, U[2]=ρu, U[3]=ρw, U[4]=E, U[5]=By, U[6]=Bx, U[7]=Bz
        ↑         ↑         ↑              ↑         ↑
       vn       vt2        vt1            bn       bt2    bt1=Bz(保持U[7])
```

**等价性检查**:
- OpenMHD访问 `vn=vy` 对应 fvm_3d的 `U[1]=ρv` ✅
- OpenMHD访问 `vt1=vz` 对应 fvm_3d的 `U[3]=ρw` ✅
- OpenMHD访问 `vt2=vx` 对应 fvm_3d的 `U[2]=ρu` ✅
- OpenMHD访问 `bn=by` 对应 fvm_3d的 `U[5]=By` ✅
- OpenMHD访问 `bt1=bz` 对应 fvm_3d的 `U[7]=Bz` ✅
- OpenMHD访问 `bt2=bx` 对应 fvm_3d的 `U[6]=Bx` ✅

### 结论

**rotate_from_normal函数在逻辑上是正确的**。

---

## 为什么测试失败？

如果`rotate_from_normal`是正确的，那么数值不稳定性的原因**不在于坐标旋转**。

### 问题所在的其他位置

根据之前的测试分析，真正的问题可能是：

1. **compute_state_L/R()中的能量公式** (最可能 - 90%)
   - Poynting通量项的符号
   - 分母`(S_L - S_M)`的数值稳定性

2. **compute_central_state()中的能量计算** (可能 - 70%)
   - v·B修正项的计算
   - 中间态密度的选择

3. **初始条件设置** (可能 - 30%)
   - Orszag-Tang初期KE计算异常
   - 能量初始化错误

4. **CFL条件或时间积分** (可能 - 50%)
   - CFL = 0.4过于保守
   - RK2不适合高度磁化流体

---

## 建议

基于本分析：

### ✅ 保留rotate_from_normal
该函数逻辑正确，不需要修改。

### 🔴 关注其他区域
1. **优先检查**: `compute_state_L()` 和 `compute_state_R()` 的能量公式
2. **次优先**: `compute_central_state()` 的Poynting通量项
3. **再次**: OpenMHD源代码中如何计算中间态能量

### 验证步骤

添加详细的调试输出到Riemann求解器：

```cpp
// 在 compute_state_L() 中
std::cout << "DEBUG State L:" << std::endl;
std::cout << "  U_L(4) [E] = " << U_L(4) << std::endl;
std::cout << "  E_Lstar = " << E_Lstar << std::endl;
std::cout << "  v·B terms = " << vdotB_L << " -> " << vdotB_Lstar << std::endl;

// 在 compute_central_state() 中
std::cout << "DEBUG Central State:" << std::endl;
std::cout << "  E_Lstar = " << E_Lstar << std::endl;
std::cout << "  v·B_Lstar = " << vdotB_Lstar << std::endl;
std::cout << "  v·B_central = " << vdotB_central << std::endl;
std::cout << "  E_central = " << E_central << std::endl;
```

---

## 总结

| 项目 | 状态 | 备注 |
|------|------|------|
| **rotate_to_normal** | ✅ 正确 | 逻辑清晰，实现无误 |
| **rotate_from_normal** | ✅ 正确 | 正确逆转坐标，与OpenMHD等价 |
| **数值不稳定性根本原因** | 🔴 未在旋转中 | 需在其他区域寻找 |
| **下一步调查方向** | 🔍 能量公式 | compute_state_L/R中的E计算 |

---

*分析日期: 2025-11-18*
*对标版本: OpenMHD 20250804, fvm_3d commit b53c6a2*
