# Day 6：干净检验完成 — 出路 2（对象错位假说）在第一层 Dyson 重求和下不成立

> 2026-09-22，Team B Phase 2 Week 2 Day 6。
> 承接 Day 4/Day 5 两次撤回：两次尝试都混淆了"换标量幅值/符号"与"加非对角结构"，
> 本文做了隔离变量的干净检验，把这个混淆消除。

## 检验设计

对角**锁定**为与 chi^s 基线完全相同的 U=1.0（不再是 Day 5 里由 Yukawa 公式
自带的负值/不同幅值），只让非对角强度 epsilon 作为唯一自变量：

```python
off_diag_shape = (Day 5 的轨道 f^xc 矩阵的非对角部分) / max(|非对角元|)   # 归一化到量级 1
fxc_test(epsilon) = U_hubbard * I + epsilon * off_diag_shape
P(epsilon) = chi0 @ inv(I - fxc_test(epsilon) @ chi0)
```

epsilon=0 时 fxc_test = U_hubbard*I，精确回到 chi^s 基线（sanity check：
epsilon=0 给出 4 个负本征值，与 chi^s 完全一致，确认代码正确）。

## 结果（单点 q=Q_AFM，chi0 固定于 k_grid=20）

```
epsilon    n_neg(P)   备注
 0.00         4       = chi^s 基线（sanity check 通过）
 0.10-2.00    4        全程不变（横跨 20 倍强度范围）
 3.00         3       非对角强度已超过对角本身（U=1.0 vs |off-diag|~3），
                        物理上不合理的参数区
-0.50~-2.00   4        负号同样不变
```

**在物理合理的非对角强度范围内（0 到与对角同量级），负本征值计数恒为 4，
与没有非对角结构时完全相同。** 只有在非对角强度远超对角本身（epsilon=3，
即非对角元比 Hubbard onsite 排斥还大 3 倍）时才勉强变化到 3——这已经
超出"顶角修正是对角相互作用的微扰"这一物理图像的合理范围，且从 4→3 只是
穿过了一个接近零的本征值（0.0419 附近），不是结构性改善。

## 结论

**出路 2（Team B 的 Hubbard-U RPA chi^s ≠ Essenberger 的 TDDFT 第一层 P，
仅因矩阵结构不同就能改善符号）在此参数化下不成立。** 单纯的第一层 Dyson
重求和加非对角轨道结构，在物理合理的耦合强度范围内，不能把 4 个负本征值
压低。

这与 Day 4/Day 5 的两次"伪成功"报告加在一起，构成一致的画面：
每次看到"改善"，独立控制检验都发现改善来自换了标量幅值/符号，
不是来自真正的对象结构差异。

## 尚存的但未测的可能性（不构成对出路 2 的证伪）

1. 本检验只测了**第一层** Dyson（chi0 → P）。Essenberger 的完整构造还有
   **第二层夹心**（P → Λ^SF via f_σ^F，Eq.31-32）。Sylvester 惯性定律笔记
   已证明：若第二层夹心矩阵可逆，符号结构精确不变——所以即使第一层给出
   全正的 P，也不需要"改善"，只要第一层本身已经翻转符号即可（第二层保序）。
   本检验说明**第一层本身在合理参数下翻转不了**，这是对第一层单独的否定，
   不是对"两层加起来"的否定——但由于第二层是符号保持算符，若第一层没做到，
   两层合起来也做不到。
2. 只测了一个特定的"非对角形状"（Day 5 的 Hund's 型矩阵，归一化到单位模）。
   没有排除某个截然不同形状的非对角结构可能在合理强度下起作用——但趋势
   （epsilon 从 0.1 到 2.0 单调、平滑、几乎不变）不支持"换个形状会突然翻转"
   的猜测。
3. 只测了 q=Q_AFM 一点。

## 严格性标注

- epsilon=0 sanity check（回到 chi^s 基线，逐位一致）：**numerical-evidence**（精确浮点核验）
- epsilon 扫描 13 点，n_neg(P) 在合理范围内恒为 4：**numerical-evidence**
  （float64，单一 q 点，单一非对角形状，13 个确定性求值，非随机采样）
- 结论「出路 2（第一层 Dyson + 轨道结构）在此参数化/此 q 点下不成立」：**numerical-evidence**
- 未覆盖：第二层夹心的完整耦合、其他非对角形状、其他 q 点——**未测**，
  已在正文列出

## 相关文件

- 前两次撤回：`day4_scalar_fxc_artifact_retraction_2026-09-22.md`、
  `day5_orbital_fxc_still_scalar_artifact_2026-09-22.md`
- Sylvester 惯性定律（第二层夹心保序）：`sylvester_inertia_blocks_naive_sandwich_2026-09-22.md`
- 本次检验直接复用 `../_scratch_ir/scdft_orbital_fxc_route2_2026-09-22.py` 中的
  `orbital_fxc_matrix()` / `tddft_dyson_proper_response()`，扫描代码内联于本笔记
  （未另存脚本；若后续需要复算，扫描逻辑见上方"检验设计"代码块，可直接粘贴执行）

## Status

- [x] 隔离变量检验：对角锁定 U=1.0，仅变 epsilon
- [x] sanity check（epsilon=0 精确回到基线）
- [x] 判决：出路 2 在此第一层构造/此参数化下不支持
- [ ] 与 team-lead 同步，决定是否转向出路 1（f_σ^F 奇异）或出路 3（频率积分语境）
