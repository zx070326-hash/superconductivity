# Day 5 订正：轨道结构 f^xc(q) 的"改善"仍是标量伪影，不是出路 2 的证据

> 2026-09-22，Team B Phase 2 Week 2 Day 5。
> **订正对象**：`scdft_orbital_fxc_route2_2026-09-22.py` 首次运行打印的
> "[ROUTE 2 SUPPORTED]"（4→1 负本征值，稳健性扫描 100% 改善）——**已独立复核撤回**。
> 真实原因：即便 f^xc(q) 有真正的非对角矩阵元素，负本征值个数**完全由对角元的
> 幅值/符号决定**，非对角结构对结果贡献为零。这与 Day 4 的伪影同形，只是伪装层级更深。

## 一、Day 5 脚本表面结果（未经控制，误导性）

`orbital_fxc_matrix(q)` 用轨道对相关的 (A_ij, kappa_ij) 矩阵构造了真正非对角、
Hermitian 的 f^xc(Q_AFM)（非对角/对角 Frobenius 范数比 0.57，非零）。

初次判决：
```
chi^s (Hubbard U=1.0):        4 个负本征值
P (轨道结构 fxc, TDDFT):      1 个负本征值   (k_grid=30)
稳健性扫描（12 组随机扰动矩阵）：100% 给出 n_neg(P) < 4
```
表面看像是"真正的对象错位 + 轨道结构"驱动了改善。

## 二、控制检验：把"是否有非对角结构"与"标量幅值/符号是否不同"分离

在同一个 chi0（k_grid=20）上，构造三个对照对象：

```python
fxc_full        = orbital_fxc_matrix(q_afm)              # 真正非对角矩阵
fxc_diag_only   = diag(diag(fxc_full))                   # 抹掉非对角，只留对角
fxc_scalar_equiv = mean(diag(fxc_full)) * I              # 单一标量 * I（无轨道差异）
```

结果：
```
fxc(Q_AFM) 对角元素: [-2.176, -2.176, -1.803, -1.600, -1.200]
mean(diag) = -1.791

P(diag-only fxc)         : 2 个负本征值
P(scalar=mean(diag)*I)   : 2 个负本征值
P(FULL orbital matrix)   : 2 个负本征值   <- 与前两者逐位一致

chi^s (Hubbard, U=+1.0)  : 4 个负本征值
```

**三个不同复杂度的对象（纯标量、对角矩阵、真正非对角矩阵）在同一 chi0 下给出
完全相同的负本征值计数（2）**。非对角结构对判决**零贡献**。

## 三、结论：Day 5 的"改善"是什么

改善（4→2 或 4→1，取决于 k_grid）完全来自：
- f^xc 对角元的**符号**（负，-1.2~-2.2）与 U=+1.0（正）**相反**
- f^xc 对角元的**幅值**（~1.2-2.2）与 U=1.0 **不同**

这正是 Day 4 已经暴露、被明确警告要避免的模式：**换了一个不同的标量（此处是
"对角平均值"），负本征值计数当然会变**，这不需要任何轨道结构、任何 TDDFT vs
Hubbard 的对象区分——单纯换个数字就能做到（Day 4 附录里 fxc_scalar 从 -3.0 到
+3.0 扫描时，n_neg 就已经从 0 到 4 全覆盖过）。

稳健性扫描"100% 改善"同样是伪影：12 组随机扰动矩阵的对角元均值都落在与
Day 4 类似的"负、大幅值"区间（因为 A_MATRIX 对角基准值全部为正、Yukawa 公式
本身自带负号），所以扫描的其实是"对角平均值恰好比 U=1.0 更负"这一件事，
被扰动的非对角结构完全是旁观者。

## 四、真正检验出路 2 需要什么

要让"非对角轨道结构"本身成为判决的自变量，必须**固定对角元（幅值+符号）
与 chi^s 所用的 U 完全一致**（比如都设为 +1.0），只改变非对角元的有无/大小，
然后看负本征值计数是否因非对角项而改变。当前两次尝试（Day 4 标量、Day 5
轨道矩阵）都没有做到这一点——两次都是"对角/标量整体换了个数"，不是
"仅改变非对角、控制对角不变"。

下一步（Day 6 计划）：
```python
fxc_test = U_hubbard * np.eye(5) + epsilon * OffDiagOnly(orbital_structure)
```
扫描 epsilon（非对角强度）从 0 开始增大，对角固定为 chi^s 基线的 U=1.0，
观察负本征值计数是否随 epsilon 单调变化——这才是"非对角结构是否重要"的
干净检验。

## 五、严格性标注

- Day 5 原始报告「ROUTE 2 SUPPORTED，稳健性 100%」：**撤回**
- 控制检验（diag-only vs scalar-equiv vs full 逐位相等）：**numerical-evidence**
  （float64，单一 q 点 Q_AFM，单一 chi0 网格 k_grid=20，三个确定性对象精确对比，非采样）
- 结论「本次轨道结构对判决零贡献」：**numerical-evidence**（非 analytic——
  未证明这对任意轨道矩阵成立，只是对这一个具体 A_MATRIX/KAPPA_MATRIX 参数化成立）
- 尚未回答：非对角结构在"对角固定、只变非对角"的干净实验下是否有效——**未测**

## Status

- [x] Day 5 上午：轨道结构 f^xc(q) 构造 + 初测（[ROUTE 2 SUPPORTED] 误报）
- [x] Day 5 下午：diag-only / scalar-equiv 控制检验，发现非对角贡献为零，**订正撤回**
- [ ] Day 6 计划：对角固定为 U=1.0、只扫描非对角强度 epsilon 的干净检验
- [ ] 与 team-lead 同步订正

## 相关文件

- 脚本：`../_scratch_ir/scdft_orbital_fxc_route2_2026-09-22.py`
- 前一次同形伪影：`day4_scalar_fxc_artifact_retraction_2026-09-22.md`
- Sylvester 惯性定律阻碍（第二层夹心）：`sylvester_inertia_blocks_naive_sandwich_2026-09-22.md`
