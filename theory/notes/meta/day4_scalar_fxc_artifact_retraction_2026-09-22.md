# Day 4 订正：「P 比 chi^s 更正」是自由参数伪影，非结构性差异

> 2026-09-22，Team B Phase 2 Week 2 Day 4。
> **订正**：`scdft_tddft_dyson_first_layer_2026-09-22.py` 首次运行报告的
> "[PARTIAL SUCCESS] P 比 chi^s 更正（4→2 负）"是**伪影**，已独立复核撤回。
> 真实原因：当前 Yukawa f^xc(q) 实现退化为**标量×单位阵**，与 Hubbard-U 公式数学等价，
> 只是选用了不同的标量值，不构成"对象错位假说"（出路 2）的证据。

## 一、sanity check：P 与 chi^s 公式等价性

当 f^xc(q) → 标量常数 × I（与我的 Yukawa 实现在轨道空间的实际结构一致）：

```
fxc(Q_AFM) = +1.0 * I  （人为设为与 U 相同）

P (TDDFT 公式, fxc=+1.0*I): [-2.618, -0.808, -0.512, -0.235, 0.066]
chi^s (Hubbard 公式, U=1.0): [-2.618, -0.808, -0.512, -0.235, 0.066]

max diff = 0.0  （逐位精确相等）
```

**结论**：TDDFT-Dyson 公式 `P = chi0/(1 - fxc·chi0)` 与 Hubbard-RPA 公式
`chi^s = chi0/(1 - U·chi0)` 在 f^xc 退化为标量的情形下**是同一个公式**，
只是把标量参数从 U 换成了 f^xc 的（此时也是标量的）取值。

## 二、扫描标量值：符号结构对参数极度敏感

```
fxc_scalar=-3.000: n_neg=0   <- 全正！
fxc_scalar=-2.254: n_neg=2   <- Day 4 报告用的值（Yukawa A=1.5,kappa=1.5 在 Q_AFM 处的输出）
fxc_scalar=-1.000: n_neg=2
fxc_scalar=-0.500: n_neg=2
fxc_scalar=+0.500: n_neg=3
fxc_scalar=+1.000: n_neg=4   <- Day 1 用的 U=1.0
fxc_scalar=+2.000: n_neg=4
fxc_scalar=+3.000: n_neg=4
```

**这就是问题所在**：负本征值个数从 0 到 4 全看你选什么标量值。
Day 4 报告的"P 比 chi^s 更正"，只是因为 Yukawa(A=1.5, κ=1.5) 在 Q_AFM 处
恰好给出了一个与 U=1.0 不同的负数（-2.254），**不是因为 f^xc 有非平凡的
q 依赖矩阵结构**——它仍然是标量×I，跟 U 是同一类对象。

## 三、真问题：我的 f^xc(q) 实现里没有真正的轨道结构

回看 `yukawa_fxc_matrix()`：

```python
fxc_q += fxc_R * phase * np.eye(5)   # <- 每个 R 项都乘 np.eye(5)
```

这意味着**所有轨道分量共享同一个傅里叶求和结果**，f^xc(q) 恒为
`(标量函数(q)) × I`，在轨道空间没有非对角元素、没有轨道选择性。

这与"出路 2"设想的机制（**q 依赖 + 轨道结构选择性增强某些通道**）完全不符——
我目前实现的"q 依赖"只是一个标量函数 g(q)，套了个恒等矩阵的壳，
本质上跟 Hubbard-U（也是标量）没有区别。

## 四、需要重新设计的检验

要真正检验出路 2，f^xc(q) 必须有**非平凡轨道结构**（非对角元素、
不同轨道对不同耦合强度），来源可能是：

1. **轨道相关的 Yukawa**：耦合强度依赖轨道对称性（如 xz/yz 与 xy 不同）
2. **真实 ALDA**：基于局域密度的 xc 核，天然含轨道混合（通过实空间波函数重叠）
3. **Hund's J 耦合**：非对角 f^xc 元素来自 Hund's 规则（轨道间自旋关联）

## 五、严格性标注（订正）

- Day 4 原始报告「P 比 chi^s 更正，出路 2 部分成立」：**撤回**
- Sanity check（fxc=+1.0*I 精确复现 chi^s）：**numerical-evidence**（float64，diff=0 精确到浮点精度）
- 标量值扫描：**numerical-evidence**（8 组值，单点 q=Q_AFM）
- 结论「当前实现的 f^xc(q) 退化为标量，不构成结构性检验」：**analytic**（由代码逐行检查确认）

## Status

- [x] Day 4 上午：TDDFT-Dyson vs Hubbard-RPA 初测（[PARTIAL SUCCESS] 报告）
- [x] Day 4 下午：sanity check 发现伪影，**订正撤回**
- [ ] 需要实现真正非对角 f^xc(q)（轨道结构）才能检验出路 2
- [ ] 与 team-lead 同步订正
