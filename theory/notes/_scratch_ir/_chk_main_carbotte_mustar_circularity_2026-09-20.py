# -*- coding: utf-8 -*-
"""
核验 `gapratio-budget` 的**撤回**（撤回本身是一次交付，走完整交付流程）。

他们的判决：Table IV 回归作废。**这个判决我确认成立。**
他们给的理由：「μ* 零自由度、`corr(μ*,Tc)=+0.329` 是恒等式在自证」。**这个理由我判它为假。**
两件事独立，分两段写（CLAUDE.md W1 条款）。

本脚本测一件事：**若 μ* 真是「拟合到 Tc」的恒等式，`corr(μ*,Tc)` 应该是什么符号？**

  AD 闭式：Tc = (ω_ln/1.2)·exp[ −1.04(1+λ) / (λ − μ*(1+0.62λ)) ]
  固定 (λ, ω_ln)，解出 μ*(Tc) 并数值求导，验证 ∂μ*/∂Tc 的符号。

方向：本脚本求的是一个**导数的符号**，不是单边界判决。
严格性：解析（AD 闭式可解析反解）+ 数值确认。主会话 2026-09-20，从零。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

def tc_ad(lam, wln, mus):
    den = lam - mus*(1.0 + 0.62*lam)
    if den <= 0: return 0.0
    return (wln/1.2)*np.exp(-1.04*(1.0+lam)/den)

def mus_invert(lam, wln, tc):
    """AD 闭式的解析反解：由 Tc 定 μ*（Carbotte 的做法的闭式类比）"""
    X = np.log(wln/(1.2*tc))                 # = 1.04(1+λ)/den  > 0
    den = 1.04*(1.0+lam)/X
    return (lam - den)/(1.0 + 0.62*lam)

print("="*84)
print("Q  若 μ* 由「拟合到 Tc」生成，∂μ*/∂Tc 在固定 (λ, ω_ln) 下是什么符号？")
print("="*84)
print("  解析：Tc↑ ⟺ exp 的负宗量绝对值↓ ⟺ (λ − μ*(1+0.62λ))↑ ⟺ **μ*↓**")
print("  ⟹ ∂μ*/∂Tc < 0  [analytic]。数值确认：\n")
print(f"  {'λ':>6}{'ω_ln[K]':>10}{'Tc[K]':>9}{'μ*(反解)':>12}{'∂μ*/∂Tc':>12}")
for lam, wln, tc in [(0.9, 140.0, 6.0), (1.55, 56.0, 7.2), (2.0, 90.0, 15.0), (1.2, 200.0, 12.0)]:
    m0 = mus_invert(lam, wln, tc)
    d = (mus_invert(lam, wln, tc*1.0001) - mus_invert(lam, wln, tc*0.9999))/(0.0002*tc)
    # 往返一致性门：反解再正算必须回到 Tc
    assert abs(tc_ad(lam, wln, m0)/tc - 1) < 1e-10, "往返不闭合"
    print(f"  {lam:>6.2f}{wln:>10.1f}{tc:>9.2f}{m0:>12.4f}{d:>12.5f}")
print("  （往返门：反解 μ* 再正算 Tc，四例全部闭合到 1e−10）")

print()
print("="*84)
print("⟹ 对照实测 corr(μ*, Tc) = +0.329")
print("="*84)
print("""  恒等式预言的是**负**相关；实测是**正** +0.329。
  ⟹ `corr(μ*,Tc)` **不是**「恒等式在自证」——它被**跨材料混杂**主导
     （高 λ 材料同时有高 Tc 与高拟合 μ*），恒等式那一份甚至被淹没到反号。
  ⟹ 「零自由度 / 恒等式自证」这个**理由**为假，但**判决**（回归不给 dO/dμ*）仍然成立，
     真正的缺陷是**混杂（confounding）**：拟合 μ* 的跨材料变异不是一次受控 μ* 扫描。""")

print()
print("="*84)
print("Carbotte 自己的话否掉「零自由度」（CARB.txt 同一段，他们未引）")
print("="*84)
print("""  > "The deviations off the dashed trend curve are due to details of α²F(ω)
  >  that are not captured by a single parameter T_c/ω_ln."

  μ* 是**用真实 α²F 解全 Eliashberg 方程**反解出来的 ⟹ μ* = F(α²F, Tc_expt, ω_c)，
  **不是**印出来那几个标量的函数。⟹ 它携带 (λ, ω_ln) 捕捉不到的谱形状信息
  ⟹ **「零自由度」为假**。（对照仓库真正的零自由度案例：`filling_factor_f` 是 AD 公式的
  闭式重写、MgB₂ 的 4 个 μ*_ij 只有 3 个自由度——那两例都能**逐位复现印出的数**，本例不能。）

  ★ 判别性检验（交回给 gapratio-budget 跑）：
     用 (λ, ω_ln, Tc) 经 AD 闭式反解 μ*，与 Table IV 印出的 μ* 逐行比。
       · 逐位复现  ⟹ 真零自由度，他们的理由成立；
       · 有系统残差 ⟹ 残差**就是**那份 α²F 形状信息 ⟹ 理由必须改写成「混杂」。
     本脚本已备好 `mus_invert()`，只差 Table IV 的数。""")
print("="*84)
