# -*- coding: utf-8 -*-
"""
数「经验 r」路线的自由度（CLAUDE.md：把一组数字当作 N 个独立检验之前先数它的自由度）。

被审对象：`isotope/gap_ratio_mainsession_verification_2026-09-20.md` §9.2
（我自己写的）把 `c_emp := O_Carb/O_ME` 称作「`c` 参数依赖的一个**独立估计**」，
并由此报 `R_real = 1.3777`、`s = 2.022`、单通道余量 `λ 6.74× / μ* 2.89×`、
以及分母被放大 `+34.9%`。

本脚本问三件事：
  (A) `R_real(r_emp)` 是哪些数的函数？ME 那一列是否整列相消？
  (B) 分母的 `+34.9%` 是什么量？
  (C) `R_Carb` 自己有多少自由度？（Carbotte 只吃 t = k_B Tc/ω_log ⟹ 会不会塌回纯 Tc 比？）

方向注记：本脚本不产生单边界判决，全部是恒等式与比值的重算。
主会话 2026-09-20，从零；数字来源只有 TC / OME 两列（与被审脚本同源，便于逐位对照）。
"""
import math, random

KB   = 0.08617333262   # meV/K
WLOG = 106.9           # meV（Einstein ⟹ ω_log = ω_E）

TC  = {'mu+': 170.807, 'mu-': 209.333, 'lam+': 213.481, 'lam-': 160.660}
OME = {'mu+': 4.5278,  'mu-': 4.7346,  'lam+': 4.8231,  'lam-': 4.4078}


def carbotte(tc, wlog=WLOG):
    t = KB * tc / wlog
    return 3.53 * (1.0 + 12.5 * t * t * math.log(1.0 / (2.0 * t)))


def budget(tc, ome, wlog=WLOG):
    """完全复刻被审脚本的定义，返回它印出的每一个量。"""
    oc = {k: carbotte(v, wlog) for k, v in tc.items()}
    half = lambda d, a, b: 0.5 * (d[a] - d[b])
    DlM, DmM = half(ome, 'lam+', 'lam-'), half(ome, 'mu+', 'mu-')
    DlC, DmC = half(oc, 'lam+', 'lam-'),  half(oc, 'mu+', 'mu-')
    O0M, O0C = 0.25 * sum(ome.values()), 0.25 * sum(oc.values())
    rl = DlC / O0C - DlM / O0M
    rm = DmC / O0C - DmM / O0M
    num = DlM + O0M * rl
    den = DmM + O0M * rm
    return dict(oc=oc, DlM=DlM, DmM=DmM, DlC=DlC, DmC=DmC, O0M=O0M, O0C=O0C,
                rl=rl, rm=rm, num=num, den=den,
                R_real=abs(num / den), R_ME=abs(DlM / DmM), R_C=abs(DlC / DmC))


b = budget(TC, OME)

print("=" * 94)
print("0. 先复刻被审脚本印出的数（逐位对照，确认我审的是同一个对象）")
print("=" * 94)
print(f"   R_ME   = {b['R_ME']:.4f}   (笔记 §4: 2.0082)")
print(f"   R_Carb = {b['R_C']:.4f}   (笔记 §9.2: 1.3777)")
print(f"   r_lam_emp = {100*b['rl']:+.2f}%   r_mu_emp = {100*b['rm']:+.2f}%   (笔记: -0.33% / -0.78%)")
print(f"   R_real(r_emp) = {b['R_real']:.4f}")

print()
print("=" * 94)
print("A. `R_real(r_emp)` 的自由度：ME 那一列是否整列相消？")
print("=" * 94)
print("   代数：r_p := D_p^C/O0^C − D_p^ME/O0^ME")
print("      ⟹ D_p^ME + O0^ME·r_p = D_p^ME + D_p^C·O0^ME/O0^C − D_p^ME = (O0^ME/O0^C)·D_p^C")
print("      ⟹ R_real = |D_l^C| / |D_m^C| = R_Carb   —— 与 D^ME、O0^ME 全部无关。")
print(f"   逐位：num = {b['num']:+.8f}   (O0M/O0C)·DlC = {b['O0M']/b['O0C']*b['DlC']:+.8f}")
print(f"         den = {b['den']:+.8f}   (O0M/O0C)·DmC = {b['O0M']/b['O0C']*b['DmC']:+.8f}")
print(f"   R_real − R_Carb = {b['R_real']-b['R_C']:+.3e}")

print()
print("   随机投毒：任意改 O_ME 那一列（保持 Tc 不动），R_real 必须一动不动。")
random.seed(20260920)
worst = 0.0
for _ in range(2000):
    ome2 = {k: v * (0.3 + 2.0 * random.random()) for k, v in OME.items()}
    b2 = budget(TC, ome2)
    worst = max(worst, abs(b2['R_real'] - b['R_C']))
    if _ < 3:
        print(f"     O_ME 乱改一次 ⟹ R_ME={b2['R_ME']:7.3f}  但 R_real={b2['R_real']:.6f}")
print(f"   2000 次随机改 O_ME：max|R_real − R_Carb| = {worst:.2e}")
print("   ⟹ 【判决】`R_real(r_emp)` 携带的 Eliashberg 信息 = 0。它**不是**独立估计，")
print("      而是「把 ME 的灵敏度整列换成 Carbotte 的灵敏度」这一操作的**定义**。")

print()
print("=" * 94)
print("B. 分母 `+34.9%` 是什么量？")
print("=" * 94)
amp = abs(b['den']) / abs(b['DmM'])
print(f"   |den|/|D_mu^ME| = {amp:.4f}  ⟹ {100*(amp-1):+.1f}%   (笔记读作「实轴修正放大了分母」)")
print(f"   而 (O0M/O0C)·|D_mu^C| / |D_mu^ME| = "
      f"{b['O0M']/b['O0C']*abs(b['DmC'])/abs(b['DmM']):.4f}")
print("   ⟹ 它**不是**实轴修正的大小，而是「Carbotte 的 μ* 灵敏度 / ME 的 μ* 灵敏度」。")
print("      量的是两套公式的分歧，不是 c 的参数依赖。")

print()
print("=" * 94)
print("C. R_Carb 自己有多少自由度？（Carbotte 只吃 t ∝ Tc）")
print("=" * 94)
dTc_l = 0.5 * (TC['lam+'] - TC['lam-'])
dTc_m = 0.5 * (TC['mu+'] - TC['mu-'])
print(f"   纯 Tc 比：|ΔTc_λ/ΔTc_μ| = {abs(dTc_l/dTc_m):.4f}")
print(f"   R_Carb                 = {b['R_C']:.4f}")
print(f"   相对差（= Carbotte 在 t 上的曲率贡献）= {100*(b['R_C']/abs(dTc_l/dTc_m)-1):+.2f}%")
print("   ⟹ R_Carb 到 0.5% 就是**纯 Tc 灵敏度之比**。Carbotte 把 2Δ/kTc 压成 t 的单变量函数,")
print("      而 t ∝ Tc ⟹ 这个『第二观测量』在 Carbotte 下**退化回 Tc 本身**，零额外判别力。")

print()
print("   ω_log 换值（μ* 不改变 ω_log，所以这只测曲率项的稳健性）：")
for wl in (80.0, 106.9, 140.0, 200.0):
    bb = budget(TC, OME, wl)
    print(f"     ω_log={wl:6.1f} meV  t={KB*TC['mu-']/wl:.4f}..{KB*TC['lam+']/wl:.4f}"
          f"   R_Carb={bb['R_C']:.4f}   vs 纯Tc比 {abs(dTc_l/dTc_m):.4f}")

print()
print("=" * 94)
print("D. 什么活着：R_ME 相对纯 Tc 的真实增益")
print("=" * 94)
print(f"   R_ME / |ΔTc_λ/ΔTc_μ| = {b['R_ME']/abs(dTc_l/dTc_m):.4f}")
print("   ⟹ 能隙比相对『只看 Tc』的净增益是这个倍数，且它全部来自 Carbotte 压掉的那部分")
print("      （本征向量形状信息）。这是第二观测量**有意义**的正面陈述，与 §9.1 的解析预算无关。")
print("=" * 94)
