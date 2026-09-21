# -*- coding: utf-8 -*-
"""
CaC₆ 是不是单带？—— 这是我把 CaC₆ 排在 MgB₂ 之前的**前提**，本脚本检验它。

数据源：Sanna et al., PRB 85, 184514 (2012) Table I（本会话从本地 PDF 文本层直接抽，
未经 agent 转手）。三带与六带两套 λ_JJ' + DOS。

三件事：
 A) 用**自洽性**定下标约定（caption 说「第一指标跑列」——但不从名字推语义，
    用「λ_iso 必须复现印出的 0.870」与「细致平衡 N_J λ_JJ' = N_J' λ_J'J」两条独立地定）；
 B) 算多带放大 ρ(λ)/λ_avg —— 它是 𝔣_scalar / 𝔣_multi 的平方，决定 CaC₆ 能否进趋势回归；
 C) 与 MgB₂ 的 1.2703 对照。

主会话 2026-09-20，从零。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

# ---- Table I 三带块，**逐字**照印刷版式抄下（行=印刷行，列=印刷列） ----
P3 = np.array([[0.301, 0.136, 0.257],
               [0.546, 0.239, 0.479],
               [0.427, 0.198, 0.367]])
N3 = np.array([0.412, 0.104, 0.249])          # states/eV/spin

P6 = np.array([[0.163, 0.126, 0.099, 0.033, 0.201, 0.046],
               [0.179, 0.140, 0.105, 0.035, 0.221, 0.050],
               [0.331, 0.245, 0.151, 0.084, 0.384, 0.096],
               [0.271, 0.202, 0.206, 0.047, 0.400, 0.080],
               [0.252, 0.194, 0.145, 0.061, 0.309, 0.073],
               [0.206, 0.157, 0.128, 0.044, 0.259, 0.060]])
N6 = np.array([0.250, 0.176, 0.075, 0.031, 0.200, 0.056])

LAM_ISO_PRINTED = 0.870
DOS_ISO_PRINTED = 0.787

def diag_report(tag, P, N, lam_print=None, dos_print=None):
    print("=" * 96)
    print(f"{tag}")
    print("=" * 96)
    print(f"  ΣN = {N.sum():.4f}" + (f"   印出的各向同性 DOS = {dos_print}"
          f"   差 {100*(N.sum()/dos_print-1):+.2f}%" if dos_print else ""))
    for tagc, LAM in [("读法 A: λ_JJ' = 印刷[J][J']（第一指标跑**行**）", P),
                      ("读法 B: λ_JJ' = 印刷[J'][J]（第一指标跑**列**，= caption）", P.T)]:
        lam_iso = float(N @ LAM @ np.ones(len(N)) / N.sum())
        # 细致平衡残差 N_J λ_JJ' - N_J' λ_J'J
        DB = (N[:, None] * LAM) - (N[None, :] * LAM.T)
        scale = np.abs(N[:, None] * LAM).max()
        print(f"\n  {tagc}")
        print(f"    λ_iso = Σ_JJ' N_J λ_JJ' / ΣN = {lam_iso:.4f}"
              + (f"   vs 印出 {lam_print}   差 {100*(lam_iso/lam_print-1):+.2f}%" if lam_print else ""))
        print(f"    细致平衡最大残差 |N_Jλ_JJ' − N_J'λ_J'J| = {np.abs(DB).max():.5f}"
              f"   (相对 {100*np.abs(DB).max()/scale:.1f}%)")
    return

diag_report("A. 三带块：定下标约定（两条独立自洽性）", P3, N3,
            LAM_ISO_PRINTED, DOS_ISO_PRINTED)
diag_report("A'. 六带块：同样两条", P6, N6, LAM_ISO_PRINTED, DOS_ISO_PRINTED)

print()
print("=" * 96)
print("B. 多带放大：ρ(λ)/λ_avg —— 决定 CaC₆ 能否进**趋势回归**")
print("=" * 96)
print("""  依据（仓库既有 [analytic]，笔记 ad75_family_independence §8.3）：
    细致平衡 ⟹ M2_avg ≤ ρ(M2) ≤ ρ(M2″_Z)，且
    𝔣_scalar/𝔣_multi = √(ρ(M2″_Z)/M2_avg) ≥ 1，逐材料不同 ⟹ 不是同一个量。
  这里用 λ 代 M2 作**一阶指示**（等价于假设 ⟨ω²⟩ 逐通道相同）。""")
print(f"\n  {'块':<26}{'λ_avg':>10}{'ρ(λ)':>10}{'ρ/λ_avg':>11}{'√(ρ/avg)':>11}")
for tag, P, N in [("CaC₆ 三带", P3, N3), ("CaC₆ 六带", P6, N6)]:
    for conv, LAM in [("(A 行)", P), ("(B 列)", P.T)]:
        avg = float(N @ LAM @ np.ones(len(N)) / N.sum())
        rho = float(max(abs(np.linalg.eigvals(LAM))))
        print(f"  {tag+conv:<26}{avg:>10.4f}{rho:>10.4f}{rho/avg:>11.4f}{np.sqrt(rho/avg):>11.4f}")

print("""
  对照 MgB₂（仓库既有）：√(ρ(M2″_Z)/M2_avg) = 1.2703 ⟹ 因此 MgB₂ 被拦在趋势回归外。""")
print("=" * 96)

# ------------------------------------------------------------------ C
print()
print("=" * 96)
print("C. 剩下的混杂量：Sanna 的 μ* = 0.21 —— 它是**拟合**出来的，且远超 AD75 的 [0.081, 0.117]")
print("=" * 96)
print("""  Sanna 原文：「μ* is adjusted to obtain the experimental Tc」
             「experimental Tc = 11.5 K is reproduced for μ* ≃ 0.21」
  ⟹ 0.21 **不是**对 CaC₆ 库仑赝势的独立测量，它吸收了该文全部模型误差。
  ⟹ 但无论它是真物理还是模型误差，只要 CaC₆ 的有效 μ* 显著高于 AD75 人群的 ~0.10，
     趋势检验就多了**第三个**混杂源（前两个：谱形状族、测量技术）。
     好消息：这一个是**可定量**的——用 AD75 自己的公式扫 μ*。""")

def tc_ad(lam, wlog, wb2, mu):
    r = wb2 / wlog
    L1 = 2.46 * (1 + 3.8 * mu); L2 = 1.82 * (1 + 6.3 * mu) * r
    f1 = (1 + (lam / L1) ** 1.5) ** (1.0 / 3.0)
    f2 = 1 + (r - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)
    return (wlog / 1.2) * np.exp(-1.04 * (1 + lam) / (lam - mu * (1 + 0.62 * lam))) * f1 * f2

lam_c = 0.870
print(f"\n  AD75 公式在 λ={lam_c} 上，Tc(μ*) / Tc(μ*=0.10) ——「库仑混杂」的大小：")
print(f"  {'ω_log/ω̄₂':>11}" + "".join(f"{'μ*=%.2f'%m:>12}" for m in [0.10, 0.13, 0.17, 0.21]))
for rat in [0.60, 0.70, 0.80]:
    wl, wb = 25.0, 25.0 / rat
    base = tc_ad(lam_c, wl, wb, 0.10)
    print(f"  {rat:>11.2f}" + "".join(f"{tc_ad(lam_c,wl,wb,m)/base:>12.4f}"
                                      for m in [0.10, 0.13, 0.17, 0.21]))
print("""
  ⟹ μ* 从 0.10 抬到 0.21，AD75 预言的 Tc **掉到 42.2%**（即 𝔣 掉 57.8%）。
     这比「预测区间 ±41%」**还大** ⟹ **μ* 是比谱形状更大的混杂源，不能忽略。**

  **对预注册判据的后果（现在写下，不是事后调整）**：
    · 若 CaC₆ 落 (i)（S 在支撑域内）且 `Tc_exp` **低于**预测区间
      ⟹ 该偏离与「CaC₆ 的有效 μ* 高于 AD75 人群」**方向一致**，
        ⟹ **不得**直接判为「AD75 的 𝔣 关系失效」——必须先扣掉 μ* 通道。
        可操作做法：报 `𝔣_obs` 与 **AD75 在 μ*=0.21 处**的预言之比，两个数都给。
    · 若 `Tc_exp` **高于**预测区间 ⟹ 与 μ* 方向**相反** ⟹ 该信号更干净、更承重。
  **方向已审：μ* 混杂只能把 CaC₆ 往低推。**""")
print("=" * 96)

# ------------------------------------------------------------------ D
print()
print("=" * 96)
print("D. ★ 订正 A 段自己的论证：λ_iso 闭合是**简并**的，判别器是细致平衡")
print("=" * 96)
print("""  A 段说「λ_iso 闭合 + 细致平衡两条独立地定下标约定」——**「λ_iso 定约定」那半句过强。**
  真实自由度：待定的是**两个**二元选择（矩阵读法 / DOS 挂首还是次指标），
  而 λ_iso 只给**一个**方程 ⟹ 它只约束二者的**组合**，四种组合塌成两个值。
  （CLAUDE.md「把一组数字当作 N 个独立检验之前，先数它的自由度」。）""")
print(f"\n  {'矩阵读法':<30}{'DOS 挂':<14}{'λ_iso':>9}{'差':>10}")
for mtag, M in [("行=首指标", P3), ("列=首指标（caption 原话）", P3.T)]:
    for wtag, val in [("**首**指标", N3 @ M @ np.ones(3) / N3.sum()),
                      ("次指标", np.ones(3) @ M @ N3 / N3.sum())]:
        print(f"  {mtag:<30}{wtag:<14}{val:>9.4f}{100*(val/LAM_ISO_PRINTED-1):>9.2f}%")
print("""
  ⟹ 两条对角线上的组合（行/首 与 列/次）**给同一个 0.8685**，λ_iso 分不开它们。

  真正的判别器是**细致平衡**——它是矩阵 + DOS 的性质，与 λ_iso 的加权约定无关：""")
for tag, M in [("行=首指标", P3), ("列=首指标（caption）", P3.T)]:
    A = N3[:, None] * M
    r = np.abs(A - A.T).max() / np.abs(A).max()
    print(f"    {tag:<24} 最大相对残差 {100*r:6.1f}%    "
          f"N₁λ₁₂ = {N3[0]*M[0,1]:.4f}  vs  N₂λ₂₁ = {N3[1]*M[1,0]:.4f}")
print("""
  ⟹ **物理首指标 = 印刷的行**；再代回 λ_iso ⟹ **N 挂首指标**。两个选择都被唯一定死。
  ⟹ caption 的 "The first index runs over the columns" 与它自己的数据矛盾。
     **两种仪器都跑过了**（CLAUDE.md 要求）：`pdftotext` 文本层的词级坐标
     （表头 1/2/3 在 x=359/395/432，行标在 x=316，DOS 行在 y=221.6）
     与 6× 页图 `_sanna_tableI_zoom.png` **逐字一致** ⟹ 不是抽取伪影，是原文笔误。""")

print()
print("  后果（数字化真正要用的那一步）：Fig.1 panel = 入射 = 首指标，颜色 = 出射 = 次指标")
wp = np.repeat(N3, 3) / np.repeat(N3, 3).sum()
wc = np.tile(N3, 3) / np.tile(N3, 3).sum()
print(f"    按 panel 加权（**正确**）: {np.round(wp, 4)}")
print(f"    按颜色加权（**错误**）  : {np.round(wc, 4)}")
print(f"    两组权重 L1 距离 = {np.abs(wp-wc).sum():.4f} ⟹ 挂错指标 ⟹ α²F_iso 形状完全不同")
print("""
  ⚠ 三带 ΣN 比印出的各向同性 DOS 少 2.80%（六带 +0.13% ⟹ 六带才是完整划分）；
    归一化权重里大部分会约掉，但**不保证完全约掉**，它会在「α²F_iso 的 λ 是否回到 0.8685」
    那一步现形 —— 照实报，不许调。
  ⚠ 九个 λ_JJ' 闭合检验**不用 DOS 权重**（每个 panel 内共享同一 N_J）
    ⟹ 验收门与本节的权重争议**正交**，可以先过门。""")
print("=" * 96)
