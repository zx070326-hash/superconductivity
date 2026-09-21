# -*- coding: utf-8 -*-
"""gapratio-noise 的 J_B 第二行：α 的**同位素约定**审计（从零，自己的实现）。

线索（由 `_chk_main_gapratio_jacobian_dof_audit` 的 AD 交叉带出）：
  · 它的 α(λ=1.876, μ*=0.13) = **0.4973**，距 μ*=0 的 0.5 只差 0.0027；
  · 但仓库自己的 H3S（λ≈2, μ*=0.13，几乎是同一个锚点）给 **α ≈ 0.36**；
  · AD 闭式在同一点给 0.479（截断重标定后 0.486）。
  ⟹ 50× 的差，不是步长（AD 上同步长差分 vs 解析只差 3.6–3.9%），只能是**约定**。

它的实现把 `w_cut = 500 meV` **绝对固定**，同位素替换只动 ω_E。
仓库 `isotope_coefficient` 有且只有这一个开关：
  mu_star_fixed=True  → μ* 数值原样保留（= 它的做法）
  mu_star_fixed=False → ω_c 随谱移动，μ*(ω_c) 按 Tolmachev 重新推导（= 物理约定）

本脚本在**同一个锚点**上把三种约定各算一遍，看 α 与 ∂α/∂μ* 各落在哪。
若 α 随约定大幅移动，则 J_B 第二行、进而 A-vs-B 的判决**依赖一个未声明的约定**。

⚠ 本文件 §3 的解读**已被** `_chk_main_alpha_jacobian_step_convergence_2026-09-20.py` **订正并取代**：
  §1/§2 的**测量**全部成立（且逐位复现 peer），但我给 `tolmachev` 贴的「物理约定」标签是**反的**。
  §1 实测 `scaled_wc` 下 α ≡ 0.500000（对每个 μ*）⟹ α 偏离 ½ **全部**来自 w_cut/ω_E 这个比值，
  而那个比值**就是** Tolmachev 对数本身 ⟹ 绝对固定 w_cut 的 `fixed_abs` 才是正确约定，不是赝物。
  真正的缺陷在别处：`fixed_abs` 的 ∂α/∂μ* **对差分步长不稳**（Δ=0.05 给 −0.154，Δ=0.02 给 −0.207）。

方向注记：全为双侧中心差分点算，无单边界放缩。

严格性：numerical-evidence。
主会话 2026-09-20，从零。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, time
from pathlib import Path
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "solvers"))
from eliashberg_iso import Spectrum, Solver, mu_star_at_cutoff   # noqa: E402

LAM, WE, MU, WCUT, NM = 1.876, 106.9, 0.13, 500.0, 400
NM_REF = 1600          # peer 的 N；§4 用它抽查 α 的 N 无关性（单点 29 s/解）


def tc(lam, we, mu, wcut, nm=NM):
    return Solver(spectrum=Spectrum.einstein(lam=lam, w0=we), mu_star=mu,
                  w_cut=wcut, n_matsubara=nm).tc()


def alpha(lam, mu, mode, we=WE, eps=0.01, nm=NM):
    """α = −dlnTc/dlnM，ω_E ∝ M^{-1/2} ⟹ α = ½·dlnTc/dlnω_E（中心差分）。

    mode='fixed_abs'   ω_c 绝对固定 500 meV、μ* 原样          ← peer 的做法
    mode='scaled_wc'   ω_c ∝ ω_E（比例固定）、μ* 原样
    mode='tolmachev'   ω_c ∝ ω_E，μ*(ω_c) 按 Tolmachev 从 500 meV 处重推 ← 物理约定
    """
    out = []
    for s in (1.0 + eps, 1.0 - eps):
        w = we * s
        if mode == "fixed_abs":
            wc, m = wcut_of(we), mu
        elif mode == "scaled_wc":
            wc, m = wcut_of(we) * s, mu
        elif mode == "tolmachev":
            wc = wcut_of(we) * s
            m = mu_star_at_cutoff(mu, wcut_of(we), wc)
        else:
            raise ValueError(mode)
        out.append(np.log(tc(lam, w, m, wc, nm)))
    # dlnTc/dlnω_E 的中心差分：分母必须是 ln((1+ε)/(1−ε))，
    # 写成 2ln(1+ε) 会在 O(ε) 上错 0.5%（足以把 μ*=0 的门从 0.5000 推到 0.5025）。
    return 0.5 * (out[0] - out[1]) / np.log((1.0 + eps) / (1.0 - eps))


def wcut_of(we):
    return WCUT


t0 = time.time()
print("=" * 96)
print("§1  锚点 λ=1.876, ω_E=106.9 meV, w_cut=500 meV, N=1600 —— 三种同位素约定下的 α")
print("=" * 96)
tc0 = tc(LAM, WE, MU, WCUT)
print(f"  Tc(锚点) = {tc0:.4f} K     （peer 报 188.58 K）")
print()
print(f"  {'μ*':>6}" + "".join(f"{m:>16}" for m in ("fixed_abs(peer)", "scaled_wc", "tolmachev")))
A = {}
for m in (0.0, 0.10, 0.13, 0.16):
    row = []
    for mode in ("fixed_abs", "scaled_wc", "tolmachev"):
        row.append(alpha(LAM, m, mode))
    A[m] = row
    print(f"  {m:>6.2f}" + "".join(f"{v:>16.6f}" for v in row))
print()
print("  门：μ*=0 三列都必须 = 0.5（Tc ∝ ω_E 严格标度）")
ok = all(abs(v - 0.5) < 2e-3 for v in A[0.0])
print(f"     {'PASS' if ok else '★ FAIL'}  实测 {['%.6f' % v for v in A[0.0]]}")

print()
print("=" * 96)
print("§2  ∂α/∂μ* —— det(J_B) 的主导项")
print("=" * 96)
DMU = 0.05
print(f"  {'约定':>16}{'α(0.13)':>12}{'∂α/∂μ*(Δ=0.05)':>18}{'∂α/∂μ*(Δ=0.02)':>18}")
G = {}
for j, mode in enumerate(("fixed_abs", "scaled_wc", "tolmachev")):
    d = []
    for h in (0.05, 0.02):
        d.append((alpha(LAM, MU + h, mode) - alpha(LAM, MU - h, mode)) / (2 * h))
    G[mode] = d[0]
    print(f"  {mode:>16}{A[0.13][j]:>12.6f}{d[0]:>18.6f}{d[1]:>18.6f}")
print(f"  peer 报：α = 0.497300，∂α/∂μ* = −0.154000")
print(f"  AD 闭式（解析，同一点）：α = 0.478504，∂α/∂μ* = −0.389014")

print()
print("=" * 96)
print("§3  判决对约定的敏感性：把 J_B 第二行换成各约定的值，重算阈值")
print("=" * 96)
dA = 90.234022                      # det(J_A)，与 peer 逐位一致（已在 DOF 审计里复核）
r1, r2 = (0.04 * 4.629) / 0.05, (0.075 * 4.629) / (0.15 * 0.4973)
GLAM = 0.006                        # ∂α/∂λ：暂用 peer 值（§4 单独审）
print(f"  σ_gap/σ_α：模型1 = {r1:.4f}，模型2 = {r2:.4f}    （σ_Tc 已证明相消）")
print(f"  {'∂α/∂μ* 来源':>22}{'值':>12}{'det(J_B)':>12}{'阈值':>10}{'模型1':>9}{'模型2':>9}  判决")
cands = [("peer (fixed_abs 报值)", -0.154)] + \
        [(f"我算 {m}", G[m]) for m in ("fixed_abs", "scaled_wc", "tolmachev")] + \
        [("AD 闭式解析", -0.389014)]
for nm, h in cands:
    d = 93.854 * h + 385.263 * GLAM
    t = abs(dA) / abs(d)
    v = "A 胜(两模型)" if t > max(r1, r2) else ("A 仅胜模型1" if t > min(r1, r2) else "★ A 败")
    print(f"  {nm:>22}{h:>12.6f}{d:>12.4f}{t:>10.4f}{t/r1:>9.3f}{t/r2:>9.3f}  {v}")

print()
print("=" * 96)
print(f"§4  N 无关性抽查：peer 用 N={NM_REF}，上表用 N={NM}")
print("=" * 96)
print(f"  Tc(N={NM})    = {tc0:.4f} K")
tcr = tc(LAM, WE, MU, WCUT, NM_REF)
print(f"  Tc(N={NM_REF})   = {tcr:.4f} K     （peer 报 188.58 K ⟹ 同一台仪器）")
a4 = alpha(LAM, MU, "fixed_abs", nm=NM)
a16 = alpha(LAM, MU, "fixed_abs", nm=NM_REF)
print(f"  α_fixed_abs(N={NM}) = {a4:.6f}   α_fixed_abs(N={NM_REF}) = {a16:.6f}"
      f"   Δ = {a16-a4:+.2e}")
print(f"  ⟹ {'α 对 N 已收敛，N=400 的表可直接用' if abs(a16-a4) < 5e-3 else '★ α 未对 N 收敛，上表不可用'}")
print(f"  ⟹ 与 peer 的 0.497300 差 {a16-0.4973:+.6f}"
      f"（{'复现' if abs(a16-0.4973) < 5e-3 else '★ 不复现'}）")

print()
print("=" * 96)
print(f"完成，用时 {time.time()-t0:.0f} 秒")
print("=" * 96)
