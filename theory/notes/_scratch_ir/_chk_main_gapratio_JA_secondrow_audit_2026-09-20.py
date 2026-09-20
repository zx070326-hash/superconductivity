# -*- coding: utf-8 -*-
"""J_A 第二行（∂gr/∂λ, ∂gr/∂μ*）—— 唯一从未被审过的两项。

gapratio-noise 的判决 `A 胜 ⟺ σ_gap/σ_α < |det J_A/det J_B|` 里，det(J_A) 完全由
这两项决定。前序核实已查完 J 的第一行（∂Tc/∂·，步长稳）与 J_B 第二行
（∂α/∂·，**因库仑斜坡扭结而逐点导数不存在**）。这两项是最后的缺口，且仓库记录
2Δ/kTc 对 N **极敏感**（λ=0.3：N=400 给 4.05、N=8000 给 3.55）⟹ 不查不能用。

四个独立的坑，逐个测，互不替代：
§1  solve_gap 的 `converged` 标志 —— peer 的脚本**没有检查它**。不收敛的解不是解。
§2  步长收敛（与 α 同一套步长），并同时印 `w_cut/(2πT)` 到整数的距离
    （斜坡在 T_low 与在 Tc 处**位置不同**，是两个独立的扭结源）
§3  N 收敛：gr 锚点与两个导数各自的 N 阶梯
§4  T_low 因子（0.03Tc 是 T→0 的代理）：gr 与两个导数对它的依赖

方向注记：全为双侧中心差分与点算，无单边界放缩。
严格性：numerical-evidence。判词全部由实测数按显式条件算出，不预写。
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
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K   # noqa: E402

LAM, WE, MU, WCUT = 1.876, 106.9, 0.13, 500.0
NM_PEER, TLOW_PEER = 1600, 0.03
PEER = {"gr": 4.629, "dgr_dlam": +0.738, "dgr_dmu": -2.068,
        "dTc_dlam": +93.854, "dTc_dmu": -385.263}
NOTCONV = []          # 记录所有未收敛的解，最后统一报


def solver(lam, mu, nm):
    return Solver(spectrum=Spectrum.einstein(lam=lam, w0=WE), mu_star=mu,
                  w_cut=WCUT, n_matsubara=nm)


def gr(lam, mu, nm=NM_PEER, tlow=TLOW_PEER):
    """2Δ(iω₀)/(k_B Tc)，与 peer 逐行同构（含它用 Δ(iω₀) 而非实轴能隙这一点）。"""
    s = solver(lam, mu, nm)
    tc = s.tc()
    wn, d, z, conv = s.solve_gap(tlow * tc)
    if not conv:
        NOTCONV.append((lam, mu, nm, tlow))
    return 2.0 * d[0] / (KB_MEV_PER_K * tc), tc, conv, WCUT / (2 * np.pi * tlow * tc * KB_MEV_PER_K)


t0 = time.time()
print("=" * 100)
print(f"§1  锚点与收敛标志（λ={LAM}, μ*={MU}, w_cut={WCUT}, N={NM_PEER}, T_low={TLOW_PEER}·Tc）")
print("=" * 100)
g0, tc0, conv0, nc_low = gr(LAM, MU)
nc_tc = WCUT / (2 * np.pi * tc0 * KB_MEV_PER_K)
print(f"  Tc       = {tc0:.4f} K        （peer 报 188.58）")
print(f"  2Δ/kTc   = {g0:.6f}          （peer 报 {PEER['gr']}）  "
      f"差 {g0-PEER['gr']:+.4f}")
print(f"  solve_gap converged = {conv0}   ← peer 的脚本未检查此标志")
print(f"  斜坡位置：w_cut/(2πTc)   = {nc_tc:.4f}（距整数 "
      f"{min(nc_tc%1, 1-nc_tc%1):.4f}）")
print(f"            w_cut/(2πT_low) = {nc_low:.4f}（距整数 "
      f"{min(nc_low%1, 1-nc_low%1):.4f}）  ← 与 Tc 处是**两个**独立扭结源")

print()
print("=" * 100)
print("§2  步长收敛（N=1600, T_low=0.03Tc）")
print("=" * 100)
print(f"  {'Δλ/λ':>9}{'∂gr/∂λ':>13}   |{'Δμ*':>9}{'∂gr/∂μ*':>13}")
row2 = {}
for rl, hm in zip((0.15, 0.05, 0.02), (0.05, 0.02, 0.01)):
    h = rl * LAM
    dl = (gr(LAM + h, MU)[0] - gr(LAM - h, MU)[0]) / (2 * h)
    dm = (gr(LAM, MU + hm)[0] - gr(LAM, MU - hm)[0]) / (2 * hm)
    row2[rl] = (dl, dm)
    print(f"  {rl:>9.3f}{dl:>13.6f}   |{hm:>9.4f}{dm:>13.6f}")
a15, afine = row2[0.15], row2[0.02]
print(f"  peer（Δλ/λ=0.15, Δμ*=0.05）：{PEER['dgr_dlam']:+.6f} / {PEER['dgr_dmu']:+.6f}")
rep_l, rep_m = abs(a15[0] - PEER["dgr_dlam"]), abs(a15[1] - PEER["dgr_dmu"])
print(f"  我复现同步长：{a15[0]:+.6f} / {a15[1]:+.6f}   逐项差 {rep_l:.4f} / {rep_m:.4f}"
      f"  ⟹ {'复现' if rep_l < 0.01 and rep_m < 0.02 else '★ 不复现'}")
drift_l = abs(afine[0] - a15[0]) / max(abs(afine[0]), 1e-12)
drift_m = abs(afine[1] - a15[1]) / max(abs(afine[1]), 1e-12)
print(f"  ★ 粗步→细步漂移：∂gr/∂λ {drift_l:.2%}，∂gr/∂μ* {drift_m:.2%}"
      f"  ⟹ {'步长稳（与 ∂α 那两项不同）' if max(drift_l, drift_m) < 0.10 else '★ 步长不稳，同 ∂α 的病'}")

print()
print("=" * 100)
print("§3  N 收敛（步长固定 Δλ/λ=0.05, Δμ*=0.02）")
print("=" * 100)
print(f"  {'N':>7}{'Tc(K)':>11}{'2Δ/kTc':>11}{'∂gr/∂λ':>13}{'∂gr/∂μ*':>13}{'det J_A':>12}{'秒':>8}")
lad = {}
for nm in (400, 800, 1600, 3200):
    t1 = time.time()
    g, tc, _, _ = gr(LAM, MU, nm)
    h, hm = 0.05 * LAM, 0.02
    dl = (gr(LAM + h, MU, nm)[0] - gr(LAM - h, MU, nm)[0]) / (2 * h)
    dm = (gr(LAM, MU + hm, nm)[0] - gr(LAM, MU - hm, nm)[0]) / (2 * hm)
    det = PEER["dTc_dlam"] * dm - PEER["dTc_dmu"] * dl
    lad[nm] = (g, dl, dm, det)
    print(f"  {nm:>7}{tc:>11.4f}{g:>11.6f}{dl:>13.6f}{dm:>13.6f}{det:>12.4f}{time.time()-t1:>8.1f}")
rel = abs(lad[3200][3] - lad[1600][3]) / abs(lad[3200][3])
print(f"  ★ det(J_A)：N=1600 → 3200 相对变化 {rel:.3%}"
      f"  ⟹ {'已收敛' if rel < 0.01 else '★ 未收敛，peer 的 N=1600 不够'}")
print(f"     （对照仓库旧战例：λ=0.3 时 gr 从 N=400 的 4.05 掉到 N=8000 的 3.55，"
      f"此处 λ={LAM} 强耦合，实测 N=400→3200 gr 变化 "
      f"{abs(lad[3200][0]-lad[400][0])/lad[3200][0]:.3%}）")

print()
print("=" * 100)
print("§4  T_low 因子（0.03Tc 是 T→0 的代理，peer 未做敏感性）")
print("=" * 100)
print(f"  {'T_low/Tc':>10}{'2Δ/kTc':>11}{'∂gr/∂λ':>13}{'∂gr/∂μ*':>13}{'det J_A':>12}")
tl = {}
for f_ in (0.05, 0.03, 0.02):
    g, _, _, _ = gr(LAM, MU, NM_PEER, f_)
    h, hm = 0.05 * LAM, 0.02
    dl = (gr(LAM + h, MU, NM_PEER, f_)[0] - gr(LAM - h, MU, NM_PEER, f_)[0]) / (2 * h)
    dm = (gr(LAM, MU + hm, NM_PEER, f_)[0] - gr(LAM, MU - hm, NM_PEER, f_)[0]) / (2 * hm)
    det = PEER["dTc_dlam"] * dm - PEER["dTc_dmu"] * dl
    tl[f_] = det
    print(f"  {f_:>10.3f}{g:>11.6f}{dl:>13.6f}{dm:>13.6f}{det:>12.4f}")
spread = (max(tl.values()) - min(tl.values())) / abs(np.mean(list(tl.values())))
print(f"  ★ det(J_A) 在 T_low∈[0.02,0.05] 上的相对跨度 {spread:.2%}")

print()
print("=" * 100)
print("§5  合并判决：阈值 |det J_A/det J_B| 的区间")
print("=" * 100)
SG1, SG2 = (0.04 * PEER["gr"]) / 0.05, (0.075 * PEER["gr"]) / (0.15 * 0.4973)
print(f"  σ_gap/σ_α：模型1 = {SG1:.4f}，模型2 = {SG2:.4f}   （σ_Tc 已证相消 [analytic]）")
# det(J_B) 的三个候选（来自 α 光滑性扫描：局部干净支 / 整段割线 / AD 闭式）
DB = {"peer（粗步）": (+0.006, -0.154),
      "局部干净支": (+0.018000, -0.101438),
      "整段割线≈周期平均": (+0.006219, -0.177161),
      "AD 闭式（换算后）": (+0.007066, -0.206062)}
detA_lo, detA_hi = min(lad[3200][3], tl[0.02], tl[0.05]), max(lad[3200][3], tl[0.02], tl[0.05])
print(f"  det(J_A) 的合理区间（N=3200 与 T_low 扫描的包络）= [{detA_lo:.3f}, {detA_hi:.3f}]")
print(f"  {'det(J_B) 来源':>22}{'det J_B':>11}{'阈值下':>10}{'阈值上':>10}{'模型1余量':>11}{'模型2余量':>11}  判决")
for nm, (gl, gm) in DB.items():
    db = PEER["dTc_dlam"] * gm - PEER["dTc_dmu"] * gl
    t_lo, t_hi = abs(detA_lo) / abs(db), abs(detA_hi) / abs(db)
    m1, m2 = min(t_lo, t_hi) / SG1, min(t_lo, t_hi) / SG2
    v = "A 胜(两模型)" if min(m1, m2) > 1 else ("A 仅胜模型1" if m1 > 1 else "★ A 败")
    print(f"  {nm:>22}{db:>11.4f}{min(t_lo,t_hi):>10.3f}{max(t_lo,t_hi):>10.3f}"
          f"{m1:>11.3f}{m2:>11.3f}  {v}")

print()
if NOTCONV:
    print(f"  ★★ solve_gap 未收敛的点共 {len(NOTCONV)} 个：{NOTCONV[:6]}")
else:
    print("  solve_gap 在本脚本所有调用上均 converged=True")
print("=" * 100)
print(f"完成，用时 {time.time()-t0:.0f} 秒")
print("=" * 100)
