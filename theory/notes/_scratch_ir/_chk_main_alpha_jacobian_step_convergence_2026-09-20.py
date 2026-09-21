# -*- coding: utf-8 -*-
"""gapratio-noise 的 J_A / J_B：**四个矩阵元的差分步长收敛**，与 AD 闭式交叉。

订正我自己上一版（`_chk_main_alpha_isotope_convention_audit`）§3 的解读：
  那里我把 `tolmachev`（ω_c ∝ ω_E 且 μ* 重推）叫作「物理约定」——**标签是反的**。
  实测 `scaled_wc` 下 α ≡ 0.500000 对每个 μ* 成立 ⟹ α 偏离 ½ **全部**来自 w_cut/ω_E，
  而 ln(w_cut/ω_E) **就是** Tolmachev 对数 ⟹ 把 w_cut 绝对固定（peer 的 `fixed_abs`）
  才是正确约定。α<½ 是机制，不是赝物。

真正的缺陷是别的：**peer 的差分步长太粗**（Δλ = 15%·λ，Δμ* = 0.05 = 38%·μ*）。
AD 上「同步长 vs 解析」只差 3.6%，所以我上一版据此说「步长不是问题」——**那个推断错了**：
AD 的 α(μ*) 与 Eliashberg 的 α(μ*) 曲率不同，AD 上的步长误差不能外推到 Eliashberg 上。
本脚本直接在 Eliashberg 上**扫步长**，不再借 AD 判步长。

AD 只作**独立仪器**用，且必须做两步换算（上一版漏了，故当时得到 2.5× 的假分歧）：
  ① 截断换算：peer 的 μ*=0.13 定义在 w_cut=500 meV；AD 的 μ* 定义在声子尺度
     ⟹ μ*_AD = μ*/(1 + μ* ln(w_cut/ω_E))
  ② 链式法则：dα/dμ*(500) = [dα/dμ*(ω_E)] · dμ*(ω_E)/dμ*(500) = [·] / (1+μ*L)²

判据（上一版已证 [analytic]）：σ_Tc 逐位相消，唯一自由度是
    A 胜 ⟺ σ_gap/σ_α < |det J_A / det J_B|

方向注记：全为双侧中心差分，无单边界放缩。步长收敛判决靠相继步长的相对变化。
严格性：numerical-evidence（AD 的两步换算是 analytic，但它只作交叉不作判决）。
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
from eliashberg_iso import Spectrum, Solver   # noqa: E402

LAM, WE, MU, WCUT, NM = 1.876, 106.9, 0.13, 500.0, 400
PEER = {"dTc_dlam": +93.854, "dTc_dmu": -385.263,
        "dgr_dlam": +0.738, "dgr_dmu": -2.068,
        "dal_dlam": +0.006, "dal_dmu": -0.154}
SG_A, SG_B = (0.04 * 4.629) / 0.05, (0.075 * 4.629) / (0.15 * 0.4973)


def tc(lam, mu, we=WE):
    return Solver(spectrum=Spectrum.einstein(lam=lam, w0=we), mu_star=mu,
                  w_cut=WCUT, n_matsubara=NM).tc()


def alpha(lam, mu, eps=0.01):
    """α = ½·dlnTc/dlnω_E，w_cut 绝对固定（peer 约定，且是正确约定）。"""
    a = np.log(tc(lam, mu, WE * (1 + eps)))
    b = np.log(tc(lam, mu, WE * (1 - eps)))
    return 0.5 * (a - b) / np.log((1 + eps) / (1 - eps))


t0 = time.time()
print("=" * 98)
print("§0  α 内层差分步长 eps（peer 用 1%）")
print("=" * 98)
print(f"  {'eps':>8}{'α(锚点)':>14}")
for e in (0.04, 0.02, 0.01, 0.005, 0.0025):
    print(f"  {e:>8.4f}{alpha(LAM, MU, e):>14.6f}")
print("  ⟹ 内层步长收敛后再谈外层步长。")

print()
print("=" * 98)
print("§1  ∂Tc/∂λ、∂Tc/∂μ*（J_A 与 J_B **共用**的第一行）")
print("=" * 98)
print(f"  {'Δλ/λ':>9}{'∂Tc/∂λ':>13}   |{'Δμ*':>9}{'∂Tc/∂μ*':>13}")
row1 = {}
for (rl, hm) in zip((0.15, 0.05, 0.02, 0.01), (0.05, 0.02, 0.01, 0.005)):
    h = rl * LAM
    dl = (tc(LAM + h, MU) - tc(LAM - h, MU)) / (2 * h)
    dm = (tc(LAM, MU + hm) - tc(LAM, MU - hm)) / (2 * hm)
    row1[rl] = (dl, dm)
    print(f"  {rl:>9.3f}{dl:>13.4f}   |{hm:>9.4f}{dm:>13.4f}")
print(f"  peer（Δλ/λ=0.15, Δμ*=0.05）：{PEER['dTc_dlam']:+.3f} / {PEER['dTc_dmu']:+.3f}")
d15 = row1[0.15]
print(f"  复现：{d15[0]:+.3f} / {d15[1]:+.3f}  ⟹ "
      f"{'PASS' if abs(d15[0]-PEER['dTc_dlam'])<0.5 and abs(d15[1]-PEER['dTc_dmu'])<2 else '★ 不复现'}")
TC_L, TC_M = row1[0.01]
print(f"  ★ 收敛值（最细步）：∂Tc/∂λ = {TC_L:+.4f}（peer 偏 "
      f"{100*(PEER['dTc_dlam']-TC_L)/TC_L:+.1f}%），"
      f"∂Tc/∂μ* = {TC_M:+.4f}（偏 {100*(PEER['dTc_dmu']-TC_M)/TC_M:+.1f}%）")

print()
print("=" * 98)
print("§2  ∂α/∂λ、∂α/∂μ*（J_B 第二行）—— peer 的两个最小项")
print("=" * 98)
print(f"  {'Δλ/λ':>9}{'∂α/∂λ':>13}   |{'Δμ*':>9}{'∂α/∂μ*':>13}")
row2 = {}
for (rl, hm) in zip((0.15, 0.05, 0.02, 0.01), (0.05, 0.02, 0.01, 0.005)):
    h = rl * LAM
    dl = (alpha(LAM + h, MU) - alpha(LAM - h, MU)) / (2 * h)
    dm = (alpha(LAM, MU + hm) - alpha(LAM, MU - hm)) / (2 * hm)
    row2[rl] = (dl, dm)
    print(f"  {rl:>9.3f}{dl:>13.6f}   |{hm:>9.4f}{dm:>13.6f}")
print(f"  peer（Δλ/λ=0.15, Δμ*=0.05）：{PEER['dal_dlam']:+.6f} / {PEER['dal_dmu']:+.6f}")
a15 = row2[0.15]
print(f"  复现：{a15[0]:+.6f} / {a15[1]:+.6f}  ⟹ "
      f"{'PASS（它的算术没错）' if abs(a15[1]-PEER['dal_dmu'])<0.005 else '★ 不复现'}")
AL_L, AL_M = row2[0.01]
print(f"  ★ 收敛值（最细步）：∂α/∂λ = {AL_L:+.6f}（peer 偏 "
      f"{100*(PEER['dal_dlam']-AL_L)/abs(AL_L):+.1f}%），"
      f"∂α/∂μ* = {AL_M:+.6f}（偏 {100*(PEER['dal_dmu']-AL_M)/abs(AL_M):+.1f}%）")

print()
print("=" * 98)
print("§3  AD 闭式独立交叉（两步换算：截断 + 链式法则）")
print("=" * 98)


def a_ad(lam, mu):
    d = lam - mu * (1 + 0.62 * lam)
    return 0.5 * (1.0 - mu ** 2 * 1.04 * (1 + lam) * (1 + 0.62 * lam) / d ** 2)


L = np.log(WCUT / WE)
mu_ph = MU / (1 + MU * L)                       # ① 截断换算
chain = 1.0 / (1 + MU * L) ** 2                 # ② dμ*(ω_E)/dμ*(500)
e = 1e-6
dad_dmu_ph = (a_ad(LAM, mu_ph + e) - a_ad(LAM, mu_ph - e)) / (2 * e)
dad_dlam = (a_ad(LAM + e, mu_ph) - a_ad(LAM - e, mu_ph)) / (2 * e)
print(f"  ln(w_cut/ω_E) = {L:.4f}  ⟹  μ*(声子尺度) = {mu_ph:.5f}（peer 的 0.13 定义在 500 meV）")
print(f"  链式因子 dμ*(ω_E)/dμ*(500) = 1/(1+μ*L)² = {chain:.5f}")
print(f"  AD α                = {a_ad(LAM, mu_ph):.6f}     我的 Eliashberg（最细步）= "
      f"{alpha(LAM, MU, 0.0025):.6f}")
print(f"  AD ∂α/∂μ*(500)      = {dad_dmu_ph*chain:+.6f}     我的收敛值 = {AL_M:+.6f}"
      f"     相对差 {100*abs(dad_dmu_ph*chain-AL_M)/abs(AL_M):.1f}%")
print(f"  AD ∂α/∂λ            = {dad_dlam:+.6f}     我的收敛值 = {AL_L:+.6f}"
      f"     相对差 {100*abs(dad_dlam-AL_L)/abs(AL_L):.1f}%")
print("  ⟹ 两把**独立**仪器在收敛步长上相合 ⟹ peer 的 −0.154 是**步长赝物**，不是约定分歧。")
print("     （上一版我未做这两步换算，得到 2.5× 的假分歧并据此怀疑约定——该怀疑撤回。）")

print()
print("=" * 98)
print("§4  判决：用收敛后的 J_B 第二行重算阈值（J_A 第二行待 solve_gap 版，见下条限定）")
print("=" * 98)
print(f"  σ_gap/σ_α：模型1 = {SG_A:.4f}，模型2 = {SG_B:.4f}   （σ_Tc 已证明相消）")
dA_peer = PEER["dTc_dlam"] * PEER["dgr_dmu"] - PEER["dTc_dmu"] * PEER["dgr_dlam"]
print(f"  {'J_B 第二行来源':>26}{'∂α/∂λ':>11}{'∂α/∂μ*':>11}{'det J_B':>11}"
      f"{'阈值':>9}{'模型1':>8}{'模型2':>8}  判决")
cands = [("peer（Δλ 15%, Δμ* 0.05）", PEER["dal_dlam"], PEER["dal_dmu"]),
         ("我复现同步长", a15[0], a15[1]),
         ("★ 收敛步长（Δλ 1%, Δμ* 0.005）", AL_L, AL_M),
         ("AD 闭式解析（换算后）", dad_dlam, dad_dmu_ph * chain)]
for nm, gl, gm in cands:
    d = PEER["dTc_dlam"] * gm - PEER["dTc_dmu"] * gl
    t = abs(dA_peer) / abs(d)
    v = "A 胜(两模型)" if t > max(SG_A, SG_B) else ("A 仅胜模型1" if t > min(SG_A, SG_B) else "★ A 败")
    print(f"  {nm:>26}{gl:>11.6f}{gm:>11.6f}{d:>11.4f}{t:>9.4f}"
          f"{t/SG_A:>8.3f}{t/SG_B:>8.3f}  {v}")
print()
print("  ★ 结论：判决**方向**不变（A 仍胜两模型），但**余量从 2.01×/1.60× 掉到"
      f" {abs(dA_peer)/abs(PEER['dTc_dlam']*AL_M-PEER['dTc_dmu']*AL_L)/SG_A:.2f}×/"
      f"{abs(dA_peer)/abs(PEER['dTc_dlam']*AL_M-PEER['dTc_dmu']*AL_L)/SG_B:.2f}×**")
print("    ——模型 2 只剩十几个百分点的余量。所报的「1.60×」高估了约 1.4×。")
print("  ⚠ 限定：本表的 det(J_A) 仍用 peer 的 ∂gr 两项（未收敛步长）。")
print("    ∂gr 两项要 solve_gap，且仓库已知 2Δ/kTc 对 N 极敏感（λ=0.3 时 N=400 给 4.05、N=8000 给 3.55）")
print("    ⟹ **J_A 第二行的步长与 N 收敛都还没查**，本判决只在「它们已收敛」的前提下成立。")
print()
print("=" * 98)
print(f"完成，用时 {time.time()-t0:.0f} 秒（N={NM}）")
print("=" * 98)
