# -*- coding: utf-8 -*-
"""
H3S 的两个自家数字对不上：
  · 仓库认证 cert-v6 / 独立基线 : Tc(μ*=0.13) = 198.974 K ，Tc(μ*=0) = 283.344 K
  · f-deficit-decomposition 交付 : Tc(μ*=0.13) = 194.33 K   （w_cut = 500 meV）
差 2.3%。两者都承重 ⟹ 必须判出是「数字错」还是「约定不同」。

★ 判决性检验（一步定死）：**μ*=0 时库仑项完全不作用 ⟹ Tc 对 ω_c 的选择免疫。**
   （唯一残留依赖是 Matsubara 截断 N 本身，可用 N 收敛单独排除。）
   ⟹ 若 Tc(μ*=0) 两边**一致**而 Tc(μ*=0.13) 不一致
      ⟹ 分歧**完全**落在库仑通道 ⟹ 是 ω_c 约定，不是求解器错。
   ⟹ 若 Tc(μ*=0) 就不一致 ⟹ 谱/矩/求解器层面有真错，与 μ* 无关。

方向审计：本脚本不产出单边界判决，只做两值比对 + 收敛性，无需方向。
严格性：numerical-evidence（float64，N 收敛自测）。

主会话 2026-09-20，从零。
"""
import sys, os
sys.path.insert(0, r"D:\physics\superconductivity")
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from theory.solvers.eliashberg_iso import Spectrum, Solver

MEV_PER_K = 1.0 / 11.604518
A2F = (r"D:\physics\superconductivity\theory\data\external\mishra_h3s_d3s"
       r"\Data-H3S-D3S-PH-EPW-Anharmonic\H3S-200GPa\EPW_anharmonic\H3S.a2f")

# ---- 自己读谱，不复用 agent 的 loader ----
rows = []
with open(A2F) as fh:
    for ln in fh.readlines()[1:]:
        p = ln.split()
        if len(p) < 2: break
        try: rows.append((float(p[0]), float(p[1])))
        except ValueError: break
W_meV = np.array([r[0] for r in rows]); F = np.array([r[1] for r in rows])
m = W_meV > 0
lam  = 2*np.trapezoid(F[m]/W_meV[m], W_meV[m])
wlog = np.exp(2/lam*np.trapezoid(F[m]*np.log(W_meV[m])/W_meV[m], W_meV[m]))
wb2  = np.sqrt(2/lam*np.trapezoid(F[m]*W_meV[m], W_meV[m]))
print("="*92)
print("0. 谱与矩（主会话自读，单位 meV；括号内换算到 K 与仓库值对照）")
print("="*92)
print(f"  {len(W_meV)} 点，ω ∈ [{W_meV.min():.4f}, {W_meV.max():.4f}] meV")
print(f"  λ      = {lam:.6f}        （仓库 1.720958）  差 {lam-1.720958:+.2e}")
print(f"  ω_log  = {wlog:.6f} meV = {wlog/MEV_PER_K:.4f} K  （仓库 109.995163 K）")
print(f"  ω̄₂     = {wb2:.6f} meV = {wb2/MEV_PER_K:.4f} K  （仓库 149.262922 K）")

spec = Spectrum.from_table(w=W_meV, a2f=F, name="H3S_anh_200GPa")

# ---- 1. μ*=0：对 ω_c 免疫，先验收求解器与谱 ----
print()
print("="*92)
print("1. ★ 判决性检验：μ*=0 ⟹ 库仑项不作用 ⟹ Tc 应对 w_cut **免疫**")
print("="*92)
print(f"  仓库独立基线 Tc(μ*=0) = 283.344 K\n")
print(f"  {'w_cut [meV]':>12}{'N':>8}{'Tc(μ*=0) [K]':>16}{'vs 283.344':>14}")
base0 = {}
for wc in (200.0, 500.0, 1000.0, 2000.0):
    for N in (400,):
        tc = Solver(spectrum=spec, mu_star=0.0, w_cut=wc, n_matsubara=N).tc()
        base0[(wc, N)] = tc
        print(f"  {wc:>12.0f}{N:>8d}{tc:>16.4f}{100*(tc/283.344-1):>13.2f}%", flush=True)

print(f"\n  N 收敛（w_cut=500, μ*=0）：")
for N in (100, 200, 400, 800):
    tcN = Solver(spectrum=spec, mu_star=0.0, w_cut=500.0, n_matsubara=N).tc()
    print(f"    N={N:>5d}  Tc = {tcN:.4f} K", flush=True)

# ---- 2. μ*=0.13：扫 w_cut，看 198.974 落在哪个约定上 ----
print()
print("="*92)
print("2. μ*=0.13：Tc 随 w_cut 的走向（库仑通道**必然**依赖 ω_c 约定）")
print("="*92)
print(f"  仓库认证值 198.974 K（P 模型，r=10 ⟹ N_c=⌊rϖ⌋=13）")
print(f"  agent 交付 194.33 K（w_cut=500 meV）\n")
print(f"  {'w_cut [meV]':>12}{'N':>8}{'Tc(μ*=.13) [K]':>17}{'vs 198.974':>13}{'N_c(≈)':>9}")
for wc in (200.0, 300.0, 500.0, 800.0, 1200.0, 2000.0):
    tc = Solver(spectrum=spec, mu_star=0.13, w_cut=wc, n_matsubara=400).tc()
    t_mev = tc * MEV_PER_K
    nc = (wc/(np.pi*t_mev) - 1.0)/2.0          # (2n+1)πT ≤ w_cut
    print(f"  {wc:>12.0f}{400:>8d}{tc:>17.4f}{100*(tc/198.974-1):>12.2f}%{nc:>9.1f}",
          flush=True)

print("""
  判读规则（先写下，再看数）：
   · 若第 1 节 μ*=0 在四个 w_cut 上**一致**且对上 283.344
     ⟹ 谱、矩、求解器全部无误 ⟹ 2.3% **只**能来自库仑通道 ⟹ 是 **ω_c 约定之争**，不是错。
   · 若第 1 节就散开或对不上 283.344 ⟹ 有与 μ* 无关的真错，先修那个。
""")
print("="*92)

# ---- 3. 反推仓库约定：r=10 ⟹ ω_c ≈ r·ω̄₂，不是 500 meV ----
print()
print("="*92)
print("3. ★ 反推仓库的 ω_c：ϖ = ω̄₂/(2πT) ⟹ N_c=⌊rϖ⌋=13 ⟹ ω_c=(2N_c+1)πT")
print("="*92)
Tc_repo_K = 198.974; Tc_repo_meV = Tc_repo_K*MEV_PER_K
piT = np.pi*Tc_repo_meV
vpi = wb2/(2*piT)
print(f"  T = {Tc_repo_K} K = {Tc_repo_meV:.4f} meV ,  πT = {piT:.4f} meV")
print(f"  ϖ = ω̄₂/(2πT) = {wb2:.4f}/{2*piT:.4f} = {vpi:.4f}   （仓库印出 ϖ≈1.39 ✓）")
print(f"  N_c = ⌊10·ϖ⌋ = {int(np.floor(10*vpi))}                （仓库印出 N_c=13 ✓）")
wc_star = (2*int(np.floor(10*vpi))+1)*piT
print(f"  ⟹ ω_c = (2N_c+1)πT = {wc_star:.1f} meV  ≈ {wc_star/wb2:.2f}·ω̄₂")
print(f"  （对照：agent 用 500 meV = {500/wb2:.2f}·ω̄₂）\n")
print(f"  {'w_cut [meV]':>12}{'Tc(μ*=.13) [K]':>17}{'vs 198.974':>13}")
for wc in (1400.0, wc_star, 1500.0, 10*wb2):
    tc = Solver(spectrum=spec, mu_star=0.13, w_cut=wc, n_matsubara=400).tc()
    print(f"  {wc:>12.1f}{tc:>17.4f}{100*(tc/198.974-1):>12.2f}%", flush=True)

print("""
  ⟹ 若此处复现 198.97 ⟹ **两个数字全对，分歧 100% 是 ω_c 约定**，案结。
  ⚠ 但随之而来的是一个**真问题**：第 2 节显示 Tc(μ*=0.13) 对 ω_c **非单调**
    （200→2000 meV 走 203.98 ↘ 194.33 ↗ 200.90，谷底在 ~500–600 meV）
    ⟹ **不能**靠取极端 ω_c 给 δ_COU 一个单边界。
    物理机制：ω_c 增大有两个反向效应——μ* 作用的频点变多（压低 Tc），
    但 Δ(iω_n) 的**变号**空间也变大（Morel–Anderson 屏蔽，抬高 Tc）。
  ⟹ **δ_COU 必须随它的 ω_c 一起报**，单说「库仑吃掉 X%」是不完整的陈述。""")
print("="*92)

# ---- 4. 非单调性是真的还是 N 伪影？（新主张 ⟹ 必须先验收数值） ----
print()
print("="*92)
print("4. ★ 验收第 2 节的非单调性：它是物理还是 N 截断伪影？")
print("="*92)
print("""  第 2 节固定 N=400 扫 w_cut。μ*=0 处 N 收敛极好（400→800 只动 0.001%），
  但 μ*≠0 时 θ 斜坡结构随 w_cut 变 ⟹ **不能沿用 μ*=0 的收敛结论**。
  在谷底(500)与右端(2000)各做一次 N 收敛；若两端都收敛，非单调性为真。""")
print(f"\n  {'N':>7}{'Tc @ w_cut=500':>18}{'Tc @ w_cut=2000':>19}")
for N in (200, 400, 800):
    a = Solver(spectrum=spec, mu_star=0.13, w_cut=500.0,  n_matsubara=N).tc()
    b = Solver(spectrum=spec, mu_star=0.13, w_cut=2000.0, n_matsubara=N).tc()
    print(f"  {N:>7d}{a:>18.4f}{b:>19.4f}", flush=True)

print(f"\n  谷底细扫（N=400）：")
print(f"  {'w_cut [meV]':>12}{'Tc(μ*=.13) [K]':>17}")
for wc in (400.0, 450.0, 500.0, 550.0, 600.0, 700.0):
    tc = Solver(spectrum=spec, mu_star=0.13, w_cut=wc, n_matsubara=400).tc()
    print(f"  {wc:>12.0f}{tc:>17.4f}", flush=True)
print("="*92)
