# -*- coding: utf-8 -*-
"""α(μ*) 与 α(λ) 到底光不光滑？—— 决定 J_B 第二行能不能用差分取到。

背景：`_chk_main_alpha_jacobian_step_convergence` 的步长扫描**不收敛**（震荡），
我在那个脚本 §3 里预写的「两把仪器相合」判词与印出的 103%/329% 自相矛盾 ⟹ 该判词作废。
不收敛有两种完全不同的成因，后果相反：
  (a) 舍入主导 —— 真导数存在，只是 h 太小；⟹ 取最优 h 或拟合即可救。
  (b) 函数本身不可微（扭结/跳变）—— ⟹ 差分取不到任何东西，J_B 第二行**不存在**。
不先分清就报数，就是把「数字错」和「证书松」混为一谈。

线索：Solver 第 565–568 行对库仑截断装了**线性斜坡** `clip((w_cut−ω_n)/(2πT)+0.5,0,1)`，
注释自陈是为了压住「某些巧合比例下 α 跳变 ~0.03」。斜坡使函数 C⁰ 但**不是 C¹**。
锚点处 2πTc ≈ 102 meV、w_cut = 500 meV ⟹ 截断内只有约 **4.9** 个松原频率，
斜坡一格就占库仑权重的 ~20% ⟹ 扭结不是小量。

本脚本不推理，直接**密采样**看曲线。
方向注记：全为点算，无单边界放缩。
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
from eliashberg_iso import Spectrum, Solver   # noqa: E402

LAM, WE, MU, WCUT, NM = 1.876, 106.9, 0.13, 500.0, 400
KB = 0.08617333262


def tc(lam, mu, we=WE):
    return Solver(spectrum=Spectrum.einstein(lam=lam, w0=we), mu_star=mu,
                  w_cut=WCUT, n_matsubara=NM).tc()


def alpha(lam, mu, eps=0.01):
    a, b = np.log(tc(lam, mu, WE * (1 + eps))), np.log(tc(lam, mu, WE * (1 - eps)))
    return 0.5 * (a - b) / np.log((1 + eps) / (1 - eps))


def ncut(t_k):
    """截断内的松原频率数 w_cut/(2πT)（斜坡中心）。"""
    return WCUT / (2 * np.pi * t_k * KB)


t0 = time.time()
print("=" * 98)
print("§0  锚点：截断里装得下几个松原频率？")
print("=" * 98)
t_anch = tc(LAM, MU)
print(f"  Tc = {t_anch:.4f} K，2πTc = {2*np.pi*t_anch*KB:.3f} meV，w_cut = {WCUT} meV")
print(f"  ⟹ w_cut/(2πTc) = {ncut(t_anch):.4f}  —— 截断内只有约 5 个频率，斜坡一格占 "
      f"{100/ncut(t_anch):.0f}% 库仑权重")
print(f"  距最近的整数边界 {min(ncut(t_anch)%1, 1-ncut(t_anch)%1):.4f} 格")

for tag, var, grid, other in (
        ("μ*", "mu", np.linspace(0.100, 0.160, 31), LAM),
        ("λ", "lam", np.linspace(1.60, 2.15, 31), MU)):
    print()
    print("=" * 98)
    print(f"§  α 对 {tag} 的密采样（31 点）")
    print("=" * 98)
    A, T = [], []
    for g in grid:
        A.append(alpha(other, g) if var == "mu" else alpha(g, other))
        T.append(tc(other, g) if var == "mu" else tc(g, other))
    A, T = np.array(A), np.array(T)
    h = grid[1] - grid[0]
    d1 = np.gradient(A, h)
    nc = ncut(T)
    print(f"  {tag:>8}{'Tc(K)':>11}{'wc/2πT':>9}{'α':>12}{'dα/d'+tag:>12}")
    for i in range(0, 31, 2):
        print(f"  {grid[i]:>8.4f}{T[i]:>11.3f}{nc[i]:>9.4f}{A[i]:>12.6f}{d1[i]:>12.5f}")
    # 光滑性诊断：二阶差分的逐点幅度 vs 其自身中位数
    d2 = np.abs(np.diff(A, 2))
    med = float(np.median(d2))
    print(f"  ★ 二阶差分 |Δ²α|：中位 {med:.3e}，最大 {d2.max():.3e}，"
          f"最大/中位 = {d2.max()/med:.1f}×  @ {tag}={grid[np.argmax(d2)+1]:.4f}")
    xs = [j for j in range(len(nc) - 1) if np.floor(nc[j]) != np.floor(nc[j + 1])]
    print(f"     wc/2πT 从 {nc.min():.3f} 到 {nc.max():.3f}；整数边界穿越 {len(xs)} 次"
          + ("".join(f"（@{tag}≈{0.5*(grid[j]+grid[j+1]):.4f}，穿 {np.floor(nc[j+1]):.0f}）"
                     for j in xs) if xs else "（无）"))
    smooth = d2.max() / med < 5.0
    print(f"     判决：{'光滑，导数可由拟合稳健取得' if smooth else '★ 不光滑（有扭结）——局部导数与周期平均导数不是同一个数'}")
    # 三个互不相同的「导数」，都印出来：局部干净支 / 整段割线 / 多项式拟合
    i0 = int(np.argmin(np.abs(grid - (MU if var == "mu" else LAM))))
    lo = max(0, i0 - 2); hi = min(len(grid) - 1, i0 + 2)
    print(f"     ① 局部干净支斜率（{tag}={grid[lo]:.4f}→{grid[hi]:.4f}，不跨边界）= "
          f"{(A[hi]-A[lo])/(grid[hi]-grid[lo]):+.6f}")
    print(f"     ② 整段割线（跨 {len(xs)} 个边界，≈周期平均）            = "
          f"{(A[-1]-A[0])/(grid[-1]-grid[0]):+.6f}")
    for deg in (2, 3, 4):
        c = np.polyfit(grid, A, deg)
        dv = np.polyval(np.polyder(c), MU if var == "mu" else LAM)
        r = A - np.polyval(c, grid)
        print(f"     ③ {deg} 次拟合 dα/d{tag} = {dv:+.6f}   残差 rms = {np.sqrt(r@r/len(r)):.2e}")

print()
print("=" * 98)
print(f"完成，用时 {time.time()-t0:.0f} 秒（N={NM}）")
print("=" * 98)
