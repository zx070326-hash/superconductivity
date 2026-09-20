# -*- coding: utf-8 -*-
"""主会话从零独立测量「理论侧曲率指数 q」——不经 AD 闭式、不经 AD75 数据。

动机（CLAUDE.md「承重结论必须从零独立重算」）：
  头条回归 `log𝔣 = a + b·logλ + c·logS + q·(logλ)²` 在 AD75 29 行上给 q = −0.525，
  但该点值**不可引用**（同材料双技术偏置的不确定度比 q 本身还宽，需 n≥24 对，现只有 2 对）。
  若**理论**能预言 q，数据就从「拟合 q」变成「被 q 检验」。
  理论侧不存在技术偏置 r₁ ⟹ 这条路绕开了那个阻断。

方法：
  振幅重标度 α²F → s·α²F ⟹ λ → sλ，而 ω_log 与 ω̄₂ **精确不变**
  （两者都是 dλ 归一的矩，s 在分子分母里约掉）⟹ 形状坐标 S := ω_log/ω̄₂ 严格钉死。
  于是 `log𝔣 vs logλ` 的曲率**纯粹**是 λ 的效应，不混形状。

谱：Mishra et al. (Figshare 31136635) 的四条真实 α²F —— H3S/D3S × 谐波/非谐 @200 GPa。
  它们的 S ∈ [0.737, 0.856] 全部落在 AD75 支撑域 [0.5446, 0.8919] **内**。
  ⚠ 不用 Einstein 谱作主线：S = 1 落在支撑域**外**，正是 CaC₆ 那条线死掉的原因。

μ*：主线 μ* = 0（此时 Tc 对 ω_c 严格免疫 [analytic]，见 h3s_mustar_cutoff_convention）。

严格性：float64 Eliashberg 数值解 ⟹ 全部结论为 numerical-evidence。
  「ω_log/ω̄₂ 在振幅重标度下不变」是 analytic，脚本里作为 1e-12 断言门检验实现。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, time
from pathlib import Path
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = Path(__file__).resolve().parents[2]          # theory/
sys.path.insert(0, str(ROOT / "solvers"))
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K   # noqa: E402

C_INF = 0.1827262477
DATA = ROOT / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"
FILES = {
    "H3S-anh": DATA / "H3S-200GPa/EPW_anharmonic/H3S.a2f",
    "D3S-anh": DATA / "D3S-200GPa/EPW_anharmonic/D3S.a2f",
    "H3S-har": DATA / "H3S-200GPa/EPW_harmonic/H3S.a2f",
    "D3S-har": DATA / "D3S-200GPa/EPW_harmonic/D3S.a2f",
}
# 仓库既有值（验收门，来自 deficit_budget_* / h3s_on_ad75_trend_closed 笔记）
REPO = {  # name: (lambda, S = w_log/w2)
    "H3S-anh": (1.7210, 0.7369), "D3S-anh": (1.7838, 0.8113),
    "H3S-har": (2.2929, 0.8287), "D3S-har": (2.2938, 0.8557),
}
WCUT = 500.0        # meV；μ*=0 时不进核，仅用于 Solver 构造
NMATS = 400


def load(path):
    w, a = [], []
    for line in path.read_text().splitlines()[1:]:
        p = line.split()
        if len(p) < 3:
            break
        try:
            x, y = float(p[0]), float(p[1])
        except ValueError:
            break
        w.append(x); a.append(max(y, 0.0))
    return np.array(w), np.array(a)


RAW = {k: load(p) for k, p in FILES.items()}

print("=" * 78)
print("门 0：谱矩复现（对照仓库既有 λ 与 S）")
print(f"{'谱':9s} {'λ_me':>8} {'λ_repo':>8} {'差%':>7}   {'S_me':>7} {'S_repo':>7} {'差%':>7}  {'ω̄₂(meV)':>9}")
g0 = True
for k, (w, a) in RAW.items():
    sp = Spectrum.from_table(w, a, name=k)
    lam, wl, w2 = sp.lambda_total(), sp.omega_log(), sp.omega_2()
    S = wl / w2
    dl = 100 * (lam / REPO[k][0] - 1); ds = 100 * (S / REPO[k][1] - 1)
    print(f"{k:9s} {lam:8.4f} {REPO[k][0]:8.4f} {dl:+7.3f}   {S:7.4f} {REPO[k][1]:7.4f} {ds:+7.3f}  {w2:9.3f}")
    if abs(dl) > 0.2 or abs(ds) > 0.2:
        g0 = False
print("门 0:", "PASS" if g0 else "*** FAIL ***")

print("\n" + "=" * 78)
print("门 1 [analytic 命题的实现检验]：振幅重标度下 ω_log 与 ω̄₂ 必须逐位不变")
w, a = RAW["H3S-anh"]
sp0 = Spectrum.from_table(w, a)
g1 = True
for s in (0.4, 1.0, 1.5, 3.0):
    sp = Spectrum.from_table(w, s * a)
    e_lam = abs(sp.lambda_total() / (s * sp0.lambda_total()) - 1)
    e_wl = abs(sp.omega_log() / sp0.omega_log() - 1)
    e_w2 = abs(sp.omega_2() / sp0.omega_2() - 1)
    print(f"  s={s:4.1f}:  |λ/(sλ₀)−1| = {e_lam:.2e}   |ω_log 变化| = {e_wl:.2e}   |ω̄₂ 变化| = {e_w2:.2e}")
    if max(e_lam, e_wl, e_w2) > 1e-12:
        g1 = False
print("门 1:", "PASS (≤1e-12)" if g1 else "*** FAIL ***")


def tc_of(sp, mu, n=NMATS, wcut=WCUT):
    return Solver(spectrum=sp, mu_star=mu, w_cut=wcut, n_matsubara=n).tc()


print("\n" + "=" * 78)
print("门 2：H3S-非谐 原谱在 μ*=0 的 Tc 与天花板（仓库：Tc(N=400)=283.3462 K，ceiling=415.21 K，𝔣=0.6824）")
lam0, w2_0 = sp0.lambda_total(), sp0.omega_2()
ceiling_K = C_INF * np.sqrt(lam0) * w2_0 / KB_MEV_PER_K
t0 = time.time()
tc0 = tc_of(sp0, 0.0)
print(f"  Tc(μ*=0, N=400) = {tc0:.4f} K   （仓库 283.3462，差 {100*(tc0/283.3462-1):+.4f}%，{time.time()-t0:.1f}s）")
print(f"  ceiling         = {ceiling_K:.4f} K   （仓库 415.21，差 {100*(ceiling_K/415.21-1):+.3f}%）")
print(f"  𝔣(μ*=0)         = {tc0/ceiling_K:.4f}   （仓库 283.3462/415.21 = {283.3462/415.21:.4f}）")

print("\n门 3：μ*=0 时 Tc 对 w_cut 严格免疫 [analytic]，实现必须逐位复现")
for wc in (200.0, 500.0, 1000.0, 2000.0):
    print(f"  w_cut={wc:7.1f} meV -> Tc = {tc_of(sp0, 0.0, wcut=wc):.6f} K")

print("\n门 4：N 阶梯（μ*=0，H3S-非谐）")
for n in (100, 200, 400, 800):
    print(f"  N={n:4d} -> Tc = {tc_of(sp0, 0.0, n=n):.4f} K")

# ------------------------------------------------------------------ 主扫描
LAM_GRID = np.exp(np.linspace(np.log(0.69), np.log(2.59), 9))   # AD75 支撑域
print("\n" + "=" * 78)
print(f"主扫描：μ*=0，N={NMATS}，λ 网格 = AD75 支撑域 [0.69, 2.59] 九点（对数等距）")
print("       振幅重标度 ⟹ 每条谱的 S 严格不动\n")

rows = []
for k, (w, a) in RAW.items():
    sp_ref = Spectrum.from_table(w, a)
    lam_ref, w2, wl = sp_ref.lambda_total(), sp_ref.omega_2(), sp_ref.omega_log()
    S = wl / w2
    print(f"--- {k}  (S = {S:.4f}, ω̄₂ = {w2:.2f} meV) ---")
    print(f"{'λ':>7} {'Tc(K)':>10} {'ceiling(K)':>11} {'𝔣':>8}")
    for lam in LAM_GRID:
        sp = Spectrum.from_table(w, a * (lam / lam_ref))
        lam_a = sp.lambda_total()
        tc = tc_of(sp, 0.0)
        ceil_K = C_INF * np.sqrt(lam_a) * w2 / KB_MEV_PER_K
        f = tc / ceil_K
        flag = "  <<< 𝔣>1 报警" if f > 1.0 else ""
        print(f"{lam_a:7.4f} {tc:10.4f} {ceil_K:11.3f} {f:8.4f}{flag}")
        rows.append((k, S, lam_a, tc, f))
    print()

# ------------------------------------------------------------------ 拟合
arr_k = np.array([r[0] for r in rows])
lam = np.array([r[2] for r in rows]); f = np.array([r[4] for r in rows])
S_all = np.array([r[1] for r in rows])
x = np.log(lam); y = np.log(f)

print("=" * 78)
print("拟合 1（逐谱，S 固定 ⟹ 无需 c 项）：log𝔣 = a + b·logλ + q·(logλ)²")
print(f"{'谱':9s} {'a':>9} {'b':>9} {'q':>9} {'R²':>8} {'rms残差':>9}")
qs = []
for k in FILES:
    m = arr_k == k
    A = np.vstack([np.ones(m.sum()), x[m], x[m] ** 2]).T
    coef, *_ = np.linalg.lstsq(A, y[m], rcond=None)
    res = y[m] - A @ coef
    r2 = 1 - res @ res / ((y[m] - y[m].mean()) @ (y[m] - y[m].mean()))
    print(f"{k:9s} {coef[0]:9.4f} {coef[1]:9.4f} {coef[2]:9.4f} {r2:8.5f} {np.sqrt((res**2).mean()):9.2e}")
    qs.append(coef[2])
qs = np.array(qs)
print(f"\n  q 的四谱范围 [{qs.min():.4f}, {qs.max():.4f}]，均值 {qs.mean():.4f}，跨度 {qs.max()-qs.min():.4f}")

print("\n拟合 2（合并四谱，带 logS 项，与头条回归同形）：log𝔣 = a + b·logλ + c·logS + q·(logλ)²")
A = np.vstack([np.ones(len(x)), x, np.log(S_all), x ** 2]).T
coef, *_ = np.linalg.lstsq(A, y, rcond=None)
res = y - A @ coef
r2 = 1 - res @ res / ((y - y.mean()) @ (y - y.mean()))
dof = len(x) - 4
s_unb = np.sqrt(res @ res / dof)
cov = s_unb ** 2 * np.linalg.inv(A.T @ A)
print(f"  a = {coef[0]:+.4f}   b = {coef[1]:+.4f} ± {np.sqrt(cov[1,1]):.4f}"
      f"   c = {coef[2]:+.4f} ± {np.sqrt(cov[2,2]):.4f}"
      f"   q = {coef[3]:+.4f} ± {np.sqrt(cov[3,3]):.4f}")
print(f"  R² = {r2:.5f}   无偏残差 s = {s_unb:.4f} (dof={dof})")
print(f"\n  对照：AD75 29 行经验回归给 (a,b,c) = (−1.0691, +0.9972, +0.7198)、q = −0.525、s = 0.0794")
print(f"  ⚠ 严格性：本脚本全部数值结论为 numerical-evidence（float64 Eliashberg + 有限 λ 网格）。")
print(f"  ⚠ 𝔣 最大值 = {f.max():.4f}" + ("  <<< 破天花板，必须追查" if f.max() > 1 else "  （≤1，天花板未破）"))
