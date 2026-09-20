# -*- coding: utf-8 -*-
"""第三步：q 的**约定敏感度**与**标度不变性**。

三个必须回答的问题（前两个是门，第三个是结果）：

门 A [analytic]：μ*=0 时各向同性 Eliashberg **无绝对能标**——把整条 α²F 的频率轴乘 c，
  Tc 也乘 c，ω̄₂ 乘 c ⟹ 𝔣 **逐位不变**。这保证「用 ω̄₂~100 meV 的氢化物谱测出的 q」
  可以拿去和「ω̄₂~5 meV 的 AD75 材料」比。若实现破坏它，整条比较就不成立。

门 B：μ*>0 时唯一的无量纲参数是 **w_cut/ω̄₂**（不是 w_cut 本身）。
  ⟹ 上一轮的「约定 II：w_cut=500 meV 绝对」对 AD75 尺度的材料是 **w_cut/ω̄₂ ≈ 100**，
  而对氢化物是 ≈ 3–5。两者不是同一个物理设定。本脚本改用比值扫描。

结果：q(μ*=0.10) 对 w_cut/ω̄₂ 的依赖有多强？两种 μ* 约定：
  (a) μ* 固定 = 0.10（朴素：把 AD75 印出的数直接当成任何截断下的 μ*）；
  (b) μ* 跑动：以 w_cut/ω̄₂ = 10 处 μ*=0.10 为参考，按 Anderson–Morel
      1/μ*(ω_c') = 1/μ*(ω_c) + ln(ω_c/ω_c') 换算（仓库 `mu_star_at_cutoff`）。

严格性：numerical-evidence（门 A 的命题本身是 analytic，这里只检验实现）。
主会话 2026-09-20，从零。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "solvers"))
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K, mu_star_at_cutoff  # noqa: E402

C_INF = 0.1827262477
DATA = ROOT / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"
FILES = {
    "H3S-anh": DATA / "H3S-200GPa/EPW_anharmonic/H3S.a2f",
    "D3S-anh": DATA / "D3S-200GPa/EPW_anharmonic/D3S.a2f",
    "H3S-har": DATA / "H3S-200GPa/EPW_harmonic/H3S.a2f",
    "D3S-har": DATA / "D3S-200GPa/EPW_harmonic/D3S.a2f",
}
NMATS = 400
LAM_GRID = np.exp(np.linspace(np.log(0.69), np.log(2.59), 9))


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

# ---------------------------------------------------------------- 门 A
print("=" * 88)
print("门 A [analytic 命题的实现检验]：μ*=0 时 𝔣 对整体频率标度 c 逐位不变")
print("  （AD75 材料 ω̄₂ ≈ 4–13 meV，氢化物 ≈ 93–149 meV，差 ~20× ⟹ 这道门决定两者可不可比）")
w, a = RAW["H3S-anh"]
base = None
for c in (1.0, 0.05, 0.2, 5.0):
    sp = Spectrum.from_table(w * c, a)          # α²F 无量纲 ⟹ 只缩频率轴，λ 不变
    tc = Solver(spectrum=sp, mu_star=0.0, w_cut=1e4 * c, n_matsubara=NMATS).tc()
    f = tc / (C_INF * np.sqrt(sp.lambda_total()) * sp.omega_2() / KB_MEV_PER_K)
    if base is None:
        base = f
    print(f"  c={c:5.2f}:  ω̄₂ = {sp.omega_2():8.3f} meV   Tc = {tc:10.4f} K   "
          f"𝔣 = {f:.10f}   |𝔣/𝔣₀−1| = {abs(f/base-1):.2e}")
print()


def scan(mu_fn, wratio):
    """mu_fn(w2) -> μ*；w_cut = wratio·ω̄₂。返回 (pooled b,c,q, err, q 范围, 𝔣max)。"""
    rows = []
    for k, (w, a) in RAW.items():
        sp_ref = Spectrum.from_table(w, a)
        lam_ref, w2, wl = sp_ref.lambda_total(), sp_ref.omega_2(), sp_ref.omega_log()
        wcut = wratio * w2
        mu = mu_fn(w2)
        for lam in LAM_GRID:
            sp = Spectrum.from_table(w, a * (lam / lam_ref))
            tc = Solver(spectrum=sp, mu_star=mu, w_cut=wcut, n_matsubara=NMATS).tc()
            assert tc > 0, f"{k} λ={lam} 给非正 Tc"
            ceil_K = C_INF * np.sqrt(sp.lambda_total()) * w2 / KB_MEV_PER_K
            rows.append((k, wl / w2, sp.lambda_total(), tc / ceil_K))
    arr_k = np.array([r[0] for r in rows])
    S = np.array([r[1] for r in rows]); lam = np.array([r[2] for r in rows])
    f = np.array([r[3] for r in rows])
    x, y = np.log(lam), np.log(f)
    qs = []
    for k in FILES:
        m = arr_k == k
        A = np.vstack([np.ones(m.sum()), x[m], x[m] ** 2]).T
        qs.append(np.linalg.lstsq(A, y[m], rcond=None)[0][2])
    A = np.vstack([np.ones(len(x)), x, np.log(S), x ** 2]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ coef
    cov = (res @ res / (len(x) - 4)) * np.linalg.inv(A.T @ A)
    return coef, np.sqrt(np.diag(cov)), (min(qs), max(qs)), f.max()


print("=" * 88)
print("结果：q(μ*≈0.10) 对无量纲截断 w_cut/ω̄₂ 的依赖")
print("对照：AD75 29 行经验回归 b = +0.9972 ± 0.041，c = +0.7198，q = −0.525")
print()
MU_REF, R_REF = 0.10, 10.0
print(f"{'w_c/ω̄₂':>8} | {'μ*(a) 固定':>10} {'b':>8} {'q':>9} {'±':>6} | "
      f"{'μ*(b) 跑动':>10} {'b':>8} {'q':>9} {'±':>6}")
for r in (3.0, 5.0, 10.0, 20.0, 50.0, 100.0):
    out = []
    for mode in ("fix", "run"):
        if mode == "fix":
            mu_fn = lambda w2: MU_REF
            mu_show = MU_REF
        else:
            # μ* 从参考截断 R_REF·ω̄₂ 跑到 r·ω̄₂；比值 ⟹ 与 ω̄₂ 无关
            mu_show = mu_star_at_cutoff(MU_REF, R_REF, r)
            mu_fn = lambda w2, m=mu_show: m
        coef, err, qr, fmax = scan(mu_fn, r)
        out.append((mu_show, coef[1], coef[3], err[3], fmax))
    a1, b1, q1, e1, fm1 = out[0]
    a2, b2, q2, e2, fm2 = out[1]
    print(f"{r:8.1f} | {a1:10.4f} {b1:+8.4f} {q1:+9.4f} {e1:6.4f} | "
          f"{a2:10.4f} {b2:+8.4f} {q2:+9.4f} {e2:6.4f}"
          + ("   <<< 𝔣>1" if max(fm1, fm2) > 1 else ""))

print()
print("读法：(a) 列把 AD75 印出的 μ*=0.10 当作任何截断下的值——**约定不自洽**，只作上下包络用；")
print("      (b) 列固定裸库仑、让 μ* 随截断跑（Anderson–Morel）——物理上自洽的那一列。")
print("严格性：全部 numerical-evidence；门 A 的不变性命题本身是 analytic。")
