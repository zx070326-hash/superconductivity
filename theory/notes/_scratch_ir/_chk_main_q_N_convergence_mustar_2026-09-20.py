# -*- coding: utf-8 -*-
"""核实 f-sample-expansion 的 q-vs-N 表，并补上它**没测**的那一档。

它交付的是 **μ*=0** 的 N 阶梯（H3S-anh，w_cut=500 meV）：
    N     100        200        400        800        1600
    q  -0.328672 -0.329274 -0.329425 -0.329463 -0.329473
    ⟹ N=400→1600 变化 -4.7e-05（0.014%）⟹ CONVERGED

**缺口**：我的生产值 `q = −0.5565` 在 **μ*=0.10** 上。库仑通道引入 `w_cut` 截断与
对 `|n| ≤ N_c` 的求和，其 N 行为与 μ*=0 **不是同一件事**——μ*=0 的收敛**不蕴含** μ*>0 的收敛。
（同族教训：`nphi` 那次正是「先测 N 收敛再判结论」。）

故本脚本：
  §1  从零重跑 μ*=0 的阶梯（对照它的五个数）
  §2  ★ 补跑 **μ*=0.10** 的阶梯——这才是我承重值所在的那一档
  §3  ★ 再补四谱合并的 μ*=0.10 阶梯（生产量本身，不是单谱）

方向注记：双侧点算，不产生单边界判决。
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
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K   # noqa: E402

C_INF = 0.1827262477
LAM_GRID = np.exp(np.linspace(np.log(0.69), np.log(2.59), 9))
NS = (100, 200, 400, 800, 1600)
DATA = ROOT / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"
FILES = {"H3S-anh": DATA / "H3S-200GPa/EPW_anharmonic/H3S.a2f",
         "D3S-anh": DATA / "D3S-200GPa/EPW_anharmonic/D3S.a2f",
         "H3S-har": DATA / "H3S-200GPa/EPW_harmonic/H3S.a2f",
         "D3S-har": DATA / "D3S-200GPa/EPW_harmonic/D3S.a2f"}


def load(p):
    w, a = [], []
    for line in p.read_text().splitlines()[1:]:
        q = line.split()
        if len(q) < 3:
            break
        try:
            x, yv = float(q[0]), float(q[1])
        except ValueError:
            break
        w.append(x); a.append(max(yv, 0.0))
    return np.array(w), np.array(a)


RAW = {k: load(p) for k, p in FILES.items()}


def fit(L, Sg, y, with_S):
    cols = [np.ones(len(L)), L] + ([Sg] if with_S else []) + [L * L]
    X = np.column_stack(cols)
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ c
    return c, float(np.sqrt(r @ r / len(y)))


def run(keys, mu, n, wcut_abs=None):
    L, S, Y = [], [], []
    for k in keys:
        w, a = RAW[k]
        sp0 = Spectrum.from_table(w, a)
        l0, w2, wl = sp0.lambda_total(), sp0.omega_2(), sp0.omega_log()
        wc = wcut_abs if wcut_abs is not None else 10.0 * w2
        for lv in LAM_GRID:
            sp = Spectrum.from_table(w, a * (lv / l0))
            t = Solver(spectrum=sp, mu_star=mu, w_cut=wc, n_matsubara=n).tc()
            assert t > 0
            f = t / (C_INF * np.sqrt(sp.lambda_total()) * w2 / KB_MEV_PER_K)
            assert f < 1.0, f"𝔣>1 报警：{f:.4f}"
            L.append(np.log(sp.lambda_total())); S.append(np.log(wl / w2)); Y.append(np.log(f))
    return np.array(L), np.array(S), np.array(Y)


def ladder(title, keys, mu, with_S, wcut_abs=None, ref=None):
    print()
    print("=" * 92)
    print(title)
    print("=" * 92)
    cols = ("a", "b", "c", "q") if with_S else ("a", "b", "q")
    print("  " + f"{'N':>6}" + "".join(f"{c:>13}" for c in cols) + f"{'rms':>12}{'秒':>8}")
    out = {}
    for n in NS:
        t0 = time.time()
        L, S, Y = run(keys, mu, n, wcut_abs)
        c, rms = fit(L, S, Y, with_S)
        out[n] = c
        print("  " + f"{n:>6}" + "".join(f"{v:>13.6f}" for v in c)
              + f"{rms:>12.3e}{time.time()-t0:>8.1f}")
    q4, q16 = out[400][-1], out[1600][-1]
    print(f"  ⟹ q(400) = {q4:.6f} ，q(1600) = {q16:.6f} ，"
          f"Δ = {q16-q4:+.2e}（相对 {abs((q16-q4)/q16):.3%}）")
    print(f"     判决：{'CONVERGED（N=400 可作生产用）' if abs((q16-q4)/q16) < 1e-2 else '★ 未收敛，N=400 不可用'}")
    if ref is not None:
        d = max(abs(out[n][-1] - ref[n]) for n in NS)
        print(f"     对照 peer 的五个数：最劣逐点差 {d:.2e} "
              f"⟹ {'复现 PASS' if d < 5e-5 else '★ 不复现，须查'}")
    return out


PEER = {100: -0.328672, 200: -0.329274, 400: -0.329425, 800: -0.329463, 1600: -0.329473}
ladder("§1  复现 peer：H3S-anh 单谱，μ*=0，w_cut=500 meV（无 logS 项，单谱无 S 变化）",
       ["H3S-anh"], 0.0, False, wcut_abs=500.0, ref=PEER)
ladder("§2  ★ peer 未测的那一档：H3S-anh 单谱，**μ*=0.10**，w_cut=10·ω̄₂",
       ["H3S-anh"], 0.10, False)
ladder("§3  ★ 生产量本身：四谱合并（含 logS 项，36 点），**μ*=0.10**，w_cut=10·ω̄₂",
       list(RAW), 0.10, True)
print()
print("=" * 92)
print("总判决")
print("=" * 92)
print("  · peer 的 μ*=0 阶梯若复现 ⟹ 它那一档的结论成立，但**它不是我承重值所在的那一档**。")
print("  · 承重的是 §3（四谱合并、μ*=0.10、含 logS）——只有它直接管 `q=−0.5565` 这个数。")
print("  · 方向：N 截断在此处**不是**单边界估计（Tc 随 N 双向漂移）⟹ 收敛性判决靠 Δq 的量级，")
print("    不靠不等号方向；这与 `nphi` 那次（N 截断本征值 fail open）是**不同**的情形。")
print("=" * 92)
