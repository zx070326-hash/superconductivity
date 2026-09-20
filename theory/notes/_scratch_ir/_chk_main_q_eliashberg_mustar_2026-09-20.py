# -*- coding: utf-8 -*-
"""接 `_chk_main_q_eliashberg_independent_2026-09-20.py`：μ*>0 轨道。

问题：μ*=0 的理论侧给 q = −0.3219 ± 0.0182，而 AD75 29 行数据给 q = −0.525。
      **最大的混杂就是 μ***——AD75 Table I 逐行印出的 μ* 在 0.092–0.111 之间
      （页图文本层可读；这里只用「≈0.10 且跨行近乎均匀」这一弱事实，不逐行引用）。
      μ* 对小 λ 的压制相对更强 ⟹ 可能同时抬高 b 与加深 q。本脚本量它。

⚠ μ*>0 时 Tc 对 ω_c **非单调**（仓库 h3s_mustar_cutoff_convention：谷底 ≈550 meV）
  ⟹ 报数必带 w_cut，且必须跑**两套约定**看 q 的约定敏感度。
  约定 I：w_cut = 10·ω̄₂（谱相对；注意 ω̄₂ 在振幅重标度下精确不变 ⟹ 整条 λ 扫描里 w_cut 恒定）
  约定 II：w_cut = 500 meV（绝对，与仓库 cert 链一致）

严格性：numerical-evidence（float64 Eliashberg + 有限 λ 网格）。
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
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K   # noqa: E402

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


def scan(mu, wcut_mode):
    """返回 (q_per_spectrum, pooled_coef, pooled_err, fmax, n_fail)。"""
    rows, nfail = [], 0
    for k, (w, a) in RAW.items():
        sp_ref = Spectrum.from_table(w, a)
        lam_ref, w2, wl = sp_ref.lambda_total(), sp_ref.omega_2(), sp_ref.omega_log()
        wcut = 10.0 * w2 if wcut_mode == "I" else 500.0
        for lam in LAM_GRID:
            sp = Spectrum.from_table(w, a * (lam / lam_ref))
            try:
                tc = Solver(spectrum=sp, mu_star=mu, w_cut=wcut,
                            n_matsubara=NMATS).tc()
            except Exception as e:                       # 不吞异常，计数并报出
                nfail += 1
                print(f"    [FAIL] {k} λ={lam:.3f} μ*={mu}: {type(e).__name__}: {e}")
                continue
            if not (tc > 0):
                nfail += 1
                continue
            ceil_K = C_INF * np.sqrt(sp.lambda_total()) * w2 / KB_MEV_PER_K
            rows.append((k, wl / w2, sp.lambda_total(), tc, tc / ceil_K, wcut))
    if len(rows) < 12:
        return None, None, None, None, nfail
    arr_k = np.array([r[0] for r in rows])
    S = np.array([r[1] for r in rows]); lam = np.array([r[2] for r in rows])
    f = np.array([r[4] for r in rows])
    x, y = np.log(lam), np.log(f)
    qs = {}
    for k in FILES:
        m = arr_k == k
        if m.sum() < 4:
            continue
        A = np.vstack([np.ones(m.sum()), x[m], x[m] ** 2]).T
        qs[k] = np.linalg.lstsq(A, y[m], rcond=None)[0][2]
    A = np.vstack([np.ones(len(x)), x, np.log(S), x ** 2]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ coef
    dof = len(x) - 4
    cov = (res @ res / dof) * np.linalg.inv(A.T @ A)
    return qs, coef, np.sqrt(np.diag(cov)), f.max(), nfail


print("=" * 86)
print("μ* 轨道：q 与 b 随 μ* 的走向（λ 网格 = AD75 支撑域 [0.69, 2.59] 九点，N=400）")
print("AD75 Table I 印出的 μ* ∈ [0.092, 0.111]（近乎均匀 ≈0.10）")
print("对照（本轮 μ*=0 已测）：b = +0.7389 ± 0.0125，q = −0.3219 ± 0.0182，c = +0.6482")
print("对照（AD75 29 行经验回归）：b = +0.9972 ± 0.041，q = −0.525，c = +0.7198")
print()
for mode, tag in (("I", "w_cut = 10·ω̄₂（谱相对）"), ("II", "w_cut = 500 meV（绝对）")):
    print("-" * 86)
    print(f"约定 {mode}：{tag}")
    print(f"{'μ*':>6} {'b':>9} {'±':>7} {'c':>9} {'q':>9} {'±':>7}   {'q 四谱范围':>20} {'𝔣max':>7} {'失败':>5}")
    for mu in (0.0, 0.05, 0.10, 0.13):
        qs, coef, err, fmax, nf = scan(mu, mode)
        if coef is None:
            print(f"{mu:6.3f}   —— 可用点不足（失败 {nf}）")
            continue
        qv = np.array(list(qs.values()))
        rng = f"[{qv.min():+.4f}, {qv.max():+.4f}]"
        alarm = " <<< 𝔣>1" if fmax > 1 else ""
        print(f"{mu:6.3f} {coef[1]:+9.4f} {err[1]:7.4f} {coef[2]:+9.4f} "
              f"{coef[3]:+9.4f} {err[3]:7.4f}   {rng:>20} {fmax:7.4f} {nf:5d}{alarm}")
    print()
print("严格性：全部为 numerical-evidence。")
