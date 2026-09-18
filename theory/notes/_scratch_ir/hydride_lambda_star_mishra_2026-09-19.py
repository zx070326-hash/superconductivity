"""
numerical-evidence 级，非区间算术。

把 lambda_lattice_width_real_material_check_2026-09-18.md 的方法（真实 alpha^2F(omega) 谱上
直接数值求解 Semenok lambda_*[P]，与 r=w2av/wlog 关联，评估 headroom=lambda_*/lambda_actual）
应用到强耦合氢化物：Mishra et al. Figshare数据集 H3S/D3S 在 200GPa/160GPa、谐波/非谐 EPW alpha^2F 谱。

数据来源（未额外下载，主仓库已有）：
  D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/
    Data-H3S-D3S-PH-EPW-Anharmonic/{H3S,D3S}-{200,160}GPa/EPW_{harmonic,anharmonic,anh}/{H3S,D3S}.a2f
  列：w[meV], a2f, integrated 2*a2f/w (累积lambda(omega))
  文件尾部另有 "Summed el-ph coupling" 自报总 lambda，可用于核验积分流程。

方法与 real_material_lambda_star_width_check_2026-09-18.py 完全一致（同一个 g(x) 核，
同一个 lambda_*[P] = 1/max_u int g(w/u) P(dw) 定义），仅改数据源与单位（meV 而非 Ry，
lambda_*[P] 是无量纲量，g(x) 的宗量 x=omega/u 与能量单位无关，只需 w 和 u 用同一单位）。
"""
import re
import numpy as np
import mpmath as mp
from scipy.interpolate import CubicSpline

mp.mp.dps = 25

def g_mp(x):
    x = mp.mpf(x)
    ix = mp.mpc(0, 1) * x
    psi1 = mp.polygamma(1, ix)
    psi2 = mp.polygamma(2, ix)
    return float(6*x**2 + 12*x**3*mp.im(psi1) + 6*x**4*mp.re(psi2))

print("building g(x) lookup table...")
xs_log = np.linspace(np.log(1e-5), np.log(200.0), 1400)
xs = np.exp(xs_log)
gs = np.array([g_mp(x) for x in xs])
g_spline = CubicSpline(xs_log, gs, extrapolate=True)

def g_fast(x):
    x = np.asarray(x, dtype=float)
    x = np.clip(x, 1e-5, 200.0)
    return g_spline(np.log(x))

def parse_epw_a2f(path):
    om, a2f = [], []
    lam_reported = None
    with open(path, errors="replace") as f:
        lines = f.readlines()
    for line in lines[1:]:
        s = line.strip()
        m = re.match(r"^Summed el-ph coupling\s+([\-0-9.Ee+]+)", s)
        if m:
            lam_reported = float(m.group(1))
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            w = float(parts[0]); a = float(parts[1])
        except ValueError:
            continue
        om.append(w); a2f.append(max(a, 0.0))
    return np.array(om), np.array(a2f), lam_reported

def moments_and_lambda_star(om, a2f, u_lo_factor=0.02, u_hi_factor=50.0, n_u=500):
    mask = om > 1e-8
    om, a2f = om[mask], a2f[mask]
    lam = 2.0 * np.trapezoid(a2f/om, om)
    w2av = np.sqrt((2.0/lam) * np.trapezoid(a2f*om, om))
    wlog = np.exp((2.0/lam) * np.trapezoid((a2f/om)*np.log(om), om))
    r = w2av/wlog
    wmin, wmax = om.min(), om.max()
    u_grid = np.exp(np.linspace(np.log(wmin*u_lo_factor), np.log(wmax*u_hi_factor), n_u))
    Fvals = np.empty(n_u)
    for i, u in enumerate(u_grid):
        Fvals[i] = (2.0/lam) * np.trapezoid((a2f/om) * g_fast(om/u), om)
    imax = np.argmax(Fvals)
    Fmax = Fvals[imax]
    lam_star = 1.0/Fmax if Fmax > 0 else np.nan
    # refine near the max with a finer local grid (edge/coarse-grid check)
    if 0 < imax < n_u-1:
        u_fine = np.exp(np.linspace(np.log(u_grid[imax-1]), np.log(u_grid[imax+1]), 200))
        Ff = np.array([(2.0/lam)*np.trapezoid((a2f/om)*g_fast(om/u), om) for u in u_fine])
        imf = np.argmax(Ff)
        Fmax2 = Ff[imf]
        if Fmax2 > Fmax:
            Fmax = Fmax2
            lam_star = 1.0/Fmax
    return dict(lam=lam, w2av=w2av, wlog=wlog, r=r, lam_star=lam_star, u_star=u_grid[imax],
                edge_hit=(imax == 0 or imax == n_u-1))

BASE = "D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic"
FILES = {
    "H3S-200GPa-har": f"{BASE}/H3S-200GPa/EPW_harmonic/H3S.a2f",
    "H3S-200GPa-anh": f"{BASE}/H3S-200GPa/EPW_anharmonic/H3S.a2f",
    "D3S-200GPa-har": f"{BASE}/D3S-200GPa/EPW_harmonic/D3S.a2f",
    "D3S-200GPa-anh": f"{BASE}/D3S-200GPa/EPW_anharmonic/D3S.a2f",
    "H3S-160GPa-anh": f"{BASE}/H3S-160GPa/EPW_anharmonic/H3S.a2f",
    "D3S-160GPa-anh": f"{BASE}/D3S-160GPa/EPW_anh/D3S.a2f",
}

print(f"\n{'case':18s} {'lambda(recomp)':>14} {'lambda(file)':>13} {'rel_err':>9} {'w2av(meV)':>10} {'wlog(meV)':>10} {'r':>7} {'lambda_*':>9} {'headroom':>9} {'edge?':>6}")
results = {}
for name, path in FILES.items():
    om, a2f, lam_rep = parse_epw_a2f(path)
    res = moments_and_lambda_star(om, a2f)
    results[name] = res
    err = abs(res['lam']-lam_rep)/lam_rep if lam_rep else float('nan')
    hr = res['lam_star']/res['lam']
    print(f"{name:18s} {res['lam']:14.4f} {lam_rep:13.4f} {err:9.4%} {res['w2av']:10.2f} {res['wlog']:10.2f} {res['r']:7.4f} {res['lam_star']:9.4f} {hr:9.4f} {str(res['edge_hit']):>6}")

print("\n=== headroom summary ===")
for name, res in results.items():
    print(f"{name:18s} headroom = lambda_*/lambda = {res['lam_star']/res['lam']:.4f}")

print("\n=== r vs headroom across the 6 spectra (direction check, n too small for correlation stats) ===")
for name, res in results.items():
    print(f"{name:18s} r={res['r']:.4f}  headroom={res['lam_star']/res['lam']:.4f}")
