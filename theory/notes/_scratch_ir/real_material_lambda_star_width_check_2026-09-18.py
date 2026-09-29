"""
数值证据（numerical-evidence 级），非区间算术。

目的：把 lambda_stability_tc_joint_ceiling_2026-09-18.md 里在"对称两点谱玩具模型"上得到的
机制 kappa_S>0（谱展宽 sigma^2 提升 Semenok 稳定性阈值 lambda_*[P]）拿到**真实材料**的
第一性原理 alpha^2F(omega) 谱上独立检验方向是否一致。

数据来源：theory/data/external/cerqueira_raw/extracted/batch-q/<material>/a2F.dos1
（Cerqueira et al. arXiv:2307.10728, Materials Cloud 3kbt5-r3n56，Espresso 计算，
batch-q 250 个材料，已在仓库内，未额外下载）。col1=omega[Ry]，col2=alpha^2F_total。

方法：
  P(domega) := 2 alpha^2F(omega) / (lambda omega) domega   (与 P3_tc_max_pressure.md /
    KAY 框架同一个测度，见 lambda_stability_tc_joint_ceiling_2026-09-18.md §3 "记号同构")
  lambda   = 2 int alpha^2F/omega domega                    (直接从 a2F.dos1 数值积分)
  w2av     = sqrt[(2/lambda) int alpha^2F * omega domega]  = sqrt(<omega^2>_P)
  wlog     = exp[(2/lambda) int (alpha^2F/omega) ln(omega) domega] = exp(<ln omega>_P)
  r        := w2av/wlog   (>=1 恒成立by Jensen; r=1 当且仅当 P 是 delta 谱即"Einstein 谱"；
                            这是本仓库对"谱展宽 sigma"的**唯一**可从两个已发表矩直接构造的
                            无量纲代理——r 不是 sigma 本身，是 sigma 的单调代理，仅在 sigma 不
                            太大时可与 toy model 的 sigma^2 参数近似对应；见笔记 §2 的严格性声明)
  lambda_*[P] = 1 / max_{u=2 pi T>0} int g(omega/u) P(domega)
     g(x) = 6x^2+12x^3 Im[psi'(ix)] + 6x^4 Re[psi''(ix)]  (Semenok, 与 P3 逐字核对过的核)

g(x) 用 mpmath 建一次查表（对数网格，dps=25），再用三次样条插值供数值积分/一维最大化重复调用，
避免 250 个材料 x 一维搜索 x 复化积分点数的组合爆炸式 mpmath 直接调用。

输出：对每个材料打印 lambda(重算) vs lambda(文件自报) 的核验、r、lambda_*，
最后给出 r vs lambda_* 的 Pearson/Spearman 相关系数（跨材料），并与 r vs lambda(实际DFT值)
的相关系数对比。
"""
import glob, os, re, sys
import numpy as np
import mpmath as mp
from scipy.interpolate import CubicSpline
from scipy.stats import pearsonr, spearmanr

mp.mp.dps = 25

# ---------- build g(x) lookup table once ----------
def g_mp(x):
    x = mp.mpf(x)
    ix = mp.mpc(0, 1) * x
    psi1 = mp.polygamma(1, ix)
    psi2 = mp.polygamma(2, ix)
    return float(6*x**2 + 12*x**3*mp.im(psi1) + 6*x**4*mp.re(psi2))

print("building g(x) lookup table...", file=sys.stderr)
xs_log = np.linspace(np.log(1e-5), np.log(200.0), 1400)
xs = np.exp(xs_log)
gs = np.array([g_mp(x) for x in xs])
g_spline = CubicSpline(xs_log, gs, extrapolate=True)

def g_fast(x):
    x = np.asarray(x, dtype=float)
    x = np.clip(x, 1e-5, 200.0)
    return g_spline(np.log(x))

print("table built.", file=sys.stderr)

# ---------- parse one material's a2F.dos1 ----------
def parse_a2f(path):
    om, a2f = [], []
    lam_reported = None
    with open(path, errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            m = re.match(r"^\s*lambda\s*=\s*([\-0-9.Ee+]+)", s)
            if m:
                lam_reported = float(m.group(1))
                continue
            parts = s.split()
            try:
                w = float(parts[0]); a = float(parts[1])
            except (ValueError, IndexError):
                continue
            om.append(w); a2f.append(a)
    om = np.array(om); a2f = np.array(a2f)
    return om, a2f, lam_reported

def moments_and_lambda_star(om, a2f, u_lo_factor=0.02, u_hi_factor=50.0, n_u=400):
    mask = om > 1e-8
    om, a2f = om[mask], a2f[mask]
    if len(om) < 10 or np.all(a2f <= 0):
        return None
    lam = 2.0 * np.trapezoid(a2f/om, om)
    if lam <= 1e-6:
        return None
    w2av = np.sqrt((2.0/lam) * np.trapezoid(a2f*om, om))
    wlog = np.exp((2.0/lam) * np.trapezoid((a2f/om)*np.log(om), om))
    r = w2av/wlog
    # lambda_*[P]: P(domega)=2 a2f/(lam*omega) domega ; integrand for fixed u:
    #   F(u) = (2/lam) * int (a2f/omega) * g(omega/u) domega
    wmin, wmax = om.min(), om.max()
    u_grid = np.exp(np.linspace(np.log(wmin*u_lo_factor), np.log(wmax*u_hi_factor), n_u))
    Fvals = np.empty(n_u)
    for i, u in enumerate(u_grid):
        Fvals[i] = (2.0/lam) * np.trapezoid((a2f/om) * g_fast(om/u), om)
    imax = np.argmax(Fvals)
    Fmax = Fvals[imax]
    lam_star = 1.0/Fmax if Fmax > 0 else np.nan
    return dict(lam=lam, w2av=w2av, wlog=wlog, r=r, lam_star=lam_star, u_star=u_grid[imax])

# ---------- loop over batch-q materials ----------
base = "D:/physics/superconductivity/theory/data/external/cerqueira_raw/extracted/batch-q"
dirs = sorted(glob.glob(os.path.join(base, "*")))
rows = []
for d in dirs:
    f1 = os.path.join(d, "a2F.dos1")
    if not os.path.exists(f1):
        continue
    name = os.path.basename(d)
    om, a2f, lam_rep = parse_a2f(f1)
    res = moments_and_lambda_star(om, a2f)
    if res is None:
        continue
    lam_err = None
    if lam_rep is not None and lam_rep > 1e-6:
        lam_err = abs(res["lam"]-lam_rep)/lam_rep
    rows.append((name, res["lam"], lam_rep, lam_err, res["w2av"], res["wlog"], res["r"], res["lam_star"]))

print(f"parsed {len(rows)} materials with usable a2F.dos1", file=sys.stderr)

# sanity: lambda recompute vs file-reported lambda
errs = [row[3] for row in rows if row[3] is not None]
print(f"\n=== lambda recompute vs a2F.dos1 self-reported 'lambda=' line ===")
print(f"n={len(errs)}  median rel err = {np.median(errs):.4%}  max rel err = {np.max(errs):.4%}")

# filter to a physically sane band: exclude lambda too small (noisy wlog/lambda_star) or too large (bad DOS)
lam_arr = np.array([row[1] for row in rows])
r_arr = np.array([row[6] for row in rows])
lamstar_arr = np.array([row[7] for row in rows])
good = (lam_arr > 0.15) & (lam_arr < 3.0) & np.isfinite(lamstar_arr) & (r_arr > 0.98) & (r_arr < 3.0)
print(f"\nn_total={len(rows)}  n_good(0.15<lambda<3, finite lambda_*, 0.98<r<3)={good.sum()}")

r_g = r_arr[good]; lamstar_g = lamstar_arr[good]; lam_g = lam_arr[good]
headroom_g = lamstar_g/lam_g  # stability margin

pear_r_lamstar = pearsonr(r_g, lamstar_g)
spear_r_lamstar = spearmanr(r_g, lamstar_g)
pear_r_headroom = pearsonr(r_g, headroom_g)
spear_r_headroom = spearmanr(r_g, headroom_g)
pear_r_lam = pearsonr(r_g, lam_g)

print(f"\n=== correlation: spectral width proxy r=w2av/wlog  vs  Semenok stability threshold lambda_*[P] ===")
print(f"Pearson(r, lambda_*)  = {pear_r_lamstar.statistic:.4f}  (p={pear_r_lamstar.pvalue:.2e})")
print(f"Spearman(r, lambda_*) = {spear_r_lamstar.statistic:.4f}  (p={spear_r_lamstar.pvalue:.2e})")
print(f"\n=== correlation: r  vs  headroom lambda_*[P]/lambda_actual (how far from instability) ===")
print(f"Pearson(r, headroom)  = {pear_r_headroom.statistic:.4f}  (p={pear_r_headroom.pvalue:.2e})")
print(f"Spearman(r, headroom) = {spear_r_headroom.statistic:.4f}  (p={spear_r_headroom.pvalue:.2e})")
print(f"\n(control) Pearson(r, lambda_actual) = {pear_r_lam.statistic:.4f}  (p={pear_r_lam.pvalue:.2e})")

# also print a handful of rows for manual spot-check
print("\n=== sample rows (first 15 of the 'good' subset) ===")
idxs = np.where(good)[0][:15]
for i in idxs:
    name, lam, lam_rep, lam_err, w2av, wlog, r, lam_star = rows[i]
    print(f"{name:35s} lambda={lam:.4f} (file {lam_rep}) r={r:.4f}  lambda_*={lam_star:.4f}  headroom={lam_star/lam:.3f}")

# bin by r into low/high width halves, compare mean lambda_* (a coarse, robust alternative to correlation)
order = np.argsort(r_g)
half = len(order)//2
lo_idx, hi_idx = order[:half], order[half:]
print(f"\n=== coarse split: low-r half (n={half}) vs high-r half (n={len(order)-half}) ===")
print(f"mean r      : low={r_g[lo_idx].mean():.4f}  high={r_g[hi_idx].mean():.4f}")
print(f"mean lambda*: low={lamstar_g[lo_idx].mean():.4f}  high={lamstar_g[hi_idx].mean():.4f}")
print(f"mean headrm : low={headroom_g[lo_idx].mean():.4f}  high={headroom_g[hi_idx].mean():.4f}")
print(f"mean lambda : low={lam_g[lo_idx].mean():.4f}  high={lam_g[hi_idx].mean():.4f}")
