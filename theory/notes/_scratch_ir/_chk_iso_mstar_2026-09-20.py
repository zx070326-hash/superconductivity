# -*- coding: utf-8 -*-
"""
alpha(M)=0 的临界质量 M*：预算 + 内插散布 + mu* 敏感度 + 残差推动。

背景：alpha := -dlnTc/dlnM = alpha_omega + alpha_lambda（精确可加分解，见
_chk_iso_alpha_decomp_2026-09-20.py）。SSCHA 行 alpha(H->D)=-0.379<0（真.反常），
alpha(D->T)=-0.148 仍 <0 但该值本身病态（领头项相消至 1.4%）。
本脚本检验：把"某一步的 alpha 值"换成"alpha(M)=0 的临界质量 M*"是否真的更良态、
是否可用（能否由三个同位素点内插稳定给出）。

数据：Errea et al. arXiv:1305.7123 Table I, SSCHA 列。
"""
import math

M = {"H": 1.007825, "D": 2.014102, "T": 3.016049}
SSCHA = {"H": (0.40, 405.0, 5.0), "D": (0.46, 304.0, 6.5), "T": (0.48, 257.0, 6.9)}
MUSTAR0 = 0.085


def tc_ad(lam, wlog, mus):
    den = lam - mus * (1.0 + 0.62 * lam)
    if den <= 1e-9:
        return float("nan")
    return (wlog / 1.2) * math.exp(-1.04 * (1.0 + lam) / den)


def dlnF_dlnlam(lam, mus):
    """d ln F / d ln lambda for F = Tc_AD/omega_log, closed form (analytic)."""
    den = lam - mus * (1.0 + 0.62 * lam)
    dden_dlam = 1.0 - 0.62 * mus
    # d(lnF)/dlam = -1.04 * [den - (1+lam)*dden_dlam] / den**2
    dlnF_dlam = -1.04 * (den - (1.0 + lam) * dden_dlam) / den ** 2
    return dlnF_dlam * lam  # convert d/dlam -> d/dlnlam = lam * d/dlam


# ---------------------------------------------------------------------------
# Step 0: budget. Compare condition number of M* vs alpha(D->T) under a
# perturbation of mu* (same perturbation, same units), and under rounding of
# Tc_ME digits (same perturbation used in the alpha-decomp script).
# ---------------------------------------------------------------------------

def alpha_of_M_loglog_quadratic(mus=MUSTAR0):
    """
    Build lambda(M), omega_log(M) as the UNIQUE quadratic in ln M interpolating
    the 3 SSCHA points exactly (Lagrange quadratic in x=ln M). Returns a
    callable alpha(M) = alpha_omega(M) + alpha_lambda(M), using the closed-form
    local log-derivatives of the quadratic (exact, no finite differences).
    """
    xs = [math.log(M[k]) for k in ("H", "D", "T")]
    lams = [SSCHA[k][0] for k in ("H", "D", "T")]
    wlogs = [SSCHA[k][1] for k in ("H", "D", "T")]
    lnlams = [math.log(v) for v in lams]
    lnws = [math.log(v) for v in wlogs]

    def lagrange_quad_and_deriv(xs, ys, x):
        # y(x) = sum_i y_i * L_i(x); dy/dx = sum_i y_i * L_i'(x)
        n = len(xs)
        y = 0.0
        dy = 0.0
        for i in range(n):
            Li = 1.0
            dLi = 0.0
            for j in range(n):
                if j == i:
                    continue
                Li *= (x - xs[j]) / (xs[i] - xs[j])
            # derivative of product via sum-of-logs trick (safe here, n=3, no zero factors at generic x)
            for j in range(n):
                if j == i:
                    continue
                term = 1.0 / (xs[i] - xs[j])
                for k in range(n):
                    if k == i or k == j:
                        continue
                    term *= (x - xs[k]) / (xs[i] - xs[k])
                dLi += term
            y += ys[i] * Li
            dy += ys[i] * dLi
        return y, dy

    def alpha_of_M(Mval):
        x = math.log(Mval)
        lnlam, dlnlam_dlnM = lagrange_quad_and_deriv(xs, lnlams, x)
        lnw, dlnw_dlnM = lagrange_quad_and_deriv(xs, lnws, x)
        lam = math.exp(lnlam)
        a_omega = -dlnw_dlnM
        a_lambda = -dlnF_dlnlam(lam, mus) * dlnlam_dlnM
        return a_omega + a_lambda, a_omega, a_lambda, lam

    return alpha_of_M


def alpha_of_M_piecewise_powerlaw(mus=MUSTAR0):
    """
    Piecewise power law: on [M_H,M_D] use the exact (H,D) log-log slope for
    BOTH lambda and omega_log (so the function matches endpoints exactly on
    each sub-interval); on [M_D,M_T] use the (D,T) slope. Discontinuous
    derivative at M_D, continuous value.
    """
    keys = ("H", "D", "T")
    segs = []
    for a, b in (("H", "D"), ("D", "T")):
        Ma, Mb = M[a], M[b]
        lam_a, w_a = SSCHA[a][0], SSCHA[a][1]
        lam_b, w_b = SSCHA[b][0], SSCHA[b][1]
        p_lam = math.log(lam_b / lam_a) / math.log(Mb / Ma)
        p_w = math.log(w_b / w_a) / math.log(Mb / Ma)
        segs.append((Ma, Mb, lam_a, w_a, p_lam, p_w))

    def alpha_of_M(Mval):
        for (Ma, Mb, lam_a, w_a, p_lam, p_w) in segs:
            if Ma - 1e-12 <= Mval <= Mb + 1e-12:
                lam = lam_a * (Mval / Ma) ** p_lam
                a_omega = -p_w
                a_lambda = -dlnF_dlnlam(lam, mus) * p_lam
                return a_omega + a_lambda, a_omega, a_lambda, lam
        # extrapolate using nearest segment
        if Mval < segs[0][0]:
            Ma, Mb, lam_a, w_a, p_lam, p_w = segs[0]
        else:
            Ma, Mb, lam_a, w_a, p_lam, p_w = segs[1]
        lam = lam_a * (Mval / Ma) ** p_lam
        a_omega = -p_w
        a_lambda = -dlnF_dlnlam(lam, mus) * p_lam
        return a_omega + a_lambda, a_omega, a_lambda, lam

    return alpha_of_M


def alpha_of_M_global_powerlaw_lsq(mus=MUSTAR0):
    """
    Global power law lambda=A*M^p, omega=B*M^q fit by least squares (log-log
    linear regression) over all 3 points. Constant exponents p,q (no
    curvature) => alpha_omega, alpha_lambda are CONSTANT in M (alpha(M) either
    has no root or is identically the same sign -- included to show what an
    over-simplified interpolation gives).
    """
    xs = [math.log(M[k]) for k in ("H", "D", "T")]
    lnlams = [math.log(SSCHA[k][0]) for k in ("H", "D", "T")]
    lnws = [math.log(SSCHA[k][1]) for k in ("H", "D", "T")]

    def lsq_slope(xs, ys):
        n = len(xs)
        xbar = sum(xs) / n
        ybar = sum(ys) / n
        num = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
        den = sum((x - xbar) ** 2 for x in xs)
        slope = num / den
        intercept = ybar - slope * xbar
        return slope, intercept

    p_lam, b_lam = lsq_slope(xs, lnlams)
    p_w, b_w = lsq_slope(xs, lnws)

    def alpha_of_M(Mval):
        x = math.log(Mval)
        lnlam = b_lam + p_lam * x
        lam = math.exp(lnlam)
        a_omega = -p_w
        a_lambda = -dlnF_dlnlam(lam, mus) * p_lam
        return a_omega + a_lambda, a_omega, a_lambda, lam

    return alpha_of_M


def find_root(alpha_of_M, lo, hi, tol=1e-8):
    flo = alpha_of_M(lo)[0]
    fhi = alpha_of_M(hi)[0]
    if flo == 0:
        return lo
    if fhi == 0:
        return hi
    if flo * fhi > 0:
        return None  # no sign change in [lo,hi]
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fmid = alpha_of_M(mid)[0]
        if fmid == 0 or (hi - lo) < tol:
            return mid
        if flo * fmid < 0:
            hi = mid
            fhi = fmid
        else:
            lo = mid
            flo = fmid
    return 0.5 * (lo + hi)


print("=" * 78)
print("STEP 0: 预算 -- M* 在 H..T 区间内是否存在，条件数是否真的更好？")
print("=" * 78)

builders = {
    "quadratic-loglog (exact interp, curved)": alpha_of_M_loglog_quadratic,
    "piecewise-powerlaw (exact interp, kinked)": alpha_of_M_piecewise_powerlaw,
    "global-powerlaw-lsq (constant exponent)": alpha_of_M_global_powerlaw_lsq,
}

roots = {}
for name, builder in builders.items():
    af = builder(MUSTAR0)
    aH, _, _, _ = af(M["H"])
    aD, _, _, _ = af(M["D"])
    aT, _, _, _ = af(M["T"])
    print(f"\n[{name}]")
    print(f"  alpha(H)={aH:+.4f}  alpha(D)={aD:+.4f}  alpha(T)={aT:+.4f}")
    r_HD = find_root(af, M["H"], M["D"])
    r_DT = find_root(af, M["D"], M["T"])
    print(f"  root in [H,D]: {r_HD}")
    print(f"  root in [D,T]: {r_DT}")
    roots[name] = (r_HD, r_DT)

print()
print("-> M_D=%.6f  M_T=%.6f  (spacing M_T-M_D=%.6f u)" % (M["D"], M["T"], M["T"] - M["D"]))

print()
print("=" * 78)
print("STEP 0b: 条件数比较 -- 扰动 mu* 由 0.085 -> 0.095 (Δmu*=0.010, 仓库常见误差棒)")
print("=" * 78)
af0 = alpha_of_M_loglog_quadratic(0.085)
af1 = alpha_of_M_loglog_quadratic(0.095)
aDT0 = af0(M["D"])[0] - 0  # value form: use direct discrete alpha(D->T) formula instead (below)


def alpha_DT_direct(mus):
    lamD, wD, _ = SSCHA["D"]
    lamT, wT, _ = SSCHA["T"]
    TcD = tc_ad(lamD, wD, mus)
    TcT = tc_ad(lamT, wT, mus)
    return -math.log(TcT / TcD) / math.log(M["T"] / M["D"])


for name, builder in builders.items():
    r0 = find_root(builder(0.085), M["D"] - 0.3, M["T"] + 0.5)
    r1 = find_root(builder(0.095), M["D"] - 0.3, M["T"] + 0.5)
    dr = None if (r0 is None or r1 is None) else (r1 - r0)
    drtxt = "n/a (no root)" if dr is None else f"{dr/0.010:.3f} u per unit mu*"
    print(f"[{name}] M*(mu*=0.085)={r0}  M*(mu*=0.095)={r1}  dM*/dmu*~{drtxt}")

a_dt0 = alpha_DT_direct(0.085)
a_dt1 = alpha_DT_direct(0.095)
print(f"\nalpha(D->T) direct-AD: mu*=0.085 -> {a_dt0:+.4f}; mu*=0.095 -> {a_dt1:+.4f}; "
      f"dalpha/dmu* ~ {(a_dt1-a_dt0)/0.010:+.3f}")

print()
print("=" * 78)
print("STEP 2: 内插散布 -- 三种内插方案给出的 M* 差多少？")
print("=" * 78)
for name, (r_HD, r_DT) in roots.items():
    print(f"[{name}]  root(H,D)={r_HD}  root(D,T)={r_DT}")

print()
print("=" * 78)
print("STEP 3: mu* 扫描 M*(mu*) in [0.06, 0.15] (quadratic-loglog 内插，主线)")
print("=" * 78)
mus_grid = [0.06 + 0.005 * i for i in range(19)]
mstars = []
for mus in mus_grid:
    af = alpha_of_M_loglog_quadratic(mus)
    r = find_root(af, M["H"], M["T"] + 1.0)
    mstars.append(r)
    print(f"  mu*={mus:.3f}  M*={r}")

valid = [(mu, r) for mu, r in zip(mus_grid, mstars) if r is not None]
if len(valid) >= 2:
    lo_mu, lo_r = valid[0]
    hi_mu, hi_r = valid[-1]
    print(f"\n  M* range over mu* in [{lo_mu},{hi_mu}]: [{min(r for _,r in valid):.4f}, {max(r for _,r in valid):.4f}] u")

print()
print("=" * 78)
print("STEP 5: 谱形状残差 (>=0.0638, alpha_ME - alpha_AD 同号负) 的推动量")
print("=" * 78)
# Apply as a constant systematic offset to alpha(M) (shifts alpha_omega down by 0.0638)
# and re-find M* under quadratic-loglog interpolation.
RESID = 0.0638


def alpha_of_M_with_offset(mus, offset):
    af = alpha_of_M_loglog_quadratic(mus)

    def shifted(Mval):
        a, aw, al, lam = af(Mval)
        return a - offset, aw - offset, al, lam

    return shifted


af_shift = alpha_of_M_with_offset(0.085, RESID)
r_shift = find_root(af_shift, M["H"], M["T"] + 1.0)
r_base = find_root(alpha_of_M_loglog_quadratic(0.085), M["H"], M["T"] + 1.0)
print(f"  M*(no offset)   = {r_base}")
print(f"  M*(offset -0.0638, worst-sign) = {r_shift}")
if r_base is not None and r_shift is not None:
    print(f"  shift = {r_shift - r_base:+.4f} u")
