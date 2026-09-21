"""
Route R' M1: freed-width trial family (tau_n = u_n a^2/(a^2+(n+1/2)^2), a=alpha*w)
+ exact reduced-resolvent R (Prop.1(ii)) instead of the Davis-Kahan relaxation
bound1 = 2x^3 E^2/g, for tightening Lemma C (H(x)=h(sqrt x) concave) on the
compact interval.

Does NOT modify any existing lemmaC_*.py baseline script; reuses O_of/top/
THETA_STAR from lemmaC_compact_probe.py (same convention as the rest of the
repo's exploration scripts).

Levels: everything here is float64 numerical-evidence unless stated otherwise.
No claim in this file is 'analytic' or 'rational-exact'.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from lemmaC_compact_probe import O_of, top, THETA_STAR

N_DEFAULT = 600


# ---------------------------------------------------------------------------
# Section 0: shared machinery
# ---------------------------------------------------------------------------

def build(w, N=N_DEFAULT):
    x = w * w
    f = lambda k: x / (k * k + x)
    phi = lambda k: k * k / (k * k + x) ** 3
    psi = lambda k: 1.0 / (k * k + x) ** 2
    A = O_of(f, N)
    Ophi = O_of(phi, N)
    Opsi = O_of(psi, N)
    h, v, l2, ev, EV = top(A)
    g = h - l2
    nphi = np.abs(np.linalg.eigvalsh(Ophi)).max()
    npsi = np.abs(np.linalg.eigvalsh(Opsi)).max()
    return dict(w=w, x=x, N=N, A=A, Ophi=Ophi, Opsi=Opsi, h=h, v=v, l2=l2,
                ev=ev, EV=EV, g=g, nphi=nphi, npsi=npsi)


def exact_true_pieces(pr):
    """(C**) pieces using the TRUE top eigenvector v (full eigh). This is the
    'oracle' ceiling: no trial-vector approximation error, only the DK-vs-
    exact-resolvent gap. numerical-evidence (float64 eigh)."""
    A, Ophi, Opsi = pr["A"], pr["Ophi"], pr["Opsi"]
    h, v, g, x, ev, EV = pr["h"], pr["v"], pr["g"], pr["x"], pr["ev"], pr["EV"]
    F = v @ (Ophi @ v)
    r = Opsi @ v
    r = r - v * (v @ r)
    E = np.linalg.norm(r)
    c = EV.T @ r
    c[-1] = 0.0
    R_exact = np.sum(c[:-1] ** 2 / (h - ev[:-1]))   # <r,(h-A)^+ r>, positive branch
    exact_pos = 2 * x * x * R_exact
    dk_pos = 2 * x * x * E * E / g
    return dict(F=F, E=E, R=R_exact, exact_pos=exact_pos, dk_pos=dk_pos,
                ratio_true=exact_pos / (2 * F), ratio_dk=dk_pos / (2 * F),
                margin_factor=dk_pos / exact_pos)  # >1: DK overestimates pos


def trial_vector(pr, alpha, c=0.5):
    N, w = pr["N"], pr["w"]
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    a = alpha * w
    t = u * a * a / (a * a + (n + c) ** 2)
    return t / np.linalg.norm(t)


def trial_dk_pieces(pr, t, g_an):
    """DK-relaxed (C**) ratio for an explicit trial vector t, exactly as in
    lemmaC_compact_probe.trial_dk / lemmaC_cert_size.dk_ratio (E, F via
    genuine quadratic forms of the psi/phi kernels, NOT the residual of A --
    that residual-as-E-proxy substitution is known unsound for w>1, see
    routeRprime_width_scan.py Section 3 / team-lead warning)."""
    A, Ophi, Opsi, x, N = pr["A"], pr["Ophi"], pr["Opsi"], pr["x"], pr["N"]
    rho = t @ (A @ t)
    res = np.linalg.norm(A @ t - rho * t)
    s_bnd = min(res / g_an, 1.0) if g_an > 0 else np.inf
    Ft = t @ (Ophi @ t)
    e_phi = np.linalg.norm(Ophi @ t - Ft * t)
    Kt = Opsi @ t
    qpsi = t @ Kt
    e_psi = np.linalg.norm(Kt - qpsi * t)
    nphi = pr["nphi"]
    npsi = pr["npsi"]
    s = s_bnd
    F_lo = (1 - s * s) * Ft - 2 * s * e_phi - s * s * nphi
    E_hi = e_psi + s * (abs(qpsi) + npsi)
    ratio_dk = x * x * E_hi ** 2 / (g_an * F_lo) if (F_lo > 0 and g_an > 0) else np.inf
    return dict(rho=rho, res=res, s_bnd=s_bnd, Ft=Ft, e_phi=e_phi, e_psi=e_psi,
                qpsi=qpsi, nphi=nphi, npsi=npsi, F_lo=F_lo, E_hi=E_hi, ratio_dk=ratio_dk)


def gap_an(w, h):
    return 0.4 * h if w <= 2 else h - THETA_STAR


CELLS = [
    (0.30, 0.42, 5.0), (0.42, 0.60, 3.6), (0.60, 0.85, 2.7), (0.85, 1.20, 2.1),
    (1.20, 1.70, 1.8), (1.70, 2.40, 1.6), (2.40, 3.40, 1.5), (3.40, 4.80, 1.42),
    (4.80, 6.80, 1.35), (6.80, 9.60, 1.30), (9.60, 13.6, 1.27), (13.6, 20.0, 1.24),
]


# ---------------------------------------------------------------------------
# Section 1: exact R vs DK margin, TRUE v, wide grid -- is 15-18% cell-wide?
# ---------------------------------------------------------------------------
print("=" * 100)
print("Section 1: exact-R (Prop.1 ii, TRUE top eigenvector) vs DK relaxation, wide grid")
print("(oracle ceiling -- no trial-vector approximation error, only DK-vs-exact-resolvent gap)")
print("=" * 100)
grid1 = [0.3, 0.42, 0.6, 0.85, 1.0, 1.2, 1.7, 2.0, 2.4, 3.0, 3.4, 4.0, 4.8, 6.0, 6.8, 8.5, 9.6, 13.6, 17.0, 20.0]
print(f"{'w':>6} {'h':>9} {'g':>9} {'F':>11} {'E':>11} {'R':>11} {'ratio_true':>10} {'ratio_dk':>9} {'margin%':>8}")
margin_by_w = {}
for w in grid1:
    pr = build(w)
    piece = exact_true_pieces(pr)
    margin_pct = (piece["margin_factor"] - 1.0) * 100.0
    margin_by_w[w] = piece["margin_factor"]
    print(f"{w:6.2f} {pr['h']:9.4f} {pr['g']:9.4f} {piece['F']:11.4e} {piece['E']:11.4e} "
          f"{piece['R']:11.4e} {piece['ratio_true']:10.4f} {piece['ratio_dk']:9.4f} {margin_pct:8.2f}")

print()
print("Interpretation: margin% = 100*(dk_pos/exact_pos - 1); this is how much slack the")
print("Davis-Kahan relaxation 2x^3E^2/g gives away relative to the exact reduced-resolvent")
print("term 2x^3R, computed from the TRUE top eigenvector (no trial-vector error).")


# ---------------------------------------------------------------------------
# Section 2: 12-cell trial family, DK version, EXACT quadratic-form E/F
# (reproduces / re-derives the claimed worst ratio 0.585, using the fixed
#  E/F machinery -- not the residual-as-E-proxy bug)
# ---------------------------------------------------------------------------
print()
print("=" * 100)
print("Section 2: 12-cell freed-width trial family, DK-relaxed ratio (cell-uniform worst)")
print("=" * 100)
print(f"{'cell':>15} {'alpha':>6} {'worst w':>8} {'sinTh_bnd':>9} {'ratio_dk':>9} {'status':>6}")
cell_worst_dk = []
for wa, wb, alpha in CELLS:
    ws = np.geomspace(wa, wb, 11)
    worst = None
    for w in ws:
        pr = build(w)
        g_an = gap_an(w, pr["h"])
        t = trial_vector(pr, alpha)
        d = trial_dk_pieces(pr, t, g_an)
        if worst is None or d["ratio_dk"] > worst[1]["ratio_dk"]:
            worst = (w, d)
    w_worst, d_worst = worst
    cell_worst_dk.append((wa, wb, alpha, w_worst, d_worst["ratio_dk"], d_worst["s_bnd"]))
    status = "PASS" if d_worst["ratio_dk"] < 1 else "FAIL"
    print(f"[{wa:5.2f},{wb:5.2f}] {alpha:6.2f} {w_worst:8.3f} {d_worst['s_bnd']:9.4f} {d_worst['ratio_dk']:9.4f} {status:>6}")
worst_overall_dk = max(c[4] for c in cell_worst_dk)
print(f"\nWorst ratio_dk over all 12 cells: {worst_overall_dk:.4f}  (report claims 0.585)")


# ---------------------------------------------------------------------------
# Section 3: REAL (derived, not assumed) upper bound on R using ONLY the
# trial vector tau (no access to the true v beyond h, which the existing
# repo convention already treats as separately certified).
#
# Derivation. T := hI - A is PSD on the WHOLE space (not just v^perp): Tv=0,
# and every other eigenvalue of T is h-lambda_i >= h-lambda_2 = g > 0.
# r := O[psi]v - v<v,O[psi]v> satisfies r _|_ v, hence r in range(T). For a
# PSD T and r in range(T), for ANY vector y (checked by hand, does NOT
# require y _|_ v -- only uses T^+T = T T^+ = orthogonal projector onto
# v^perp and r _|_ v):
#     R := <r,T^+ r> = <y,Ty> + 2<y,rho> + <rho,T^+ rho>,   rho := r - Ty.
# Since <rho,T^+rho> <= ||rho||^2/g (T^+ has operator norm 1/g on v^perp,
# and rho _|_ v because r _|_ v and Ty's v-component cancels in the
# algebra), this is a valid identity + one relaxation.
#
# We pick y _|_ tau (computable) and bound the O(s) leakage from tau != v
# uniformly over the whole sin(angle) <= s cone:
#   <y,r(v)> <= |<y,O[psi]tau>| + 2 s npsi ||y||   =: B_r(y)   (safe/upper)
#   <y,Ty>, ||Ty||                                  exact (no v needed)
#   ||rho|| = ||r(v)-Ty|| <= ||r(v)||+||Ty|| <= E_hi+||Ty||   (triangle;
#       E_hi is the PRE-EXISTING tau-only DK bound on ||r(v)||)
# giving the one-sided (upper-bound direction throughout) estimate
#   R <= UB_combined(y) := -<y,Ty> + 2 B_r(y) + (E_hi+||Ty||)^2 / g_an.
# Check y=0: UB_combined(0) = E_hi^2/g_an, EXACTLY the existing DK bound.
# So minimizing UB_combined over any finite candidate set of y (always
# including y=0) can only match or IMPROVE on DK -- never worse, and stays
# a valid upper bound for every y tried (safe to numerically minimize).
#
# y ranges over a 1- or 2-vector Krylov-type space:
#   y0 := e_psi_vec/||e_psi_vec||  (= O[psi]tau - tau*qpsi, unit, _|_ tau)
#   y1 := Gram-Schmidt(T y0) against {tau, y0}
# with coefficients optimized numerically (grid + local refine).
#
# CAVEAT (explicit): this still uses the EXACT true h, exactly as the
# pre-existing DK machinery already does via g_an = gap_an(w, pr['h']).
# It does NOT need the true eigenvector v anywhere except through the
# pre-existing (already tau-only) E_hi bound. So it genuinely removes the
# "needs a linear solve/full diagonalization at v" objection for the
# resolvent term -- h was already assumed separately certified before this
# script even started.
# ---------------------------------------------------------------------------
print()
print("=" * 100)
print("Section 3: REAL trial-vector-only upper bound on R via 1-2 vector Krylov space")
print("=" * 100)


def trial_combined(pr, t, g_an):
    """Real (derived) upper bound on the (C**) ratio using only the trial
    vector t (plus the pre-existing E_hi/F_lo tau-only DK machinery, plus a
    small Krylov subspace for the resolvent term). Returns ratio_dk (m=0,
    identical to trial_dk_pieces), ratio_m1, ratio_m2 -- valid upper bounds
    satisfying ratio_m2 <= ratio_m1 <= ratio_dk by construction (minimizing
    over a strictly larger candidate set each time)."""
    A, x, h = pr["A"], pr["x"], pr["h"]
    Opsi = pr["Opsi"]
    npsi = pr["npsi"]
    d = trial_dk_pieces(pr, t, g_an)
    s, F_lo, E_hi = d["s_bnd"], d["F_lo"], d["E_hi"]
    ratio_dk = d["ratio_dk"]
    if F_lo <= 0 or g_an <= 0:
        return dict(ratio_dk=ratio_dk, ratio_m1=ratio_dk, ratio_m2=ratio_dk, s=s)

    Kt = Opsi @ t
    qpsi = t @ Kt
    e_psi_vec = Kt - qpsi * t
    e_psi = np.linalg.norm(e_psi_vec)
    if e_psi < 1e-14:
        return dict(ratio_dk=ratio_dk, ratio_m1=ratio_dk, ratio_m2=ratio_dk, s=s)

    y0 = e_psi_vec / e_psi                       # unit, exactly perp t
    Ty0 = h * y0 - A @ y0
    Kt0 = e_psi                                  # = y0 . Kt  (exact, y0 perp t)

    def UB_m1(c):
        c = abs(c)
        yTy = c * c * (y0 @ Ty0)
        cross = c * Kt0 + 2 * s * npsi * c        # Kt0>0 by construction
        Tynorm = c * np.linalg.norm(Ty0)
        rho_bound = E_hi + Tynorm
        return -yTy + 2 * cross + rho_bound ** 2 / g_an

    cs = np.concatenate(([0.0], np.geomspace(1e-6, 50.0, 300)))
    vals = np.array([UB_m1(c) for c in cs])
    j = int(np.argmin(vals))
    from scipy.optimize import minimize_scalar
    lo = cs[max(j - 1, 0)]
    hi = cs[min(j + 1, len(cs) - 1)]
    if hi <= lo:
        hi = lo + 1e-6
    ropt = minimize_scalar(UB_m1, bounds=(lo, hi + 1e-9), method="bounded")
    ub_m1 = min(vals[0], vals[j], float(ropt.fun))
    ratio_m1 = min(x * x * ub_m1 / (g_an * F_lo), ratio_dk)

    # m=2: second Krylov direction, Gram-Schmidt vs {t, y0}
    y1_raw = Ty0 - (Ty0 @ t) * t - (Ty0 @ y0) * y0
    n1 = np.linalg.norm(y1_raw)
    ratio_m2 = ratio_m1
    if n1 > 1e-10:
        y1 = y1_raw / n1
        Ty1 = h * y1 - A @ y1
        Kt1 = y1 @ Kt
        G00 = float(y0 @ Ty0)
        G01 = float(y0 @ Ty1)
        G11 = float(y1 @ Ty1)
        TG00 = float(Ty0 @ Ty0)
        TG01 = float(Ty0 @ Ty1)
        TG11 = float(Ty1 @ Ty1)

        def UB_m2(c0, c1):
            yTy = c0 * c0 * G00 + 2 * c0 * c1 * G01 + c1 * c1 * G11
            crossval = c0 * Kt0 + c1 * Kt1
            cross = np.abs(crossval) + 2 * s * npsi * np.hypot(c0, c1)
            Ty_sq = c0 * c0 * TG00 + 2 * c0 * c1 * TG01 + c1 * c1 * TG11
            Tynorm = np.sqrt(np.maximum(Ty_sq, 0.0))
            rho_bound = E_hi + Tynorm
            return -yTy + 2 * cross + rho_bound ** 2 / g_an

        c0_seed = cs[j]
        c0_grid = np.concatenate(([0.0, c0_seed], np.geomspace(1e-6, 30.0, 50)))
        c1_grid = np.concatenate(([0.0], np.geomspace(1e-6, 30.0, 50), -np.geomspace(1e-6, 30.0, 50)))
        best = ub_m1
        best_pt = (c0_seed, 0.0)
        for c0v in c0_grid:
            row = UB_m2(c0v, c1_grid)
            jm = int(np.argmin(row))
            if row[jm] < best:
                best = row[jm]
                best_pt = (c0v, c1_grid[jm])
        from scipy.optimize import minimize
        r2 = minimize(lambda p: UB_m2(p[0], p[1]), x0=np.array(best_pt), method="Nelder-Mead",
                       options=dict(xatol=1e-9, fatol=1e-13, maxiter=800))
        best = min(best, float(r2.fun))
        ratio_m2 = min(x * x * best / (g_an * F_lo), ratio_m1)

    return dict(ratio_dk=ratio_dk, ratio_m1=ratio_m1, ratio_m2=ratio_m2, s=s)


print("-- sanity check: feed the TRUE top eigenvector v as the 'trial' vector (s=0), --")
print("-- compare ratio_dk/ratio_m2 against Section 1's margin_factor --")
for w in (0.3, 1.0, 2.0, 8.5, 20.0):
    pr = build(w)
    g_an = gap_an(w, pr["h"])
    tv = pr["v"].copy()
    d = trial_combined(pr, tv, g_an)
    mf = margin_by_w.get(w)
    line = (f"  w={w:5.2f}  ratio_dk={d['ratio_dk']:.6e}  ratio_m1={d['ratio_m1']:.6e}  "
            f"ratio_m2={d['ratio_m2']:.6e}  ratio_dk/ratio_m2={d['ratio_dk']/d['ratio_m2']:.4f}")
    if mf:
        line += f"   (Sec.1 margin_factor={mf:.4f})"
    print(line)
print()

print(f"{'cell':>15} {'alpha':>6} {'worst w':>8} {'ratio_dk':>9} {'ratio_m1':>9} {'ratio_m2':>9} {'status':>6}")
cell_worst_real = []
for wa, wb, alpha, w_worst, r_dk, s_bnd in cell_worst_dk:
    # re-scan the whole cell under the REAL bound (the worst point may shift
    # once the numerator changes -- do not just reuse Section-2's worst w)
    ws = np.geomspace(wa, wb, 11)
    worst = None
    for w in ws:
        pr = build(w)
        g_an = gap_an(w, pr["h"])
        t = trial_vector(pr, alpha)
        d = trial_combined(pr, t, g_an)
        if worst is None or d["ratio_m2"] > worst[1]["ratio_m2"]:
            worst = (w, d)
    w2, d2 = worst
    cell_worst_real.append((wa, wb, alpha, w2, d2["ratio_dk"], d2["ratio_m1"], d2["ratio_m2"]))
    status = "PASS" if d2["ratio_m2"] < 1 else "FAIL"
    print(f"[{wa:5.2f},{wb:5.2f}] {alpha:6.2f} {w2:8.3f} {d2['ratio_dk']:9.4f} {d2['ratio_m1']:9.4f} {d2['ratio_m2']:9.4f} {status:>6}")

worst_overall_m1 = max(c[5] for c in cell_worst_real)
worst_overall_m2 = max(c[6] for c in cell_worst_real)
print(f"\nWorst REAL ratio over all 12 cells: m=1 -> {worst_overall_m1:.4f}   m=2 -> {worst_overall_m2:.4f}")
print(f"(reference: DK-only worst was {worst_overall_dk:.4f}; the EARLIER, now-superseded,")
print(" unverified projection was 0.4801)")
print("This is a genuinely derived, one-sided (upper-bound-direction) result, not an assumption")
print("that the true-v margin transfers -- see the derivation comment above for the exact chain")
print("of inequalities (all in the safe/PASS-supporting direction per CLAUDE.md's one-sided rule).")
print("STILL NOT a certificate: everything here is float64 (no interval arithmetic), and the")
print("Krylov coefficient optimization (grid + Nelder-Mead) is a numerical search, not a")
print("closed-form optimum -- a certificate would need either a closed-form optimal y or an")
print("interval-sound outer bound on the numerically-found minimum.")


# ---------------------------------------------------------------------------
# Section 4: CERT-B size re-estimate using the REAL m=2 bound (not the
# earlier, now-superseded, unverified projection).
# ---------------------------------------------------------------------------
print()
print("=" * 100)
print("Section 4: CERT-B size re-estimate, DK baseline vs REAL m=2 bound")
print("=" * 100)


def real_ratio_trunc(t, pr, g_an, m2=True):
    t = t / np.linalg.norm(t)
    d = trial_combined(pr, t, g_an)
    return (d["ratio_m2"] if m2 else d["ratio_dk"]), t


def min_support_real(pr, g_an, thresh, m2=True):
    v, N = pr["v"], pr["N"]
    for Ns in [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384]:
        t = np.zeros(N)
        t[:Ns] = v[:Ns]
        r, _ = real_ratio_trunc(t, pr, g_an, m2)
        if r < thresh:
            return Ns, r
    return N, None


def cell_halfwidth_real(pr, t_fixed, m2=True):
    x0 = pr["x"]
    dmax = 0.0
    for dd in [0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]:
        ok = True
        for xp in (x0 * (1 + dd), x0 * (1 - dd)):
            wp = np.sqrt(xp)
            prp = build(wp)
            g_an = gap_an(wp, prp["h"])
            r, _ = real_ratio_trunc(t_fixed, prp, g_an, m2)
            if not (r < 1.0):
                ok = False
                break
        if ok:
            dmax = dd
        else:
            break
    return dmax


def rescan_real(w_list, m2, label):
    print(f"\n-- {label} --")
    print(f"{'w':>5} {'h':>7} {'g_an':>7} | {'Ns':>4} {'ratio':>9} | {'half-width d':>12} {'cells/e-fold':>12}")
    rows = []
    for w in w_list:
        pr = build(w)
        g_an = gap_an(w, pr["h"])
        Ns, r = min_support_real(pr, g_an, 1.0, m2)
        v, N = pr["v"], pr["N"]
        t = np.zeros(N)
        t[:Ns] = v[:Ns]
        t = t / np.linalg.norm(t)
        dmax = cell_halfwidth_real(pr, t, m2)
        cells_per_efold = 1.0 / (2 * dmax) if dmax > 0 else np.inf
        rows.append((w, Ns, r, dmax, cells_per_efold))
        rv = r if r is not None else float("nan")
        print(f"{w:5.2f} {pr['h']:7.4f} {g_an:7.4f} | {Ns:4d} {rv:9.4f} | {dmax:12.3f} {cells_per_efold:12.1f}")
    return rows


w_scan = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.5]
w_scan_wide = w_scan + [11.0, 14.0, 17.0, 20.0]

rows_dk2 = rescan_real(w_scan, False, "DK baseline (re-run for apples-to-apples), (0.3, 8.5]")
rows_real = rescan_real(w_scan, True, "REAL m=2 bound, (0.3, 8.5]")
rows_real_wide = rescan_real(w_scan_wide, True, "REAL m=2 bound, (0.3, 20]")


def grid_summary(rows, lo, hi):
    """Piecewise-constant Riemann estimate of cells over log(w) -- crude,
    NOT a rigorous uniform-cell construction. numerical-evidence only."""
    ws = np.array([r[0] for r in rows])
    dens = np.array([r[4] for r in rows])
    lws = np.log(ws)
    order = np.argsort(lws)
    lws, dens = lws[order], dens[order]
    llo, lhi = np.log(lo), np.log(hi)
    total = 0.0
    edges = np.concatenate(([llo], (lws[:-1] + lws[1:]) / 2, [lhi]))
    for i in range(len(lws)):
        left = max(edges[i], llo)
        right = min(edges[i + 1], lhi)
        if right > left:
            total += dens[i] * (right - left)
    max_Ns = max(r[1] for r in rows)
    return total, max_Ns


ncells_dk2, ns_dk2 = grid_summary(rows_dk2, 0.3, 8.5)
ncells_real, ns_real = grid_summary(rows_real, 0.3, 8.5)
ncells_real_wide, ns_real_wide = grid_summary(rows_real_wide, 0.3, 20.0)

# Rational-constant bookkeeping. DK version certifies per cell: alpha + 7
# trial quantities (rho, res, Ft, e_phi, e_psi, nphi, npsi) = 8.
# REAL m=2 version needs, per cell, additionally: the 2 optimized Krylov
# coefficients (c0, c1) PLUS the 8 extra scalars that enter UB_m2 and are
# not already among the DK 7 (Kt0, Kt1, G00, G01, G11, TG00, TG01, TG11).
# Once (c0,c1) are fixed (certified) numbers, those 8 Gram-type scalars are
# themselves single certified evaluations, not free search variables -- so
# 2 + 8 = 10 extra rationals per cell on top of the DK 8.
RAT_PER_CELL_DK = 8
RAT_PER_CELL_REAL = 8 + 10

print()
print("=" * 100)
print("CERT-B size table (rough Riemann estimate over log w, NOT a rigorous cell count)")
print("=" * 100)
print(f"{'version':>28} {'domain':>12} {'~cells':>7} {'max Ns':>7} {'rat/cell':>8} {'~#rationals':>11}")
for label, ncells, ns, rat, dom in [
    ("DK baseline", ncells_dk2, ns_dk2, RAT_PER_CELL_DK, "(0.3,8.5]"),
    ("REAL m=2 bound", ncells_real, ns_real, RAT_PER_CELL_REAL, "(0.3,8.5]"),
    ("REAL m=2 bound (wide)", ncells_real_wide, ns_real_wide, RAT_PER_CELL_REAL, "(0.3,20]"),
]:
    print(f"{label:>28} {dom:>12} {ncells:7.1f} {ns:7d} {rat:8d} {ncells * rat:11.1f}")

print()
print("NOTE: this is a Riemann-sum extrapolation from 10-14 sample points, not an actual")
print("uniform-cell construction; treat the cell counts as order-of-magnitude, not final.")
print("Net effect: cell count changes by REAL/DK ratio below; per-cell rational count goes UP")
print("(8 -> 18) because of the 2 Krylov coefficients + 8 extra Gram scalars.")
if ncells_dk2 > 0:
    print(f"cells: DK={ncells_dk2:.1f} -> REAL={ncells_real:.1f}  (x{ncells_real / ncells_dk2:.3f})")
    print(f"rationals: DK={ncells_dk2 * RAT_PER_CELL_DK:.0f} -> REAL={ncells_real * RAT_PER_CELL_REAL:.0f}  "
          f"(x{(ncells_real * RAT_PER_CELL_REAL) / (ncells_dk2 * RAT_PER_CELL_DK):.3f})")
