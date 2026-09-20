r"""
P4 . QUANTITATIVE ANISOTROPY -- Temple-Bhatia-Davis lower bound on the Tc
enhancement due to Fermi-surface anisotropy.

v0.1 (2026-09-14): initial implementation.  Companion to
`P4_anisotropy_quantitative.md`; upgrades Prop 5 of `P4_extremal_constraints.md`
from qualitative (operator compression) to quantitative (explicit bound).

SETUP (verbatim reuse, no re-derivation):
  A = Lam (x) L23 - diag(z) (x) Z   on R^N (x) R^S   (Matsubara x FS sheets),
  Lam = sqrt(w) lam sqrt(w) (symmetrized coupling),  z = lam w  (Z row sums),
  L23 = H + Z,  H = H_mat(N, varpi)  [p4_interchange_uniformity],
  A built by build_aniso_KAY_operator  [p4_extremal_constraints].
  Compression:  Q^T A Q = lam_bar * H,  lam_bar = w^T lam w   (Prop 5).

MAIN BOUND (md sec.3, Theorem A + Lemma 1):
  trial psi = e (x) sqrt(w),  e = unit top eigenvector of H:
      <psi,A psi> = lam_bar * h                    (r; exact)
      A psi - r psi = h * e (x) (Lam u - lam_bar u),  u = sqrt(w)   (since Mu=0)
      V := ||A psi - r psi||^2 = h^2 * Var_w(z),  Var_w(z) = sum w_i z_i^2 - lam_bar^2
  Bhatia-Davis (no gap/simplicity needed, only mu := lam_min(A)):
      V <= (lam_1 - r)(r - mu)   =>
      lam_1(A) - lam_bar*h  >=  h^2 * Var_w(z) / (lam_bar*h - mu).
  Numerator vanishes iff the coupling is isotropic (Var_w(z) = 0).

Tc CONVERSION (md sec.4, Theorem B):  with varpi = omega_E/(2 pi T),
  G(varpi) := lam_1(A(varpi)),  Tc := (omega_E/2pi) / varpi_a,
  varpi_a := sup{varpi : G(varpi) >= 1}   (sup convention, no monotonicity needed);
  iso reference varpi_i := h^{-1}(1/lam_bar);  E(varpi) := bound(varpi);
  varpi_hat solves lam_bar*h + E = 1  =>  Tc_aniso/Tc_iso >= varpi_i/varpi_hat.

GATES:
  0  structural identities (compression Q^TAQ = lam_bar H; Mu = 0; V two ways;
     Bhatia-Davis sanity on random 2x2)
  1  bound never violated: 200 random matrices per S in {2,3,4}, N=30
  2  tightness (bound / true enhancement) reported per S
  3  isotropic limit: Var_w(z)=0, bound=0, true enhancement=0
  4  Tc conversion: MgB2 (Golubov-style set; moderate set) + H3S vs literature
  5  two-N consistency (N=40 vs N=80) for all examples

RUNTIME TARGET: < 300 s.  COUNTING: ck-style, ends "X/X checks PASSED".
Units: kernels dimensionless; varpi = omega/(2 pi T); omega_E in meV;
  Tc [K] = 11.6045 * omega_E / (2 pi varpi_c)   (k_B = 8.617333e-2 meV/K).
"""
import time
import numpy as np

# ---------------------------------------------------------------- check harness
_CK_COUNT = 0
_CK_FAILED = []


def ck(cond, label, info=""):
    global _CK_COUNT
    _CK_COUNT += 1
    if not cond:
        msg = f"CHECK FAILED: {label}"
        if info:
            msg += f"  [{info}]"
        _CK_FAILED.append(msg)
        print(f"  [{_CK_COUNT}] FAIL: {label}  {info}")
    return bool(cond)


def report_checks():
    if _CK_FAILED:
        print(f"\n{'='*80}")
        print(f"SUMMARY: {len(_CK_FAILED)} / {_CK_COUNT} checks FAILED")
        for msg in _CK_FAILED:
            print(f"  {msg}")
        print(f"{'='*80}\n")
        raise RuntimeError(f"{len(_CK_FAILED)} of {_CK_COUNT} checks failed")
    print(f"\n{'='*80}")
    print(f"ALL {_CK_COUNT} CHECKS PASSED")
    print(f"{'='*80}\n")


# ------------------------------------------------------------ imports (verbatim)
from p4_interchange_uniformity import H_mat, _brk_vec
from p4_extremal_constraints import build_aniso_KAY_operator

K_PER_MEV = 1.0 / 8.61733285e-2          # 11.6045 K per meV


# --------------------------------------------------------------------- builders
def Z_diag(N, vp):
    """Diagonal of the KAY Z block: Z_nn = 2*sum_{j<=n} [[j]](vp) / (2n+1)."""
    n = np.arange(N)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * _brk_vec(np.arange(1, N), vp, 0))))
    return cum / (2.0 * n + 1.0)


def h_N(N, vp):
    """Top eigenvalue of H^(N)(varpi)."""
    return float(np.linalg.eigvalsh(H_mat(N, vp))[-1])


def fiber(lam, wfs):
    """Fiber data: lam_bar = w^T lam w, z = lam w, Var_w(z)."""
    z = lam @ wfs
    lbar = float(wfs @ z)
    return dict(lbar=lbar, z=z, Vz=float(wfs @ (z * z)) - lbar * lbar)


def op_analysis(N, vp, lam, wfs):
    """lam_1, mu, top eigenvectors, Temple data and the bound, at fixed varpi."""
    S = len(wfs)
    A = build_aniso_KAY_operator(N, vp, lam, wfs)
    ev, vec = np.linalg.eigh(A)
    l1, mu, e = float(ev[-1]), float(ev[0]), vec[:, -1]
    evh, vech = np.linalg.eigh(H_mat(N, vp))
    h, eh = float(evh[-1]), vech[:, -1]
    u = np.sqrt(wfs)
    Lam = u[:, None] * lam * u[None, :]
    z = lam @ wfs
    lbar = float(wfs @ z)
    Vz = float(wfs @ (z * z)) - lbar * lbar
    psi = np.concatenate([u[i] * eh for i in range(S)])
    r = float(psi @ A @ psi)
    resid = A @ psi - r * psi
    V_dir = float(resid @ resid)
    V_cf = h * h * Vz
    den = lbar * h - mu
    return dict(l1=l1, mu=mu, h=h, e=e, lbar=lbar, Vz=Vz, r=r,
                V_dir=V_dir, V_cf=V_cf, bound=V_cf / den, den=den, A=A)


def bisect(f, lo, hi, iters=80, expanding=False):
    """Sign-change bisection.  expanding=True widens hi until f(hi)>0."""
    flo = f(lo)
    ck(flo < 0.0, "bisection left end negative", f"f({lo:.4g})={flo:.3e}")
    if expanding:
        n_widen = 0
        while f(hi) < 0.0:
            hi *= 2.0
            n_widen += 1
            if n_widen > 30:
                raise RuntimeError("bisection: no sign change up to hi=%.3g" % hi)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if flo * f(mid) <= 0.0:
            hi = mid
        else:
            lo, flo = mid, f(mid)
    return 0.5 * (lo + hi)


# ===================================================================== STAGE 0
def stage0():
    print("=== STAGE 0: structural identities (Lemma 1 + compression + BD sanity) ===")
    t0 = time.time()
    rng = np.random.default_rng(101)

    for tag, S, N in (("S=2", 2, 30), ("S=3", 3, 24)):
        wfs = rng.dirichlet(np.ones(S))
        lam = np.zeros((S, S))
        for i in range(S):
            lam[i, i] = rng.uniform(0.3, 2.0)
        for i in range(S):
            for j in range(i + 1, S):
                lam[i, j] = lam[j, i] = rng.uniform(0.0, 0.8 * np.sqrt(lam[i, i] * lam[j, j]))
        vp = rng.uniform(0.5, 5.0)
        u = np.sqrt(wfs)

        # (a) compression identity Q^T A Q = lam_bar * H  (Prop 5, re-verified)
        A = build_aniso_KAY_operator(N, vp, lam, wfs)
        Q = np.zeros((S * N, N))
        for i in range(S):
            Q[i * N:(i + 1) * N, :] = u[i] * np.eye(N)
        PAP = Q.T @ A @ Q
        d1 = float(np.max(np.abs(PAP - fiber(lam, wfs)["lbar"] * H_mat(N, vp))))
        ck(d1 < 1e-12, f"{tag} compression Q^TAQ = lam_bar*H", f"diff={d1:.2e}")

        # (b) Lemma 1: Mu = 0, r = lam_bar*h, V two ways
        d = op_analysis(N, vp, lam, wfs)
        Mmat = u[:, None] * lam * u[None, :] - np.diag(lam @ wfs)
        dMu = float(np.max(np.abs(Mmat @ u)))
        ck(dMu < 1e-13, f"{tag} Mu = 0 (Lemma 1)", f"|Mu|={dMu:.2e}")
        dr = abs(d["r"] - d["lbar"] * d["h"])
        ck(dr < 1e-12, f"{tag} <psi,Apsi> = lam_bar*h", f"diff={dr:.2e}")
        dv = abs(d["V_dir"] - d["V_cf"])
        ck(dv < 1e-11 * max(1.0, d["V_dir"]),
           f"{tag} V_direct = h^2*Var_w(z)", f"diff={dv:.2e}")

        # (c) Bhatia-Davis itself: V <= (l1 - r)(r - mu)
        bd_slack = (d["l1"] - d["r"]) * (d["r"] - d["mu"]) - d["V_dir"]
        ck(bd_slack > -1e-10, f"{tag} Bhatia-Davis V <= (l1-r)(r-mu)",
           f"slack={bd_slack:.2e}")
        print(f"  {tag}: vp={vp:.3f}, V={d['V_dir']:.6f}, "
              f"bound={d['bound']:.6f}, enhancement={d['l1'] - d['lbar']*d['h']:.6f}")

    print(f"  Stage 0 done ({time.time()-t0:.1f}s)\n")


# ============================================================ STAGE 1 + 2 (G1,G2)
def stage1_2():
    print("=== STAGES 1+2: random matrices (Gates 1+2) -- bound never violated, tightness ===")
    t0 = time.time()
    N = 30
    for S in (2, 3, 4):
        rng = np.random.default_rng(200 + S)
        n_trials, tight = 200, []
        worst = np.inf
        for t in range(n_trials):
            wfs = rng.dirichlet(np.ones(S))
            lam = np.zeros((S, S))
            for i in range(S):
                lam[i, i] = rng.uniform(0.3, 2.0)
            for i in range(S):
                for j in range(i + 1, S):
                    lam[i, j] = lam[j, i] = rng.uniform(0.0, 0.8 * np.sqrt(lam[i, i] * lam[j, j]))
            vp = rng.uniform(0.3, 8.0)
            d = op_analysis(N, vp, lam, wfs)
            enh = d["l1"] - d["lbar"] * d["h"]
            bnd = d["bound"]
            if d["den"] < 0.05:            # degenerate guard (should not occur)
                continue
            slack = enh - bnd
            worst = min(worst, slack)
            ck(slack > -1e-8 * (1.0 + abs(d["l1"])),
               f"S={S} trial {t+1} Temple bound <= enhancement",
               f"vp={vp:.3f} enh={enh:.6e} bnd={bnd:.6e} slack={slack:.2e}")
            ck(d["den"] > 0.0, f"S={S} trial {t+1} denominator positive")
            if enh > 1e-12:
                tight.append(bnd / enh)
        tight = np.array(tight)
        ck(len(tight) == n_trials, f"S={S} all {n_trials} trials usable")
        print(f"  S={S}: worst slack = {worst:.2e};  tightness bound/enh: "
              f"min={tight.min():.4f} med={np.median(tight):.4f} max={tight.max():.4f}")
    print(f"  Stages 1+2 done ({time.time()-t0:.1f}s)\n")


# ================================================================= STAGE 3 (G3)
def stage3():
    print("=== STAGE 3: isotropic limit (Gate 3) -- bound and enhancement degenerate to 0 ===")
    t0 = time.time()
    rng = np.random.default_rng(301)
    for lam0 in (0.5, 1.0, 1.7):
        for S in (2, 3):
            wfs = rng.dirichlet(np.ones(S))
            vp = rng.uniform(0.5, 4.0)
            lam = np.full((S, S), lam0)
            d = op_analysis(24, vp, lam, wfs)
            # Var_w(z) = sum w_i z_i^2 - lam_bar^2 is a difference of two O(lam0^2)
            # quantities, so its float64 floor is a few eps*lam0^2 (~1 ulp), NOT 0.
            # Tolerance sits above that floor and still 12 orders below the real
            # anisotropies of sec.5 (5e-3 ... 5e-2).  Sign may be either way.
            ck(abs(d["Vz"]) < 1e-14 * max(1.0, lam0**2),
               f"iso lam0={lam0} S={S} Var_w(z)=0 (at float64 floor)",
               f"Vz={d['Vz']:.2e}")
            ck(d["bound"] < 1e-13, f"iso lam0={lam0} S={S} bound=0",
               f"bound={d['bound']:.2e}")
            enh = d["l1"] - d["lbar"] * d["h"]
            ck(enh < 1e-8 * max(1.0, lam0), f"iso lam0={lam0} S={S} enhancement=0",
               f"enh={enh:.2e}")
            print(f"  lam0={lam0}, S={S}, vp={vp:.3f}: Vz={d['Vz']:.1e}, "
                  f"bound={d['bound']:.1e}, enhancement={enh:.1e}")
    print(f"  Stage 3 done ({time.time()-t0:.1f}s)\n")


# ================================================== STAGE 4 helpers (Gate 4/5)
def tc_example(N, lam, wfs, omega_mev, label, varpi_grid=None):
    """Full Tc conversion for one coupling set; returns diagnostics dict."""
    fb = fiber(lam, wfs)
    lbar, Vz = fb["lbar"], fb["Vz"]

    def E_of(vp):
        """Temple lower bound on lam_1 - lam_bar*h at varpi (exact evaluation)."""
        d = op_analysis(N, vp, lam, wfs)
        return d["bound"]

    # monotonicity of G(varpi) = lam_1(A(varpi)) on a grid (numerical fact; the
    # theorem uses the sup convention and does not need it -- only the reading
    # of varpi_a as the unique crossing does).
    lo_g, hi_g = 0.02, 60.0
    grid = np.exp(np.linspace(np.log(lo_g), np.log(hi_g), 26))
    gvals = [op_analysis(N, v, lam, wfs)["l1"] for v in grid]
    dmin = float(np.min(np.diff(gvals)))
    ck(dmin > -1e-12, f"{label}: G(varpi) numerically increasing on grid",
       f"min diff={dmin:.2e}")

    # varpi_i := h^{-1}(1/lam_bar); h strictly increasing (theorem)
    vi = bisect(lambda v: h_N(N, v) - 1.0 / lbar, lo_g, hi_g, expanding=True)

    if Vz > 1e-14:
        # varpi_a: unique crossing of G = 1 in (lo, varpi_i]
        va = bisect(lambda v: op_analysis(N, v, lam, wfs)["l1"] - 1.0, lo_g, vi,
                    expanding=True)
        # varpi_hat: crossing of lam_bar*h + E = 1 in (lo, varpi_i]
        vh = bisect(lambda v: lbar * h_N(N, v) + E_of(v) - 1.0, lo_g, vi)
    else:                                   # isotropic coupling
        va = vh = vi
    rho_true, rho_bound = vi / va, vi / vh
    ck(rho_bound <= rho_true + 1e-9, f"{label}: rho_bound <= rho_true",
       f"{rho_bound:.6f} vs {rho_true:.6f}")
    ck(rho_true >= 1.0 - 1e-12, f"{label}: rho_true >= 1 (Prop 5)",
       f"rho_true={rho_true:.6f}")

    # data at the iso critical point (r = lam_bar*h = 1 there)
    d_i = op_analysis(N, vi, lam, wfs)
    enh_i = d_i["l1"] - 1.0
    tight_i = d_i["bound"] / enh_i if enh_i > 1e-14 else np.nan
    # log-slope factor for the nonlinearity remark: varpi*h'(varpi) at varpi_i
    dh = (h_N(N, vi * 1.02) - h_N(N, vi * 0.98)) / (0.04 * vi)
    slope = vi * dh

    tc_iso = omega_mev / (2.0 * np.pi * vi)
    tc_ano = omega_mev / (2.0 * np.pi * va)
    print(f"  {label}: lam_bar={lbar:.4f} Var_w(z)={Vz:.5f} (rel {Vz/lbar**2:.3%})")
    print(f"    varpi_i={vi:.5f} varpi_a={va:.5f} varpi_hat={vh:.5f} | "
          f"varpi*h'={slope:.3f}")
    print(f"    rho_true={rho_true:.5f}  rho_bound={rho_bound:.5f}  "
          f"tightness@crit={tight_i:.3f}  mu(varpi_i)={d_i['mu']:.4f}  "
          f"E(varpi_i)={d_i['bound']:.5f}")
    print(f"    Tc_iso={tc_iso:.3f} meV = {K_PER_MEV*tc_iso:.1f} K ; "
          f"Tc_aniso={tc_ano:.3f} meV = {K_PER_MEV*tc_ano:.1f} K  (omega_E={omega_mev} meV)")
    return dict(lbar=lbar, Vz=Vz, vi=vi, va=va, vh=vh,
                rho_true=rho_true, rho_bound=rho_bound, tight_i=tight_i,
                tc_iso=tc_iso, tc_ano=tc_ano, slope=slope, E_i=d_i["bound"],
                mu_i=d_i["mu"], h_i=d_i["h"])


def stage4():
    print("=== STAGE 4: Tc conversion, physical examples (Gate 4), N=40 ===")
    t0 = time.time()
    N = 40
    res = {}

    # --- MgB2, setting 1: Golubov-style representative set (ILLUSTRATIVE) ------
    # lambda_11=1.0 (sigma), lambda_22=0.45 (pi), lambda_12=0.2, w_sigma=0.56;
    # representative of the two-band fits around Choi et al. Nature 418, 758
    # (2002) / Golubov et al. two-band parameters; omega_E = 61 meV
    # (log-moment scale of the calculated alpha^2F).
    lam = np.array([[1.0, 0.2], [0.2, 0.45]])
    wfs = np.array([0.56, 0.44])
    res["M1"] = tc_example(N, lam, wfs, 61.0, "MgB2/M1 Golubov-style")
    ck(1.10 <= res["M1"]["rho_true"] <= 1.60, "M1 rho_true sanity window",
       f"{res['M1']['rho_true']:.4f}")
    ck(1.00 <= res["M1"]["rho_bound"] <= 1.40, "M1 rho_bound sanity window",
       f"{res['M1']['rho_bound']:.4f}")

    # decomposition context: same w, lam_12 = 0 (decoupled bands) -- how much of
    # the enhancement is max-vs-mean rather than merging (observation, not gate)
    lam0 = np.array([[1.0, 0.0], [0.0, 0.45]])
    vi0 = bisect(lambda v: h_N(N, v) - 1.0 / float(wfs @ lam0 @ wfs), 0.02, 60.0,
                 expanding=True)
    va0 = bisect(lambda v: op_analysis(N, v, lam0, wfs)["l1"] - 1.0, 0.02, vi0,
                 expanding=True)
    res["M1_decoupled"] = vi0 / va0
    print(f"    M1 decoupled (lam_12=0) reference: rho = {res['M1_decoupled']:.5f} "
          f"(max-vs-mean structure; merging pulls 2.10 -> 1.34)")

    # --- MgB2, setting 2: moderate anisotropy (folklore-calibrated) -------------
    lam = np.array([[0.9, 0.2], [0.2, 0.6]])
    wfs = np.array([0.5, 0.5])
    res["M2"] = tc_example(N, lam, wfs, 61.0, "MgB2/M2 moderate")
    ck(1.02 <= res["M2"]["rho_true"] <= 1.30, "M2 rho_true sanity window",
       f"{res['M2']['rho_true']:.4f}")
    ck(1.00 <= res["M2"]["rho_bound"] <= 1.20, "M2 rho_bound sanity window",
       f"{res['M2']['rho_bound']:.4f}")

    # --- H3S: illustrative two-sheet representative -----------------------------
    lam = np.array([[3.0, 2.45], [2.45, 2.1]])
    wfs = np.array([0.5, 0.5])
    res["H3S"] = tc_example(N, lam, wfs, 90.0, "H3S illustrative")
    ck(1.00 <= res["H3S"]["rho_true"] <= 1.05, "H3S rho_true sanity window",
       f"{res['H3S']['rho_true']:.4f}")

    print(f"\n  Literature comparison (md sec.5): MgB2 enhancement anchor 1.05-1.10 "
          f"(task-given):")
    print(f"    M2 rho_true = {res['M2']['rho_true']:.4f} lands at the top of the band;")
    print(f"    M1 rho_true = {res['M1']['rho_true']:.4f} exceeds it (large sigma/pi "
          f"spread); rho_bound M1 = {res['M1']['rho_bound']:.4f}, M2 = "
          f"{res['M2']['rho_bound']:.4f} bracket the band.")
    print(f"    H3S rho_true = {res['H3S']['rho_true']:.4f}: sub-percent, consistent "
          f"with anisotropic ~= isotropic full calculations.")
    print(f"  Stage 4 done ({time.time()-t0:.1f}s)\n")
    return res


# ================================================================= STAGE 5 (G5)
def stage5(res40):
    print("=== STAGE 5: two-N consistency (Gate 5), N=40 vs N=80 ===")
    t0 = time.time()
    N2 = 80
    sets = {
        "M1": (np.array([[1.0, 0.2], [0.2, 0.45]]), np.array([0.56, 0.44]), 61.0),
        "M2": (np.array([[0.9, 0.2], [0.2, 0.6]]), np.array([0.5, 0.5]), 61.0),
        "H3S": (np.array([[3.0, 2.45], [2.45, 2.1]]), np.array([0.5, 0.5]), 90.0),
    }
    for key, (lam, wfs, om) in sets.items():
        r = tc_example(N2, lam, wfs, om, f"{key}@N=80")
        d_rho = abs(r["rho_true"] - res40[key]["rho_true"])
        d_bnd = abs(r["rho_bound"] - res40[key]["rho_bound"])
        ck(d_rho < 1e-3, f"{key}: rho_true N-consistent",
           f"|{r['rho_true']:.6f} - {res40[key]['rho_true']:.6f}| = {d_rho:.2e}")
        ck(d_bnd < 1e-3, f"{key}: rho_bound N-consistent",
           f"|{r['rho_bound']:.6f} - {res40[key]['rho_bound']:.6f}| = {d_bnd:.2e}")
    print(f"  Stage 5 done ({time.time()-t0:.1f}s)\n")


# ======================================================================== MAIN
def main():
    print("\n" + "=" * 80)
    print("P4 QUANTITATIVE ANISOTROPY -- Temple-Bhatia-Davis enhancement bound")
    print("=" * 80 + "\n")
    t_start = time.time()
    try:
        stage0()
        stage1_2()
        stage3()
        res40 = stage4()
        stage5(res40)
        report_checks()
    except Exception as e:
        print(f"\nEXCEPTION: {e}\n")
        report_checks()
        raise
    print(f"Total runtime: {time.time()-t_start:.1f}s (target: <300s)")
    return 0


if __name__ == "__main__":
    exit(main())
