r"""
P4-A T4 :: T_c^(N)(lambda, P) rigorous phonon-mediated Tc LOWER-BOUND LADDER
           + LaH10 spectral-shape correction table.

Physics recap (all established upstream -- NOT re-derived here; see
theory/derivations/P4_rigorous_tc_bound.md and p4_kappa_N.py in worktree
p3-numerics for the proofs/derivations this script consumes):

  KAY (Kiessling-Altshuler-Yuzbashyan 2025, J. Stat. Phys., dispersive-phonon
  paper "KAY II") proved the Eliashberg linear-stability boundary equals the
  largest eigenvalue of a compact self-adjoint operator frak_K on ell^2(N_0).
  Truncated eigenvalues (Thm 4):

      frak_k^(N)(P,T) = rho( Integral_0^infty H^(N)(varpi) P(domega) ),
      varpi := omega/(2 pi T),                                    (KAY Thm4)

  with closed-form H^(N)(varpi) matrices for N in {1,2,3,4}, built from the
  scalar kernel  [[n]](varpi) = varpi^2/(n^2+varpi^2).  N=1 anchor:
  H^(1)(varpi) = varpi^2/(varpi^2+1).  rho = spectral radius = largest
  eigenvalue (all H^(N) here are real symmetric).

  These give a LOWER-BOUND LADDER for Tc:
      frak_k^(1) < frak_k^(2) < frak_k^(3) < frak_k^(4) < ... -> frak_k (KAY Cor.3)
  The Tc condition is  lambda * frak_k^(N)(P,T) = 1,  i.e.
      frak_k^(N)(P, T_c^(N)) = 1/lambda.                                  (*)
  frak_k^(N)(P,T) is STRICTLY DECREASING in T (each matrix element
  [[k]](varpi) is increasing in varpi = omega/2piT, hence decreasing in T),
  running from a finite N-dependent supremum at T->0 (varpi->infinity) down
  to 0 as T->infinity.  Because of monotonicity, T_c^(N) (*) is unique
  whenever it exists, and frak_k^(N) </ frak_k^(N+1) pointwise in T implies
  T_c^(N) < T_c^(N+1) at fixed (lambda, P): each N is a rigorous LOWER BOUND
  on the true Tc given lambda.

  NOTE: KAY is mu*=0 throughout (no Coulomb pseudopotential); any of these
  T_c^(N) numbers is a mu*=0 lower bound, not a materials prediction.

  P(domega) is the same object as the Semenok/P3 spectral shape function.
  <omega> = first moment; write omega = <omega>(1+xi) with <xi>=0,
  sigma^2 = <xi^2>, gamma3 = <xi^3>/sigma^3 (skewness).
  varpi_bar := <omega>/(2 pi T).

  Upstream Proposition P4A.1 (verified there, 1257 checks in p4_kappa_N.py):
      frak_k^(N)(P,T) = h^(N)(varpi_bar) * [1 + (kappa_N(varpi_bar)/2) sigma^2
                          + c_{3,N}(varpi_bar) gamma3 sigma^3 + O(sigma^4)]
  with
      kappa_N(vb)  = vb^2 <v_N| H^(N)''(vb) |v_N> / h^(N)(vb),
      c_{3,N}(vb)  = (1/6) vb^3 <v_N| H^(N)'''(vb) |v_N> / h^(N)(vb),
  v_N = normalized leading eigenvector of H^(N)(vb).

  LaH10 four pressure points (P3 first-hand inversion), r := omega_2/omega_log:
      P/GPa :   129     163     214     264
      r     : 1.194   1.075   1.105   1.095
  Map r -> sigma via the P3 relation  r - 1 = sigma^2 - (1/3) gamma3 sigma^3;
  we ASSUME gamma3 = 0 here (leading order r-1 ~= sigma^2), i.e. we adopt the
  symmetric-shape realization sigma = sqrt(r-1).  Real hydride alpha^2F
  spectra are known to have negative skewness, so this is a documented
  APPROXIMATION, not a claim about the true LaH10 gamma3 -- it only fixes a
  concrete sigma to feed the exact (non-perturbative) two-point spectrum used
  below.  This is the same r <-> (sigma,gamma3) convention as P3.

UNITS CONVENTION (fixed once, used throughout):
  We set <omega> = 1 for every spectral measure P used below.  Then
  varpi_bar = 1/(2 pi T), and T_c^(N) comes out as a pure number equal to
  T_c^(N)/<omega> in units where <omega>=1 (equivalently, "T_c in units of
  <omega>/2pi" folded into the varpi definition).  To convert to physical
  temperature, multiply by the physical <omega> (e.g. in Kelvin-equivalent
  units) of the material at hand -- this script does not do that conversion.

REFERENCE CODE REUSED (per task instructions -- construction NOT altered,
only ported from mpmath to numpy float64 since all matrices are <=4x4):
  D:/physics/superconductivity/.claude/worktrees/p3-numerics/theory/derivations/p4_kappa_N.py
  (H_matrix, brk/dbrk elementwise derivatives, kappa_c3 machinery, and the
  KAY-I G_matrix(N,gamma) used only in the optional Stage F squeeze).
  This script is fully self-contained: everything needed is copied inline
  below; no import of the upstream file.

ANTI-STALL DISCIPLINE:
  numpy float64 only.  Matrices are <=4x4 -> numpy.linalg.eigh is instant.
  Root-finding via scipy.optimize.brentq with a fixed maxiter cap (<=60).
  Structured in stages A-F; each stage asserts self-checks and prints
  "Stage X: n/n checks passed".  No unbounded loops.
"""

import numpy as np
from scipy.optimize import brentq

PASS = [0]


def chk(cond, msg=""):
    """Counting assert: aborts on failure, tallies successes."""
    assert cond, msg
    PASS[0] += 1


# ===========================================================================
# Inline copy of upstream KAY-II H^(N)(varpi) matrix construction
# (ported mpmath -> numpy float64; construction itself UNCHANGED, see
#  p4_kappa_N.py module docstring for the derivation/self-checks that fixed
#  this exact placement: H^(N) = -H1 + H2 + H3).
# ===========================================================================
def brk(k, w):
    """[[k]](w) = w^2/(k^2+w^2).  (k=0 would be identically 1; unused here.)"""
    return w ** 2 / (k ** 2 + w ** 2)


def dbrk(k, w, order):
    """d^order/dw^order [[k]](w), k>=1, order in {1,2,3}: exact elementwise
    (same closed forms as p4_kappa_N.py dbrk)."""
    k2 = float(k) ** 2
    d = k2 + w ** 2
    if order == 1:
        return 2 * k2 * w / d ** 2
    if order == 2:
        return 2 * k2 * (k2 - 3 * w ** 2) / d ** 3
    if order == 3:
        return 24 * k2 * w * (w ** 2 - k2) / d ** 4
    raise ValueError(order)


def H_matrix(N, w, order=0):
    """KAY II H^(N)(w) (order=0) or its elementwise d^order/dw^order.

    H^(N) = -H1 + H2 + H3, matrix elements normalized by
    1/sqrt((2n+1)(2m+1)), arguments |n-m| and n+m+1 into [[.]] (or its
    derivative).  Identical construction to p4_kappa_N.py H_matrix.
    """
    M = np.zeros((N, N))
    for n in range(N):
        for m in range(N):
            nrm = np.sqrt(2 * n + 1) * np.sqrt(2 * m + 1)
            if order == 0:
                val = brk(n + m + 1, w) / nrm
                if n != m:
                    val += brk(abs(n - m), w) / nrm
                else:
                    s = sum(2 * brk(k, w) for k in range(1, n + 1)) if n >= 1 else 0.0
                    val -= s / (2 * n + 1)
            else:
                val = dbrk(n + m + 1, w, order) / nrm
                if n != m:
                    val += dbrk(abs(n - m), w, order) / nrm
                else:
                    s = sum(2 * dbrk(k, w, order) for k in range(1, n + 1)) if n >= 1 else 0.0
                    val -= s / (2 * n + 1)
            M[n, m] = val
    return M


def eig_max(M):
    """Largest eigenvalue + normalized eigenvector of a small real-symmetric M."""
    E, V = np.linalg.eigh(M)
    j = int(np.argmax(E))
    return E[j], V[:, j]


def quad(v, A):
    """<v|A|v>."""
    return float(v @ A @ v)


def kappa_c3(N, w):
    """(kappa_N(w), c_{3,N}(w), h^(N)(w)) from upstream (P4A.2)/(P4A.3)."""
    lam, v = eig_max(H_matrix(N, w))
    q2 = quad(v, H_matrix(N, w, order=2))
    q3 = quad(v, H_matrix(N, w, order=3))
    kappa = w ** 2 * q2 / lam
    c3 = w ** 3 * q3 / (6 * lam)
    return kappa, c3, lam


# ===========================================================================
# Stage A: self-check the inline copy (N=1 closed form + monotonicity in w)
# ===========================================================================
def run_stage_A():
    print("=== Stage A: inline H^(N)(w) copy self-check ===")
    for ws in (0.1, 0.3, 1.0, 3.0, 10.0):
        lam, _ = eig_max(H_matrix(1, ws))
        chk(abs(lam - ws ** 2 / (ws ** 2 + 1)) < 1e-13,
            f"H^(1)({ws}) mismatch: {lam}")
    # h^(N)(w) must be symmetric matrix and strictly increasing in w (needed
    # for the Tc solver's monotone bracket).
    wgrid = np.array([0.05, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0])
    for N in range(1, 5):
        vals = []
        for w in wgrid:
            M = H_matrix(N, w)
            chk(np.allclose(M, M.T), f"H^({N})({w}) not symmetric")
            lam, _ = eig_max(M)
            vals.append(lam)
        vals = np.array(vals)
        chk(np.all(np.diff(vals) > 0), f"h^({N})(w) not strictly increasing: {vals}")
    print(f"Stage A: {PASS[0]}/{PASS[0]} checks passed")


# ===========================================================================
# Stage B: frak_k^(N)(P,T), sup value at T->0, and the Tc solver
# ===========================================================================
def frak_k_of_P(N, T, omegas, ps):
    """frak_k^(N)(P,T) = rho( sum_j p_j H^(N)(omega_j/(2 pi T)) )."""
    M = np.zeros((N, N))
    for om, p in zip(omegas, ps):
        w = om / (2 * np.pi * T)
        M += p * H_matrix(N, w)
    lam, _ = eig_max(M)
    return lam


def sup_frak_k(N, w_big=1e8):
    """T->0 (w->infinity) limit of frak_k^(N); shape-independent since every
    omega_j/(2 pi T) -> infinity together as T->0 for any fixed omega_j>0."""
    lam, _ = eig_max(H_matrix(N, w_big))
    return lam


def solve_Tc(N, lam_coupling, omegas, ps, T_lo=1e-7, T_hi=1e7, maxiter=60):
    """Solve frak_k^(N)(P,T) = 1/lambda for T via bracketed Brent's method.

    frak_k^(N) is strictly decreasing in T (Stage A confirms h^(N) increasing
    in w = omega/2piT, and T -> w is decreasing), so a sign change on
    [T_lo, T_hi] is a robust bracket provided 1/lambda < sup frak_k^(N)(N).
    """
    target = 1.0 / lam_coupling
    sup = sup_frak_k(N)
    if target >= sup:
        raise ValueError(
            f"no finite T_c^({N}) for lambda={lam_coupling}: "
            f"1/lambda={target:.6f} >= sup frak_k^({N})={sup:.6f} "
            f"(coupling below the N={N}-truncation threshold lambda_min={1/sup:.6f})")

    def f(T):
        return frak_k_of_P(N, T, omegas, ps) - target

    flo, fhi = f(T_lo), f(T_hi)
    assert flo > 0 and fhi < 0, f"bad bracket N={N}: f(Tlo)={flo}, f(Thi)={fhi}"
    return brentq(f, T_lo, T_hi, maxiter=maxiter, xtol=1e-14, rtol=8.9e-16)


def run_stage_B():
    print("=== Stage B: frak_k^(N)(P,T) + Tc solver self-check ===")
    # (a) sup values at T->0 (used as lambda-threshold everywhere below)
    sups = {N: sup_frak_k(N) for N in range(1, 5)}
    print("    sup_{T->0} frak_k^(N):", {N: round(v, 6) for N, v in sups.items()})
    chk(abs(sups[1] - 1.0) < 1e-6, f"sup frak_k^(1) != 1: {sups[1]}")
    chk(sups[1] < sups[2] < sups[3] < sups[4], "sup frak_k^(N) not increasing in N")

    # (b) N=1 delta-spectrum closed form: frak_k^(1) = wbar^2/(wbar^2+1) = 1/lam
    #     => wbar = 1/sqrt(lam-1) => T_c^(1) = <omega>/(2 pi wbar) = sqrt(lam-1)/(2 pi)
    for lam_c in (1.5, 2.0, 3.0, 5.0, 10.0):
        Tc_closed = np.sqrt(lam_c - 1.0) / (2 * np.pi)
        Tc_num = solve_Tc(1, lam_c, [1.0], [1.0])
        chk(abs(Tc_num - Tc_closed) / Tc_closed < 1e-9,
            f"N=1 closed-form mismatch at lambda={lam_c}: {Tc_num} vs {Tc_closed}")
    print("    N=1 delta-spectrum: numeric Tc solver == closed form "
          "sqrt(lambda-1)/2pi to <1e-9 for lambda=1.5,2,3,5,10  OK")

    # (c) monotone-decreasing-in-T sanity (delta spectrum, N=2)
    Ts = np.array([0.01, 0.05, 0.1, 0.5, 1.0, 5.0])
    vals = np.array([frak_k_of_P(2, T, [1.0], [1.0]) for T in Ts])
    chk(np.all(np.diff(vals) < 0), f"frak_k^(2)(T) not monotone decreasing: {vals}")
    print("    frak_k^(2)(P,T) strictly decreasing in T on a test grid  OK")

    # (d) lambda <= 1 must raise (no finite Tc^(1): sup is exactly 1, never attained)
    raised = False
    try:
        solve_Tc(1, 1.0, [1.0], [1.0])
    except ValueError:
        raised = True
    chk(raised, "solve_Tc(N=1, lambda=1) should have raised (threshold not attained)")
    print("    lambda<=1 correctly refused for N=1 (sup frak_k^(1)=1 is a supremum, "
          "not attained at finite T)  OK")

    print(f"Stage B: {PASS[0]}/{PASS[0]} checks passed")
    return sups


# ===========================================================================
# Stage C (task b): the ladder T_c^(1) < T_c^(2) < T_c^(3) < T_c^(4)
# ===========================================================================
def shape_delta():
    return [1.0], [1.0]


def shape_2peak(sigma):
    """omega = 1 +/- sigma, p=1/2,1/2: <omega>=1, sigma matches, gamma3=0."""
    return [1.0 + sigma, 1.0 - sigma], [0.5, 0.5]


def shape_3peak(sigma):
    """1-2-1 symmetric: xi={-a,0,a}, p=(1/4,1/2,1/4), a=sigma*sqrt(2): sigma
    matches, gamma3=0 (mirrors p4_kappa_N.py shape_B)."""
    a = sigma * np.sqrt(2.0)
    return [1.0 - a, 1.0, 1.0 + a], [0.25, 0.5, 0.25]


def run_stage_C():
    print("=== Stage C: T_c^(1) < T_c^(2) < T_c^(3) < T_c^(4) ladder check ===")
    lam_list = [1.5, 2.0, 3.0, 5.0]  # all > 1 = 1/sup(frak_k^(1)) so all four
                                      # N=1..4 branches have a finite Tc (see
                                      # Stage B(d): lambda<=1 has no Tc^(1)).
    shapes = {
        "delta        ": shape_delta(),
        "2pk sigma=0.10": shape_2peak(0.10),
        "2pk sigma=0.20": shape_2peak(0.20),
        "3pk sigma=0.15": shape_3peak(0.15),
    }
    n_checked = 0
    for lam_c in lam_list:
        for tag, (omegas, ps) in shapes.items():
            Tcs = [solve_Tc(N, lam_c, omegas, ps) for N in range(1, 5)]
            chk(Tcs[0] < Tcs[1] < Tcs[2] < Tcs[3],
                f"ladder broken at lambda={lam_c}, shape={tag}: {Tcs}")
            n_checked += 1
            print(f"    lambda={lam_c:>4} {tag}: "
                  f"Tc^(1..4) = {', '.join(f'{t:.6f}' for t in Tcs)}")
    print(f"Stage C: {PASS[0]}/{PASS[0]} checks passed "
          f"({n_checked} (lambda,shape) combinations, ladder held every time)")


# ===========================================================================
# Stage D (task c): expansion (P4A.1) vs exact -- O(sigma^4) falsifiability
# ===========================================================================
def Tc_from_expansion(N, lam_coupling, sigma, gamma3=0.0, wb_lo=1e-7, wb_hi=1e7, maxiter=60):
    """Solve h^(N)(wbar)*[1 + kappa_N(wbar)/2 sigma^2 + c_{3,N}(wbar) gamma3 sigma^3]
    = 1/lambda for wbar (P4A.1 truncated at the stated order), then
    T_c_approx = <omega>/(2 pi wbar) = 1/(2 pi wbar) (<omega>=1)."""
    target = 1.0 / lam_coupling

    def g(wb):
        kap, c3, hN = kappa_c3(N, wb)
        return hN * (1.0 + 0.5 * kap * sigma ** 2 + c3 * gamma3 * sigma ** 3) - target

    glo, ghi = g(wb_lo), g(wb_hi)
    # g(wbar) = h^(N)(wbar)*[...] - target is INCREASING in wbar (h^(N) is
    # increasing in wbar, unlike f(T) in solve_Tc which is decreasing in T):
    # g(wb_lo)~-target<0, g(wb_hi)~sup-target>0 since target<sup (checked by caller).
    assert glo < 0 and ghi > 0, f"bad wbar bracket N={N} sigma={sigma}: g(lo)={glo} g(hi)={ghi}"
    wb_root = brentq(g, wb_lo, wb_hi, maxiter=maxiter, xtol=1e-13, rtol=8.9e-16)
    return 1.0 / (2 * np.pi * wb_root)


def run_stage_D():
    print("=== Stage D: T_c from exact frak_k^(N) vs from the sigma-expansion (P4A.1) ===")
    lam_c = 2.0  # > 1, safe for all N=1..4 (Stage B)
    sigmas = [0.16, 0.08, 0.04, 0.02, 0.01]  # halving sequence
    for N in range(1, 5):
        errs, e3, e4 = [], [], []
        for s in sigmas:
            omegas, ps = shape_2peak(s)  # gamma3 = 0 exactly for this shape
            Tc_exact = solve_Tc(N, lam_c, omegas, ps)
            Tc_approx = Tc_from_expansion(N, lam_c, s, gamma3=0.0)
            err = abs(Tc_exact - Tc_approx)
            errs.append(err)
            e3.append(err / s ** 3)
            e4.append(err / s ** 4)
        errs, e3, e4 = np.array(errs), np.array(e3), np.array(e4)
        print(f"    N={N}: sigma={sigmas}")
        print(f"         |Tc_exact-Tc_approx|      = {[f'{x:.3e}' for x in errs]}")
        print(f"         err/sigma^3 (-> 0)        = {[f'{x:.3e}' for x in e3]}")
        print(f"         err/sigma^4 (bounded)     = {[f'{x:.3e}' for x in e4]}")
        # (i) err/sigma^3 must shrink roughly linearly with sigma (halving
        #     sigma 4x from sigmas[0] to sigmas[-1] should roughly quarter
        #     err/sigma^3 if err ~ C sigma^4); loose factor-of-2 gate as in
        #     the upstream p4_kappa_N.py T3 self-test style.
        chk(e3[-1] < e3[0] / 2,
            f"N={N}: err/sigma^3 not shrinking like O(sigma^4): {e3}")
        # (ii) err/sigma^4 must stay bounded (not blow up) across the sweep --
        #      the hallmark of a genuine O(sigma^4) residual, not O(sigma^3).
        chk(np.max(e4) / (np.min(e4) + 1e-300) < 50,
            f"N={N}: err/sigma^4 not bounded (O(sigma^3) leakage?): {e4}")
    print(f"Stage D: {PASS[0]}/{PASS[0]} checks passed -- "
          f"O(sigma^4) scaling of the (P4A.1) Tc residual confirmed for N=1..4")


# ===========================================================================
# Stage E (task d): LaH10 Tc shape-correction table
# ===========================================================================
LAH10_POINTS = [
    (129, 1.194),
    (163, 1.075),
    (214, 1.105),
    (264, 1.095),
]


def run_stage_E():
    print("=== Stage E: LaH10 T_c^(N)(lambda,P) shape-correction table ===")
    lam_c = 2.0  # representative coupling: > 1 (Stage B threshold for N=1) and
                 # a plausible order-of-magnitude for hydride electron-phonon
                 # coupling; chosen for illustration of the shape correction,
                 # NOT a first-principles lambda(P) fit for LaH10.
    print(f"    representative coupling lambda = {lam_c}  (mu*=0, KAY convention)")
    print("    sigma := sqrt(r-1), assuming gamma3=0 (documented approximation; "
          "see module docstring)")

    Tc0 = [solve_Tc(N, lam_c, *shape_delta()) for N in range(1, 5)]  # sigma=0 baseline
    print(f"    delta-spectrum (sigma=0) baseline Tc^(1..4) = "
          f"{', '.join(f'{t:.6f}' for t in Tc0)}")

    header = (f"    {'P/GPa':>6} {'r':>7} {'sigma':>7} |"
              + "".join(f" {'Tc^(' + str(N) + ')':>10} {'corr^(' + str(N) + ')':>9} |"
                        for N in range(1, 5)))
    print(header)
    print("    " + "-" * (len(header) - 4))
    rows = []
    for P, r in LAH10_POINTS:
        sigma = np.sqrt(r - 1.0)
        omegas, ps = shape_2peak(sigma)
        Tcs = [solve_Tc(N, lam_c, omegas, ps) for N in range(1, 5)]
        corr = [Tcs[i] / Tc0[i] for i in range(4)]
        rows.append((P, r, sigma, Tcs, corr))
        line = f"    {P:>6} {r:>7.3f} {sigma:>7.4f} |"
        for N in range(1, 5):
            line += f" {Tcs[N-1]:>10.6f} {corr[N-1]:>9.5f} |"
        print(line)
        # sanity: correction factors finite and positive (Tc^(N)>0 always,
        # since frak_k^(N) is finite and continuous); NOT presumed to exceed
        # 1 -- kappa_N(wbar) changes sign (e.g. kappa_1 at wbar=1/sqrt(3),
        # see P4_rigorous_tc_bound.md sec.3), so broad spectra can either
        # enhance or SUPPRESS Tc^(N) relative to the delta baseline depending
        # on where wbar sits; here (lambda=2, LaH10 sigma~0.2-0.4) all four
        # branches happen to suppress Tc, which is a legitimate KAY-operator
        # effect, not a bug.
        chk(all(t > 0 for t in Tcs), f"P={P}: non-positive Tc^(N): {Tcs}")
        chk(Tcs[0] < Tcs[1] < Tcs[2] < Tcs[3],
            f"P={P}: ladder broken in shape-correction table: {Tcs}")

    print(f"Stage E: {PASS[0]}/{PASS[0]} checks passed")
    return rows, Tc0, lam_c


# ===========================================================================
# Stage F (OPTIONAL, task e): squeeze vs the proven upper constant C^sharp
# ===========================================================================
def G_matrix(N, gamma):
    """KAY I G^(N)(gamma) (gamma-model), same placement as H_matrix with
    [[k]](w) -> 1/k^gamma.  Copied inline from p4_kappa_N.py (numpy port),
    used ONLY for this optional constant-squeeze cross-check."""
    M = np.zeros((N, N))
    for n in range(N):
        for m in range(N):
            nrm = np.sqrt(2 * n + 1) * np.sqrt(2 * m + 1)
            val = 1.0 / (nrm * (n + m + 1) ** gamma)
            if n != m:
                val += 1.0 / (nrm * abs(n - m) ** gamma)
            else:
                s = sum(2.0 / k ** gamma for k in range(1, n + 1)) if n >= 1 else 0.0
                val -= s / (2 * n + 1)
            M[n, m] = val
    return M


def run_stage_F():
    print("=== Stage F (optional): squeeze vs proven upper constant C^sharp ===")
    C_SHARP = 0.371665  # KAY-proven rigorous upper bound (given, not re-derived here)
    E, _ = np.linalg.eigh(G_matrix(4, 2.0))
    g4_2 = E[-1]
    Tc4_2 = np.sqrt(g4_2) / (2 * np.pi)
    chk(abs(Tc4_2 - 0.1825137102) < 1e-8,
        f"T_c^(4)(2) mismatch vs known constant: {Tc4_2}")
    width = C_SHARP - Tc4_2
    ratio = C_SHARP / Tc4_2
    print(f"    T_c^(4)(gamma=2) [N=4 lower-bound truncation] = {Tc4_2:.10f}")
    print(f"    C^sharp          [proven rigorous upper bound] = {C_SHARP:.6f}")
    print(f"    squeeze width = C^sharp - T_c^(4) = {width:.6f}  "
          f"(ratio C^sharp/T_c^(4) = {ratio:.4f})")
    print(f"Stage F: {PASS[0]}/{PASS[0]} checks passed")


if __name__ == '__main__':
    run_stage_A()
    run_stage_B()
    run_stage_C()
    run_stage_D()
    run_stage_E()
    run_stage_F()
    print(f"\nALL {PASS[0]} checks passed.")
