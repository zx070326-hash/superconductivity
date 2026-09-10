r"""
P4-D . CERTIFICATION ROUND for the Coulomb pseudopotential inside the KAY pencil.

RUN:  python theory/derivations/p4d_coulomb_certify.py            (all stages)
      python theory/derivations/p4d_coulomb_certify.py A C D      (a subset)

Stages: A (mu*=0 self-check gate), B (operator-order layer), C (limit gap
delta), D (uniform gap g_alt at the nodes), E (interval closure), F (the N<N_c
segment), G (physics anchors).

================================================================================
WHAT IS CERTIFIED (read this first)
================================================================================
`P4D_coulomb_operator.md` v0.3 fixes the Tc criterion as loss of positive
definiteness of the OPERATOR PENCIL

    I - lambda K_ph(wbar) + mu* K_c(N_c),                                    (P)
    K_c = 2|u><u|,   u_n = (2n+1)^{-1/2} 1_{n<N_c},   N_c = floor(w_c/2piT),

with K_ph = H(wbar) the certified KAY II operator of `p4_interchange_uniformity`
and NO -D_{N_c} term (that variant was falsified; see the note's Sec. 1).
Dividing (P) by lambda > 0, the pencil is positive definite iff

    lam_1( M(wbar, mutilde, N_c) ) < 1/lambda,
    M := H(wbar) - mutilde K_c(N_c),        mutilde := mu*/lambda,           (M)

so the ONE-PARAMETER self-adjoint family M is the object that has to carry the
whole certificate machine.  ** mutilde = mu*/lambda is the certified parameter,
not mu* itself. **  The requested physical values mu* in {0.1, 0.13} are
reported at lambda = 1 (i.e. mutilde = mu*).  COVERING LEMMA (reviewer
2026-09-10, `notes/2026-09-10_P4D_review.md` sec.4.2): K_c >= 0 gives
M(mutilde) >= M(mutilde_0) for mutilde <= mutilde_0, hence lam_1(M(mutilde,N_c))
>= lam_1(M(mutilde_0,N_c)) for every N_c, every Rayleigh/2x2 trial bound is
monotone the same way, and alpha is mutilde-free -- so EVERY certificate proved
at a node mutilde_0 (delta, g_alt, interval closure) holds verbatim for all
mutilde in [0, mutilde_0].  The nodes therefore certify the whole range below
them; MUTILDES_EXT adds the lambda<1 corner node mutilde=0.26 (mu*=0.13,
lambda=0.5) to the interval closure of Stage E.

Three certificates are re-established for M, in the same discipline as the
mu*=0 file (float64 `eigvalsh` only as a TRIAL-VECTOR / SEARCH device, never
inside a bound; mpmath dps=30 for every quantity that enters a bound):

  delta(wbar, mutilde, N_c) = mu - alpha   <=  lam_1(M) - lam_2(M)      [Stage C]
  g_alt(wbar, mutilde)      = mu_2 - alpha <=  inf_{N>=2} gap_N(M)      [Stage D]
  the same over the whole INTERVAL of wbar                              [Stage E]

THE ONE STRUCTURAL GIFT (used everywhere, and it is exact):
    K_c = 2|u><u| >= 0   ==>   M = H - mutilde K_c  <=  H   in operator order
    ==>  lam_k(M) <= lam_k(H) for EVERY k (min-max), in particular
         lam_2(M) <= lam_2(H) <= alpha(wbar)
    ==>  the ENTIRE alpha side of the mu*=0 certificate (the Cholesky-deflation
         head bound + the Schur ||B||,||D|| envelopes) is inherited VERBATIM,
         uniformly in mutilde >= 0 and in N_c.  Nothing is recomputed, and no
         new tail analysis is needed: only the lam_1 side moves.
Combined with Poincare separation (lam_2(P_N M P_N) <= lam_2(M)) this makes the
alpha side of g_alt free as well.

WHAT IS *NOT* CERTIFIED, AND WHY (mapped, not hidden)
-----------------------------------------------------
Taking the infimum over ALL N_c at FIXED wbar, the pencil gap genuinely
COLLAPSES at the small-wbar (high-temperature) end: at wbar = 0.3, mutilde = 0.1
the worst window is N_c = 3, where lam_1(M) ~ 1e-4 and gap ~ 6e-5.  This is not
a certificate artifact -- Stage C prints the float64 truth next to the bound --
it is the physical statement that at T ~ 0.53 w_E a Coulomb repulsion of this
strength very nearly cancels the phonon attraction, so the top eigenvalue itself
dives into the cluster at 0.  Consequently the INTERVAL theorem is certified on
a SUB-interval [wbar_min(mutilde), 5], with wbar_min located by the certificate
itself (Stage E) and the collapse boundary of the true operator reported
alongside.  The mu*=0 theorem on all of [0.3,5] is recovered exactly as
mutilde -> 0.

NUMERICAL DISCIPLINE.  Same as the mu*=0 file, plus: the exhaustive N_c scan of
Stage C is EXHAUSTIVE, not a sample -- a trial vector supported in [0,Nv) sees
K_c(N_c) only through u_{min(N_c,Nv)}, so scanning N_c = 0..Nv covers every
integer window N_c >= 0.
"""

import io
import os
import re
import sys
import time
from contextlib import redirect_stdout

import numpy as np
from mpmath import mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import p4_interchange_uniformity as piu                              # noqa: E402
from p4_interchange_uniformity import (H_mat, H_matrix, ROUND_PAD,   # noqa: E402
                                       NV_RAY, WBARS, alpha_certificate_H,
                                       alpha_cell_ub, schur_full_seg, wbar_grid)

mp.dps = 30

PASS = [0]
_LAST = [0]

MUTILDES = (0.1, 0.13)      # the requested physical values, at lambda = 1
MUTILDES_EXT = (0.26,)      # lambda<1 corner (mu*=0.13 at lambda=0.5); Stage E interval
                            # closure only (its nodal wbar=1 is inside the collapsed
                            # prefix, so the nodal wbar>=1 assertions of C/D do not apply)


def chk(name, ok, detail=""):
    if not ok:
        raise AssertionError(f"CHECK FAILED: {name} {detail}")
    PASS[0] += 1


def _summary(tag, note=""):
    k = PASS[0] - _LAST[0]
    _LAST[0] = PASS[0]
    print(f"  [{tag}] {k} checks passed. {note}")


# ============================================================================
# THE COULOMB OPERATOR AND THE PENCIL
# ============================================================================

def u_vec(N, Nc):
    """u_n = (2n+1)^{-1/2} 1_{n<N_c}, zero-padded to length N (float64)."""
    u = np.zeros(N)
    k = max(0, min(int(Nc), N))
    u[:k] = 1.0 / np.sqrt(2.0 * np.arange(k) + 1.0)
    return u


def u_vec_mp(N, Nc):
    """The same vector in mpmath (exact to dps): used in every BOUND."""
    k = max(0, min(int(Nc), N))
    return [mp.one / mp.sqrt(2 * n + 1) if n < k else mp.zero for n in range(N)]


def K_c(N, Nc):
    """K_c = 2|u><u| restricted to the leading N x N block.

    Note the EXACT compression identity used all over Stage F:
        P_N K_c(N_c) P_N = K_c(min(N, N_c)),
    i.e. truncating the Matsubara index at N simply shrinks the Coulomb window
    to min(N, N_c).  It is a leading principal block for every N -- there is no
    N-dependence hidden in the entry formulas, exactly as for H."""
    u = u_vec(N, Nc)
    return 2.0 * np.outer(u, u)


def M_mat(N, w, order=0, mutilde=0.0, Nc=0):
    """The pencil operator M = H(wbar) - mutilde K_c(N_c) on the N x N block.

    `order` is passed through to H_mat (the Coulomb term is wbar-INDEPENDENT, so
    every wbar-derivative of M equals the corresponding derivative of H -- this
    is what makes the whole Stage-1/Stage-4 envelope machinery of the mu*=0 file
    apply to M verbatim, and it is asserted in Stage B).

    At mutilde = 0 this returns H_mat(N, w, order) BIT-FOR-BIT: the subtraction
    of 0.0 * K_c is exact in IEEE-754 (0.0*x is +-0.0 and y-(+-0.0) == y).
    Stage A routes ALL 1668 checks of the mu*=0 certificate through this builder
    to prove exactly that."""
    A = H_mat(N, w, order)
    if order != 0:
        return A - 0.0 * K_c(N, Nc)      # derivative of a constant is 0
    return A - mutilde * K_c(N, Nc)


def S_of(Nc):
    """S(N_c) = sum_{n<N_c} 1/(2n+1) = ||u||^2 = (1/2) ln N_c + 0.98176 + o(1).
    ||K_c|| = 2S is the Morel-Anderson logarithm (K_c is rank one, so 2S is its
    only nonzero eigenvalue)."""
    return float(np.sum(1.0 / (2.0 * np.arange(max(0, int(Nc))) + 1.0)))


# ============================================================================
# STAGE A -- mu* = 0 SELF-CHECK GATE (non-negotiable; the whole round's floor)
#
# The gate is NOT "reproduce a few numbers".  `p4_interchange_uniformity`'s
# module-global `H_mat` is REBOUND to the pencil builder at mutilde = 0 and the
# complete 6-stage certificate is re-run.  Every one of the 1668 checks --
# Stage 0 builders, Stage 1 Schur/HS envelopes, Stage 2 Cholesky deflation,
# Stage 3 Poincare, Stage 3B the 246-cell interval closure, Stage 4 derivative
# norms, Stage 5 the sigma^4 assembly -- is therefore evaluated on matrices
# produced by THIS file's pencil code path.  The printed output is then diffed
# line-by-line against an UNPATCHED run of the same module in the same process
# (timings stripped), so agreement is BIT-level on every printed digit.
#
# Both runs go through one `_run_all` helper, which resets every module-level
# accumulator first, so the two runs cannot contaminate each other.  (An earlier
# version spawned the baseline as a subprocess; that doubles the peak footprint
# and died with a numpy allocation failure on this box.  Same-process is both
# lighter and a strictly tighter comparison -- identical interpreter, identical
# libm, so any difference really is the builder.)
# ============================================================================

_TIMING = re.compile(r"\(\d+(\.\d+)?s\)")


def _strip(txt):
    return [_TIMING.sub("(TIME)", ln.rstrip()) for ln in txt.splitlines()]


def _run_all():
    """Run the whole mu*=0 certificate with a clean slate, capturing its output.
    Returns (printed_text, n_checks)."""
    piu.PASS[0] = 0
    piu._STAGE0[0] = 0
    for d in (piu.DELTA, piu.ALPHA, piu.MU2, piu.GALT, piu.GSTAR, piu.N0,
              piu.BINF, piu.SMAX):
        d.clear()
    piu.G_INTERVAL[0] = None
    piu._ENV_CACHE.clear()
    buf = io.StringIO()
    with redirect_stdout(buf):
        piu.run_stage0()
        piu.run_stage1()
        piu.run_stage2()
        piu.run_stage3()
        piu.run_stage3b()
        piu.run_stage4()
        piu.run_stage5()
    return buf.getvalue(), piu.PASS[0]


def run_stage_A():
    print("=== STAGE A: mu*=0 SELF-CHECK GATE (non-negotiable) ===")
    t0 = time.time()

    # A1. bit-for-bit identity of the builders, all orders, before anything else
    for N in (1, 5, 17, 40, 120):
        for w in (0.3, 0.5, 1.0 / np.sqrt(3.0), 1.0, 2.0, 3.0, 5.0):
            for order in (0, 1, 2, 3, 4, 5):
                for Nc in (0, 1, 7, N, 10 * N):
                    d = np.max(np.abs(M_mat(N, w, order, 0.0, Nc) - H_mat(N, w, order)))
                    chk("M_mat(mutilde=0) == H_mat BIT-FOR-BIT", d == 0.0,
                        f"N={N} w={w} order={order} Nc={Nc} d={d:.3e}")
            d2 = np.max(np.abs(M_mat(N, w, 0, 0.0, N) - H_matrix(N, w)))
            chk("M_mat(mutilde=0) == H_matrix (verbatim T5 builder)", d2 <= 1e-15,
                f"N={N} w={w} d={d2:.3e}")
    print(f"  A1: pencil builder identical to both certified phonon builders "
          f"({PASS[0]} bit-level checks)")

    # A2. the UNPATCHED certificate (baseline)
    print("  A2: running the UNPATCHED certificate (baseline) ...")
    base_txt, n_base = _run_all()
    chk("baseline total is the certified 1668", n_base == 1668, f"got {n_base}")
    print(f"      baseline: ALL {n_base} checks passed")

    # A3. the SAME certificate with H_mat rebound to the pencil at mu*=0
    print("  A3: re-running ALL of it with H_mat := M_mat(., mutilde=0) ...")
    orig_H_mat = piu.H_mat
    calls = [0]

    def patched(N, w, order=0):
        calls[0] += 1
        # Nc deliberately NONZERO and varying, so the Coulomb branch is really
        # exercised (and must still be exactly annihilated by mutilde = 0).
        return M_mat(N, w, order, 0.0, Nc=(N // 2) + 1)

    piu.H_mat = patched
    try:
        pat_txt, n_patched = _run_all()
    finally:
        piu.H_mat = orig_H_mat

    chk("patched re-run passes the SAME number of checks (1668)",
        n_patched == n_base == 1668, f"patched={n_patched} baseline={n_base}")
    # Every matrix the certificate touches came from the pencil builder.  The
    # count is deterministic (2873 on the v0.6 certificate); the assertion is
    # kept loose enough to survive an added check upstream but tight enough that
    # a silently-bypassed patch (which would give 0) cannot slip through.
    chk("the pencil builder was actually the one used", calls[0] > 2000,
        f"only {calls[0]} calls")
    print(f"      patched : {n_patched} checks passed, "
          f"{calls[0]} pencil-builder calls")

    # A4. line-by-line diff of the printed certificate (timings stripped)
    b_lines, p_lines = _strip(base_txt), _strip(pat_txt)
    chk("same number of printed certificate lines", len(b_lines) == len(p_lines),
        f"{len(b_lines)} vs {len(p_lines)}")
    ndig = 0
    for i, (a, b) in enumerate(zip(b_lines, p_lines)):
        chk("printed certificate line identical", a == b,
            f"line {i}:\n  base: {a!r}\n  penc: {b!r}")
        ndig += sum(c.isdigit() for c in a)
    print(f"      diff: {len(b_lines)} printed lines identical, {ndig} digits, "
          f"0 differences")
    print("  ==> GATE A PASSED: at mu*=0 the pencil reproduces the 1668-check "
          "certificate DIGIT-FOR-DIGIT.")
    _summary("Stage A", f"({time.time()-t0:.1f}s)")
    return n_base


# ============================================================================
# STAGE B -- the operator-order layer (Layer 1 of the note's Sec. 3)
#
# Everything here is an exact algebraic fact about K_c, asserted numerically:
#   (i)   K_c is rank one and PSD; its only nonzero eigenvalue is 2S(N_c),
#         and 2S = ln N_c + 1.9635 + o(1) (the Morel-Anderson logarithm).
#   (ii)  M(mutilde) <= M(mutilde') in operator order whenever mutilde >=
#         mutilde' >= 0  ==>  lam_k(M) is NONINCREASING in mutilde for every k.
#         With the Tc condition lam_1(M) = 1/lambda this is precisely the note's
#         Layer-1 claim "mu* lowers Tc", obtained with no certificate at all.
#   (iii) lam_2(M) <= lam_2(H) <= alpha(wbar): the alpha side is INHERITED.
#   (iv)  d^j/dwbar^j M = d^j/dwbar^j H for j >= 1 (K_c is wbar-independent), so
#         Stage 1's envelopes and Stage 4's B_2^inf/B_4^inf transfer verbatim.
#   (v)   NO monotonicity in N_c.  K_c(N_c+1) - K_c(N_c) is an indefinite rank-2
#         perturbation, so lam_1(M) is NOT monotone in the window size -- this is
#         why Stage C scans N_c exhaustively instead of bounding an endpoint.
# ============================================================================

def run_stage_B():
    print("\n=== STAGE B: operator-order layer (Layer 1: mu* lowers Tc, no "
          "certificate needed) ===")
    t0 = time.time()
    # Matrix-level structure (rank one / PSD / spectrum) on windows small enough
    # to materialize; the Morel-Anderson log is then continued analytically,
    # since 2S(N_c) needs no matrix at all.
    for Nc in (1, 2, 3, 10, 100, 1000):
        N = max(2 * Nc, 4)
        Kc = K_c(N, Nc)
        ev = np.linalg.eigvalsh(Kc)
        S = S_of(Nc)
        chk("K_c is PSD", ev[0] > -1e-12 * max(1.0, ev[-1]), f"Nc={Nc} {ev[0]:.3e}")
        chk("K_c is rank one", abs(ev[-2]) < 1e-12 * max(1.0, ev[-1]), f"Nc={Nc}")
        chk("the single nonzero eigenvalue of K_c is 2S(N_c)",
            abs(ev[-1] - 2 * S) < 1e-11 * max(1.0, S), f"Nc={Nc}")
        chk("K_c is symmetric", np.max(np.abs(Kc - Kc.T)) == 0.0, f"Nc={Nc}")
    for Nc in (10, 100, 1000, 10000, 100000):
        chk("2S = ln N_c + 1.9635 (Morel-Anderson log)",
            abs(2 * S_of(Nc) - (np.log(Nc) + 1.9635)) < 0.03,
            f"Nc={Nc} 2S={2*S_of(Nc):.6f}")
    print(f"  {'N_c':>7}  {'||K_c||=2S':>11}  {'ln N_c + 1.9635':>16}  {'rank':>5}")
    for Nc in (10, 100, 1000):
        print(f"  {Nc:>7}  {2*S_of(Nc):>11.6f}  {np.log(Nc)+1.9635:>16.6f}  "
              f"{np.linalg.matrix_rank(K_c(2*Nc, Nc)):>5}")
    for Nc in (10000, 100000):
        print(f"  {Nc:>7}  {2*S_of(Nc):>11.6f}  {np.log(Nc)+1.9635:>16.6f}  "
              f"{1:>5}")

    # (ii) monotone decrease of every eigenvalue in mutilde  =>  Layer 1
    for w in WBARS:
        for Nc in (1, 4, 17, 200):
            N = 200
            prev = None
            for mt in (0.0, 0.05, 0.1, 0.13, 0.2, 0.5, 1.0):
                E = np.linalg.eigvalsh(M_mat(N, w, 0, mt, Nc))
                if prev is not None:
                    chk("lam_k(M) nonincreasing in mutilde (Layer 1)",
                        np.all(E <= prev + 1e-12),
                        f"w={w} Nc={Nc} mt={mt} worst={np.max(E-prev):.3e}")
                prev = E

    # (iii) the inherited alpha side
    print(f"  {'wbar':>7}  {'alpha (mu*=0 cert)':>18}  "
          f"{'max_over(mutilde,N_c) lam_2(M) at N=400':>40}")
    for w in WBARS:
        alpha = float(alpha_certificate_H(w)[0])
        worst = -np.inf
        for mt in (0.0,) + MUTILDES + (0.5,):
            for Nc in (0, 1, 3, 8, 40, 200, 400):
                E = np.linalg.eigvalsh(M_mat(400, w, 0, mt, Nc))
                chk("lam_2(M) <= lam_2(H) (operator order)",
                    E[-2] <= float(np.linalg.eigvalsh(H_mat(400, w))[-2]) + 1e-12,
                    f"w={w} mt={mt} Nc={Nc}")
                chk("lam_2(M) <= alpha: the alpha side is INHERITED",
                    E[-2] <= alpha, f"w={w} mt={mt} Nc={Nc} lam2={E[-2]:.6e}")
                worst = max(worst, E[-2])
        print(f"  {w:>7.4f}  {alpha:>18.7f}  {worst:>40.7e}")

    # (iv) derivatives of M == derivatives of H
    for w in (0.3, 1.0, 5.0):
        for order in (1, 2, 3, 4, 5):
            for mt in MUTILDES:
                d = np.max(np.abs(M_mat(60, w, order, mt, 30) - H_mat(60, w, order)))
                chk("d^j M / dwbar^j == d^j H / dwbar^j (K_c is wbar-independent)",
                    d == 0.0, f"w={w} order={order}")

    # (v) NO operator order in N_c -- documented as an assertion, since it is the
    #     reason Stage C must scan the window exhaustively.
    found = 0
    for Nc in (2, 3, 5, 9, 33):
        N = 4 * Nc
        D = K_c(N, Nc + 1) - K_c(N, Nc)
        ev = np.linalg.eigvalsh(D)
        if ev[0] < -1e-12 and ev[-1] > 1e-12:
            found += 1
    chk("K_c(N_c+1) - K_c(N_c) is INDEFINITE for every window tested",
        found == 5, f"only {found}/5 indefinite")
    print("  K_c(N_c+1)-K_c(N_c) is indefinite (rank 2, one sign each) => lam_1(M)")
    print("  is NOT monotone in N_c; Stage C therefore scans N_c EXHAUSTIVELY.")
    _summary("Stage B", f"({time.time()-t0:.1f}s)")


# ============================================================================
# STAGE C -- PILLAR 2 for the pencil: delta(wbar, mutilde, N_c) = mu - alpha.
#
# mu <= lam_1(M):  Rayleigh quotient of a finitely-supported trial vector, in
#   mpmath.  The trial vector is the float64 top eigenvector of M^(Nv) -- it is
#   only a TRIAL vector (any vector gives a valid, weaker bound), so no float64
#   eigenvalue enters the certificate.  Because v is supported in [0, Nv),
#
#       <v, M(N_c) v> = <v, H v> - 2 mutilde <u_{min(N_c,Nv)}, v>^2
#
#   EXACTLY, for every N_c (the part of u beyond Nv never couples to v).  Hence
#   scanning N_c = 0 .. Nv is EXHAUSTIVE over all integer windows N_c >= 0.
#   The Coulomb part is evaluated with u in exact mpmath.
# alpha >= lam_2(M):  INHERITED from the mu*=0 certificate (Stage B (iii)).
# ============================================================================

def rayleigh_mu_pencil(w, mutilde, Nc, Nv=NV_RAY, Am=None, A_np=None):
    """Certified lower bound mu <= lam_1(M(w, mutilde, N_c)) (mpmath, exact on
    the float64 matrix entries, minus ROUND_PAD for entry rounding)."""
    if A_np is None:
        A_np = H_mat(Nv, w)
    M = A_np - mutilde * K_c(Nv, Nc)
    v = np.linalg.eigh(M)[1][:, -1]
    vm = [mp.mpf(float(x)) for x in v]
    if Am is None:
        Am = _mp_H(w, Nv)
    num = mp.zero
    for i in range(Nv):
        if vm[i] == 0:
            continue
        row = Am[i]
        acc = mp.zero
        for j in range(Nv):
            acc += row[j] * vm[j]
        num += vm[i] * acc
    den = sum(x * x for x in vm)
    um = u_vec_mp(Nv, Nc)
    uv = sum(um[i] * vm[i] for i in range(Nv))
    num -= 2 * mp.mpf(float(mutilde)) * uv * uv
    return num / den - ROUND_PAD


_MPH = {}


def _mp_H(w, Nv):
    key = (float(w), Nv)
    if key not in _MPH:
        if len(_MPH) > 3:
            _MPH.clear()
        A = H_mat(Nv, w)
        _MPH[key] = [[mp.mpf(float(A[i, j])) for j in range(Nv)] for i in range(Nv)]
    return _MPH[key]


DELTA_P = {}


def run_stage_C():
    print("\n=== STAGE C: PILLAR 2 for the pencil -- delta = mu - alpha, "
          "N_c scanned EXHAUSTIVELY ===")
    t0 = time.time()
    print(f"  Rayleigh support Nv={NV_RAY} (=> N_c = 0..{NV_RAY} covers EVERY "
          f"integer window); alpha inherited from the mu*=0 certificate.")
    print(f"  {'wbar':>7}  {'mutilde':>8}  {'alpha':>9}  {'worst N_c':>9}  "
          f"{'mu (bound)':>11}  {'delta (bound)':>13}  {'true gap @N=400':>15}  ok")
    for mt in (0.0,) + MUTILDES:
        for w in WBARS:
            alpha = alpha_certificate_H(w)[0]
            A_np = H_mat(NV_RAY, w)
            Am = _mp_H(w, NV_RAY)
            best = None
            for Nc in range(0, NV_RAY + 1):
                mu = rayleigh_mu_pencil(w, mt, Nc, NV_RAY, Am, A_np)
                if best is None or mu < best[0]:
                    best = (mu, Nc)
            mu, Nc_w = best
            delta = float(mu - alpha)
            DELTA_P[(w, mt)] = (delta, Nc_w, float(mu), float(alpha))
            # float64 truth at the worst window (EVIDENCE, never in a bound)
            E = np.linalg.eigvalsh(M_mat(400, w, 0, mt, Nc_w))
            chk("mu <= lam_1(M) numerically", float(mu) <= E[-1] + 1e-9,
                f"w={w} mt={mt} Nc={Nc_w}")
            chk("delta <= true gap numerically", delta <= E[-1] - E[-2] + 1e-9,
                f"w={w} mt={mt}")
            if mt == 0.0:
                # Non-vacuous regardless of whether Stage A has run: recompute
                # the mu*=0 file's own Rayleigh bound and demand agreement.
                chk("mutilde=0 reproduces the certified delta(wbar) exactly",
                    abs(float(mu) - float(piu.rayleigh_mu(w))) < 1e-15, f"w={w}")
                chk("mutilde=0 selects the empty window N_c=0", Nc_w == 0, f"w={w}")
            ok = "YES" if delta > 0 else "no"
            print(f"  {w:>7.4f}  {mt:>8.3f}  {float(alpha):>9.6f}  {Nc_w:>9}  "
                  f"{float(mu):>11.7f}  {delta:>13.6f}  {E[-1]-E[-2]:>15.6e}  {ok}")
    # the mu*=0 row must be positive everywhere (that is the old theorem)
    for w in WBARS:
        chk("delta > 0 at mutilde = 0 (the mu*=0 theorem, recovered)",
            DELTA_P[(w, 0.0)][0] > 0, f"w={w}")
    print("\n  Positivity boundary in wbar (worst-case over ALL windows N_c):")
    for mt in MUTILDES:
        good = [w for w in WBARS if DELTA_P[(w, mt)][0] > 0]
        bad = [w for w in WBARS if DELTA_P[(w, mt)][0] <= 0]
        print(f"    mutilde={mt}: delta>0 at wbar in {good};  FAILS at {bad}")
        chk("delta > 0 at every node with wbar >= 1", all(DELTA_P[(w, mt)][0] > 0
                                                          for w in WBARS if w >= 1.0),
            f"mt={mt}")
    print("  (the small-wbar failures are REAL, not slack: the 'true gap' column")
    print("   collapses to ~1e-4 there too -- see the module docstring.)")
    _summary("Stage C", f"({time.time()-t0:.1f}s)")


# ============================================================================
# STAGE D -- PILLAR 3 for the pencil: inf_{N>=2} gap_N >= g_alt = mu_2 - alpha.
#
# Both ingredients of the mu*=0 Poincare route survive, and neither needs any
# positivity (Krein-Rutman) argument -- which is exactly why our route, and not
# KAY's, carries the Coulomb term:
#   (i)  POINCARE SEPARATION on the compression M_N = P_N M P_N:
#          lam_2(M_N) <= lam_2(M) <= lam_2(H) <= alpha.
#        The middle step is the operator order of Stage B (iii); the outer steps
#        are unchanged from mu*=0.  So the alpha side is free AND N_c-uniform.
#   (ii) CAUCHY INTERLACING on nested principal blocks: lam_1(M_N) is
#        nondecreasing in N.  This is where the compression identity
#        P_N K_c(N_c) P_N = K_c(min(N,N_c)) matters: M_N really is the leading
#        principal block of the single fixed operator M, so interlacing applies
#        for EVERY N, including N < N_c (Stage F).
#   ==>  inf_{N>=2} gap_N >= lam_1(M_2) - alpha, and lam_1(M_2) is bounded below
#        in CLOSED FORM on the 2x2 block, where the Coulomb term only sees
#        min(N_c,2) in {0,1,2}: three cases, exhaustive in N_c, take the min.
# ============================================================================

U2_CASES = {0: (mp.zero, mp.zero),
            1: (mp.one, mp.zero),
            2: (mp.one, mp.one / mp.sqrt(3))}


def mu2_pencil_cell(a, b, mutilde):
    """RIGOROUS lower bound on lam_1(M^(2)(w)) valid for EVERY w in [a,b] and
    EVERY window N_c, in closed form with ZERO Lipschitz slack.

    On the 2x2 block, H^(2)(w) = [[ [[1]], ([[1]]+[[2]])/sqrt3 ],
                                  [ .    , ([[3]]-2[[1]])/3    ]]
    and K_c^(2) = 2 u u^T with u = (u_0, u_1) determined by min(N_c,2) alone.
    Fix v = the top eigenvector of M^(2)(a) for that case (any fixed vector is
    admissible).  Then

        Q(w) = k1 [[1]](w) + k2 [[2]](w) + k3 [[3]](w) - 2 mutilde <u,v>^2/||v||^2

    with the SAME k1,k2,k3 as the mu*=0 file, and the Coulomb term is
    w-INDEPENDENT (K_c does not depend on wbar).  Every [[k]](w) increases in w,
    so k1,k2,k3 >= 0 implies Q is increasing and Q(w) >= Q(a) exactly.

    ** THE k-SIGN TEST IS LOAD-BEARING, NOT COSMETIC. **  Unlike the mu*=0 case
    (where k1 ~ 1.34..1.44 always), the repulsion tips the 2x2 top eigenvector to
    a sign change at small wbar, k2 goes negative, Q stops being increasing, and
    the returned value is then NOT a valid bound on the cell.  An adversarial
    sweep (`_scratch_ir/p4d_adversarial.py`, 400 random cells) found 276
    violations of the cell bound -- and EVERY ONE of them had a failing k-sign,
    with zero violations among the 1364 k-sign-valid cells tested.  So callers
    MUST check the returned signs; Stage E does, and excludes those cells from
    the certified interval (that exclusion is what sets wbar_min).

    NOTE the degenerate case a = b: there is no cell to cross, the returned value
    is the EXACT 2x2 top eigenvalue minus ROUND_PAD, and it is a valid bound
    whatever the k-signs do.  That is why Stage D reports the nodal g_alt and the
    k-signs in SEPARATE columns -- only Stage E's interval closure needs both.

    ** k-SIGNS ARE TAKEN OVER ALL THREE CASES (reviewer item 11, 2026-09-10). **
    min_case lam_case(a) is a bound on the whole cell for every N_c only if
    EVERY case's own Q_case is increasing, i.e. every case's k-signs hold -- not
    just the minimising case's.  The returned k1,k2,k3 are therefore the
    componentwise MINIMUM over the three cases (conservative: all-nonnegative
    iff all nine signs hold).  The reviewer's 1781-cell sweep found no cell
    where the minimising case passes and another fails, so no number changes.
    Returns (bound, k1, k2, k3, worst_case_Nc)."""
    aw = mp.mpf(float(a))
    b1 = aw * aw / (1 + aw * aw)
    b2 = aw * aw / (4 + aw * aw)
    b3 = aw * aw / (9 + aw * aw)
    h00, h01, h11 = b1, (b1 + b2) / mp.sqrt(3), (b3 - 2 * b1) / 3
    mt = mp.mpf(float(mutilde))
    best = None
    kmin = None
    for case, (u0, u1) in U2_CASES.items():
        m00 = h00 - 2 * mt * u0 * u0
        m01 = h01 - 2 * mt * u0 * u1
        m11 = h11 - 2 * mt * u1 * u1
        tr, det = m00 + m11, m00 * m11 - m01 * m01
        lam = (tr + mp.sqrt(tr * tr - 4 * det)) / 2
        # eigenvector of the 2x2 (m01 != 0 for every case we use)
        v0, v1 = (m01, lam - m00) if m01 != 0 else (mp.one, mp.zero)
        nn = v0 * v0 + v1 * v1
        k1 = (v0 * v0 + 2 * v0 * v1 / mp.sqrt(3) - 2 * v1 * v1 / 3) / nn
        k2 = (2 * v0 * v1 / mp.sqrt(3)) / nn
        k3 = (v1 * v1 / 3) / nn
        kmin = (k1, k2, k3) if kmin is None else (min(kmin[0], k1),
                                                  min(kmin[1], k2),
                                                  min(kmin[2], k3))
        if best is None or lam < best[0]:
            best = (lam - ROUND_PAD, case)
    return (best[0], kmin[0], kmin[1], kmin[2], best[1])


GALT_P = {}


def run_stage_D():
    print("\n=== STAGE D: PILLAR 3 for the pencil -- inf_{N>=2} gap_N >= "
          "mu_2 - alpha ===")
    t0 = time.time()
    print("  (i)  Poincare: lam_2(M_N) <= lam_2(M) <= lam_2(H) <= alpha  [alpha FREE]")
    print("  (ii) Cauchy:   lam_1(M_N) nondecreasing in N  =>  >= lam_1(M_2) >= mu_2")
    print("  mu_2 exact on the 2x2, min over the three window cases min(N_c,2)"
          " in {0,1,2}")
    print("  (the k-sign column is NOT needed for these NODAL values -- a=b is the")
    print("   exact 2x2 eigenvalue.  It is reported because Stage E's CELLS need it.)")
    print(f"  {'wbar':>7}  {'mutilde':>8}  {'mu_2':>11}  {'case':>4}  {'alpha':>9}  "
          f"{'g_alt':>10}  {'min_N min_Nc gap_N (true)':>25}  {'ok':>3}  {'k>=0':>5}")
    for mt in (0.0,) + MUTILDES:
        for w in WBARS:
            mu2, k1, k2, k3, case = mu2_pencil_cell(w, w, mt)
            alpha = alpha_certificate_H(w)[0]
            galt = float(mu2 - alpha)
            GALT_P[(w, mt)] = (galt, float(mu2), case, bool(k1 >= 0 and k2 >= 0
                                                            and k3 >= 0))
            # float64 truth: min over N and over ALL windows (EVIDENCE ONLY)
            worst = np.inf
            for N in (2, 3, 4, 6, 8, 12, 16, 24, 32, 64, 140, 400):
                A = H_mat(N, w)
                for Nc in range(0, N + 1):
                    E = np.linalg.eigvalsh(A - mt * K_c(N, Nc))
                    worst = min(worst, E[-1] - E[-2])
            chk("mu_2 == lam_1(M^(2)) of the worst window (numeric)",
                abs(float(mu2 + ROUND_PAD)
                    - min(float(np.linalg.eigvalsh(M_mat(2, w, 0, mt, c))[-1])
                          for c in (0, 1, 2))) < 1e-11, f"w={w} mt={mt}")
            if galt > 0:
                chk("g_alt <= every directly computed gap_N over every window",
                    worst >= galt - 1e-12, f"w={w} mt={mt} true={worst:.6f} "
                                           f"g_alt={galt:.6f}")
            ok = "YES" if galt > 0 else "no"
            print(f"  {w:>7.4f}  {mt:>8.3f}  {float(mu2):>11.7f}  {case:>4}  "
                  f"{float(alpha):>9.6f}  {galt:>10.6f}  {worst:>25.6e}  {ok:>3}  "
                  f"{('yes' if GALT_P[(w, mt)][3] else 'NO'):>5}")
    for w in WBARS:
        chk("g_alt > 0 at mutilde=0 (the mu*=0 theorem, recovered)",
            GALT_P[(w, 0.0)][0] > 0, f"w={w}")
        # Non-vacuous standalone: compare against the mu*=0 file's own mu_2.
        chk("mutilde=0 mu_2 matches the certified mu2_cell_lb exactly",
            abs(GALT_P[(w, 0.0)][1] - float(piu.mu2_cell_lb(w, w)[0])) < 1e-15,
            f"w={w}")
    for mt in MUTILDES:
        for w in WBARS:
            if w >= 1.0:
                chk("g_alt > 0 at every node wbar >= 1",
                    GALT_P[(w, mt)][0] > 0,
                    f"w={w} mt={mt} g_alt={GALT_P[(w, mt)][0]:.6f}")
    print("  ANSWER (ii) of the task: g_alt(0.1) and g_alt(0.13) at the nodes are "
          "the 'g_alt' column;")
    print("  the worst node with wbar>=1 is wbar=1: "
          f"g_alt(0.1)={GALT_P[(1.0, 0.1)][0]:.6f}, "
          f"g_alt(0.13)={GALT_P[(1.0, 0.13)][0]:.6f}.")
    _summary("Stage D", f"({time.time()-t0:.1f}s)")


# ============================================================================
# STAGE E -- the INTERVAL closure for the pencil (Stage 3B of the mu*=0 file).
#
# Every ingredient of the 246-cell closure is reused UNCHANGED except mu_2:
#   ||B||, ||D||, lam_2(head A), the Lipschitz slack L = schur_full_seg(a,b,1)
#     -- all belong to the alpha side, which Stage B (iii) shows is inherited
#        verbatim (and is N_c- and mutilde-uniform).  `alpha_cell_ub` is called
#        as-is, with the SAME Cholesky c values as the mu*=0 run.
#   mu_2 -- replaced by `mu2_pencil_cell`, which is still EXACT per cell (zero
#        Lipschitz slack) because K_c carries no wbar dependence.
# The output is the largest wbar_min on the certified grid such that the closure
# holds on [wbar_min, 5]; at mutilde = 0 it must return the whole [0.3, 5] and
# the published infimum 0.093225.
# ============================================================================

G_INT = {}


def run_stage_E():
    print("\n=== STAGE E: INTERVAL closure for the pencil on [wbar_min, 5] ===")
    t0 = time.time()
    grid = wbar_grid()
    print(f"  grid: {len(grid)} geometric points on [0.3,5] (the certified mu*=0 "
          f"grid, unchanged)")
    cs = []
    for w in grid:
        _, ok, c, _, _ = alpha_certificate_H(float(w))
        chk("Cholesky SPD on the refined grid (unchanged from mu*=0)", ok, f"w={w}")
        cs.append(c)
    Ls, als = [], []
    for i in range(len(grid) - 1):
        a, b = float(grid[i]), float(grid[i + 1])
        L = schur_full_seg(a, b, 1)
        Ls.append(L)
        als.append(alpha_cell_ub(a, b, cs[i], cs[i + 1], L))

    print(f"  {'mutilde':>8}  {'wbar_min':>9}  {'cells kept':>10}  "
          f"{'inf g_alt on [wbar_min,5]':>25}  {'worst cell':>20}")
    for mt in (0.0,) + MUTILDES + MUTILDES_EXT:
        rows = []
        for i in range(len(grid) - 1):
            a, b = float(grid[i]), float(grid[i + 1])
            mu2, k1, k2, k3, case = mu2_pencil_cell(a, b, mt)
            g = float(mu2 - als[i])
            ksig = bool(k1 >= 0 and k2 >= 0 and k3 >= 0)
            rows.append((a, b, float(mu2), float(als[i]), g, ksig, case))
        # largest suffix of cells on which BOTH g>0 and the k-signs hold
        start = None
        for i in range(len(rows)):
            if all(r[4] > 0 and r[5] for r in rows[i:]):
                start = i
                break
        chk("a nonempty certified suffix exists", start is not None, f"mt={mt}")
        kept = rows[start:]
        ginf = min(r[4] for r in kept)
        worst = min(kept, key=lambda r: r[4])
        G_INT[mt] = (kept[0][0], ginf, len(kept))
        chk("interval infimum positive", ginf > 0, f"mt={mt} {ginf:.6f}")
        for r in kept:
            chk("cell g_alt > 0 on the certified suffix", r[4] > 0,
                f"mt={mt} cell [{r[0]:.4f},{r[1]:.4f}]")
            chk("cell mu_2 monotonicity k-signs hold", r[5],
                f"mt={mt} cell [{r[0]:.4f},{r[1]:.4f}]")
        print(f"  {mt:>8.3f}  {kept[0][0]:>9.4f}  {len(kept):>4}/{len(rows):<5}  "
              f"{ginf:>25.6f}  [{worst[0]:.4f},{worst[1]:.4f}]")

        # SOUNDNESS of the cell bound on the cells we KEEP.  The k-sign test is
        # what makes `mu2_pencil_cell` a valid segment bound, and an adversarial
        # sweep showed the bound really does fail without it (see that function's
        # docstring).  So verify directly that on every RETAINED cell the bound
        # sits below the true 2x2 top eigenvalue at interior points and at both
        # endpoints, for every window case.  (Evidence, float64: it audits the
        # certificate rather than entering it.)
        for r in kept[::7]:
            a, b = r[0], r[1]
            for wq in (a, 0.5 * (a + b), b):
                for Nc in (0, 1, 2, 9, 300):
                    v = float(np.linalg.eigvalsh(M_mat(2, wq, 0, mt, Nc))[-1])
                    chk("retained cell: mu_2 bound <= true lam_1(M^(2)(w))",
                        v >= r[2] - 1e-13,
                        f"mt={mt} cell [{a:.4f},{b:.4f}] w={wq:.4f} Nc={Nc}: "
                        f"true={v:.9f} bound={r[2]:.9f}")
        # and that the DISCARDED prefix really is discarded for a reason
        if start > 0:
            drop = rows[:start]
            chk("the discarded prefix fails g_alt>0 or the k-signs (not both ok)",
                any((not r[5]) or r[4] <= 0 for r in drop), f"mt={mt}")
        if mt == 0.0:
            chk("mutilde=0 recovers the FULL interval [0.3,5]",
                abs(kept[0][0] - 0.3) < 1e-12 and len(kept) == len(rows),
                f"start={kept[0][0]}")
            chk("mutilde=0 recovers the published infimum 0.093225",
                abs(ginf - 0.093225) < 1e-5, f"{ginf:.6f}")
    # COVERING LEMMA consistency: wbar_min must be nondecreasing in mutilde
    # (a node's interval is inherited by every smaller mutilde, so a larger
    # mutilde can never certify a LARGER interval).
    mts = (0.0,) + MUTILDES + MUTILDES_EXT
    for m1, m2 in zip(mts[:-1], mts[1:]):
        chk("covering lemma: wbar_min nondecreasing in mutilde",
            G_INT[m1][0] <= G_INT[m2][0] + 1e-12, f"{m1}->{m2}")
    chk("lambda<1 corner node mutilde=0.26 certifies [wbar_min,5] with wbar_min<1.17",
        G_INT[0.26][0] < 1.17, f"wbar_min={G_INT[0.26][0]:.4f}")
    # the k-signs never bind: the first DROPPED cell fails on g<=0 with k-signs ok
    print("  ==> mu*=0 row reproduces the v0.6 theorem (inf gap >= 0.093225 on the")
    print("      whole [0.3,5]); the mutilde>0 rows are the honest sub-intervals,")
    print("      each inherited by every smaller mutilde (covering lemma).")
    _summary("Stage E", f"({time.time()-t0:.1f}s)")


# ============================================================================
# STAGE F -- THE N < N_c SEGMENT (task item (iii)).
#
# The note's Sec. 2.3 flags "N < N_c needs separate treatment: the truncation
# cuts into the Coulomb window and UNDER-estimates the repulsion, so KAY Thm 4's
# variational ladder needs N >= N_c".  The resolution, and it is exact:
#
#   P_N K_c(N_c) P_N = K_c(min(N, N_c))            (asserted below, entrywise)
#
# so for N >= N_c the compression is the identity map on K_c (the note's clean
# case), and for N < N_c it shrinks the window.  EITHER WAY,
#
#   M_N := P_N M P_N  is the leading principal N x N block of the SINGLE fixed
#   self-adjoint operator M = H - mutilde K_c,
#
# because neither H's nor K_c's entry formulas depend on the truncation size.
# Cauchy interlacing on nested principal blocks therefore applies for EVERY
# N >= 1, with no N >= N_c proviso, and gives
#
#   lam_1(M_N) NONDECREASING in N, lam_1(M_N) -> lam_1(M),                  (F)
#
# i.e. the variational ladder is intact through the window boundary.  The reason
# the "under-estimated repulsion" objection does not bite is that the truncation
# under-estimates the ATTRACTION by exactly as much: it is the ASSEMBLED pencil,
# not K_ph and K_c separately, that the interlacing theorem is applied to.  What
# would break is a ladder applied to the two pieces separately, or a claim that
# the truncated problem is "Eliashberg with cutoff w_c" -- for N < N_c it is
# instead the same theory with the SMALLER cutoff 2 pi T N.  That reading is
# stated, and the monotone ladder is verified numerically ACROSS N = N_c below.
#
# Consequence for the certificates: nothing in Stage C/D/E needs an N >= N_c
# restriction.  Stage D's infimum sits at N = 2, the deepest possible truncation,
# far INSIDE the Coulomb window for every physical N_c -- and it is certified
# there in closed form.
# ============================================================================

def run_stage_F():
    print("\n=== STAGE F: the N < N_c segment -- compression identity + ladder "
          "across the window ===")
    t0 = time.time()
    # (a) the compression identity, entrywise exact
    for Nc in (1, 2, 3, 7, 16, 33, 128):
        for N in (1, 2, 3, 5, 8, 17, 40, 129, 300):
            lhs = K_c(N, Nc)                              # = P_N K_c(N_c) P_N
            rhs = K_c(N, min(N, Nc))
            chk("P_N K_c(N_c) P_N == K_c(min(N,N_c)) exactly",
                np.max(np.abs(lhs - rhs)) == 0.0, f"N={N} Nc={Nc}")
            big = K_c(max(N, 4 * Nc), Nc)[:N, :N]
            chk("K_c has truncation-independent entries (leading principal block)",
                np.max(np.abs(big - lhs)) == 0.0, f"N={N} Nc={Nc}")
    print("  (a) P_N K_c(N_c) P_N = K_c(min(N,N_c)) verified entrywise exactly:")
    print("      N >= N_c is the identity (the note's clean case); N < N_c merely")
    print("      shrinks the window to N, i.e. the SAME theory at cutoff 2 pi T N.")

    # (b) the monotone ladder ACROSS N = N_c, which is where the objection lived
    print("  (b) Cauchy ladder lam_1(M_N) nondecreasing in N, checked THROUGH the")
    print("      window boundary N = N_c (marked *):")
    for mt in MUTILDES:
        for w in (0.577, 1.0, 2.0, 5.0):
            for Nc in (4, 16, 64):
                prev = -np.inf
                seq = []
                for N in range(1, 4 * Nc + 1):
                    l1 = float(np.linalg.eigvalsh(M_mat(N, w, 0, mt, Nc))[-1])
                    chk("lam_1(M_N) nondecreasing in N across N=N_c (Cauchy)",
                        l1 >= prev - 1e-13,
                        f"w={w} mt={mt} Nc={Nc} N={N}: {l1} < {prev}")
                    prev = l1
                    if N in (1, 2, Nc - 1, Nc, Nc + 1, 2 * Nc, 4 * Nc):
                        seq.append((N, l1))
                l_inf = float(np.linalg.eigvalsh(M_mat(8 * Nc, w, 0, mt, Nc))[-1])
                chk("the ladder converges from BELOW to the limit lam_1(M)",
                    prev <= l_inf + 1e-12, f"w={w} mt={mt} Nc={Nc}")
                if (w, Nc) == (1.0, 16):
                    print(f"      wbar={w} mutilde={mt} N_c={Nc}:  " +
                          "  ".join(f"N={N}{'*' if N == Nc else ''}:{v:.6f}"
                                    for N, v in seq))
    # (c) the interlacing bound is BINDING at N=2, deep inside the window.
    #     g_alt is recomputed here from `mu2_pencil_cell` rather than read out of
    #     Stage D's table, so Stage F can be run on its own.
    for mt in MUTILDES:
        for w in WBARS:
            if w < 1.0:
                continue
            g = float(mu2_pencil_cell(w, w, mt)[0] - alpha_certificate_H(w)[0])
            for Nc in (2, 8, 64, 400):
                for N in (2, 3, 5, 9, 33, 129):
                    E = np.linalg.eigvalsh(M_mat(N, w, 0, mt, Nc))
                    chk("g_alt bounds gap_N for N < N_c as well as N >= N_c",
                        E[-1] - E[-2] >= g - 1e-12,
                        f"w={w} mt={mt} N={N} Nc={Nc}: gap={E[-1]-E[-2]:.6f} "
                        f"g_alt={g:.6f}")
    print("  (c) g_alt bounds gap_N at every N tested, INCLUDING N < N_c: the")
    print("      infimum sits at N=2 (deepest truncation, deep inside every")
    print("      physical window) and is certified there in closed form.")
    print("  ==> NO N >= N_c proviso is needed anywhere in Stages C/D/E.")
    _summary("Stage F", f"({time.time()-t0:.1f}s)")


# ============================================================================
# STAGE G -- PHYSICS ANCHORS (task item (iii)): no regression after embedding.
#
# G1 (mu* = 0, Marsiglio-Carbotte eq. 94).  The Tc solver built on the PENCIL
#    must give Tc/(sqrt(lambda) w_E) -> sqrt(frak_g(2))/2pi = 0.18272 = M-C's
#    0.183, confirming that the Tc-condition code path (not just the operator) is
#    right.
# G2 (mu* != 0, TOLMACHEV).  At fixed Tc, scanning the window N_c and solving for
#    the mu* that keeps Tc fixed must give the parameter-free Morel-Anderson-
#    Tolmachev slope  d(1/mu*)/d ln w_c = -1.  Slope -1 is locked jointly by the
#    prefactor 2 in 2|u><u| and the 1/sqrt(2n+1) weight (prefactor 1 -> slope -2).
#    Closed form, no bisection: 1 is an eigenvalue of lambda H - 2 mu* |u><u| iff
#        1 = 2 mu* u^T (lambda H - I)^{-1} u                                  (S)
#    -- one linear solve per window; every solution is verified to be the TOP
#    eigenvalue, so no interior crossing is mistaken for the Tc condition.
# ============================================================================

_HC = {}


def _H(N, w):
    key = (N, float(w))
    if key not in _HC:
        if len(_HC) > 4:
            _HC.clear()
        _HC[key] = H_mat(N, w)
    return _HC[key]


def solve_Tc(lam, N=800):
    """mu*=0 Tc: solve rho(lambda H(1/(2 pi T))) = 1.  w_E = 1.  rho is decreasing
    in T, so bisection is safe."""
    def f(T):
        return lam * float(np.linalg.eigvalsh(_H(N, 1.0 / (2 * np.pi * T)))[-1]) - 1.0
    lo, hi = 1e-3, 1e4
    assert f(lo) > 0 > f(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-12 * hi:
            break
    return 0.5 * (lo + hi)


def mustar_fixed_T(N, w, lam, Nc):
    """The mu* that puts the TOP eigenvalue of lambda H - mu* K_c at 1, via (S)."""
    B = lam * _H(N, w) - np.eye(N)
    u = u_vec(N, Nc)
    d = 2.0 * float(u @ np.linalg.solve(B, u))
    return 1.0 / d if d > 0 else None


def run_stage_G():
    print("\n=== STAGE G: physics anchors -- no regression after embedding K_c ===")
    t0 = time.time()
    # G1: strong-coupling limit at mu* = 0
    N = 800
    gN = float(np.linalg.eigvalsh(piu.G_matrix(N, 2))[-1])
    target = np.sqrt(gN) / (2 * np.pi)
    print(f"  G1 (mu*=0): Tc/(sqrt(lam) w_E) -> sqrt(g^({N})(2))/2pi = "
          f"{target:.7f}   [M-C eq.94: 0.183]")
    print(f"     {'lambda':>9}  {'Tc/w_E':>12}  {'Tc/(sqrt(lam) w_E)':>19}")
    prev = None
    for lam in (2.0, 10.0, 100.0, 1000.0, 10000.0):
        Tc = solve_Tc(lam, N)
        r = Tc / np.sqrt(lam)
        chk("Tc/sqrt(lam) stays below the KAY constant", r < target + 1e-9,
            f"lam={lam}")
        if prev is not None:
            chk("Tc/sqrt(lam) increases toward the constant", r > prev, f"lam={lam}")
        prev = r
        print(f"     {lam:>9.0f}  {Tc:>12.6f}  {r:>19.7f}")
    chk("strong-coupling limit reproduces sqrt(g(2))/2pi", abs(prev - target) < 2e-4,
        f"{prev:.7f} vs {target:.7f}")
    chk("strong-coupling limit reproduces M-C eq.94 constant 0.183",
        abs(prev - 0.183) < 5e-4, f"{prev:.7f}")
    print(f"     => {prev:.7f}  (M-C eq.94: 0.183)  OK -- the Tc code path on the")
    print(f"        pencil is unchanged at mu*=0.")

    # G2: Tolmachev logarithm with the Coulomb term switched ON
    print("  G2 (mu*!=0): Tolmachev slope d(1/mu*)/d ln w_c -> -1 "
          "(w_c = 2 pi Tc N_c at fixed Tc)")
    for wbar, Ncs, fac in ((2.0, [16, 32, 64, 128, 256, 512], 4),
                           (6.0, [32, 64, 128, 256, 512], 4)):
        Nc_ref, = (Ncs[len(Ncs) // 2],)
        Nref = fac * Nc_ref
        lam0 = 1.0 / float(np.linalg.eigvalsh(_H(Nref, wbar))[-1])
        lo, hi, f = lam0 * (1 + 1e-9), None, 1.02
        while f < 4.0:
            m = mustar_fixed_T(Nref, wbar, lam0 * f, Nc_ref)
            if m is not None and m >= 0.1:
                hi = lam0 * f
                break
            if m is not None:
                lo = lam0 * f
            f *= 1.02
        chk("lambda bracket calibrating mu*(N_c_ref)=0.1 found", hi is not None,
            f"wbar={wbar}")
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            m = mustar_fixed_T(Nref, wbar, mid, Nc_ref)
            if m is None or m < 0.1:
                lo = mid
            else:
                hi = mid
        lam = hi
        chk("lambda calibration lands on mu* = 0.100",
            abs(mustar_fixed_T(Nref, wbar, lam, Nc_ref) - 0.1) < 1e-6, f"wbar={wbar}")
        print(f"     --- varpi={wbar}  (Tc/w_E={1/(2*np.pi*wbar):.5f}), "
              f"lambda={lam:.6f} calibrated so mu*(N_c={Nc_ref})=0.100 ---")
        print(f"     {'N_c':>6}  {'w_c/w_E':>8}  {'mu*':>10}  {'1/mu*':>9}")
        xs, ys = [], []
        for Nc in Ncs:
            Nt = fac * Nc
            m = mustar_fixed_T(Nt, wbar, lam, Nc)
            chk("secular solution exists", m is not None, f"wbar={wbar} Nc={Nc}")
            rho = float(np.linalg.eigvalsh(lam * _H(Nt, wbar)
                                           - m * K_c(Nt, Nc))[-1])
            chk("the secular mu* puts the TOP eigenvalue exactly at 1",
                abs(rho - 1.0) < 1e-9, f"Nc={Nc} rho={rho:.12f}")
            m2 = mustar_fixed_T(int(1.5 * Nt), wbar, lam, Nc)
            chk("mu* converged in the truncation N",
                abs(m2 - m) < 2e-3 * m, f"Nc={Nc}: {m:.8f} vs {m2:.8f}")
            xs.append(np.log(Nc))
            ys.append(1.0 / m)
            print(f"     {Nc:>6}  {2*np.pi*Nc/(2*np.pi*wbar):>8.2f}  {m:>10.6f}  "
                  f"{1/m:>9.4f}")
        xs, ys = np.array(xs), np.array(ys)
        loc = np.diff(ys) / np.diff(xs)
        A = np.vstack([xs[-3:], np.ones(3)]).T
        s_hi = float(np.linalg.lstsq(A, ys[-3:], rcond=None)[0][0])
        print(f"     local slopes: " + " ".join(f"{v:+.5f}" for v in loc))
        print(f"     fit over the largest 3 windows: {s_hi:+.5f}   [prediction -1]")
        chk("Tolmachev slope -> -1 at the largest window",
            abs(loc[-1] + 1.0) < 0.01, f"wbar={wbar} {loc[-1]:.5f}")
        chk("Tolmachev slope -> -1 (fit over the largest 3 windows)",
            abs(s_hi + 1.0) < 0.02, f"wbar={wbar} {s_hi:.5f}")
        chk("local slopes converge MONOTONICALLY to -1",
            all(abs(loc[i + 1] + 1.0) < abs(loc[i] + 1.0) for i in range(len(loc) - 1)),
            f"wbar={wbar} {loc}")
        chk("the prefactor-1 alternative (slope -2) is EXCLUDED",
            abs(s_hi + 2.0) > 0.9, f"wbar={wbar} {s_hi:.5f}")
    print("     => the embedding is anchor-clean: mu*=0 strong coupling unchanged,")
    print("        mu*!=0 Tolmachev logarithm reproduced to <1e-2.")
    _summary("Stage G", f"({time.time()-t0:.1f}s)")


# ============================================================================

def run_verdict():
    print("\n=== VERDICT ===")
    print("  (i) GATE A: at mu*=0 the pencil reproduces the 1668-check certificate")
    print("      digit-for-digit (H_mat rebound to the pencil builder, whole file")
    print("      re-run, printed output diffed line by line).")
    print(f"  (ii) pencil gap certificates ({'mutilde = mu*/lambda' })")
    print(f"      {'mutilde':>8}  {'delta(wbar=1)':>13}  {'g_alt(wbar=1)':>13}  "
          f"{'interval [wbar_min,5]':>21}  {'inf g_alt':>10}")
    for mt in (0.0,) + MUTILDES + MUTILDES_EXT:
        wmin, ginf, ncell = G_INT[mt]
        d1 = (f"{DELTA_P[(1.0, mt)][0]:>13.6f}" if (1.0, mt) in DELTA_P
              else f"{'n/a (wbar=1 < wbar_min)':>13}")
        g1 = (f"{GALT_P[(1.0, mt)][0]:>13.6f}" if (1.0, mt) in GALT_P
              else f"{'n/a':>13}")
        print(f"      {mt:>8.3f}  {d1}  {g1}  "
              f"{'[%.4f, 5]' % wmin:>21}  {ginf:>10.6f}")
    print("      (each row holds verbatim for every smaller mutilde: covering lemma)")
    print("  (iii) N < N_c: no separate treatment NEEDED.  P_N K_c(N_c) P_N =")
    print("      K_c(min(N,N_c)) exactly, so M_N is a leading principal block of the")
    print("      one fixed operator M and Cauchy interlacing holds for every N.  The")
    print("      'under-estimated repulsion' is exactly offset by the equally")
    print("      under-estimated attraction; what N < N_c changes is the READING of")
    print("      the truncated model (cutoff 2 pi T N, not w_c), not the ladder.")
    print("  (iv) anchors: M-C eq.94 constant 0.183 (mu*=0) and the Tolmachev slope")
    print("      -1 (mu*!=0) both reproduced through the pencil code path.")
    print("  (v) HONEST BOUNDARY: taking the infimum over ALL windows N_c, the gap")
    print("      genuinely collapses at the small-wbar end (lam_1(M) ~ 1e-4 at")
    print("      wbar=0.3, N_c=3, mutilde=0.1), so the interval theorem is certified")
    print("      on a sub-interval only.  The mu*=0 theorem on all of [0.3,5] with")
    print("      inf gap >= 0.093225 is recovered exactly as mutilde -> 0.")


if __name__ == "__main__":
    want = [s.upper() for s in sys.argv[1:]] or ["A", "B", "C", "D", "E", "F", "G"]
    T0 = time.time()
    if "A" in want:
        run_stage_A()
    if "B" in want:
        run_stage_B()
    if "C" in want:
        run_stage_C()
    if "D" in want:
        run_stage_D()
    if "E" in want:
        run_stage_E()
    if "F" in want:
        run_stage_F()
    if "G" in want:
        run_stage_G()
    if set(want) >= {"C", "D", "E"}:
        run_verdict()
    print(f"\nALL {PASS[0]} checks passed  ({time.time()-T0:.1f}s).")
