r"""
P4 uniqueness x P4-D Coulomb pencil: is there a mu*>0 monotonicity statement?

CONTEXT (read `P4_uniqueness_strategy.md` and `P4D_coulomb_operator.md` first).
The mu*=0 proof shows T -> lam1(K(P,T)) strictly decreasing via a cone C =
{v >= 0 nonincreasing} that the top eigenvector of K(P,T) always sits in
(Lemma A/B, semigroup cone invariance).  The Coulomb pencil

    M(wbar, mutilde, N_c) = H(wbar) - mutilde * K_c(N_c),   K_c = 2|u><u|,
    u_n = (2n+1)^{-1/2} 1_{n<N_c},   mutilde = mu*/lambda in [0, 0.33],

kicks the top eigenvector out of C in 58/75 sampled cases (`uniq_probe3.py`),
so the cone-invariance route does not transfer as is.  This script asks, with
fresh numerics (N=600, float64, code independent of `uniq_probe3.py` and of
`p4d_coulomb_certify.py` -- it only IMPORTS `H_mat`, never their logic):

  (1a) fixed window N_c: is wbar -> lam1(M) strictly increasing, for all
       N_c in a grid subset of {0..300} and mutilde in
       {0.05, 0.1, 0.13, 0.26, 0.33}?
  (1b) physical stepping N_c(wbar) = floor(r*wbar) or floor(r*wbar+1/2),
       r in {3, 10, 30}: where are the decreasing steps, how big, and do they
       cause multiple crossings of 1/lambda for lambda in [0.3, 5]?
  (1c) multi-atom P (2-5 Einstein atoms): at fixed N_c, is
       s -> lam1(sum_i p_i H(c_i s) - mutilde K_c) increasing (s = 1/T-like
       scale; increasing s = decreasing T)?
  (2)  secular-equation structure: derive lam1(M) as the largest root of
           1 = 2 mutilde <u, (H-lam)^{-1} u>            for lam in (lam2(H), lam1(H))
       (rank-one downward perturbation interlacing forces the root into this
       interval); use it as a FAST reusable computational engine (one eigh(H)
       per wbar serves every (mutilde, N_c) pair via a cached cumulative
       projection of u), and test whether the pencil's top eigenvector lies in
       the window-shifted cone C' = {v: v_n >= 0 nonincreasing for n >= N_c};
       then the Hellmann-Feynman sign test dlam1/dwbar = <v, H'(wbar) v>, which
       is what actually carries the fixed-window monotonicity (see FINDINGS 2).

RUN:  python theory/derivations/p4_uniqueness_coulomb_probe.py [stage ...]
      stages: 0 1a 1b 1c 2   (default: all)

Only import from certified files (`p4_interchange_uniformity.H_mat`); no
existing script is modified.

============================================================================
FINDINGS
============================================================================
N=600 float64 unless noted; logs and scratch probes in `theory/notes/_scratch_ir/`.

(1a) FIXED WINDOW: MONOTONE EXCEPT ONE SLIVER.  The first pass used an N_c grid
     that jumped 3 -> 5 -> 7 -> 10 and reported zero violations over 67830
     comparisons.  That was a sampling artifact: N_c=4 was never tested, and it
     is the one place monotonicity fails.  Sweeping ALL N_c in 1..40 against
     mutilde in {0.05, 0.1, 0.13, 0.2, 0.26, 0.33} (240 pairs) gives exactly ONE
     non-monotone pair:
         N_c = 4, mutilde in ~[0.322, 0.334], wbar in [0.71, 0.77],
     where wbar -> lam1(M) dips by 1.60e-6 (local max 3.00732e-4 at wbar=0.7075,
     local min 2.99138e-4 at wbar=0.7775) and lam1 = 1/lambda therefore has
     3 ROOTS -- explicitly counted -- but only for
         lambda in (3325.2, 3342.9).
     Monotone again at mutilde <= 0.320 and at mutilde >= 0.335, and for every
     other N_c.  So the failure is real but confined to lam1 ~ 3e-4, i.e. three
     orders of magnitude past the physical lambda <= 5; it lives in the regime
     where the Coulomb term has all but annihilated the top eigenvalue.
     Everywhere else in the box, lam1 = 1/lambda has AT MOST ONE root.

(2)  WHY (1a) nearly holds -- and why it is not a one-liner.  K_c is
     wbar-independent, so Hellmann-Feynman gives dlam1/dwbar = <v, H'(wbar) v>,
     which was positive in every sampled case EXCEPT the N_c=4 sliver above, and
     matches a central difference to rel 3e-7 there (so the negative sign is not
     float noise -- it matters, the raw value is -3.25e-5).  Both easy routes to
     the positive cases are closed:
       * H'(wbar) is strongly INDEFINITE -- lam_min/lam_max(H') ~ -1 at every
         wbar, N-independent to 6 digits (e.g. -0.514 / +0.828 at wbar=1).  So
         positivity is NOT "H' >= 0"; v must be localized somehow.
       * The cone C = {v >= 0 nonincreasing} that carried the mu*=0 proof FAILS
         for the pencil (267/320 cases), as does the window-shifted C' (204/320).
     Margin behaviour: restricted to mutilde <= 0.32 the normalized margin
     <v,H'v>/||H'|| stays positive, with the thinnest values at small N_c; in the
     wbar -> 0 corner it has a floor ~4.50e-4 (N-converged: 4.23e-4 at N=100 ->
     4.4975e-4 at N=600, attained at N_c=3).  The raw margin collapses to
     ~2.4e-6 there only because ||H'|| ~ wbar.  At mutilde=0.33 the margin
     crosses zero at N_c=4.
     => Any proof must either restrict mutilde below ~0.32 or carry the
     lambda restriction; and it still needs a localization of v on which H' is
     positive, since the old cone is not it.  That is the open gap.

(1b) STEPPED WINDOW N_c(wbar) = floor(r*wbar): monotonicity can FAIL, and
     whether that costs uniqueness depends on r.
       r=3  : all 24 N_c-steps are DOWNWARD (drops to -1.22e-1 for mutilde
              <= 0.13), and this DOES break uniqueness -- the exact root count
              finds 2 ROOTS at mutilde=0.10, lambda=1.0.
       r=10 : only 3 of 80 steps downward, all at wbar <= 0.1.  r=30: 3 of 240,
              all at wbar <= 0.033.  No multiple roots at any lambda in [0.5,3].
     The "0 roots at lambda=1.0, mutilde=0.13, r=10" line in
     uniq_coulomb_crossings.log is a RED HERRING for uniqueness: the step at
     wbar=2.2 is UPWARD (+1.77e-3), so lam1 jumps OVER the level
     (b_21=0.99838 < 1 < a_22=1.00015).  The level is never attained, but lam1
     is still increasing -- Tc is then a sup, not a second root.  The original
     classifier reported this identically to the dangerous DOWNWARD skips; the
     two are now distinguished.  Gap width ~ 1/r, so it closes as r -> infinity.

     r-THRESHOLD (log `uniq_coulomb_rscan.log`, script
     `_scratch_ir/p4_uniq_coulomb_rscan.py`, N=600, 1042 eigh, 71 s on
     petlove-train).  The gap between "r=3 broken" and "r>=10 clean" is closed
     by scanning r in {3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10} against the full
     physical grid lambda in {0.5, 0.7, 1.0, 1.5, 2.0} x mutilde in
     {0.05, 0.10, 0.13, 0.26} (20 combinations each), counting crossings per
     plateau with 4 interior samples plus the open right limit, so an
     in-plateau dip would be caught too (none occurs: every plateau is
     increasing at every scanned r, so ALL multiplicity here comes from the
     stepping, not from the fixed-window function):

         r = 3    : 2 ROOTS at mutilde=0.10, lambda in {1.0, 1.5}   -- FAILS
         r >= 3.5 : exactly one T_c for all 20 combinations, at every scanned r.

     So the threshold is r = 3.5 among the scanned values; the true threshold
     lies in (3, 3.5].  For a quotable hypothesis, r >= 4 is the safe integer
     statement.  Caveats: verified on this 20-point grid, wbar in [0.02, 8],
     N=600, and only at the listed r -- it is numerical evidence for the
     window-regularity hypothesis, not a proof of it.

     Two kinds of "0 crossings" appear and NEITHER is a uniqueness failure:
       * SUP (r = 4.5, 5, 7, 8, 9, 10, one (lambda,mutilde) each): an UPWARD gap
         steps over the level, exactly the r=10 red herring above.  lam1 is
         still globally increasing, {T : lambda*lam1 >= 1} is still an interval,
         so Tc = sup of that set is well defined and unique.  These are
         scattered and r-dependent because the gap positions move with r.
       * RANGE (mutilde=0.26, lambda=0.5, at EVERY r including r=3): the level
         1/lambda = 2.0 simply is not reached anywhere on wbar <= 8.  This is
         window truncation, not stepping: lam1 keeps rising and crosses 2.0 near
         wbar ~ 15 (1.516 at wbar=8 -> 2.018 at wbar=16, r=3).  Being
         r-independent it carries no information about the threshold.
     Because RANGE fires at every r, a "strictly one attained crossing at every
     (lambda,mutilde)" threshold does not exist on this window -- that is an
     artifact of wbar <= 8 and of Tc-as-sup, and hides no multiplicity.

(1c) MULTI-ATOM P: monotonicity is FALSE; uniqueness survives at physical lambda.
     One of 5 sampled P's is non-monotone in the scale s:
         P = 0.8 delta(w - 0.3 s) + 0.2 delta(w - 4.0 s),
     and only for N_c in [5, 12] AND mutilde >= ~0.17 (monotone at 0.16, fails
     at 0.18): s -> lam1(M) rises to s~0.45, FALLS on (0.45, 1.33), then rises.
     This is real, not an engine artifact -- the secular engine agrees with a
     dense eigvalsh(M) to 1e-16 and lam1(H)-lam2(H) = 0.41 (no degeneracy).
     HF pins the cause (FD agrees to 9 digits): at s=0.81 the c=4.0 component
     contributes -9.66e-2 against the c=0.3 component's +5.04e-2.  A mixture
     forces one v to sit in the NEGATIVE cone of H'(c_2 s) while sitting in the
     positive cone of H'(c_1 s) -- something a single atom never has to do.
     The dip does yield 3 ROOTS, but only for
         lambda in (24.3, 75.0) at mutilde=0.26,  (26.0, 1738) at mutilde=0.33,
     i.e. only where the level 1/lambda < ~0.041 -- FAR outside the physical
     lambda in [0.3, 5] (level >= 0.2).  The other 4 P's (2-atom symmetric,
     3-atom, 4/5-atom random, c-spread up to 15.4) are monotone in all 15
     (N_c, mutilde) series each.

NET EVIDENCE ON UNIQUENESS (for Section 3)
  "T -> lam1 decreasing" is FALSE in general once the Coulomb term is on.  Three
  independent failures were found, and they are NOT on the same footing --
    (i)   fixed N_c: one sliver, N_c=4 and mutilde ~0.325, giving 3 roots at
          lambda ~ 3.3e3  -- unphysical by three orders of magnitude;
    (ii)  asymmetric 2-atom P with mutilde >= ~0.17: 3 roots at
          lambda in (24, 75)  -- unphysical by a factor ~5;
    (iii) COARSE WINDOW STEPPING r=3: 2 roots at lambda = 1.0, mutilde = 0.10
          -- squarely INSIDE the physical range.  This one is not a curiosity.
          It is also SHARP in r: the r-scan above gives uniqueness across the
          whole physical grid at every scanned r >= 3.5, so r=3 is the only
          scanned window coarse enough to break it.
  So the honest statement is conditional, not absolute:
    * At fixed N_c, or with a sufficiently fine window (r >= 3.5 on the scanned
      grid; safe integer form r >= 4), uniqueness of the root of lam1 = 1/lambda
      is SUPPORTED throughout mutilde in [0, 0.33], lambda in [0.3, 5].
    * It is NOT supported for a coarsely stepped window.  Uniqueness therefore
      requires a regularity hypothesis on N_c(wbar) -- a smooth/continuous
      window, or r bounded below -- and this hypothesis is doing real work, not
      bookkeeping.  Section 3 should state it rather than inherit monotonicity.
  The cone/monotonicity route that settled mu*=0 does not extend as is: the cone
  is not invariant for the pencil (267/320), and H' is indefinite, so even the
  cases that DO work are unproved -- only verified.

SCOPE REDUCTIONS (shared machine, BLAS pinned to 1 thread)
  Stage 1c runs at N=600 with ds=0.04.  The mutilde-threshold, N_c-range and
  lambda-window refinements quoted in (1c) were run at N=400 with ds=0.01 and
  0.004 (`probe_1c_window.py`), cross-checked against N=600/ds=0.01 on the
  (N_c=10, mutilde=0.26) series (`probe_1c_multiroot.py`).  Stage-2 cone counts
  use the 8-point wbar grid; the HF sweep uses a 14-point grid extended down to
  wbar=0.002 (thin-margin corner) and including wbar=0.74 (the sign failure).
  The dense N_c=1..40 sweep behind (1a), the dip resolution and the mutilde
  threshold are in `probe_nc4_dip.py`, `probe_nc4_extent.py`, `probe_hf_thin.py`.
  No claim here rests on a single grid resolution; the two that matter most --
  the r=3 two-root case and the N_c=4 dip -- were both confirmed by dense
  eigvalsh(M) and by finite differences at three step sizes, independently of
  the secular engine.

N=1200 CONFIRMATION (petlove-train, `p4_uniq_coulomb_recheck_N1200.py`, log
  `uniq_coulomb_recheck_N1200.log`).  Both load-bearing counterexamples are
  N-stable, so neither is a truncation artefact:
    * N_c=4 dip: depth 1.5946e-6 (N=600) -> 1.5917e-6 (N=1200), a 0.2% change;
      local max/min move by < 1e-8; the 3-root lambda window is
      (3325.22, 3342.94) -> (3325.16, 3342.85); <v,H'v> at wbar=0.74 is
      -3.2488e-5 -> -3.2448e-5 (FD-confirmed at h=1e-4,1e-5,1e-6 at both N);
      the mutilde window is identical, [0.325, 0.330], monotone at 0.320 and
      0.335.  The near-degeneracy that made this suspect (gap of M ~3e-4) is
      itself stable, 3.06e-4 -> 3.01e-4 -- it is a genuine feature, not a
      symptom of under-resolution.
    * r=3 two-root case: plateau endpoints agree to 9-10 digits between N=600
      and N=1200; the downward jump at wbar=2.0 is -8.514512e-3 at both.

CAVEAT ON NEGATIVE RESULTS
  "No multiple root found" is grid evidence, not a proof: every failure above was
  found only after refining a grid that had previously looked clean (N_c=4 was
  missed by a grid that skipped it; the thin HF margin was missed by a coarse
  wbar grid).  Treat the supported cases as unfalsified, not established.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p4_interchange_uniformity import H_mat        # noqa: E402  (certified builder, read-only)

PASS = [0]
_LAST = [0]


def chk(name, ok, detail=""):
    if not ok:
        raise AssertionError(f"CHECK FAILED: {name} {detail}")
    PASS[0] += 1


def _summary(tag, note=""):
    k = PASS[0] - _LAST[0]
    _LAST[0] = PASS[0]
    print(f"  [{tag}] {k} checks passed. {note}")


# ============================================================================
# THE PENCIL, BUILT INDEPENDENTLY OF p4d_coulomb_certify.py
# ============================================================================

def u_full_vec(N):
    n = np.arange(N)
    return 1.0 / np.sqrt(2.0 * n + 1.0)


def K_c_local(N, Nc):
    u = np.zeros(N)
    k = max(0, min(int(Nc), N))
    u[:k] = 1.0 / np.sqrt(2.0 * np.arange(k) + 1.0)
    return 2.0 * np.outer(u, u)


def M_direct(N, w, mutilde, Nc):
    """Direct dense construction, used ONLY for the Stage-0 cross-check and
    small spot probes (never inside the fast scan)."""
    return H_mat(N, w) - mutilde * K_c_local(N, Nc)


# ---- fast engine: one eigh(H(wbar)) serves every (mutilde, N_c) -----------

def eig_H(N, w):
    """Ascending eigendecomposition of H(wbar) plus the cumulative projection
    cumU[k-1, :] = sum_{n<k} u_full[n] * V[n, :] = <e_i, (I placeholder)>,
    i.e. cumU[k-1, :] is V^T u(N_c=k) for EVERY k at once."""
    Hm = H_mat(N, w)
    evals, V = np.linalg.eigh(Hm)              # ascending
    uf = u_full_vec(N)
    cumU = np.cumsum(V * uf[:, None], axis=0)  # cumU[k-1,:] = V^T u(Nc=k)
    return evals, V, cumU


def lam1_secular_batch(evals, cumU, Nc_list, mutilde, iters=60):
    """lam1(M) for a BATCH of N_c values at fixed mutilde, given the cached
    (evals, cumU) of H(wbar).  Returns an array aligned with Nc_list.

    Derivation (see module docstring / the .md write-up): write
    M = H - mutilde*2|u><u|.  If (H-lam)v = 2*mutilde*c*u with c = <u,v> != 0,
    then v = 2*mutilde*c*(H-lam)^{-1}u, and taking <u,.> gives the SECULAR
    EQUATION  1 = 2*mutilde*<u,(H-lam)^{-1}u> = 2*mutilde*sum_i uhat_i^2/(lam_i-lam).
    As lam -> lam1(H)^- the i=top term -> +inf; as lam -> lam2(H)^+ the
    i=second term -> -inf; the sum is STRICTLY increasing in lam wherever
    finite (derivative = 2*mutilde*sum uhat_i^2/(lam_i-lam)^2 > 0) so there is
    EXACTLY one root in (lam2(H), lam1(H)), found here by bisection.
    Rank-one downward-perturbation interlacing (Weyl) guarantees
    lam2(H) <= lam1(M) <= lam1(H), i.e. the root IS lam1(M)."""
    Nc_arr = np.asarray(Nc_list, dtype=int)
    N = evals.size
    out = np.empty(Nc_arr.shape, dtype=float)
    zero_mask = (Nc_arr <= 0) | (mutilde == 0.0)
    out[zero_mask] = evals[-1]
    nz = ~zero_mask
    if not np.any(nz):
        return out
    idx = np.clip(Nc_arr[nz], 1, N) - 1
    Uhat = cumU[idx, :]                      # (m, N)
    lam1H, lam2H = evals[-1], evals[-2]
    span = lam1H - lam2H
    eps = 1e-13 * max(1.0, abs(span))
    lo = np.full(Uhat.shape[0], lam2H + eps)
    hi = np.full(Uhat.shape[0], lam1H - eps)

    def g(lam_vec):
        return 2.0 * mutilde * np.sum(Uhat**2 / (evals[None, :] - lam_vec[:, None]), axis=1)

    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        gm = g(mid)
        below = gm < 1.0
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)
    out[nz] = 0.5 * (lo + hi)
    return out


def lam1_secular_one(N, w, mutilde, Nc):
    evals, V, cumU = eig_H(N, w)
    return float(lam1_secular_batch(evals, cumU, [Nc], mutilde)[0])


# ============================================================================
# STAGE 0: self-check of the secular-equation engine against eigvalsh(M)
# ============================================================================

def run_stage0():
    print("=== STAGE 0: secular-equation engine self-check ===")
    t0 = time.time()
    N = 300
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(40):
        w = float(rng.uniform(0.1, 6.0))
        mutilde = float(rng.choice([0.05, 0.1, 0.13, 0.26, 0.33]))
        Nc = int(rng.integers(0, 250))
        lam1_direct = float(np.linalg.eigvalsh(M_direct(N, w, mutilde, Nc))[-1])
        lam1_fast = lam1_secular_one(N, w, mutilde, Nc)
        d = abs(lam1_direct - lam1_fast)
        worst = max(worst, d)
        chk("secular engine == eigvalsh(M)", d < 1e-8,
            f"w={w:.4f} mt={mutilde} Nc={Nc} direct={lam1_direct:.10f} fast={lam1_fast:.10f} d={d:.2e}")
    print(f"  40 random (wbar,mutilde,N_c), N=300: max|secular-eigvalsh| = {worst:.2e}")

    # sign check: interlacing lam2(H) <= lam1(M) <= lam1(H)
    for _ in range(20):
        w = float(rng.uniform(0.1, 6.0))
        mutilde = float(rng.uniform(0.01, 0.4))
        Nc = int(rng.integers(1, 250))
        evals, V, cumU = eig_H(N, w)
        lam1M = float(lam1_secular_batch(evals, cumU, [Nc], mutilde)[0])
        chk("lam2(H) <= lam1(M) <= lam1(H)", evals[-2] - 1e-9 <= lam1M <= evals[-1] + 1e-9,
            f"w={w} mt={mutilde} Nc={Nc} lam2H={evals[-2]} lam1M={lam1M} lam1H={evals[-1]}")
    _summary("Stage 0", f"({time.time()-t0:.1f}s)")


# ============================================================================
# STAGE 1a: fixed-window monotonicity in wbar
# ============================================================================

N_MAIN = 600
# N_c grid.  4, 6, 8, 9, 12, 25, 40 were ADDED after the first pass: the original
# grid jumped 3 -> 5 -> 7 -> 10, and N_c=4 is exactly where the only fixed-window
# monotonicity failure lives (see EXPECTED_1A_FAILURES and FINDINGS 1a).
NC_FIXED = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250, 300]
MU_LIST = [0.05, 0.1, 0.13, 0.26, 0.33]
R_LIST = [3, 10, 30]
WMAX = 8.0
DW = 0.01

# The one (N_c, mutilde) pair at which wbar -> lam1(M) is NOT increasing.  Pinned
# as an exact expectation so the probe is a regression test in both directions:
# a new failure OR the disappearance of this one both trip the check.
EXPECTED_1A_FAILURES = {(4, 0.33)}


def _wgrid():
    return np.arange(0.02, WMAX + 1e-9, DW)


def run_stage1a_and_1b():
    """Single pass over a common wbar grid; produces BOTH the fixed-N_c
    series (1a) and the stepped-N_c series (1b), reusing one eigh(H(wbar))
    per grid point for everything (the expensive part)."""
    print("=== STAGE 1a/1b: monotonicity scan, N=%d, dw=%.3f, wbar in [0.02,%.1f] ===" % (N_MAIN, DW, WMAX))
    t0 = time.time()
    ws = _wgrid()
    nW = ws.size
    fixed = {(Nc, mt): np.empty(nW) for Nc in NC_FIXED for mt in MU_LIST}
    variants = {}
    for r in R_LIST:
        variants[("floor", r)] = np.empty(nW, dtype=int)
        variants[("round", r)] = np.empty(nW, dtype=int)
    stepped = {(v, r, mt): np.empty(nW) for v in ("floor", "round") for r in R_LIST for mt in MU_LIST}

    for i, w in enumerate(ws):
        evals, V, cumU = eig_H(N_MAIN, w)
        for r in R_LIST:
            variants[("floor", r)][i] = int(np.floor(r * w))
            variants[("round", r)][i] = int(np.floor(r * w + 0.5))
        for mt in MU_LIST:
            lam = lam1_secular_batch(evals, cumU, NC_FIXED, mt)
            for j, Nc in enumerate(NC_FIXED):
                fixed[(Nc, mt)][i] = lam[j]
            for v in ("floor", "round"):
                Ncs = [variants[(v, r)][i] for r in R_LIST]
                lamv = lam1_secular_batch(evals, cumU, Ncs, mt)
                for j, r in enumerate(R_LIST):
                    stepped[(v, r, mt)][i] = lamv[j]
        if i % 200 == 0:
            print(f"    ... wbar={w:.3f} ({i+1}/{nW}, {time.time()-t0:.0f}s)")

    print(f"  grid computed in {time.time()-t0:.0f}s")

    # ---- 1a report: strictly increasing at fixed N_c? ----
    print("\n  -- 1a: fixed-N_c monotonicity --")
    worst_decr = {}
    found = set()
    for (Nc, mt), arr in fixed.items():
        d = np.diff(arr)
        neg = d < -1e-11
        if np.any(neg):
            found.add((Nc, mt))
            # lam1 levels ON the decreasing steps -- these, not the global range
            # of lam1, are the levels whose root count the dip can corrupt.
            lv = np.concatenate([arr[:-1][neg], arr[1:][neg]])
            worst_decr[(Nc, mt)] = (float(d[neg].min()), float(ws[:-1][neg][np.argmin(d[neg])]),
                                    float(lv.min()), float(lv.max()))
    print(f"  {len(fixed)} (N_c, mutilde) combos tested x {nW-1} consecutive pairs each"
          f" = {(len(fixed))*(nW-1)} monotonicity comparisons")
    print(f"  non-monotone combos (decreasing step beyond tol 1e-11): {len(found)}")
    for key, (dmin, wat, lmin, lmax) in sorted(worst_decr.items()):
        print(f"    Nc={key[0]} mt={key[1]}: worst decrease {dmin:.3e} at wbar~{wat:.3f};"
              f" lam1 on the decreasing steps in [{lmin:.6e}, {lmax:.6e}]"
              f" => only levels with lambda ~ {1.0/lmax:.0f}..{1.0/lmin:.0f} are affected")
    # Monotonicity is NOT universal; assert the exact failure set instead, so that
    # both a NEW failure and the disappearance of the known one are caught.
    chk("fixed-N_c monotonicity fails exactly on EXPECTED_1A_FAILURES",
        found == EXPECTED_1A_FAILURES,
        f"found={sorted(found)} expected={sorted(EXPECTED_1A_FAILURES)}")
    for key in EXPECTED_1A_FAILURES:
        _, _, lmin, lmax = worst_decr[key]
        chk(f"known failure (Nc={key[0]}, mt={key[1]}) only affects unphysically large lambda",
            1.0 / lmax > 100.0,
            f"lam1 up to {lmax:.6e} on the decreasing steps => lambda down to"
            f" {1.0/lmax:.1f}, not > 100")

    return ws, fixed, variants, stepped


def analyze_stage1b(ws, variants, stepped):
    print("\n  -- 1b: physical stepping N_c(wbar) = floor(r*wbar) or floor(r*wbar+1/2) --")
    lam_grid = np.array([0.3, 0.5, 0.7, 1.0, 1.3, 1.5, 2.0, 3.0, 5.0])
    for v in ("floor", "round"):
        for r in R_LIST:
            Ncs = variants[(v, r)]
            for mt in MU_LIST:
                arr = stepped[(v, r, mt)]
                d = np.diff(arr)
                dNc = np.diff(Ncs)
                neg = d < -1e-10
                n_neg = int(np.sum(neg))
                worst = float(d.min()) if n_neg else 0.0
                worst_w = float(ws[:-1][np.argmin(d)]) if n_neg else float("nan")
                # crossing counts of level 1/lambda
                crossings = {}
                for lam in lam_grid:
                    level = 1.0 / lam
                    s = np.sign(arr - level)
                    s[s == 0] = 1.0
                    ncross = int(np.sum(np.diff(s) != 0))
                    crossings[lam] = ncross
                cr_str = " ".join(f"{lam}:{crossings[lam]}" for lam in lam_grid)
                print(f"    variant={v:5s} r={r:2d} mt={mt:5.2f}: decreasing steps on grid = {n_neg:4d}"
                      f"  worst drop {worst:+.3e} at wbar~{worst_w:.3f}  |  crossings[lambda]: {cr_str}")


# exact +-eps jump sizes at the true jump locations (small supplementary stage)

def run_stage1b_exact(r, variant, mutilde_list, wmax=8.0, N=N_MAIN, eps=1e-6):
    print(f"  -- exact jump sizes, r={r}, variant={variant} --")
    t0 = time.time()
    jmax = int(np.floor(r * wmax)) + 1
    results = {mt: [] for mt in mutilde_list}
    for j in range(1, jmax + 1):
        wj = (j / r) if variant == "floor" else ((j - 0.5) / r)
        if wj <= 0.02 or wj > wmax:
            continue
        w_below, w_above = wj - eps, wj + eps
        Nc_below = int(np.floor(r * w_below)) if variant == "floor" else int(np.floor(r * w_below + 0.5))
        Nc_above = int(np.floor(r * w_above)) if variant == "floor" else int(np.floor(r * w_above + 0.5))
        if Nc_above <= Nc_below:
            continue
        evals_b, V_b, cumU_b = eig_H(N, w_below)
        evals_a, V_a, cumU_a = eig_H(N, w_above)
        for mt in mutilde_list:
            fb = float(lam1_secular_batch(evals_b, cumU_b, [Nc_below], mt)[0])
            fa = float(lam1_secular_batch(evals_a, cumU_a, [Nc_above], mt)[0])
            results[mt].append((wj, Nc_below, Nc_above, fa - fb))
    for mt in mutilde_list:
        drops = [x for x in results[mt] if x[3] < 0]
        print(f"    mt={mt}: {len(results[mt])} jumps checked, {len(drops)} DECREASING"
              f" (largest drop {min((x[3] for x in drops), default=0.0):.3e}"
              f" at wbar={min(drops, key=lambda x: x[3])[0]:.4f})" if drops else
              f"    mt={mt}: {len(results[mt])} jumps checked, 0 decreasing")
    print(f"    ({time.time()-t0:.0f}s)")
    return results


# ---- EXACT root count of lam1 = 1/lambda along the physical path -----------
# Stage 1a says lam1(., N_c) is continuous and strictly increasing on each
# N_c-constant interval, so on interval j the attained set is exactly the open
# range (a_j, b_j) of its endpoint values and a level is hit AT MOST ONCE
# there.  Hence  #roots = #{j : level in (a_j, b_j)}  is EXACT -- independent of
# any wbar grid (the dw=0.01 scan above can only lower-bound it).  A level that
# lies in no range but between b_j and a_{j+1} is SKIPPED by a downward jump:
# the Tc equation then has no solution and Tc must be read as a sup.

def run_stage1b_crossings(r, variant, mutilde_list, lam_levels,
                          wlo=0.02, whi=8.0, N=N_MAIN, eps=1e-7):
    print(f"  -- EXACT root count of lam1=1/lambda, r={r}, variant={variant} --")
    t0 = time.time()
    jmax = int(np.floor(r * whi)) + 2
    # ranges[mt][j] = (wl, wr, a_j, b_j) for the N_c=j plateau
    ranges = {mt: [] for mt in mutilde_list}
    for j in range(0, jmax + 1):
        if variant == "floor":
            wl, wr = j / r, (j + 1) / r
        else:
            wl, wr = (j - 0.5) / r, (j + 0.5) / r
        wl, wr = max(wl, wlo), min(wr, whi)
        if wr - wl <= 2 * eps:
            continue
        ev_l, _, cu_l = eig_H(N, wl + eps)
        ev_r, _, cu_r = eig_H(N, wr - eps)
        for mt in mutilde_list:
            a = float(lam1_secular_batch(ev_l, cu_l, [j], mt)[0])
            b = float(lam1_secular_batch(ev_r, cu_r, [j], mt)[0])
            chk(f"plateau j={j} r={r} {variant} mt={mt}: increasing across plateau",
                b > a - 1e-12, f"a={a:.12f} b={b:.12f}")
            ranges[mt].append((wl, wr, a, b, j))
    for mt in mutilde_list:
        rg = ranges[mt]
        for lam in lam_levels:
            level = 1.0 / lam
            hits = [(wl, wr, jj) for (wl, wr, a, b, jj) in rg if a < level < b]
            # A level can be missed by the range of lam1 in TWO inequivalent ways
            # at an N_c step, and the distinction decides uniqueness:
            #   up_gap   b_j < level < a_{j+1}: lam1 jumps UP over the level.  The
            #            level is simply not attained -- lam1 is still globally
            #            increasing, Tc is a sup, uniqueness is NOT threatened.
            #   down_jmp a_{j+1} < level < b_j: lam1 jumps DOWN over the level.
            #            This is the dangerous one: it can hand the level back a
            #            second time and create extra roots.
            up_gap, down_jmp = [], []
            for (p, q) in zip(rg[:-1], rg[1:]):
                if p[3] < level < q[2]:
                    up_gap.append((p[4], q[4], p[1]))
                if q[2] < level < p[3]:
                    down_jmp.append((p[4], q[4], p[1]))
            flag = ""
            if len(hits) == 0:
                flag = ("  <== NO ROOT (skipped by UPWARD gap; lam1 still increasing,"
                        " Tc is a sup)" if up_gap and not down_jmp else "  <== NO ROOT")
            elif len(hits) > 1:
                flag = f"  <== {len(hits)} ROOTS (Tc not unique)"
            sk = (f"  DOWNWARD-jump skip at wbar={[f'{x[2]:.4f}' for x in down_jmp]}"
                  if down_jmp else "")
            ug = (f"  upward-gap at wbar={[f'{x[2]:.4f}' for x in up_gap]}"
                  if up_gap and len(hits) == 0 else "")
            print(f"    mt={mt:5.2f} lambda={lam:4.1f} (level 1/lam={level:.4f}):"
                  f" {len(hits)} root(s) at wbar in {[f'({a:.4f},{b:.4f})' for a, b, _ in hits][:6]}"
                  f"{flag}{sk}{ug}")
    print(f"    ({time.time()-t0:.0f}s)")
    return ranges


# ============================================================================
# STAGE 1c: multi-atom P, T-monotonicity at fixed N_c
# ============================================================================

def H_multi(N, s, comps, order=0):
    """H for P = sum_i p_i delta(w - c_i*<w>); comps = list of (p_i, c_i)."""
    out = np.zeros((N, N))
    for p, c in comps:
        out += p * H_mat(N, c * s, order)
    return out


def run_stage1c():
    print("=== STAGE 1c: multi-atom P, fixed N_c, T-monotonicity ===")
    t0 = time.time()
    N = N_MAIN
    rng = np.random.default_rng(7)
    Ps = [
        ("2-atom symmetric", [(0.5, 0.4), (0.5, 2.2)]),
        ("2-atom asymmetric", [(0.8, 0.3), (0.2, 4.0)]),
        ("3-atom", [(0.3, 0.3), (0.4, 1.0), (0.3, 3.5)]),
    ]
    for k in (4, 5):
        cs = np.exp(rng.uniform(np.log(0.15), np.log(6.0), k))
        ps = rng.dirichlet(0.8 * np.ones(k))
        cs = cs / np.sum(ps * cs)     # normalize so <w>=1 (mean frequency 1)
        Ps.append((f"{k}-atom random", list(zip(ps.tolist(), cs.tolist()))))

    ss = np.arange(0.05, 6.0 + 1e-9, 0.04)
    Nc_list = [0, 3, 10, 50, 150]
    mt_list = [0.1, 0.13, 0.26]
    n_violations = 0
    n_series = 0
    violations = []
    uf = u_full_vec(N)
    for name, comps in Ps:
        # ONE eigh per (P, s): the cached (evals, cumU) serve every (N_c, mutilde)
        lam = {(Nc, mt): np.empty(ss.size) for Nc in Nc_list for mt in mt_list}
        for i, s in enumerate(ss):
            Hm = H_multi(N, s, comps)
            evals, V = np.linalg.eigh(Hm)
            cumU = np.cumsum(V * uf[:, None], axis=0)
            for mt in mt_list:
                vals = lam1_secular_batch(evals, cumU, Nc_list, mt)
                for j, Nc in enumerate(Nc_list):
                    lam[(Nc, mt)][i] = vals[j]
        for Nc in Nc_list:
            for mt in mt_list:
                n_series += 1
                arr = lam[(Nc, mt)]
                d = np.diff(arr)
                neg = d < -1e-10
                if np.any(neg):
                    n_violations += 1
                    # shape of the failure: how many sign changes does the
                    # derivative have, and over which levels 1/lambda does the
                    # resulting dip hand back multiple roots?
                    sgn = np.sign(d)
                    sgn[sgn == 0] = 1.0
                    flips = np.where(np.diff(sgn) != 0)[0] + 1
                    msg = (f"    NON-MONOTONE: {name} Nc={Nc} mt={mt}:"
                           f" worst drop {d.min():.3e} at s~{ss[:-1][np.argmin(d)]:.3f},"
                           f" {int(neg.sum())}/{d.size} decreasing steps")
                    if len(flips) >= 2:
                        L_hi, L_lo = float(arr[flips[0]]), float(arr[flips[1]])
                        if L_lo < L_hi:
                            msg += (f"\n      local max lam1={L_hi:.6f} @s={ss[flips[0]]:.3f},"
                                    f" local min lam1={L_lo:.6f} @s={ss[flips[1]]:.3f}"
                                    f" => 3 roots for lambda in ({1.0/L_hi:.2f}, {1.0/L_lo:.2f})")
                    violations.append(msg)
                    print(msg)
    print(f"  {n_series} (P, N_c, mutilde) series x {ss.size-1} steps each tested;"
          f" NON-MONOTONE series: {n_violations}")
    # Engine-level invariant that MUST hold even where monotonicity fails:
    # interlacing lam2(H) <= lam1(M) <= lam1(H) (checked inside the loop above
    # would cost another eigh; assert the weaker finiteness/ordering here).
    chk("stage 1c produced finite lam1 everywhere",
        all(np.all(np.isfinite(v)) for v in lam.values()))
    _summary("Stage 1c", f"({time.time()-t0:.0f}s); monotonicity is NOT universal -- see above")
    return violations


# ============================================================================
# STAGE 2: cone-membership probe (C and the window-shifted C')
# ============================================================================

def run_stage2_cone():
    print("=== STAGE 2: cone membership of the pencil's top eigenvector ===")
    t0 = time.time()
    N = N_MAIN
    ws = [0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
    mts = [0.05, 0.1, 0.13, 0.26, 0.33]
    Ncs = [1, 3, 5, 10, 20, 50, 100, 200]
    tot = 0
    vio_C = 0
    vio_Cprime = 0
    vio_pos = 0
    examples = []
    for w in ws:
        Hm = H_mat(N, w)
        for Nc in Ncs:
            Kc = K_c_local(N, Nc)
            for mt in mts:
                tot += 1
                M = Hm - mt * Kc
                evals, V = np.linalg.eigh(M)
                v = V[:, -1] * np.sign(V[np.argmax(np.abs(V[:, -1])), -1])
                pos_ok = v.min() > -1e-10
                C_ok = pos_ok and np.all(np.diff(v) <= 1e-11)
                tail = v[Nc:] if Nc < N else np.array([])
                Cprime_ok = (tail.size < 2) or (
                    np.all(tail >= -1e-10) and np.all(np.diff(tail) <= 1e-11))
                if not pos_ok:
                    vio_pos += 1
                if not C_ok:
                    vio_C += 1
                if not Cprime_ok:
                    vio_Cprime += 1
                    if len(examples) < 6:
                        examples.append((w, Nc, mt))
    print(f"  {tot} (wbar,N_c,mutilde) cases (N={N}):")
    print(f"    violates v>=0 : {vio_pos}/{tot}")
    print(f"    violates cone C (v>=0 nonincreasing everywhere): {vio_C}/{tot}")
    print(f"    violates window-shifted C' (v nonincr. for n>=N_c only): {vio_Cprime}/{tot}")
    if examples:
        print(f"    C' counterexamples (wbar,Nc,mt): {examples}")
    _summary("Stage 2 (cone probe)", f"({time.time()-t0:.0f}s)")

    # ---- Hellmann-Feynman: dlam1/dwbar = <v, H'(wbar) v> at FIXED N_c --------
    # This is the positive result of the whole script.  Two facts make it
    # non-trivial rather than a one-liner:
    #   (i)  H'(wbar) is strongly INDEFINITE -- lam_min(H') ~ -lam_max(H') at
    #        every wbar and every N tested -- so "<v,H'v> > 0 for all v" is
    #        FALSE and some localization of v is genuinely needed;
    #   (ii) the cone C (and C') that carried the mu*=0 proof fails for the
    #        pencil (counts above), so that localization is not the old one.
    # Nevertheless <v,H'v> > 0 held in every case sampled.  The margin is
    # reported NORMALIZED by ||H'||, because ||H'|| ~ wbar as wbar -> 0 and the
    # raw number therefore collapses for trivial scaling reasons; the
    # normalized margin has a floor ~4.5e-4 attained at N_c=3, mutilde=0.33.
    print("  -- Hellmann-Feynman sign check <v, H'(wbar) v> at the pencil top eigenvector --")
    # 0.74 is in the grid deliberately: it is the centre of the ONE region where
    # the sign fails (N_c=4, mutilde~0.325-0.33).  Without it the sweep reports
    # uniform positivity and the probe would assert something false.
    ws_hf = [0.002, 0.005, 0.02, 0.1, 0.3, 0.5, 0.7, 0.74, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0]
    Ncs_hf = [1, 2, 3, 4, 5, 7, 10, 20, 50, 100, 200, 300]
    hf_neg_keys = []
    hf_tot = 0
    worst_abs = (np.inf, None)
    worst_nrm = (np.inf, None)
    worst_nrm_safe = (np.inf, None)          # restricted to mutilde <= 0.32
    Hp_indef = 0
    for w in ws_hf:
        Hm = H_mat(N, w)
        Hp = H_mat(N, w, order=1)
        eHp = np.linalg.eigvalsh(Hp)
        nrmHp = float(np.max(np.abs(eHp)))
        if eHp[0] < 0:
            Hp_indef += 1
        for Nc in Ncs_hf:
            Kc = K_c_local(N, Nc)
            for mt in mts:
                hf_tot += 1
                evals, V = np.linalg.eigh(Hm - mt * Kc)
                v = V[:, -1]
                hf = float(v @ Hp @ v)
                if hf <= 0:
                    hf_neg_keys.append((w, Nc, mt, hf))
                if hf < worst_abs[0]:
                    worst_abs = (hf, (w, Nc, mt))
                if hf / nrmHp < worst_nrm[0]:
                    worst_nrm = (hf / nrmHp, (w, Nc, mt))
                if mt <= 0.32 and hf / nrmHp < worst_nrm_safe[0]:
                    worst_nrm_safe = (hf / nrmHp, (w, Nc, mt))
        chk(f"H'(wbar={w}) is indefinite (so <v,H'v>>0 needs v-localization)", eHp[0] < 0,
            f"lam_min(H')={eHp[0]:.6e} lam_max(H')={eHp[-1]:.6e}")
    print(f"    H'(wbar) indefinite (lam_min<0) at {Hp_indef}/{len(ws_hf)} sampled wbar"
          f"  ==> monotonicity does NOT follow from H'>=0")
    print(f"    {hf_tot} cases: <v,H'(wbar)v> <= 0 in {len(hf_neg_keys)} cases")
    for w, Nc, mt, hf in hf_neg_keys:
        print(f"      NEGATIVE: wbar={w} Nc={Nc} mt={mt}  <v,H'v>={hf:+.6e}")
    print(f"    worst raw margin      {worst_abs[0]:+.6e} at (wbar,Nc,mt)={worst_abs[1]}")
    print(f"    worst margin/||H'||   {worst_nrm[0]:+.6e} at (wbar,Nc,mt)={worst_nrm[1]}")
    print(f"    worst margin/||H'|| restricted to mutilde<=0.32:"
          f" {worst_nrm_safe[0]:+.6e} at (wbar,Nc,mt)={worst_nrm_safe[1]}")
    # The TRUE statements (a universal "hf > 0" is FALSE -- see FINDINGS 1a/2):
    chk("every HF sign failure has N_c=4 and mutilde=0.33",
        all(Nc == 4 and mt == 0.33 for _, Nc, mt, _ in hf_neg_keys),
        f"unexpected failures {[(w, Nc, mt) for w, Nc, mt, _ in hf_neg_keys if not (Nc == 4 and mt == 0.33)]}")
    chk("HF sign is positive throughout mutilde <= 0.32", worst_nrm_safe[0] > 0.0,
        f"worst {worst_nrm_safe[0]:.3e} at {worst_nrm_safe[1]}")
    chk("the documented N_c=4 failure is still reproduced at wbar=0.74",
        any(abs(w - 0.74) < 1e-12 and Nc == 4 and mt == 0.33 for w, Nc, mt, _ in hf_neg_keys),
        "expected a negative HF at (wbar=0.74, N_c=4, mutilde=0.33)")

    # The HF identity itself, cross-checked against a central difference so the
    # sign cannot be a float64 artifact (it matters: the margin is ~1e-5 here).
    w, Nc, mt = worst_nrm[1]
    h = 1e-6
    Kc = K_c_local(N, Nc)
    l_p = float(np.linalg.eigvalsh(H_mat(N, w + h) - mt * Kc)[-1])
    l_m = float(np.linalg.eigvalsh(H_mat(N, w - h) - mt * Kc)[-1])
    fd = (l_p - l_m) / (2.0 * h)
    evals, V = np.linalg.eigh(H_mat(N, w) - mt * Kc)
    hf = float(V[:, -1] @ H_mat(N, w, order=1) @ V[:, -1])
    rel = abs(fd - hf) / max(abs(hf), 1e-300)
    print(f"    HF identity at the worst point: <v,H'v>={hf:+.10e} vs central FD={fd:+.10e}"
          f" (rel {rel:.2e})")
    chk("HF identity matches central finite difference at the worst point", rel < 1e-5,
        f"hf={hf:.10e} fd={fd:.10e} rel={rel:.3e}")
    _summary("Stage 2 (HF sign)", f"({time.time()-t0:.0f}s)")


# ============================================================================
if __name__ == "__main__":
    args = sys.argv[1:] or ["0", "1a", "1b", "1c", "2"]
    T0 = time.time()
    if "0" in args:
        run_stage0()
    ws = fixed = variants = stepped = None
    if "1a" in args or "1b" in args:
        ws, fixed, variants, stepped = run_stage1a_and_1b()
    if "1b" in args:
        analyze_stage1b(ws, variants, stepped)
        run_stage1b_exact(3, "floor", [0.05, 0.1, 0.13, 0.26, 0.33])
        run_stage1b_exact(3, "round", [0.1, 0.13])
        run_stage1b_exact(10, "floor", [0.1, 0.13])
        run_stage1b_exact(30, "floor", [0.1, 0.13])
        LEVELS = [0.5, 0.7, 1.0, 1.3, 2.0, 3.0]
        run_stage1b_crossings(3, "floor", [0.1, 0.13, 0.26], LEVELS)
        run_stage1b_crossings(5, "floor", [0.1, 0.13], LEVELS)
        run_stage1b_crossings(10, "floor", [0.1, 0.13], LEVELS)
    if "1c" in args:
        run_stage1c()
    if "2" in args:
        run_stage2_cone()
    print(f"\nTOTAL: {PASS[0]} checks passed. ({time.time()-T0:.0f}s)")
