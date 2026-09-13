r"""
P4-D . UNIQUENESS OF Tc WITH A COULOMB PSEUDOPOTENTIAL (mu* > 0)
================================================================================

RUN:  python theory/derivations/p4d_coulomb_uniqueness.py           (all stages)
      python theory/derivations/p4d_coulomb_uniqueness.py S0 S2     (a subset)
      PN=500 python theory/derivations/p4d_coulomb_uniqueness.py    (override N)

Companion note: `P4D_coulomb_uniqueness.md`.
Base files (NOT modified, only imported): `p4_interchange_uniformity.py`
(H_mat, _brk_vec), `p4_uniqueness_lemmas.py` (the mu*=0 Lemmas A/B),
`p4d_coulomb_certify.py` (the pencil certificate machine).

--------------------------------------------------------------------------------
0. WHAT THIS FILE IS ABOUT
--------------------------------------------------------------------------------
`P4_uniqueness_strategy.md` proves KAY II Conjecture 2 (uniqueness of Tc) for
mu* = 0.  Its engine is a CONE argument:

    C = {v : v_n >= 0, v nonincreasing},   v = sum_L a_L e_L,  e_L = 1_{[0,L)},
    Lemma A  : Q_{L,L'}(delta_k) >= 0                 ==> <v,O[f]v> > 0 on C\{0}
    Lemma B  : W_{L,L'}(K) >= 0 for L != L'           ==> exp(t O[F]) C subset C
    ==> the Perron vector of K_ph lies in C ==> d lam_1/dT < 0 ==> Tc unique.

With mu* > 0 the pencil operator is

    M(wbar, mutilde, N_c) = A - mutilde * K_c,     A := K_ph(wbar) = H(wbar),
    K_c = 2|u_c><u_c|,  u_c[n] = (2n+1)^{-1/2} 1_{n<N_c},  N_c = floor(w_c/2piT),
    mutilde = mu*/lambda,   Tc boundary:  lambda * lam_1(M) = 1.

The Coulomb term is a NEGATIVE rank-one perturbation, so in the step (Abel)
basis it contributes ONE negative generator and the cone argument breaks: the
top eigenvector of M leaves C.  This file quantifies that breakdown and
replaces the cone argument by a RANK-ONE RESOLVENT / PERTURBATION argument that
survives it, with an explicit threshold mutilde_*(wbar, N_c).

MAIN RESULT PROVED HERE (analytic; constants float64, see HONEST BOUNDARIES)
--------------------------------------------------------------------------------
    THEOREM (within-window strict monotonicity).
    Fix N_c >= 1 and r = w_c/<w>, and let  W(N_c) := [N_c/r, (N_c+1)/r)  be the
    wbar-window on which N_c(wbar) = floor(r wbar) is constant.  For each wbar
    let mu_1 > mu_2 >= ... be the eigenvalues of A = K_ph(wbar) with orthonormal
    eigenvectors psi_j, and set

        gam  = mu_1 - mu_2                      (phonon spectral gap)
        c_1  = <u_c, psi_1>  ( > 0 ),   |u_perp| = ( |u_c|^2 - c_1^2 )^{1/2}
        q_1  = <psi_1, O[f] psi_1>  ( > 0  by Lemma A, since psi_1 in C )
        b    = || P_perp O[f] psi_1 ||
        m    = max(0, -lam_min(O[f]))
        t_max= ( -b + sqrt(b^2 + m q_1) ) / m           (= +inf if m = 0)
        mutilde_*(wbar,N_c) := t_max gam / ( 2 c_1 |u_perp| + 2 t_max c_1^2 )

    If
        mutilde  <  inf_{wbar in W(N_c)}  mutilde_*(wbar, N_c)

    then for EVERY wbar in W(N_c), <vhat, O[f] vhat> > 0 for the top eigenvector
    vhat of M, hence  T |--> lam_1(M(T))  is STRICTLY DECREASING on the whole
    N_c-window, hence  lambda * lam_1(M(T)) = 1  has AT MOST ONE root there.

    MIND THE QUANTIFIER.  Lemmas 1-4 are evaluated at ONE wbar and deliver the
    sign of T d_T lam_1 at THAT wbar only; window-wise monotonicity therefore
    needs the hypothesis at every wbar of the window, i.e. the INF over the
    window -- not the value at the nodal point wbar_c.  mutilde_* is NOT
    constant on a window: it dips by up to -9.2% below its wbar_c value
    (lambda = 2, r = 10 anchor).  Stage S3 tabulates the window inf next to the
    nodal value; all eleven anchors still pass.  (Reviewer note
    `theory/notes/2026-09-12_coulomb_uniqueness_review.md`, CHALLENGED 1.)

    Every physical anchor of `P4D_coulomb_operator.md` satisfies
    mutilde < inf_W mutilde_*  (Stage S3), including the lambda = 0.4, r = 3
    corner with mutilde = 0.325.

    The threshold correctly EXCLUDES the known counterexample
    (N_c = 4, wbar ~ 0.74, mutilde in [0.322, 0.334]):  mutilde_*(0.74,4)
    = 0.1249 < 0.322  (Stage S5).

WHY THE CONE ROUTE IS DEAD (Stage S1, honest negative result)
--------------------------------------------------------------------------------
The Kamke (quasi-monotonicity) coefficients of K_c have the exact closed form

    Wc_{L,L'} = 2 * S(min(L',N_c)) * ( u_{L-1} 1_{L-1<N_c} - u_L 1_{L<N_c} ),
    S(j) = sum_{i<j} u_i,   u_i = (2i+1)^{-1/2},      (L != L'),

which is >= 0 EVERYWHERE, with support exactly the rows L <= N_c (it vanishes
identically for L > N_c).  Because the Coulomb term enters with a MINUS sign,
M = A - mutilde K_c has Kamke coefficients  W^F - mutilde Wc,  so the single
positive generator Wc is exactly the obstruction: Lemma B survives iff

    mutilde <= mu_Metz(wbar,N_c) := min_{L != L', Wc > 0} W^F_{L,L'} / Wc_{L,L'}

(this is the precise, exact form of "tolerating one negative Abel coefficient":
the bad generator is the rank-one block u_c u_c^T and its Kamke footprint is
the closed form above, supported on N_c rows only)

and mu_Metz decays like  wbar^2 / N_c^3  (1.5e-3 at N_c=10, 2.7e-5 at N_c=40 for
wbar=1) -- two to four orders of magnitude BELOW the physical mutilde.  The
"generalized Lemma B tolerating one negative Abel coefficient" therefore exists
with an explicit epsilon, but the epsilon is physically useless.  The winning
route is NOT cone invariance.

THE BACKUP "OSCILLATION CONE" IS ALSO REFUTED (Stage S1c)
--------------------------------------------------------------------------------
The top eigenvector of M has EXACTLY TWO sign changes (+ for small n, - across
n ~ N_c, + in the tail), N-independent.  A Gantmacher-Krein cone of "at most
one sign change in consecutive differences" is the wrong object.

HONEST BOUNDARIES  (H1-H6; SAME numbering as `P4D_coulomb_uniqueness.md` Sec.7,
also restated verbatim by run_verdict() at the end of the run)
--------------------------------------------------------------------------------
(H1) gam, q_1, b, m, c_1, |u_perp| are computed in float64 `eigh` on the
     N-truncation.  The THEOREM is analytic given these six numbers; making it
     a certificate in this repo's sense needs interval enclosures for them
     (cf. the Tier1/Tier2 split of `p4_cii_m3_interval.py`).  The reported
     mutilde_* is an OVERESTIMATE: Stage S6 shows mutilde_*(N) decreasing
     MONOTONICALLY in N on all six cells tested, so the finite-N value is an
     upper bound for the N -> infinity limit.  Convergence is 6-7 digits only
     on the small cells; on the load-bearing anchors it is 2.1-3.5 digits
     (rel(N=400,N=1600) = 3.0e-4 at (13.5,135), 4.5e-4 at (15.6,46), 7.1e-3 at
     (51.2,153), the last still drifting, channel gam).  The margins survive
     (the binding anchor keeps 7.3%), but "6 significant digits from N = 200"
     is FALSE off the small cells.  Stage S7 spot-checks lam_1 with mpmath
     dps = 40 at N = 60.  Which constants the margin actually rests on is the
     subject of Stage S8: gam (elasticity exactly +1) and c_1 (-1.47..-1.72)
     are load-bearing, m is comfortably slack (flip needs +101%).
(H2) The theorem is a WITHIN-WINDOW statement (N_c fixed, hypothesis at the
     window INF).  Across N_c steps lam_1 jumps; Stage S4 gives a rigorous
     two-sided Rayleigh sandwich for the jump and explicit lambda-intervals
     where GLOBAL uniqueness genuinely FAILS (r = w_c/<w> = 3).  For r = 10 all
     steps on the true path are upward on [0.34, 6] (step 0.004, S4b) and on
     (6.0, 8.2] (step 0.02, S4d), so there global uniqueness holds FOR EVERY
     lambda WHOSE ROOT LIES IN THAT wbar-RANGE -- not for every lambda: below
     wbar = 0.34 there are three genuine downward steps (wbar = 0.105, 0.205,
     0.305; the first is N_c 0 -> 1, lam_1 dropping 1.3e-2 -> 1.0e-4), giving
     non-uniqueness for lambda in roughly (77, 1.1e4).  Unphysical, but real
     (Stage S4c).  This is not a gap in the theorem; it is a real feature of
     the stepped model, and Tc must be defined as a sup.
(H3) Between mutilde_* (provable) and mutilde_crit (the sharp numerical
     threshold where <vhat,O[f]vhat> first vanishes) there is a factor
     ~1.4-3 of unproved territory (the adversarial review measured
     min mutilde_crit/mutilde_* = 1.4025 over 224 (wbar,N_c) cells, at
     (wbar=20,N_c=50)).  Stage S2 tabulates both.  The known
     counterexample (N_c = 4, wbar ~ 0.74, mutilde in [0.322,0.334]) lies in
     that band and is a genuine failure of monotonicity -- so the band cannot
     be filled in trivially.  Correctly excluded by the theorem (Stage S5).
(H4) UPSTREAM DEPENDENCY, and one ingredient this file ADDS.  q_1 > 0 rests on
     Lemma A of `p4_uniqueness_lemmas.py` plus psi_1 in C, i.e. on the mu*=0
     theorem -- so this file is a genuine EXTENSION of
     `P4_uniqueness_strategy.md`, not a replacement.  BUT: that upstream proof
     states four times over that it needs NEITHER simplicity of the top
     eigenvalue NOR Krein-Rutman; it only proves E_1 (intersect) C != {0}.
     This file needs MORE: lam_1(A) SIMPLE (so that psi_1, q_1, and the
     Hellmann-Feynman step of Lemma 4 are well defined at all).  That extra
     ingredient is supplied here, not upstream, and it is supplied by exactly
     the positivity argument upstream avoided: every off-diagonal entry of
     A = K_ph(wbar) is u_n u_m [F(|n-m|) + F(n+m+1)] with F(k) = [[k]](wbar) > 0
     for every k >= 0, hence STRICTLY POSITIVE, hence Perron-Frobenius on the
     N-truncation gives lam_1(A) simple with psi_1 > 0 entrywise.  Recorded
     honestly: THIS THEOREM RE-INTRODUCES THE PERRON-FROBENIUS POSITIVITY THAT
     THE mu*=0 THEOREM DELIBERATELY AVOIDED.  Stage S0(e) verifies all the
     inputs.  (Reviewer note, CHALLENGED 3.)
(H5) N -> infinity tail.  All quantities are computed on P_N M P_N.  The mu*=0
     certificate machinery (`p4d_coulomb_certify.py`, alpha side + Poincare
     separation) controls the lam_2 tail, but the N -> infinity limits of q_1,
     b, m are not certified.  The Perron-Frobenius step of (H4) is also a
     finite-matrix statement; the operator version would need the compact
     Krein-Rutman analogue.
(H6) EINSTEIN-ONLY.  All six constants are evaluated on the SINGLE-ATOM
     EINSTEIN kernel A = K_ph(wbar) = H(wbar).  The mu*=0 theorem being
     extended holds for EVERY P in the class Pcal (<w^2> < infinity).  Lemmas
     1-4 only need A compact self-adjoint with simple top eigenvalue and psi_1
     in C, so the machine extends to multi-atom P -- but the constants have
     never been computed there, and the sibling probe
     `p4_uniqueness_coulomb_probe.py` (Stage 1c) finds a genuine WITHIN-WINDOW
     non-monotonicity at the PHYSICAL mutilde = 0.26 (lambda = 0.5,
     mu* = 0.13 -- the value of two anchors here) for the 2-atom
     P = 0.8 delta(w - 0.3 s) + 0.2 delta(w - 4.0 s) with N_c in [5,12]:
     lam_1 decreases by 1.8e-3 over s in (0.45,1.33).  So the restriction to
     Einstein is NOT cosmetic.  (Reviewer note, DOWNGRADE 6.)
"""

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p4_interchange_uniformity import H_mat, _brk_vec  # noqa: E402

PASS = [0]
N = int(os.environ.get("PN", "300"))


def ck(cond, msg):
    if not cond:
        raise AssertionError("FAILED: " + msg)
    PASS[0] += 1


def head(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ------------------------------------------------------------------ primitives
def u_vec(NN=None):
    NN = NN or N
    return 1.0 / np.sqrt(2.0 * np.arange(NN) + 1.0)


def uc_vec(Nc, NN=None):
    NN = NN or N
    u = u_vec(NN).copy()
    u[min(Nc, NN):] = 0.0
    return u


def K_c(Nc, NN=None):
    """K_c = 2 |u_c><u_c|  (the falsified -D variant is NOT included; see the note)."""
    u = uc_vec(Nc, NN)
    return 2.0 * np.outer(u, u)


def O_step(K, NN=None):
    """O[1_{k<=K}] on the NN x NN block; entries are N-independent."""
    NN = NN or N
    n = np.arange(NN)
    u = u_vec(NN)
    off = (np.abs(n[:, None] - n[None, :]) <= K).astype(float)
    np.fill_diagonal(off, 0.0)
    M = (off + ((n[:, None] + n[None, :] + 1) <= K)) * np.outer(u, u)
    M[np.diag_indices(NN)] -= 2.0 * np.minimum(n, K) / (2.0 * n + 1.0)
    return M


def kamke(M):
    """W[L-1, L'-1] = (M e_{L'})_{L-1} - (M e_{L'})_L,  e_{L'} = 1_{[0,L')}."""
    Y = np.cumsum(M, axis=1)
    return Y[:-1, :] - Y[1:, :]


def O_f(wbar, NN=None):
    """O[f] with f(k) = wbar * d/dwbar [[k]](wbar) >= 0;  -T d K/dT = O[f]."""
    return wbar * H_mat(NN or N, wbar, 1)


def split_L1_L23(wbar, NN=None):
    """O[f] = -L1 + (L2+L3), the decomposition of `P4_uniqueness_strategy.md` Sec.1."""
    NN = NN or N
    n = np.arange(NN)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)
    off = wbar * _brk_vec(np.abs(n[:, None] - n[None, :]), wbar, 1)
    np.fill_diagonal(off, 0.0)
    L23 = (off + wbar * _brk_vec(n[:, None] + n[None, :] + 1, wbar, 1)) * nrm
    fk = wbar * _brk_vec(np.arange(1, NN), wbar, 1)
    L1 = np.diag(np.concatenate(([0.0], np.cumsum(2.0 * fk))) / (2.0 * n + 1.0))
    return L1, L23


# --------------------------------------------------------------- the Cell class
class Cell:
    """All threshold ingredients for one (wbar, N_c) window."""

    def __init__(self, wbar, Nc, NN=None):
        NN = NN or N
        self.wbar, self.Nc, self.N = wbar, Nc, NN
        A = H_mat(NN, wbar, 0)
        Of = O_f(wbar, NN)
        e, V = np.linalg.eigh(A)
        self.e = e
        self.G = V.T @ Of @ V                      # O[f] in the A-eigenbasis
        self.c = V.T @ uc_vec(Nc, NN)              # <u_c, psi_j>
        self.mu1, self.mu2 = float(e[-1]), float(e[-2])
        self.gam = self.mu1 - self.mu2
        self.q1 = float(self.G[-1, -1])
        self.b = float(np.linalg.norm(self.G[:-1, -1]))
        self.m = max(0.0, -float(np.linalg.eigvalsh(Of)[0]))
        self.tmax = ((-self.b + np.sqrt(self.b ** 2 + self.m * self.q1)) / self.m
                     if self.m > 0 else np.inf)
        self.c1 = abs(float(self.c[-1]))
        self.up = float(np.linalg.norm(self.c[:-1]))
        self.thr = (self.tmax * self.gam
                    / (2 * self.c1 * self.up + 2 * self.tmax * self.c1 ** 2))

    # -- Lemma 1: secular equation --------------------------------------------
    def lam_of_mut(self, mut):
        """lam_1(M) = unique root in (mu_2, mu_1) of  2 mut sum_j c_j^2/(mu_j-lam) = 1."""
        c2 = self.c ** 2

        def f(l):
            return 2 * mut * float(np.sum(c2 / (self.e - l))) - 1.0

        sc = max(1.0, abs(self.mu1))
        lo, hi = self.mu2 + 1e-15 * sc, self.mu1 - 1e-15 * sc
        if f(lo) > 0:            # root pushed below mu_2 by a huge mut: dense fallback
            return float(np.linalg.eigvalsh(
                np.diag(self.e) - 2 * mut * np.outer(self.c, self.c))[-1])
        for _ in range(300):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    def vhat_coeffs(self, mut):
        """coefficients of vhat = (A-lam)^{-1} u_c in the A-eigenbasis, normalised."""
        lam = self.lam_of_mut(mut)
        y = self.c / (self.e - lam)
        nn = float(np.linalg.norm(y))
        return (y / nn if nn > 0 and np.isfinite(nn) else y), lam

    def q_of_mut(self, mut):
        """the monotonicity functional <vhat, O[f] vhat>  (>0 <=> d lam_1/dT < 0)."""
        y, lam = self.vhat_coeffs(mut)
        if not np.all(np.isfinite(y)) or np.linalg.norm(y) == 0:
            return np.nan
        return float(y @ (self.G @ y))

    def kappa(self, mut):
        return 2 * mut * self.c1 ** 2 / self.gam

    def tan_theta(self, mut):
        y, _ = self.vhat_coeffs(mut)
        c = abs(float(y[-1]))
        return np.sqrt(max(0.0, 1.0 - c * c)) / max(c, 1e-300)

    def mut_crit(self, hi=3.0):
        """sharp NUMERICAL threshold: the first mut>0 where q(vhat) <= 0."""
        grid = np.concatenate([np.linspace(1e-4, 0.6, 3000), np.linspace(0.6, hi, 1200)])
        prev = 1e-4
        for mut in grid:
            q = self.q_of_mut(mut)
            if not np.isfinite(q) or q <= 0:
                a_, b_ = prev, mut
                for _ in range(60):
                    mid = 0.5 * (a_ + b_)
                    qq = self.q_of_mut(mid)
                    if np.isfinite(qq) and qq > 0:
                        a_ = mid
                    else:
                        b_ = mid
                return a_
            prev = mut
        return np.inf


def top_pencil(A, uc, mut):
    e, V = np.linalg.eigh(A - 2.0 * mut * np.outer(uc, uc))
    v = V[:, -1]
    return float(e[-1]), float(e[-2]), v * (1.0 if v[0] >= 0 else -1.0)


# ============================================================== STAGE S0 (gate)
def run_S0():
    head("S0.  GATE + exact structure identities for K_c in the step (Abel) basis")

    # --- GATE: at mutilde = 0 we must reproduce the mu*=0 theory verbatim ------
    print("\n(a) GATE  mutilde = 0  reproduces `p4_uniqueness_lemmas.py`")
    NN = 200
    Om = O_step(10 ** 6, NN)                      # O[1]
    u = u_vec(NN)
    ck(np.abs(Om - (2 * np.outer(u, u) - np.eye(NN))).max() < 1e-13,
       "O[1] = 2 u u^T - I")
    print("    O[1] = 2 u u^T - I                                  max err "
          f"{np.abs(Om - (2*np.outer(u,u)-np.eye(NN))).max():.2e}")
    # Lemma B (Metzler) for the phonon kernel alone
    WF = kamke(H_mat(NN, 1.0, 0))
    L = np.arange(1, NN)[:, None]
    Lp = np.arange(1, NN + 1)[None, :]
    offd = (L != Lp)
    ck(WF[offd].min() > -1e-14, "Lemma B: W^F >= 0 off-diagonal (mutilde = 0)")
    print(f"    Lemma B off-diagonal min W^F (wbar=1)               {WF[offd].min():+.3e}")
    # Cor A': psi_1 in C and <psi_1, O[f] psi_1> > 0
    for wbar in (0.5, 1.0, 5.0, 20.0):
        C0 = Cell(wbar, 1, 200)
        e, V = np.linalg.eigh(H_mat(200, wbar, 0))
        p1 = V[:, -1] * (1.0 if V[0, -1] >= 0 else -1.0)
        ck(p1.min() > -1e-12 and (p1[:-1] - p1[1:]).min() > -1e-12,
           f"psi_1 in C at wbar={wbar}")
        ck(C0.q1 > 0, f"Cor A': q_1 > 0 at wbar={wbar}")
        print(f"    wbar={wbar:>5}: psi_1 in C = True,  q_1 = {C0.q1:.6e} > 0")

    # --- (b) the task brief's identification is NOT exact ----------------------
    print("\n(b) is  K_c == O[1_{k<N_c}] ?   (the natural guess -- it is FALSE)")
    for Nc in (3, 5, 10, 20):
        d = np.abs(K_c(Nc, NN) - O_step(Nc - 1, NN)).max()
        ck(d > 0.5, f"K_c != O[1_(k<Nc)] at Nc={Nc}")
        print(f"    Nc={Nc:>3d}:  max|K_c - O[1_(k<=Nc-1)]| = {d:.4f}   (NOT zero)")

    # --- (c) the exact identity ----------------------------------------------
    print("\n(c) EXACT:  K_c = P_{N_c} ( O[1] + I ) P_{N_c}")
    print("    PROOF (two lines, no numerics):  O[1] = 2 u u^T - I  (S0a), hence")
    print("    P(O[1]+I)P = 2 (Pu)(Pu)^T = 2 u_c u_c^T = K_c, because Pu = u_c by")
    print("    the definition u_c[n] = u_n 1_{n<N_c}.  The numbers below are")
    print("    BOOKKEEPING (they confirm the code, not the identity).")
    for Nc in (3, 5, 10, 20):
        P = np.zeros((NN, NN))
        P[np.arange(Nc), np.arange(Nc)] = 1.0
        d = np.abs(K_c(Nc, NN) - P @ (Om + np.eye(NN)) @ P).max()
        ck(d < 1e-13, f"K_c = P(O[1]+I)P at Nc={Nc}")
        print(f"    Nc={Nc:>3d}:  max|K_c - P(O[1]+I)P| = {d:.3e}")
        # the intermediate step of the proof, also as bookkeeping
        Pu = np.zeros(NN)
        Pu[:Nc] = u[:Nc]
        ck(np.abs(P @ (Om + np.eye(NN)) @ P - 2 * np.outer(Pu, Pu)).max() < 1e-13,
           f"P(O[1]+I)P = 2 (Pu)(Pu)^T at Nc={Nc}")

    # --- (d) closed form of the Kamke coefficients of K_c ---------------------
    print("\n(d) CLOSED FORM  Wc_{L,L'} = 2 S(min(L',N_c)) (u_{L-1}1_{L-1<N_c}"
          " - u_L 1_{L<N_c})")
    ue = u_vec(NN + 1)
    cs = np.concatenate([[0.0], np.cumsum(ue[:NN])])
    for Nc in (1, 2, 3, 5, 10, 37):
        W = kamke(K_c(Nc, NN))
        Sc = cs[np.minimum(Lp, Nc)]
        Wcf = 2.0 * Sc * (ue[L - 1] * (L - 1 < Nc) - ue[L] * (L < Nc))
        d = np.abs((W - Wcf)[offd]).max()
        pos = np.unique(np.where((W > 1e-15) & offd)[0]) + 1
        ck(d < 1e-13, f"Kamke closed form for K_c at Nc={Nc}")
        ck(W[offd].min() > -1e-15, f"Wc >= 0 off-diagonal at Nc={Nc}")
        ck(pos.size == 0 or pos.max() <= Nc,
           f"support of Wc is rows L <= N_c at Nc={Nc}")
        print(f"    Nc={Nc:>3d}: max|W-closed| = {d:.2e}  min Wc(offdiag) = "
              f"{W[offd].min():+.2e}  support rows L in [1,{pos.max() if pos.size else 0}]"
              f"  (<= N_c = {Nc})")
    print("\n    ==> Wc >= 0 everywhere, supported on the N_c rows L <= N_c.")
    print("        M = A - mutilde K_c therefore has Kamke matrix W^F - mutilde Wc:")
    print("        ONE positive generator block is subtracted, and that is the exact")
    print("        form of the 'one negative Abel coefficient' obstruction.")

    # --- (e) the ingredient this file ADDS to the upstream mu*=0 theorem ------
    print("\n(e) SIMPLICITY OF lam_1(A) -- the extra hypothesis this file needs (H4).")
    print("    Upstream (`P4_uniqueness_strategy.md`) proves only  E_1 (cap) C != {0}")
    print("    and states four times that it needs NEITHER simplicity of the top")
    print("    eigenvalue NOR Krein-Rutman.  Lemma 1 (simple root), q_1 = <psi_1,")
    print("    O[f] psi_1> and the Hellmann-Feynman step of Lemma 4 all need psi_1 to")
    print("    be a well-defined single vector, i.e. lam_1(A) SIMPLE.  Proof:")
    print("      A_{nm} = u_n u_m [ F(|n-m|) + F(n+m+1) ],  F(k) = [[k]](wbar)")
    print("             = wbar^2/(k^2+wbar^2) > 0  for every k >= 0,")
    print("      so every OFF-DIAGONAL entry of A is STRICTLY POSITIVE.  Then for any")
    print("      psi,  <|psi|,A|psi|> - <psi,A psi> = sum_{n!=m} A_{nm}(|psi_n psi_m|")
    print("      - psi_n psi_m) >= 0, the diagonal cancelling.  If psi is a top")
    print("      eigenvector then |psi| is too and equality holds, forcing psi not to")
    print("      change sign; a nonnegative top eigenvector has FULL support (if")
    print("      psi_n = 0 then 0 = sum_{m!=n} A_{nm} psi_m > 0, absurd); two")
    print("      orthogonal full-support nonnegative vectors cannot exist.  Hence")
    print("      lam_1(A) is SIMPLE and psi_1 > 0 entrywise.  This IS the")
    print("      Perron-Frobenius positivity upstream deliberately avoided --")
    print("      recorded as (H4).  Consequence (the short proof of c_1 > 0):")
    print("      psi_1 in C ==> psi_1 nonneg nonincreasing ==> psi_1[0] > 0, and")
    print("      u_c[0] = 1, so  c_1 = <u_c,psi_1> >= psi_1[0] > 0.")
    print(f"\n    {'wbar':>7} {'min offdiag A':>14} {'gam = mu_1-mu_2':>16} "
          f"{'min psi_1':>11} {'psi_1[0]':>9} {'c_1 (N_c=10)':>13}")
    for wbar in (0.05, 0.5, 1.0, 5.0, 13.5, 51.2):
        A = H_mat(NN, wbar, 0)
        offA = A.copy()
        np.fill_diagonal(offA, np.inf)
        ee, VV = np.linalg.eigh(A)
        p1 = VV[:, -1] * (1.0 if VV[0, -1] >= 0 else -1.0)
        gg = float(ee[-1] - ee[-2])
        c1 = float(uc_vec(10, NN) @ p1)
        ck(offA.min() > 0.0, f"A off-diagonal strictly positive at wbar={wbar}")
        ck(gg > 0.0, f"lam_1(A) simple (gap > 0) at wbar={wbar}")
        ck(p1.min() > 0.0, f"psi_1 > 0 entrywise at wbar={wbar}")
        ck(c1 >= p1[0] - 1e-14, f"c_1 >= psi_1[0] > 0 at wbar={wbar}")
        print(f"    {wbar:>7.2f} {offA.min():>14.4e} {gg:>16.6e} {p1.min():>11.3e} "
              f"{p1[0]:>9.4f} {c1:>13.6f}")
    print("    (all three inputs hold at every wbar tested; the closed form of")
    print("     A_{nm} makes the first one analytic, not numerical.)")



# =============================================================== STAGE S1 (cone)
def run_S1():
    head("S1.  The cone route: generalized Lemma B with one negative generator "
         "-- and why it is quantitatively dead")

    NN = 200
    L = np.arange(1, NN)[:, None]
    Lp = np.arange(1, NN + 1)[None, :]
    offd = L != Lp

    print("\n(a) GENERALIZED LEMMA B (exact).  O[F] - mutilde K_c is Metzler in the")
    print("    step basis  <=>  mutilde <= mu_Metz := min_{L!=L', Wc>0} W^F/Wc.")
    print("    Wc >= 0 is supported on the N_c rows L <= N_c (S0d), so the min runs")
    print("    over that block only and the eps of the 'one bad generator' lemma is")
    print("    EXPLICIT and closed-form.\n")
    print(f"    {'wbar':>6} {'N_c':>5} {'mu_Metz':>12} {'argmin (L,L*)':>15} "
          f"{'mu_Metz*N_c^3/wbar^2':>21}")
    tab = {}
    for wbar in (0.5, 1.0, 2.0, 5.0):
        for Nc in (2, 3, 5, 10, 20, 40, 80):
            WF = kamke(H_mat(NN, wbar, 0))
            Wc = kamke(K_c(Nc, NN))
            msk = offd & (Wc > 1e-14)
            r = np.where(msk, WF / np.where(Wc > 1e-14, Wc, 1.0), np.inf)
            i = np.unravel_index(np.argmin(r), r.shape)
            mm = float(r[i])
            tab[(wbar, Nc)] = mm
            print(f"    {wbar:>6.2f} {Nc:>5d} {mm:>12.5e} "
                  f"{str((int(i[0])+1, int(i[1])+1)):>15} {mm*Nc**3/wbar**2:>21.4f}")
        print()
    for wbar in (0.5, 1.0, 2.0, 5.0):
        ck(tab[(wbar, 40)] < tab[(wbar, 10)],
           f"mu_Metz decreasing in N_c at wbar={wbar}")
        ratio = tab[(wbar, 10)] / tab[(wbar, 40)]
        ck(20.0 < ratio < 200.0, f"mu_Metz ~ N_c^-3 scaling at wbar={wbar}")
        print(f"    wbar={wbar:>5}:  mu_Metz(10)/mu_Metz(40) = {ratio:7.2f} "
              f"(N_c^3 would give 64)")
    ck(tab[(1.0, 10)] < 0.01, "mu_Metz(wbar=1,N_c=10) far below physical mutilde")
    print("\n    VERDICT: mu_Metz ~ wbar^2 / N_c^3, e.g. 1.5e-3 at (1,10) and 2.7e-5")
    print("    at (1,40), against physical mutilde in [0.04, 0.33].  The generalized")
    print("    Lemma B is TRUE with an explicit eps but USELESS.  ==> abandon the cone.")

    print("\n(b) how badly does the top eigenvector leave C?  (cone violations)")
    print(f"    {'wbar':>6} {'N_c':>5} {'mutilde':>8} {'#(v<0)':>7} {'#(increases)':>13} "
          f"{'#sign changes':>14} {'v_{N_c}/max|v|':>15}")
    for (wbar, Nc) in ((1.0, 10), (2.0, 20), (0.74, 4), (5.0, 40)):
        A = H_mat(NN, wbar, 0)
        for mut in (0.0, 0.065, 0.13, 0.26, 0.33):
            _, _, v = top_pencil(A, uc_vec(Nc, NN), mut)
            a = np.abs(v).max()
            neg = int((v < -1e-11 * a).sum())
            inc = int(((v[:-1] - v[1:]) < -1e-11 * a).sum())
            nz = v[np.abs(v) > 1e-12 * a]
            sc = int((np.sign(nz[:-1]) != np.sign(nz[1:])).sum())
            print(f"    {wbar:>6.2f} {Nc:>5d} {mut:>8.3f} {neg:>7d} {inc:>13d} "
                  f"{sc:>14d} {v[min(Nc,NN-1)]/a:>+15.4e}")
            if mut == 0.0:
                ck(neg == 0 and inc == 0, f"mutilde=0: v in C at ({wbar},{Nc})")
            if mut >= 0.13:
                ck(neg > 0, f"mutilde={mut}: v leaves C at ({wbar},{Nc})")
        print()

    print("(c) BACKUP ROUTE REFUTED: the number of sign changes is 2, not <= 1,")
    print("    and it is N-independent.")
    for NN2 in (200, 400, 800):
        A = H_mat(NN2, 1.0, 0)
        _, _, v = top_pencil(A, uc_vec(10, NN2), 0.26)
        a = np.abs(v).max()
        nz = v[np.abs(v) > 1e-12 * a]
        sc = int((np.sign(nz[:-1]) != np.sign(nz[1:])).sum())
        ck(sc == 2, f"exactly 2 sign changes at N={NN2}")
        print(f"    N = {NN2:>4d}:  sign changes of the top eigenvector = {sc}")
    print("    ==> a Gantmacher-Krein 'at most one sign change' cone cannot work.")


# ============================================== STAGE S2 (the rank-one theorem)
def run_S2():
    head("S2.  The rank-one resolvent route: Lemmas 1-3 and the threshold "
         "mutilde_*(wbar,N_c)")

    print("\n(a) LEMMA 1 (secular reduction).  HYPOTHESIS: c_1 != 0 AND")
    print("    kappa := 2 mutilde c_1^2 / gam < 1.  Then lam := lam_1(M) is the unique")
    print("    root in (mu_2, mu_1) of  phi(lam) := 2 mutilde sum_j c_j^2/(mu_j-lam)")
    print("    = 1; it is SIMPLE; the top eigenvector is EXACTLY (A-lam)^{-1} u_c;")
    print("    and mu_1 - lam <= 2 mutilde c_1^2.")
    print("    WHY kappa < 1 IS NEEDED (reviewer CHALLENGED 2): the old proof said")
    print("    'phi runs from -infinity to +infinity on (mu_2,mu_1)', which needs")
    print("    c_2 != 0 -- never argued, and min|c_2| = 1.9e-3 is attained at exactly")
    print("    the counterexample cell (0.74,4).  With kappa < 1 no c_2 is needed:")
    print("      lam_0 := mu_1 - 2 mutilde c_1^2 > mu_2   (this IS kappa < 1), and")
    print("      phi(lam_0) = 1 - 2 mutilde sum_{j>=2} c_j^2/(lam_0-mu_j) <= 1,")
    print("      phi(mu_1^-) = +infinity,  phi increasing  ==>  the root lies in")
    print("      [lam_0, mu_1) (subset of (mu_2,mu_1)), which also re-proves (iii).")
    print("      Every 'inherited' eigenvalue mu_j (those with c_j = 0) is <= mu_2 <")
    print("      lam, so the root IS lam_1(M) and it is simple.")
    print("    The theorem's hypothesis mutilde < mutilde_* IMPLIES kappa < 1 (see")
    print("    (e) / run_verdict), so the patch costs nothing -- but Lemma 1 is a")
    print("    stand-alone statement and must carry the hypothesis itself.")
    print(f"    {'wbar':>6} {'N_c':>5} {'mutilde':>8} {'kappa':>7} "
          f"{'|lam_sec-lam_eigh|':>19} {'|v_sec-v_eigh|':>15} {'mu_1-lam':>10} "
          f"{'2 mut c_1^2':>12} {'phi(lam_0)':>11} {'lam_1-lam_2(M)':>15}")
    for (wbar, Nc) in ((1.0, 10), (0.74, 4), (2.0, 20), (5.0, 40)):
        C = Cell(wbar, Nc)
        A = H_mat(N, wbar, 0)
        for mut in (0.065, 0.13, 0.26, 0.33):
            lam_e, _, v_e = top_pencil(A, uc_vec(Nc), mut)
            y, lam_s = C.vhat_coeffs(mut)
            e, V = np.linalg.eigh(A)
            v_s = V @ y
            if float(v_s @ v_e) < 0:
                v_s = -v_s
            dl, dv = abs(lam_s - lam_e), float(np.linalg.norm(v_s - v_e))
            kap = C.kappa(mut)
            lam0 = C.mu1 - 2 * mut * C.c1 ** 2
            phi0 = 2 * mut * float(np.sum(C.c ** 2 / (C.e - lam0)))
            eM = np.linalg.eigvalsh(A - 2.0 * mut * np.outer(uc_vec(Nc),
                                                             uc_vec(Nc)))
            gapM = float(eM[-1] - eM[-2])
            tag = "" if kap < 1.0 else "   (kappa >= 1: Lemma 1 makes no claim)"
            # (i)(ii)(iii) are verified numerically everywhere -- they happen to
            # survive kappa >= 1 as long as the root stays in (mu_2,mu_1) -- but
            # the LEMMA only claims them under kappa < 1.
            ck(dl < 1e-11, f"Lemma 1 secular root ({wbar},{Nc},{mut})")
            ck(dv < 1e-7, f"Lemma 1 eigenvector formula ({wbar},{Nc},{mut})")
            ck(C.mu1 - lam_s <= 2 * mut * C.c1 ** 2 + 1e-12,
               f"Lemma 1 bound mu_1-lam <= 2 mut c_1^2 ({wbar},{Nc},{mut})")
            if kap < 1.0:
                ck(lam0 > C.mu2, f"Lemma 1: kappa<1 gives lam_0 > mu_2 "
                                 f"({wbar},{Nc},{mut})")
                ck(phi0 <= 1.0 + 1e-12,
                   f"Lemma 1: phi(lam_0) <= 1 ({wbar},{Nc},{mut})")
                ck(lam0 - 1e-12 <= lam_s < C.mu1,
                   f"Lemma 1: root in [lam_0,mu_1) ({wbar},{Nc},{mut})")
                ck(gapM > 0.0, f"Lemma 1: lam_1(M) simple ({wbar},{Nc},{mut})")
            print(f"    {wbar:>6.2f} {Nc:>5d} {mut:>8.3f} {kap:>7.4f} {dl:>19.3e} "
                  f"{dv:>15.3e} {C.mu1-lam_s:>10.3e} {2*mut*C.c1**2:>12.3e} "
                  f"{phi0:>11.6f} {gapM:>15.3e}{tag}")


    print("\n    how small can c_2 get?  (the quantity the OLD proof silently needed)")
    worst_c2, arg_c2 = np.inf, None
    for wbar in (0.5, 0.74, 1.0, 2.0, 5.0, 10.0):
        for Nc in (2, 4, 10, 20, 40, 80):
            if Nc > N // 2:
                continue
            C = Cell(wbar, Nc)
            a = abs(float(C.c[-2]))
            if a < worst_c2:
                worst_c2, arg_c2 = a, (wbar, Nc)
    print(f"    min |c_2| over the (wbar,N_c) grid = {worst_c2:.4e} at {arg_c2}")
    ck(worst_c2 < 1e-2, "min |c_2| is tiny somewhere on the grid")
    ck(arg_c2 == (0.74, 4), "the min |c_2| cell is the counterexample cell (0.74,4)")
    print("    ==> c_2 is NOT safely bounded away from 0, and its worst cell is")
    print("        exactly the counterexample cell.  The kappa<1 route avoids c_2.")

    print("\n(b) LEMMA 2 (angle bound).  If kappa := 2 mutilde c_1^2/gam < 1 then")

    print("        tan theta <= 2 mutilde c_1 |u_perp| / (gam - 2 mutilde c_1^2).")
    print(f"    {'wbar':>6} {'N_c':>5} {'mutilde':>8} {'kappa':>8} {'tan theta':>11} "
          f"{'bound':>11} {'ok':>5}")
    for (wbar, Nc) in ((1.0, 10), (0.74, 4), (2.0, 20), (5.0, 40)):
        C = Cell(wbar, Nc)
        for mut in (0.065, 0.13, 0.26):
            kap = C.kappa(mut)
            t = C.tan_theta(mut)
            if kap >= 1.0:
                print(f"    {wbar:>6.2f} {Nc:>5d} {mut:>8.3f} {kap:>8.4f} {t:>11.4f} "
                      f"{'--':>11} {'n/a':>5}   (kappa >= 1: lemma inapplicable)")
                continue
            bd = 2 * mut * C.c1 * C.up / (C.gam - 2 * mut * C.c1 ** 2)
            ck(t <= bd * (1 + 1e-9) + 1e-12,
               f"Lemma 2 angle bound ({wbar},{Nc},{mut})")
            print(f"    {wbar:>6.2f} {Nc:>5d} {mut:>8.3f} {kap:>8.4f} {t:>11.4f} "
                  f"{bd:>11.4f} {'OK':>5}")

    print("\n(c) LEMMA 3 (robust Cor A').  For vhat = cos th psi_1 + sin th phi,")
    print("        <vhat,O[f]vhat> >= cos^2 th [ q_1 - 2 b tan th - m tan^2 th ] > 0")
    print("    whenever tan th < t_max = (-b + sqrt(b^2 + m q_1))/m.")
    print(f"    {'wbar':>6} {'N_c':>5} {'mutilde':>8} {'q(vhat)':>12} {'lower bd':>12} "
          f"{'tan th':>9} {'t_max':>9} {'ok':>5}")
    for (wbar, Nc) in ((1.0, 10), (0.74, 4), (2.0, 20), (5.0, 40)):
        C = Cell(wbar, Nc)
        for mut in (0.065, 0.13, 0.26):
            t = C.tan_theta(mut)
            cos2 = 1.0 / (1.0 + t * t)
            lb = cos2 * (C.q1 - 2 * C.b * t - C.m * t * t)
            q = C.q_of_mut(mut)
            ck(q >= lb - 1e-12, f"Lemma 3 lower bound ({wbar},{Nc},{mut})")
            print(f"    {wbar:>6.2f} {Nc:>5d} {mut:>8.3f} {q:>12.4e} {lb:>12.4e} "
                  f"{t:>9.4f} {C.tmax:>9.4f} {'OK':>5}")

    print("\n(d) LEMMA 4 (Hellmann-Feynman).  T d lam_1(M)/dT = -<vhat,O[f]vhat>")
    h = 1e-6
    for (wbar, mut, Nc) in ((1.0, 0.13, 10), (2.0, 0.26, 20), (0.74, 0.13, 4)):
        uc = uc_vec(Nc)
        l0 = top_pencil(H_mat(N, wbar - h, 0), uc, mut)[0]
        l1 = top_pencil(H_mat(N, wbar + h, 0), uc, mut)[0]
        fd = (l1 - l0) / (2 * h) * wbar          # wbar d/dwbar = -T d/dT
        C = Cell(wbar, Nc)
        q = C.q_of_mut(mut)
        ck(abs(fd - q) <= 1e-6 * max(1.0, abs(q)), f"HF at ({wbar},{mut},{Nc})")
        print(f"    wbar={wbar:>5} mut={mut:<5} N_c={Nc:<3d}:  FD = {fd:.10f}   "
              f"q(vhat) = {q:.10f}   rel {abs(fd-q)/abs(q):.2e}")

    print("\n(e) THEOREM: the threshold mutilde_*  vs the sharp numerical mutilde_crit")
    print(f"    {'wbar':>6} {'N_c':>5} {'mut_*(PROVED)':>14} {'mut_crit(num)':>14} "
          f"{'gam':>9} {'c_1':>8} {'|u_perp|':>9} {'t_max':>8} {'kappa@crit':>11}")
    worst = np.inf
    for wbar in (0.5, 0.74, 1.0, 2.0, 5.0, 10.0):
        for Nc in (2, 4, 10, 20, 40, 80):
            if Nc > N // 2:
                continue
            C = Cell(wbar, Nc)
            mc = C.mut_crit()
            ck(C.thr < mc * (1 + 1e-9),
               f"mut_* <= mut_crit at ({wbar},{Nc})")   # the theorem is not vacuous
            kap = C.kappa(mc) if np.isfinite(mc) else np.inf
            worst = min(worst, C.thr)
            print(f"    {wbar:>6.2f} {Nc:>5d} {C.thr:>14.5f} {mc:>14.5f} "
                  f"{C.gam:>9.5f} {C.c1:>8.4f} {C.up:>9.4f} {C.tmax:>8.4f} "
                  f"{kap:>11.4f}")
        print()
    print(f"    worst mutilde_* over this grid = {worst:.5f}")
    print("    NOTE kappa@mut_crit is O(1) throughout: kappa = 2 mutilde c_1^2/gam is")
    print("    the true control parameter, and the obstruction is a NEAR-CROSSING of")
    print("    lam_1(M) with mu_2(A), not a gradual cone exit.")
    print("    NOTE mutilde_* decays only LOGARITHMICALLY in N_c (c_1^2 ~ ln N_c),")
    print("    against N_c^{-3} for the cone route -- that is the whole gain.")


# ================================================ STAGE S3 (physical anchors)
ANCHORS = [(3.0, 0.13, 0.60, 10), (2.0, 0.13, 0.98, 10), (2.0, 0.10, 0.98, 10),
           (1.0, 0.13, 2.20, 10), (1.0, 0.10, 2.00, 10), (0.7, 0.13, 4.56, 10),
           (0.5, 0.13, 13.5, 10), (0.5, 0.13, 15.6, 3), (0.4, 0.13, 51.2, 3),
           (2.0, 0.13, 3.00, 3), (1.0, 0.13, 2.50, 3)]

WIN_PTS = 9        # sample points per window, left endpoint included


def window_inf(Nc, r, npts=WIN_PTS, NN=None):
    """min of mutilde_*(., N_c) over an `npts`-sample of [N_c/r, (N_c+1)/r).

    This is an UPPER BOUND for the true inf (a sample, not an enclosure) -- see
    (H1)/Sec.8 item 2 of the note.  The sample is nondecreasing in wbar on every
    window tested, so the argmin is the LEFT endpoint and the bound is tight.
    """
    lo, hi = Nc / r, (Nc + 1) / r
    ws = np.linspace(lo, hi, npts + 1)[:-1]
    vals = np.array([Cell(w, Nc, NN).thr for w in ws])
    i = int(np.argmin(vals))
    mono = bool(np.all(np.diff(vals) > -1e-12))
    return float(vals[i]), float(ws[i]), mono, (lo, hi)


def run_S3():
    head("S3.  Physical coverage: every anchor of `P4D_coulomb_operator.md`, "
         "with the WINDOW INF")

    print("\n    anchors (lambda, mu*, wbar_c, r) from the note's physics table;")
    print("    N_c = floor(r * wbar_c),  mutilde = mu*/lambda.")
    print("    The theorem's hypothesis is  mutilde < inf over the WINDOW")
    print("    W(N_c) = [N_c/r, (N_c+1)/r)  of mutilde_*(., N_c) -- NOT the value at")
    print(f"    the nodal point wbar_c.  Sampled at {WIN_PTS} points of each window")
    print("    (left endpoint included, right endpoint excluded).  mutilde_* is")
    print("    nondecreasing in wbar on every window here, so the argmin is the left")
    print("    endpoint; the sample is still a SAMPLE, not an enclosure (H1).\n")
    print(f"    {'lambda':>7} {'mu*':>6} {'r':>4} {'wbar_c':>8} {'N_c':>5} "
          f"{'mutilde':>8} {'mut_*(wbar_c)':>14} {'inf_W mut_*':>12} {'argmin':>8} "
          f"{'drop%':>7} {'mut_crit':>10} {'kappa':>8} {'verdict':>9}")
    nproved = 0
    worst_margin, worst_drop = np.inf, 0.0
    for (lam, mus, wc, r) in ANCHORS:
        mut = mus / lam
        Nc = max(1, int(np.floor(r * wc)))
        if Nc > N // 2:
            print(f"    {lam:>7.2f} {mus:>6.2f} {r:>4d} {wc:>8.2f} {Nc:>5d} "
                  f"{mut:>8.4f}   skipped (N_c > N/2; raise PN)")
            continue
        C = Cell(wc, Nc)
        mc = C.mut_crit()
        inf_w, arg_w, mono, (wlo, whi) = window_inf(Nc, r)
        ok = "PROVED" if mut < inf_w else ("nodal" if mut < C.thr
                                          else ("numeric" if mut < mc else "FAILS"))
        if ok == "PROVED":
            nproved += 1
            ck(mut < inf_w,
               f"anchor lambda={lam} mu*={mus} r={r}: mutilde < inf_W mut_*")
            ck(mut < C.thr, f"anchor lambda={lam} mu*={mus} r={r} proved at wbar_c")
        ck(inf_w <= C.thr * (1 + 1e-9),
           f"window inf <= nodal value at lambda={lam} r={r}")
        ck(mono, f"mut_*(.,N_c) nondecreasing in wbar on W({Nc}) (r={r})")
        worst_margin = min(worst_margin, inf_w / mut)
        worst_drop = min(worst_drop, inf_w / C.thr - 1.0)
        print(f"    {lam:>7.2f} {mus:>6.2f} {r:>4d} {wc:>8.2f} {Nc:>5d} {mut:>8.4f} "
              f"{C.thr:>14.5f} {inf_w:>12.5f} {arg_w:>8.3f} "
              f"{100*(inf_w/C.thr-1):>+7.2f} {mc:>10.5f} {C.kappa(mut):>8.4f} "
              f"{ok:>9}")
    print(f"\n    {nproved} anchors PROVED by the WINDOW-INF hypothesis")
    print(f"    (mutilde < inf_W mutilde_*).  Worst margin inf_W mut_*/mutilde = "
          f"{worst_margin:.3f};")
    print(f"    largest within-window drop of mut_* below its nodal value = "
          f"{100*worst_drop:+.2f}%.")
    print("    All have kappa < 0.5, i.e. they sit at half the Lemma-2 radius or less.")
    print("    Physical ranges: lambda >= 2 (hydrides) => mutilde <= 0.065; lambda = 1")
    print("    => 0.13; the lambda < 1 corner reaches mutilde = 0.325 and is still")
    print("    covered.")
    print("    CAVEAT (H1): these mut_* are finite-N OVERESTIMATES (S6); the margins")
    print("    above absorb the drift, the tightest case losing 0.8% more.")

    print(f"\n    REFINE: the {WIN_PTS}-point sample could miss an interior dip, so")
    print("    the three windows that matter most (the largest within-window drop,")
    print("    and the two tightest margins) are re-sampled at 33 points:")
    print(f"    {'r':>4} {'N_c':>5} {'window':>20} {'inf_W mut_* (33pt)':>19} "
          f"{'argmin':>9} {'= left end':>11} {'nondecr':>9}")
    for (r, Nc) in ((10, 9), (10, 135), (3, 153)):
        if Nc > N // 2:
            print(f"    {r:>4d} {Nc:>5d}   skipped (N_c > N/2; raise PN)")
            continue
        inf_w, arg_w, mono, (wlo, whi) = window_inf(Nc, r, npts=33)
        at_left = abs(arg_w - wlo) < 1e-12
        ck(mono, f"33-point refine: mut_* nondecreasing on W({Nc}) (r={r})")
        ck(at_left, f"33-point refine: argmin of mut_* is the left end of W({Nc})")
        print(f"    {r:>4d} {Nc:>5d} {'[%.4f,%.4f)' % (wlo, whi):>20} "
              f"{inf_w:>19.6f} {arg_w:>9.4f} {str(at_left):>11} {str(mono):>9}")
    print("    This still does not make the sample an enclosure; it only rules out")
    print("    dips on a 3x finer grid (H1).")



# ================================================ STAGE S4 (N_c step jumps)
def run_S4():
    head("S4.  Across N_c windows: rigorous jump sandwich, and where global "
         "uniqueness genuinely FAILS")

    print("\n(a) TWO-SIDED RAYLEIGH SANDWICH for the jump lam_1(N_c)-lam_1(N_c+1).")
    print("    With w = top vector at N_c, v = top vector at N_c+1, and")
    print("    u := (2N_c+1)^{-1/2}:")
    print("        2 mut u v_{N_c} ( 2<u_c^{(N_c)},v> + u v_{N_c} )")
    print("            <= lam_1(N_c) - lam_1(N_c+1) <=")
    print("        2 mut u w_{N_c} ( 2<u_c^{(N_c)},w> + u w_{N_c} )")
    print("    (both sides are Rayleigh quotients of M(N_c+1) resp. M(N_c)).")
    print(f"\n    {'r':>3} {'mut':>6} {'N_c':>5} {'wbar':>8} {'jump(+=up)':>12} "
          f"{'lower bd':>11} {'upper bd':>11} {'v_{N_c}':>11}")
    nup = ndn = 0
    for r in (10, 3):
        for mut in (0.13,):
            for Nc in range(1, 16):
                wb = (Nc + 0.5) / r
                A = H_mat(N, wb, 0)
                u0, u1 = uc_vec(Nc), uc_vec(Nc + 1)
                l0, _, w = top_pencil(A, u0, mut)
                l1, _, vs = top_pencil(A, u1, mut)
                uN = 1.0 / np.sqrt(2.0 * Nc + 1.0)
                lo = 2 * mut * uN * vs[Nc] * (2 * float(u0 @ vs) + uN * vs[Nc])
                hi = 2 * mut * uN * w[Nc] * (2 * float(u0 @ w) + uN * w[Nc])
                ck(lo - 1e-10 <= l0 - l1 <= hi + 1e-10,
                   f"jump sandwich r={r} mut={mut} N_c={Nc}")
                if l1 > l0:
                    nup += 1
                else:
                    ndn += 1
                if Nc <= 8:
                    print(f"    {r:>3d} {mut:>6.3f} {Nc:>3d}->{Nc+1:<2d} {wb:>8.3f} "
                          f"{l1-l0:>+12.3e} {lo:>11.3e} {hi:>11.3e} {vs[Nc]:>+11.3e}")
            print()
    print(f"    sandwich verified on all {nup+ndn} sampled steps "
          f"({nup} upward, {ndn} downward).")

    print("(b) explicit NON-UNIQUENESS windows on the true stepped path "
          "N_c(wbar) = floor(r wbar),  wbar in [0.34, 6), step 0.004")
    for r in (10, 3):
        for mut in (0.065, 0.10, 0.13, 0.26):
            ws = np.arange(0.34, 6.0, 0.004)
            vals = np.empty(len(ws))
            for i, wb in enumerate(ws):
                Nc = max(0, int(np.floor(r * wb)))
                vals[i] = top_pencil(H_mat(N, wb, 0), uc_vec(Nc), mut)[0]
            run = np.maximum.accumulate(vals)
            bad = np.where(vals < run - 1e-12)[0]
            if len(bad) == 0:
                print(f"    r={r:>2d} mut={mut:<5}: lam_1 nondecreasing on the WHOLE "
                      f"scanned range [0.34,6)")
                print(f"    {'':>18}  ==> unique Tc for every lambda WHOSE ROOT "
                      f"lam_1 = 1/lambda")
                print(f"    {'':>18}      IS ATTAINED IN [0.34,6), i.e. lambda in "
                      f"[{1/vals.max():.4f},{1/vals.min():.4f}].")
                print(f"    {'':>18}      NOT 'for every lambda': see (c) for the "
                      f"real dips at wbar < 0.34.")
                continue

            segs, i = [], 0
            while i < len(bad):
                j = i
                while j + 1 < len(bad) and bad[j + 1] == bad[j] + 1:
                    j += 1
                s = bad[i:j + 1]
                segs.append((ws[s[0]], ws[s[-1]], run[s[0]], vals[s].min()))
                i = j + 1
            print(f"    r={r:>2d} mut={mut:<5}: {len(segs)} downward stretches; "
                  f"first 4 bad lambda-intervals:")
            for (a, b_, hi, lo) in segs[:4]:
                print(f"        wbar in [{a:.3f},{b_:.3f}]: lam_1 dips "
                      f"{hi:.6f} -> {lo:.6f}  ==> NON-UNIQUE for lambda in "
                      f"({1/hi:.5f}, {1/lo:.5f})")
            hit = [s for s in segs if 1.0 / s[2] <= 1.0 <= 1.0 / s[3]]
            if hit:
                a, b_, hi, lo = hit[0]
                print(f"        *** lambda = 1 IS non-unique here: wbar in "
                      f"[{a:.3f},{b_:.3f}], lam_1 {hi:.6f} -> {lo:.6f}")
    print("\n(c) THE LOW-wbar CAVEAT: r = 10 is NOT monotone below wbar = 0.34.")
    print("    True stepped path on [0.05, 0.40], step 0.005.  The N_c = 0 -> 1 step")
    print("    is a genuine, large DOWNWARD jump (at N_c = 0 there is no Coulomb term")
    print("    at all, so lam_1 = mu_1(A); switching on K_c crushes it).")
    print(f"    {'mut':>6} {'#down steps':>12} {'wbar of down steps':>24} "
          f"{'first dip lam_1':>22} {'bad lambda interval':>24}")
    for mut in (0.065, 0.10, 0.13, 0.26):
        ws = np.arange(0.05, 0.4005, 0.005)
        vals = np.array([top_pencil(H_mat(N, wb, 0),
                                    uc_vec(max(0, int(np.floor(10 * wb)))), mut)[0]
                         for wb in ws])
        d = np.diff(vals)
        dn = np.where(d < -1e-12)[0]
        ck(len(dn) > 0, f"r=10 has downward steps below wbar=0.34 (mut={mut})")
        ck(bool(np.all(ws[dn + 1] < 0.34)),
           f"all r=10 downward steps are below wbar=0.34 (mut={mut})")
        i = int(dn[0])
        print(f"    {mut:>6.3f} {len(dn):>12d} "
              f"{str([round(float(x),3) for x in ws[dn+1]]):>24} "
              f"{vals[i]:>10.6f} -> {vals[i+1]:>9.3e} "
              f"({1/vals[i]:>7.1f}, {1/vals[i+1]:>8.2e})")
    print("    ==> 'unique Tc for EVERY lambda' is FALSE as stated; the true claim is")
    print("    'for every lambda whose root lies in the scanned wbar-range'.  The bad")
    print("    lambda are ~77 ... 1e4, i.e. unphysical, but the quantifier must be")
    print("    written (reviewer DOWNGRADE 7).")

    print("\n(d) the docstring (H2) claim reaches wbar = 8: verify (6.0, 8.2] too")
    print("    (step 0.02; the N_c steps of r = 10 are 0.1 apart, so every step is")
    print("    straddled).")
    for mut in (0.065, 0.13, 0.26):
        ws = np.arange(6.0, 8.205, 0.02)
        vals = np.array([top_pencil(H_mat(N, wb, 0),
                                    uc_vec(int(np.floor(10 * wb))), mut)[0]
                         for wb in ws])
        d = np.diff(vals)
        dn = int((d < -1e-12).sum())
        ck(dn == 0, f"r=10 monotone on (6.0,8.2] at mut={mut}")
        print(f"    r=10 mut={mut:<6}: {dn} downward steps on (6.0,8.2]; lam_1 in "
              f"[{vals.min():.6f}, {vals.max():.6f}]")

    print("\n    ==> For r = w_c/<w> >= 10 the stepped path is monotone on")
    print("    [0.34, 8.2] and Tc is unique for every lambda whose root lies there.")
    print("    Below wbar = 0.34 it is not (part c).  For r = 3 it is not, at physical")
    print("    mutilde -- the failure is the DISCONTINUITY of N_c(T), not within-window")
    print("    behaviour, so Tc must be defined as a sup (cf.")
    print("    `P4D_coulomb_operator.md` Layer 1).")



# =============================================== STAGE S5 (the counterexample)
def run_S5():
    head("S5.  HARD TEST: the known within-window counterexample "
         "(N_c = 4, wbar ~ 0.74, mutilde in [0.322,0.334])")

    print("\n    Independently confirmed at N = 600 and N = 1200 by the sibling probe")
    print("    `p4_uniqueness_coulomb_probe.py`: <v,H'v> = -3.24e-5 there, and the")
    print("    non-monotone mutilde window is EXACTLY [~0.322, ~0.334], for N_c = 4")
    print("    ONLY (N_c = 1..3 and 5..40 are monotone at mutilde = 0.33).")
    print("\n(a) our machinery must REPRODUCE the failure ...")
    print(f"    {'wbar':>7} {'mutilde':>9} {'q(vhat)':>13} {'monotone?':>10}")
    Cbad = Cell(0.74, 4)
    for mut in (0.300, 0.320, 0.325, 0.330, 0.335, 0.340):
        q = Cbad.q_of_mut(mut)
        print(f"    {0.74:>7.2f} {mut:>9.3f} {q:>13.4e} "
              f"{'yes' if q > 0 else 'NO':>10}")
    q325, q330 = Cbad.q_of_mut(0.325), Cbad.q_of_mut(0.330)
    q320, q340 = Cbad.q_of_mut(0.320), Cbad.q_of_mut(0.340)
    ck(q320 > 0 and q340 > 0, "monotone outside the bad window at wbar=0.74,N_c=4")
    ck(q325 <= 0 or q330 <= 0, "non-monotone inside [0.322,0.334] reproduced")
    mc = Cbad.mut_crit()
    print(f"\n    mutilde_crit(0.74, 4) = {mc:.5f}   (sibling: window starts ~0.322)")
    ck(0.31 < mc < 0.335, "mut_crit reproduces the sibling window start")

    print("\n(b) ... and our THEOREM must NOT claim it.")
    print(f"    mutilde_*(0.74, 4) = {Cbad.thr:.5f}  <  0.322   ==> the bad point is")
    print("    outside the theorem's hypothesis.  EXCLUDED CORRECTLY.")
    ck(Cbad.thr < 0.322, "the threshold excludes the counterexample")

    print("\n(c) the MECHANISM: kappa = 2 mutilde c_1^2 / gam at the bad point")
    for mut in (0.13, 0.26, 0.33):
        print(f"    mutilde = {mut:>5.3f}:  kappa = {Cbad.kappa(mut):.4f}"
              + ("   <-- Lemma 2 hypothesis kappa<1 VIOLATED" if Cbad.kappa(mut) >= 1
                 else ""))
    lam = Cbad.lam_of_mut(0.33)
    print(f"    lam_1(M) = {lam:.6e},  mu_2(A) = {Cbad.mu2:.6e},  "
          f"lam_1(M)-mu_2 = {lam-Cbad.mu2:.3e}")
    ck(Cbad.kappa(0.33) > 1.5, "kappa >> 1 at the counterexample")
    print("    ==> lam_1(M) is driven down to within ~3e-4 of mu_2(A): an isolated")
    print("    NEAR-CROSSING.  tan theta blows up, Lemma 3's quadratic goes negative.")
    print("    The bad-Abel-coefficient picture does NOT single out N_c = 4; the")
    print("    near-crossing picture DOES (it needs 2 mutilde c_1^2 ~ gam, and the")
    print("    ratio c_1^2/gam peaks at small N_c and small wbar).")

    print("\n(d) sweep: is the theorem's hypothesis violated ONLY where it must be?")
    print(f"    {'wbar':>7} {'N_c':>5} {'mut_*':>9} {'mut_crit':>10} {'ratio':>7}")
    for wbar in (0.74,):
        for Nc in range(1, 13):
            C = Cell(wbar, Nc)
            mc = C.mut_crit()
            ck(C.thr <= mc * (1 + 1e-9), f"mut_* <= mut_crit at N_c={Nc}")
            print(f"    {wbar:>7.2f} {Nc:>5d} {C.thr:>9.5f} {mc:>10.5f} "
                  f"{mc/C.thr:>7.2f}")
    print("    ratio mut_crit/mut_* is the unproved margin (H3): factor ~1.4-3")
    print("    (adversarial review: min = 1.4025 over 224 cells, at (20,50)).")


# =============================================== STAGE S6 (N-convergence)
def run_S6():
    head("S6.  N-truncation: stability of the six threshold constants, on the "
         "LOAD-BEARING cells too")

    print("\n    Two SMALL cells (where 6-7 digits really do hold) and four cells that")
    print("    carry physical anchors (where they do NOT).  All numbers float64 eigh.")
    print(f"\n    {'wbar':>6} {'N_c':>5} {'N':>6} {'gam':>11} {'q_1':>11} {'b':>11} "
          f"{'m':>11} {'c_1':>9} {'t_max':>9} {'mut_*':>10}")
    cells = [(1.0, 10, (100, 200, 400, 800)),
             (0.74, 4, (100, 200, 400, 800)),
             (4.56, 45, (200, 400, 800, 1600)),
             (13.5, 135, (400, 800, 1600)),
             (15.6, 46, (200, 400, 800, 1600)),
             (51.2, 153, (400, 800, 1600))]
    summary = []
    for (wbar, Nc, Ns) in cells:
        vals = []
        for NN in Ns:
            C = Cell(wbar, Nc, NN)
            vals.append(C.thr)
            print(f"    {wbar:>6.2f} {Nc:>5d} {NN:>6d} {C.gam:>11.7f} {C.q1:>11.7f} "
                  f"{C.b:>11.7f} {C.m:>11.7f} {C.c1:>9.6f} {C.tmax:>9.6f} "
                  f"{C.thr:>10.6f}")
        # (i) the drift is ONE-SIDED: mut_*(N) decreases, so finite N OVERESTIMATES
        for i in range(len(vals) - 1):
            ck(vals[i + 1] <= vals[i] + 1e-12,
               f"mut_*(N) nonincreasing in N at ({wbar},{Nc}), "
               f"N={Ns[i]}->{Ns[i+1]}")
        rel = abs(vals[-1] - vals[-2]) / max(1e-30, abs(vals[-1]))
        dig = -np.log10(max(rel, 1e-16))
        summary.append((wbar, Nc, Ns, vals, rel, dig))
        print(f"    {'':>12} -> monotone DECREASING; rel(N={Ns[-2]},N={Ns[-1]}) = "
              f"{rel:.1e}  ({dig:.1f} digits)")
        print()

    print("    HONEST SUMMARY (H1), replacing the old '6 significant digits from")
    print("    N = 200 on' claim, which is true only on the two small cells:")
    print(f"    {'wbar':>7} {'N_c':>5} {'mut_* at largest N':>19} "
          f"{'rel(last two N)':>16} {'digits':>7} {'verdict':>18}")
    for (wbar, Nc, Ns, vals, rel, dig) in summary:
        vd = ("6+ digits" if rel < 1e-6 else
              ("4-5 digits" if rel < 1e-4 else "2-3 digits ONLY"))
        print(f"    {wbar:>7.2f} {Nc:>5d} {vals[-1]:>19.6f} {rel:>16.1e} "
              f"{dig:>7.1f} {vd:>18}")
    small = [s for s in summary if s[0] < 2.0]
    heavy = [s for s in summary if s[0] >= 10.0]
    ck(all(s[4] < 1e-6 for s in small), "6+ digits on the small cells")
    ck(any(s[4] > 1e-4 for s in heavy),
       "the load-bearing cells do NOT reach 6 digits (the old claim was false)")
    # the drift channel on the worst cell is gam
    print("\n    the drift channel on the worst cell (51.2,153) is gam:", end=" ")
    for NN in (400, 800, 1600):
        print(f"N={NN}: {Cell(51.2, 153, NN).gam:.7f}", end="  ")
    print()
    ck(Cell(51.2, 153, 1600).gam < Cell(51.2, 153, 400).gam,
       "gam(N) decreasing on (51.2,153)")
    print("    CONCLUSION: the reported mutilde_* are OVERESTIMATES (mut_*(N) falls")
    print("    monotonically in N on all six cells).  6-7 digits only on small cells;")
    print("    on the load-bearing anchors it is 2-3 digits.  The anchors' margins")
    print("    absorb it: the tightest anchor (lambda=0.5, r=10) has 7.3% of room")
    print("    (S8) against a 0.03% N-drift; the lambda=0.4 corner has 19.4% against")
    print("    0.8%.  This is EVIDENCE of truncation-independence, not a proof (H1).")



# =============================================== STAGE S7 (high precision)
def run_S7():
    head("S7.  mpmath spot-check of the secular reduction (dps = 40)")
    try:
        import mpmath as mp
    except ImportError:
        print("    mpmath unavailable -- skipped.")
        return
    mp.mp.dps = 40
    NN = 60
    for (wbar, Nc, mut) in ((1.0, 10, 0.13), (0.74, 4, 0.33)):
        A = H_mat(NN, wbar, 0)
        uc = uc_vec(Nc, NN)
        Am = mp.matrix(NN, NN)
        for i in range(NN):
            for j in range(NN):
                Am[i, j] = mp.mpf(float(A[i, j])) - 2 * mp.mpf(mut) * \
                    mp.mpf(float(uc[i])) * mp.mpf(float(uc[j]))
        ev = mp.eigsy(Am, eigvals_only=True)
        lam_mp = max(ev)
        lam_f = top_pencil(A, uc, mut)[0]
        rel = abs(float(lam_mp) - lam_f) / max(1e-300, abs(lam_f))
        ck(rel < 1e-10, f"mpmath lam_1 spot-check ({wbar},{Nc},{mut})")
        print(f"    wbar={wbar:<5} N_c={Nc:<3d} mut={mut:<5} N={NN}:  "
              f"lam_1 = {mp.nstr(lam_mp, 20)}   float64 rel.err {rel:.2e}")


# ========================================== STAGE S8 (intervalisation budget)
CONST_NAMES = ("gam", "c_1", "|u_perp|", "q_1", "b", "m")
# direction in which an ERROR in that constant LOWERS mutilde_*
CONST_SIGN = (-1, +1, +1, -1, +1, +1)


def thr_from(gam, c1, up, q1, b, m):
    """mutilde_* as an explicit function of the six constants (no linear algebra)."""
    if m > 0:
        t = (-b + np.sqrt(b * b + m * q1)) / m
    else:
        return gam / (2 * c1 * up) if c1 * up > 0 else np.inf
    return t * gam / (2 * c1 * up + 2 * t * c1 ** 2)


def _consts(C):
    return [C.gam, C.c1, C.up, C.q1, C.b, C.m]


def elasticity(C, i, h=1e-7):
    """d ln mutilde_* / d ln X_i, by central difference on the closed form."""
    v = _consts(C)
    base = thr_from(*v)
    vp, vm = list(v), list(v)
    vp[i] = v[i] * (1 + h)
    vm[i] = v[i] * (1 - h)
    return (thr_from(*vp) - thr_from(*vm)) / (2 * h * base)


def flip_error(C, i, mut):
    """signed relative error in X_i that drags mutilde_* down to mut (bisection)."""
    v = _consts(C)
    s = CONST_SIGN[i]
    lo, hi = 0.0, (1.0 - 1e-15 if s < 0 else 1e9)
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        vp = list(v)
        vp[i] = v[i] * (1 + s * mid)
        if thr_from(*vp) > mut:
            lo = mid
        else:
            hi = mid
    return s * 0.5 * (lo + hi)


def run_S8():
    head("S8.  WHICH constants are load-bearing: elasticity and flip budget "
         "(the real intervalisation priority)")

    print("\n(a) mutilde_* = t_max gam / (2 c_1 |u_perp| + 2 t_max c_1^2) with")
    print("    t_max = (-b + sqrt(b^2+m q_1))/m is an explicit function of the six")
    print("    constants.  gam enters LINEARLY and nowhere else, so")
    print("        d ln mutilde_* / d ln gam = 1   EXACTLY,")
    print("    hence any certified lower bound gam_cert simply rescales the")
    print("    threshold:   mutilde_*(cert) = mutilde_* * (gam_cert/gam).")
    Cb = None
    for (lam, mus, wc, r) in ANCHORS:
        Nc = max(1, int(np.floor(r * wc)))
        if Nc > N // 2:
            continue
        C = Cell(wc, Nc)
        ck(abs(elasticity(C, 0) - 1.0) < 1e-6,
           f"d ln mut_*/d ln gam = 1 at (lambda={lam},r={r})")
        ck(abs(thr_from(C.gam * 0.5, C.c1, C.up, C.q1, C.b, C.m)
               - 0.5 * C.thr) < 1e-12 * max(1.0, C.thr),
           f"mut_* is exactly linear in gam at (lambda={lam},r={r})")
        ck(abs(thr_from(*_consts(C)) - C.thr) < 1e-12 * max(1.0, C.thr),
           f"closed form reproduces Cell.thr at (lambda={lam},r={r})")
    print("    (verified at every anchor: elasticity 1 to 1e-6 and exact rescaling.)")

    print("\n(b) FLIP BUDGET: the smallest relative error in one constant (others")
    print("    exact) that would drag mutilde_* down to the anchor's own mutilde.")
    print(f"    {'lambda':>7} {'r':>3} {'wbar_c':>7} {'mutilde':>8} {'mut_*':>8} "
          + " ".join(f"{n:>10}" for n in CONST_NAMES))
    binding, bflip = None, None
    for (lam, mus, wc, r) in ANCHORS:
        mut = mus / lam
        Nc = max(1, int(np.floor(r * wc)))
        if Nc > N // 2:
            continue
        C = Cell(wc, Nc)
        fl = [flip_error(C, i, mut) for i in range(6)]
        ck(abs(fl[0] - (mut / C.thr - 1.0)) < 1e-6,
           f"gam flip = mutilde/mut_* - 1 at (lambda={lam},r={r})")
        if binding is None or abs(fl[0]) < abs(bflip[0]):
            binding, bflip = (lam, mus, wc, r, Nc, mut, C), fl
        print(f"    {lam:>7.2f} {r:>3d} {wc:>7.2f} {mut:>8.4f} {C.thr:>8.5f} "
              + " ".join(f"{100*x:>+9.1f}%" for x in fl))
    (lam, mus, wc, r, Nc, mut, C) = binding
    print(f"\n    BINDING anchor: lambda={lam}, mu*={mus}, r={r}, wbar_c={wc}, "
          f"N_c={Nc}, mutilde={mut}")
    print(f"    {'constant':>10} {'value':>12} {'elasticity':>12} "
          f"{'flip needed':>12} {'verdict':>20}")
    fl = bflip
    order = sorted(range(6), key=lambda i: abs(fl[i]))
    for i in order:
        a = abs(fl[i])
        vd = ("LOAD-BEARING" if a < 0.10 else
              ("moderate" if a < 0.20 else
               ("slack" if a < 0.50 else "COMFORTABLY SLACK")))
        print(f"    {CONST_NAMES[i]:>10} {_consts(C)[i]:>12.6f} "
              f"{elasticity(C, i):>+12.4f} {100*fl[i]:>+11.1f}% {vd:>20}")
    ck(set(order[:2]) == {0, 1},
       "gam and c_1 are the two most load-bearing constants")
    ck(order[-1] == 5, "m is the LEAST load-bearing of the six")
    ck(abs(fl[5]) > 5 * abs(fl[0]),
       "m needs > 5x the relative error of gam to flip the anchor")
    print("    ==> the note's long footnote about m having no cheap analytic bound is")
    print("    INVERSELY proportional to m's importance: m is the slackest of the six.")
    print("    Intervalise gam and c_1 FIRST (reviewer DOWNGRADE 5).")

    print("\n(c) what gam_cert/gam does each anchor need?  (elasticity 1 ==> the")
    print("    requirement is exactly mutilde/mutilde_*, and with the window inf of")
    print("    S3 it is mutilde/inf_W mutilde_*.)")
    print(f"    {'lambda':>7} {'r':>3} {'wbar_c':>7} {'gam':>9} "
          f"{'need at wbar_c':>15} {'need with inf_W':>16}")
    need_max = 0.0
    n_skip = 0
    for (lam2, mus2, wc2, r2) in ANCHORS:
        mut2 = mus2 / lam2
        Nc2 = max(1, int(np.floor(r2 * wc2)))
        if Nc2 > N // 2:
            n_skip += 1
            continue
        C2 = Cell(wc2, Nc2)
        inf_w = window_inf(Nc2, r2)[0]
        need_max = max(need_max, mut2 / inf_w)
        print(f"    {lam2:>7.2f} {r2:>3d} {wc2:>7.2f} {C2.gam:>9.5f} "
              f"{mut2/C2.thr:>15.4f} {mut2/inf_w:>16.4f}")
    # what the EXISTING delta = mu - alpha machinery delivers (published nodes,
    # `p4d_coulomb_certify.py` Stage C/D; the delta VALUES are quoted, the true
    # gam they are compared against is recomputed here)
    DELTA = ((1.0, 0.657276), (5.0, 1.803759))
    print("\n    against what the EXISTING gam-certificate machinery delivers")
    print("    (delta = mu - alpha at the published mutilde=0 nodes of")
    print("     `p4d_coulomb_certify.py` Stage C/D -- delta QUOTED, gam recomputed):")
    ratios = []
    for (wb, dl) in DELTA:
        ee = np.linalg.eigvalsh(H_mat(N, wb, 0))
        gtrue = float(ee[-1] - ee[-2])
        ratios.append(dl / gtrue)
        ck(dl < gtrue, f"delta is a valid lower bound for gam at wbar={wb}")
        print(f"      wbar={wb:>5}:  delta = {dl:.6f}  vs true gam = {gtrue:.6f}  "
              f"==> gam_cert/gam = {dl/gtrue:.4f}")
    print("      and it DEGRADES with wbar (inf g_alt on [20,60] is 1.076 against a")
    print("      true gam of ~3.5 there), while the binding anchors sit at")
    print("      wbar = 13.5 and 51.2.")
    best_avail = min(ratios)
    print(f"\n    worst anchor requirement = {need_max:.4f};  the existing mechanism's")
    print(f"    tightness at its highest published node (wbar=5) = {best_avail:.4f}.")
    ck(need_max > 0.70, "even the light anchors need gam_cert/gam > 0.70")
    if n_skip == 0:
        ck(need_max > best_avail,
           "the existing delta machinery is NOT tight enough for the binding anchor")
        ck(need_max > 0.90, "the binding anchor needs gam_cert/gam > 0.90")
    else:
        print(f"    [N = {N} skips {n_skip} heavy anchor(s) (N_c > N/2), so need_max is")
        print("     UNDERSTATED here; the two budget assertions run only at PN >= 400,")
        print("     where the binding anchor (lambda=0.5, r=10, wbar=13.5) needs 0.927.]")
    print("    ==> Sec.8 item 1 ('reuse the delta/g_alt mechanism') does NOT close at")
    print("    the present precision: a gam envelope 2-3x tighter is needed, and that")
    print("    is new work, not mechanical labour.  (The extrapolation of delta to")
    print("    wbar = 13.5/51.2 is the reviewer's PLAUSIBLE rating, not a recomputed")
    print("    certificate.)  Note also that this budget does NOT yet subtract the")
    print("    within-window inf (up to 9.2%, S3) or the N -> infinity drift (up to")
    print("    0.8%, S6) -- those eat into the SAME margin.")


# ==================================================================== verdict
def run_verdict():
    head("VERDICT")
    print("""
PROVED (analytic argument, constants in float64 -- see H1):
  * THEOREM.  Fix N_c and r, let W(N_c) = [N_c/r,(N_c+1)/r).  If
        mutilde < inf_{wbar in W(N_c)} mutilde_*(wbar,N_c),
        mutilde_* = t_max gam / (2 c_1 |u_perp| + 2 t_max c_1^2),
    then T |--> lam_1(M(T)) is strictly decreasing on the whole window, so
    lambda lam_1(M(T)) = 1 has at most one solution inside it.  Proof = Lemma 1
    (rank-one secular reduction under c_1 != 0 AND kappa = 2 mut c_1^2/gam < 1:
    lam_1 simple and in [mu_1 - 2 mut c_1^2, mu_1), v = (A-lam)^{-1}u_c exactly,
    mu_1-lam <= 2 mut c_1^2 -- NO hypothesis on c_2 is needed)
    + Lemma 2 (angle bound, also under kappa < 1)
    + Lemma 3 (robust Cor A': q_1 - 2 b t - m t^2 > 0 for t < t_max, with q_1 > 0
      supplied by Lemma A of the mu*=0 theory since psi_1 in C)
    + Lemma 4 (Hellmann-Feynman, legitimate because lam_1(M) is simple).
    The hypothesis mutilde < mutilde_* is equivalent to
    2 mut c_1 |u_perp| < t_max (gam - 2 mut c_1^2) and therefore IMPLIES
    kappa < 1, so Lemmas 1 and 2 are applicable and the assembly is not circular.
    MIND THE QUANTIFIER: the inf over the window, not the nodal value -- mut_*
    dips up to 9.2% inside a window (S3).
  * All eleven physical anchors of `P4D_coulomb_operator.md` satisfy the
    WINDOW-INF hypothesis, including the lambda = 0.4, r = 3 corner
    (mutilde = 0.325 < inf_W mutilde_* = 0.4027).  So the ENTIRE physical
    mutilde range [0, 0.33] is covered at the physical (wbar_c, N_c) pairs.
    Tightest margin: lambda = 0.5, r = 10 (inf_W mut_*/mutilde = 1.079).
  * GENERALIZED LEMMA B is true with the explicit epsilon
        mutilde <= mu_Metz(wbar,N_c) = min_{L != L', Wc > 0} W^F_{L,L'}/Wc_{L,L'},
    where Wc >= 0 has the closed form Wc_{L,L'} = 2 S(min(L',N_c)) (u_{L-1}
    1_{L-1<N_c} - u_L 1_{L<N_c}) supported on the N_c rows L <= N_c;  but
    mu_Metz ~ wbar^2/N_c^3 is 2-4 orders below physical, so this route is
    abandoned.

REFUTED:
  * K_c = O[1_{k<N_c}]  (the natural guess).  Exact identity instead:
    K_c = P_{N_c}(O[1]+I)P_{N_c}  -- a two-line consequence of O[1] = 2uu^T - I
    (S0c), numerics being bookkeeping only.
  * The 'at most one sign change' oscillation cone: the top eigenvector has
    exactly TWO sign changes, N-independently.

OPEN / NOT CLAIMED  (H1-H6, same numbering as the note's Sec.7):
  (H1) interval enclosures for gam, q_1, b, m, c_1, |u_perp| (float64 today).
       The reported mutilde_* are OVERESTIMATES (mut_*(N) falls monotonically in
       N, S6) and are 2-3 digit converged on the load-bearing anchors, not 6.
       S8: gam (elasticity exactly 1) and c_1 (-1.47..-1.72) are load-bearing,
       m is comfortably slack (needs +101% to flip) -- intervalise gam and c_1
       first.  The existing delta/g_alt gam machinery (0.919 at wbar=1, 0.819 at
       wbar=5, degrading) is NOT tight enough for the binding anchor's 0.927.
  (H2) uniqueness ACROSS N_c steps: genuinely false for r = 3 at mutilde = 0.10
       on explicit lambda-intervals; for r = 10 all steps are upward on
       [0.34, 8.2] (S4b/S4d) -- so uniqueness there is for every lambda WHOSE
       ROOT LIES IN THAT RANGE, not for every lambda: below wbar = 0.34 there
       are three real downward steps (S4c), the worst being N_c 0 -> 1.
  (H3) the band mutilde_* < mutilde < mutilde_crit (factor 1.4-3, global min
       1.4025 over 224 cells) is numerically
       fine but not proved; the counterexample (N_c=4, wbar~0.74, mutilde in
       [0.322,0.334]) lies in it and is a genuine failure of monotonicity, so
       the band cannot be filled in trivially.  Correctly excluded (S5).
  (H4) upstream dependency AND one ingredient added here: lam_1(A) SIMPLE, which
       `P4_uniqueness_strategy.md` explicitly declines to provide (it needs
       neither simplicity nor Krein-Rutman, proving only E_1 (cap) C != {0}).
       Supplied in S0(e) from the STRICT POSITIVITY of A's off-diagonal entries,
       i.e. by Perron-Frobenius -- this theorem RE-INTRODUCES the positivity the
       mu*=0 theorem deliberately avoided.
  (H5) N -> infinity tail: the limits of q_1, b, m are not certified, and the
       Perron-Frobenius step of (H4) is a finite-matrix statement.
  (H6) EINSTEIN-ONLY: the six constants are only ever evaluated on
       A = K_ph(wbar) = H(wbar), a single-atom Einstein kernel, whereas the
       mu*=0 theorem being extended holds for every P in Pcal.  The sibling
       probe finds a genuine within-window non-monotonicity at the PHYSICAL
       mutilde = 0.26 for a 2-atom P (0.8 delta(w-0.3s) + 0.2 delta(w-4.0s),
       N_c in [5,12], lam_1 down 1.8e-3 over s in (0.45,1.33)).
""")


STAGES = {"S0": run_S0, "S1": run_S1, "S2": run_S2, "S3": run_S3,
          "S4": run_S4, "S5": run_S5, "S6": run_S6, "S7": run_S7,
          "S8": run_S8}

if __name__ == "__main__":
    want = [s.upper() for s in sys.argv[1:]] or list(STAGES)
    t0 = time.time()
    print(f"p4d_coulomb_uniqueness.py   N = {N}   stages = {want}")
    for s in want:
        STAGES[s]()
    if set(want) >= {"S2", "S3", "S5"}:
        run_verdict()
    print(f"\nALL {PASS[0]} checks passed  ({time.time()-t0:.1f}s).")
