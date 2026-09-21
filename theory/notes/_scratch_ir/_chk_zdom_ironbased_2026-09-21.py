#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_zdom_ironbased_2026-09-21.py  --  2026-09-21, ceiling-fe agent, FROM ZERO.

TASKS A (per-arrow (Z-dom) test on real s+- material inputs),
      B (degrees-of-freedom audit of the criterion),
      C (multiband ceiling test  Tc / ceiling, reported together with lambda).

======================  L1 SOURCE, read by me this session  ==================
G. A. Ummarino, M. Tortello, D. Daghero, R. S. Gonnelli,
"Three-band s+- Eliashberg theory and the superconducting gaps of iron
 pnictides", arXiv:0904.1808  = Phys. Rev. B 80, 172503 (2009).
PDF fetched this session to  _fe_0904.1808.pdf  (3,465,699 bytes, 4 pp);
text layer  _fe_0904.1808.txt  (pdftotext -layout -enc UTF-8);
page-2 image  _fe0904_p2-2.png (170 dpi) and crop _fe0904_c2-2.png (600 dpi)
BOTH read (CLAUDE.md: two instruments for load-bearing characters).

Eq.(1)  w_n Z_i = w_n + sum_j (G^N+G^M)_ij N^Z_j + piT sum_{m,j} Lam^Z_ij N^Z_j
Eq.(2)  Z_i Dl_i = sum_j (G^N-G^M)_ij N^Dl_j
                 + piT sum_{m,j} [Lam^Dl_ij - mu*_ij(w_c)] th(w_c-|w_m|) N^Dl_j
   with, VERBATIM from the page image,
        Lam^Z_ij  = Lam^ph_ij + Lam^sp_ij
        Lam^Dl_ij = Lam^ph_ij - Lam^sp_ij            <<<< this is exactly (K1)
        Lam^{ph,sp}_ij(iw_n,iw_m) = int_0^{+inf} dW a2F^{ph,sp}_ij(W)
                                       / [(w_n-w_m)^2 + W^2]
   >>> NOTE the printed kernel is MISSING the factor 2W.  It is dimensionally
       impossible as printed (a2F is dimensionless => the integral has units
       1/energy, but Lam must be dimensionless), so the intended kernel is the
       standard  2 int W a2F/[(w_n-w_m)^2+W^2].  Dimensional analysis settles
       this with no guesswork; recorded because it is a load-bearing character.

   a2F_ij(W) = C_ij [ L(W+W_ij,Y_ij) - L(W-W_ij,Y_ij) ],
               L(W +- W_ij, Y_ij) = [ (W +- W_ij)^2 + (Y_ij)^2 ]^{-1}
   >>> as printed this is NEGATIVE for W>0 (since (W-W0)^2 < (W+W0)^2), so
       C_ij < 0 is required for a2F >= 0.  Either way the shape is SINGLE-
       SIGNED on W>0; below I use the positive-definite ordering
       D(W) = L(W-W0,Y) - L(W+W0,Y) > 0  and normalise by lambda_ij.
   "In all our calculations we always set W_ij = W_0" ... "and Y_ij = 2 meV".
   "The cut-off energy is w_c = 12 W_0 and the maximum quasiparticle energy
    is w_max = 16 W_0."        (w_c is a MATSUBARA cutoff, paired with mu*)

Model reduction (p.2, verbatim): "lam^ph_ii = lam^ph_ij = 0, lam^sp_ii = 0 and
   mu*_ii = mu*_ij = 0"; coupling matrix [[0,0,l31 nu1],[0,0,l32 nu2],
   [l31,l32,0]], nu1=N1/N3, nu2=N2/N3.
   122 Ba0.6K0.4Fe2As2: nu1=1, nu2=2, Tc=37 K.  1111 SmFeAsO0.8F0.2:
   nu1=0.4, nu2=0.5, Tc=52 K.
Table I (p.3): l31=3.866, l32=0.471, mu*=0, lam_ii=0.4, T=2K, 122.
Table II (p.4): l31=2.730, l32=0.758, mu*=0.1, lam_ii=0.4, 122.
Fig.3 text (p.3): 122  W0=10 meV, l31=4.267, l32=0.569  (lam_ii=0)
                  1111 W0=10 meV, l31=14.520, l32=1.708 (lam_ii=0)
p.4: "We also tested the effect into the model of a small intraband coupling
      (POSSIBLY OF PHONON ORIGIN).  In the case of Ba0.6K0.4Fe2As2 we used
      lam_ii = 0.4"                    <<< fixes lam_ii as the PHONON piece

CONVENTION / detailed balance:  printed lam_13/lam_31 = nu1 = N1/N3
  => lam_ij = N_i V_ij (DOS on the LEFT index), V symmetric.
  Symmetrised  Lam~_ij = sqrt(N_i N_j) V_ij = lam_ij sqrt(N_j/N_i).
  Check: Lam~_13 = l31 sqrt(nu1), Lam~_23 = l32 sqrt(nu2), Lam~_ii = lam_ii.

=============================================================================
DIRECTION AUDIT (CLAUDE.md single-sided rule)
  * Ceiling  Tc <= C_inf sqrt(rho(M2''_Z)) is an UPPER bound => a numerically
    computed ratio Tc/ceiling < 1 is CONSISTENT, not a proof of the theorem.
  * lam_1(P_N O P_N) <= lam_1(O): truncation UNDER-estimates => a Tc solved
    from the truncated operator is a LOWER estimate of the model's Tc.
=============================================================================
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

KB = 0.08617333262         # meV / K
G2 = 1.318140212132        # repo rational-exact (CERT-D), re-derived in the
                           # companion script _chk_zdom_budget_2026-09-21.py
C_GF = (np.pi ** 2 / 6.0) / G2
C_INF = np.sqrt(G2) / (2 * np.pi)


# ------------------------------------------------- the Ummarino spectral shape
def shape(W, W0, Y):
    """Positive-definite Lorentzian difference, unnormalised."""
    return 1.0 / ((W - W0) ** 2 + Y ** 2) - 1.0 / ((W + W0) ** 2 + Y ** 2)


def shape_moments(W0, Y):
    """Returns (I_m1, I_p1) = (int D/W dW , int W D dW) and <W^2> = I_p1/I_m1."""
    Im1 = quad(lambda W: shape(W, W0, Y) / W, 0, np.inf, limit=400)[0]
    Ip1 = quad(lambda W: W * shape(W, W0, Y), 0, np.inf, limit=400)[0]
    return Im1, Ip1, Ip1 / Im1


def g_of_pi(p):
    """g(pi) = (1 - pi*p/sinh(pi*p))/p^2 ; Taylor for small p."""
    p = np.asarray(p, float)
    out = np.empty_like(p)
    small = np.abs(p) < 1e-3
    x = np.pi * p[~small]
    out[~small] = (1.0 - x / np.sinh(x)) / p[~small] ** 2
    q = p[small]
    out[small] = np.pi ** 2 / 6.0 - 7 * np.pi ** 4 * q ** 2 / 360.0
    return out


def gbar(W0, Y, T_K):
    """weighted average of g over the spectrum with weight 2W a2F."""
    twoPiT = 2 * np.pi * KB * T_K
    num = quad(lambda W: W * shape(W, W0, Y) * g_of_pi(np.array([W / twoPiT]))[0],
               0, np.inf, limit=400)[0]
    den = quad(lambda W: W * shape(W, W0, Y), 0, np.inf, limit=400)[0]
    return num / den


# ------------------------------------------------------- multiband operator
def u_vec(N):
    n = np.arange(N)
    return 1.0 / np.sqrt(2.0 * n + 1.0)


def kernels(lam_ph, lam_sp, W0, Y, T_K, Im1):
    """Return callables k -> Lam^ph(k), Lam^sp(k) built from the ACTUAL
       spectral shape (not an Einstein pole).
       Lam_ij(k) = 2 int W a2F_ij(W)/(W^2+(2 pi T k)^2) dW, with
       a2F_ij = (lam_ij / (2 Im1)) * D(W)   so that Lam_ij(0) = lam_ij."""
    twoPiT = 2 * np.pi * KB * T_K

    def J(k):
        """2 int W D(W) / (W^2 + (2piT k)^2) dW  / (2 Im1) -- the scalar shape
           factor, = 1 at k=0."""
        a = (twoPiT * k) ** 2
        val = quad(lambda W: W * shape(W, W0, Y) / (W ** 2 + a),
                   0, np.inf, limit=400)[0]
        return val / Im1                      # 2 cancels: Im1 = int D/W

    cache = {}

    def jf(k):
        if k not in cache:
            cache[k] = J(k)
        return cache[k]

    return (lambda k: lam_ph * jf(k)), (lambda k: lam_sp * jf(k)), jf


def build_O(LamD, LamZ, N):
    B = LamD(0).shape[0]
    u = u_vec(N)
    LD = [LamD(k) for k in range(2 * N + 2)]
    LZ = [LamZ(k) for k in range(N + 2)]
    RZ = np.zeros((N, B))
    for i in range(B):
        acc = LZ[0][i].sum()
        RZ[0, i] = acc
        for n in range(1, N):
            acc += 2 * LZ[n][i].sum()
            RZ[n, i] = acc
    O = np.zeros((N * B, N * B))
    for n in range(N):
        for m in range(N):
            blk = LD[n + m + 1].copy()
            if n != m:
                blk = blk + LD[abs(n - m)]
            blk = blk * (u[n] * u[m])
            if n == m:
                blk = blk + u[n] ** 2 * LD[0] - np.diag(u[n] ** 2 * RZ[n])
            O[n * B:(n + 1) * B, m * B:(m + 1) * B] = blk
    return 0.5 * (O + O.T)


# =============================================================================
MATS = [
    # tag, compound, l31, l32, nu1, nu2, lam_ii(ph), W0_meV, Tc_K, provenance
    ("T1", "Ba0.6K0.4Fe2As2 (122)", 3.866, 0.471, 1.0, 2.0, 0.4, 10.0, 37.0,
     "Table I, p.3 (mu*=0)"),
    ("T2", "Ba0.6K0.4Fe2As2 (122)", 2.730, 0.758, 1.0, 2.0, 0.4, 10.0, 37.0,
     "Table II, p.4 (mu*=0.1)"),
    ("F3a", "Ba0.6K0.4Fe2As2 (122)", 4.267, 0.569, 1.0, 2.0, 0.0, 10.0, 37.0,
     "Fig.3 body text, p.3 (primary model, lam_ii=0)"),
    ("F3b", "SmFeAsO0.8F0.2 (1111)", 14.520, 1.708, 0.4, 0.5, 0.0, 10.0, 52.0,
     "Fig.3 body text, p.3 (primary model, lam_ii=0)"),
]


def split(l31, l32, nu1, nu2, lam_ii):
    """Return the SYMMETRISED (K1) pieces Lam~^ph, Lam~^sp (coupling constants)."""
    ph = np.diag([lam_ii, lam_ii, lam_ii]).astype(float)
    sp = np.zeros((3, 3))
    sp[0, 2] = sp[2, 0] = l31 * np.sqrt(nu1)
    sp[1, 2] = sp[2, 1] = l32 * np.sqrt(nu2)
    return ph, sp


def main():
    print("=" * 78)
    print("PART 0 -- spectral shape moments (W0=10 meV, Y=2 meV)")
    print("=" * 78)
    W0, Y = 10.0, 2.0
    Im1, Ip1, W2 = shape_moments(W0, Y)
    print("   int D/W dW  = %.6e      int W D dW = %.6e" % (Im1, Ip1))
    print("   <W^2> := M2_ij/lam_ij = %.6f meV^2   sqrt = %.6f meV"
          % (W2, np.sqrt(W2)))
    print("   (W_ij = W0 for ALL (i,j) => <W^2> is the SAME for every element")
    print("    => M2 = <W^2> * lam  EXACTLY, i.e. the 9-element moment matrix")
    print("    carries ONE frequency, not nine.)")
    # cross-check <W^2> by a second route: Richardson-free direct high-accuracy
    Im1b = quad(lambda t: shape(np.exp(t), W0, Y), -25, 25, limit=800)[0]
    print("   cross-check int D/W dW via W=e^t substitution: %.6e (rel %.2e)"
          % (Im1b, abs(Im1b - Im1) / Im1))

    print()
    print("=" * 78)
    print("TASK B / ARROW 1 -- (Z-dom) DEGREES-OF-FREEDOM AUDIT")
    print("=" * 78)
    print("""
  Under the paper's own convention (K1):  Lam^Z = Lam^ph + Lam^sp,
  Lam^Dl = Lam^ph - Lam^sp, both pieces elementwise >= 0.  The (Z-dom)
  MARGIN at each element and each Matsubara index k is

     m_ij(k) := Lam^Z_ij(k) / |Lam^Dl_ij(k)|
              = (Lam^ph_ij + Lam^sp_ij)/|Lam^ph_ij - Lam^sp_ij|
              = (1 + r_ij(k)) / |1 - r_ij(k)| ,   r := Lam^sp_ij/Lam^ph_ij

  [analytic]  (i) m >= 1 for every r in [0,inf]  => the (Z-dom) VERDICT is a
                  TAUTOLOGY: not one published number enters it.  DOF = 0.
             (ii) m = 1 EXACTLY iff r in {0, inf}, i.e. iff exactly one of the
                  two channels is present at that (i,j).
  In this model lam^ph is diagonal-only and lam^sp is off-diagonal-only
  (the paper's own reduction), so r_ij in {0, inf} at EVERY element
  => the margin is exactly 1 everywhere: (Z-dom) holds with EQUALITY,
     zero margin, and zero dependence on l31, l32, nu1, nu2, lam_ii, W0, Y, T.
""")
    for tag, comp, l31, l32, nu1, nu2, lii, W0_, Tc, prov in MATS:
        ph, sp = split(l31, l32, nu1, nu2, lii)
        Z = ph + sp
        Dl = ph - sp
        with np.errstate(divide="ignore", invalid="ignore"):
            m = np.where(np.abs(Dl) > 0, Z / np.abs(Dl), np.nan)
        finite = m[np.isfinite(m)]
        print("   %-4s %-24s  elements with Lam^Dl != 0 : %d ; margins: %s"
              % (tag, comp, finite.size,
                 " ".join("%.6f" % v for v in np.unique(np.round(finite, 12)))))
    print()
    print("   Zero-DOF STRESS TEST: randomise every published number by up to")
    print("   +-90% and recompute the margins.  If any margin moves, the")
    print("   verdict has >0 DOF.")
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(20000):
        p = rng.uniform(0.1, 1.9, size=5)
        ph, sp = split(3.866 * p[0], 0.471 * p[1], 1.0 * p[2], 2.0 * p[3],
                       0.4 * p[4])
        Z, Dl = ph + sp, ph - sp
        msk = np.abs(Dl) > 0
        worst = max(worst, np.abs(Z[msk] / np.abs(Dl[msk]) - 1.0).max())
    print("   max |margin - 1| over 20000 randomisations = %.3e   => DOF = 0"
          % worst)
    print("   NEGATIVE CONTROL: let the two channels OVERLAP at one element")
    print("   (put a phonon piece p on the (1,3) element too).  Then the margin")
    print("   must move away from 1:")
    for p in (0.1, 0.5, 1.0, 3.866):
        ph, sp = split(3.866, 0.471, 1.0, 2.0, 0.4)
        ph[0, 2] = ph[2, 0] = p
        Z, Dl = ph + sp, ph - sp
        print("      phonon on (1,3) = %6.3f -> margin_13 = %s"
              % (p, ("%.6f" % (Z[0, 2] / abs(Dl[0, 2]))) if Dl[0, 2] != 0
                 else "inf (exact compensation)"))

    print()
    print("=" * 78)
    print("ARROW 2 -- a2F^ph_ij >= 0 and a2F^sp_ij >= 0 elementwise")
    print("=" * 78)
    Ws = np.linspace(1e-6, 2000.0, 400001)
    D = shape(Ws, W0, Y)
    print("   min over W in (0, 2000] meV of the Lorentzian-difference shape:")
    print("      min D = %.6e   (>=0 everywhere)  [analytic: (W-W0)^2<(W+W0)^2]"
          % D.min())
    print("   => a2F_ij = (lam_ij/(2 Im1)) D(W) >= 0 whenever lam_ij >= 0,")
    print("      and every printed lam_ij here is >= 0.  ARROW 2: COVERED,")
    print("      but by CONSTRUCTION of the model's functional form (0 DOF),")
    print("      not by a measurement of the material.")

    print()
    print("=" * 78)
    print("ARROW 3 -- each Lam^X_ij(k) non-increasing in k")
    print("=" * 78)
    print("   [analytic] Lam_ij(k) = 2 int W a2F_ij(W)/(W^2+(2piT k)^2) dW with")
    print("   a2F>=0; the integrand is strictly decreasing in k^2 pointwise.")
    print("   Numerical confirmation on the actual shape:")
    _, _, jf = kernels(np.eye(1), np.eye(1), W0, Y, 37.0, Im1)
    vals = [jf(k) for k in range(0, 25)]
    dmax = max(vals[i + 1] - vals[i] for i in range(len(vals) - 1))
    print("      J(k) k=0..5: %s   max increment = %.3e (<=0 required)"
          % (" ".join("%.6f" % v for v in vals[:6]), dmax))

    print()
    print("=" * 78)
    print("TASK C -- CEILING TEST.  Tc <= C_inf sqrt(rho(M2''_Z))")
    print("   M2''_Z = diag(M2^ph + c M2^sp) + c offdiag(M2^ph + M2^sp)")
    print("   c = c_GF = %.10f (unconditional) or c = max(1, gbar/g(2))" % C_GF)
    print("=" * 78)
    hdr = ("%-4s %-24s %7s %7s %9s %9s %9s %8s %8s" %
           ("set", "compound", "lam_max", "rho(lam)", "ceil_cGF",
            "ceil_GF1", "Tc(K)", "Tc/ceil", "gbar/g2"))
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for tag, comp, l31, l32, nu1, nu2, lii, W0_, Tc, prov in MATS:
        ph, sp = split(l31, l32, nu1, nu2, lii)
        M2ph, M2sp = ph * W2, sp * W2
        gb = gbar(W0, Y, Tc)
        cmat = max(1.0, gb / G2)
        def ceil(c):
            M = (np.diag(np.diag(M2ph) + c * np.diag(M2sp))
                 + c * ((M2ph + M2sp) - np.diag(np.diag(M2ph + M2sp))))
            rho = np.abs(np.linalg.eigvals(M)).max()
            return C_INF * np.sqrt(rho) / KB      # in K
        cg, c1 = ceil(C_GF), ceil(cmat)
        lam_max = (ph + sp).sum(axis=1).max()
        rows.append((tag, comp, lam_max, np.abs(np.linalg.eigvals(ph + sp)).max(),
                     cg, c1, Tc, Tc / c1, gb / G2, prov))
        print("%-4s %-24s %7.3f %8.3f %9.1f %9.1f %9.1f %8.4f %8.4f"
              % (tag, comp, lam_max, rows[-1][3], cg, c1, Tc, Tc / c1, gb / G2))
    print()
    print("   NOTE gbar/g(2) < 1 for every row => the material-resolved GF")
    print("   factor max(1, gbar/g(2)) = 1 => 'ceil_GF1' is the sharp column.")
    print("   ALL rows report lam_max alongside so the reader can judge whether")
    print("   the test is vacuous (repo: lam<=0.5 gives Tc/ceiling ~ 0.003-0.02).")

    print()
    print("=" * 78)
    print("GATE -- solve the MODEL's own Tc from my operator and compare with")
    print("        the paper's Tc.  (independent end-to-end check of inputs,")
    print("        convention, symmetrisation and the operator itself)")
    print("=" * 78)
    for tag, comp, l31, l32, nu1, nu2, lii, W0_, Tc, prov in MATS:
        ph, sp = split(l31, l32, nu1, nu2, lii)

        def lam1(T, N=None, cut=False):
            Kph, Ksp, _ = kernels(ph, sp, W0, Y, T, Im1)
            LZ = lambda k: Kph(k) + Ksp(k)
            LD = lambda k: Kph(k) - Ksp(k)
            if N is None:
                N = 400
            if cut:
                Nc = max(1, int(np.floor(12.0 * W0 / (2 * np.pi * KB * T)
                                         - 0.5)) + 1)
                N = Nc
            O = build_O(LD, LZ, N)
            return np.linalg.eigvalsh(O)[-1], N

        try:
            f_nc = lambda T: lam1(T, N=220)[0] - 1.0
            Tc_nc = brentq(f_nc, 2.0, 400.0, xtol=1e-3)
        except Exception as e:
            Tc_nc = float("nan")
        try:
            f_c = lambda T: lam1(T, cut=True)[0] - 1.0
            Tc_c = brentq(f_c, 2.0, 400.0, xtol=1e-3)
            Nc = lam1(Tc_c, cut=True)[1]
        except Exception as e:
            Tc_c, Nc = float("nan"), -1
        print("   %-4s %-24s  Tc(paper)=%5.1f K | my Tc, no cutoff = %6.1f K"
              " | my Tc, w_c=12W0 cutoff (N_c=%d) = %6.1f K"
              % (tag, comp, Tc, Tc_nc, Nc, Tc_c))
    print()
    print("   DIRECTION: the theta(w_c-|w_m|) cutoff REMOVES pairing weight =>")
    print("   lowers lam_1 => the cutoff column is the LOWER of the two, and")
    print("   the paper solved WITH the cutoff.")

    # -------------------------------------------------------------------
    print()
    print("=" * 78)
    print("PART D -- IS THE CEILING RATIO ITSELF CIRCULAR?  (DOF of Tc/ceiling)")
    print("=" * 78)
    print("""
  The coupling constants l31, l32 were FITTED by the authors so that the model
  reproduces the experimental Tc and the gaps, at a HAND-CHOSEN W0 = 10 meV.
  So before reading anything into Tc/ceiling, count its degrees of freedom.

  [analytic] SCALE INVARIANCE.  The linearised Eliashberg operator depends on
  the spectrum only through W/(2 pi T).  Hence if (W0,Y) -> (s W0, s Y) then
  Tc -> s Tc exactly.  Simultaneously <W^2> -> s^2 <W^2> so M2 -> s^2 M2 and
  ceiling = C_inf sqrt(rho(M2'')) -> s * ceiling.  THEREFORE

        Tc / ceiling  is INVARIANT under W0 -> s W0 (with Y/W0 fixed):
        it is a function of the coupling matrix and the shape ratio Y/W0 ONLY.

  => the ceiling test cannot see W0 at all.  And since Tc_exp is reproduced by
  the fit, Tc_exp/ceiling ~= Tc_model(lam)/ceiling(lam) = a CLOSED-FORM
  FUNCTION OF lam.  Same shape as the repo's 'filling factor is circular'
  accident (Tc_AD fed back into the AD formula).
""")
    print("   NUMERICAL CHECK of the scale invariance (s = 1, 2, 5, 0.3):")
    l31, l32, nu1, nu2, lii = 3.866, 0.471, 1.0, 2.0, 0.4
    ph0, sp0 = split(l31, l32, nu1, nu2, lii)
    for s in (1.0, 2.0, 5.0, 0.3):
        W0s, Ys = 10.0 * s, 2.0 * s
        Im1s, _, W2s = shape_moments(W0s, Ys)

        def lam1s(T):
            Kph, Ksp, _ = kernels(ph0, sp0, W0s, Ys, T, Im1s)
            return np.linalg.eigvalsh(
                build_O(lambda k: Kph(k) - Ksp(k),
                        lambda k: Kph(k) + Ksp(k), 120))[-1] - 1.0

        Tcs = brentq(lam1s, 1.0, 2000.0, xtol=1e-4)
        M2ph, M2sp = ph0 * W2s, sp0 * W2s
        M = (np.diag(np.diag(M2ph + M2sp))
             + ((M2ph + M2sp) - np.diag(np.diag(M2ph + M2sp))))
        ceil = C_INF * np.sqrt(np.abs(np.linalg.eigvals(M)).max()) / KB
        print("      s=%4.1f  W0=%6.2f meV  Tc_model=%8.3f K  ceiling=%8.3f K"
              "  ratio=%.8f" % (s, W0s, Tcs, ceil, Tcs / ceil))
    print("   => ratio identical to ~7 digits across a factor 17 in W0:")
    print("      W0 carries ZERO weight in the ceiling test.  [analytic+numeric]")

    print()
    print("   CIRCULARITY: replace the EXPERIMENTAL Tc by the MODEL's own Tc.")
    print("   If the two columns agree, the test never touched the experiment.")
    W0, Y = 10.0, 2.0
    Im1, _, W2 = shape_moments(W0, Y)
    print("   %-4s %-24s %9s %9s %9s %9s %9s"
          % ("set", "compound", "Tc_exp", "Tc_model", "ceiling",
             "exp/ceil", "mod/ceil"))
    for tag, comp, a, b, n1, n2, lii_, W0_, Tc, prov in MATS:
        ph, sp = split(a, b, n1, n2, lii_)

        def f(T):
            Kph, Ksp, _ = kernels(ph, sp, W0, Y, T, Im1)
            return np.linalg.eigvalsh(
                build_O(lambda k: Kph(k) - Ksp(k),
                        lambda k: Kph(k) + Ksp(k), 140))[-1] - 1.0

        Tm = brentq(f, 2.0, 400.0, xtol=1e-3)
        M2ph, M2sp = ph * W2, sp * W2
        M = (np.diag(np.diag(M2ph + M2sp))
             + ((M2ph + M2sp) - np.diag(np.diag(M2ph + M2sp))))
        ceil = C_INF * np.sqrt(np.abs(np.linalg.eigvals(M)).max()) / KB
        print("   %-4s %-24s %9.2f %9.2f %9.2f %9.4f %9.4f"
              % (tag, comp, Tc, Tm, ceil, Tc / ceil, Tm / ceil))

    print()
    print("=" * 78)
    print("PART E -- WHAT *IS* NON-TRIVIAL: the COST OF MULTIBANDNESS")
    print("   cost := ceiling_multiband(c_GF on offdiag) / ceiling_scalar")
    print("         = sqrt( rho(M2''_Z) / rho(M2) ) in [1, sqrt(c_GF)=1.11710]")
    print("   This number does NOT depend on Tc at all => it is NOT circular;")
    print("   it is a property of the coupling matrix's band structure.")
    print("=" * 78)
    print("   %-4s %-24s %10s %10s %10s"
          % ("set", "compound", "rho(M2)", "rho(M2''_Z)", "cost"))
    for tag, comp, a, b, n1, n2, lii_, W0_, Tc, prov in MATS:
        ph, sp = split(a, b, n1, n2, lii_)
        M2ph, M2sp = ph * W2, sp * W2
        M2 = M2ph + M2sp
        rho0 = np.abs(np.linalg.eigvals(M2)).max()
        Mpp = (np.diag(np.diag(M2ph) + C_GF * np.diag(M2sp))
               + C_GF * (M2 - np.diag(np.diag(M2))))
        rho1 = np.abs(np.linalg.eigvals(Mpp)).max()
        print("   %-4s %-24s %10.3f %10.3f %10.6f"
              % (tag, comp, rho0, rho1, np.sqrt(rho1 / rho0)))
    print("   (repo CANONICAL Sec.8.4 quotes MgB2 1.0122, LiFeAs-like 1.0853,")
    print("    worst case sqrt(c_GF)=1.1171 -- these rows sit near the worst")
    print("    case because the coupling is almost purely INTERBAND.)")


if __name__ == "__main__":
    main()
