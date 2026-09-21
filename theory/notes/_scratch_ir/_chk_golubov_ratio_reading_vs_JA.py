#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_golubov_ratio_reading_vs_JA.py -- 2026-09-18, main session, round 11, from zero.

source-hunt's round 11 puts a THIRD-PARTY bare-Coulomb ratio next to Golubov's
printed "2.23/2.48/1" and observes that the two competing readings of that
sentence disagree by 27% in opposite places.  It explicitly declined to
adjudicate.  This script adjudicates.

THE TWO READINGS of  "The sigsig-, pipi- and sigpi-values were in the ratio
2.23/2.48/1"  (Golubov cond-mat/0111262v2 printed p.5):

  (A) ROW-WISE  (my round-9 verdict):  2.23 = mu_ss/mu_sp ,  2.48 = mu_pp/mu_ps
  (B) SINGLE REFERENCE (literal grammar): both referred to the same mu_sp.

L1 anchors used, all read by me from the text layer in earlier rounds:
  Golubov p.4 : N_s = 0.300, N_p = 0.410 states/cell*eV
  Golubov p.5 : mu*_ss=0.210, mu*_sp=0.095, mu*_ps=0.069, mu*_pp=0.172
  Golubov p.5 : "This allows one to express mu*_ij(omega_c) by these ratios
                 and one single free parameter"
  Mazin-Antropov p.11 : "<<V_C>>_sigsig : <<V_C>>_pipi : <<V_C>>_sigpi
                 to be approximately 3 : 2.5 : 1"  [Jepsen-Andersen, TB-LMTO]
  Mazin-Antropov p.11 : mu_ij "is formally defined as <<V_C>>_ij N_j"
  Mazin-Antropov p.11 : "assuming mu_ss = mu_pp = alpha mu_sp ... alpha ~ 2.3"

THE DECISIVE TEST.  Reading (B) can only survive if the printed mu* is NOT
elementwise proportional to mu -- i.e. if a matrix Morel-Anderson
renormalisation sits between them.  So: CAN the MA map carry a bare V with
Jepsen-Andersen ratios onto the V* ratios that the printed mu* encodes?
If no, (B) is dead on the paper's own numbers and (A) stands.
"""

import numpy as np
from itertools import product

N_S, N_P = 0.300, 0.410
MUSTAR = np.array([[0.210, 0.095],
                   [0.069, 0.172]])
NB = np.array([N_S, N_P])


def K_of(M):
    """strip the mu_ij = N_j K_ij convention -> the symmetric kernel K."""
    return M / NB[None, :]


def ratios(K):
    """(K_ss/K_sp, K_pp/K_sp) using the SYMMETRISED off-diagonal."""
    off = 0.5 * (K[0, 1] + K[1, 0])
    return K[0, 0] / off, K[1, 1] / off


def ma_map(V, v, DL):
    """bare symmetric V (scale v) -> renormalised symmetric kernel K*.
       mu^s = D_N^{1/2} (v V) D_N^{1/2};  (mu^s)* = ((mu^s)^{-1} + D_L)^{-1};
       K* = D_N^{-1/2} (mu^s)* D_N^{-1/2}."""
    Dh = np.diag(np.sqrt(NB))
    Dhi = np.diag(1.0 / np.sqrt(NB))
    mus = Dh @ (v * V) @ Dh
    star = np.linalg.inv(np.linalg.inv(mus) + DL)
    return Dhi @ star @ Dhi


def main():
    print("=" * 76)
    print("0.  what the PRINTED mu* encodes, convention stripped")
    print("=" * 76)
    K = K_of(MUSTAR)
    print("    K = mu*_ij / N_j  = [[%.5f, %.5f], [%.5f, %.5f]]"
          % (K[0, 0], K[0, 1], K[1, 0], K[1, 1]))
    print("    symmetry check |K_sp - K_ps| / K_sp = %.4f %%"
          % (100 * abs(K[0, 1] - K[1, 0]) / K[0, 1]))
    tgt = ratios(K)
    print("    => V*-ratios encoded by the printed mu* : %.4f : %.4f : 1"
          % tgt)

    print("\n" + "=" * 76)
    print("1.  the two readings, taken at face value (mu* proportional to mu)")
    print("=" * 76)
    print("    the sentence says mu* is fixed by 'these ratios and ONE single")
    print("    free parameter' => mu*_ij = p * (ratio pattern) => the printed")
    print("    mu* must itself display the quoted ratios.  Test both readings:")
    rA = (MUSTAR[0, 0] / MUSTAR[0, 1], MUSTAR[1, 1] / MUSTAR[1, 0])
    rB = (MUSTAR[0, 0] / MUSTAR[0, 1], MUSTAR[1, 1] / MUSTAR[0, 1])
    for tag, (x, y) in (("(A) row-wise       ", rA),
                        ("(B) single reference", rB)):
        print("      %s : printed gives %.4f / %.4f  vs quoted 2.23 / 2.48"
              "   -> %+.2f %% / %+.2f %%"
              % (tag, x, y, 100 * (x / 2.23 - 1), 100 * (y / 2.48 - 1)))
    print("    => on the paper's OWN four numbers, (A) matches to <1 %, (B)")
    print("       misses mu*_pp by 27 %.  This is the round-9 verdict, and")
    print("       source-hunt independently reproduced it (its item 1).")

    print("\n" + "=" * 76)
    print("2.  the Jepsen-Andersen cross-check, at the V level")
    print("=" * 76)
    f = N_S / N_P
    print("    mu_ij = <<V_C>>_ij N_j  =>  mu_ss/mu_sp = (N_s/N_p) V_ss/V_sp")
    print("                                mu_pp/mu_ps = (N_p/N_s) V_pp/V_sp")
    print("                                mu_pp/mu_sp =           V_pp/V_sp")
    print("    N_s/N_p = %.5f" % f)
    VA = (2.23 / f, 2.48 / (1 / f))
    VB = (2.23 / f, 2.48)
    for tag, V in (("(A) row-wise        ", VA),
                   ("(B) single reference", VB)):
        print("      %s -> V_ss:V_pp:V_sp = %.4f : %.4f : 1   vs JA 3 : 2.5 : 1"
              "   (%+.1f %% / %+.1f %%)"
              % (tag, V[0], V[1], 100 * (V[0] / 3 - 1), 100 * (V[1] / 2.5 - 1)))
    print("    => the discriminant is V_pp/V_sp alone, exactly as source-hunt")
    print("       says.  (B) agrees with JA to 0.8 %; (A) misses by 27.3 %.")

    print("\n    consistency of JA with Mazin-Antropov's own alpha ~ 2.3:")
    print("      single-reference at the mu level: mu_ss/mu_sp = %.4f,"
          % (f * 3))
    print("      mu_pp/mu_sp = %.4f  -> common approximation alpha ~ %.2f"
          % (2.5, (f * 3 + 2.5) / 2))
    print("      MA writes 'mu_ss = mu_pp = alpha mu_sp' -- literally single")
    print("      reference -- and quotes alpha ~ 2.3.  So MA's own arithmetic")
    print("      IS reading (B).  [this is MA's convention, not Golubov's]")

    print("\n" + "=" * 76)
    print("3.  ** THE DECISIVE TEST **  can a matrix MA renormalisation carry")
    print("    a Jepsen-Andersen bare V onto the printed mu* ?")
    print("=" * 76)
    print("""    Reading (B) survives ONLY if mu* is not proportional to mu, i.e.
    only if the MA map sits in between.  Target V*-ratios %.4f : %.4f : 1.
    Scan the bare V over BOTH candidate starting points, the scale v over 6
    decades, and D_L over scalar AND band-indexed diagonals.""" % tgt)
    starts = {
        "JA  3 : 2.5 : 1": np.array([[3.0, 1.0], [1.0, 2.5]]),
        "(B) 3.048 : 2.48 : 1": np.array([[VB[0], 1.0], [1.0, VB[1]]]),
        "(A) 3.048 : 1.815 : 1": np.array([[VA[0], 1.0], [1.0, VA[1]]]),
    }
    vs = np.geomspace(1e-3, 1e3, 241)
    ls = np.geomspace(1e-3, 1e3, 121)
    for tag, V in starts.items():
        best = (np.inf, None)
        for v in vs:
            for l1, l2 in product(ls, ls):
                DL = np.diag([l1, l2])
                try:
                    Ks = ma_map(V, v, DL)
                except np.linalg.LinAlgError:
                    continue
                if Ks[0, 0] <= 0 or Ks[1, 1] <= 0 or Ks[0, 1] <= 0:
                    continue
                r = ratios(Ks)
                err = max(abs(r[0] / tgt[0] - 1), abs(r[1] / tgt[1] - 1))
                if err < best[0]:
                    best = (err, (v, l1, l2, r))
        v, l1, l2, r = best[1]
        print("    start %-22s best mismatch %7.3f %%   at v=%.3g L=(%.3g,%.3g)"
              % (tag, 100 * best[0], v, l1, l2))
        print("      %-22s reached V*-ratios %.4f : %.4f : 1"
              % ("", r[0], r[1]))

    print("""
    NOTE on the scan: D_L >= 0 is the PHYSICAL branch (omega_c below the
    bandwidth).  The scan above allows any positive (l1,l2), i.e. it is
    GENEROUS to reading (B): if even the unconstrained scan cannot reach the
    target, the constrained one certainly cannot.  Direction audit: this is a
    best-case search, so a large residual is a genuine NO, while a small
    residual would only mean 'not excluded'.""")

    print("\n" + "=" * 76)
    print("4.  why the map cannot get there -- the structural reason")
    print("=" * 76)
    print("""    (mu^s)* = ((mu^s)^{-1} + D_L)^{-1} suppresses the LARGER
    eigenvalue of mu^s more than the smaller one (x -> x/(1+Lx) for scalar L,
    and monotonically in the resolvent order generally).  Starting from a V
    whose two diagonal entries are 3 and 2.5 (nearly isotropic), the map
    cannot INCREASE the sigma/pi diagonal asymmetry from 3/2.5 = 1.20 to the
    target 3.032/1.817 = 1.67 while holding V_ss/V_sp fixed at ~3.03.""")
    for tag, V in (("JA 3:2.5:1", np.array([[3.0, 1.0], [1.0, 2.5]])),):
        print("    diagonal asymmetry V_ss/V_pp along the map (v=1):")
        for L in (0.0, 0.1, 0.3, 1.0, 3.0, 10.0):
            Ks = ma_map(V, 1.0, L * np.eye(2))
            print("      L = %5.2f :  V*_ss/V*_pp = %.4f   (target %.4f)"
                  % (L, Ks[0, 0] / Ks[1, 1], K[0, 0] / K[1, 1]))


if __name__ == "__main__":
    main()
