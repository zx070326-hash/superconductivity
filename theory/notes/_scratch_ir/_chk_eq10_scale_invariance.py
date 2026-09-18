#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_eq10_scale_invariance.py -- 2026-09-18, main session, round 12, from zero.

source-hunt's round 12 brought Ummarino-Gonnelli-Massidda-Bianconi
cond-mat/0310284v3 Eq.(10), which prints the WHOLE mu* matrix and settles the
2.23/2.48 reading in favour of (A).  I re-pulled the PDF and read the text
layer myself (source-hunt read page images -- two independent extraction
paths).  My L1, `pdftotext -layout -enc UTF-8`, p.3:

    mu*(p) = | mu*_ss  mu*_sp |
             | mu*_ps  mu*_pp |

           = mu(w_c) NN_tot(EF,p) x | 2.23/NN_s(EF,p)   1/NN_s(EF,p)  |     (10)
                                    | 1/NN_p(EF,p)      2.48/NN_p(EF,p) |

  -> ROW 1 carries NN_s in BOTH denominators, ROW 2 carries NN_p in both.
     That is per-row normalisation = reading (A), printed as a matrix.

  p.4 verbatim: "The numbers 2.23 and 2.48 in the Coulomb matrix have been
  calculated in the MgB2 case and we suppose that, as first approximation,
  they do not depend on x ... In conclusion we put the cut-off energy
  w_c = 700 meV.  So we have now only one free parameter, mu(w_c,x) that must
  be determined by reproducing exactly the experimental critical temperatures."
  p.5 verbatim: "we find mu(w_c, x = 0) = 0.0326".

  p.2 verbatim: "as shown in ref. 14, lam_ij(x)/lam_ji(x) = mu*_ij(x)/mu*_ji(x)
  = NN_j(w=EF,x)/NN_i(w=EF,x)"      [14] = Adrian et al.

** THE POINT THIS SCRIPT MAKES, WHICH IS A CORRECTION TO MY OWN SS10.1. **
The criterion  mu*_ss mu*_pp >= mu*_sp mu*_ps  applied to Eq.(10) is
INVARIANT under the overall scale mu(w_c) NN_tot AND under both DOS values:
every one of them cancels.  The margin is exactly 2.23 x 2.48.
=> the criterion never touches the FITTED content of this mu*.  It tests the
   first-principles ratio pair from refs [15]=Golubov, [16]=Brinkman -- i.e.
   exactly the part I had already declared uninformative (Bochner).
=> "fitted mu*_ij is PSD" is NOT closed by MgB2.  Re-opened.
"""

from fractions import Fraction as F

# ---- Golubov cond-mat/0111262v2 p.5, four printed numbers (my earlier L1) ----
MU_SS, MU_SP = F("0.210"), F("0.095")
MU_PS, MU_PP = F("0.069"), F("0.172")
N_S, N_P = F("0.300"), F("0.410")

# ---- Ummarino cond-mat/0310284v3 Eq.(10), the two printed numbers ----
A, B = F("2.23"), F("2.48")
MU_WC_X0 = F("0.0326")          # p.5, for x = 0


def margin(m11, m12, m21, m22):
    return (m11 * m22) / (m12 * m21)


def main():
    print("=" * 74)
    print("1.  margin from Eq.(10) -- exact rational, DOS-free, scale-free")
    print("=" * 74)
    # mu*_ij = c * P_ij with c = mu(w_c) NN_tot and
    #   P = [[A/N_s, 1/N_s], [1/N_p, B/N_p]]
    for c in (F(1), F(7, 3), MU_WC_X0, F(10**6)):
        for ns, npi in ((N_S, N_P), (F(1), F(1)), (F("0.17"), F("2.9"))):
            P11, P12 = c * A / ns, c * F(1) / ns
            P21, P22 = c * F(1) / npi, c * B / npi
            m = margin(P11, P12, P21, P22)
            assert m == A * B, (c, ns, npi, m)
    print("    margin = A*B exactly, for every (c, N_s, N_p) tried")
    print("    A*B = %s = %.6f" % (A * B, float(A * B)))
    print("    => the ONE fitted number mu(w_c) cancels; both DOS cancel.")
    print("    => the criterion is blind to everything this paper fitted.")

    print("\n" + "=" * 74)
    print("2.  margin from Golubov's four printed numbers (independent path)")
    print("=" * 74)
    mg = margin(MU_SS, MU_SP, MU_PS, MU_PP)
    print("    margin = %s = %.6f" % (mg, float(mg)))
    print("    Eq.(10) form gives %.6f   ->  gap %.3f %%"
          % (float(A * B), 100 * (float(mg) / float(A * B) - 1)))

    print("\n" + "=" * 74)
    print("3.  can ONE scale reproduce all four printed numbers?  (coherence)")
    print("=" * 74)
    # fix c from mu*_ss, then predict the other three
    c = MU_SS * N_S / A
    pred = {"mu*_ss": c * A / N_S, "mu*_sp": c / N_S,
            "mu*_ps": c / N_P, "mu*_pp": c * B / N_P}
    obs = {"mu*_ss": MU_SS, "mu*_sp": MU_SP, "mu*_ps": MU_PS, "mu*_pp": MU_PP}
    print("    c := mu*_ss N_s / 2.23 = %.6f" % float(c))
    worst = 0.0
    for k in ("mu*_ss", "mu*_sp", "mu*_ps", "mu*_pp"):
        d = 100 * (float(pred[k]) / float(obs[k]) - 1)
        worst = max(worst, abs(d))
        print("      %-7s predicted %.6f   printed %.3f   %+.2f %%"
              % (k, float(pred[k]), float(obs[k]), d))
    print("    worst deviation %.2f %%  (mu* printed to 3 decimals"
          " -> rounding alone is up to ~0.7 %%)" % worst)
    print("    => Golubov's FOUR numbers carry, to <1 %%, only THREE pieces of")
    print("       information: 2.23, 2.48 and one scale.  Reading (A) is not")
    print("       merely the better fit -- it is the generating form.")

    print("\n" + "=" * 74)
    print("4.  the two readings against Eq.(10), for the record")
    print("=" * 74)
    print("      (A) row-wise        : mu*_ss/mu*_sp = %.4f  vs 2.23  (%+.2f %%)"
          % (float(MU_SS / MU_SP), 100 * (float(MU_SS / MU_SP) / 2.23 - 1)))
    print("                            mu*_pp/mu*_ps = %.4f  vs 2.48  (%+.2f %%)"
          % (float(MU_PP / MU_PS), 100 * (float(MU_PP / MU_PS) / 2.48 - 1)))
    print("      (B) single ref      : mu*_pp/mu*_sp = %.4f  vs 2.48  (%+.2f %%)"
          % (float(MU_PP / MU_SP), 100 * (float(MU_PP / MU_SP) / 2.48 - 1)))
    print("    Eq.(10) prints the matrix, so this is now a consistency check,")
    print("    not the evidence.  The evidence is the shared row denominator.")

    print("\n" + "=" * 74)
    print("5.  NOT load-bearing: does the scale match mu(w_c,x=0)=0.0326 ?")
    print("=" * 74)
    ntot = c / MU_WC_X0
    print("    c / mu(w_c) = NN_tot would be %.4f states/cell.eV" % float(ntot))
    print("    N_s + N_p   = %.4f" % float(N_S + N_P))
    print("    ratio %.3f -- NOT a contradiction: Golubov's mu* and Ummarino's"
          % float(ntot / (N_S + N_P)))
    print("    w_c=700 meV need not be the same cutoff, and NN_tot need not be")
    print("    N_s+N_p.  UNCHECKED, flagged, used for nothing.")


if __name__ == "__main__":
    main()
