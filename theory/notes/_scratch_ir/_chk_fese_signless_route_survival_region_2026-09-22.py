"""
team-lead, 2026-09-22.  Independent re-verification of essenberger's
_chk_fese_signless_amplitude_test_2026-09-22.py  (verdict: UNDETERMINED),
plus the one thing it did not compute: the exact survival region of the
sign-agnostic route in terms of the UNPRINTED SF diagonals.

Inputs (I read them myself from a fresh 130 dpi render of 1411.2121 p.3,
_fe1411_tl_p3-3.png; all match the agent's transcription):
    phonon  lam_++=0.028  lam_--=0.029  lam_+-=0.042  lam_-+=0.048   (Fig.3a labels)
    SF      lam_max = 0.07 / 0.16 / 0.37 / 0.48  at alpha/alpha_c = 0.90/0.92/0.94/0.98
            (Fig.3c inset labels);  text: "In the limit alpha->alpha_c ... since (in
            this limit) SF are dominated by off diagonal components".

Facts used (analytic):
  (F1) For a 2x2 entrywise-nonnegative matrix [[a,b],[c,d]],
         lam_max = (a+d)/2 + sqrt(((a-d)/2)^2 + bc)   <=>   bc = (lam_max - a)(lam_max - d).
       So lam_max alone bounds bc ABOVE (bc <= lam_max^2, at a=d=0) but not below
       (bc -> 0 as a -> lam_max).  This is the agent's counterexample [[0,100],[0,0]]
       said exactly: lam_max = 0 there because it is nilpotent, and the formula gives
       bc = (0-0)(0-0) = 0 consistently -- the entry 100 is invisible to lam_max.
  (F2) The consumption point of condition (ii) is the SYMMETRISED kernel
       (w_i = sqrt(N_i)), whose off-diagonal is the geometric mean sqrt(lam_ij lam_ji)
       (detailed balance N_i lam_ij = N_j lam_ji).  So the sign-agnostic magnitude
       route needs   sqrt(lam^ph_+- lam^ph_-+)  >=  sqrt(b c)   i.e.
                     lam^ph_+- * lam^ph_-+  >=  (lam_max - a)(lam_max - d).
       Its FAILURE is a sound FAIL of the sign-agnostic ROUTE (a sufficient condition
       for (ii)), NOT of (ii) itself: (ii) can still hold if the SF off-diagonal has
       the right sign -- which is Eq.(37)'s content again.
"""
from fractions import Fraction as Fr
from math import sqrt

ph = {"++": Fr(28, 1000), "--": Fr(29, 1000), "+-": Fr(42, 1000), "-+": Fr(48, 1000)}
lam_max = {"0.90": Fr(7, 100), "0.92": Fr(16, 100), "0.94": Fr(37, 100), "0.98": Fr(48, 100)}

T = ph["+-"] * ph["-+"]          # exact: 0.002016
print(f"threshold  lam^ph_+- * lam^ph_-+ = {T}  = {float(T):.6f}   (sqrt = {sqrt(T):.5f})")
print(f"detailed-balance ratio r = lam_-+/lam_+- = {ph['-+']/ph['+-']}  (= N-ratio if both channels obey it)\n")

for conv, fac in (("as printed", Fr(1)), ("/1.5", Fr(2, 3)), ("x1.5", Fr(3, 2))):
    print(f"--- convention {conv} ---")
    for key, lm in lam_max.items():
        L = lm * fac
        # survival needs (L-a)(L-d) <= T.  Symmetric-diagonal case a=d: a >= L - sqrt(T)
        a_sym = float(L) - sqrt(T)
        # prose case a=d=0 (off-diagonal dominated): bc = L^2
        ratio_prose = float(L) / sqrt(T)
        print(f"  alpha/alpha_c={key}: lam_max={float(L):.4f} | survive only if (L-a)(L-d) <= {float(T):.6f}"
              f" | a=d case: diag >= {max(a_sym,0):.4f} = {100*max(a_sym,0)/float(L):.1f}% of lam_max"
              f" | prose (a=d=0): sqrt(bc)/sqrt(ph) = {ratio_prose:.2f}x -> {'FAIL' if ratio_prose>1 else 'pass'}")
    print()

print("Reading: at alpha/alpha_c = 0.98 the sign-agnostic route survives only if the SF matrix is")
print(">= 90.6% diagonal (a=d case) -- i.e. the OPPOSITE of the paper's own 'dominated by off")
print("diagonal components'.  Under that prose the route fails by ~10.7x (as printed) and by")
print("7.1x / 16x under the two 1.5x conventions.  At alpha/alpha_c = 0.90 (lam_max = 0.07) the")
print("route is genuinely undetermined: a 36% diagonal share already saves it.")
