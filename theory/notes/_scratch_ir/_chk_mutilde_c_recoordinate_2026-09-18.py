"""
Pure arithmetic (no operator work): re-express the repo's measured weak-coupling
extremality boundary in the variable  mutilde = mu*(omega_c)/lambda.

Source of the input numbers (5 rows, r = omega_c/omega_rms = 5):
  theory/notes/coulomb/mustar_extremality_failure_map_2026-09-17.md, section 4.1
  columns: mu*_rms | bare mu | s_c = lambda_c/mu*_rms bracket | lambda_c
Convention of that note: mu*(omega) = mu / (1 + mu ln(omega_c/omega)),
so the "bare mu" column IS mu*(omega_c), and the operator uses
mutilde = mu*(omega_c)/lambda.   (This is checked below, column mu*_rms_check.)
"""
from fractions import Fraction
import math

ROWS = [  # (mu*_rms, bare mu, s_c_lo, s_c_hi, lambda_c)  -- transcribed from the note
    (0.10,   0.1192, 2.1500, 2.1625, 0.2156),
    (0.13,   0.1644, 2.2500, 2.2625, 0.2933),
    (0.15,   0.1977, 2.3250, 2.3375, 0.3497),
    (0.18,   0.2534, 2.4406, 2.4625, 0.4413),
    (0.2023, 0.3000, 2.5250, 2.5500, 0.5133),
]
r = 5.0
lnr = math.log(r)

print(f"{'mu*_rms':>8} {'mu(wc)':>8} {'mu*_rms_check':>14} {'s_c':>8} "
      f"{'lambda_c':>9} {'mutilde_c':>10} {'mutilde_c bracket':>22}")
mts = []
for mus_rms, mu, slo, shi, lc in ROWS:
    check = mu / (1 + mu * lnr)
    mt = mu / lc
    # bracket from the s_c bracket: lambda_c = s_c * mu*_rms
    lo, hi = slo * mus_rms, shi * mus_rms
    mt_hi, mt_lo = mu / lo, mu / hi          # mutilde decreasing in lambda_c
    mts.append(mt)
    print(f"{mus_rms:8.4f} {mu:8.4f} {check:14.6f} "
          f"{(slo+shi)/2:8.4f} {lc:9.4f} {mt:10.4f} "
          f"[{mt_lo:.4f}, {mt_hi:.4f}]")

print()
print(f"  spread of lambda_c   : {max(l for *_ , l in ROWS)/min(l for *_, l in ROWS):.3f}x "
      f"(+{100*(max(l for *_, l in ROWS)/min(l for *_, l in ROWS)-1):.0f}%)")
scs = [(a + b) / 2 for _, _, a, b, _ in ROWS]
print(f"  spread of s_c        : {max(scs)/min(scs):.3f}x (+{100*(max(scs)/min(scs)-1):.1f}%)")
print(f"  spread of mutilde_c  : {max(mts)/min(mts):.3f}x (+{100*(max(mts)/min(mts)-1):.1f}%)")
print()
print("  leading-order prediction  lambda < 2 mu*_rms  <=>  mutilde > (1 + mu ln r)/2 :")
print(f"  {'mu(wc)':>8} {'mutilde_c LO':>13} {'mutilde_c measured':>19} {'meas/LO':>9}")
for (mus_rms, mu, slo, shi, lc), mt in zip(ROWS, mts):
    lo_pred = (1 + mu * lnr) / 2
    print(f"  {mu:8.4f} {lo_pred:13.4f} {mt:19.4f} {mt/lo_pred:9.3f}")
