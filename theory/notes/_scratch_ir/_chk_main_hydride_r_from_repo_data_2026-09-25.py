"""Hydride r = Tc_exp / ceiling from numbers already in the repo (2026-09-25).

ceiling[K] = C_inf * sqrt(M2),  M2 = lambda * (w2bar)^2,  C_inf = 0.1827262477 (analytic).
Inputs are transcribed from (a) deficit_budget_h3s_2026-09-17.md table (cross-check of the
unit chain against the note's own ceiling column) and (b) omega2_data_sources_2026-09-17.md
Part 2 (Ding et al. arXiv:2507.03383 Table I, page-image reading recorded there).
Rows in (b) are NOT re-read from the PDF by this script; see the note for what is/isn't verified.
Run: OMP_NUM_THREADS=1 python this_file.py
"""
import math

CINF = 0.1827262477
MEV_TO_K = 11.604518  # e/k_B in K per meV

# (a) H3S/D3S: lambda, sqrt<w^2> [meV], ceiling column of the note [K], Tc_exp [K]
h3s = [
    ("H3S-200 anh", 1.721, 149.26, 415.2, 203),
    ("D3S-200 anh", 1.784, 105.43, 298.6, 150),
    ("H3S-160 anh", 1.876, 137.87, 400.4, 203),
]
print("== (a) unit-chain cross-check against note's ceiling column ==")
for name, lam, w_meV, ceil_note, tc_exp in h3s:
    ceil = CINF * math.sqrt(lam) * w_meV * MEV_TO_K
    print(f"{name}: ceiling={ceil:.2f} K (note {ceil_note}; rel diff {ceil/ceil_note-1:+.2e}) "
          f" r=Tc_exp/ceiling={tc_exp/ceil:.4f}")

# (b) Ding et al Table I: lambda(0), w2bar [K], Tc [K], Tc with renormalised mu*(N) [K] or None
yh = [
    ("YH6 DFPT",  2.534, 1399, 270, 248),
    ("YH6 SSCHA", 1.779, 1662, 245, 218),
    ("YH6 SPIA",  1.712, 1738, 212, None),
    ("YH9 DFPT",  3.038, 1308, 252, 236),
    ("YH9 SPIA",  1.749, 1905, 240, None),
]
print("\n== (b) YH6/YH9 from Ding Table I (Tc here is the CALCULATED Tc, mu*=0.11, not experiment) ==")
for name, lam, w2, tc, tcr in yh:
    ceil = CINF * math.sqrt(lam) * w2
    s = f"{name}: ceiling={ceil:.1f} K  Tc_calc/ceiling={tc/ceil:.4f}"
    if tcr:
        s += f"  (renorm {tcr/ceil:.4f})"
    print(s)

# M2 invariance across methods for YH6 (Hopfield sum rule check on printed numbers)
print("\nYH6 ceiling spread across DFPT/SSCHA/SPIA:")
c = [CINF*math.sqrt(l)*w for _, l, w, _, _ in yh[:3]]
print([f"{x:.1f}" for x in c], f"spread={(max(c)-min(c))/min(c):.3%}")
