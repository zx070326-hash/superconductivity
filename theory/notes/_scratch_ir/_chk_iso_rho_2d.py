"""2D dependence table  (Om_ac_max, w_B) -> alpha_budget interval.

Why 2D and not two marginalised 1D legs (team-lead, and I agree):
beta_lo is ATTAINED and climbs monotonically with w_B, and its two entry points
BOTH carry w_B --
    beta_lo = 1/(1 + (w_A Om^2)/(w_B G^2)),   G = (wlog / Om^{w_A})^{1/w_B}
so "Om_ac_max unknown" and "w_B unknown" live on ONE constraint surface.
Marginalising each into a +-interval and summing would over- or under-state the
width depending on which side of the surface the truth sits.

Structure worth reading off: the UPPER end s_hi = t^(2 w_B) depends on w_B ALONE
(beta <= 1 is free), so only the LOWER end is two-dimensional.

Run: OMP_NUM_THREADS=1 python _chk_iso_rho_2d.py
"""
import numpy as np

M = {"H": 1.00782503, "D": 2.01410178}
T = np.sqrt(M["H"] / M["D"])
WLOG = 405.0
SHIFT = 0.5 / np.log(M["D"] / M["H"])
A_BASE = -0.3992
EXP = (-0.2898, -0.3223)          # isotope_curated.md refs [7], [8]


def beta_lo(om, w_B):
    w_A = 1.0 - w_B
    G = (WLOG / om**w_A) ** (1.0 / w_B)
    return 1.0 / (1.0 + (w_A * om**2) / (w_B * G**2))


def alpha_of(w_B, beta):
    s = T ** (2.0 * w_B) / (1.0 - beta * (1.0 - T**2))
    return A_BASE + SHIFT * np.log(s)


OMS = (150.0, 190.0, 230.0, 270.0, 310.0)
WBS = (0.750, 0.790, 0.830, 0.870, 0.9023)

print("=== alpha_budget interval  [lower, upper]  on the (Om_ac_max, w_B) surface ===")
print("  upper end depends on w_B ALONE (beta<=1 free); only the lower end is 2D.")
print()
hdr = f"{'Om\\w_B':>8} |" + "".join(f"{w:^19.4f}" for w in WBS)
print(hdr)
print("-" * len(hdr))
for om in OMS:
    row = f"{om:8.0f} |"
    for w_B in WBS:
        lo = alpha_of(w_B, beta_lo(om, w_B))
        hi = alpha_of(w_B, 1.0)
        row += f"  [{lo:+.3f},{hi:+.3f}]  "
    print(row)

print()
print("=== which cells COVER an experimental value? ===")
print(f"  experiment: {EXP[0]:+.4f} (ref 7)   {EXP[1]:+.4f} (ref 8)")
print()
hdr2 = f"{'Om\\w_B':>8} |" + "".join(f"{w:^12.4f}" for w in WBS)
print(hdr2)
print("-" * len(hdr2))
for om in OMS:
    row = f"{om:8.0f} |"
    for w_B in WBS:
        lo, hi = alpha_of(w_B, beta_lo(om, w_B)), alpha_of(w_B, 1.0)
        hits = [n for n, e in zip(("7", "8"), EXP) if lo <= e <= hi]
        row += f"{(','.join(hits) if hits else '-'):^12}"
    print(row)

print()
print("=== why marginalising is wrong: 2D envelope vs product of 1D legs ===")
lo2d = min(alpha_of(w, beta_lo(o, w)) for o in OMS for w in WBS)
hi2d = max(alpha_of(w, 1.0) for w in WBS)
# 1D legs, each taken with the OTHER input frozen at a mid value
om_mid, wb_mid = 230.0, 0.830
leg_om = [alpha_of(wb_mid, beta_lo(o, wb_mid)) for o in OMS]
leg_wb = [alpha_of(w, beta_lo(om_mid, w)) for w in WBS]
half_om = 0.5 * (max(leg_om) - min(leg_om))
half_wb = 0.5 * (max(leg_wb) - min(leg_wb))
centre = alpha_of(wb_mid, beta_lo(om_mid, wb_mid))
print(f"  true 2D envelope          : [{lo2d:+.4f}, {hi2d:+.4f}]   width {hi2d-lo2d:.4f}")
print(f"  1D legs about the midpoint: Om leg half-width {half_om:.4f}, "
      f"w_B leg half-width {half_wb:.4f}")
print(f"    linear sum   : {centre:+.4f} +- {half_om+half_wb:.4f} -> "
      f"[{centre-half_om-half_wb:+.4f}, {centre+half_om+half_wb:+.4f}]")
print(f"    quadrature   : {centre:+.4f} +- {np.hypot(half_om,half_wb):.4f} -> "
      f"[{centre-np.hypot(half_om,half_wb):+.4f}, {centre+np.hypot(half_om,half_wb):+.4f}]")
print("  ==> both 1D composites MISS the true envelope on the upper side")
print("      (the upper end is set by beta=1 at low w_B, a corner no 1D leg visits).")
