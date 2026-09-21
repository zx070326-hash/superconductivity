"""BUDGET PROBE (not the job): intrinsic slack of an analytic two-sided bound on

    s_partsub := rho_H / rho_D   under strictly HARMONIC partial substitution
                 (P-weights frozen, substituted group B rigidly scaled w -> t w,
                  unsubstituted group A untouched),   t = sqrt(M_H/M_D).

Consumption point (team-lead): the term 0.7221*ln(s) in alpha_budget.
Required precision: s to about +-0.05 (0.025 in alpha vs 0.08-0.11 to adjudicate).
Making it tighter than that buys nothing.

Question this script answers: BEFORE doing the work, how loose is the estimator
itself -- the slack I cannot remove by trying harder?

Run: OMP_NUM_THREADS=1 python _chk_iso_rho_budget.py
"""
import numpy as np

M = {"H": 1.00782503, "D": 2.01410178, "T": 3.01604928}
t_HD = np.sqrt(M["H"] / M["D"])
WLOG_H = 405.0          # Errea 2013 Table I, SSCHA, cm^-1  (L1)
WOP = {"H": 0.750, "D": 0.783, "T": 0.792}   # lambda share of H-character optical block


def rho_of_P(w, p):
    w, p = np.asarray(w, float), np.asarray(p, float)
    assert abs(p.sum() - 1.0) < 1e-12
    return float(p @ w**2) / np.exp(2.0 * float(p @ np.log(w)))


# ------------------------------------------------------------------ step 1
print("=== 1. Is there an estimator at all?  Candidate CLOSED FORM ===")
print("""
  Under frozen weights + rigid scaling of B:
     <w^2>_D = a2 + t^2 b2      (a2 := \\int_A w^2 P,  b2 := \\int_B w^2 P)
     ln wlog_H - ln wlog_D = w_B ln(1/t)      EXACTLY  (only B's logs shift)
  ==>  s = (<w^2>_H/<w^2>_D) * t^(2 w_B)
         = t^(2 w_B) / (1 - beta (1 - t^2)),    beta := b2/(a2+b2)
  i.e. s is an IDENTITY in (w_B, beta, t) -- no inequality, no two-block needed.
""")


def s_closed(w_B, beta, t):
    return t ** (2.0 * w_B) / (1.0 - beta * (1.0 - t**2))


# verify the identity against brute force on GENERAL (non-two-block) spectra
rng = np.random.default_rng(7)
worst = 0.0
for _ in range(3000):
    nA, nB = rng.integers(1, 6), rng.integers(1, 6)
    wA = rng.uniform(40.0, 260.0, nA)          # unsubstituted (Pd-character)
    wB = rng.uniform(200.0, 1200.0, nB)        # substituted (H-character)
    p = rng.dirichlet(np.ones(nA + nB))
    pA, pB = p[:nA], p[nA:]
    t = rng.choice([t_HD, np.sqrt(M["D"] / M["T"]), 0.5, 0.9])
    s_bf = rho_of_P(np.r_[wA, wB], p) / rho_of_P(np.r_[wA, t * wB], p)
    a2, b2 = float(pA @ wA**2), float(pB @ wB**2)
    s_cf = s_closed(pB.sum(), b2 / (a2 + b2), t)
    worst = max(worst, abs(s_cf / s_bf - 1.0))
print(f"  identity vs brute force, 3000 random GENERAL spectra: max rel dev = {worst:.2e}")
print("  ==> ESTIMATOR SLACK = 0.  It is an identity, not a bound.")
print("      All remaining width comes from not knowing (w_B, beta).")

# ------------------------------------------------------------------ step 2
print()
print("=== 2. Sensitivity: how well must (w_B, beta) be known? ===")
w_B = WOP["H"]
beta0 = 0.98415       # from the two-block Omega_ac=130 point, for reference
print(f"  d ln s / d beta = (1-t^2)/(1-beta(1-t^2)) "
      f"= {(1-t_HD**2)/(1-beta0*(1-t_HD**2)):.4f}   (at beta={beta0})")
print(f"  d ln s / d w_B  = 2 ln t = {2*np.log(t_HD):.4f}")
need = 0.05
print(f"  ==> for s to +-{need}: need beta to +-{need/ (s_closed(w_B,beta0,t_HD)*(1-t_HD**2)/(1-beta0*(1-t_HD**2))):.3f}"
      f"   and  w_B to +-{need/(s_closed(w_B,beta0,t_HD)*abs(2*np.log(t_HD))):.3f}")

# ------------------------------------------------------------------ step 3
print()
print("=== 3. The FREE endpoint: beta <= 1 is automatic. What does it alone give? ===")
for wb in (0.75, 0.79, 0.90):
    lo, hi = s_closed(wb, 0.0, t_HD), s_closed(wb, 1.0, t_HD)
    print(f"  w_B={wb:.2f}:  beta in [0,1] -> s in [{lo:.4f}, {hi:.4f}]  (width {hi-lo:.4f})")
print("  ==> the free upper end s <= t^(-2 w_A) is TIGHT-ish; the free lower end is useless.")
print("      So the whole job reduces to: LOWER-bound beta.")

# ------------------------------------------------------------------ step 4
print()
print("=== 4. Lower-bounding beta with data already in hand (no alpha^2F, no lit) ===")
print("""  a2 = w_A <w^2>_A <= w_A * Om_ac_max^2                       [needs Pd band top]
  b2 = w_B <w^2>_B >= w_B * G_op^2         (Jensen: <w^2> >= exp 2<ln w>)
  G_op = (wlog / G_ac^{w_A})^{1/w_B},  G_ac <= Om_ac_max            [same input]
  ==> beta >= 1/(1 + (w_A Om_ac_max^2)/(w_B G_op_min^2))
""")
print(f"{'Om_ac_max':>10} | {'G_op_min':>9} {'beta_lo':>8} | {'s_lo':>7} {'s_hi':>7} {'width':>7} "
      f"| {'alpha_lo':>9} {'alpha_hi':>9}")
w_A = 1.0 - w_B
SHIFT = 0.5 / np.log(M["D"] / M["H"])
for om in (150.0, 190.0, 220.0, 260.0, 300.0):
    G_op_min = (WLOG_H / om**w_A) ** (1.0 / w_B)
    beta_lo = 1.0 / (1.0 + (w_A * om**2) / (w_B * G_op_min**2))
    s_lo, s_hi = s_closed(w_B, beta_lo, t_HD), s_closed(w_B, 1.0, t_HD)
    print(f"{om:10.0f} | {G_op_min:9.1f} {beta_lo:8.5f} | {s_lo:7.4f} {s_hi:7.4f} "
          f"{s_hi-s_lo:7.4f} | {-0.3992+SHIFT*np.log(s_lo):+9.4f} {-0.3992+SHIFT*np.log(s_hi):+9.4f}")

print()
print("  same, at w_B = 0.79 (top of Errea's 75-79% quote):")
w_B2, w_A2 = 0.79, 0.21
for om in (190.0, 220.0, 260.0):
    G_op_min = (WLOG_H / om**w_A2) ** (1.0 / w_B2)
    beta_lo = 1.0 / (1.0 + (w_A2 * om**2) / (w_B2 * G_op_min**2))
    s_lo, s_hi = s_closed(w_B2, beta_lo, t_HD), s_closed(w_B2, 1.0, t_HD)
    print(f"{om:10.0f} | {G_op_min:9.1f} {beta_lo:8.5f} | {s_lo:7.4f} {s_hi:7.4f} {s_hi-s_lo:7.4f}")

# ------------------------------------------------------------------ step 5
print()
print("=== 5. Cross-check: does the bound contain the two-block values 1.14-1.18? ===")
for om_ac, wop_h in ((100.0, 0.750), (130.0, 0.750), (160.0, 0.750), (190.0, 0.750)):
    Wop_H = float(np.exp((np.log(WLOG_H) - (1 - wop_h) * np.log(om_ac)) / wop_h))
    a2, b2 = (1 - wop_h) * om_ac**2, wop_h * Wop_H**2
    beta = b2 / (a2 + b2)
    print(f"  Om_ac={om_ac:5.0f}  Om_op={Wop_H:6.1f}  beta={beta:.5f}  "
          f"s={s_closed(wop_h, beta, t_HD):.4f}")

# ------------------------------------------------------------------ step 6
print()
print("=== 6. w_B is NOT free either: the two harmonic steps pin it ===")
HARM = {"H": 205.0, "D": 150.0, "T": 125.0}
for A, B in (("H", "D"), ("D", "T")):
    print(f"  {A}->{B}: w_B = 2 ln({HARM[A]/HARM[B]:.4f}) / ln({M[B]/M[A]:.4f}) = "
          f"{2*np.log(HARM[A]/HARM[B])/np.log(M[B]/M[A]):.4f}")
print("  (SSCHA row quotes 0.75-0.79; the harmonic row reads out ~0.90 -- an open 0.11 tension)")
