"""rho := <w^2>/wlog^2 under isotope substitution.

Three questions, all on data already in hand (no alpha^2F, no literature):

  (A) Does rho cancel exactly under a UNIFORM dilation of alpha^2F?  (team-lead's lemma)
  (B) Is uniform dilation actually the right hypothesis for PdH?      (only H is substituted)
  (C) Decompose the two-block s = rho_H/rho_D into
      [partial-substitution, present already harmonically] x [anharmonic shape change].

Conventions (Allen-Dynes):
  lam        = 2 \int a2F/w dw
  P(dw)     := (2/lam) (a2F/w) dw          -- a probability measure
  <w^2>      = \int w^2 P(dw)
  ln wlog    = \int ln w P(dw)
  rho        = <w^2> / wlog^2 >= 1 (Jensen, applied to ln)

Run: OMP_NUM_THREADS=1 python _chk_iso_rho.py
"""
import numpy as np

M = {"H": 1.00782503, "D": 2.01410178, "T": 3.01604928}

# Errea 2013 Table I as transcribed in the deliverable (cm^-1 for wlog).
SSCHA_LAM = {"H": 0.40, "D": 0.46, "T": 0.48}
SSCHA_WLOG = {"H": 405.0, "D": 304.0, "T": 257.0}
HARM_LAM = {"H": 1.55, "D": 1.55, "T": 1.55}
HARM_WLOG = {"H": 205.0, "D": 150.0, "T": 125.0}

# Errea's lambda share of the H-character optical block (75-79%, isotope-dependent).
SSCHA_WOP = {"H": 0.750, "D": 0.783, "T": 0.792}


def rho_of_P(w, p):
    """rho for a discrete measure P = sum_i p_i delta(w_i). p must sum to 1."""
    p = np.asarray(p, float)
    w = np.asarray(w, float)
    assert abs(p.sum() - 1.0) < 1e-12
    w2 = float(p @ w**2)
    lnwl = float(p @ np.log(w))
    return w2 / np.exp(2.0 * lnwl)


# ---------------------------------------------------------------- (A) the lemma
print("=== (A) rho under UNIFORM dilation w -> c w (all modes) ===")
rng = np.random.default_rng(20260918)
worst = 0.0
for trial in range(6):
    n = rng.integers(2, 12)
    w = np.sort(rng.uniform(50.0, 900.0, n))
    p = rng.dirichlet(np.ones(n))
    for c in (0.70738, 0.5, 1.7, 0.31):
        d = abs(rho_of_P(c * w, p) / rho_of_P(w, p) - 1.0)
        worst = max(worst, d)
print(f"  6 random spectra x 4 dilation factors: max |rho(cw)/rho(w) - 1| = {worst:.3e}")
print("  (analytic: <w^2> -> c^2<w^2>, wlog -> c wlog, lam invariant ==> rho invariant)")

# ---------------------------------------------------------------- (B) hypothesis test
print()
print("=== (B) Is PdH's HARMONIC row consistent with uniform dilation? ===")
print("  uniform dilation would force wlog_B/wlog_A = sqrt(M_A/M_B) EXACTLY.")
for A, B in (("H", "D"), ("D", "T")):
    obs = HARM_WLOG[A] / HARM_WLOG[B]
    uni = np.sqrt(M[B] / M[A])
    # w_op implied if ONLY the H-sublattice block scales, P-weights fixed:
    #   ln wlog_A - ln wlog_B = w_op * ln(sqrt(M_B/M_A))
    wop = np.log(obs) / np.log(uni)
    print(f"  {A}->{B}: observed {obs:.4f}   uniform-dilation {uni:.4f}   "
          f"dev {100*(obs/uni-1):+.2f}%   implied w_op = {wop:.4f}")
print("  SSCHA row, same test (anharmonic, so deviation is expected to be larger):")
for A, B in (("H", "D"), ("D", "T")):
    obs = SSCHA_WLOG[A] / SSCHA_WLOG[B]
    uni = np.sqrt(M[B] / M[A])
    print(f"  {A}->{B}: observed {obs:.4f}   uniform-dilation {uni:.4f}   dev {100*(obs/uni-1):+.2f}%")

# ---------------------------------------------------------------- (C) decomposition
print()
print("=== (C) two-block s = rho_H/rho_D, and its harmonic/anharmonic split ===")


def solve_wop_freq(wlog, w_op, W_ac):
    """Omega_op from  ln wlog = w_ac ln W_ac + w_op ln W_op."""
    w_ac = 1.0 - w_op
    return float(np.exp((np.log(wlog) - w_ac * np.log(W_ac)) / w_op))


print(f"{'W_ac':>6} | {'Wop_H':>8} {'Wop_D':>8} {'ratio':>7} | {'rho_H':>7} {'rho_D':>7} "
      f"{'s_full':>7} | {'s_partsub':>9} {'s_anh':>7}")
for W_ac in (100.0, 130.0, 160.0, 190.0):
    Wop_H = solve_wop_freq(SSCHA_WLOG["H"], SSCHA_WOP["H"], W_ac)
    Wop_D = solve_wop_freq(SSCHA_WLOG["D"], SSCHA_WOP["D"], W_ac)
    rH = rho_of_P([W_ac, Wop_H], [1 - SSCHA_WOP["H"], SSCHA_WOP["H"]])
    rD = rho_of_P([W_ac, Wop_D], [1 - SSCHA_WOP["D"], SSCHA_WOP["D"]])
    s_full = rH / rD

    # Control: SAME H spectrum, but D obtained by the strictly HARMONIC partial
    # substitution -- optical block scaled by c, acoustic block fixed, weights frozen.
    c = np.sqrt(M["H"] / M["D"])
    rD_h = rho_of_P([W_ac, c * Wop_H], [1 - SSCHA_WOP["H"], SSCHA_WOP["H"]])
    s_partsub = rH / rD_h
    print(f"{W_ac:6.0f} | {Wop_H:8.1f} {Wop_D:8.1f} {Wop_H/Wop_D:7.4f} | "
          f"{rH:7.4f} {rD:7.4f} {s_full:7.4f} | {s_partsub:9.4f} {s_full/s_partsub:7.4f}")
print(f"  (uniform-dilation ratio for reference: sqrt(M_D/M_H) = {np.sqrt(M['D']/M['H']):.4f})")

# Same control on the harmonic row, where lam is mass-independent by construction.
print()
print("  HARMONIC row (lam=1.55 for all three, implied w_op from (B)):")
wop_h = np.log(HARM_WLOG["H"] / HARM_WLOG["D"]) / np.log(np.sqrt(M["D"] / M["H"]))
for W_ac in (100.0, 130.0, 160.0):
    Wop_H = solve_wop_freq(HARM_WLOG["H"], wop_h, W_ac)
    c = np.sqrt(M["H"] / M["D"])
    rH = rho_of_P([W_ac, Wop_H], [1 - wop_h, wop_h])
    rD = rho_of_P([W_ac, c * Wop_H], [1 - wop_h, wop_h])
    print(f"    W_ac={W_ac:5.0f}  W_op_H={Wop_H:7.1f}  R=W_op/W_ac={Wop_H/W_ac:5.2f}  "
          f"rho_H={rH:7.4f} rho_D={rD:7.4f}  s={rH/rD:7.4f}")

# ---------------------------------------------------------------- (D) sign of s
print()
print("=== (D) sign of s under partial substitution (frozen weights) ===")
print("  claim: d ln rho / d ln t = 2[ t^2 b2/(a2+t^2 b2) - w_B ] > 0  iff  t^2<w^2>_B > <w^2>_A")
bad = 0
for trial in range(4000):
    W_ac = rng.uniform(50.0, 250.0)
    W_op = rng.uniform(260.0, 900.0)
    w_op = rng.uniform(0.3, 0.95)
    t = np.sqrt(M["H"] / rng.choice([M["D"], M["T"]]))
    rH = rho_of_P([W_ac, W_op], [1 - w_op, w_op])
    rX = rho_of_P([W_ac, t * W_op], [1 - w_op, w_op])
    cond = (t * W_op) ** 2 > W_ac**2
    if cond and not (rH > rX):
        bad += 1
print(f"  4000 random two-block partial substitutions satisfying the condition: "
      f"violations of s>1 = {bad}")

# alpha shift induced by s, at the deliverable's numbers
lnMD = np.log(M["D"] / M["H"])
print()
print("=== (E) what s does to alpha ===")
print(f"  alpha shift = 0.5*ln(s)/ln(M_D/M_H) = {0.5/lnMD:.4f} * ln s")
for s in (1.0, 1.032, 1.10, 1.1667, 1.2536, 1.30):
    print(f"    s={s:6.4f}  ->  alpha_budget = {-0.3992 + 0.5*np.log(s)/lnMD:+.4f}")
