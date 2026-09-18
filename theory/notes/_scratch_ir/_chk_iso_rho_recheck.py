"""INDEPENDENT recomputation of the s-interval table, second path.

Path 1 (previous, _chk_iso_rho_budget.py) used a CLOSED FORMULA for beta_lo.
Path 2 (here) does NOT use that formula.  It searches directly over admissible
multi-mode spectra and minimises beta numerically, so any algebra error in path 1
shows up as a disagreement.

Admissible spectrum, given (w_B, Om_ac_max, wlog):
   group A (unsubstituted, Pd-character): any modes with  w <= Om_ac_max, total P-weight 1-w_B
   group B (substituted,  H-character):   any modes,                      total P-weight w_B
   constraint:  sum p_i ln w_i = ln wlog      (this is the measured wlog)
and  s = t^(2 w_B) / (1 - beta (1-t^2)),  beta = b2/(a2+b2)   [identity, already verified]

Run: OMP_NUM_THREADS=1 python _chk_iso_rho_recheck.py
"""
import numpy as np
from scipy.optimize import minimize

M = {"H": 1.00782503, "D": 2.01410178}
T_HD = np.sqrt(M["H"] / M["D"])
WLOG = 405.0
SHIFT = 0.5 / np.log(M["D"] / M["H"])
A_BASE = -0.3992


def s_of_beta(w_B, beta, t=T_HD):
    return t ** (2.0 * w_B) / (1.0 - beta * (1.0 - t**2))


def beta_of_spectrum(wA, pA, wB, pB):
    a2 = float(np.asarray(pA) @ np.asarray(wA) ** 2)
    b2 = float(np.asarray(pB) @ np.asarray(wB) ** 2)
    return b2 / (a2 + b2)


def min_beta_numeric(w_B, om_max, nA=3, nB=3, seed=0):
    """Minimise beta over admissible spectra by direct search. No closed form used.

    The wlog constraint is enforced EXACTLY by construction: after choosing the
    A-group and the SHAPE of the B-group, the whole B-group is rescaled by c with
        ln c = [ln wlog - pA.ln fA - pB.ln fB] / w_B
    which hits ln wlog exactly, so no penalty term and no residual slack.
    """
    rng = np.random.default_rng(seed)
    w_A = 1.0 - w_B
    best = (np.inf, None)

    def softmax(z, tot):
        z = z - z.max()
        e = np.exp(z)
        return tot * e / e.sum()

    def unpack(z):
        fA = om_max / (1.0 + np.exp(-np.clip(z[:nA], -40, 40)))
        pA = softmax(z[nA:2 * nA], w_A)
        fB = np.exp(np.clip(z[2 * nA:2 * nA + nB], -20, 20))
        pB = softmax(z[2 * nA + nB:], w_B)
        lnc = (np.log(WLOG) - float(pA @ np.log(fA)) - float(pB @ np.log(fB))) / w_B
        return fA, pA, np.exp(lnc) * fB, pB

    def obj(z):
        fA, pA, fB, pB = unpack(z)
        return beta_of_spectrum(fA, pA, fB, pB)

    for _ in range(40):
        z0 = np.r_[rng.normal(2.0, 1.5, nA), rng.normal(0, 1.0, nA),
                   rng.normal(0.0, 0.5, nB), rng.normal(0, 1.0, nB)]
        res = minimize(obj, z0, method="Nelder-Mead",
                       options={"maxiter": 40000, "maxfev": 40000,
                                "xatol": 1e-12, "fatol": 1e-14})
        fA, pA, fB, pB = unpack(res.x)
        lw = float(pA @ np.log(fA) + pB @ np.log(fB))
        assert abs(lw - np.log(WLOG)) < 1e-9, "wlog constraint broken by construction"
        b = beta_of_spectrum(fA, pA, fB, pB)
        if b < best[0]:
            best = (b, (fA, pA, fB, pB))
    return best


print("=== A. beta_lo by DIRECT SEARCH over admissible spectra (no closed form) ===")
print(f"  Om_ac_max = 190 cm^-1, wlog = {WLOG} cm^-1")
print(f"{'w_B':>8} | {'beta_lo(search)':>15} | minimiser structure (top modes, weight)")
found = {}
for w_B in (0.750, 0.790, 0.9023):
    b, sp = min_beta_numeric(w_B, 190.0, seed=11)
    found[w_B] = b
    fA, pA, fB, pB = sp
    iA, iB = np.argsort(-pA), np.argsort(-pB)
    desc = ("A:" + ",".join(f"{fA[i]:.0f}({pA[i]:.2f})" for i in iA[:2])
            + "  B:" + ",".join(f"{fB[i]:.0f}({pB[i]:.2f})" for i in iB[:2]))
    print(f"{w_B:8.4f} | {b:15.5f} | {desc}")

print()
print("=== B. Table recomputed from the SEARCHED beta_lo ===")
print(f"{'w_B':>8} | {'beta_lo':>8} | {'s_lo':>7} {'s_hi':>7} {'width':>7} | "
      f"{'alpha_lo':>9} {'alpha_hi':>9} {'a_width':>8}")
for w_B in (0.750, 0.790, 0.9023):
    b = found[w_B]
    s_lo, s_hi = s_of_beta(w_B, b), s_of_beta(w_B, 1.0)
    a_lo, a_hi = A_BASE + SHIFT * np.log(s_lo), A_BASE + SHIFT * np.log(s_hi)
    print(f"{w_B:8.4f} | {b:8.5f} | {s_lo:7.4f} {s_hi:7.4f} {s_hi-s_lo:7.4f} | "
          f"{a_lo:+9.4f} {a_hi:+9.4f} {a_hi-a_lo:8.4f}")

print()
print("=== C. Disputed point: should beta's range be FROZEN across w_B? ===")
print("  team-lead scanned beta in [0.958, 1.000] for EVERY row.")
print("  0.958 is beta_lo evaluated at w_B=0.750.  Admissible at w_B=0.9023?")
b_min = found[0.9023]
print(f"  beta_lo(w_B=0.9023) by direct search = {b_min:.5f}")
print(f"  beta=0.958 at w_B=0.9023 -> "
      f"{'reachable' if 0.958 >= b_min else 'BELOW the attainable minimum (not a spectrum)'}")
print("  what Om_ac_max would be needed to push beta_lo down to 0.958 at w_B=0.9023:")
for om in (190.0, 260.0, 340.0, 420.0, 500.0):
    w_A = 1 - 0.9023
    G = (WLOG / om**w_A) ** (1.0 / 0.9023)
    bb = 1.0 / (1.0 + (w_A * om**2) / (0.9023 * G**2))
    print(f"    Om_ac_max={om:6.0f} cm^-1 -> beta_lo={bb:.5f}")

print()
print("=== D. Direction audit of the disagreement ===")
print("  verdict the table supports: 'end-to-end width >> 0.05' == target NOT reachable")
lo_all = min(s_of_beta(w, found[w]) for w in found)
hi_all = max(s_of_beta(w, 1.0) for w in found)
tl_lo, tl_hi = s_of_beta(0.9023, 0.958), s_of_beta(0.750, 1.0)
print(f"  mine (per-row self-consistent beta_lo): [{lo_all:.4f}, {hi_all:.4f}] width {hi_all-lo_all:.4f}")
print(f"  team-lead (beta>=0.958 everywhere)    : [{tl_lo:.4f}, {tl_hi:.4f}] width {tl_hi-tl_lo:.4f}")
print("  both >> 0.05 ==> STOP verdict unchanged either way.")
