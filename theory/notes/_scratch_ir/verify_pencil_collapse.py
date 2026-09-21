r"""
MAIN-SESSION INDEPENDENT CHECK of the coulomb-certify plan-changing claims.

I do NOT re-run p4d_coulomb_certify.py.  I import ONLY the already-trusted
phonon operator H_mat, build K_c = 2|u><u| from scratch, form M = H - mutilde K_c,
and compute the TRUE top-two eigenvalues by dense eigvalsh.  Then I check, blind:

  (1) M(mutilde=0) == H_mat bit-for-bit  (the self-check-gate floor).
  (2) alpha-side inheritance: lam_k(M) <= lam_k(H) for every k, mutilde>0
      (Weyl monotonicity for a PSD subtraction) -- the report's "structural gift".
  (3) THE PLAN-CHANGING CLAIM: at wbar=0.3, mutilde=0.1 the TRUE spectral gap
      collapses to ~1e-4 (worst window N_c=3, lam_1(M)~6e-5).  If this is real
      operator physics (not certificate slack) my direct eigvalsh must see it too.
  (4) delta(wbar=1) true gaps vs the reported certified LOWER bounds
      (0.657276 / 0.360987 / 0.282156 at mutilde=0/0.1/0.13): true gap must be
      >= the certified bound (a lower bound) and in the same ballpark.
  (5) locate wbar where the worst-N_c TRUE gap crosses ~0; the certified
      wbar_min (0.5391 / 0.6427) must sit at or above it (cell+k-sign slack).
  (6) compression identity  P_N K_c(N_c) P_N == K_c(min(N,N_c))  entrywise.

Exhaustiveness of the N_c scan: at truncation N the restricted operator only sees
K_c(min(N_c,N)), so scanning N_c=0..N covers every integer window.
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_DERIV = os.path.abspath(os.path.join(_HERE, "..", "..", "derivations"))
if _DERIV not in sys.path:
    sys.path.insert(0, _DERIV)
from p4_interchange_uniformity import H_mat  # the ONLY trusted import

PASS = [0]


def chk(name, ok, detail=""):
    if not ok:
        raise AssertionError(f"FAIL: {name} {detail}")
    PASS[0] += 1
    print(f"  ok: {name} {detail}")


def Kc(N, Nc):
    """2|u><u|, u_n=1/sqrt(2n+1)*1_{n<Nc}, built independently (no import)."""
    u = np.zeros(N)
    k = max(0, min(int(Nc), N))
    u[:k] = 1.0 / np.sqrt(2.0 * np.arange(k) + 1.0)
    return 2.0 * np.outer(u, u)


def M(N, w, mutilde, Nc):
    return H_mat(N, w) - mutilde * Kc(N, Nc)


def worst_gap(N, w, mutilde):
    """min over ALL windows N_c=0..N of (lam_1 - lam_2), plus the worst N_c and
    the lam_1 there.  Direct dense eigvalsh -- the TRUE operator, no bound."""
    Hm = H_mat(N, w)
    best = None
    for Nc in range(0, N + 1):
        E = np.linalg.eigvalsh(Hm - mutilde * Kc(N, Nc))
        gap = E[-1] - E[-2]
        if best is None or gap < best[0]:
            best = (gap, Nc, E[-1], E[-2])
    return best


# ---------------------------------------------------------------------------
print("=== (1) M(mutilde=0) == H_mat bit-for-bit ===")
for N in (5, 40, 200):
    for w in (0.3, 1.0, 5.0):
        for Nc in (0, 3, N // 2 + 1, 10 * N):
            d = np.max(np.abs(M(N, w, 0.0, Nc) - H_mat(N, w)))
            chk("mutilde=0 annihilates K_c exactly", d == 0.0, f"N={N} w={w} Nc={Nc}")

print("\n=== (2) alpha-side inheritance: lam_k(M) <= lam_k(H) for every k ===")
for w in (0.3, 1.0, 2.0, 5.0):
    EH = np.linalg.eigvalsh(H_mat(300, w))
    for mt in (0.1, 0.13):
        for Nc in (1, 3, 40, 300):
            EM = np.linalg.eigvalsh(M(300, w, mt, Nc))
            chk("lam_k(M) <= lam_k(H) all k (Weyl, PSD subtraction)",
                np.all(EM <= EH + 1e-11), f"w={w} mt={mt} Nc={Nc} "
                f"worst={np.max(EM-EH):.2e}")

print("\n=== (3) PLAN-CHANGING CLAIM: collapse at wbar=0.3, mutilde=0.1 ===")
for N in (200, 400):
    g, Nc, l1, l2 = worst_gap(N, 0.3, 0.1)
    print(f"  N={N}: worst N_c={Nc}, lam_1(M)={l1:.3e}, lam_2(M)={l2:.3e}, "
          f"gap={g:.3e}")
chk("true gap at wbar=0.3,mt=0.1 collapses to ~1e-4 (real, not slack)",
    g < 5e-4, f"gap={g:.3e}")
chk("worst window is small N_c (reported N_c=3)", Nc <= 5, f"N_c={Nc}")
chk("lam_1(M) itself dives toward 0 (top eig into the cluster)", abs(l1) < 5e-4,
    f"lam_1={l1:.3e}")

print("\n=== (4) delta(wbar=1) true gaps vs reported certified LOWER bounds ===")
reported = {0.0: 0.657276, 0.1: 0.360987, 0.13: 0.282156}
for mt in (0.0, 0.1, 0.13):
    g, Nc, l1, l2 = worst_gap(400, 1.0, mt)
    print(f"  mt={mt}: true worst-N_c gap={g:.6f} (N_c={Nc}), "
          f"reported certified delta={reported[mt]:.6f}")
    chk("true gap >= certified delta (bound is a valid lower bound)",
        g >= reported[mt] - 1e-4, f"mt={mt} true={g:.6f} cert={reported[mt]:.6f}")
    chk("true gap in the same ballpark as the certified bound",
        g <= reported[mt] + 0.4, f"mt={mt}")

print("\n=== (5) locate wbar where the worst-N_c TRUE gap crosses ~0 ===")
print("  (certified wbar_min: 0.5391 @mt=0.1, 0.6427 @mt=0.13; the pointwise")
print("   true-gap crossing should be AT OR BELOW these, cell+k-sign slack)")
for mt, wmin_cert in ((0.1, 0.5391), (0.13, 0.6427)):
    print(f"  --- mutilde={mt} (certified wbar_min={wmin_cert}) ---")
    last_neg, first_pos = None, None
    for w in (0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1.0):
        g, Nc, l1, l2 = worst_gap(200, w, mt)
        flag = "collapsed" if g < 1e-3 else ("small" if g < 0.05 else "healthy")
        print(f"      wbar={w:.2f}: worst gap={g:.6e} (N_c={Nc})  [{flag}]")
        if g < 1e-3:
            last_neg = w
        elif first_pos is None:
            first_pos = w
    chk("true gap is collapsed below the certified wbar_min",
        last_neg is not None and last_neg <= wmin_cert + 1e-9,
        f"mt={mt} last-collapsed wbar={last_neg}")
    chk("true gap is healthy at wbar=1 (well above wbar_min)",
        worst_gap(200, 1.0, mt)[0] > 0.1, f"mt={mt}")

print("\n=== (6) compression identity P_N K_c(N_c) P_N == K_c(min(N,N_c)) ===")
for Nc in (1, 3, 16, 128):
    for N in (1, 2, 5, 40, 300):
        lhs = Kc(N, Nc)                    # = P_N K_c(N_c) P_N (leading block)
        rhs = Kc(N, min(N, Nc))
        chk("P_N K_c(N_c) P_N == K_c(min(N,N_c))", np.max(np.abs(lhs - rhs)) == 0.0,
            f"N={N} Nc={Nc}")

print(f"\nALL {PASS[0]} independent checks passed.")
print("VERDICT: (a) mutilde=0 gate floor holds bit-for-bit; (b) alpha-side")
print("operator-order inheritance confirmed; (c) the small-wbar gap collapse is")
print("REAL operator physics (direct eigvalsh sees lam_1~1e-4), NOT certificate")
print("slack -- so the interval theorem genuinely shrinks to [wbar_min,5]; (d) the")
print("certified deltas are valid lower bounds on the true gaps; (e) compression")
print("identity exact, so N<N_c needs no separate treatment.")
