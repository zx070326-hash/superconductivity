"""[numerical-evidence] Scoping P_beta: where does it actually BIND, and at which beta?

lam2bound is attacking "prove rowsum_beta(n) - Delta_n is eventually decreasing,
for all beta, for all varpi".  Before spending on that, check how much of that
generality the 12-cell criterion actually consumes:

  (i) On which cells is min(P_beta, Lam_full, 0.6h) attained by P_beta?
      On the rest, S~/Lam_full governs and its finitization is ALREADY closed
      (Lemma F), so P_beta's sup never needs to be evaluated there.
 (ii) Lemma Lambda holds for EACH beta separately.  So it is sound to FIX one
      rational beta per cell (the one achieving the min) -- the "for all beta"
      quantifier is a convenience of the driver, not a requirement of the proof.
(iii) Over the varpi range where P_beta binds, where is the argmax n, and how
      much headroom is there between the peak and the tail?  Headroom is what a
      CRUDE analytic tail bound needs; it does not need to be tight.

Reuses _chk_grid_coarsen.peel_all_betas / BETAS unchanged, and _chk_padR_full for
Lam_full.  Run: OMP_NUM_THREADS=1 python _chk_Pbeta_scope.py
"""
import os, io, contextlib, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, HERE)


def _load(name, alias):
    sp = importlib.util.spec_from_file_location(alias, os.path.join(HERE, name))
    m = importlib.util.module_from_spec(sp)
    with contextlib.redirect_stdout(io.StringIO()):
        sp.loader.exec_module(m)
    return m


gc = _load("_chk_grid_coarsen.py", "gc")
pf = _load("_chk_padR_full.py", "pf")
C_R = 1.60
CELLS = [(0.30, 0.42), (0.42, 0.60), (0.60, 0.85), (0.85, 1.20), (1.20, 1.70),
         (1.70, 2.40), (2.40, 3.40), (3.40, 4.80), (4.80, 6.80), (6.80, 9.60),
         (9.60, 13.60), (13.60, 20.00)]


def row_profile(w, beta, NROW=400, MT=20000):
    """V_beta(n) for n=1..NROW, independently re-implemented (not gc's blocked loop).

    V_beta(n) = (2n+1)^{(beta-1)/2} * sum_{m>=1} K(n,m) (2m+1)^{-(1+beta)/2}
                - l_G(n),
    K(n,m) = f(|n-m|) + f(n+m+1), except K(n,n) = f(2n+1)  (no Toeplitz f(0)).
    No tail pad here: this probe is about WHERE the max sits, not its value.
    """
    x = w * w
    f = lambda t: x / (t * t + x)
    m = np.arange(1, MT, dtype=float)
    wm = (2.0 * m + 1.0) ** (-(1.0 + beta) / 2.0)
    out = np.empty(NROW)
    fk = f(np.arange(0, NROW + 2, dtype=float)); fk[0] = 0.0
    cum = np.cumsum(fk)
    for i in range(1, NROW + 1):
        K = f(np.abs(i - m)) + f(i + m + 1.0)
        K[i - 1] = f(2.0 * i + 1.0)
        out[i - 1] = ((2.0 * i + 1.0) ** ((beta - 1.0) / 2.0) * (K @ wm)
                      - 2.0 * cum[i] / (2.0 * i + 1.0))
    return out


if __name__ == "__main__":
    print("== (i)+(ii) which bound binds per cell, and at which beta ==")
    print(f"{'cell':>15} {'P_beta':>9} {'beta*':>6} {'Lam_full':>9} {'binds':>9}")
    binds = []
    for lo, hi in CELLS:
        pb_w, lf_w, bstar = -np.inf, -np.inf, None
        for w in np.geomspace(lo, hi, 11):
            pb = gc.peel_all_betas(w)
            j = int(np.argmin(pb))
            if pb[j] > pb_w:
                pb_w, bstar = pb[j], gc.BETAS[j]
            lf_w = max(lf_w, pf.lam_full(w, C_R))
        who = "P_beta" if pb_w < lf_w else "Lam_full"
        binds.append((lo, hi, who, bstar))
        print(f"[{lo:5.2f},{hi:5.2f}] {pb_w:9.4f} {bstar:6.2f} {lf_w:9.4f} {who:>9}")

    pcells = [c for c in binds if c[2] == "P_beta"]
    print(f"\n   P_beta binds on {len(pcells)}/{len(CELLS)} cells: "
          + ", ".join(f"[{a:.2f},{b:.2f}](b={bb:g})" for a, b, _, bb in pcells))
    print("   => outside those, Lemma F already closes the finitization.")

    print("\n== (iii) peak location and tail headroom where P_beta binds ==")
    print(f"{'varpi':>6} {'beta':>5} {'argmax n':>9} {'peak':>10} {'V(50)':>10}"
          f" {'V(100)':>10} {'V(400)':>10} {'peak/V(50)':>11} {'mono from':>10}")
    for lo, hi, _, bstar in pcells:
        for w in (lo, np.sqrt(lo * hi), hi):
            v = row_profile(w, float(bstar))
            k = int(np.argmax(v)) + 1
            d = np.diff(v[k - 1:])
            n0 = k + (int(np.argmax(d[::-1] > 0)) if (d > 0).any() else 0)
            n0 = k if not (d > 0).any() else k + int(np.max(np.where(d > 0)[0])) + 1
            print(f"{w:6.3f} {bstar:5.2f} {k:9d} {v[k-1]:10.5f} {v[49]:10.3e}"
                  f" {v[99]:10.3e} {v[399]:10.3e} {v[k-1]/v[49]:11.1f} {n0:10d}")
