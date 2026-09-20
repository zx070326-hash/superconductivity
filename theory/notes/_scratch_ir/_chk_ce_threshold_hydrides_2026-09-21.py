"""INDEPENDENT recomputation of the CE' extremality threshold for the hydride rows.

This does NOT cite the repo's 14/14 anchor table; every constant below is rebuilt
from the single primitive H_mat.  Theorem CE' (P4D_coulomb_extremal_threshold.md
sec.6), restated:

  A_P := O[F_P] = sum_i w_i H_mat(N, wbar_i)      (O is linear in the kernel)
  (mu_j, psi_j) := eigenpairs of A_P, mu_1 > mu_2 >= ...
  u_c[n] := (2n+1)^{-1/2} 1_{n<Nc},   K_c := 2 u_c u_c^T,   c_j := <u_c, psi_j>
  lambda := lam_1(A_P - mutilde K_c)                     (the secular root)
  d(k) := [[k]](wbar_rms) - F_P(k) >= 0                  (kernel-level Jensen)
  O[d] = H_mat(N, wbar_rms) - sum_i w_i H_mat(N, wbar_i)
  q1 := <psi_1, O[d] psi_1>,  b := ||P_perp O[d] psi_1||,  m := max(0,-lam_min(O[d]))
  t^d := (-b + sqrt(b^2 + m q1))/m          (= inf if m = 0)
  tanTheta_P(mutilde) := (mu_1 - lambda) * ( sum_{j>=2} c_j^2/(lambda-mu_j)^2 )^{1/2} / c_1

  THEOREM CE':  tanTheta_P(mutilde) < t^d   ==>   Q(P) <= Q(E),
  i.e. the Einstein spectrum at wbar_rms maximizes lam_1 AT THE SAME mutilde and Nc.

  mutilde := mu* / lam       <-- NOT mu*.  (repo memory: this exact slip cost 8x.)

PROBE FAMILY (declared, because the threshold is a functional of P):
  symmetric two-delta,  weights (1/2,1/2),  frequencies wbar_rms * s / kappa and
  wbar_rms * s * kappa  with  s = sqrt(2/(kappa^2 + kappa^-2))  so that the rms is
  EXACTLY wbar_rms for every kappa.  kappa = 1 is the Einstein point (d == 0,
  threshold degenerate); kappa > 1 is a genuine spread.

DIRECTION AUDIT: CE' is a SUFFICIENT condition.  tanTheta < t^d => PASS.
tanTheta >= t^d proves NOTHING (the extremality may still hold); it must be
reported as UNDECIDED, never as a failure of extremality.

LEVEL: numerical-evidence (float64, N truncation, Perron-Frobenius used as in the
repo's (H1)-(H5)).
"""
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DERIV = os.path.normpath(os.path.join(HERE, "..", "..", "derivations"))
sys.path.insert(0, DERIV)
sys.path.insert(0, HERE)
from p4_interchange_uniformity import H_mat  # noqa: E402
import importlib.util as _ilu  # noqa: E402


def _load(stem):
    spec = _ilu.spec_from_file_location(
        stem.replace("-", "_"), os.path.join(HERE, stem + ".py"))
    mod = _ilu.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


_E = _load("_phiE_engine_2026-09-21")
N_DEF = 400


def two_delta(wrms, kappa):
    s = np.sqrt(2.0 / (kappa**2 + kappa**-2))
    return np.array([0.5, 0.5]), np.array([wrms * s / kappa, wrms * s * kappa])


def ce_prime(wrms, kappa, Nc, mu_t, N=N_DEF):
    """Return dict with tanTheta_P, t^d and the six constants, all rebuilt here."""
    wts, wbs = two_delta(wrms, kappa)
    assert abs(np.sqrt(np.dot(wts, wbs**2)) / wrms - 1.0) < 1e-13
    A = sum(w * H_mat(N, float(x)) for w, x in zip(wts, wbs))
    Od = H_mat(N, float(wrms)) - A
    mu, psi = np.linalg.eigh(A)
    mu = mu[::-1]
    psi = psi[:, ::-1]
    p1 = psi[:, 0]
    if p1.sum() < 0:
        p1 = -p1
        psi[:, 0] = p1
    u = np.zeros(N)
    k = max(0, min(int(Nc), N))
    u[:k] = 1.0 / np.sqrt(2.0 * np.arange(k) + 1.0)
    c = psi.T @ u
    # six constants on A_P and the kernel d
    q1 = float(p1 @ (Od @ p1))
    y = Od @ p1
    b = float(np.linalg.norm(y - (p1 @ y) * p1))
    m = max(0.0, -float(np.linalg.eigvalsh(Od)[0]))
    td = np.inf if m == 0.0 else (-b + np.sqrt(b * b + m * q1)) / m
    # secular root lambda = lam_1(A - mu_t K_c)
    M = A - mu_t * 2.0 * np.outer(u, u)
    lam = float(np.linalg.eigvalsh(M)[-1])
    den = (lam - mu[1:]) ** 2
    tanT = (mu[0] - lam) * np.sqrt(float(np.sum(c[1:] ** 2 / den))) / c[0]
    return dict(tanTheta=float(tanT), td=float(td), q1=q1, b=b, m=m,
                gamma=float(mu[0] - mu[1]), c1=float(c[0]), lam=lam,
                ucnorm=float(np.linalg.norm(u)))


def mutilde_star(wrms, kappa, Nc, N=N_DEF, hi=8.0):
    """solve tanTheta_P(mutilde) = t^d (tanTheta strictly increasing in mutilde)."""
    r0 = ce_prime(wrms, kappa, Nc, 0.0, N)
    if not np.isfinite(r0['td']):
        return np.inf, r0
    if r0['tanTheta'] >= r0['td']:
        return 0.0, r0
    lo = 0.0
    if ce_prime(wrms, kappa, Nc, hi, N)['tanTheta'] < r0['td']:
        return hi, r0
    for _ in range(46):
        mid = 0.5 * (lo + hi)
        if ce_prime(wrms, kappa, Nc, mid, N)['tanTheta'] < r0['td']:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), r0


# ------------------------------------------------------------------- materials
# lam, wbar2 (meV) or None, mu*, label.  wbar2 from theory/data/hydrides_theory.csv
MATS = [
    ("H3S 200 harm",  2.2929, 129.32, 0.13),
    ("H3S 200 anh",   1.7210, 149.26, 0.13),
    ("LaH10 129",     3.6200, None,   0.10),
    ("LaH10 163",     2.6700, None,   0.10),
    ("LaH10 214",     2.0600, None,   0.10),
    ("LaH10 264",     1.7300, None,   0.10),
    ("MgB2 (scalar)", 0.7400, None,   0.10),
]


def main():
    t0 = time.time()
    print("=" * 86)
    print("CE' THRESHOLD, independently rebuilt   [numerical-evidence]")
    print("=" * 86)
    print("mutilde := mu*/lam.  Probe P: symmetric two-delta, rms pinned to wbar_rms.")
    print()
    for r in (3.0, 10.0):
        print(f"--- r = omega_c/omega_E = {r:g} ---")
        print(f"{'material':<15}{'lam':>6}{'mu*':>6}{'mutilde':>9}{'wbar*':>8}"
              f"{'Nc':>4}{'kappa':>7}{'tanTh':>9}{'t^d':>9}{'mu~_*':>9}"
              f"{'margin':>8}  verdict")
        print("-" * 104)
        for name, lam, wb2, mus in MATS:
            mu_t = mus / lam
            w, nc, _ = _E.wstar_sup(lam, mus, r, N=300)
            if w is None:
                print(f"{name:<15}{lam:6.3f}{mus:6.2f}  no root")
                continue
            for kappa in (1.05, 1.2, 1.5, 2.0):
                ms, r0 = mutilde_star(w, kappa, nc)
                rr = ce_prime(w, kappa, nc, mu_t)
                ok = rr['tanTheta'] < rr['td']
                print(f"{name:<15}{lam:6.3f}{mus:6.2f}{mu_t:9.4f}{w:8.4f}"
                      f"{nc:4d}{kappa:7.2f}{rr['tanTheta']:9.4f}{rr['td']:9.4f}"
                      f"{ms:9.4f}{ms/mu_t:8.3f}x  "
                      f"{'PASS (extremality certified)' if ok else 'UNDECIDED'}")
        print()

    print("N-convergence of the decision quantities (H3S 200 harm, r=10, kappa=1.2):")
    lam, mus, r = 2.2929, 0.13, 10.0
    w, nc, _ = _E.wstar_sup(lam, mus, r, N=300)
    for N in (200, 300, 400, 600):
        rr = ce_prime(w, 1.2, nc, mus / lam, N=N)
        print(f"  N={N:4d}  tanTheta={rr['tanTheta']:.8f}  t^d={rr['td']:.8f}  "
              f"q1={rr['q1']:.8f}  b={rr['b']:.8f}  m={rr['m']:.8f}")
    print()
    print("TESTED RANGES (this is the premise's honest scope; the threshold clause of")
    print("the conditional theorem is numerical-evidence, not certificate):")
    print("  spectral family : symmetric two-delta, kappa in {1.05,1.2,1.5,2.0}")
    print("                    (rms pinned exactly; kappa=1 is the degenerate point)")
    print("  N_c             : the REALIZED window at each material's own root")
    print("  mutilde         : mu*/lam in [0.028, 0.135] for the rows above")
    print("  truncation      : N = 400 (convergence shown above)")
    print()
    print("DIRECTION: 'UNDECIDED' rows are NOT counterexamples to extremality; CE' is")
    print("only sufficient.  A ceiling that relies on extremality may be quoted ONLY")
    print("on PASS rows.")
    print(f"wall {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
