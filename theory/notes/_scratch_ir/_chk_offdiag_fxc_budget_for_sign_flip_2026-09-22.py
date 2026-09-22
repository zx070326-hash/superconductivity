"""
Budget-only task (team-lead, 2026-09-22): given the REAL Im chi^s used by
essenberger's pipeline, how strong must an ORBITAL OFF-DIAGONAL f^xc_{ll'}
be (relative to the diagonal, which is currently implicitly 1) to move
|D_ij|/|Lambda_ij| from its baseline median 0.998 down to a level where the
overall sign of Lambda_ij could plausibly flip?

Model: f^xc_{ll} = 1 (diagonal, unchanged), f^xc_{ll'} = kappa for l != l'
(uniform off-diagonal strength, single scalar knob -- this is the minimal
model needed to answer "how strong", not a full orbital-resolved fit).

Then  Lambda_ij(kappa) = D_ij + kappa * R_ij
  where D_ij  = sum_l            w_l(i) w_l(j) Im chi^s_ll(q_ij)      (l==l', kappa-independent)
        R_ij  = sum_{l!=l'}      w_l(i) w_l'(j) Im chi^s_ll'(q_ij)    (l!=l', baseline weight 1,
                                                                        i.e. R_ij = Lambda_ij(1) - D_ij)

Sign of Lambda_ij(kappa) can flip only once kappa*R_ij is comparable in
magnitude to -D_ij (opposite sign), i.e. at kappa* = -D_ij/R_ij (if that is
positive; if D_ij and R_ij already have the same sign, NO finite positive
kappa flips it -- increasing kappa only reinforces the sign, and that pair
should be excluded / reported separately).

Reuses the exact same pipeline (H, mu, chi0 wrapper, FS point set, seed) as
_chk_lambdaSF_positivity_is_perron_not_material_2026-09-22.py so the D_ij,
R_ij, |D|/|Lambda| numbers reproduce the 0.998 / 0.959 medians before any
new claim is trusted.
"""
import numpy as np, importlib.util, os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
def load(f):
    s = importlib.util.spec_from_file_location(f, os.path.join(SCRIPT_DIR, f))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m

N = 5
gate = load("_chk_gauge_invariance_gate_dynamic_2026-09-22.py")
ham  = load("_graser1003_hamiltonian_2026-09-22.py")
H_func = ham.build_H_graser
mu = float(np.load(os.path.join(SCRIPT_DIR, "_graser_mu_value.npy")))

def chi0(q, w):
    return -gate.chi0_dynamic(q, w, H_func, mu, T=0.02, eta=0.02, k_grid=16, norb=N, rephase=False)
def dress(c, U=0.3):
    return c @ np.linalg.inv(np.eye(N) - (U*np.eye(N)) @ c)

ks = np.linspace(-np.pi, np.pi, 36, endpoint=False)
pts = []
for kx in ks:
    for ky in ks:
        e, U = np.linalg.eigh(H_func(kx, ky))
        for n in range(N):
            if abs(e[n]-mu) < 0.03:
                pts.append((kx, ky, np.abs(U[:, n])**2))
sub = np.random.default_rng(20260922).choice(len(pts), size=25, replace=False)
pts = [pts[i] for i in sub]

for omega in (0.05, 0.2):
    cache = {}
    def chiq(qx, qy):
        k = (round(qx,6), round(qy,6))
        if k not in cache: cache[k] = dress(chi0((qx,qy), omega))
        return cache[k]

    Ds, Rs, ratios0 = [], [], []
    kappa_star = []      # kappa needed to flip, only where D,R have opposite sign
    same_sign_count = 0
    total = 0
    for i in range(len(pts)):
        for j in range(i+1, len(pts)):
            M = chiq(pts[i][0]-pts[j][0], pts[i][1]-pts[j][1]).imag
            wi, wj = pts[i][2], pts[j][2]
            D = float(np.sum(wi*np.diag(M)*wj))
            full1 = float(wi @ M @ wj)     # Lambda_ij at kappa=1 (baseline, current pipeline)
            R = full1 - D
            total += 1
            Ds.append(D); Rs.append(R); ratios0.append(abs(D)/max(abs(full1),1e-30))
            if D == 0 or R == 0:
                continue
            if np.sign(D) == np.sign(R):
                same_sign_count += 1
                continue
            kappa_star.append(abs(D / R))

    ratios0 = np.array(ratios0)
    kappa_star = np.array(kappa_star)
    print(f"=== omega={omega} ===")
    print(f"  pairs: {total}   |D|/|Lambda(kappa=1)| median = {np.median(ratios0):.4f}  (paper quotes 0.998/0.959)")
    print(f"  pairs where D,R SAME sign (no finite kappa flips them): {same_sign_count}/{total} = {100*same_sign_count/total:.1f}%")
    print(f"  pairs where D,R OPPOSITE sign (flippable): {len(kappa_star)}/{total}")
    if len(kappa_star):
        print(f"  kappa* = |D/R| needed to flip:  median={np.median(kappa_star):.3g}"
              f"  min={np.min(kappa_star):.3g}  10th pct={np.percentile(kappa_star,10):.3g}"
              f"  90th pct={np.percentile(kappa_star,90):.3g}  max={np.max(kappa_star):.3g}")
    print()
