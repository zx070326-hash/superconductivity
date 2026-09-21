import os, sys, time
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "derivations"))
from p4_interchange_uniformity import H_mat
from p4d_coulomb_uniqueness import uc_vec, top_pencil
from p4d_coulomb_extremal_threshold import CellExtremal

N = 500
T_fixed = 1.0

# Wider search: vary the atom split and separation more aggressively,
# push mut higher, and use more extreme two-atom spectra (larger separation,
# more skewed weight) to hunt for a genuine Q(P) > Q(E) violation.
rng = np.random.default_rng(12345)

best_ratio = 0.0
best_info = None
n_tested = 0
n_violations = 0
n_below_thr_violations = 0

t0 = time.time()
for trial in range(400):
    p1 = rng.uniform(0.05, 0.95)
    weights = np.array([p1, 1 - p1])
    sep = rng.uniform(2.0, 40.0)
    base = rng.uniform(0.05, 2.0)
    omegas = np.array([base, base * sep])
    Nc = int(rng.choice([2, 3, 4, 5, 6, 8, 12, 16]))
    mut = float(rng.choice([0.10, 0.13, 0.20, 0.26, 0.33, 0.40, 0.50]))

    C = CellExtremal(weights, omegas, T_fixed, Nc, N)
    Q_P = C.Q_P(mut)
    Q_E = C.Q_E(mut)
    ratio = Q_P / Q_E if Q_E > 1e-15 else np.nan
    n_tested += 1

    if np.isfinite(ratio) and ratio > best_ratio:
        best_ratio = ratio
        best_info = (p1, sep, base, Nc, mut, C.thr, ratio)

    if ratio > 1 + 1e-5:
        n_violations += 1
        below = mut < C.thr
        if below:
            n_below_thr_violations += 1
            print(f"FALSE NEGATIVE: p1={p1:.3f} sep={sep:.2f} base={base:.3f} "
                  f"Nc={Nc} mut={mut} thr={C.thr:.4f} ratio={ratio:.5f}")

print(f"\nscan done in {time.time()-t0:.1f}s, {n_tested} trials")
print(f"violations found: {n_violations}  (false negatives below threshold: {n_below_thr_violations})")
print(f"\nbest (highest) Q(P)/Q(E) ratio found: {best_ratio:.6f}")
if best_info:
    p1, sep, base, Nc, mut, thr, ratio = best_info
    print(f"  at p1={p1:.3f}, sep={sep:.2f}, base={base:.3f}, Nc={Nc}, mut={mut}, thr={thr:.4f}, ratio={ratio:.6f}")
