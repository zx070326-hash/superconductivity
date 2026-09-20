import os, sys, time
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "derivations"))
from p4_interchange_uniformity import H_mat
from p4d_coulomb_uniqueness import uc_vec, top_pencil
from p4d_coulomb_extremal_threshold import (
    random_spectrum, spectrum_to_K, spectrum_wrms, O_d_for_spectrum,
    CellExtremal
)

N = 600

def H_multi(NN, s, comps, order=0):
    out = np.zeros((NN, NN))
    for p, c in comps:
        out += p * H_mat(NN, c * s, order)
    return out

# Test extremality Q(P)<=Q(E) on the same 2-atom danger spectrum
comps = [(0.8, 0.3), (0.2, 4.0)]
ss = np.arange(0.05, 3.0 + 1e-9, 0.02)
Nc_list = [5, 8, 12]
mut_list = [0.13, 0.26, 0.33]
T_fixed = 1.0

n_tested = 0
n_violations = 0
n_below_threshold = 0
n_correctly_excluded = 0

t0 = time.time()
for i, s in enumerate(ss):
    weights = np.array([p for p, c in comps])
    omegas = np.array([c * s * (2 * np.pi * T_fixed) for p, c in comps])
    
    for Nc in Nc_list:
        C = CellExtremal(weights, omegas, T_fixed, Nc, N)
        
        for mt in mut_list:
            n_tested += 1
            Q_P = C.Q_P(mt)
            Q_E = C.Q_E(mt)
            violation = Q_P > Q_E * (1 + 1e-5)
            below = mt < C.thr
            
            if violation:
                n_violations += 1
                if below:
                    print(f"VIOLATION BELOW THRESHOLD: s={s:.3f} Nc={Nc} mut={mt:.3f}")
                    print(f"  Q(P)/Q(E) = {Q_P/Q_E:.6f}, thr = {C.thr:.5f}")
                else:
                    n_correctly_excluded += 1
            
            if below:
                n_below_threshold += 1

print(f"\nscan done in {time.time()-t0:.1f}s")
print(f"\n{n_tested} points tested")
print(f"{n_below_threshold} points below threshold (mut < mutilde_*^CE)")
print(f"{n_violations} violations Q(P)>Q(E) found")
print(f"{n_correctly_excluded} violations correctly excluded by threshold")
print(f"\nVIOLATIONS BELOW THRESHOLD (false negatives): {n_violations - n_correctly_excluded}")
if n_violations - n_correctly_excluded == 0:
    print("PASS: All Q(P)>Q(E) violations correctly excluded by mutilde_*^CE(P,T,Nc)")
else:
    print("FAIL: Theorem CE does NOT extend to this multi-atom spectrum")
