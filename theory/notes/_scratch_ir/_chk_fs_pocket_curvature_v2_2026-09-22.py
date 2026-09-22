"""v2: check curvature of EVERY band (not just nearest-to-mu) at each high-symmetry
point, to find which bands actually have an extremum poking through mu there."""
import importlib.util, os
import numpy as np

_here = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    "graser1003_hamiltonian", os.path.join(_here, "_graser1003_hamiltonian_2026-09-22.py"))
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
build_H_graser = _mod.build_H_graser
NBAND = 5

def bands(kx, ky):
    return np.linalg.eigvalsh(build_H_graser(kx, ky, 0.0))

mu = 0.0109  # from v1, N=128, 6 e/Fe

pts = {"Gamma": (0.0, 0.0), "X=(pi,0)": (np.pi, 0.0), "Y=(0,pi)": (0.0, np.pi),
       "M=(pi,pi)": (np.pi, np.pi)}
h = 0.02
for name, (kx0, ky0) in pts.items():
    e0 = bands(kx0, ky0)
    print(f"\n{name}: bands (E-mu) = {np.round(e0-mu,4)}")
    for iband in range(NBAND):
        def eb(kx, ky, n=iband):
            return bands(kx, ky)[n]
        e_val = eb(kx0, ky0)
        d2x = (eb(kx0+h, ky0) - 2*e_val + eb(kx0-h, ky0)) / h**2
        d2y = (eb(kx0, ky0+h) - 2*e_val + eb(kx0, ky0-h)) / h**2
        curv = d2x + d2y
        tag = ""
        if abs(e_val - mu) < 1.2:
            if e_val > mu and curv < 0:
                tag = "  <- HOLE-pocket candidate (max above mu)"
            elif e_val < mu and curv > 0:
                tag = "  <- ELECTRON-pocket candidate (min below mu)"
        print(f"  band#{iband}: E-mu={e_val-mu:+.4f}  curv={curv:+.4f}{tag}")
