"""Is the ceiling's anharmonicity-invariance the Hopfield sum rule?
lam*<om^2> = 2 int a2F(om) om dom =: eta/M  -- purely electronic, no phonon freqs.
Check on first-hand data: conserved under harmonic->anharmonic? scales as 1/M?
Computed two ways (via lam,<om2> and by direct quadrature of 2 int a2F om dom).
numerical-evidence (float64)."""
import os, math, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_fill_factor import moments, ROOT, MEV_K
from fill_factor_probe import CINF

vals = {}
print(f"{'dataset':>26} {'lam':>7} {'sq<om2>':>9} {'lam<om2>':>11} {'2int a2F om':>12} {'rel':>9} {'ceil[K]':>8}")
for tag in ("H3S-200GPa", "D3S-200GPa", "H3S-160GPa", "D3S-160GPa"):
    for sub in ("EPW_anharmonic", "EPW_harmonic"):
        gl = [x for x in glob.glob(f"{ROOT}/{tag}/{sub}/*.a2f") if not x.endswith(".a2fv")]
        if not gl: continue
        lam, om2, olog, w, a2f = moments(gl[0])
        A = lam*om2                                    # via the two moments
        B = 2.0*np.trapezoid(a2f*w, w)                 # direct sum rule integral
        vals[(tag, sub)] = A
        print(f"{tag+'/'+sub:>26} {lam:>7.4f} {math.sqrt(om2):>9.3f} {A:>11.1f} {B:>12.1f} "
              f"{abs(A-B)/B:>9.1e} {CINF*math.sqrt(A)*MEV_K:>8.2f}")

print("\nanharmonic vs harmonic (eta/M must be phonon-spectrum blind):")
for tag in ("H3S-200GPa", "D3S-200GPa", "H3S-160GPa"):
    k1, k2 = (tag,"EPW_anharmonic"), (tag,"EPW_harmonic")
    if k1 in vals and k2 in vals:
        print(f"  {tag:>12}: anh={vals[k1]:.1f}  har={vals[k2]:.1f}  ratio={vals[k1]/vals[k2]:.5f}")

print("\nisotope scaling (eta/M ~ 1/M; M_D/M_H = 2 on the H sublattice only):")
for sub in ("EPW_anharmonic", "EPW_harmonic"):
    k1, k2 = ("H3S-200GPa",sub), ("D3S-200GPa",sub)
    if k1 in vals and k2 in vals:
        print(f"  {sub:>15}: (eta/M)_H/(eta/M)_D = {vals[k1]/vals[k2]:.4f}   (pure H-sublattice: 2.0)")
