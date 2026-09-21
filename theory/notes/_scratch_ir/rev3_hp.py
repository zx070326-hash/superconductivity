"""rev3 item 1/2: weighted-row-0 L (server recompute) vs float64 max ||H'(w)|| at N=3000 on the three spot cells."""
import numpy as np, sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from rev2_truth_lib import H
for a, b, Lw in ((20.0, 20.1152, 1.2262), (34.9065, 35.1075, 0.9626), (59.8827, 60.0, 0.7515)):
    m = max(np.max(np.abs(np.linalg.eigvalsh(H(3000, w, 1)))) for w in np.linspace(a, b, 5))
    print(f"[{a},{b}] L_w0={Lw}  max||H'|| N=3000 = {m:.5f}  ratio {Lw/m:.1f}x  {'OK' if Lw >= m else 'FATAL'}", flush=True)
