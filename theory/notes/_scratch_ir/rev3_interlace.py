"""rev3 item 4: lam1(H_N) vs N at the X5 nodes -- how much room does the Rayleigh mu at support Nv have
against lam1 at N=3000?  (float64 truth, self-written H)"""
import numpy as np, sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from rev2_truth_lib import H
for w, Nv in ((30.0, 960), (40.0, 1280), (60.0, 1920)):
    vals = {N: np.linalg.eigvalsh(H(N, w))[-1] for N in (Nv // 2, Nv, 2 * Nv, 3000, 3840)}
    print(f"w={w:>4.0f} Nv={Nv}: " + "  ".join(f"lam1(N={N})={v:.9f}" for N, v in vals.items()) + f"  | lam1(3000)-lam1(Nv)={vals[3000]-vals[Nv]:+.2e}", flush=True)
