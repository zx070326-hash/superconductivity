"""rev3 item 6/1: the committed CFG default is kna=14 but the log was run with kna=16.
What would the envelope terms be at the committed default on the binding cells (no Cholesky,
assume c=1e-7 as at kna=16)?  Also K_* from the FULL-precision server node values."""
import sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
from mpmath import mp
from p4_interchange_uniformity import spec_seg, schur_coupling, schur_tailblock
mp.dps = 30
for a, b, Na16, Lw, slack, al16 in ((59.5399, 59.8827, 960, 0.7515, 0.1297, 0.590), (59.8827, 60.0, 960, 0.7515, 0.0441, 0.504), (20.0, 20.1152, 340, 1.2262, 0.0706, 0.454)):
    sp = spec_seg(a, b, 0)
    for Na in (Na16, int(-(-14 * b // 20)) * 20):
        nB = float(schur_coupling(Na, b, 0, spec=sp)); nD = float(schur_tailblock(Na, b, 0, spec=sp))
        al = 1e-7 + slack + nB + nD + 1e-12
        print(f"cell [{a},{b}] Na={Na}: ||B||={nB:.4f} ||D||={nD:.4f} slack={slack} -> alpha_cell~{al:.4f}  g(0.33)~{0.7856 if a>50 else 0.7772:.4f}-{al:.4f}={ (0.7856 if a>50 else 0.7772)-al:+.4f}   (log al at kna=16: {al16})", flush=True)
print("K_* from full-precision node values (server node_job rerun):")
for w, g, B2, B4 in ((30.0, 1.6624286882410773 - 0.4246470263874016, 0.04878288142105583, 0.009924257737949577), (60.0, 1.6656031892425494 - 0.45955797242667923, 0.023698692501203204, 0.001265490547282798)):
    K = w**4 * B4 * 4 / 24 + 0.5 * w**4 * B2**2 / g; s = min(0.25, 0.95 * (g / (2 * w**2 * B2))**0.5)
    print(f"  w={w}: g_alt={g:.5f} K_*={K:.4e} s_max={s:.4f}")
