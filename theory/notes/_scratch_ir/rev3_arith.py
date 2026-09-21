"""rev3 items 2 + 4c: grid arithmetic, head sizes (kna 14 vs 16), K_*/s_max from logged numbers."""
import numpy as np, sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
from p4_interval_extension import BASE_RATIO, ext_grid
def grid_on(lo, hi, refine):
    r = BASE_RATIO**(1.0/refine); g = [lo]
    while g[-1]*r < hi-1e-12: g.append(g[-1]*r)
    g.append(hi); return [float(x) for x in g]
def head_size(b, kna):
    if b <= 5+1e-12: return 70
    if b <= 10+1e-12: return 140
    if b <= 20+1e-12: return 280
    return int(np.ceil(kna*b/20.0)*20)
g = grid_on(20, 60, 2); r = BASE_RATIO**0.5
print(f"ratio={r:.5f} BASE_RATIO={BASE_RATIO:.6f}; [20,60]: {len(g)} points {len(g)-1} cells; g[0]={g[0]} g[-1]={g[-1]}; last full cell [{g[-3]:.4f},{g[-2]:.4f}] h={g[-2]-g[-3]:.4f}; final cell [{g[-2]:.4f},{g[-1]:.4f}] h={g[-1]-g[-2]:.4f}")
print(f"  cell 192 in log [20.000,20.115]: g[1]={g[1]:.4f}; cell 2 [59.540,59.883]: g[-3]={g[-3]:.4f} g[-2]={g[-2]:.4f}")
print(f"  grid points sitting exactly on a Na band boundary (16b/20 integer): {[x for x in g if abs(16*x/20-round(16*x/20))<1e-9]}")
g520 = grid_on(5, 20, 1); e = ext_grid(5.0, 20.0)
print(f"[5,20] refine=1: {len(g520)-1} cells; identical to X1 ext_grid: {len(g520)==len(e) and np.allclose(g520, e)}; ends {g520[0]},{g520[-1]}")
print("head sizes kna=16 (log) vs kna=14 (file default CFG):")
for b in (20.115, 30.0, 40.0, 59.883, 60.0):
    print(f"  b={b}: Na16={head_size(b,16)}  Na14={head_size(b,14)}")
print("X5 K_*/s_max from logged g_alt/B2/B4 (4-digit rounded):")
for w, g_, B2, B4, K, s in ((30, 1.23778, 0.0488, 0.0099, 2.1184e3, 0.1128), (40, 1.22569, 0.0360, 0.0042, 3.1549e3, 0.0980), (60, 1.20605, 0.0237, 0.0013, 5.7510e3, 0.0799)):
    Kc = w**4*B4*4/24 + 0.5*w**4*B2**2/g_; sc = min(0.25, 0.95*np.sqrt(g_/(2*w**2*B2)))
    print(f"  w={w}: K recomputed={Kc:.4e} (log {K:.4e}, rel {abs(Kc-K)/K:.1e}); s_max={sc:.4f} (log {s}); cond w^2 B2 s^2={w**2*B2*sc**2:.4f} < g/2={g_/2:.4f}: {w**2*B2*sc**2 < g_/2}")
    print(f"     K_tail={w**4*B4*4/24:.4e}  K_2nd={0.5*w**4*B2**2/g_:.4e}  (B4 log rounding +-5e-5 -> K_tail +-{w**4*5e-5*4/24:.1e})")
