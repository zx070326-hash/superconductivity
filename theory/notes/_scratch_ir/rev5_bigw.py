"""rev5_bigw.py -- item 6: does (x^2 E^2/g)/F saturate ~0.21 as wbar -> inf?  N ~ 40 wbar, dense eigh."""
import sys, time, numpy as np
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from rev5_truth import Hb, Phi, Psi, fd_hpp
t0 = time.time()
print(f"{'w':>6} {'N':>5} {'lam1':>10} {'lam2':>10} {'g':>8} {'w^4 F':>10} {'w^4 E^2':>10} {'w^4 E':>10} {'x^2E^2/g /F':>12} {'exact pos/F':>12} {'x^2 H''':>10} {'x^2 H''(fd)':>12}")
for w, N in ((44.0, 1800), (60.0, 2400), (60.0, 3600), (100.0, 4000), (130.0, 5200)):
    x = w * w
    e, V = np.linalg.eigh(Hb(N, w)); v = V[:, -1]
    F = float(v @ Phi(N, x) @ v)
    y = Psi(N, x) @ v; cj = V.T @ y; E = float(np.linalg.norm(y - float(v @ y) * v))
    pos_exact = x * x * float(np.sum(cj[:-1] ** 2 / (e[-1] - e[:-1])))
    g = e[-1] - e[-2]
    hpp = -2 * F + 2 * pos_exact
    fd = fd_hpp(N, w, 2e-3) if N <= 4000 else float('nan')
    print(f"{w:6.1f} {N:5d} {e[-1]:10.6f} {e[-2]:10.4e} {g:8.4f} {w**4*F:10.5f} {w**4*E*E:10.5f} {w**4*E:10.5f} {x*x*E*E/g/F:12.4f} {pos_exact/F:12.4f} {x*x*hpp:10.5f} {x*x*fd:12.5f}   [{time.time()-t0:.0f}s]", flush=True)
