"""4a: exact concavity criterion H'' = -2 Psi + R < 0 via the identity
      R = 2 < y_perp, (h-A)^{-1} y_perp >,  y_perp = (I - v v^T) O[f'] v,
valid for mu~ > 0 (at mu~ = 0 the j=1 term is a removable 0/0: h = mu_1, a_1 = 0).
(a) validation: mu~ -> 0 must reproduce the proven mu*=0 concavity H'' < 0;
(b) physical boundary mu~_c(varpi, r) with N_c = floor(r*varpi).
numerical-evidence (float64, N as stated)."""
import os, math, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concavity_exact_resolvent import build

_cache = {}


def hpp(N, x, Nc, mut):
    key = (N, x)
    if key not in _cache:
        n, u, A, Mp, Opsi = build(N, x)
        _cache[key] = (n, u, A, Mp, Opsi) + eigh(A)
    n, u, A, Mp, Opsi, mu, PSI = _cache[key]
    uc = np.where(n < Nc, u, 0.0)
    M = A - 2.0*mut*np.outer(uc, uc) if mut > 0 else A
    lam, W = eigh(M)
    h = lam[-1]; v = W[:, -1]
    Psi = float(v @ (Opsi @ v))
    y = Mp @ v; yp = y - float(v @ y)*v
    a = PSI.T @ yp
    d = h - mu
    if mut == 0.0:                       # removable: drop the j=1 (0/0) term
        R = float(np.sum(2.0*a[:-1]**2/d[:-1]))
    else:
        R = float(np.sum(2.0*a**2/d))
    return -2.0*Psi + R, R, Psi


print("(a) mu~ -> 0 validation: the PROVEN mu*=0 concavity (H'' < 0) must be recovered")
print(f"{'x':>6} {'Nc':>3} {'mu~':>8} {'H\"':>12} {'R':>11} {'2Psi':>9}")
bad = 0
for x in (0.04, 0.09, 0.25, 1.0, 4.0):
    for mut in (0.0, 1e-5, 1e-4, 1e-3):
        H, R, P = hpp(400, x, 6, mut)
        flag = "  <-- H''>0" if H > 0 else ""
        bad += H > 0
        print(f"{x:>6} {6:>3} {mut:>8.0e} {H:>12.4e} {R:>11.3e} {2*P:>9.4f}{flag}")
    print()
print(f"validation: {bad} sign violations at mu~ <= 1e-3 (expect 0)\n")

print("(b) concavity boundary mu~_c(varpi, r): H'' < 0 for mu~ < mu~_c")
print("    ('all' = concave for every mu~ in [1e-3, 1.0]; 'none' = fails already at 1e-3)")
VPS = [0.3, 0.4, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 4.5, 6.0]
GRID = np.geomspace(1e-3, 1.0, 61)
for r in (3, 10):
    print(f"\n  r = {r}   (N_c = floor(r*varpi))")
    print(f"  {'varpi':>6} {'Nc':>4} {'N':>5} {'mu~_c':>9} {'H\" @1e-3':>11} "
          f"{'H\" @0.13':>11} {'H\" @0.5':>11}  {'signs over grid':>16}")
    for vp in VPS:
        x = vp*vp; Nc = max(1, int(r*vp)); N = max(400, 8*Nc)
        sg = np.array([hpp(N, x, Nc, m)[0] > 0 for m in GRID])
        if not sg.any():
            mc = "all"
        elif sg[0]:
            mc = "none"
        else:
            mc = f"{GRID[np.argmax(sg)]:.4f}"
        nfail = int(sg.sum())
        blocks = "".join("X" if s else "." for s in sg[::4])
        print(f"  {vp:>6.1f} {Nc:>4} {N:>5} {mc:>9} "
              f"{hpp(N,x,Nc,1e-3)[0]:>11.3e} {hpp(N,x,Nc,0.13)[0]:>11.3e} "
              f"{hpp(N,x,Nc,0.5)[0]:>11.3e}  {blocks:>16} ({nfail}/61)")
