"""Is lambda_1(M) monotone DECREASING in N_c?  v3's docstring assumes it
("larger N_c -> stronger Coulomb -> lower Tc") to justify its directional rounding.
Checked at FIXED T, so no bisection is involved.
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
from scipy.linalg import eigh
from _main_indep_baseline import read_a2f, P_measure, F_T, O_of_f, A2F

w,a2f,_ = read_a2f(A2F); lam,p,om2 = P_measure(w,a2f)
T, N, mus = 17.146237, 100, 0.13
F = F_T(2*N, w, p, T); n,u,L1,A = O_of_f(F, N)
vp = math.sqrt(om2)/(2*math.pi*T)
print(f"T={T} meV  varpi={vp:.6f}  r*varpi={10*vp:.6f}  (floor=13, ceil=14)")
_,V = eigh(A)
v = V[:,-1]
print(f"\ntop eigenvector of A: v[11..16] = {np.array2string(v[11:17], precision=4)}")
print(f"{'N_c':>5}{'lam1(M)':>14}{'d lam1':>12}{'|<u_c,v>|':>12}")
prev=None
for Nc in range(10, 19):
    uc = np.where(n < Nc, u, 0.0)
    M = A - 2.0*(mus/lam)*np.outer(uc, uc)
    l1 = eigh(M, eigvals_only=True)[-1]
    d = "" if prev is None else f"{l1-prev:+12.3e}"
    print(f"{Nc:>5}{l1:14.9f}{d:>12}{abs(uc@v):12.6f}")
    prev=l1
print("\nv3 docstring assumes lam1 DECREASES with N_c. "
      "Any '+' in the d-lam1 column falsifies that.")

print("\n=== full N_c sweep (same T, N=100) ===")
vals=[]
for Nc in list(range(0, 101, 1)):
    uc = np.where(n < Nc, u, 0.0)
    M = A - 2.0*(mus/lam)*np.outer(uc, uc)
    vals.append(eigh(M, eigvals_only=True)[-1])
vals=np.array(vals)
d=np.diff(vals)
print(f"lam1(N_c=0)  = {vals[0]:.9f}   (no Coulomb = lam1(A))")
print(f"lam1 min     = {vals.min():.9f} at N_c={vals.argmin()}")
print(f"lam1(N_c=100)= {vals[-1]:.9f}")
print(f"increments: {np.sum(d<0)} negative (down-jumps), {np.sum(d>0)} positive (up-jumps)")
print("first 12 increments:", np.array2string(d[:12], precision=3))
print("N_c where the sign of d(lam1) flips:",
      [i+1 for i in range(len(d)-1) if (d[i]>0) != (d[i+1]>0)])
