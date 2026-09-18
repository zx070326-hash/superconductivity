"""Can (C**) close with ANALYTIC gap + norm bounds instead of certificates?

H''(x) = <v,O[fxx]v> + 2 sum_{j>1} <v,O[fx]w_j>^2/(h-lam_j)
      <= <v,O[fxx]v> + 2||(1-P1)O[fx]v||^2 / g,   g := h - lam_2.

Analytic ingredients available:
  (G) gap: 2026-09-13_gap_lorentz_probe.md gives lam_2 <= sqrt(tr B^2 - h_lo^2)
      with B from the general-kernel rearrangement identity -> lam_2/h <= 0.27..0.58
      on varpi in [0.2,2]; route R gives h - Theta* for varpi >= 8 (Theta*=0.6349).
  (N) ||O[fx]|| via Schur/Holmgren closed form.
Question: with g replaced by the analytic lower bound, does term1 + term2_bound < 0?
Report the slack for each route so we know which range needs what.
"""
import numpy as np, sys
sys.path.insert(0, r"D:/physics/superconductivity/theory/derivations")

def O_of_kernel(N, fvals):
    n = np.arange(N); u = 1.0/np.sqrt(2*n+1.0)
    i = n[:,None]; j = n[None,:]
    idx_d = np.abs(i-j)-1
    A = np.where(i!=j, fvals[np.clip(idx_d,0,None)], 0.0)
    A = A + fvals[i+j]
    A = A*u[:,None]*u[None,:]
    csum = np.concatenate(([0.0], np.cumsum(fvals)))
    A[n,n] -= 2.0*u**2*csum[n]
    return A

def schur(A):
    r = np.abs(A).sum(1)
    return np.sqrt(r.max()*r.max())

THETA = 0.6349
rows=[]
for x in [0.04,0.1,0.25,1.0,4.0,9.0,25.0,100.0]:
    N = max(160, int(40*np.sqrt(x))+160); K=2*N+2
    k = np.arange(1,K+1,dtype=float)
    A   = O_of_kernel(N, x/(k**2+x))
    Ax  = O_of_kernel(N, k**2/(k**2+x)**2)
    Axx = O_of_kernel(N, -2.0*k**2/(k**2+x)**3)
    w,V = np.linalg.eigh(A); o=np.argsort(w)[::-1]; w=w[o]; V=V[:,o]
    h=w[0]; l2=w[1]; v=V[:,0]; g=h-l2
    t1 = v@(Axx@v)
    r  = Ax@v; r1 = r - (v@r)*v
    num=(V.T@(Ax@v))**2
    t2_true = 2.0*np.sum(num[1:]/(h-w[1:]))
    t2_perp = 2.0*np.linalg.norm(r1)**2/g          # collapse to gap
    t2_nrm  = 2.0*schur(Ax)**2/g                   # collapse to norm
    g_gap   = h*(1-0.58) if x<=4 else max(h-THETA,1e-9)   # analytic gap route
    t2_ag   = 2.0*np.linalg.norm(r1)**2/g_gap
    t2_full = 2.0*schur(Ax)**2/g_gap
    rows.append((x,h,g,g_gap,t1,t2_true,t2_perp,t2_ag,t2_full,t1+t2_ag,t1+t2_full))

print(f"{'x':>7} {'h':>7} {'g_true':>7} {'g_anal':>7} {'term1':>10} {'t2_true':>8} {'t2_perp':>8} {'t2:perp+ag':>10} {'t2:norm+ag':>10} {'H\"(perp,ag)':>12} {'H\"(norm,ag)':>12}")
for r in rows:
    print(f"{r[0]:7.2f} {r[1]:7.4f} {r[2]:7.4f} {r[3]:7.4f} {r[4]:10.5f} {r[5]:8.5f} {r[6]:8.5f} {r[7]:10.5f} {r[8]:10.4f} {r[9]:12.5f} {r[10]:12.4f}")
print("\nCLOSES (<0) with perp-norm + analytic gap:", all(r[9]<0 for r in rows))
print("CLOSES (<0) with Schur norm + analytic gap:", all(r[10]<0 for r in rows))
