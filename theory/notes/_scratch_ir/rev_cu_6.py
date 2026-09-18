"""rev_cu_6: hand-trace of one anchor, Lemma-2/3 chain, mut_crit/mut_* margin scan,
mu_Metz spot table, sign-change count, and the Sec.6 r=3 non-uniqueness numbers."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "D:/physics/superconductivity/theory/derivations")
import numpy as np
from rev_cu_lib import RCell, uc, myH, myOf, Ocal, brk

print("(1) HAND TRACE  wbar=1, N_c=10, mutilde=0.1   (N=300, all numbers recomputed)")
C = RCell(1.0, 10, 300)
tmax = (-C.b + np.sqrt(C.b**2 + C.m*C.q1))/C.m
thr = tmax*C.gam/(2*C.c1*C.up + 2*tmax*C.c1**2)
mut = 0.1
kap = 2*mut*C.c1**2/C.gam
tbd = 2*mut*C.c1*C.up/(C.gam - 2*mut*C.c1**2)
lb  = (1/(1+tbd**2))*(C.q1 - 2*C.b*tbd - C.m*tbd**2)
tact = C.tan_theta_dense(mut); qact = C.q_dense(mut)
print(f"   gam={C.gam:.6f} c1={C.c1:.6f} |u_perp|={C.up:.6f} q1={C.q1:.6f} b={C.b:.6f} m={C.m:.6f}")
print(f"   t_max={tmax:.6f}  mut_*={thr:.6f}  (mut=0.1 < mut_*  => hypothesis holds)")
print(f"   kappa={kap:.6f} (<1 OK)   tan_theta bound={tbd:.6f}  actual tan_theta={tact:.6f}"
      f"   bound holds: {tact<=tbd}")
print(f"   Lemma-3 lower bound={lb:.6f}   actual <vhat,O[f]vhat>={qact:.6f}"
      f"   bound holds: {qact>=lb}")
h = 1e-6
def lam1(w, Nc, mt, NN=300):
    A = myH(NN, w); u = uc(NN, Nc)
    return float(np.linalg.eigvalsh(A - 2*mt*np.outer(u, u))[-1])
fd = (lam1(1.0+h,10,mut) - lam1(1.0-h,10,mut))/(2*h)*1.0
print(f"   Hellmann-Feynman: wbar*dlam1/dwbar (central FD) = {fd:.8f}  vs  q = {qact:.8f}"
      f"   rel {abs(fd-qact)/abs(qact):.2e}   (so T dlam1/dT = -q < 0)")

print("\n(2) margin scan: min over the grid of mut_crit/mut_*   (ratio<=1 would be FATAL)")
def q_of(C, Nc, mt):
    """secular engine (verified == dense eigh to 1e-16 in rev_cu_5b)"""
    c2 = C.c**2
    f = lambda l: 2*mt*float(np.sum(c2/(C.e - l))) - 1.0
    lo, hi = C.mu2 + 1e-18, C.mu1 - 1e-18
    if f(lo) > 0:
        u = uc(C.N, Nc)
        e, V = np.linalg.eigh(C.A - 2*mt*np.outer(u, u))
        return float(V[:, -1] @ (C.Of @ V[:, -1]))
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if f(mid) > 0: hi = mid
        else: lo = mid
    y = C.c/(C.e - 0.5*(lo+hi)); y = y/np.linalg.norm(y)
    return float(y @ (C.G @ y))
NN = 250
worst = (np.inf, None)
t0 = time.time()
for w in (0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.74,0.9,1.0,1.5,2.0,3.0,5.0,10.0,20.0):
    for Nc in (1,2,3,4,5,6,8,10,15,20,30,50,80,120):
        C = RCell(w, Nc, NN)
        grid = np.concatenate([np.linspace(1e-4, 1.0, 120), np.linspace(1.0, 3.0, 40)])
        mc = np.inf; prev = 1e-4
        for mt in grid:
            if q_of(C, Nc, mt) <= 0:
                a_, b_ = prev, mt
                for _ in range(45):
                    mid = 0.5*(a_+b_)
                    if q_of(C, Nc, mid) > 0: a_ = mid
                    else: b_ = mid
                mc = a_; break
            prev = mt
        ratio = mc/C.thr
        if ratio < worst[0]: worst = (ratio, (w, Nc, C.thr, mc))
print(f"   min mut_crit/mut_* = {worst[0]:.4f} at (wbar,N_c)={worst[1][:2]}"
      f"  mut_*={worst[1][2]:.5f} mut_crit={worst[1][3]:.5f}   ({time.time()-t0:.0f}s)")

print("\n(3) mu_Metz spot check (wbar=1)  and sign-change count of the pencil top vector")
NN = 200
def kamke(M):
    Y = np.cumsum(M, axis=1); return Y[:-1,:] - Y[1:,:]
L = np.arange(1, NN)[:, None]; Lp = np.arange(1, NN+1)[None, :]; offd = (L != Lp)
WF = kamke(myH(NN, 1.0))
for Nc in (10, 40, 80):
    Wc = kamke(2*np.outer(uc(NN,Nc), uc(NN,Nc)))
    msk = offd & (Wc > 1e-14)
    r = np.where(msk, WF/np.where(Wc>1e-14, Wc, 1.0), np.inf)
    print(f"   wbar=1 N_c={Nc:>3d}: mu_Metz={r.min():.4e}   (note: 1.55e-3 / 2.65e-5 / 3.40e-6)")
for NN2 in (200, 400, 800):
    A = myH(NN2, 1.0); u = uc(NN2, 10)
    e, V = np.linalg.eigh(A - 2*0.26*np.outer(u,u)); v = V[:,-1]
    a = np.abs(v).max(); nz = v[np.abs(v) > 1e-12*a]
    print(f"   N={NN2}: sign changes of top eigenvector (wbar=1,Nc=10,mut=0.26) = "
          f"{int((np.sign(nz[:-1])!=np.sign(nz[1:])).sum())}")

print("\n(4) Sec.6: r=3, mutilde=0.10 -- the lambda=1 non-uniqueness window")
NN = 400
for w in (1.996, 2.000, 2.004, 2.008, 2.012, 2.016, 2.020, 2.100, 2.200, 2.333, 2.400):
    Nc = int(np.floor(3*w))
    print(f"   wbar={w:.3f} N_c={Nc:>3d}  lam1={lam1(w, Nc, 0.10, NN):.6f}")
