"""rev_cu_1: sanity of the independent build + reproduce the note's tables."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "D:/physics/superconductivity/theory/derivations")
import numpy as np
from rev_cu_lib import myH, myOf, uc, RCell, Ocal, brk
from p4_interchange_uniformity import H_mat

N = 300
print("(0) my builder vs certified H_mat")
for w in (0.5, 0.74, 1.0, 2.0, 5.0, 51.2):
    d = np.abs(myH(N, w) - H_mat(N, w, 0)).max()
    d1 = np.abs(myOf(N, w) - w * H_mat(N, w, 1)).max()
    print(f"   w={w:>6}: max|myH-H_mat|={d:.2e}   max|myOf - w*H'|={d1:.2e}")

print("\n(1) K_c = P(O[1]+I)P ?   and Kamke closed form")
O1 = Ocal(N, np.ones(2 * N + 2))
u = 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)
print(f"   max|O[1] - (2uu^T - I)| = {np.abs(O1 - (2*np.outer(u,u)-np.eye(N))).max():.2e}")
for Nc in (3, 5, 10, 20, 37):
    P = np.diag((np.arange(N) < Nc).astype(float))
    Kc = 2 * np.outer(uc(N, Nc), uc(N, Nc))
    print(f"   Nc={Nc:>3d}: max|K_c - P(O[1]+I)P| = {np.abs(Kc - P@(O1+np.eye(N))@P).max():.2e}"
          f"   max|K_c - O[1_(k<Nc)]| = {np.abs(Kc - Ocal(N, (np.arange(2*N+2)<Nc).astype(float))).max():.4f}")
    # Kamke closed form
    Y = np.cumsum(Kc, axis=1); W = Y[:-1, :] - Y[1:, :]
    ue = 1.0/np.sqrt(2.0*np.arange(N+1)+1.0)
    cs = np.concatenate([[0.0], np.cumsum(ue[:N])])
    L = np.arange(1, N)[:, None]; Lp = np.arange(1, N+1)[None, :]
    Wcf = 2.0*cs[np.minimum(Lp, Nc)]*(ue[L-1]*(L-1 < Nc) - ue[L]*(L < Nc))
    offd = (L != Lp)
    print(f"        Kamke: max|W - closed| (offdiag) = {np.abs((W-Wcf)[offd]).max():.2e}"
          f"  min W offdiag = {W[offd].min():+.2e}"
          f"  max row with W>1e-15: {(np.where((W>1e-15)&offd)[0].max()+1) if np.any((W>1e-15)&offd) else 0}")

print("\n(2) threshold table, N=300  (compare to the note Sec.3 table, N=200)")
print(f"   {'w':>6} {'Nc':>4} {'mu*_thr':>10} {'gam':>10} {'c1':>8} {'|uperp|':>8} "
      f"{'tmax':>8} {'q1':>10} {'b':>10} {'m':>10} {'c2':>9} {'mu2-mu3':>9}")
rows = [(0.74,2),(0.74,4),(0.74,10),(1.0,10),(1.0,80),(2.0,10),(5.0,20),(5.0,40),(10.0,40)]
for (w, Nc) in rows:
    C = RCell(w, Nc, N)
    print(f"   {w:>6.2f} {Nc:>4d} {C.thr:>10.5f} {C.gam:>10.5f} {C.c1:>8.4f} {C.up:>8.4f} "
          f"{C.tmax:>8.4f} {C.q1:>10.5f} {C.b:>10.5f} {C.m:>10.5f} {C.c2:>9.2e} {C.mu2-C.mu3:>9.2e}")

print("\n(3) anchors  (N=300)")
anchors = [(3.0,0.13,0.60,10),(2.0,0.13,0.98,10),(2.0,0.10,0.98,10),(1.0,0.13,2.20,10),
           (1.0,0.10,2.00,10),(0.7,0.13,4.56,10),(0.5,0.13,13.5,10),(0.5,0.13,15.6,3),
           (0.4,0.13,51.2,3),(2.0,0.13,3.00,3),(1.0,0.13,2.50,3)]
print(f"   {'lam':>5} {'mu*':>5} {'r':>3} {'w_c':>7} {'Nc':>5} {'mut':>7} {'mut_*':>9} "
      f"{'kappa':>7} {'q(vhat)':>11} {'verdict':>8}")
for (lam, mus, wc, r) in anchors:
    mut = mus/lam; Nc = max(1, int(np.floor(r*wc)))
    if Nc > N//2:
        print(f"   {lam:>5} {mus:>5} {r:>3d} {wc:>7.2f} {Nc:>5d} {mut:>7.4f}   SKIP (Nc>N/2)")
        continue
    C = RCell(wc, Nc, N)
    print(f"   {lam:>5} {mus:>5} {r:>3d} {wc:>7.2f} {Nc:>5d} {mut:>7.4f} {C.thr:>9.5f} "
          f"{C.kappa(mut):>7.4f} {C.q_dense(mut):>11.4e} "
          f"{'PROVED' if mut < C.thr else 'FAILS':>8}")
