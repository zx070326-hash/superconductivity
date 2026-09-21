"""rev3 item 3: mutilde=0.33 pencil.  (a) true gap of M(w,0.33,N_c) at N=1500 over ALL N_c
(rank-one secular sweep, self-written, cross-checked vs direct eigvalsh);  (b) retained-cell
2x2 bound (mu2_pencil_cell, imported = the claim) vs self-written 2x2 truth at 11 interior
points x 3 window cases, on cells sampled from the three grids used by X4/X4L."""
import numpy as np, time, sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
from rev2_truth_lib import H, brk
from mpmath import mp
from p4d_coulomb_certify import mu2_pencil_cell
from p4_interval_extension import BASE_RATIO
import p4_interchange_uniformity as piu
mp.dps = 30
t0 = time.time(); MT = 0.33; RHO = 2*MT; N = 1500
def u_vec(N, Nc):
    u = np.zeros(N); k = min(Nc, N); u[:k] = 1/np.sqrt(2*np.arange(k)+1.0); return u
def all_windows(N, w):
    d, Q = np.linalg.eigh(H(N, w)); d = d[::-1]; Q = Q[:, ::-1]
    Upref = np.cumsum(Q*(1/np.sqrt(2*np.arange(N)+1.0))[:,None], axis=0)
    Z2 = np.vstack([np.zeros(N), Upref**2])
    out = np.empty((N+1, 2))
    for k in (0, 1):
        lo = np.full(N+1, d[k+1])+1e-13; hi = np.full(N+1, d[k])-1e-13
        for _ in range(70):
            mid = 0.5*(lo+hi); f = 1-RHO*np.sum(Z2/(d[None,:]-mid[:,None]), axis=1)
            up = f > 0; lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
        out[:, k] = 0.5*(lo+hi)
    out[0] = d[:2]; return out, d
print("(a) true gap over ALL N_c=0..1500, mutilde=0.33, N=1500")
print(f"{'w':>6} {'gap(H)':>8} {'inf_Nc gap(M)':>14} {'argmin':>6} {'gap@3w':>8} {'gap@10w':>8} {'xcheck':>8}")
for w in (1.5, 1.5666, 1.7, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 45.0, 60.0):
    out, d = all_windows(N, w); gap = out[:,0]-out[:,1]; j = int(np.argmin(gap))
    xc = 0.0
    for Nc in (1, 2, int(3*w), int(10*w), 1500):
        u = u_vec(N, Nc); E = np.linalg.eigvalsh(H(N, w)-RHO*np.outer(u, u)); xc = max(xc, abs((E[-1]-E[-2])-gap[Nc]))
    print(f"{w:>6.3f} {d[0]-d[1]:>8.4f} {gap.min():>14.5f} {j:>6} {gap[int(3*w)]:>8.4f} {gap[int(10*w)]:>8.4f} {xc:>8.1e}", flush=True)
print("\n(b) retained cells: mu2_pencil_cell(a,b,0.33) vs true lam1(M^(2)(w)) (self-written), 11 interior w x 3 cases")
def lam1_M2(w, case):
    b1, b2, b3 = brk(1, w), brk(2, w), brk(3, w)
    h = np.array([[b1, (b1+b2)/np.sqrt(3)], [(b1+b2)/np.sqrt(3), (b3-2*b1)/3]])
    u = {0: np.array([0., 0.]), 1: np.array([1., 0.]), 2: np.array([1., 1/np.sqrt(3)])}[case]
    return np.linalg.eigvalsh(h-2*MT*np.outer(u, u))[-1]
def grid_on(lo, hi, refine):
    r = BASE_RATIO**(1.0/refine); g = [lo]
    while g[-1]*r < hi-1e-12: g.append(g[-1]*r)
    g.append(hi); return [float(x) for x in g]
g1 = grid_on(20.0, 60.0, 2); g2 = grid_on(5.0, 20.0, 1); g3 = [float(x) for x in piu.wbar_grid()]
cells = [(g1[i], g1[i+1]) for i in range(0, len(g1)-1, 8)] + [(g1[-2], g1[-1])] \
      + [(g2[i], g2[i+1]) for i in range(0, len(g2)-1, 10)] + [(g2[-2], g2[-1])] \
      + [(g3[i], g3[i+1]) for i in range(len(g3)-1) if g3[i] >= 1.5666-1e-9][::8]
nviol = 0; worst = (np.inf, None); nks_fail = 0
for a, b in cells:
    m, k1, k2, k3, case = mu2_pencil_cell(a, b, MT); m = float(m); ks = k1 >= 0 and k2 >= 0 and k3 >= 0
    nks_fail += (not ks)
    for w in np.linspace(a, b, 11):
        for c in (0, 1, 2):
            v = lam1_M2(w, c) - m
            nviol += (v < 0)
            if v < worst[0]: worst = (v, (a, b, w, c))
print(f"  cells={len(cells)} points={len(cells)*33} violations={nviol} k-sign fails={nks_fail}")
print(f"  min(true lam1(M2) - bound) = {worst[0]:.2e} at cell [{worst[1][0]:.4f},{worst[1][1]:.4f}] w={worst[1][2]:.4f} case={worst[1][3]}")
print(f"({time.time()-t0:.0f}s)")
