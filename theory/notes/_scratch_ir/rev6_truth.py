# rev6: own float64 truth for M1 (items 3, 5, 7). Own builder from the KAY II formulas; script imported only for its bounds().
import sys, os, time, numpy as np
from scipy.sparse.linalg import eigsh, LinearOperator
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
import p4_cii_m1 as m1
import p4_interchange_uniformity as piu
t0 = time.time()
def opK(N, f):           # f: array indexed by integer k >= 0 (f[0] unused), length >= 2N
    n = np.arange(N); sq = np.sqrt(2.0*n+1.0)
    D = np.abs(n[:,None]-n[None,:]); S = n[:,None]+n[None,:]+1
    M = (np.where(D>0, f[D], 0.0) + f[S])/(sq[:,None]*sq[None,:])
    cs = np.concatenate(([0.0], np.cumsum(2.0*f[1:N])))
    M[np.diag_indices(N)] -= cs/(2.0*n+1.0)
    return M
def tab(N, fun, w):
    k = np.arange(0, 2*N+2, dtype=float); f = np.zeros_like(k); f[1:] = fun(k[1:]/w); return f
Gf = lambda y: 1/(1+y*y); PHf = lambda y: y*y/(1+y*y)**3; KAf = lambda y: y*y/(1+y*y)**2
def RG(s,t): return 2*s*s*t*t*(5+s*s+t*t)/((1+(s-t)**2)*(1+(s+t)**2)*(1+s*s)*(1+t*t))
# one-time check of own builder against certified H_mat
for w in (3.0, 45.0):
    print(f"[0] own builder vs piu.H_mat N=60 w={w}: {np.max(np.abs(opK(60, tab(60,Gf,w)) - piu.H_mat(60,w,0))):.2e}")
m1.SPL = m1.build_profiles(); print(f"[0] profiles built ({time.time()-t0:.0f}s)")
def topk(M, k=3):
    if M.shape[0] <= 2000:
        e = np.linalg.eigh(M); return e[0][::-1][:k], e[1][:, ::-1][:, :k]
    e, V = eigsh(M, k=k, which='LA', tol=1e-12); o = np.argsort(e)[::-1]; return e[o], V[:, o]
print("\n[7] float64 truth (N=40w) vs M1 bounds.  side: OK means bound on the correct side")
hdr = " w   |   h      rho_lo  |  lam2    Theta*  ||Rop||  lam1(Rop-Dg) supSchur | ||r||  C_r  | sin   s_DK | F(v)   F_lo  | E(v)   E_hi | g    g_lo | ratio_true ratio_M1"
print(hdr)
rows = {}
for w in (8.0, 20.0, 45.0, 200.0):
    N = int(40*w); n = np.arange(N); wn = 1/(2*n+1.0); sn = (n+0.5)/w; p = np.sqrt(wn)
    H = opK(N, tab(N, Gf, w)); OP = opK(N, tab(N, PHf, w)); OK = opK(N, tab(N, KAf, w))
    e, V = topk(H, 3); h, l2 = e[0], e[1]; v = V[:,0]*np.sign(V[0,0])
    tau = p*Gf(sn); B = tau@tau; th = tau/np.sqrt(B); rho = th@H@th
    r = H@th - rho*th
    S, T = np.meshgrid(sn, sn, indexing='ij'); Rop = np.outer(p,p)*RG(S,T)
    Dg = wn + 2*np.concatenate(([0.0], np.cumsum(Gf(np.arange(1,N)/w))))/(2*n+1.0)
    M1 = 2*np.outer(tau,tau) + Rop - np.diag(Dg)
    idn = np.max(np.abs(M1-H))
    nR = topk(Rop,1)[0][0]; lRD = topk(Rop-np.diag(Dg),1)[0][0]
    schur = np.max(np.sum(RG(S,T)*wn[None,:], axis=1))
    Fv = v@OP@v; y = OK@v; y -= (v@y)*v; Ev = np.linalg.norm(y)
    sin_ = np.sqrt(max(0, 1-(th@v)**2))
    # true operator norms of O[phi], K (largest |eig|)
    nOP = np.max(np.abs(np.linalg.eigvalsh(OP))) if N<=2000 else abs(eigsh(OP,k=1,which='LM')[0][0])
    nOK = np.max(np.abs(np.linalg.eigvalsh(OK))) if N<=2000 else abs(eigsh(OK,k=1,which='LM')[0][0])
    b = m1.bounds(w)
    g = h-l2
    rows[w] = dict(h=h,l2=l2,rho=rho,B=B,nr=np.linalg.norm(r),nR=nR,lRD=lRD,schur=schur,Fv=Fv,Ev=Ev,sin=sin_,g=g,nOP=nOP,nOK=nOK,b=b,idn=idn)
    ok = lambda c: "OK" if c else "VIOL"
    print(f"{w:4g} | {h:.4f} {b['rho']:.4f} | {l2:.4f} {m1.THETA_STAR:.4f} {nR:.4f} {lRD:.4f} {schur:.4f} | {np.linalg.norm(r):.4f} {b['Cr']:.4f} | {sin_:.4f} {b['s']:.4f} | {Fv:.5f} {b['F_lo']:.5f} | {Ev:.5f} {b['E_hi']:.5f} | {g:.3f} {b['g']:.3f} | {Ev**2/(g*Fv):.4f} {b['ratio']:.4f}")
    print(f"      (M1) identity max|.|={idn:.1e}; sides: rho_lo<=h {ok(b['rho']<=h)}, rho_lo<=rho {ok(b['rho']<=rho)}, lam2<=||Rop|| {ok(l2<=nR)}, ||Rop||<=Theta* {ok(nR<=m1.THETA_STAR)}, "
          f"Cr>=||r|| {ok(b['Cr']>=np.linalg.norm(r))}, s>=sin {ok(b['s']>=sin_)}, F_lo<=F(v) {ok(b['F_lo']<=Fv)}, E_hi>=E(v) {ok(b['E_hi']>=Ev)}, g_lo<=g {ok(b['g']<=g)}, "
          f"NP>=||O[phi]|| {ok(b['NP']>=nOP)} ({b['NP']:.4f} vs {nOP:.4f}), NK>=||K|| {ok(b['NK']>=nOK)} ({b['NK']:.4f} vs {nOK:.4f})   [{time.time()-t0:.0f}s]")
# [3] where does the Theta* looseness come from: Perron weight vs p weight
w = 45.0; N = int(200*w); n = np.arange(N); wn = 1/(2*n+1.0); sn=(n+0.5)/w; p=np.sqrt(wn)
S,T = np.meshgrid(sn,sn,indexing='ij'); Rop = np.outer(p,p)*RG(S,T)
lam, vec = eigsh(Rop, k=1, which='LA'); vec = np.abs(vec[:,0])
print(f"\n[3] N=200w, w=45: ||Rop||={lam[0]:.5f}; p-weighted Schur sup={np.max((Rop@p)/p):.5f} (=Theta(s_n) max); Perron-weighted Schur sup={np.max((Rop@vec)/vec):.5f}; row-sum profile Theta(s) at s=0.5,1,1.5,2,3: {[round(float(np.sum(RG(s,sn)*wn)),4) for s in (0.5,1,1.5,2,3)]}")
print(f"    Perron vector shape: v_n*sqrt(2n+1) at s=0.1,0.5,1,1.5,3,6: {[round(float(vec[np.argmin(np.abs(sn-s))]/p[np.argmin(np.abs(sn-s))]/ (vec[0]/p[0])),3) for s in (0.1,0.5,1,1.5,3,6)]}")
# [6] vacuity: two-dim variational at w=45 with z = u 1_[300,2900) (main session's probe)
w=45.0; N=int(40*w)
H = opK(N, tab(N,Gf,w)); e,V = np.linalg.eigh(H); h=e[-1]; v=V[:,-1]*np.sign(V[0,-1])
u = 1/np.sqrt(2*np.arange(N)+1.0)
for (a,bb) in [(300,2900),(90,1800)]:
    bb=min(bb,N); z=np.zeros(N); z[a:bb]=u[a:bb]; z-= (v@z)*v; z/=np.linalg.norm(z)
    uv,uz = u@v, u@z; Hz=H@z
    Mm = np.array([[2*uv*uv-h, 2*uv*uz - v@Hz],[2*uv*uz - v@Hz, 2*uz*uz - z@Hz]])
    print(f"[6] w=45 N={N} z=u*1_[{a},{bb}): <z,Hz>={z@Hz:.4f} <u,z>={uz:.3f} <u,v>={uv:.3f}  h+lam_min(M) <= {h+np.linalg.eigvalsh(Mm)[0]:.4f}  (h={h:.4f}, true gap {h-e[-2]:.4f})")
print(f"[done {time.time()-t0:.0f}s]")
