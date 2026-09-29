import sys
import numpy as np
sys.path.insert(0, "D:/physics/superconductivity/.claude/worktrees/agent-a8932ec4c14045b28/theory/derivations")
from p4_concavity_q_route import O_ker, kernels, top_eigen

# Analyticize the Bhatia-Davis inputs bmax,bmin,q for Var_v(B)<=BD.
# Gershgorin (row abs-sum) gives analytic bmax_hi, bmin_lo for O[f^2].
# q_lo: q = h - x<v,Dv> >= h_lo - x*lam1(O[f']); lam1(O[f']) via Gershgorin.
#   BD_analytic = (bmax_hi - q_lo)(q_lo - bmin_lo) if q_lo >= midpoint else BDmid.
THETA=0.6349
N=400
def gersh(M):
    d=np.diag(M).copy(); r=np.sum(np.abs(M),axis=1)-np.abs(d)
    return float(np.max(d+r)), float(np.min(d-r))   # lam1_hi, lammin_lo
TRIAL={"e0":np.array([1.0]),"inv_n":1.0/(np.arange(40)+1.0),
       "inv_sqrt":1.0/np.sqrt(2*np.arange(40)+1.0)}
def rq(w,x,KX):
    M=len(w); k=np.arange(KX,dtype=float); fk=x/(k*k+x); fk[0]=1.0
    n=np.arange(M); u=1.0/np.sqrt(2.0*n+1.0)
    Ai=np.abs(n[:,None]-n[None,:]); S=n[:,None]+n[None,:]+1
    off=fk[np.minimum(Ai,KX-1)].copy(); np.fill_diagonal(off,0.0)
    Mm=(off+fk[np.minimum(S,KX-1)])*(np.outer(u,u))
    cum=np.concatenate(([0.0],np.cumsum(2.0*fk[1:M]))); Mm[np.diag_indices(M)]-=cum/(2.0*n+1.0)
    return float(w@(Mm@w))/float(w@w)
def gdom(x,Nloc,KX=60000):
    n=np.arange(Nloc); u2=1.0/(2.0*n+1.0); u=np.sqrt(u2)
    k=np.arange(KX,dtype=float); fk=x/(k*k+x); fk[0]=0.0
    csum=np.cumsum(fk); total=csum[-1]+x/(KX-0.5); T=total-csum[:Nloc]
    S=n[:,None]+n[None,:]+1; G=(u[:,None]*u[None,:])*fk[np.minimum(S,KX-1)]
    G[np.diag_indices(Nloc)]=u2*T+u2*fk[np.minimum(2*n+1,KX-1)]
    return float(np.sum(G*G)), float(u2[0]*T[0]+u2[0]*fk[1])

hdr=("varpi  bmaxT bmaxH  bminT bminL  qT    qlo   mid   Ok? BDt   BDan  an/T  ratio")
print(hdr); print("-"*len(hdr))
worst=0.0; wv=0
for varpi in [0.2,0.5,0.8,1.0,1.2,1.35,1.5,1.6,1.8,2.0,2.5,3.0,5.0,8.0,20.0,40.0]:
    x=varpi**2
    ker=kernels(x,KX=8*N); A=O_ker(N,ker["f"]); B=O_ker(N,ker["f2"]); Dm=O_ker(N,ker["df"])
    h,v,l2=top_eigen(A); gt=h-l2
    E=v@(O_ker(N,ker["g"])@v)
    q=float(v@(B@v)); Bv=B@v; Var=float(Bv@Bv-q*q)
    eB=np.linalg.eigvalsh(B); bmaxT=float(eB[-1]); bminT=float(eB[0])
    BDt=(bmaxT-q)*(q-bminT)
    bmaxH,bminL=gersh(B)
    dmaxH,_=gersh(Dm)
    hlo=max(rq(w,x,8*N) for w in TRIAL.values())
    qlo=hlo - x*dmaxH
    mid=(bmaxH+bminL)/2.0
    ok = qlo>=mid
    if ok:
        BDan=(bmaxH-qlo)*(qlo-bminL)
    else:
        BDan=((bmaxH-bminL)/2.0)**2
    trG2,g00=gdom(x,N); l2hi=min(np.sqrt(max(trG2-g00*g00,0.0)),THETA)
    glo=hlo-l2hi
    ratio=(2.0/x)*BDan/glo/E if glo>0 else np.inf
    if np.isfinite(ratio) and ratio>worst: worst=ratio; wv=varpi
    print(f"{varpi:5.2f} {bmaxT:5.3f} {bmaxH:5.3f}  {bminT:+5.3f} {bminL:+5.3f} "
          f"{q:5.3f} {qlo:5.3f} {mid:5.3f} {'Y' if ok else 'N'}  {BDt:.4f} {BDan:.4f} "
          f"{BDan/Var:5.2f} {ratio:6.3f}")
print()
print(f"worst full-analytic ratio = {worst:.3f} at varpi={wv}")
