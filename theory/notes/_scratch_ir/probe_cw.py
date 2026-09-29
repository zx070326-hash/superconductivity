import sys
import numpy as np
sys.path.insert(0, "D:/physics/superconductivity/.claude/worktrees/agent-a8932ec4c14045b28/theory/derivations")
from p4_concavity_q_route import O_ker, kernels, top_eigen

# Collatz-Wielandt UPPER bound for top eigenvalue of an ENTRYWISE-POSITIVE
# dominating operator Gdom[phi] >= O[phi]:  lam1 <= max_n (Gdom[phi] w)_n / w_n.
# Test tightness for beta_max (phi=f^2) and for lam1(O[f']) used in q_lo.
THETA=0.6349; N=400
def Gdom_ker(phi, Nloc, KX):
    # phi: array phi[0..KX-1], phi[k]=phi(k); build dominating operator (drop Dirichlet)
    n=np.arange(Nloc); u2=1.0/(2.0*n+1.0); u=np.sqrt(u2)
    csum=np.cumsum(np.concatenate(([0.0],phi[1:KX])))  # csum[j]=sum_{k=1}^{j} phi(k)
    tot=csum[-1]
    T=tot-csum[:Nloc]                     # T(n)=sum_{k>n} phi(k)
    S=n[:,None]+n[None,:]+1
    G=(u[:,None]*u[None,:])*phi[np.minimum(S,KX-1)]
    G[np.diag_indices(Nloc)]=u2*T+u2*phi[np.minimum(2*n+1,KX-1)]
    return G
def CW(G, w):
    Gw=G@w
    return float(np.max(Gw/w))
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

# CW test vectors (positive), moderate length
wcw = 1.0/(np.arange(N)+1.0)
wcw2 = 1.0/np.sqrt(2*np.arange(N)+1.0)
print("varpi  bmaxT  bmaxCW  cw/T   l1DfT  l1DfCW  h'true  x*l1Df  qtrue  qlo")
for varpi in [0.5,1.0,1.5,2.0,3.0,5.0,8.0,20.0]:
    x=varpi**2; KX=8*N
    ker=kernels(x,KX=KX); A=O_ker(N,ker["f"]); B=O_ker(N,ker["f2"]); Dm=O_ker(N,ker["df"])
    h,v,l2=top_eigen(A)
    q=float(v@(B@v)); hp=float(v@(Dm@v))
    bmaxT=float(np.linalg.eigvalsh(B)[-1])
    l1DfT=float(np.linalg.eigvalsh(Dm)[-1])
    # dominating operators for f^2 and f'
    k=np.arange(KX,dtype=float)
    f2=(x/(k*k+x))**2; f2[0]=0.0
    df=k*k/((k*k+x)**2); df[0]=0.0
    Gf2=Gdom_ker(f2,N,KX); Gdf=Gdom_ker(df,N,KX)
    bmaxCW=min(CW(Gf2,wcw),CW(Gf2,wcw2))
    l1DfCW=min(CW(Gdf,wcw),CW(Gdf,wcw2))
    hlo=max(rq(w,x,KX) for w in TRIAL.values())
    qlo=hlo-x*l1DfCW
    print(f"{varpi:5.2f} {bmaxT:6.4f} {bmaxCW:6.4f} {bmaxCW/bmaxT:5.3f}  "
          f"{l1DfT:6.4f} {l1DfCW:6.4f}  {hp:6.4f} {x*l1DfCW:6.4f} {q:6.4f} {qlo:+6.4f}")
