"""Route: keep v. r_perp = (1/x)||(I-P1) O[f^2] v||.
Use an explicit TRIAL vector t (normalized) for v and control the error by the
gap via Davis-Kahan: ||v - t|| <= sqrt(2)*||(A-h)t||/g  (sin-theta).
Then r_perp <= (1/x)[ ||(I-P1)O[f^2]t|| + ||O[f^2]||*||v-t|| ].
Candidate trial vectors (route R in P4_main_theorems already uses one for varpi>=8):
  t1: u_n/||u||_N truncated  (u_n = 1/sqrt(2n+1))  -- the tau-direction
  t2: tau_n = p_n G(s_n) from (M1), p_n = u_n/||u||, G(s)=s^2/(1+s^2), s_n=(n+1/2)/varpi
Report overlap |<t,v>| and the resulting bound vs truth.
"""
import numpy as np
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
print(f"{'x':>7} {'varpi':>6} {'|<t_tau,v>|':>11} {'|<t_u,v>|':>10} {'r_perp':>9} {'bnd_tau':>9} {'ratio':>6}")
for x in [0.04,0.25,1.0,4.0,9.0,25.0,100.0]:
    vp=np.sqrt(x); N=max(200,int(40*vp)+200); K=2*N+2
    k=np.arange(1,K+1,dtype=float); f=x/(k**2+x)
    A=O_of_kernel(N,f); Af2=O_of_kernel(N,f**2); Ax=O_of_kernel(N,k**2/(k**2+x)**2)
    w,V=np.linalg.eigh(A); o=np.argsort(w)[::-1]; w=w[o]; V=V[:,o]
    h=w[0]; v=V[:,0]; g=h-w[1]
    if v.sum()<0: v=-v
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0)
    t_u=u/np.linalg.norm(u)
    s=(n+0.5)/vp; G=s**2/(1+s**2); tau=t_u*G; t_tau=tau/np.linalg.norm(tau)
    r=Ax@v; rp=np.linalg.norm(r-(v@r)*v)
    # bound with trial tau + Davis-Kahan
    res=np.linalg.norm((A-h*np.eye(N))@t_tau)
    dv=np.sqrt(2)*res/g
    nb=np.abs(Af2).sum(1).max()
    bnd=(np.linalg.norm(Af2@t_tau-(t_tau@(Af2@t_tau))*t_tau)+nb*dv)/x
    print(f"{x:7.2f} {vp:6.2f} {abs(t_tau@v):11.6f} {abs(t_u@v):10.6f} {rp:9.5f} {bnd:9.5f} {bnd/rp:6.2f}")
