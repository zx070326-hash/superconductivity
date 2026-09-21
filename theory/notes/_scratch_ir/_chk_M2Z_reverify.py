"""[numerical-evidence] W1 re-verification of the CORRECTION itself.
Sample DELIBERATELY at the degenerate set M2^Delta ~ 0 that both random scans missed."""
import os
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
G2=1.318140212132; PI26=np.pi**2/6
rng=np.random.default_rng(20260918)

def build(LZf,LDf,N):
    B=LZf.shape[0]; n=np.arange(N); u=1/np.sqrt(2*n+1.)
    ai=np.abs(n[:,None]-n[None,:]); si=n[:,None]+n[None,:]+1
    SZ=LZf.sum(axis=1); M=np.zeros((B*N,B*N))
    for i in range(B):
        for j in range(B):
            fv=LDf[i,j]; T=fv[ai].copy(); np.fill_diagonal(T,0.)
            blk=(T+fv[si])*np.outer(u,u)
            if i==j:
                blk=blk-np.diag(2.*np.concatenate(([0.],np.cumsum(SZ[i][1:N])))/(2*n+1.))\
                       -np.diag(u*u*(SZ[i][0]-fv[0]))
            else: blk=blk+np.diag(u*u*fv[0])
            M[i*N:(i+1)*N,j*N:(j+1)*N]=blk
    return M

def trial(B,N,r_target):
    """r_target = M2^sf/M2^ph ratio; r_target=1 => M2^Delta = 0 exactly."""
    K=2*N+2; kk=np.arange(K,dtype=float)
    wph=10**rng.uniform(-1.3,0.8,(B,B)); wph=(wph+wph.T)/2
    wsf=10**rng.uniform(-1.3,0.8,(B,B)); wsf=(wsf+wsf.T)/2
    cph=rng.uniform(0,1.2,(B,B)); cph=(cph+cph.T)/2
    M2ph=cph*wph**2
    M2sf=r_target*M2ph                      # force near-cancellation
    csf=M2sf/wsf**2
    Lph=np.empty((B,B,K)); Lsf=np.empty((B,B,K))
    for i in range(B):
        for j in range(B):
            Lph[i,j]=cph[i,j]*wph[i,j]**2/(kk**2+wph[i,j]**2)
            Lsf[i,j]=csf[i,j]*wsf[i,j]**2/(kk**2+wsf[i,j]**2)
    lam1=np.linalg.eigvalsh(build(Lph+Lsf,Lph-Lsf,N))[-1]
    if lam1<=0: return None
    M2D=np.abs(M2ph-M2sf)
    M2Zp=np.diag(np.diag(M2ph+PI26/G2*M2sf))+ (M2ph+M2sf)*PI26/G2*(1-np.eye(B))
    rho=lambda A: np.abs(np.linalg.eigvals(A)).max()
    return lam1/(G2*rho(M2D)+1e-300), lam1/(G2*rho(M2Zp))

for r in (1.0, 0.999, 0.99, 0.9, 0.5):
    for B in (2,3):
        res=[t for _ in range(150) if (t:=trial(B,90,r))]
        a=np.array(res)
        print(f"r=M2sf/M2ph={r:<6}  B={B}   rho(|M2^D|) ver: {int((a[:,0]>1).sum())}/{len(a)} viol, worst {a[:,0].max():.4g}"
              f"   |  my Sec8.3 rho(M2''_Z): {int((a[:,1]>1).sum())}/{len(a)} viol, worst {a[:,1].max():.4f}")

# N-stability of the sharpest surviving form at the hardest setting
print()
for N in (90,180,360):
    rng2=np.random.default_rng(7); 
    globals()['rng']=rng2
    a=np.array([t for _ in range(60) if (t:=trial(3,N,1.0))])
    print(f"  N={N:4d}  Sec8.3 worst ratio = {a[:,1].max():.6f}   ({int((a[:,1]>1).sum())}/{len(a)} viol)")
