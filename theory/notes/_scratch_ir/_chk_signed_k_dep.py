"""[numerical-evidence] The arrow my Sec 6 never tested: Lambda^Delta with a
k-DEPENDENT sign (omega_ph != omega_sf).  My step-2 (|O[-f]|=|O[f]|) assumes the
sign is k-independent; R2's 0/1500 used ONE varpi for both channels."""
import os
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
G2 = 1.318140212132
rng = np.random.default_rng(20260918)

def build(LZf, LDf, N):
    """LZf/LDf: (B,B,K) kernel arrays, K >= 2N+2.  General k-dependent kernels."""
    B = LZf.shape[0]; n = np.arange(N); u = 1/np.sqrt(2*n+1.)
    ai = np.abs(n[:,None]-n[None,:]); si = n[:,None]+n[None,:]+1
    SZ = LZf.sum(axis=1)                                  # SZ[i,k]
    M = np.zeros((B*N,B*N))
    for i in range(B):
        for j in range(B):
            fv = LDf[i,j]
            T = fv[ai].copy(); np.fill_diagonal(T,0.0)
            blk = (T + fv[si])*np.outer(u,u)
            if i==j:
                d = 2.0*np.concatenate(([0.],np.cumsum(SZ[i][1:N])))/(2*n+1.)
                blk = blk - np.diag(d) - np.diag(u*u*(SZ[i][0]-fv[0]))
            else:
                blk = blk + np.diag(u*u*fv[0])
            M[i*N:(i+1)*N, j*N:(j+1)*N] = blk
    return M

def M2mat(Kf, kk):
    """M2_ij/kappa for f(k)=sum c*w^2/(k^2+w^2): the moment is c*w^2 = f(k)k^2 as k->inf.
    Compute exactly from the construction instead (passed in separately)."""
    raise NotImplementedError

K = 2*200+2; kk = np.arange(K, dtype=float)
def eins(c,w): 
    v = c*w*w/(kk*kk+w*w); return v, c*w*w      # kernel array, and M2/kappa

B=2; N=90
nvA=nvZ=0; wA=wZ=-np.inf; nmono=0; ex=None
for trial in range(800):
    # phonon and spin-fluct blocks, DIFFERENT frequencies -> k-dependent sign
    cph = rng.uniform(0,1.2,(B,B)); cph=(cph+cph.T)/2
    csf = rng.uniform(0,1.2,(B,B)); csf=(csf+csf.T)/2
    wph = 10**rng.uniform(-1,0.6); wsf = 10**rng.uniform(-1,0.6)
    LZf = np.zeros((B,B,K)); LDf = np.zeros((B,B,K))
    M2D = np.zeros((B,B)); M2Z = np.zeros((B,B))
    for i in range(B):
        for j in range(B):
            a,ma = eins(cph[i,j],wph); b,mb = eins(csf[i,j],wsf)
            LZf[i,j]=a+b; LDf[i,j]=a-b; M2D[i,j]=ma-mb; M2Z[i,j]=ma+mb
    # is |Lambda^Delta(k)| non-increasing on k>=1 ?
    mono = all(np.all(np.diff(np.abs(LDf[i,j,1:60]))<=1e-15) for i in range(B) for j in range(B))
    nmono += (not mono)
    lam1 = np.linalg.eigvalsh(build(LZf,LDf,N))[-1]
    rA = lam1/(G2*np.abs(np.linalg.eigvals(np.abs(M2D))).max())
    rZ = lam1/(G2*np.abs(np.linalg.eigvals(M2Z)).max())
    nvA += rA>1.0; nvZ += rZ>1.0
    if rA>wA: wA,ex = rA,(cph,csf,wph,wsf,mono)
    wZ = max(wZ,rZ)
print(f"k-dependent sign, {800} trials  (|Lam^D| non-monotone in {nmono}/800)")
print(f"  ceiling with rho(|M2^Delta|) : {nvA}/800 VIOLATE, worst {wA:.4f}")
print(f"  ceiling with rho(M2^Z)       : {nvZ}/800 violate, worst {wZ:.4f}")
if nvA:
    cph,csf,wph,wsf,mono = ex
    print(f"  worst at w_ph={wph:.4f} w_sf={wsf:.4f}  |Lam^D| monotone? {mono}")
    for NN in (90,180,360):
        LZf=np.zeros((B,B,2*NN+2)); LDf=np.zeros((B,B,2*NN+2))
        kk2=np.arange(2*NN+2,dtype=float); M2D=np.zeros((B,B))
        for i in range(B):
            for j in range(B):
                a=cph[i,j]*wph**2/(kk2**2+wph**2); b=csf[i,j]*wsf**2/(kk2**2+wsf**2)
                LZf[i,j]=a+b; LDf[i,j]=a-b; M2D[i,j]=cph[i,j]*wph**2-csf[i,j]*wsf**2
        l=np.linalg.eigvalsh(build(LZf,LDf,NN))[-1]
        print(f"    N={NN:4d} -> ratio {l/(G2*np.abs(np.linalg.eigvals(np.abs(M2D))).max()):.6f}")
