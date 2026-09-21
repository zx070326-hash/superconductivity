"""[numerical-evidence] Independent check of universality-lift's M2^Delta=0 counterexample,
using MY builder (gate-0'd against the scalar O), plus my own diagnosis."""
import os
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
G2 = 1.318140212132; PI26 = np.pi**2/6

def sO(fv, N):
    n = np.arange(N); u = 1/np.sqrt(2*n+1.)
    T = fv[np.abs(n[:,None]-n[None,:])].copy(); np.fill_diagonal(T,0.0)
    M = (T + fv[n[:,None]+n[None,:]+1])*np.outer(u,u)
    M[n,n] -= 2.0*np.concatenate(([0.],np.cumsum(fv[1:N])))/(2*n+1.)
    return M

def build1(fZ, fD, N):
    """B=1 multiband builder: L1 driven by fZ, pairing blocks by fD."""
    n = np.arange(N); u = 1/np.sqrt(2*n+1.)
    T = fD[np.abs(n[:,None]-n[None,:])].copy(); np.fill_diagonal(T,0.0)
    M = (T + fD[n[:,None]+n[None,:]+1])*np.outer(u,u)
    M[n,n] -= 2.0*np.concatenate(([0.],np.cumsum(fZ[1:N])))/(2*n+1.)
    M[n,n] -= u*u*(fZ[0]-fD[0])
    return M

wph, lph = 0.1, 10.0
M2ph = lph*wph**2                    # = 0.1  (units of kappa)
wsf = 5.0; lsf = M2ph/wsf**2         # = 0.004
M2sf = lsf*wsf**2
print(f"M2^ph/k = {M2ph:.6f}   M2^sf/k = {M2sf:.6f}   M2^Delta/k = {M2ph-M2sf:.3e}")

for N in (100,200,400,800,1600):
    kk = np.arange(2*N+2, dtype=float)
    Lph = lph*wph**2/(kk**2+wph**2); Lsf = lsf*wsf**2/(kk**2+wsf**2)
    lam1 = np.linalg.eigvalsh(build1(Lph+Lsf, Lph-Lsf, N))[-1]
    if N in (100,1600) or N==400:
        print(f"  N={N:5d}  lam1(O) = {lam1:.6f}")
lam1 = lam1

N=800; kk=np.arange(2*N+2,dtype=float)
Lph=lph*wph**2/(kk**2+wph**2); Lsf=lsf*wsf**2/(kk**2+wsf**2); LD=Lph-Lsf
print(f"\nCEILING on rho(|M2^Delta|)  = G2*{abs(M2ph-M2sf):.3e} = {G2*abs(M2ph-M2sf):.3e}"
      f"   vs lam1 = {lam1:.6f}   -> VIOLATED, ratio = {lam1/(G2*abs(M2ph-M2sf)+1e-300):.3e}")
print(f"  my Sec 8.3 bound  [g2*M2ph + (pi^2/6)*M2sf] = {G2*M2ph+PI26*M2sf:.6f}"
      f"   -> ratio {lam1/(G2*M2ph+PI26*M2sf):.4f}   (holds)")
print(f"  univ-lift M2^G (c=3)      [g2*(M2ph+3*M2sf)] = {G2*(M2ph+3*M2sf):.6f}"
      f"   -> ratio {lam1/(G2*(M2ph+3*M2sf)):.4f}   (holds, looser)")

# --- is the hypothesis of my Sec 6.2 theorem satisfied?  |Lam^D| nonneg non-increasing?
d = np.diff(np.abs(LD[1:200]))
print(f"\nHYPOTHESIS CHECK of my Sec 6.2 thm:  Lam^Delta(k)>0 for all k>=1? "
      f"{bool(np.all(LD[1:200]>0))};  |Lam^Delta| non-increasing? {bool(np.all(d<=1e-18))}"
      f"   ==> premises HOLD, conclusion FALSE  ==> my Sec 6.2 theorem is FALSIFIED")

# --- MY diagnosis: Lam^Delta is positive & non-increasing but leaves the class P.
print("\nDIAGNOSIS: tail exponent of Lam^Delta (class P requires exactly k^-2):")
for k in (50,200,1000,5000):
    print(f"   k={k:5d}   k^2*Lam^D(k) = {k**2*(lph*wph**2/(k**2+wph**2)-lsf*wsf**2/(k**2+wsf**2)):.6e}"
          f"    k^4*Lam^D(k) = {k**4*(lph*wph**2/(k**2+wph**2)-lsf*wsf**2/(k**2+wsf**2)):.6e}")
print("   -> k^2*f -> 0 while k^4*f -> const  ==>  f ~ k^-4  ==>  f NOT in P")
print("      (same fact as the repo's 'gamma>2 kernels are outside Theorem E class':")
print("       a positive Einstein mixture has tail exactly k^-2, never faster.)")
print("   ==> Lemma M / GF bound lam1 by the moment of the POSITIVE representing measure;")
print("       the DIFFERENCE of two such moments is not a moment of anything nonneg.")

# --- does GF itself survive on this f?  (GF only needs f>=0 non-increasing)
def twoC(fv, Nth=1<<16):
    a=np.zeros(Nth); K=min(len(fv)-1,Nth//2-1); a[1:K+1]=fv[1:K+1]
    return -2.0*np.fft.rfft(a).real.min()
lminD = np.linalg.eigvalsh(sO(np.r_[0.,LD[1:]], 400))[0]
print(f"\nGF on f=Lam^Delta directly (f>=0 non-incr, so GF applies): "
      f"-lmin={-lminD:.6f} <= max(f(1),2C_f)={max(LD[1],twoC(np.r_[0.,LD[1:]])):.6f}  "
      f"{-lminD <= max(LD[1],twoC(np.r_[0.,LD[1:]]))}")
