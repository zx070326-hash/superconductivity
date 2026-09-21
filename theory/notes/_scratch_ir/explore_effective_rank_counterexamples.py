"""
Counterexample / generalization search for effective-rank universality.
Builds the generic kernel operator O[phi] (same formula family as
P4_concentration_bound.md Sec.0) for various toy phi(k), and checks
whether the spectral gap ratio (lam1-lam2)/lam1 stays close to 1,
or collapses (candidate counterexample to "effective rank always ~1").

numerical-evidence only, float64, not committed.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

def O_kernel(N, phi_func):
    """Generic operator (O[phi])_{nm} = u_n u_m [phi(|n-m|) 1_{n!=m} + phi(n+m+1)]
       - delta_nm * (2/(2n+1)) * sum_{k=1}^n phi(k)"""
    n = np.arange(N)
    u = 1.0/np.sqrt(2*n+1)
    kmax = 2*N+2
    phivals = np.array([phi_func(k) for k in range(0, kmax+1)])
    nn, mm = np.meshgrid(n, n, indexing='ij')
    absdiff = np.abs(nn-mm)
    off = np.where(nn != mm, phivals[absdiff], 0.0)
    summ = phivals[nn+mm+1]
    diagcorr = np.zeros(N)
    csum = np.cumsum(phivals[1:N+1])
    for i in range(1, N):
        diagcorr[i] = (2.0/(2*i+1))*csum[i-1]
    O = np.outer(u,u)*(off+summ) - np.diag(diagcorr)
    return O

def gap_stats(H):
    evals = np.sort(np.linalg.eigvalsh(H))[::-1]
    lam1, lam2, lam3 = evals[0], evals[1], evals[2]
    gap = lam1-lam2
    return lam1, lam2, lam3, gap, (gap/lam1 if lam1 != 0 else np.nan)

if __name__ == "__main__":
    N = 500
    xs = [1, 5, 20, 60]

    print("=== baseline Lorentz phi(k)=x/(k^2+x) (should match H(varpi), varpi=sqrt(x)) ===")
    for x in xs:
        phi = lambda k, x=x: x/(k*k+x)
        O = O_kernel(N, phi)
        lam1,lam2,lam3,gap,ratio = gap_stats(O)
        print(f"x={x:6.1f} lam1={lam1:8.4f} lam2={lam2:8.4f} lam3={lam3:8.4f} gap/lam1={ratio:8.5f}")

    print()
    print("=== oscillatory kernel phi(k)=cos(k)*x/(k^2+x)  [sign-alternating, long range] ===")
    for x in xs:
        phi = lambda k, x=x: np.cos(k)*x/(k*k+x)
        O = O_kernel(N, phi)
        lam1,lam2,lam3,gap,ratio = gap_stats(O)
        print(f"x={x:6.1f} lam1={lam1:8.4f} lam2={lam2:8.4f} lam3={lam3:8.4f} gap/lam1={ratio:8.5f}")

    print()
    print("=== heavy-tail kernel phi(k)=x/(k^2+x)^0.6  [slower decay, still even] ===")
    for x in xs:
        phi = lambda k, x=x: x/((k*k+x)**0.6)
        O = O_kernel(N, phi)
        lam1,lam2,lam3,gap,ratio = gap_stats(O)
        print(f"x={x:6.1f} lam1={lam1:8.4f} lam2={lam2:8.4f} lam3={lam3:8.4f} gap/lam1={ratio:8.5f}")

    print()
    print("=== power-law tail kernel phi(k)=1/(1+k)^1.2 [marginal, long range, no x scale] ===")
    for p in [1.05, 1.2, 1.5, 2.0]:
        phi = lambda k, p=p: 1.0/(1.0+k)**p
        O = O_kernel(N, phi)
        lam1,lam2,lam3,gap,ratio = gap_stats(O)
        print(f"p={p:5.2f} lam1={lam1:8.4f} lam2={lam2:8.4f} lam3={lam3:8.4f} gap/lam1={ratio:8.5f}")

    print()
    print("=== bimodal kernel: sum of two Lorentzians at different scales (two-peak spectrum) ===")
    for x1, x2, wgt in [(1,1,0.5), (0.5, 50, 0.5), (1, 100, 0.5), (0.5, 200, 0.9), (0.5,200,0.1)]:
        phi = lambda k, x1=x1, x2=x2, wgt=wgt: wgt*x1/(k*k+x1) + (1-wgt)*x2/(k*k+x2)
        O = O_kernel(N, phi)
        lam1,lam2,lam3,gap,ratio = gap_stats(O)
        print(f"x1={x1:6.1f} x2={x2:6.1f} w={wgt:4.2f} lam1={lam1:8.4f} lam2={lam2:8.4f} lam3={lam3:8.4f} gap/lam1={ratio:8.5f}")

    print()
    print("=== rank-boost test: TWO widely-separated Lorentz peaks with COMPARABLE weight & similar x (near-degenerate design attempt) ===")
    for x1, x2, wgt in [(3,3.01,0.5), (3, 3.1, 0.5), (3,4,0.5), (3,10,0.5)]:
        phi = lambda k, x1=x1, x2=x2, wgt=wgt: wgt*x1/(k*k+x1) + (1-wgt)*x2/(k*k+x2)
        O = O_kernel(N, phi)
        lam1,lam2,lam3,gap,ratio = gap_stats(O)
        print(f"x1={x1:6.2f} x2={x2:6.2f} w={wgt:4.2f} lam1={lam1:8.4f} lam2={lam2:8.4f} lam3={lam3:8.4f} gap/lam1={ratio:8.5f}")

if __name__ == "__main__":
    print()
    print("=== Mexican-hat kernel: positive core minus broader negative shoulder (still even, partly negative) ===")
    N=500
    for a,b,c in [(1.0,0.3,4.0),(1.0,0.5,4.0),(1.0,0.8,4.0),(1.0,0.95,4.0)]:
        phi = lambda k,a=a,b=b,c=c: a/(1+k*k) - b*a/(1+(k/c)**2)
        O = O_kernel(N, phi)
        lam1,lam2,lam3,gap,ratio = gap_stats(O)
        negfrac = phi(0) 
        print(f"core-a={a} shoulder-b={b} c={c}  lam1={lam1:8.4f} lam2={lam2:8.4f} gap/lam1={ratio:8.5f}")

    print()
    print("=== extreme far-separated bimodal, equal weight ===")
    for x1,x2 in [(0.3, 1000.0),(0.1,5000.0)]:
        phi = lambda k,x1=x1,x2=x2: 0.5*x1/(k*k+x1) + 0.5*x2/(k*k+x2)
        O = O_kernel(N, phi)
        lam1,lam2,lam3,gap,ratio = gap_stats(O)
        print(f"x1={x1} x2={x2}  lam1={lam1:8.4f} lam2={lam2:8.4f} gap/lam1={ratio:8.5f}")

    print()
    print("=== mild oscillation: phi(k) = (1+0.3cos(0.5k)) * x/(k^2+x)  [perturbed-positive, not sign changing] ===")
    for x in [1,5,20,60]:
        phi = lambda k,x=x: (1+0.3*np.cos(0.5*k))*x/(k*k+x)
        O = O_kernel(N, phi)
        lam1,lam2,lam3,gap,ratio = gap_stats(O)
        print(f"x={x}  lam1={lam1:8.4f} lam2={lam2:8.4f} gap/lam1={ratio:8.5f}  (phi always positive: {np.all(np.array([phi(k) for k in range(200)])>0)})")
