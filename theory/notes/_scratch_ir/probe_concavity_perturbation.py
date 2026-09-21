"""Probe: is H(x)=h(sqrt x) concavity readable from 2nd-order perturbation theory?

H''(x) = <v, O[f_xx] v>  +  2 * sum_{j>=2} <v, O[f_x] w_j>^2 / (lam_1 - lam_j)
where f(k;x) = x/(k^2+x)  (so H(x)=lam_1(O[f(.;x)]), x=varpi^2),
f_x = k^2/(k^2+x)^2,  f_xx = -2k^2/(k^2+x)^3.
Term 1 <= 0 is NOT automatic; term 2 >= 0 always (lam_1 top). So concavity needs
term1 to beat term2. Print both, plus the two structural quantities that could
make it analytic: A := <v,O[f_xx]v>  and  B := 2*sum_j ... .
"""
import numpy as np, sys
sys.path.insert(0, r"D:/physics/superconductivity/theory/derivations")
from p4_interchange_uniformity import H_mat

def O_of_kernel(N, fvals):
    """O[f] for f given as array fvals[k-1], k=1..2N."""
    n = np.arange(N); u = 1.0/np.sqrt(2*n+1.0)
    A = np.zeros((N,N))
    # L2 + L3
    for i in range(N):
        for j in range(N):
            s = 0.0
            if i != j: s += fvals[abs(i-j)-1]
            s += fvals[i+j]           # k = i+j+1
            A[i,j] = s*u[i]*u[j]
    # -L1 diagonal
    for i in range(N):
        A[i,i] -= 2.0*u[i]**2*fvals[:i].sum() if i>0 else 0.0
    return A

def parts(N, x):
    K = 2*N+2
    k = np.arange(1, K+1, dtype=float)
    f   = x/(k**2+x)
    fx  = k**2/(k**2+x)**2
    fxx = -2.0*k**2/(k**2+x)**3
    A = O_of_kernel(N, f)
    w, V = np.linalg.eigh(A)
    i1 = np.argmax(w); l1 = w[i1]; v = V[:,i1]
    Ax  = O_of_kernel(N, fx)
    Axx = O_of_kernel(N, fxx)
    term1 = v @ (Axx @ v)
    num = (V.T @ (Ax @ v))**2
    denom = l1 - w
    mask = np.ones(N, bool); mask[i1] = False
    term2 = 2.0*np.sum(num[mask]/denom[mask])
    return l1, term1, term2, term1+term2

print(f"{'x':>9} {'varpi':>7} {'H':>10} {'term1(<=0?)':>13} {'term2(>=0)':>11} {'H\"':>11} {'ratio |t1|/t2':>13}")
for x in [0.01,0.04,0.1,0.25,1.0,4.0,9.0,25.0,100.0,400.0,1600.0]:
    N = max(120, int(30*np.sqrt(x))+120)
    l1,t1,t2,h2 = parts(N, x)
    print(f"{x:9.3f} {np.sqrt(x):7.3f} {l1:10.5f} {t1:13.6f} {t2:11.6f} {h2:11.6f} {abs(t1)/t2:13.4f}")
