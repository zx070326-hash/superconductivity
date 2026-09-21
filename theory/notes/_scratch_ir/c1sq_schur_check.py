# 独立验证：朴素 Schur 界 vs 真实 ||C||  (float64, numerical-evidence)
# C_{nm} = u_n u_m ( f(|n-m|) + f(n+m+1) ),  n<=N<m
import numpy as np

def f(k, w):  return w*w/(w*w + k*k)
def u(n):     return 1.0/np.sqrt(2.0*n+1.0)

def true_Cnorm(w, N, M):
    n = np.arange(0, N+1)
    m = np.arange(N+1, M+1)
    K = f(np.abs(n[:,None]-m[None,:]), w) + f(n[:,None]+m[None,:]+1.0, w)
    C = u(n)[:,None]*u(m)[None,:]*K
    return np.linalg.norm(C, 2), C

def schur_bound_finite(C):
    R = C.sum(axis=1).max()     # sup row sum
    Kc = C.sum(axis=0).max()    # sup col sum
    return np.sqrt(R*Kc), R, Kc

def old_bound(w, N):
    fh0 = np.pi*w/np.tanh(np.pi*w)
    fhpi = np.pi*w/np.sinh(np.pi*w)
    return u(N+1)*(fh0-fhpi)/2.0

def my_estimate(w, N):
    # u_n u_{N+1} * 2*sum_{j>=1} f(j),  worst n=N  ->  ~ pi*w/(2N+3) for large w
    S = 2.0*sum(f(j,w) for j in range(1, 200000))
    return u(N)*u(N+1)*S

for w in (1.0, 5.0, 30.0, 100.0):
    for N in (500, 1000):
        M = max(20*N, N+40000)
        tn, C = true_Cnorm(w, N, M)
        sb, R, Kc = schur_bound_finite(C)
        print(f"w={w:6.1f} N={N:5d} M={M:6d} | true={tn:.6f}  schur(finite)={sb:.6f}"
              f"  ratio={sb/tn:6.3f} | old={old_bound(w,N):.6f} ({old_bound(w,N)/tn:7.1f}x)"
              f" | my_est={my_estimate(w,N):.6f}  R={R:.6f} K={Kc:.6f}")
