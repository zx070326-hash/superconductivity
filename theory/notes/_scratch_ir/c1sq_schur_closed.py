# 我给 agent 的粗估「最坏在 n=N」是否成立？ + 闭式行/列和上界是否真包住
import numpy as np
def f(k,w): return w*w/(w*w+k*k)
def u(n):   return 1.0/np.sqrt(2.0*n+1.0)

def tail_sum(J, w, JMAX=400000):     # sum_{j>=J} f(j)
    j = np.arange(J, JMAX); return f(j,w).sum()

for w in (1.,5.,30.,100.):
    for N in (500,1000):
        M = max(20*N, N+40000)
        n = np.arange(0,N+1); m = np.arange(N+1,M+1)
        K = f(np.abs(n[:,None]-m[None,:]),w)+f(n[:,None]+m[None,:]+1.,w)
        C = u(n)[:,None]*u(m)[None,:]*K
        R_exact = C.sum(axis=1); argn = int(R_exact.argmax())
        # 闭式上界: u_m <= u_{N+1};  sum_m f(m-n) = sum_{j>=N+1-n} f(j);  sum_m f(n+m+1) = sum_{j>=N+n+2} f(j)
        R_cf = np.array([u(nn)*u(N+1)*(tail_sum(N+1-nn,w)+tail_sum(N+nn+2,w)) for nn in n])
        ok = bool((R_cf >= R_exact - 1e-12).all())
        print(f"w={w:6.1f} N={N:5d} | argmax_n(R)={argn:4d} (n=N? {argn==N})"
              f" R_max={R_exact.max():.6f}  closed-form_max={R_cf.max():.6f}"
              f"  encloses_all_rows={ok}")
