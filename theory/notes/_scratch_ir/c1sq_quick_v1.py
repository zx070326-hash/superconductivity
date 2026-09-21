import numpy as np

def build_float(varpi, N):
    dim = N+1
    def f(k): return varpi**2/(varpi**2+k**2)
    def u(n): return 1/np.sqrt(2*n+1)
    a = np.zeros(dim)
    running=0.0
    for n in range(1,dim):
        running += f(n)
        a[n] = 2*running/(2*n+1)
    H = np.zeros((dim,dim))
    for n in range(dim):
        un=u(n)
        for m in range(dim):
            um=u(m)
            if n==m:
                H[n,m] = -a[n]+f(2*n+1)*un*un
            else:
                H[n,m] = f(abs(n-m))*un*um + f(n+m+1)*un*um
    return H

for varpi in [0.5,1.0,2.0,3.0,5.0,10.0,30.0]:
    H = build_float(varpi,500)
    evals,evecs = np.linalg.eigh(H)
    v1=evecs[:,-1]
    print(varpi, "v1[0]=",v1[0], "c1sq_approx(N=500)=", v1[0]**2)
