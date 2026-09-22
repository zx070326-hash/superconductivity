import numpy as np

def check(N, c, u, c_off):
    chi0 = c*np.eye(N) + c_off*(np.ones((N,N))-np.eye(N))
    U = u*np.eye(N)
    rpa = chi0 @ np.linalg.inv(np.eye(N) - U@chi0)
    def D(M):
        return min(M[i,i]-sum(abs(M[i,j]) for j in range(N) if j!=i) for i in range(N))
    return D(chi0), D(rpa), D(rpa)/D(chi0)

for N in [2,3,4,5,6,7,8,9,10]:
    c=1.0; u=0.3; c_off=0.5*c/(N-1)
    d0,drpa,ratio = check(N,c,u,c_off)
    print(N, d0, drpa, ratio)
