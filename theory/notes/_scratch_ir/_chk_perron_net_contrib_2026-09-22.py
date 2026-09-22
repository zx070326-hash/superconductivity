import numpy as np

for N in [2,3,4,5,6]:
    v = np.ones(N)/np.sqrt(N)
    diag = v[0]**2
    offdiag_sum = sum(v[0]*v[b] for b in range(1,N))
    net = diag - offdiag_sum
    print(N, diag, offdiag_sum, net)
