import numpy as np
np.set_printoptions(linewidth=150, suppress=True)
N = np.array([0.556, 0.646, 0.616, 0.370])   # Table 1, states/spin/eV
rows = {"sf":(0.00,1.78,0.66,0.45,0.52), "sf,?":(2.10,1.15,0.80,0.45,0.30),
        "sf,ph(1)":(0.86,1.06,0.79,0.42,0.30), "sf,ph(2)":(0.90,1.15,0.80,0.45,0.30)}
def build(l11,l13,l23,l14,l24):
    L=np.zeros((4,4)); L[0,0]=l11
    L[0,2],L[1,2],L[0,3],L[1,3]=l13,l23,l14,l24
    for (i,j) in [(0,2),(1,2),(0,3),(1,3)]: L[j,i]=L[i,j]*N[i]/N[j]
    return L
print("nu_ij=N_i/N_j: nu13=%.4f nu14=%.4f nu23=%.4f nu24=%.4f  (paper: .9019/1.5010/1.0483/1.7447)"
      % (N[0]/N[2],N[0]/N[3],N[1]/N[2],N[1]/N[3]))
for nm,p in rows.items():
    L=build(*p); K=L/N[None,:]; Ks=(K+K.T)/2; e=np.linalg.eigvalsh(Ks)
    viol=[(i+1,j+1,round(L[i,i]*L[j,j]-L[i,j]*L[j,i],4)) for i in range(4) for j in range(i+1,4)
          if L[i,i]*L[j,j]-L[i,j]*L[j,i] < -1e-12]
    print(f"\n{nm}: K symmetric? max|K-K^T|=%.2e | ALL lambda_ij >= 0 ? %s"
          % (np.abs(K-K.T).max(), bool((L>=-1e-15).all())))
    print("   eig(K_sym)=",e.round(4)," PSD?",e.min()>=-1e-9,
          " lam_max=%.4f rho=%.4f equal?%s"%(e.max(),abs(e).max(),abs(e.max()-abs(e).max())<1e-9))
    print("   2x2 principal-minor violations:",viol)
