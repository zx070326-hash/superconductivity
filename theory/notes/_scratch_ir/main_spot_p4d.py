import numpy as np, sys
sys.path.insert(0,'.')
from p4_interchange_uniformity import H_mat
from p4d_coulomb_certify import mu2_pencil_cell, ROUND_PAD
def Kc(N,Nc):
    u=np.zeros(N); k=min(Nc,N); u[:k]=1/np.sqrt(2*np.arange(k)+1); return 2*np.outer(u,u)
N=400
print("== (iii) kappa: pencil vs phonon, and FD sigma^2 coefficient of the pencil ==")
for w,mt,Nc in ((1.0,0.1,10),(0.7,0.1,3),(1.0,0.0,0)):
    H=H_mat(N,w); H2=H_mat(N,w,2); M=H-mt*Kc(N,Nc)
    eM,VM=np.linalg.eigh(M); eH,VH=np.linalg.eigh(H)
    kM=w*w*(VM[:,-1]@H2@VM[:,-1])/eM[-1]; kH=w*w*(VH[:,-1]@H2@VH[:,-1])/eH[-1]
    def k_of(s):  # 2-pk symmetric, p=1/2 each
        A=0.5*(H_mat(N,w*(1+s))+H_mat(N,w*(1-s)))-mt*Kc(N,Nc); return np.linalg.eigvalsh(A)[-1]
    k0=k_of(0.0)
    def kap(s): return 2*(k_of(s)-k0)/(k0*s*s)
    kFD=(4*kap(0.01)-kap(0.02))/3
    print(f" w={w} mt={mt} Nc={Nc}: kappa_M={kM:.5f} kappa_H={kH:.5f} ratio={kM/kH:.3f}  kappa_FD={kFD:.5f} FD/M={kFD/kM:.4f}")
print("== (iv) covering lemma: lam1(M(mt)) >= lam1(M(mt0)) for mt<=mt0, every Nc ==")
viol=0
for w in (0.6,1.0,2.0):
    for Nc in range(0,60):
        l=[np.linalg.eigvalsh(H_mat(300,w)-mt*Kc(300,Nc))[-1] for mt in (0.0,0.05,0.1,0.13,0.26)]
        if any(l[i]<l[i+1]-1e-13 for i in range(4)): viol+=1
print(" violations:",viol)
print("== (ii) first dropped cell for mt=0.1: is k-sign ok, is g<=0? ==")
from p4_interchange_uniformity import wbar_grid, alpha_certificate_H, alpha_cell_ub, schur_full_seg
grid=wbar_grid(); 
i=[k for k in range(len(grid)-1) if abs(float(grid[k+1])-0.5391)<2e-3][0]
a,b=float(grid[i]),float(grid[i+1])
ca=alpha_certificate_H(a)[2]; cb=alpha_certificate_H(b)[2]
al=alpha_cell_ub(a,b,ca,cb,schur_full_seg(a,b,1))
mu2,k1,k2,k3,case=mu2_pencil_cell(a,b,0.1)
print(f" cell [{a:.4f},{b:.4f}] mu2={float(mu2):.5f} alpha_cell={float(al):.5f} g={float(mu2-al):.5f} ksign_ok={bool(k1>=0 and k2>=0 and k3>=0)} case={case}")
tl2=max(np.linalg.eigvalsh(H_mat(400,a)-0.1*Kc(400,Nc))[-2] for Nc in range(0,60))
tl1=min(np.linalg.eigvalsh(H_mat(400,a)-0.1*Kc(400,Nc))[-1] for Nc in range(0,60))
print(f" true at a: min_Nc lam1={tl1:.5f}  max_Nc lam2={tl2:.2e}  (alpha_cell/lam2 ~ {float(al)/tl2:.0f}x)")
