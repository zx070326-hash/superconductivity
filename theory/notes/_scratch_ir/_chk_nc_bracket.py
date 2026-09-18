"""Verify certv2's claim that the 0.089 K width on the mu*=0.13 Tc interval comes
from the N_c = floor/ceil bracketing.  Main session, own implementation.
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
from scipy.linalg import eigh
from _main_indep_baseline import read_a2f, P_measure, F_T, O_of_f, A2F

K_PER_MEV = 11.604518

def lam1_mode(T, w, p, lam, om2, mus, r, N, mode):
    F = F_T(2*N, w, p, T)
    n, u, L1, M = O_of_f(F, N)
    vp = math.sqrt(om2)/(2*math.pi*T)
    Nc = int(math.floor(r*vp)) if mode=='floor' else int(math.ceil(r*vp))
    if mus > 0 and Nc > 0:
        uc = np.where(n < Nc, u, 0.0)
        M = M - 2.0*(mus/lam)*np.outer(uc, uc)
    return eigh(M, eigvals_only=True)[-1], Nc, vp

def Tc_mode(w,p,lam,om2,mus,r,N,mode,lo=1.0,hi=40.0):
    f=lambda T: lam*lam1_mode(T,w,p,lam,om2,mus,r,N,mode)[0]-1.0
    assert f(lo)>0>f(hi)
    for _ in range(200):
        m=0.5*(lo+hi)
        if f(m)>0: lo=m
        else: hi=m
    return 0.5*(lo+hi)

w,a2f,_ = read_a2f(A2F); lam,p,om2 = P_measure(w,a2f)
for N in (40,100):
    print(f"\nN={N}  (H3S-200 anharmonic, mu*=0.13, r=10)")
    out={}
    for mode in ('floor','ceil'):
        T=Tc_mode(w,p,lam,om2,0.13,10.0,N,mode)
        l1,Nc,vp=lam1_mode(T,w,p,lam,om2,0.13,10.0,N,mode)
        out[mode]=T*K_PER_MEV
        print(f"  N_c={mode:5s}: T={T:.6f} meV  Tc={T*K_PER_MEV:.4f} K"
              f"   varpi={vp:.4f}  r*varpi={10*vp:.4f}  N_c={Nc}")
    print(f"  bracket width = {abs(out['floor']-out['ceil']):.4f} K"
          f"   (certv2 reported 0.089 / 0.090 K)")
    print(f"  ordering: ceil (stronger Coulomb) must give the LOWER Tc -> "
          f"{'OK' if out['ceil'] < out['floor'] else '*** REVERSED ***'}")
