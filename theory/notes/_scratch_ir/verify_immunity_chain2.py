import numpy as np, math, os
os.environ["OMP_NUM_THREADS"]="1"
print("(i') q series with enough terms")
K=120
from fractions import Fraction
import mpmath as mp
mp.mp.dps=60
c=[2*mp.pi-mp.pi**3/6]+[mp.pi**(2*j-1)*((2*j+1)*2*j-mp.pi**2)/mp.factorial(2*j+1) for j in range(2,K)]
print("   all c_j>0 for j=1..%d: %s ; c_1=%.6f"%(K-1,all(x>0 for x in c),c[0]))
for w in (1e-3,0.4646,1.0,2.0,10.0,20.0):
    s=sum(c[j-1]*mp.mpf(w)**(2*j+1) for j in range(1,K)); ex=mp.pi*w*(1+w*w)-(1-w*w)*mp.sinh(mp.pi*w)
    print("   w=%7.4f  series=%s exact=%s relerr=%.2e"%(w,mp.nstr(s,12),mp.nstr(ex,12),float(abs(s-ex)/abs(ex))))
print("\n(v) Toeplitz+Hankel:  T_f + H_f  >=  fhat(pi) I   for Lorentz f (min eig minus fhat(pi))")
for N in (200,400):
  for w in (0.2,1.0,5.0,50.0):
    k=np.arange(N); n,m=np.meshgrid(k,k,indexing='ij'); f=lambda q: w*w/(q*q+w*w)
    S=f(np.abs(n-m))+f(np.minimum(n+m+1,10**9))
    fh=np.pi*w/np.sinh(np.pi*w)
    print("   N=%d w=%5.1f  lmin(T+H)=%.6f  fhat(pi)=%.6f  diff=%+.3e"%(N,w,np.linalg.eigvalsh(S)[0],fh,np.linalg.eigvalsh(S)[0]-fh))
