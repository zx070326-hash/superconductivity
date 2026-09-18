import numpy as np, os
os.environ["OMP_NUM_THREADS"]="1"
# A) pointwise positivity of psi(w) = 2w^2/(w^2+1) + pi*w/sinh(pi*w) - 1
def psi(w): 
    x=np.pi*w
    s=np.where(x<1e-8, 1.0, x/np.sinh(np.clip(x,0,700)))
    return 2*w*w/(w*w+1) + s - 1
w=np.concatenate([np.linspace(1e-6,3,300001), np.logspace(np.log10(3),3,2000)])
p=psi(w); i=np.argmin(p)
print("A) min over w>0 of psi(w) = %.6e at w=%.6f ;  psi(0+)=%.3e  psi(1)=%.4f  psi(inf)->1"%(p[i],w[i],psi(1e-6),psi(1.0)))
print("   small-w expansion coeff: 2 - pi^2/6 = %.6f ; psi(w)/w^2 at w=1e-3: %.6f"%(2-np.pi**2/6, psi(1e-3)/1e-6))
# B) c1 = psi_1(0) and gap ratio of the pure-attraction operator H(w)
def Hmat(N,w):
    u=1/np.sqrt(2*np.arange(N)+1); k=np.arange(N); n,m=np.meshgrid(k,k,indexing='ij')
    f=lambda q: w*w/(q*q+w*w)
    L1=np.diag(np.concatenate(([0.0],2*np.cumsum(f(np.arange(1,N)))))/(2*k+1))
    return -L1+np.where(n!=m,f(np.abs(n-m)),0.0)*np.outer(u,u)+f(n+m+1)*np.outer(u,u)
print("\nB) pure-attraction H(w):  mu1, mu2, c1^2=psi1(0)^2, c1^2*gamma, and threshold mu1_max s.t. mu1-c1^2*gamma<1")
for N in (400,800):
  print("  N=%d"%N)
  for w in (0.5,1,2,3,5,10,30,100):
    H=Hmat(N,w); ev,V=np.linalg.eigh(H); mu1,mu2=ev[-1],ev[-2]; c1=V[0,-1]
    g=mu1-mu2; 
    print("   w=%6.1f mu1=%8.4f mu2=%+9.2e  c1^2=%.4f  gamma=%.4f  c1^2*gamma=%.4f  (mu1-c1^2 g)/mu1=%.4f"%(w,mu1,mu2,c1*c1,g,c1*c1*g,(mu1-c1*c1*g)/mu1))
