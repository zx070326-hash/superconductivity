import numpy as np, os
os.environ["OMP_NUM_THREADS"]="1"
rng=np.random.default_rng(20260916)
N=300; k=np.arange(N); u=1/np.sqrt(2*k+1); nn,mm=np.meshgrid(k,k,indexing='ij')
def L1(fv):  # fv[j]=f(j), j=0..N-1
    return np.diag(2*np.concatenate(([0.0],np.cumsum(fv[1:])))/(2*k+1))
def P(fv):
    T=np.where(nn!=mm,fv[np.abs(nn-mm)],0.0); H=fv[np.minimum(nn+mm+1,N-1)]
    return (T+H)*np.outer(u,u)
def O(fv): return -L1(fv)+P(fv)
def bos(P_w,P_a):   # f(k)=sum_i a_i w_i^2/(w_i^2+k^2)
    return sum(a*w*w/(w*w+k.astype(float)**2) for w,a in zip(P_w,P_a))
def fhat_pi(P_w,P_a): return sum(a*np.pi*w/np.sinh(np.pi*w) for w,a in zip(P_w,P_a))

print("(i) q(w)=pi w(1+w^2)-(1-w^2)sinh(pi w) >=0 via nonneg odd series coefficients")
import math
c=[2*np.pi-np.pi**3/6]+[np.pi**(2*j-1)*((2*j+1)*2*j-np.pi**2)/math.factorial(2*j+1) for j in range(2,12)]
print("   c_1..c_5 =", ["%.4e"%x for x in c[:5]], " all>0:", all(x>0 for x in c))
ws=np.array([1e-3,0.1,0.4646,0.8,1.0,2.0,10.0])
qser=sum(c[j-1]*ws**(2*j+1) for j in range(1,12))
qex=np.pi*ws*(1+ws**2)-(1-ws**2)*np.sinh(np.pi*ws)
print("   series vs exact q:", np.max(np.abs(qser-qex)/(1+np.abs(qex))))

print("\n(ii) Lemma Z:  O'[f]=L1+P >= -f(0) e0e0^T  for random bosonic f (min eig of O'+f(0)e0e0^T)")
worst=1e9
for t in range(12):
    nw=rng.integers(1,5); W=10**rng.uniform(-1.2,1.7,nw); A=rng.uniform(0,1,nw)
    fv=bos(W,A); Op=L1(fv)+P(fv); E=np.zeros((N,N)); E[0,0]=fv[0]
    lm=np.linalg.eigvalsh(Op+E)[0]; worst=min(worst,lm)
    # also check the sharpened diagonal criterion 2f(1)+fhat(pi) >= f(0)
    crit=2*fv[1]+fhat_pi(W,A)-fv[0]
    print("   trial%2d nw=%d  f(0)=%.3f  2f(1)+fhat(pi)-f(0)=%+.4f  lmin(O'+f(0)e0e0^T)=%+.3e  lmin(O')=%+.4f"%(t,nw,fv[0],crit,lm,np.linalg.eigvalsh(Op)[0]))
print("   worst lmin =",worst)

print("\n(iii) THEOREM: lambda2(M) <= (1-c1^2)mu1 + c1^2 mu2 ,  M=-L1[phiZ]+P[f+ - f-], phiZ partial sums >= f+ + f-")
print("     [also reports actual lambda2(M) and lambda1(M)]")
bad=0
for t in range(20):
    nw=rng.integers(1,4); Wp=10**rng.uniform(-0.5,1.5,nw); Ap=rng.uniform(0,1.2,nw); fp=bos(Wp,Ap)
    nm_=rng.integers(1,4); Wm=10**rng.uniform(-1,2,nm_); Am=rng.uniform(0,2.0,nm_); fm=bos(Wm,Am)
    extra=rng.uniform(0,1.0)*bos(10**rng.uniform(-1,1.5,1),rng.uniform(0,1,1))  # extra nonneg Z weight
    phiZ=fp+fm+extra
    M=-L1(phiZ)+P(fp)-P(fm)
    Hp=O(fp); ev,V=np.linalg.eigh(Hp); mu1,mu2,c1=ev[-1],ev[-2],V[0,-1]
    maj=(1-c1*c1)*mu1+c1*c1*mu2
    e=np.linalg.eigvalsh(M); ok=e[-2]<=maj+1e-10
    bad+= (not ok)
    print("   t%2d mu1=%7.3f c1^2=%.3f  maj=%8.4f | l1(M)=%8.4f l2(M)=%+9.4f  ok=%s"%(t,mu1,c1*c1,maj,e[-1],e[-2],ok))
print("   violations:",bad)

print("\n(iv) rank-one secular bound: lambda2(H+rho e0e0^T) <= mu1-c1^2*gamma for all rho")
for w in (1.0,5.0,30.0):
    fv=bos([w],[1.0]); Hp=O(fv); ev,V=np.linalg.eigh(Hp); mu1,mu2,c1=ev[-1],ev[-2],V[0,-1]
    bnd=mu1-c1*c1*(mu1-mu2); row=[]
    for rho in (0.1,1,10,1e3,1e6):
        E=np.zeros((N,N)); E[0,0]=rho; row.append(np.linalg.eigvalsh(Hp+E)[-2])
    print("   w=%5.1f bound=%.5f  l2 at rho=0.1,1,10,1e3,1e6: "%(w,bnd)+" ".join("%.5f"%x for x in row)+"  ok=%s"%all(x<=bnd+1e-9 for x in row))
