"""有限 N 是 theta_max^inf 的下估（Cauchy 交错，需 exact S 且 N>=N_c）
 ==> N 可以取小，仍然声音。测最小可用 N，并验证 Rayleigh 商是有理函数。"""
import os; os.environ["OMP_NUM_THREADS"]="1"
import numpy as np
from scipy.linalg import eigh
from fractions import Fraction as F

def parts(w, nd):
    n = np.arange(nd); u = 1.0/np.sqrt(2.0*n+1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    G = w*w/(w*w + np.arange(2*nd+2)**2)
    Gc = w*w/(w*w + np.arange(nd)**2)
    S = Gc[0] + 2.0*(np.cumsum(Gc)-Gc[0])
    return (G[np.abs(i-j)]+G[i+j+1])*np.outer(u,u) - np.diag(u*u*S), u

def theta(w, mus, r, nd):
    K,u = parts(w,nd); Nc=int(np.floor(r*w)); B=np.eye(nd)
    if mus>0 and Nc>0:
        uc=np.where(np.arange(nd)<Nc,u,0.0); B=B+2.0*mus*np.outer(uc,uc)
    val,vec = eigh(K,B,subset_by_index=[nd-1,nd-1])
    return float(val[0]), vec[:,0]

W, MUS, R = 1.2, 0.13, 10.0
print("N-截断亏损（闭式 S，下估）：theta_max^N 相对 N=400 的亏损")
ref,_ = theta(W,MUS,R,400)
for nd in (16,24,32,48,64,100,200,400):
    t,_ = theta(W,MUS,R,nd)
    print(f"   N={nd:4d}  theta={t:.12f}  亏损={1-t/ref:.3e}  c={t/(np.pi*W):.9f}")

# ── 有理性验证：把 u 吸进向量 (y_n=u_n v_n) 后所有系数为有理 ──
print()
print("有理化检验（w=6/5, mu*=13/100, N_c=12）：Rayleigh 商是有理数？")
ND = 32
Wf, MUf = F(6,5), F(13,100); Nc = 12
def Gf(k): return Wf*Wf/(Wf*Wf + F(k*k))
Sf = []; acc = Gf(0)
for n in range(ND):
    Sf.append(Gf(0) + 2*sum(Gf(j) for j in range(1,n+1)))
_, vec = theta(W,MUS,R,ND)
DEN = 10**6
y = [F(round(float(vec[k])/np.sqrt(2*k+1)*DEN), DEN) for k in range(ND)]   # y_n = u_n v_n
num = sum((Gf(abs(i-j))+Gf(i+j+1))*y[i]*y[j] for i in range(ND) for j in range(ND))
num -= sum(Sf[i]*y[i]*y[i] for i in range(ND))
den = sum(y[i]*y[i]*F(2*i+1) for i in range(ND)) + 2*MUf*(sum(y[:Nc]))**2
rq = num/den
print(f"   Rayleigh(有理) = {float(rq):.12f}")
print(f"   float64 theta  = {theta(W,MUS,R,ND)[0]:.12f}")
print(f"   相对亏损       = {1-float(rq)/theta(W,MUS,R,ND)[0]:.3e}   (>=0 才声音)")
print(f"   c 下界(有理 theta / (pi w)) = {float(rq)/(np.pi*W):.9f}")
print(f"   有理数个数: y {ND} 个 + G {2*ND+2} 个 + S {ND} 个 = {ND+2*ND+2+ND}")
print(f"   分子分母位数: num {len(str(num.numerator))}/{len(str(num.denominator))}, "
      f"den {len(str(den.numerator))}/{len(str(den.denominator))}")
