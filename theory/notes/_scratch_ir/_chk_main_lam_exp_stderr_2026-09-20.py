# -*- coding: utf-8 -*-
"""0.9972(真实) vs 0.9053(AD 世界) 的 10% 差距是不是噪声？给系数标准误。"""
import numpy as np
exec(open("theory/notes/_scratch_ir/_chk_main_lambda_exponent_artifact_2026-09-20.py",
          encoding="utf-8").read().split('print("=" * 78)')[0])

C=C_INF
def ad(l,w,b,m):
    L1=2.46*(1+3.8*m); r=b/w; L2=1.82*(1+6.3*m)*r
    return (1+(l/L1)**1.5)**(1/3)*(1+(r-1)*l**2/(l**2+L2**2))*(w/1.20)*np.exp(-1.04*(1+l)/(l-m*(1+0.62*l)))
f=lambda t,l,b: t/(C*np.sqrt(l)*b)
A=np.column_stack([np.ones(29),np.log(lam),np.log(wlog/wb2)])

def fit(y):
    b,*_=np.linalg.lstsq(A,y,rcond=None); r=y-A@b
    s2=float(r@r)/(29-3); cov=s2*np.linalg.inv(A.T@A)
    return b,np.sqrt(np.diag(cov)),float(r@r)

y_obs=np.log(f(Tc,lam,wb2)); y_ad=np.log(f(ad(lam,wlog,wb2,mus),lam,wb2))
bo,so,_=fit(y_obs); ba,sa,_=fit(y_ad)
print(f"真实数据 : lam {bo[1]:.4f} +- {so[1]:.4f}   shape {bo[2]:.4f} +- {so[2]:.4f}")
print(f"AD 世界  : lam {ba[1]:.4f} +- {sa[1]:.4f}   shape {ba[2]:.4f} +- {sa[2]:.4f}")
d=bo[1]-ba[1]; sd=np.hypot(so[1],sa[1])
print(f"lam  差 {d:+.4f}   合并标准误 {sd:.4f}   => {d/sd:+.2f} sigma")
d2=bo[2]-ba[2]; sd2=np.hypot(so[2],sa[2])
print(f"shape差 {d2:+.4f}   合并标准误 {sd2:.4f}   => {d2/sd2:+.2f} sigma")
print()
print("配对检验（两个回归共用同一 A，差值 = 残差 y_obs-y_ad 对同一 A 的回归）：")
bd,sdd,_=fit(y_obs-y_ad)
print(f"  Delta-lam 指数 {bd[1]:+.4f} +- {sdd[1]:.4f}  => {bd[1]/sdd[1]:+.2f} sigma")
print(f"  Delta-shape指数 {bd[2]:+.4f} +- {sdd[2]:.4f}  => {bd[2]/sdd[2]:+.2f} sigma")
print("  ^^^ 配对才是对的检验：两个 y 共享 X，误差高度相关，独立标准误会高估不确定度")
print()
print("H0: 'lam 指数 = 1' 在真实数据上：", f"{(bo[1]-1)/so[1]:+.2f} sigma  => 不可拒绝" if abs((bo[1]-1)/so[1])<2 else "拒绝")
print("H0: 'lam 指数 = 1' 在 AD 世界上：", f"{(ba[1]-1)/sa[1]:+.2f} sigma")
