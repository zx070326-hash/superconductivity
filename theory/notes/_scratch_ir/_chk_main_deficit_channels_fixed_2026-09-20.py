# -*- coding: utf-8 -*-
"""重算 deficit_decomposition_analysis_2026-09-20.md 的任务 B/C，
   用 **AD75 Eq.(35) 原式的 f2**（agent 版多了一个开方）。"""
import numpy as np
exec(open("theory/notes/_scratch_ir/_chk_main_lambda_exponent_artifact_2026-09-20.py",
          encoding="utf-8").read().split('print("=" * 78)')[0])
C=C_INF
def ad(l,w,b,m):
    l=np.asarray(l,float); w=np.asarray(w,float); b=np.asarray(b,float); m=np.asarray(m,float)
    L1=2.46*(1+3.8*m); r=b/w; L2=1.82*(1+6.3*m)*r
    f1=(1+(l/L1)**1.5)**(1/3); f2=1+(r-1)*l**2/(l**2+L2**2)
    return f1*f2*(w/1.20)*np.exp(-1.04*(1+l)/(l-m*(1+0.62*l)))

print("=== 0. AD 的 lam->inf 渐近常数 vs 我们的 C_inf ===")
for L in (1e3,1e5,1e7):
    v=ad(L,100.,100.,0.0)/(np.sqrt(L)*100.)
    print(f"   lam={L:.0e}  Tc/(sqrt(lam)*w) = {float(v):.7f}")
print(f"   解析极限 e^-1.04/(1.20*sqrt(2.46)) = {np.exp(-1.04)/(1.20*np.sqrt(2.46)):.7f}")
print(f"   我们的 C_inf                      = {C:.7f}   "
      f"=> AD 渐近高 {100*(np.exp(-1.04)/(1.20*np.sqrt(2.46))/C-1):+.2f}%")
print("   ^^^ AD 自报 Eq.(22) 常数 0.183；其 f1 形式实际给 0.1878 => 分解里 delta_FIN 可为负")

ceil = C*np.sqrt(lam)*wb2
TcE0 = ad(lam, wb2, wb2, 0.0)          # Einstein: wlog=wbar2 => f2=1 (**不受 bug 影响**)
TcP0 = ad(lam, wlog, wb2, 0.0)
TcP  = ad(lam, wlog, wb2, mus)
dFIN = 1-TcE0/ceil; dBRD = 1-TcP0/TcE0; dCOU = 1-TcP/TcP0; dRES = 1-Tc/TcP
print()
print("=== B. 四通道（AD75 原式 f2）===")
print(f"{'通道':<22}{'均值':>9}{'std':>8}{'min':>9}{'max':>9}   agent 报的均值")
for nm,v,a in (("delta_FIN 有限 lambda",dFIN,"33.7%"),("delta_BRD 谱展宽",dBRD,"6.1%"),
               ("delta_COU 库仑",dCOU,"37.3%"),("delta_RES AD 残差",dRES,"11.5%")):
    print(f"{nm:<22}{v.mean()*100:8.1f}%{v.std()*100:7.1f}%{v.min()*100:8.1f}%{v.max()*100:8.1f}%   {a}")
order=sorted([("库仑",dCOU.mean()),("有限lambda",dFIN.mean()),("谱展宽",dBRD.mean()),("AD残差",abs(dRES).mean())],key=lambda t:-t[1])
print("   排序:", " > ".join(f"{k} {v*100:.1f}%" for k,v in order))
print("   agent 排序: 库仑 37.3% > 有限lambda 33.7% > AD残差 11.5% > 谱展宽 6.1%")
print()
print("=== C. 可达上限（AD75 原式 f2, shape=0.892 最佳观测, mu*=0.10）===")
for L in (2.59,5,10,30,100):
    w=100.; b=w/0.892
    print(f"   lam={L:>6}  f = {float(ad(L,w,b,0.10)/(C*np.sqrt(L)*b)):.4f}"
          f"   (mu*=0: {float(ad(L,w,b,0.0)/(C*np.sqrt(L)*b)):.4f})")
print("   agent 报: lam=10 -> 0.909, lam=30 -> 0.953, lam=100 -> 0.972; 上限 0.85-0.90")
