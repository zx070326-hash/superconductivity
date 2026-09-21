# -*- coding: utf-8 -*-
"""曲率 q=-0.525 (t=-7.01) 的稳健性：在我据此订正仓库「λ指数=1」之前，先排除三种伪影。
   R1 高λ端离群点 a-Pb0.45Bi0.55 (残差 -2.68s) 单点驱动？
   R2 形状变量 S=log(ω_log/ω̄₂) 与 λ 的非线性相关被二次项吸收？
   R3 只在合金上成立、单质上没有？（单质 8 行）
   R4 留一元素出：q 的符号与显著性稳不稳？
方向注记：双侧检验，不产生单边界判决。主会话 2026-09-20，从零。"""
import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, re
C=0.1827262477
T=[("Pb",56,65,1.55,7.2),("In",68,89,0.805,3.40),("Sn",99,121,0.72,3.75),("Hg",29,49,1.6,4.19),
   ("Tl",52,64,0.795,2.36),("Ta",132,148,0.69,4.48),("a-Ga",55,101,1.62,8.56),("b-Ga",87,129,0.97,5.90),
   ("Tl0.9Bi0.1",48,62,0.78,2.30),("Pb0.4Tl0.6",48,62,1.15,4.60),("Pb0.6Tl0.4",50,62,1.38,5.90),
   ("Pb0.8Tl0.2",50,61,1.53,6.80),("Pb0.6Tl0.2Bi0.2",48,58,1.81,7.26),("Pb0.9Bi0.1",50,60,1.66,7.65),
   ("Pb0.8Bi0.2",46,57,1.88,7.95),("Pb0.7Bi0.3",47,57,2.01,8.45),("Pb0.65Bi0.35",45,55,2.13,8.95),
   ("In0.9Tl0.1",63,86,0.85,3.28),("In0.73Tl0.27",55,77,0.93,3.36),("In0.67Tl0.33",57,79,0.90,3.26),
   ("In0.57Tl0.43",53,74,0.85,2.60),("In0.5Tl0.5",53,73,0.83,2.52),("In0.27Tl0.73",42,63,1.09,3.64),
   ("In0.17Tl0.83",45,63,0.98,3.19),("In0.07Tl0.93",49,63,0.89,2.77),("In2Bi",46,67,1.40,5.6),
   ("Sb2Tl7",37,58,1.43,5.2),("Bi2Tl",47,59,1.63,6.4),("a-Pb0.45Bi0.55",29,47,2.59,7.0)]
nm=np.array([r[0] for r in T]); wl=np.array([float(r[1]) for r in T]); w2=np.array([float(r[2]) for r in T])
lam=np.array([float(r[3]) for r in T]); tc=np.array([float(r[4]) for r in T])
fr=tc/(C*np.sqrt(lam)*w2); L=np.log(lam); S=np.log(wl/w2); y=np.log(fr)
def fit(X,yy):
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b; d=len(yy)-X.shape[1]; s2=float(r@r/d)
    return b,np.sqrt(np.diag(s2*np.linalg.inv(X.T@X))),np.sqrt(s2)
def quad(m):
    X=np.column_stack([np.ones(m.sum()),L[m],S[m],L[m]**2]); return fit(X,y[m])
print("="*90); print("R1. 删掉高λ端离群点 a-Pb0.45Bi0.55（残差 -2.68s，λ=2.59）"); print("="*90)
b,s,_=quad(np.ones(29,bool)); print(f"   全 29 行 : q = {b[3]:+.4f} ± {s[3]:.4f}   t = {b[3]/s[3]:+.2f}")
for drop in ["a-Pb0.45Bi0.55","Pb0.65Bi0.35","Ta","Sn"]:
    m=nm!=drop; b,s,_=quad(m); print(f"   删 {drop:<16}: q = {b[3]:+.4f} ± {s[3]:.4f}   t = {b[3]/s[3]:+.2f}   (n={m.sum()})")
print(); print("="*90); print("R2. 形状变量 S 是否本身随 λ 弯曲？（若是，L² 可能在替 S 干活）"); print("="*90)
bS,sS,_=fit(np.column_stack([np.ones(29),L,L**2]),S)
print(f"   S ~ 1 + logλ + (logλ)² :  线性 {bS[1]:+.4f}±{sS[1]:.4f}   二次 {bS[2]:+.4f}±{sS[2]:.4f} (t={bS[2]/sS[2]:+.2f})")
print(f"   corr(L, S) = {np.corrcoef(L,S)[0,1]:+.4f}")
b,s,_=fit(np.column_stack([np.ones(29),L,L**2]),y)
print(f"   **不含 S** 的二次拟合: q = {b[2]:+.4f} ± {s[2]:.4f}  t = {b[2]/s[2]:+.2f}  ⟹ 曲率{'不依赖 S' if b[2]/s[2]<-3 else '可能来自 S'}")
b,s,_=fit(np.column_stack([np.ones(29),L,S,L**2,S**2,L*S]),y)
print(f"   全二次（含 S², L·S）    : q_LL = {b[3]:+.4f} ± {s[3]:.4f}  t = {b[3]/s[3]:+.2f}")
print(); print("="*90); print("R3. 只在合金上成立？单质 8 行（无 S 项，n 太小）"); print("="*90)
pure=np.array([not re.search(r"[0-9]",n) for n in nm]); alloy=~pure
for tag,m in [("单质 8 行",pure),("合金 21 行",alloy)]:
    b,s,_=fit(np.column_stack([np.ones(m.sum()),L[m],L[m]**2]),y[m])
    print(f"   {tag:<12} n={m.sum():<3} λ∈[{lam[m].min():.2f},{lam[m].max():.2f}]  q = {b[2]:+.4f} ± {s[2]:.4f}  t = {b[2]/s[2]:+.2f}")
print(); print("="*90); print("R4. 留一元素出：q 的符号/显著性"); print("="*90)
print(f"{'删去':<6}{'n':>5}{'q':>11}{'±':>9}{'t':>8}")
ok=True
for e in ["Pb","Bi","Tl","In","Sn","Hg","Ta","Ga","Sb"]:
    m=np.array([not re.search(e+r"(?![a-z])",re.sub(r"^[ab]-","",n)) for n in nm])
    if m.sum()<8: print(f"{e:<6}{m.sum():>5}   (太少，跳过)"); continue
    b,s,_=quad(m); t=b[3]/s[3]; ok &= (b[3]<0 and t<-2)
    print(f"{e:<6}{m.sum():>5}{b[3]:>11.4f}{s[3]:>9.4f}{t:>8.2f}")
print(); print("="*90); print(f"判决：q 在全部扰动下{'保持负号且 |t|>2 ⟹ 曲率稳健' if ok else '有子样本失去显著性 ⟹ 见上表'}"); print("="*90)
