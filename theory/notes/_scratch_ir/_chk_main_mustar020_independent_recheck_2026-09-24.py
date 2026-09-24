"""mustar020 (mu*=0.20) 交付的从零复核。不 import 该 agent 任何模块。
阈值用 w-参数化引理反解（lam*(w)=1/theta_max(w)，一次广义 eigh、零二分），
与它的二分求根是不同算法。S_n 用闭式（本轮 leavens 复核确立）。"""
import os; os.environ["OMP_NUM_THREADS"]="1"
import numpy as np
from scipy.linalg import eigh

C_INF, R, N = 0.1827262477, 10.0, 800

def parts(w, nd):
    n = np.arange(nd); u = 1.0/np.sqrt(2.0*n+1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    G  = w*w/(w*w + np.arange(2*nd+2)**2)
    Gc = w*w/(w*w + np.arange(nd)**2)
    S  = Gc[0] + 2.0*(np.cumsum(Gc)-Gc[0])
    return (G[np.abs(i-j)]+G[i+j+1])*np.outer(u,u) - np.diag(u*u*S), u

def Bmat(w, mus, nd):
    _, u = parts(w, nd); Nc = int(np.floor(R*w)); B = np.eye(nd)
    if mus > 0 and Nc > 0:
        uc = np.where(np.arange(nd) < Nc, u, 0.0); B = B + 2.0*mus*np.outer(uc, uc)
    return B, Nc

def lam_star(w, mus, nd):                      # lam*(w) = 1/theta_max(w)
    K, _ = parts(w, nd); B, Nc = Bmat(w, mus, nd)
    th = float(eigh(K, B, subset_by_index=[nd-1, nd-1], eigvals_only=True)[0])
    return 1.0/th, Nc

def lam1(w, lam, mus, nd):                     # lambda_1(lam*Khat - 2mu* uc uc^T)
    K, u = parts(w, nd); M = lam*K; Nc = int(np.floor(R*w))
    if mus > 0 and Nc > 0:
        uc = np.where(np.arange(nd) < Nc, u, 0.0); M = M - 2.0*mus*np.outer(uc, uc)
    return float(eigh(M, subset_by_index=[nd-1, nd-1], eigvals_only=True)[0])

def wE(lam, mus, nd, lo=0.05, hi=14.0, nb=90): # lam*(w)=lam 的根
    a, b = lo, hi
    for _ in range(nb):
        m = 0.5*(a+b)
        if lam_star(m, mus, nd)[0] < lam: b = m
        else: a = m
    return 0.5*(a+b)

C = lambda w, lam: 1.0/(2.0*np.pi*w*np.sqrt(lam))
MUS = 0.20
ROWS = [  # lam, wE_rep, w_up, Nc_up, worst, C_up_rep, w_lo, Nc_lo, mu_lo, C_lo_rep
 (2.0,    1.089141700, 1.078250, 10, 0.999838910, 0.104372372, 1.100033, 11, 1.021190594, 0.102305592),
 (3.2953, 0.724463583, 0.717219,  7, 0.999935012, 0.122242127, 0.731708,  7, 1.020081705, 0.119821491),
 (4.0,    0.634149822, 0.627808,  6, 0.999654397, 0.126754407, 0.640491,  6, 1.021883712, 0.124244418),
 (6.0,    0.490832248, 0.485924,  4, 0.999919035, 0.133713798, 0.495741,  4, 1.025180308, 0.131066000),
]
MU13 = {2.0:(0.114376567,0.116687205), 3.2953:(0.131240147,0.133361200),
        4.0:(0.135034160,0.137762123), 6.0:(0.142376142,0.144677171)}

print("="*104)
print("M1 Einstein 阈值（我的 w-参数化 vs 它报的）  M2 算术 C=1/(2*pi*w*sqrt(lam))  M5 N_c=floor(10w)")
bad = 0
for lam, wr, wu, ncu, wo, cur, wl, ncl, ml, clr in ROWS:
    wm = wE(lam, MUS, N)
    cu, cl, ce = C(wu, lam), C(wl, lam), C(wm, lam)
    d1, d2, d3 = abs(wm-wr)/wm, abs(cu-cur)/cu, abs(cl-clr)/cl
    ju, jl = int(np.floor(R*wu)), int(np.floor(R*wl))
    ok = (d1 < 3e-6) and (d2 < 3e-6) and (d3 < 3e-6) and ju == ncu and jl == ncl
    bad += not ok
    print(f" lam={lam:<7} wE 相对差={d1:.2e}  C_up 差={d2:.2e}  C_lo 差={d3:.2e}  "
          f"N_c {ju}/{jl} (报 {ncu}/{ncl}) {'OK' if ok else '** 不符 **'}")
print(f"  ==> M1/M2/M5: {4-bad}/4")

print()
print("M3 方向对撞：上侧须 lam1(E,w_up) <= worst；下侧须 mu_lo <= lam1(E,w_lo)")
b3 = 0
for lam, wr, wu, ncu, wo, cur, wl, ncl, ml, clr in ROWS:
    a, b = lam1(wu, lam, MUS, N), lam1(wl, lam, MUS, N)
    ou, ol = (a <= wo + 1e-12), (ml <= b + 1e-12)
    b3 += (not ou) + (not ol)
    print(f" lam={lam:<7} lam1(w_up)={a:.9f} <= worst={wo:.9f} {'OK' if ou else '**反**'}"
          f"   |  mu_lo={ml:.9f} <= lam1(w_lo)={b:.9f} {'OK' if ol else '**反**'}"
          f"   (上侧余 {wo-a:+.2e}, 下侧余 {b-ml:+.2e})")
print(f"  ==> M3: {8-b3}/8")

print()
print("M4 包含性 C_lower < C_E < C_upper   M6 收紧深度 f=0.010 核对   M7 C_upper < C_INF")
b4 = 0
for lam, wr, wu, ncu, wo, cur, wl, ncl, ml, clr in ROWS:
    wm = wE(lam, MUS, N); ce = C(wm, lam)
    fu, fl = 1.0 - wu/wm, wl/wm - 1.0
    ok = (clr < ce < cur) and abs(fu-0.010) < 2e-5 and abs(fl-0.010) < 2e-5 and cur < C_INF
    b4 += not ok
    print(f" lam={lam:<7} {clr:.9f} < C_E={ce:.9f} < {cur:.9f}   f_up={fu:.6f} f_lo={fl:.6f}"
          f"   C_up/C_INF={cur/C_INF:.5f} {'OK' if ok else '** 不符 **'}")
print(f"  ==> M4/M6/M7: {4-bad*0+0-b4}/4")

print()
print("M9 ★ N_c 台阶是否落在 [w_up, w_lo] 内？(它的 GATE J 只探了锚点 +-1e-6，未覆盖此段)")
print("    若跨台阶，须查 lam1(w)-1 在该区间的变号次数（赝根风险）")
for lam, wr, wu, ncu, wo, cur, wl, ncl, ml, clr in ROWS:
    steps = [k/R for k in range(int(np.floor(R*wu)), int(np.floor(R*wl))+2) if wu < k/R <= wl]
    ws = np.linspace(wu, wl, 121)
    v = np.array([lam1(w, lam, MUS, N) - 1.0 for w in ws])
    ncross = int(np.sum(np.sign(v[:-1]) != np.sign(v[1:])))
    ndec = int(np.sum(np.diff(v) < 0))
    print(f" lam={lam:<7} 区间内台阶 {steps if steps else '无'}   lam1-1 变号 {ncross} 次, "
          f"下降步 {ndec}/120  {'(单调升, 单根 OK)' if ncross==1 and ndec==0 else '** 需人工看 **'}")

print()
print("M8 ★ 跨 mu* 分离（第②项的真判据）：须 C_upper(0.20) < C_lower(0.13)")
b8 = 0
for lam, wr, wu, ncu, wo, cur, wl, ncl, ml, clr in ROWS:
    lo13, up13 = MU13[lam]
    ok = cur < lo13; b8 += not ok
    print(f" lam={lam:<7} C_up(0.20)={cur:.9f} < C_lo(0.13)={lo13:.9f} {'OK' if ok else '**重叠**'}"
          f"   余量 {lo13-cur:+.3e}  (锐化 {(1-cur/lo13)*100:.2f}%)")
print(f"  ==> M8: {4-b8}/4")
