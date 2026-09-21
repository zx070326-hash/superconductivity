# -*- coding: utf-8 -*-
"""Table II 的 lambda 反解——改用**完整** Eq.(34)（含 f1,f2），见 AD75 Eq.(35)-(38)。"""
import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.optimize import brentq

T2 = [  # name, wlog, w2, Tc, c6(mu*=.10), c7(mu*=.13)
    ("V",            222, 238,  5.30, 0.60, 0.67),
    ("Nb",           166, 183,  9.22, 0.85, 0.94),
    ("Mo",           236, 251,  0.92, 0.39, 0.45),
    ("Ta",           142, 153,  4.48, 0.67, 0.74),
    ("Ta(tunnel)",   132, 148,  4.48, 0.69, 0.76),
    ("W",            208, 220, 0.012, 0.25, 0.28),
    ("Nb3Sn(tunnel)",125, 163,  18.1, 1.67, 1.83),
]

def tc34(lam, wlog, w2, mus, use_f=True):
    den = lam - mus*(1.0 + 0.62*lam)
    if den <= 0: return 0.0
    if use_f:
        L1 = 2.46*(1.0 + 3.8*mus)
        r  = w2/wlog
        L2 = 1.82*(1.0 + 6.3*mus)*r
        f1 = (1.0 + (lam/L1)**1.5)**(1.0/3.0)
        f2 = 1.0 + (r - 1.0)*lam**2/(lam**2 + L2**2)
    else:
        f1 = f2 = 1.0
    return f1*f2*wlog/1.20*np.exp(-1.04*(1.0+lam)/den)

def inv(tc, wlog, w2, mus, use_f):
    g = lambda l: tc34(l, wlog, w2, mus, use_f) - tc
    a = mus/(1.0-0.62*mus) + 1e-9
    while tc34(a, wlog, w2, mus, use_f) >= tc and a < 60: a *= 1.0002
    return brentq(g, a+1e-9, 60.0, xtol=1e-13, rtol=1e-15)

print("="*104)
print("Table II lambda 反解：完整 Eq.(34)（含 f1 f2） vs f1=f2=1")
print("="*104)
print(f"  {'材料':>14} | {'列6印':>6}{'完整':>9}{'差':>8}{'(f=1)':>9} | {'列7印':>6}{'完整':>9}{'差':>8}{'(f=1)':>9}")
n6 = n7 = 0
for nm, wlog, w2, tc, c6, c7 in T2:
    F6, F7 = inv(tc, wlog, w2, .10, True),  inv(tc, wlog, w2, .13, True)
    N6, N7 = inv(tc, wlog, w2, .10, False), inv(tc, wlog, w2, .13, False)
    n6 += abs(F6-c6) <= 0.005; n7 += abs(F7-c7) <= 0.005
    print(f"  {nm:>14} | {c6:>6.2f}{F6:>9.4f}{F6-c6:>+8.4f}{N6:>9.4f} |"
          f" {c7:>6.2f}{F7:>9.4f}{F7-c7:>+8.4f}{N7:>9.4f}")
print(f"\n  完整 Eq.(34) 复现：列6 {n6}/7   列7 {n7}/7     判据 |差| <= 0.005（印数分辨力，2 位小数）")
print("\n  ★ 注意 Nb3Sn（lambda=1.67>1，f1f2 **不**可忽略）：f=1 给 2.02（错 21%），完整式给回印出值")
print("     ⟹ 这本身就证明 AD75 用的是**完整 Eq.(34) 反解**，不是别的来源。")

# ---- 印数舍入的区间检验：wlog 与 Tc 都按印出位数取半宽区间，看印出 lambda 是否被包住 ----
print()
print("="*104)
print("印数舍入区间检验：wlog / Tc 按印出位数取 +-半个末位，反解 lambda 的区间是否包住印出值")
print("="*104)
def halfwidth(x, s):
    """s = 印出的末位步长"""
    return s/2.0
ROWS = [  # name, wlog, step_wlog, w2, Tc, step_Tc, c6, c7
    ("V",            222, 1, 238,  5.30, 0.01, 0.60, 0.67),
    ("Nb",           166, 1, 183,  9.22, 0.01, 0.85, 0.94),
    ("Mo",           236, 1, 251,  0.92, 0.01, 0.39, 0.45),
    ("Ta",           142, 1, 153,  4.48, 0.01, 0.67, 0.74),
    ("Ta(tunnel)",   132, 1, 148,  4.48, 0.01, 0.69, 0.76),
    ("W",            208, 1, 220, 0.012, 0.001,0.25, 0.28),
    ("Nb3Sn(tunnel)",125, 1, 163,  18.1, 0.1,  1.67, 1.83),
]
tot = hit6 = hit7 = 0
for nm, wl, sw, w2, tc, st, c6, c7 in ROWS:
    tot += 1
    out = [nm]
    for mus, cp in ((0.10, c6), (0.13, c7)):
        vals = [inv(tc+dt, wl+dw, w2, mus, True)
                for dw in (-sw/2, sw/2) for dt in (-st/2, st/2)]
        lo, hi = min(vals), max(vals)
        # 印出值本身代表区间 [cp-0.005, cp+0.005]
        ok = not (hi < cp-0.005 or lo > cp+0.005)
        if mus == 0.10: hit6 += ok
        else:           hit7 += ok
        out.append(f"[{lo:.4f},{hi:.4f}] vs {cp:.2f} {'OK' if ok else 'MISS'}")
    print(f"  {out[0]:>14} | mu*=.10 {out[1]:<28} | mu*=.13 {out[2]}")
print(f"\n  区间相交：列6 {hit6}/{tot}   列7 {hit7}/{tot}")
print("  ⟹ Table II 的 lambda 列在**印数分辨力内完全由 (Tc, wlog, w2, mu*) 决定**，")
print("     零材料信息 ⟹ 与 `filling_factor_f` 同形的**循环量**。[numerical-evidence]")
