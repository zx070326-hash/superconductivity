# -*- coding: utf-8 -*-
"""GATE Z2 — 交付日志里的 `mu_lo` 与我的复现差 ~1e−07（符号混杂）到底是谁的问题。

假设：日志用 `w̄={w:.6f}` 打印，**只有 6 位小数**；真正代入的是 `w = wE*(1+f)`。
把 w 从日志的 6 位值换成由 wE_rep(9 位) 重建的 `wE*1.01`，若差值塌下去 ⟹
「复现失败」是**我把日志的四舍五入值当输入**造成的，不是交付方的数。
（方向：dλ₁/dw>0，w 少 4.5e−07 相对 ⟹ mu_lo 少 ~1e−07，量级与符号都可预测。）
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util, math
from mpmath import iv
_d = os.path.dirname(os.path.abspath(__file__))
def _load(n):
    s = importlib.util.spec_from_file_location(n, os.path.join(_d, n + ".py"))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
BB = _load("_chk_main_box_branch_bound_certificate_2026-09-23"); KT = BB.KT
QS = _load("_chk_main_qsigma_coords_and_deficit_identity_2026-09-23")
R, NV, PAD, MUS, F = 10.0, 80, 200, 0.20, 0.010
# (lam, wE_rep 9位, 日志 w̄ 6位, 日志 mu_lo)
ROWS = [(2.0, 1.089141700, 1.100033, 1.021190594),
        (3.2953, 0.724463583, 0.731708, 1.020081705),
        (4.0, 0.634149822, 0.640491, 1.021883712),
        (6.0, 0.490832248, 0.495741, 1.025180308)]
def mulo(lam, w):
    Nc = int(math.floor(R*w)); Lr = NV + PAD
    v = BB.vstar(lam, w, 0.5, 0.5, MUS, Nc, NV)
    g = BB.G_box(lam, w, (0.5,)*4, Lr + NV + 2)
    return BB.kt_box(g, iv.mpf(lam), iv.mpf(lam)*iv.mpf(w)*iv.mpf(w), MUS, Nc, v, NV, Lr, -1e9)[1], Nc
print(f"{'lam':>7} {'w(日志6位)':>12} {'w=wE*1.01':>12} | {'Δ用日志w':>12} {'Δ用重建w':>12}  收缩")
for lam, wE, wlog, ref in ROWS:
    w2 = wE*(1.0+F)
    m1,_ = mulo(lam, wlog); m2,_ = mulo(lam, w2)
    d1, d2 = m1-ref, m2-ref
    print(f"{lam:7.4f} {wlog:12.6f} {w2:12.9f} | {d1:+12.3e} {d2:+12.3e}  {abs(d1)/max(abs(d2),1e-18):8.1f}x")
