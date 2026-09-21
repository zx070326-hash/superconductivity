"""[numerical-evidence] P_beta 的「区间化预算」：max_{n<=300} V_beta(n) 允许被外向放大多少倍，
12 格判据仍然通过？ 用 kappa*P_beta 模拟外向区间包封（包封 >= float 值，故 kappa>=1）。

方向审计：lam2 <= min(...)；把 P_beta 抬高 = 上界变松 = gap 变小 = ratio 变大 = 判据变难。
所以「kappa 处 PASS」蕴含「任何宽度因子 <= kappa 的外向包封 PASS」。单边方向正确（只产 PASS）。

链条 A（我要的）: min(kappa*P_beta, Lam_full(C_R=6), 0.6h)   —— Theta* 已被踢出
链条 B（参考）  : min(kappa*P_beta, Lam_full, Theta*, 0.6h)
"""
import os, io, contextlib, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    m = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(m)
    return m
gc = _load("_chk_grid_coarsen")
af = _load("_chk_12cell_analytic_form")
THETA_STAR = gc.THETA_STAR
CELLS = af.CELLS
C_R = 6.0                      # analytic (partial fractions), 本轮已从零重验 sup=3sqrt(3)<6
V2G_OVER_16 = af.V2G_OVER_16   # pad_D 常数 = V2[G]/16 = (3 sqrt3 /2)/16

def lam_full(w, NMAX=200000):
    n = np.arange(NMAX, dtype=float); s = (n + 0.5)/w
    padD = V2G_OVER_16/(w*w*s)
    padR = (1.0/(8.0*w*w))*C_R/s
    return float((af.theta(s) - np.arctan(s)/s + padD + padR).max())

# 预计算每格采样点的 h, P_beta, Lam_full（与 kappa 无关）
pre = []
for lo, hi, al in CELLS:
    row = []
    for w in np.geomspace(lo, hi, 11):
        bs = gc.base(w)
        row.append((bs, al, gc.peel_all_betas(w).min(), lam_full(w), 0.6*bs["h"]))
    pre.append(((lo, hi, al), row))

def worst(kappa, with_theta):
    W, tab = -np.inf, []
    for (lo, hi, al), row in pre:
        rc, who = -np.inf, ""
        for bs, al_, Pb, L, h06 in row:
            cands = {"P_b": kappa*Pb, "Lam": L, "0.6h": h06}
            if with_theta:
                cands["Th*"] = THETA_STAR
            k = min(cands, key=cands.get)
            r = gc.ratio(bs, al_, cands[k])
            if r > rc: rc, who = r, k
        tab.append(((lo, hi), rc, who)); W = max(W, rc)
    return W, tab

print(" kappa |  chainA(no Theta*)  |  chainB(with Theta*)")
for kap in [1.0, 1.5, 2, 3, 5, 8, 12, 20, 40, 100, 1e3, 1e6, np.inf]:
    a, _ = worst(kap, False); b, _ = worst(kap, True)
    print(f"{kap:7g} | {a:10.4f} {'PASS' if a<1 else 'FAIL':>5}   | {b:10.4f} {'PASS' if b<1 else 'FAIL':>5}")

print("\n-- 逐格归属 (chainA) --")
for kap in [1.0, 12.0, np.inf]:
    W, tab = worst(kap, False)
    print(f"kappa={kap:g}  worst={W:.4f}")
    print("   " + "  ".join(f"[{lo:.2f},{hi:.2f}]{r:.3f}/{w}" for (lo,hi), r, w in tab))
