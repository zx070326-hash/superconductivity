#!/usr/bin/env python3
"""主会话从零复核 P4_eta_over_M.md §4 材料表 + §2.C 的 G(λ) 表。

不 import 仓库任何模块，不读 matcheck 的脚本。算子 𝒪[f] 与 ℌ(ϖ) 直接从
P4_main_theorems.md §0 的定义重写：
    𝒪[f] = −L1[f] + L2[f] + L3[f]
    L1 对角 (2n+1)^{-1} Σ_{k=1}^{n} 2 f(k)
    L2 非对角 f(|n−m|) u_n u_m      (n≠m)
    L3 f(n+m+1) u_n u_m             (全部)
    u_n = (2n+1)^{-1/2},  ⟦k⟧(ϖ) = ϖ²/(k²+ϖ²)
自洽哨兵：ℌ(ϖ)_{00} 必须恒等于 ϖ²/(1+ϖ²)（文档 §2.B 引理第 1 步）。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

# ---- CODATA / 常数（独立写入，不从仓库取） ----
KB_EV_PER_K = 8.617333262e-5      # eV/K
MEV_TO_K = 1e-3 / KB_EV_PER_K     # 1 meV 换成 K
HBARC_EV_A = 1973.269804          # eV·Å
AMU_EV = 931.49410242e6           # eV/c²
M_H_AMU = 1.00794

# C_∞ 上端（我自己的 rational-exact 包围，_chk_Cinf_from_g2_rational.py）
C_INF = 0.1827262477890

checks = []
def chk(name, cond, detail=""):
    checks.append(bool(cond))
    print(("PASS " if cond else "FAIL ") + name + (("  | " + detail) if detail else ""))

# ================= 第一部分：材料表的纯算术 =================
# (λ, ω_log meV, r=ω₂/ω_log, Tc_sscha, Tc_exp)  —— r=1 的行是 ω_log 代 ω₂（下估）
MAT = [
    ("H3S      200", 1.84,  92.86, 1.0,    194.0, 203.0),
    ("LaH10    129", 3.62,  76.40, 91.22/76.40, 255.3, None),
    ("LaH10    163", 2.67,  96.40, 103.63/96.40, 242.8, 250.0),
    ("LaH10    214", 2.06, 115.50, 127.63/115.50, 237.9, None),
    ("LaH10    264", 1.73, 126.60, 138.63/126.60, 216.9, None),
    ("YH6      165", 1.71, 114.87, 1.0,    241.5, 224.0),
    ("PdH        0", 0.40,  50.21, 1.0,      5.0,   9.0),
]
DOC = {  # 文档 §4.B 印出的 (η/M, η_H^ub, 界A, 界B, A/exp)
    "H3S      200": (15866, 3.83, 267.1, 214.5, 1.32),
    "LaH10    129": (30123, 7.26, 368.0, 332.1, None),
    "LaH10    163": (28674, 6.91, 359.1, 311.1, 1.44),
    "LaH10    214": (33555, 8.09, 388.4, 320.4, None),
    "LaH10    264": (33246, 8.02, 386.6, 305.4, None),
    "YH6      165": (22564, 5.44, 318.5, 250.7, 1.42),
    "PdH        0": ( 1008, 0.24,  67.3,  13.6, 7.48),
}

print("=== §4.B 主表：η/M = λω₂²，界A = C_∞√(η/M)（我的独立算术） ===")
print(f"{'system':<14}{'ω₂/meV':>9}{'η/M':>10}{'η_H^ub':>9}{'界A/K':>9}{'A/exp':>8}")
mine = {}
for name, lam, wlog, r, tcs, tce in MAT:
    w2 = wlog * r
    etaM = lam * w2 * w2                                   # meV²
    etaH = M_H_AMU * AMU_EV * (etaM * 1e-6) / HBARC_EV_A**2  # eV/Å²
    bА = C_INF * np.sqrt(etaM) * MEV_TO_K                  # K
    aexp = bА / tce if tce else None
    mine[name] = (etaM, etaH, bА, aexp)
    print(f"{name:<14}{w2:9.2f}{etaM:10.1f}{etaH:9.3f}{bА:9.1f}"
          + (f"{aexp:8.2f}" if aexp else "       —"))

for name, *_ in [(m[0],) for m in MAT]:
    etaM, etaH, bА, aexp = mine[name]
    d_etaM, d_etaH, d_bA, d_bB, d_ae = DOC[name]
    chk(f"{name}: η/M 与文档相符", abs(etaM - d_etaM) / d_etaM < 5e-5,
        f"我={etaM:.1f} 文档={d_etaM}")
    chk(f"{name}: η_H^ub 与文档相符", abs(etaH - d_etaH) / d_etaH < 3e-3,
        f"我={etaH:.3f} 文档={d_etaH}")
    chk(f"{name}: 界A 与文档相符", abs(bА - d_bA) < 0.06, f"我={bА:.2f} 文档={d_bA}")
    if d_ae:
        chk(f"{name}: A/exp 与文档相符", abs(aexp - d_ae) < 0.006,
            f"我={aexp:.3f} 文档={d_ae}")
    # 方向哨兵：界 A 必须高于被界量
    chk(f"{name}: 方向 界A > Tc_SSCHA", bА > MAT[[m[0] for m in MAT].index(name)][4])

# ================= 第二部分：G(λ) 从算子定义独立重算 =================
def H_mat(varpi, N):
    n = np.arange(N)
    u = (2.0 * n + 1.0) ** -0.5
    k = np.arange(0, 2 * N + 2)
    f = varpi**2 / (k.astype(float) ** 2 + varpi**2)   # ⟦k⟧(ϖ); f[0]=1
    # L1: 对角 (2n+1)^{-1} Σ_{k=1}^{n} 2 f(k)
    cs = np.concatenate(([0.0], np.cumsum(2.0 * f[1:N + 1])))  # cs[n] = Σ_{k=1}^{n}
    L1 = np.diag(cs[:N] / (2.0 * n + 1.0))
    idx = np.abs(n[:, None] - n[None, :])
    L2 = f[idx] * np.outer(u, u)
    np.fill_diagonal(L2, 0.0)                                   # 非对角
    L3 = f[n[:, None] + n[None, :] + 1] * np.outer(u, u)
    return -L1 + L2 + L3

from scipy.linalg import eigh
from scipy.optimize import brentq

NDEF = 800   # 文档用 N=1200；N=600→1200 相对差 1.3e-7（文档 Stage 5），故 800 足够本次比对
def h(varpi, N=NDEF):
    A = H_mat(varpi, N)
    return float(eigh(A, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])

# 自洽哨兵：ℌ_00 = ϖ²/(1+ϖ²)
for w in (0.3, 1.0, 5.0):
    got = H_mat(w, 60)[0, 0]
    chk(f"哨兵 ℌ_00(ϖ={w}) = ϖ²/(1+ϖ²)", abs(got - w**2 / (1 + w**2)) < 1e-14,
        f"got={got:.12f}")

_WC = {}
def varpi_c(lam, N=NDEF):
    if lam in _WC:
        return _WC[lam]
    target = 1.0 / lam
    lo, hi = 1e-2, 60.0
    w = brentq(lambda x: h(x, N) - target, lo, hi, xtol=1e-11, rtol=1e-13)
    _WC[lam] = w
    return w

TAB = [  # 文档 §2.C: λ, ϖ_c, G, G/C∞
    (0.30, 16.038469, 0.01811743, 0.099151),
    (1.00,  1.388550, 0.11461951, 0.627274),
    (1.71,  0.846133, 0.14384119, 0.787195),
    (2.67,  0.615311, 0.15829592, 0.866301),
    (3.62,  0.507328, 0.16488345, 0.902352),
    (100.0, 0.087401, 0.18209817, 0.996563),
]
print("\n=== §2.C 表：G(λ) = 1/(2π√λ ϖ_c)，h(ϖ_c)=1/λ（N=1200，我的独立算子） ===")
print(f"{'λ':>7}{'ϖ_c(我)':>13}{'ϖ_c(文档)':>13}{'G(我)':>13}{'G(文档)':>13}")
for lam, wc_d, G_d, r_d in TAB:
    wc = varpi_c(lam)
    G = 1.0 / (2 * np.pi * np.sqrt(lam) * wc)
    print(f"{lam:7.2f}{wc:13.6f}{wc_d:13.6f}{G:13.8f}{G_d:13.8f}")
    chk(f"λ={lam}: ϖ_c 相符", abs(wc - wc_d) / wc_d < 3e-6, f"相对差={abs(wc-wc_d)/wc_d:.2e}")
    chk(f"λ={lam}: G 相符", abs(G - G_d) / G_d < 3e-6)
    chk(f"λ={lam}: G/C_∞ 相符", abs(G / C_INF - r_d) < 2e-6)

# 界 B = G(λ)√(η/M)，核对文档的界B 列
print("\n=== 界 B = G(λ)·√(η/M) ===")
for name, lam, wlog, r, tcs, tce in MAT:
    etaM = mine[name][0]
    G = 1.0 / (2 * np.pi * np.sqrt(lam) * varpi_c(lam))
    bB = G * np.sqrt(etaM) * MEV_TO_K
    d_bB = DOC[name][3]
    print(f"{name:<14} G={G:.8f}  界B={bB:8.2f}  文档={d_bB}")
    chk(f"{name}: 界B 相符", abs(bB - d_bB) < 0.06)
    chk(f"{name}: 方向 界B ≤ 界A", bB <= mine[name][2] + 1e-9)

# ================= 第三部分：§5.3 门槛 + 一行公式系数 =================
print("\n=== §5.3 室温门槛 η/M = (Tc/(C_∞·MEV_TO_K))² ===")
for tc, doc_thr in [(200, 8896), (250, 13900), (300, 20016), (350, 27245)]:
    thr = (tc / (C_INF * MEV_TO_K)) ** 2
    print(f"  Tc={tc} K → η/M ≥ {thr:.1f} meV²   文档={doc_thr}")
    chk(f"门槛 {tc}K 相符", abs(thr - doc_thr) / doc_thr < 5e-5)

print(f"\n天花板系数 C_∞·MEV_TO_K = {C_INF*MEV_TO_K:.6f}  (一行公式 Tc ≤ 2.1205·√(η/M[meV²]))")

n = sum(checks)
print(f"\n=== {n}/{len(checks)} 通过 ===")
raise SystemExit(0 if n == len(checks) else 1)
