# -*- coding: utf-8 -*-
"""
主会话 2026-09-25：「加一个谱形因子能不能把天花板压紧？」——预算（先于任何认证）。

━━━ 问题 ━━━
仓库的天花板 = sup_{ν∈𝒞̄} Tc[ν]，𝒞̄ = {ν=α²F dω ≥ 0 : ∫(2/ω)dν ≤ λ, ∫2ω dν ≤ M2}。
算子 M[ν] 只经 G(k)=∫ g_k dν (k≥1) 依赖 ν，g_k(ω)=2ω/(ω²+(2πTk)²) ∈ C₀((0,∞))。
加第三条线性约束 ∫h dν = c（「新因子」）后极点 ≤3 原子（Richter–Rogosinski）——
**但只有当逃逸质量（软端 ω→0 / 硬端 ω→∞）不能在对 G 贡献 0 的同时满足该约束时，
新约束才会真正排除 Einstein**。两个候选：

  (F1) AD75 的 f₂ 因子 ρ_log := ω_log/ω̄₂，ln ω_log = (2/λ)∫ ln ω · dν/ω。
       预言（analytic 草证，本脚本 GATE L 数值核）：软端逃逸 —— 在 Einstein 上加
       λ 份额 ε、频率 ω0=exp(ln ρ_log/ε) 的原子，ε→0 时 ln ω_log 恰满足约束、
       对 G(k≥1) 贡献 ~λε ω0²/k² → 0 ⟹ sup = Einstein，**因子一分钱不买**。
  (F2) ρ₁ := ω̄₁/ω̄₂ = ⟨ω⟩/√⟨ω²⟩，ω̄₁=(2/λ)∫dν（面积 A=λω̄₁/2）。
       软逃逸质量 m=λεω0/2→0、硬逃逸 m≤M2/(2Ω)→0 ⟹ 两端逃逸对 ∫dν 都贡献 0
       ⟹ A 在 𝒞̄ 上弱*连续 ⟹ A<A_E 真排除 Einstein ⟹ 预言**严格压紧**。

约束用「A ≤ A0」（松弛版）：集合随 ρ₁ 嵌套 ⟹ C(λ,ρ₁) 对 ρ₁ 单调不减（analytic），
于是只需 ρ₁ 的**上界**即可用（方向干净）。极点：A 活跃 ⟹ ≤3 原子；不活跃 ⟹ ≤2 原子。

单位：ω̄₂≡1，t:=2πT/ω̄₂；C := kTc/√(λ⟨ω²⟩) = t*/(2π√λ)。μ*=0（约定无关：
μ*≥0 只压低 Tc，μ*=0 的 sup 对一切 μ*≥0 都是上界）。

严格性：F1 逃逸构造、F2 连续性、嵌套单调 = analytic（一行论证，见笔记）；
其余一切 = numerical-evidence（float64；有限 N 截断是 λ₁ 的下估 ⟹ fail-open，GATE N 量它）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
import math
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq, minimize
from multiprocessing import Pool

NMAT = int(os.environ.get("NMAT", "128"))


def _mk(N):
    n = np.arange(N); u = 1/np.sqrt(2*n+1)
    I, J = np.meshgrid(n, n, indexing="ij")
    return n, u, I, J, np.outer(u, u), np.arange(2*N+2, dtype=float)


_C = {}


def lam1(lam, a, w, t, N=None):
    N = N or NMAT
    if N not in _C: _C[N] = _mk(N)
    n, u, I, J, UU, k = _C[N]
    a = np.asarray(a, float); w = np.asarray(w, float)/t
    g = lam*(a[:, None]*w[:, None]**2/(w[:, None]**2 + k[None, :]**2)).sum(0)
    S = g[0] + 2*np.concatenate([[0.], np.cumsum(g[1:N])])
    K = (g[np.abs(I-J)] + g[I+J+1])*UU - np.diag(u*u*S)
    return eigh(K, eigvals_only=True, subset_by_index=[N-1, N-1])[0]


def tstar(lam, a, w, N=None):
    f = lambda t: lam1(lam, a, w, t, N) - 1.0
    lo, hi = 0.02, 2.0
    while f(hi) > 0: hi *= 2
    if f(lo) < 0: return 0.0
    return brentq(f, lo, hi, xtol=1e-12, rtol=1e-12)


def Cof(lam, t): return t/(2*math.pi*math.sqrt(lam))


# ───────── 家族 ─────────
def two_atom(w1, w2):
    """λ、M2 恰满足（ω̄₂=1）：a1 w1² + a2 w2² = 1，a1+a2=1，w1<1<w2。"""
    a1 = (w2*w2-1)/(w2*w2-w1*w1); return np.array([a1, 1-a1]), np.array([w1, w2])


def three_atom(w, rho1):
    V = np.vstack([np.ones(3), w, w*w]); b = np.array([1.0, rho1, 1.0])
    try: a = np.linalg.solve(V, b)
    except np.linalg.LinAlgError: return None
    return a if np.all(a >= -1e-14) else None


# ───────── 任务 ─────────
# 松弛集 𝒞̄_ρ = {Σa≤1, Σa w²≤1, Σa w≤ρ1}（λ、M2、A 三条都取 ≤；紧化理由见 09-23 笔记 §10）。
# 极点：k 原子、恰 k 条约束活跃（k=1,2,3），其余约束须满足。逐面枚举，不漏面。
ROWS = {"lam": lambda w: np.ones_like(w), "M2": lambda w: w*w, "A": lambda w: w}
RHS = {"lam": 1.0, "M2": 1.0}
FACES = {1: [("lam",), ("M2",), ("A",)], 2: [("lam", "M2"), ("lam", "A"), ("M2", "A")], 3: [("lam", "M2", "A")]}


def solve_face(w, act, rho1):
    V = np.vstack([ROWS[c](w) for c in act]); b = np.array([RHS.get(c, rho1) for c in act])
    try: a = np.linalg.solve(V, b)
    except np.linalg.LinAlgError: return None
    if np.any(a < -1e-14) or not np.all(np.isfinite(a)): return None
    a = np.clip(a, 0, None)
    for c in ("lam", "M2", "A"):
        if c not in act and ROWS[c](w) @ a > RHS.get(c, rho1) + 1e-12: return None
    return a


def job_sup(args):
    lam, rho1, ng = args
    ws = np.exp(np.linspace(math.log(0.03), math.log(30.), ng))
    best = (-1.0, None)
    import itertools
    for k in (1, 2, 3):
        for idx in itertools.combinations(range(ng), k):
            w = ws[list(idx)]
            for act in FACES[k]:
                a = solve_face(w, act, rho1)
                if a is None or a.sum() == 0: continue
                t = tstar(lam, a, w)
                if t > best[0]: best = (t, (f"{k}{''.join(c[0] for c in act)}", a, w, act))
    t0, (kind, a0, w0, act) = best
    LO, HI = math.log(0.01), math.log(100.)
    def neg(x):
        w = np.sort(np.exp(np.clip(x, LO, HI))); a = solve_face(w, act, rho1)
        return 1.0 if a is None or a.sum() == 0 else -tstar(lam, a, w)
    r = minimize(neg, np.log(w0), method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-13, maxiter=600))
    if -r.fun > t0:
        w = np.sort(np.exp(np.clip(r.x, LO, HI))); best = (-r.fun, (kind, solve_face(w, act, rho1), w, act))
    t, (kind, a, w, act) = best
    return lam, rho1, (t, (kind, a, w))


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("="*96); print("GATE 0  Einstein 基线 vs 仓库认证括号（μ*=0）"); print("="*96)
    CERT = {2.0: (0.148241256, 0.151236028), 3.2953: (0.161456556, 0.164718304),
            4.0: (0.164972121, 0.168304891), 6.0: (0.170378794, 0.173820789)}
    CE = {}
    for lam, (lo, hi) in CERT.items():
        c = Cof(lam, tstar(lam, [1.], [1.])); CE[lam] = c
        print(f"  λ={lam:<7} C_E={c:.9f}  认证括号 [{lo},{hi}]  {'IN' if lo < c < hi else '**OUT**'}")

    print("\n" + "="*96); print("GATE N  截断 N 收敛（λ=1, Einstein 与一个三原子成员）"); print("="*96)
    for N in (64, 128, 256):
        print(f"  N={N:4d}  t*_E={tstar(1.0,[1.],[1.],N):.12f}  "
              f"t*_3={tstar(1.0,[0.3,0.4,0.3],[0.3,0.9,1.4],N):.12f}")

    print("\n" + "="*96); print("GATE L  (F1) ω_log 因子：软端逃逸序列（λ=1.5，目标 ρ_log=0.80）"); print("="*96)
    lam = 1.5; tE = tstar(lam, [1.], [1.])
    for eps in (0.2, 0.1, 0.05, 0.02, 0.01, 0.005):
        w0 = math.exp(math.log(0.8)/eps)
        # 调 Einstein 频率使 M2 恰 = 1：(1-ε)wE² + ε w0² = 1
        wE = math.sqrt((1 - eps*w0*w0)/(1-eps))
        a = np.array([1-eps, eps]); w = np.array([wE, w0])
        rlog = math.exp(a @ np.log(w)); t = tstar(lam, a, w)
        print(f"  ε={eps:<6} ω0={w0:.3e}  ρ_log={rlog:.6f}  M2={a@w**2:.12f}  Tc/Tc_E={t/tE:.9f}")
    print("  ⟹ 预言：Tc/Tc_E → 1（ρ_log 恒 0.8 而上界不降）")

    print("\n" + "="*96); print("GATE F2  ρ₁=⟨ω⟩/√⟨ω²⟩ 压紧曲线 C_≤(λ,ρ₁)/C_E(λ)"); print("="*96)
    lams = [0.8, 1.2, 1.6, 2.0, 3.0]
    rhos = [1.0, 0.99, 0.97, 0.95, 0.92, 0.88, 0.84, 0.80]
    ng = int(os.environ.get("NG", "34"))
    tasks = [(l, r, ng) for l in lams for r in rhos]
    with Pool(int(os.environ.get("NPROC", "18"))) as P:
        res = P.map(job_sup, tasks)
    tab = {}
    for lam, rho1, (t, (kind, a, w)) in res:
        tab[(lam, rho1)] = (t, kind, a, w)
    for lam in lams:
        tE = tstar(lam, [1.], [1.])
        print(f"  ── λ={lam}  C_E={Cof(lam,tE):.6f}")
        for rho1 in rhos:
            t, kind, a, w = tab[(lam, rho1)]
            print(f"     ρ₁={rho1:.2f}  C/C_E={t/tE:.6f}  argmax={kind}-atom a={a} w={w}")
    np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "_out_shape_factor_budget_2026-09-25.npy"),
            np.array([(l, r, tab[(l, r)][0]) for l in lams for r in rhos]))


if __name__ == "__main__":
    main()
