# -*- coding: utf-8 -*-
"""
主会话从零：**两矩包络**门禁。μ\*=0（刻意——避开库仑截断约定）。

问题（仓库台账记的空位）：Leavens 1975 用 0 阶矩 A、LAD/本仓库用 1 阶矩 M。
两界只在**不同** λ 饱和（Leavens 在 λ≈1.143、我们在 λ→∞）⟹ 中间区两者都松。
    Env(A, M) := sup{ Tc : ∫α²F dω = A , ∫ω·α²F dω = M }  ≤  min(两界)
问：中间区 Env 比 min 低多少？低得不够多 ⟹ 该路线不承重，关线。

━━ 标度化（先做，省掉一个自由度）━━
ω→sω 下：A→sA, M→s²M, λ→λ, Tc→sTc  ⟹  唯一无量纲量是
    m := M/A²   （与 memory 里 Cerqueira 那条「m:=M₁/A² 不是形状量」是同一个 m）
两界写成 Tc/A 的函数：
    Leavens : Tc/A ≤ c0            , c0 = 0.2308678（本脚本自己重算，不转引）
    本仓库  : Tc/A ≤ C_∞·√(2m)     , C_∞ = 0.1827262477
    交叉点  : m* = (c0/C_∞)²/2

━━ 方向审计（承重，单边界）━━
本脚本扫的是 **2 原子谱子族**，它给出的是 Env 的**下估**：
      sup_{2原子} Tc  ≤  Env  ≤  min(两界)
⟹ **只有「sup_{2原子} ≈ min(两界)」这一侧是声音结论**（夹逼），判决 = 「包络不承重，关线」。
   反过来「sup_{2原子} ≪ min」**不**证明 Env ≪ min（真极值谱可能不是 2 原子）。
   那一侧只算线索，要判 PASS 必须另证「Tc 在矩体上的极值点是 2 原子」（= Tc 对 α²F 的凸性，未证）。
   本脚本**不**输出 PASS。

严格性：numerical-evidence（float64；已做 N 收敛与两条独立锚点）。
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.linalg import eigh

C_INF = 0.1827262477            # 本仓库认证常数（CERT-D 上端）

# ───────────────────────────────────────────────────────────────────
# 1. 线性化 Eliashberg 算子：直接用 λ(k)，不拆成 λ×F(k)
#    （算子对 λ(k) 线性 ⟹ 原 Einstein 脚本的 λ·M(F) = M(λF) = M(G)）
#    来源：P4_main_theorems.md:17，与 _main_indep_einstein.py 同约定
# ───────────────────────────────────────────────────────────────────
def lam_of_k(kmax, T, atoms):
    """G[k] = λ(k) = 2 Σ_j a_j ω_j /(ω_j² + (2πTk)²).   G[0] = 2Σ a_j/ω_j = λ."""
    k = np.arange(kmax, dtype=float)
    x2 = (2.0 * np.pi * T * k) ** 2
    G = np.zeros(kmax)
    for a, w in atoms:
        G += 2.0 * a * w / (w * w + x2)
    return G


def top_eig(T, atoms, N):
    G = lam_of_k(2 * N + 2, T, atoms)
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    # L1[nn] = (Σ_{k=1}^{nn} 2G[k])/(2nn+1)，含 k=nn；L1[0]=0
    L1 = np.zeros(N)
    L1[1:] = 2.0 * np.cumsum(G[1:N]) / (2.0 * n[1:] + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    M = (G[np.abs(i - j)] + G[i + j + 1]) * u[i] * u[j] - np.diag(L1)
    M[n, n] -= G[0] * u * u                     # 去掉向量化多放的对角 F[0]u²
    # 只要最大本征值（比全谱分解快一个量级）
    return eigh(M, eigvals_only=True, subset_by_index=[N - 1, N - 1])[0]


def Tc(atoms, N=600, lo=1e-5, hi=3.0, iters=44):
    f = lambda T: top_eig(T, atoms, N) - 1.0
    if not (f(lo) > 0 > f(hi)):
        return np.nan
    for _ in range(iters):                      # 44 次几何二分 ⟹ 相对精度 ~1e-13
        mid = np.sqrt(lo * hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


# ───────────────────────────────────────────────────────────────────
# 2. 锚点 1：Leavens 常数 c0 = max_m Tc_E/A，Einstein 族 A=1 ⟹ ω_E=m, λ=2/m
#    Leavens 1975 印 0.2309 与 ω_E/A=1.750。这是第三方数、不是我的输入。
# ───────────────────────────────────────────────────────────────────
def tc_einstein_over_A(m, N=400):
    return Tc([(1.0, m)], N=N)          # A=a=1, M=a·ω=m ⟹ ω_E=m


print("=" * 86, flush=True)
print("锚点 1（外部校验）：c0 := max_m Tc_E/A   —— Leavens 1975 印 0.2309 @ ω_E/A=1.750")
print("=" * 86, flush=True)
a_, b_ = 1.40, 2.10                     # 黄金分割，缓存端点避免重复求值
fa_, fb_ = None, None
m1, m2 = a_ + 0.382 * (b_ - a_), a_ + 0.618 * (b_ - a_)
f1, f2 = tc_einstein_over_A(m1), tc_einstein_over_A(m2)
for _ in range(28):
    if f1 > f2:
        b_, m2, f2 = m2, m1, f1
        m1 = a_ + 0.382 * (b_ - a_); f1 = tc_einstein_over_A(m1)
    else:
        a_, m1, f1 = m1, m2, f2
        m2 = a_ + 0.618 * (b_ - a_); f2 = tc_einstein_over_A(m2)
m_opt = 0.5 * (a_ + b_)
c0 = tc_einstein_over_A(m_opt, N=800)
print(f"  我的 c0      = {c0:.7f}      Leavens 印 0.2309   相对差 {abs(c0-0.2309)/0.2309*100:.4f}%")
print(f"  我的 ω_E/A   = {m_opt:.5f}       Leavens 印 1.750    相对差 {abs(m_opt-1.750)/1.750*100:.4f}%")
print(f"  对应 λ_opt   = {2.0/m_opt:.5f}")
ok = abs(c0 - 0.2309) / 0.2309 < 2e-3
print(f"  ⟹ 求解器盲测{'通过' if ok else '**失败，禁止继续**'}", flush=True)
assert ok, "求解器未通过 Leavens 盲测，后续全部数字不可用"

print()
print("锚点 2（内部自洽）：λ→∞（m→0）时 Tc_E/(C_∞√(2m)) → 1 ?", flush=True)
for m, NN in [(0.20, 500), (0.10, 700), (0.05, 1000), (0.02, 1600)]:
    t = tc_einstein_over_A(m, N=NN)
    print(f"   m={m:5.2f}  λ={2/m:7.1f}  N={NN:5d}  Tc/A={t:.6f}"
          f"   Tc/(C_∞√2m)={t/(C_INF*np.sqrt(2*m)):.5f}", flush=True)

m_star = (c0 / C_INF) ** 2 / 2.0
print(f"\n  两界交叉点 m* = (c0/C_∞)²/2 = {m_star:.5f}   （此处 min 界 = c0 = {c0:.5f}）", flush=True)

# ───────────────────────────────────────────────────────────────────
# 3. 主计算：2 原子谱在固定 (A=1, M=m) 上的 sup Tc/A
#    参数化：a = 低频原子权重，ω₁ = 低频频率，ω₂ = (m−aω₁)/(1−a)
#    策略：粗网格（低 N）定位 ⟹ Nelder–Mead 细化（高 N），避免全网格高 N
# ───────────────────────────────────────────────────────────────────
from scipy.optimize import minimize


def tc_pair(a, w1, m, N):
    if not (0.0 < a < 1.0):
        return np.nan
    w2 = (m - a * w1) / (1.0 - a)
    if w1 <= 1e-6 or w2 <= 1e-6:
        return np.nan
    return Tc([(a, w1), (1.0 - a, w2)], N=N)


def sup_two_atom(m, N_coarse=250, N_fine=600, na=13, nw=13):
    best = (-1.0, None)
    for a in np.linspace(0.05, 0.95, na):
        for w1 in np.linspace(0.02 * m, min(0.98 * m / a, 6.0 * m), nw):
            t = tc_pair(a, w1, m, N_coarse)
            if np.isfinite(t) and t > best[0]:
                best = (t, (a, w1))
    if best[1] is None:
        return np.nan, None
    neg = lambda p: -(tc_pair(p[0], p[1], m, N_fine) if np.isfinite(
        tc_pair(p[0], p[1], m, N_fine)) else -1e9)
    r = minimize(neg, np.array(best[1]), method="Nelder-Mead",
                 options=dict(xatol=1e-5, fatol=1e-9, maxiter=120))
    a, w1 = r.x
    t = tc_pair(a, w1, m, N_fine)
    if not np.isfinite(t) or t < tc_pair(best[1][0], best[1][1], m, N_fine):
        a, w1 = best[1]
        t = tc_pair(a, w1, m, N_fine)
    return t, (a, w1, (m - a * w1) / (1.0 - a))


print()
print("=" * 86, flush=True)
print("主计算：固定 (A=1, M=m) 上 2 原子谱的 sup Tc/A，对比 min(两界)")
print("  ⚠ 方向：sup_{2原子} 是 Env 的**下估** ⟹ 只有「≈min」侧可判（关线），「≪min」不判 PASS")
print("=" * 86, flush=True)
print(f"  {'m':>7}{'Einstein':>11}{'sup 2原子':>11}{'min(两界)':>11}{'Env/min':>10}"
      f"{'距 min':>9}{'距 Einst':>10}   最优 (a, ω₁, ω₂)", flush=True)
for m in [0.30, 0.50, m_star, 1.00, 1.30, m_opt, 2.20]:
    ours = C_INF * np.sqrt(2.0 * m)
    mn = min(c0, ours)
    te = tc_einstein_over_A(m, N=600)
    ts, arg = sup_two_atom(m)
    tag = f"({arg[0]:.3f}, {arg[1]:.4f}, {arg[2]:.4f})" if arg else "-"
    print(f"  {m:7.4f}{te:11.6f}{ts:11.6f}{mn:11.6f}{ts/mn:10.4f}"
          f"{(1-ts/mn)*100:8.2f}%{(ts/te-1)*100:+9.2f}%   {tag}", flush=True)

print()
print("  解读：")
print("   · 'Env/min' 接近 1  ⟹ 两矩包络买不到东西，**关线**（这一侧是声音的）")
print("   · 'Env/min' 明显 <1 ⟹ 只是线索；要变成界还欠「Tc 在矩体上极值点是 2 原子」的证明")
print("   · '距 Einst' > 0 ⟹ Einstein 在固定 (A,M) 下**不**极值（与猜想 (FA) 被证伪同向）")

# ───────────────────────────────────────────────────────────────────
# 4. N 收敛（承重：低频原子逼近松原截断的 λ_c(N) 地板）
# ───────────────────────────────────────────────────────────────────
print()
print("=" * 86, flush=True)
print("N 收敛检验（承重量：m* 处的 sup 2 原子）")
print("=" * 86, flush=True)
for N in [300, 600, 1200]:
    t, arg = sup_two_atom(m_star, N_coarse=200, N_fine=N, na=9, nw=9)
    print(f"   N={N:5d}   sup_2atom(m*) = {t:.6f}   @ {arg}", flush=True)
print("=" * 86)
