# -*- coding: utf-8 -*-
"""
主会话 2026-09-25：接 `_chk_main_shape_factor_budget_2026-09-25.py`。
那里的面枚举 + 局部精修在全部 40 格上给出的 argmax **都是单原子**（其余权重=0）。
本脚本做三件承重的事：

【S】单原子闭式。𝒞̄_ρ = {ν≥0: λ(ν)≤λ, M2(ν)≤M2, A(ν)≤A0}，ρ₁ := 2A0/(λ ω̄₂) = ⟨ω⟩/√⟨ω²⟩。
     单原子 (λ份额 a, 频率 w)（ω̄₂≡1）：a≤1, a w²≤1, a w≤ρ₁ ⟹ w ≤ min(ρ₁/a, 1/√a)。
     Einstein 的 Tc 对 w 线性：t*(aλ, w) = w·τ(aλ)，τ(x):=t*(x,1)。
     ⟹ t_1(λ,ρ₁) = max_{a∈(0,1]} min(ρ₁/a, 1/√a)·τ(aλ)      （一维，τ 单调）
【R】随机反证：真正的 2/3 原子成员（非退化权重）+ 多起点 Nelder–Mead，
     问有没有成员超过 t_1。**判据**：max ratio ≤ 1+1e−9 ⟹ 单原子是 sup（numerical-evidence）。
     抽样族必须给「非单原子」正质量——三原子权重全部 ≥ 0.02 才计入（否则就是在测单原子自己）。
【L】同信息类比较（记忆 same-information-class-comparison-vs-leavens）：𝒞̄_ρ ⊂ Leavens 集
     {A(ν)≤A0}，且 ⊂ 𝒞̄（M2 天花板的集）⟹ 新界 ≤ min(Leavens, Einstein)。
     **只有严格 < min 才是新信息**；报 new / min(...)。
     Leavens 1975：kTc ≤ 0.2309·A（μ*=0，其本身只是数值局部极大，见记忆）。
     单位换算：A0 = λρ₁ω̄₂/2 ⟹ t_L = 2π·0.2309·λρ₁/2。

严格性：【S】的归约（给定单原子）= analytic（Einstein 标度 + 约束代数）；
「单原子是 sup」= numerical-evidence（【R】）；【L】= numerical-evidence（Leavens 常数本身非严格）。
μ*=0 ⟹ 对一切 μ*≥0 都是上界（库仑 PSD，仓库 analytic）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import math
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from multiprocessing import Pool

_d = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("B", os.path.join(_d, "_chk_main_shape_factor_budget_2026-09-25.py"))
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
tstar, Cof = B.tstar, B.Cof

LEAVENS = 0.2309


def t1(lam, rho1):
    f = lambda a: -min(rho1/a, 1/math.sqrt(a))*tstar(a*lam, [1.], [1.])
    # 分两段：a≤ρ₁²（M2 活跃）与 a≥ρ₁²（A 活跃），各自有界标量极小
    best = (-f(1.0), 1.0)
    for lo, hi in ((1e-3, rho1*rho1), (rho1*rho1, 1.0)):
        if hi - lo < 1e-12: continue
        r = minimize_scalar(f, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-10))
        if -r.fun > best[0]: best = (-r.fun, r.x)
        for a in np.linspace(lo, hi, 41):
            v = -f(a)
            if v > best[0]: best = (v, a)
    return best


def feasible(a, w, rho1):
    return (a.min() >= 0 and a.sum() <= 1+1e-12 and a @ w**2 <= 1+1e-12 and a @ w <= rho1+1e-12)


def project(a, w, rho1):
    """把任意正 (a,w) 缩放进 𝒞̄_ρ（整体乘 s≤1 保方向）——得到可行成员。"""
    s = min(1.0, 1/a.sum(), 1/(a @ w**2), rho1/(a @ w))
    return a*s


def job_R(args):
    lam, rho1, seed, nrand, nstart = args
    rng = np.random.default_rng(seed)
    tt, _ = t1(lam, rho1)
    worst = 0.0; nvalid = 0
    for _ in range(nrand):
        k = rng.integers(2, 4)
        w = np.exp(rng.uniform(math.log(0.05), math.log(8.0), k))
        a = rng.dirichlet(np.ones(k))
        if a.min() < 0.02: continue
        a = project(a, w, rho1); nvalid += 1
        worst = max(worst, tstar(lam, a, w)/tt)
    # 多起点局部优化（参数：log a、log w，罚投影）
    for _ in range(nstart):
        k = 3
        x0 = np.concatenate([rng.normal(-1, 0.5, k), rng.uniform(math.log(0.2), math.log(3.0), k)])
        def neg(x):
            a = np.exp(np.clip(x[:k], -30, 0)); w = np.exp(np.clip(x[k:], -5, 5))
            a = project(a, w, rho1)
            return -tstar(lam, a, w)
        r = minimize(neg, x0, method="Nelder-Mead", options=dict(maxiter=1500, xatol=1e-9, fatol=1e-13))
        worst = max(worst, -r.fun/tt)
    return lam, rho1, worst, nvalid


def main():
    lams = [0.6, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0]
    rhos = [0.99, 0.95, 0.90, 0.85, 0.80]
    print("="*100); print("【S】单原子闭式 t_1(λ,ρ₁) vs 预算脚本 argmax（交叉核）"); print("="*100)
    prev = np.load(os.path.join(_d, "_out_shape_factor_budget_2026-09-25.npy"))
    mx = 0
    for l, r, t in prev:
        tt, a = t1(l, r); mx = max(mx, abs(t/tt - 1))
    print(f"  40 格：|t_budget/t_1 − 1| 最大 = {mx:.2e}（t_1 ≥ t_budget 应成立；网格只会更低）")
    worst_sign = max(t/t1(l, r)[0] for l, r, t in prev)
    print(f"  max t_budget/t_1 = {worst_sign:.12f}  （>1+1e−9 即 t_1 漏了面）")

    print("\n" + "="*100); print("【R】随机非单原子反证（权重 ≥0.02 的 2/3 原子 + 多起点 NM）"); print("="*100)
    tasks = [(l, r, 1000*i+j, 1500, 6) for i, l in enumerate(lams) for j, r in enumerate(rhos)]
    with Pool(18) as P: res = P.map(job_R, tasks)
    gw = max(x[2] for x in res); nv = sum(x[3] for x in res)
    for l, r, w, n in res:
        if w > 1 - 1e-4: print(f"  λ={l} ρ₁={r}: max ratio {w:.12f}（n={n}）")
    print(f"  总计 {nv} 个非退化随机成员 + {len(tasks)*6} 次 NM：max t/t_1 = {gw:.12f}"
          f"  ⟹ {'单原子是 sup（numerical-evidence）' if gw <= 1+1e-9 else '**有成员超过 t_1**'}")

    print("\n" + "="*100); print("【L】同信息类：new / min(Leavens·ρ₁, Einstein-M2)"); print("="*100)
    print("   λ     ρ₁   C_new/C_E   C_new/C_Leav   C_new/min   argmin(旧)   a*")
    rows = []
    for l in lams:
        tE = tstar(l, [1.], [1.])
        for r in [1.0] + rhos:
            tt, a = t1(l, r); tL = 2*math.pi*LEAVENS*l*r/2
            m = min(tE, tL)
            rows.append((l, r, tt/tE, tt/tL, tt/m))
            print(f"  {l:4.1f}  {r:.2f}   {tt/tE:.6f}    {tt/tL:.6f}      {tt/m:.6f}    "
                  f"{'Einstein' if tE < tL else 'Leavens '}    {a:.4f}")
    bad = [x for x in rows if x[4] > 1+1e-9]
    print(f"  一致性（new ≤ min 必须成立，因集合包含）：违反 {len(bad)}；"
          f"（Leavens 常数只是数值局部极大 ⟹ 若违反，先怀疑 0.2309，不是本脚本）")
    gain = [x for x in rows if x[1] < 1]
    print(f"  ρ₁<1 行里 new/min 的范围：[{min(x[4] for x in gain):.4f}, {max(x[4] for x in gain):.4f}]")


if __name__ == "__main__":
    main()
