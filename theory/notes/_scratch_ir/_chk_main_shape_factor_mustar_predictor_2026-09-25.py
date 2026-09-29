# -*- coding: utf-8 -*-
"""
主会话 2026-09-25：第三块——把 μ* 也作为「因子」加进三约束天花板，问它能不能变成**预测器**。

判据（先写死，再看数）：若 r := Tc_exp / Tc_B(λ, ω̄₁, ω̄₂, μ*) 在 29 行上的**相对散布**
明显小于 r_∞（头条）、r_E（λ 分辨）、r_B0（三约束 μ*=0），则「天花板 × 常数」就是一台
计算器，且它**永远不会越过物理上界**。散布用 CV = std/mean；对照是 AD75 公式自己
（|Tc_AD/Tc − 1| 的均值 3.85%，仓库已核）。

库仑约定：AD75 Table I 印 μ*(ω_ph) 与 ω_ph ⟹ ω_c := ω_ph（逐行，与 ν 无关 ⟹ 仿射性保持，
凸极点归约成立）。N_c = ⌊ω_c/(2πT)⌋ 随 T 变 ⟹ 按 sup 规则取根：从高 T 向下扫，取第一个
λ₁≥1 的区间再二分（记忆 coulomb-saturation-creates-a-pseudo-root：二分起点决定落哪个根）。

单原子归约在 μ*>0 下**重新验证**（不沿用 μ*=0 的【R】）：【R′】同形随机反证，缩小规模。

严格性：全部 numerical-evidence（float64）。μ* 是**拟合量**（McMillan–Rowell 反演同时拟合 λ、μ*）
⟹ 本表只量「散布」，不产生任何 PASS/FAIL；μ*>0 的天花板以「μ* 真值 ≥ 印刷值」为条件。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import math
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq, minimize_scalar
from multiprocessing import Pool

_d = os.path.dirname(os.path.abspath(__file__))


def _load(fn, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, fn))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


R = _load("_chk_main_shape_factor_real_data_tightness_2026-09-25.py", "R")
T1 = R.T1
NM = 160
n = np.arange(NM); u = 1/np.sqrt(2*n+1)
I, J = np.meshgrid(n, n, indexing="ij"); UU = np.outer(u, u); kk = np.arange(2*NM+2, dtype=float)


def lam1(lam, a, w, t, mus, wc):
    a = np.asarray(a, float); w = np.asarray(w, float)/t
    g = lam*(a[:, None]*w[:, None]**2/(w[:, None]**2 + kk[None, :]**2)).sum(0)
    S = g[0] + 2*np.concatenate([[0.], np.cumsum(g[1:NM])])
    K = (g[np.abs(I-J)] + g[I+J+1])*UU - np.diag(u*u*S)
    Nc = int(math.floor(wc/t))           # ω_c/(2πT) in units ω̄₂≡1
    uc = np.where(n < Nc, u, 0.0)
    K = K - 2*mus*np.outer(uc, uc)
    return eigh(K, eigvals_only=True, subset_by_index=[NM-1, NM-1])[0]


def tsup(lam, a, w, mus, wc):
    """sup{t : λ₁(t) ≥ 1}：几何下扫定位最高穿越，再二分。"""
    f = lambda t: lam1(lam, a, w, t, mus, wc) - 1.0
    ts = np.exp(np.linspace(math.log(3.0), math.log(0.01), 260))
    prev = ts[0]
    if f(prev) >= 0: return float("nan")
    for t in ts[1:]:
        if f(t) >= 0:
            return brentq(f, t, prev, xtol=1e-11, rtol=1e-11)
        prev = t
    return 0.0


def a_of_w(w, rho1):
    """单原子极点曲线（ω̄₂≡1）：λ 活跃 w≤ρ₁ ⟹ a=1；A 活跃 ρ₁≤w≤1/ρ₁ ⟹ a=ρ₁/w；M2 活跃 w≥1/ρ₁ ⟹ a=1/w²。
    μ*>0 时 Tc 对 w 不必单调（ω_c 固定），故三支**全扫**，不只取端点。"""
    return 1.0 if w <= rho1 else (rho1/w if w <= 1/rho1 else 1/(w*w))


def tB(lam, rho1, mus, wc):
    g = lambda lw: -math.exp(lw)*tsup(a_of_w(math.exp(lw), rho1)*lam, [1.], [1.], mus, wc/math.exp(lw))
    lws = np.linspace(math.log(0.3), math.log(3.0), 61)
    vals = [-g(x) for x in lws]
    k = int(np.argmax(vals)); best = (vals[k], math.exp(lws[k]))
    r = minimize_scalar(g, bounds=(lws[max(k-1, 0)], lws[min(k+1, len(lws)-1)]), method="bounded",
                        options=dict(xatol=1e-8))
    if -r.fun > best[0]: best = (-r.fun, math.exp(r.x))
    return best


def job_row(args):
    m, wl, w1, w2, lam, wph, mus, tc, eta = args
    rho = min(w1/w2, 1.0); wc = wph/w2
    tb0, _ = tB(lam, rho, 0.0, wc)
    tbm, a = tB(lam, rho, mus, wc)
    tEm = tsup(lam, [1.], [1.], mus, wc)
    to_K = w2/(2*math.pi)
    return m, lam, mus, tc/(0.1827262477*math.sqrt(lam)*w2), tc/(tb0*to_K), tc/(tEm*to_K), tc/(tbm*to_K), a


def job_R(args):
    lam, rho1, mus, wc, seed = args
    rng = np.random.default_rng(seed)
    tt, _ = tB(lam, rho1, mus, wc)
    worst = 0.0
    for _ in range(400):
        k = rng.integers(2, 4)
        w = np.exp(rng.uniform(math.log(0.1), math.log(5.0), k))
        a = rng.dirichlet(np.ones(k))
        if a.min() < 0.02: continue
        s = min(1.0, 1/a.sum(), 1/(a @ w**2), rho1/(a @ w)); a = a*s
        worst = max(worst, tsup(lam, a, w, mus, wc)/tt)
    return lam, rho1, worst


def main():
    print("="*100); print("【R′】μ*>0 下单原子是否仍是三约束集的 sup（μ*=0.10，ω_c=3ω̄₂ 与 1.8ω̄₂）"); print("="*100)
    tasks = [(l, r, 0.10, wc, 7*i+j) for i, (l, wc) in enumerate([(0.8, 3.0), (1.5, 1.8), (2.0, 3.0), (1.2, 1.8)])
             for j, r in enumerate([0.95, 0.85])]
    with Pool(18) as P:
        rr = P.map(job_R, tasks)
        for l, r, w in rr: print(f"  λ={l} ρ₁={r}: max t/t_B = {w:.10f}")
        gw = max(x[2] for x in rr)
        print(f"  ⟹ max = {gw:.10f}  {'单原子仍是 sup' if gw <= 1+1e-8 else '**单原子不是 sup —— 下表作废**'}")
        res = P.map(job_row, T1)
    print("\n" + "="*100); print("【P】29 行 AD75：四种归一的 r 与其相对散布"); print("="*100)
    print(f"{'材料':18s} {'λ':>5s} {'μ*':>6s} {'r_∞':>7s} {'r_B0':>7s} {'r_E(μ*)':>8s} {'r_B(μ*)':>8s}  w*/ω̄₂")
    for m, lam, mus, ri, rb0, rem, rbm, a in res:
        print(f"{m:18s} {lam:5.2f} {mus:6.3f} {ri:7.4f} {rb0:7.4f} {rem:8.4f} {rbm:8.4f}  {a:.3f}")
    cols = {"r_∞": 3, "r_B0 (三约束 μ*=0)": 4, "r_E(μ*) (λ,M2,μ*)": 5, "r_B(μ*) (三约束+μ*)": 6}
    print()
    for name, i in cols.items():
        v = np.array([x[i] for x in res])
        print(f"  {name:22s} mean {v.mean():.4f}  CV {v.std(ddof=1)/v.mean()*100:5.2f}%  "
              f"range [{v.min():.4f}, {v.max():.4f}]  max>1? {'**YES**' if v.max() > 1 else 'no'}")
    v = np.array([x[6] for x in res])
    pred = v.mean()
    err = np.abs(v/pred - 1)
    print(f"\n  计算器 Tc ≈ {pred:.4f}·Tc_B(λ,ω̄₁,ω̄₂,μ*)：|误差| 均值 {err.mean()*100:.2f}%  最大 {err.max()*100:.2f}%"
          f"   （对照：AD75 公式 3.85% / 12.62%；**该常数由同一 29 行拟合 ⟹ 样本内**）")
    # 留一交叉验证（诚实的样本外）
    loo = [abs(v[i]/np.delete(v, i).mean() - 1) for i in range(len(v))]
    print(f"  留一交叉验证 |误差| 均值 {np.mean(loo)*100:.2f}%  最大 {np.max(loo)*100:.2f}%")


if __name__ == "__main__":
    main()
