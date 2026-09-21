# -*- coding: utf-8 -*-
"""
主会话对 alpha-crossover 交付（alpha_sign_change_critical_mass_2026-09-20.md）的从零复核。
两个承重问题：
  (A) "M* 的存在性依赖内插方案" 是否成立？
      检验：有限差 alpha(D->T) 恰是 alpha(M) 在 [M_D,M_T] 上的对数平均
            alpha_fd = ∫ alpha(M) dlnM / ln(M_T/M_D)。
      凡是逐点经过数据的 C^1 内插都必须复现 alpha_fd；不复现的方案不是"另一种内插"，
      是没过数据点的拟合。三方案逐个算对数平均并与 alpha_fd 比对。
  (B) 把观测量从 alpha(D->T) 换成 M*（alpha 的零点位置）能否改善判别性？
      在单根处 dM*/dp = -(∂α/∂p)/(∂α/∂M)，分母对每个扰动 p 相同
      ⟹ |机制灵敏度|/|μ* 灵敏度| 在两种参数化下逐字相等（analytic）。数值验证之。

数据：Errea et al. arXiv:1305.7123 Table I（SSCHA 列）。独立重写，不 import 对方脚本。
跑： OMP_NUM_THREADS=1 PYTHONIOENCODING=utf-8 python _chk_iso_mstar_reverify_2026-09-20.py
"""
import math

M = {"H": 1.007825, "D": 2.014102, "T": 3.016049}
SSCHA = {"H": (0.40, 405.0), "D": (0.46, 304.0), "T": (0.48, 257.0)}
MU0 = 0.085
KEYS = ("H", "D", "T")
X = [math.log(M[k]) for k in KEYS]
LNLAM = [math.log(SSCHA[k][0]) for k in KEYS]
LNW = [math.log(SSCHA[k][1]) for k in KEYS]


def lnF(lam, mus):
    den = lam - mus * (1.0 + 0.62 * lam)
    return -math.log(1.2) - 1.04 * (1.0 + lam) / den


def dlnF_dlnlam(lam, mus):
    den = lam - mus * (1.0 + 0.62 * lam)
    return lam * (-1.04) * (den - (1.0 + lam) * (1.0 - 0.62 * mus)) / den ** 2


def alpha_fd(a, b, mus):
    """有限差 alpha = -Δln Tc/Δln M，Tc = ω_log·F。"""
    la, wa = SSCHA[a]
    lb, wb = SSCHA[b]
    lnTa = math.log(wa) + lnF(la, mus)
    lnTb = math.log(wb) + lnF(lb, mus)
    return -(lnTb - lnTa) / (math.log(M[b]) - math.log(M[a]))


# ---- 三种 lnλ(x), lnω(x) 的构造（x = ln M），各返回 (值, 导数) ----
def quad_lagrange(ys):
    def f(x):
        y = dy = 0.0
        for i in range(3):
            j, k = [t for t in range(3) if t != i]
            Li = (x - X[j]) * (x - X[k]) / ((X[i] - X[j]) * (X[i] - X[k]))
            dLi = ((x - X[k]) + (x - X[j])) / ((X[i] - X[j]) * (X[i] - X[k]))
            y += ys[i] * Li
            dy += ys[i] * dLi
        return y, dy
    return f


def piecewise(ys):
    def f(x):
        i = 0 if x <= X[1] else 1
        s = (ys[i + 1] - ys[i]) / (X[i + 1] - X[i])
        return ys[i] + s * (x - X[i]), s
    return f


def lsq_line(ys):
    xb = sum(X) / 3
    yb = sum(ys) / 3
    s = sum((x - xb) * (y - yb) for x, y in zip(X, ys)) / sum((x - xb) ** 2 for x in X)
    b = yb - s * xb
    return lambda x: (b + s * x, s)


def make_alpha(scheme, mus, offset=0.0):
    fl, fw = scheme(LNLAM), scheme(LNW)

    def alpha(Mv):
        x = math.log(Mv)
        lnlam, dlam = fl(x)
        _, dw = fw(x)
        return -dw - dlnF_dlnlam(math.exp(lnlam), mus) * dlam - offset
    return alpha


def log_average(alpha, a, b, n=4000):
    """∫ alpha dlnM / ln(Mb/Ma)，Simpson。"""
    xa, xb = math.log(M[a]), math.log(M[b])
    h = (xb - xa) / n
    # 端点向区间内侧偏 1e-12，避免分段方案在结点 D 处取到相邻段的导数
    s = alpha(M[a] * (1 + 1e-12)) + alpha(M[b] * (1 - 1e-12))
    for i in range(1, n):
        s += (4 if i % 2 else 2) * alpha(math.exp(xa + i * h))
    return s * h / 3 / (xb - xa)


def root(alpha, lo, hi):
    flo, fhi = alpha(lo), alpha(hi)
    if flo * fhi > 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if alpha(mid) * flo < 0:
            hi = mid
        else:
            lo, flo = mid, alpha(mid)
    return 0.5 * (lo + hi)


SCHEMES = {"quadratic-loglog": quad_lagrange, "piecewise-powerlaw": piecewise,
           "global-powerlaw-lsq": lsq_line}

print("=== (A) 存在性：有限差 alpha 是 alpha(M) 的对数平均 ===")
for a, b in (("H", "D"), ("D", "T")):
    print(f"  有限差 alpha_AD({a}->{b}) = {alpha_fd(a, b, MU0):+.4f}")
fd_DT = alpha_fd("D", "T", MU0)
fd_HD = alpha_fd("H", "D", MU0)
for name, sch in SCHEMES.items():
    al = make_alpha(sch, MU0)
    # 各方案是否经过数据点
    fl, fw = sch(LNLAM), sch(LNW)
    miss = max(max(abs(fl(x)[0] - y) for x, y in zip(X, LNLAM)),
               max(abs(fw(x)[0] - y) for x, y in zip(X, LNW)))
    r = root(al, M["H"], M["T"])
    print(f"  [{name:20s}] 数据点最大偏离 {miss:.2e} | 对数平均 [H,D]={log_average(al,'H','D'):+.4f} "
          f"[D,T]={log_average(al,'D','T'):+.4f} | 根 M*={'无' if r is None else f'{r:.4f}'}")
print(f"  ⟹ 经过数据点的两方案都复现 [D,T] 对数平均 = {fd_DT:+.4f}（正）与 [H,D] = {fd_HD:+.4f}（负）；")
print("     lsq 方案不过数据点，其'无根'量的是拟合对数据的偏离，不是 M* 的性质。")
print("     MVT：对数平均在 [H,D] 为负、在 [D,T] 为正 ⟹ 任何 C^1 数据内插在 (M_H,M_T) 内必有零点。")

print("\n=== (A′) quadratic-loglog 下 M*(μ*) 全程存在 ===")
for mus in (0.06, 0.08, 0.085, 0.10, 0.12, 0.15):
    r = root(make_alpha(quad_lagrange, mus), M["H"], M["T"] + 1.0)
    print(f"  μ*={mus:.3f}  M*={'无' if r is None else f'{r:.4f}'}")

print("\n=== (B) 判别性比值在换参数化下不变 ===")
DMU, DMECH = 0.09, 0.0638   # μ* 物理带宽；谐波对照行定位的谱形状残差
h = 1e-6
a0 = make_alpha(quad_lagrange, MU0)
Ms = root(a0, M["H"], M["T"] + 1.0)
dA_dmu = (make_alpha(quad_lagrange, MU0 + h)(Ms) - make_alpha(quad_lagrange, MU0 - h)(Ms)) / (2 * h)
dA_dmech = (make_alpha(quad_lagrange, MU0, offset=h)(Ms) - make_alpha(quad_lagrange, MU0, offset=-h)(Ms)) / (2 * h)
dA_dM = (a0(Ms + h) - a0(Ms - h)) / (2 * h)
Ms_mu = root(make_alpha(quad_lagrange, MU0 + h), M["H"], M["T"] + 1.0)
Ms_mech = root(make_alpha(quad_lagrange, MU0, offset=h), M["H"], M["T"] + 1.0)
dMs_dmu = (Ms_mu - Ms) / h
dMs_dmech = (Ms_mech - Ms) / h
print(f"  M* = {Ms:.4f} u   ∂α/∂M|M* = {dA_dM:+.4f} /u")
print(f"  ∂α/∂μ* = {dA_dmu:+.4f}   ∂α/∂δ_mech = {dA_dmech:+.4f}")
print(f"  dM*/dμ* = {dMs_dmu:+.4f} u   dM*/dδ_mech = {dMs_dmech:+.4f} u")
print(f"  隐函数检验 dM*/dμ* 应 = -(∂α/∂μ*)/(∂α/∂M) = {-dA_dmu / dA_dM:+.4f}  "
      f"dM*/dδ 应 = {-dA_dmech / dA_dM:+.4f}")
R_alpha = abs(dA_dmech * DMECH) / abs(dA_dmu * DMU)
R_M = abs(dMs_dmech * DMECH) / abs(dMs_dmu * DMU)
print(f"  判别性比值 |机制|/|μ*|：观测量=α(M*) → {R_alpha:.4f}；观测量=M* → {R_M:.4f}；商 = {R_M / R_alpha:.6f}")
print("  ⟹ 商恒为 1（analytic：两者分母 ∂α/∂M 相同，约掉）。换参数化只改条件数，不改判别性。")
