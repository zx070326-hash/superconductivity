"""
d 维 KAY 算子构造（d=1,2,3）—— numerical-evidence 级探索脚本。

物理设定（见最终报告的详细推导，此处只放实现）：
  标准（d=3, KAY II 原版）算子 O[f] 作用在权重 u_n = (2n+1)^{-1/2} 上，
  u_n^2 = 1/(2n+1) 来自「费米面局部平坦态密度」假设下的能量壳积分
    chi_3(w_n) = N(0) * Integral d eps / (w_n^2+eps^2) = N(0) pi/|w_n|.
  一般 d 维、在费米面存在 van Hove 型奇异性的态密度 N_d(eps) ~ N(0)*(W/|eps|)^{(2-d)/2}
  （d=2 退化为对数：N_2(eps)=N(0)[1+c ln(W/|eps|)]），代入同一个能量壳积分给出
  一般的 chi_d(w_n)，本脚本把 tau_n^(d) := chi_d(w_n)/chi_3-normalization 代回
  同一个 O[f] 算子结构（f 本身——声子 Lorentzian 核——维度无关，只有「能量壳权重」
  u_n -> sqrt(tau_n^(d)) 变了）。

  tau_n^(3) = 1/(2n+1)                                   [flat DOS, 原始 KAY]
  tau_n^(2) = (1/(2n+1)) * [1 + kappa2 * ln(Lambda/(2n+1))]   [2D van Hove log]
  tau_n^(1) = kappa1 * Lambda^{1/2} * (2n+1)^{-3/2}           [1D nesting/band-edge 1/sqrt(eps)]

  Lambda := W/(pi T)（能带宽度截断，约定取 Lambda=1000，仅在 M=200 截断内 2n+1<=399<<Lambda，
  对数项恒正）；kappa1=kappa2=1（O(1) 建模常数，绝对值可吸收进耦合常数 lambda，
  这里只关心 n 依赖的形状，不做绝对定标）。

矩阵元（对齐 P4_main_theorems.md §0 的 O[f] 定义，把 u_n -> sqrt(tau_n^(d)) 代入）：
  K_{nm}^(d) = [f(|n-m|)(1-delta_nm) + f(n+m+1)] * sqrt(tau_n^(d) tau_m^(d))
               - delta_nm * tau_n^(d) * sum_{k=1}^n 2 f(k)
  f(k) := ⟦k⟧(w) = w^2/(k^2+w^2)   （Einstein 声子 Lorentzian 核，维度无关）

h^(d)(w) := lambda_max( K^(d)(w) )，M=200 截断，mpmath dps=30，mp.eigsy。

严格性级别：整条 tau_n^(d) 推导是 **numerical-evidence** 级物理建模（态密度形式是标准
van Hove 分类的简化模型，非从头计算；Lambda, kappa1, kappa2 是未定标的 O(1) 建模选择）。
d=3 一行是精确复现已有 analytic 结果（KAY II 原始算子），用作基线/自检。
"""
import time
from mpmath import mp, mpf, log, matrix, eigsy

mp.dps = 30

M = 200
LAMBDA = mpf(1000)     # 能带宽度截断 W/(pi T)，需 >> 2n+1 (n<M) 才能保持 d=2 对数项恒正
KAPPA1 = mpf(1)
KAPPA2 = mpf(1)


def f_lorentz(k, w):
    """Einstein 声子 Lorentzian 核 f(k)=w^2/(k^2+w^2)，维度无关。k 可以是 0。"""
    if k == 0:
        return mpf(1)
    k = mpf(k)
    return w * w / (k * k + w * w)


def tau_n(n, d):
    """d 维费米面能量壳权重 tau_n^(d)，见模块 docstring 推导。"""
    n = mpf(n)
    denom = 2 * n + 1
    if d == 3:
        return 1 / denom
    if d == 2:
        return (1 / denom) * (1 + KAPPA2 * log(LAMBDA / denom))
    if d == 1:
        return KAPPA1 * (LAMBDA ** mpf('0.5')) * denom ** mpf('-1.5')
    raise ValueError(d)


def build_K(M, w, d):
    w = mpf(w)
    tau = [tau_n(n, d) for n in range(M)]
    sqrt_tau = [t ** mpf('0.5') for t in tau]
    # cumulative sum of 2 f(k), k=1..n
    cum = [mpf(0)] * (M + 1)
    for k in range(1, M + 1):
        cum[k] = cum[k - 1] + 2 * f_lorentz(k, w)
    K = matrix(M, M)
    for n in range(M):
        for m in range(n, M):
            if n == m:
                val = f_lorentz(2 * n + 1, w) * tau[n] - tau[n] * cum[n]
            else:
                val = (f_lorentz(abs(n - m), w) + f_lorentz(n + m + 1, w)) * sqrt_tau[n] * sqrt_tau[m]
            K[n, m] = val
            K[m, n] = val
    return K


def h_of(w, d, M=M):
    K = build_K(M, w, d)
    E, _ = eigsy(K, eigvals_only=False)
    return max(E)


if __name__ == '__main__':
    varpis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    results = {1: [], 2: [], 3: []}
    t0 = time.time()
    for d in (3, 2, 1):
        for w in varpis:
            t1 = time.time()
            h = h_of(w, d)
            dt = time.time() - t1
            results[d].append((w, h))
            print(f"d={d} varpi={w:2d}  h={mp.nstr(h, 12)}   ({dt:.1f}s)", flush=True)
    print("TOTAL", time.time() - t0, "s")
    print()
    print("=== SUMMARY TABLE ===")
    print(f"{'varpi':>6} {'h^(1)':>16} {'h^(2)':>16} {'h^(3)':>16}")
    for i, w in enumerate(varpis):
        row = [mp.nstr(results[d][i][1], 10) for d in (1, 2, 3)]
        print(f"{w:>6} {row[0]:>16} {row[1]:>16} {row[2]:>16}")
