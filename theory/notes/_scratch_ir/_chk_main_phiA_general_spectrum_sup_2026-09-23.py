# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）：**路线 1 的决定性探针** —— (Phi-A) 天花板在一般谱上的 sup 落在哪里？

背景：`ceiling/phiA_roomT_threshold_rescale_correction_2026-09-23.md` 定下
  lam_ceil(Einstein) = 3.2330 > lam_area = 2.5382，
但 lam_ceil 只对 **Einstein 族**成立，lam_area 是**一般谱**排除 ⟹ 不是同一作用域。
要把 (Phi-A) 抬到一般谱，需要

  (i)  mu1 <= mu1^E   （固定 (lam, wbar2)；AD1975 mu*=0 极值性，analytic）
  (ii) Delta >= 一致下界（未做；预算：须保住 Einstein Delta 的 >=76%）

**本脚本判 (i)+(ii) 的合体形式**：直接问
      sup_{谱 in 族} [ (Phi-A) ceiling ]  是否在 Einstein 处取到？
若是，Einstein 的数就直接界住整族，路线 1 成立（剩下只是证明）。
若否，量出超出多少 —— 这是路线 1 的**价码**，也是「该不该派活」的预算。

⚠ 已知风险（memory `mustar-resolved-ceiling-PhiA-2026-09-21`）：
  猜想 (FA)「固定 (lam, A) 时 Einstein 对 (Phi-A) 天花板极值」**已在 lam=3 被三反例证伪**。
  本脚本换的是**另一组固定量 (lam, wbar2)**（= 固定 eta/M，室温判据用的那一组），
  两者不是同一命题 —— 必须重新判，不得从 (lam,A) 的证伪推 (lam,wbar2)。

━━━━━━━━━ 族的参数化（两原子，固定 lam 与 wbar2） ━━━━━━━━━
alpha^2F = sum_j (lam_j w_j / 2) delta(w - w_j),  j=1,2
  lam  = lam_1 + lam_2                       （固定）
  wbar2^2 = (lam_1 w_1^2 + lam_2 w_2^2)/lam  （固定；= <w^2>，eta/M = lam wbar2^2/HB^2）
参数化：lam_2 = x*lam, lam_1 = (1-x)*lam, rho = w_2/w_1
        ⟹ w_1 = wbar2 / sqrt((1-x) + x rho^2),  w_2 = rho w_1.
Einstein = 退化成员（rho=1，或 x->0/1）—— **连续型退化集**，必须**定向**扫到
（CLAUDE.md：连续型退化集随机抽样永远打不到）。本脚本的网格显式含 rho=1 与 rho->1。

核：G(k) = sum_j lam_j wbar_j^2/(wbar_j^2 + k^2),  wbar_j := w_j/(2 pi T)。

━━━━━━━━━ 库仑截断约定（本脚本的关键选择，必须写明） ━━━━━━━━━
Einstein 情形仓库用 N_c = floor(r * wbar)（即 w_c = r * w_E）。一般谱没有唯一的 w_E，
本脚本取 **w_c = r * wbar2** ⟹ N_c = floor(r * wbar2/(2 pi T)) 对**全族相同**。
好处：把 memory `threshold-statistic-needs-calibrated-baseline` 记的 floor() 整数台阶
从族内比较中**完全移除**（同一个 N_c，不是每个成员一个）⟹ 比较是干净的。

严格性：算子形式与齐次性 analytic（见 `_chk_main_leavens_c_and_phiA_rescale_2026-09-23.py`）；
        本脚本一切数值为 **numerical-evidence**（float64, eigh, N=300）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh

HB = 64.65415134      # meV / sqrt(eV/A^2/amu)
KB = 8.617333262e-2   # meV/K


def Kfull(lams, wbars, N):
    """一般（多原子）谱的 lam=真实值 的对称化声子核。"""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    k = np.arange(2 * N + 2)
    G = np.zeros(2 * N + 2)
    for L, W in zip(lams, wbars):
        G += L * W * W / (W * W + k * k)
    Kp = G[np.abs(i - j)]
    Km = G[i + j + 1]
    S = (Kp - Km).sum(axis=1)
    return (Kp + Km) * np.outer(u, u) - np.diag(u * u * S), u


def bundle(lams, wbars, mus, Nc, N):
    """返回 (lam1_exact, ceiling, Delta, gamma, mu1, c1s, PA_ok)。"""
    K, u = Kfull(lams, wbars, N)
    uc = np.where(np.arange(N) < Nc, u, 0.0)
    M = K - 2.0 * mus * np.outer(uc, uc)
    l1 = float(eigh(M, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])
    ev, V = eigh(K)
    m1, m2 = ev[-1], ev[-2]
    c1s = float(V[:, -1] @ uc) ** 2
    un = float(uc @ uc)
    gam = m1 - m2
    den = gam - 2.0 * mus * un
    if den <= 0.0:
        return l1, None, None, gam, m1, c1s, False
    sp = 2.0 * mus * (un - c1s) / den
    d = 2.0 * mus * c1s / (1.0 + sp)
    return l1, m1 - d, d, gam, m1, c1s, True


def members(lam, wbarE, x, rho):
    """固定 (lam, wbar2) 的两原子成员；wbarE := wbar2/(2 pi T)。返回 (lams, wbars)。"""
    w1 = wbarE / np.sqrt((1.0 - x) + x * rho * rho)
    return np.array([(1.0 - x) * lam, x * lam]), np.array([w1, rho * w1])


def einstein_root(lam, mus, r, N, lo=0.05, hi=14.0, ngrid=30, nbis=40):
    """解 lam1(M_Einstein(wbarE)) = 1，返回 wbarE*（= wbar2/(2 pi Tc)）。"""
    def f(w):
        Nc = int(np.floor(r * w))
        l1, *_ = bundle(np.array([lam]), np.array([w]), mus, Nc, N)
        return l1 - 1.0
    ws = np.linspace(lo, hi, ngrid)
    v = np.array([f(w) for w in ws])
    idx = [t for t in range(len(ws) - 1) if np.sign(v[t]) != np.sign(v[t + 1])]
    if not idx:
        return None
    a, b, fa = ws[idx[0]], ws[idx[0] + 1], v[idx[0]]
    for _ in range(nbis):
        m = 0.5 * (a + b)
        fm = f(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def main():
    N, mus, r = 300, 0.13, 10.0
    print("=" * 96)
    print(f"参数: N={N}, mu*={mus}, r={r} (w_c = r*wbar2, 全族同一 N_c)")
    print("固定量: (lam, wbar2)  <=>  固定 (lam, eta/M)  —— 室温判据用的那一组")
    print()

    for lam in (2.0, 3.2953, 4.5):
        wE = einstein_root(lam, mus, r, N)
        if wE is None:
            print(f"lam={lam}: Einstein 无根，跳过")
            continue
        Nc = int(np.floor(r * wE))
        l1E, cE, dE, gE, m1E, c1E, okE = bundle(
            np.array([lam]), np.array([wE]), mus, Nc, N)
        print("=" * 96)
        print(f"lam = {lam:.4f}   Einstein 根 wbarE* = {wE:.6f}   N_c = {Nc}   "
              f"(核验 lam1_E = {l1E:.9f})")
        print(f"  Einstein:  ceiling = {cE:.6f}   Delta = {dE:.6f}   "
              f"gamma = {gE:.6f}   mu1 = {m1E:.6f}   c1^2 = {c1E:.6f}")
        print()
        print(f"{'x':>6}{'rho':>8}{'lam1':>12}{'ceiling':>12}{'Delta':>11}"
              f"{'gamma':>11}{'mu1':>11}{'D/D_E':>9}{'g/g_E':>9}{'ceil-ceilE':>12}")
        best = (-1e9, None)
        bestl1 = (-1e9, None)
        npa = 0
        rows = 0
        for x in (0.1, 0.25, 0.5, 0.75, 0.9):
            for rho in (1.0, 1.05, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0):
                lams, wbars = members(lam, wE, x, rho)
                l1, c, d, g, m1, c1s, ok = bundle(lams, wbars, mus, Nc, N)
                rows += 1
                if not ok:
                    npa += 1
                    print(f"{x:6.2f}{rho:8.2f}{l1:12.6f}{'(PA)失效':>12}"
                          f"{'-':>11}{g:11.6f}{m1:11.6f}{'-':>9}"
                          f"{g/gE:9.4f}{'-':>12}")
                    continue
                if c > best[0]:
                    best = (c, (x, rho))
                if l1 > bestl1[0]:
                    bestl1 = (l1, (x, rho))
                print(f"{x:6.2f}{rho:8.2f}{l1:12.6f}{c:12.6f}{d:11.6f}"
                      f"{g:11.6f}{m1:11.6f}{d/dE:9.4f}{g/gE:9.4f}{c-cE:12.3e}")
        print(f"  -> (PA) 失效 {npa}/{rows}")
        print(f"  -> max ceiling = {best[0]:.6f} at (x,rho)={best[1]}   "
              f"vs Einstein {cE:.6f}   超出 {best[0]-cE:+.3e}  "
              f"{'<<< Einstein 不是族内 sup' if best[0] > cE + 1e-9 else '(Einstein 是 sup)'}")
        print(f"  -> max lam1    = {bestl1[0]:.6f} at (x,rho)={bestl1[1]}   "
              f"vs Einstein {l1E:.6f}   超出 {bestl1[0]-l1E:+.3e}  "
              f"{'<<< Einstein 不是精确极值' if bestl1[0] > l1E + 1e-9 else '(Einstein 是精确极值)'}")
        print()


if __name__ == "__main__":
    main()
