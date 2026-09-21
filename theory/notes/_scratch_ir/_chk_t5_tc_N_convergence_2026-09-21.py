#!/usr/bin/env python
"""N 收敛: t5 的 Tc 是不是固定 N=400 的截断赝物?

背景: t5 的 CSV 里 211/248 行 Tc_real=0 (二分两端同号, 返回下端 1e-12)。
分界干净: 失败全在 lam<=0.138, 成功全在 lam>=0.1442。
猜测: 固定 N 时可取的库珀对数只有 ln(N) (窗口 [nu_0, nu_N] 根本够不到声子标度 omega),
      不是 ln(omega/T)。==> 过零条件被截断决定, 而不是被 Tc 决定。

判据 (收敛性, 不是单边界): 若 Tc(N) 随 N 单调下降且不收敛 ==> 每个印出的 Tc 都是 N 的函数。
方向审计: Tc 被**低估**时, 上界判据 `Tc <= ceiling` 只会更容易 PASS
         ==> 「0/249 违反」正是单边界方向规则禁止据以下结论的那种 PASS。

等级 numerical-evidence。
"""
import os, math
import numpy as np
from scipy.linalg import eigvalsh

RY_TO_MEV = 13605.693122994
MEV_TO_K = 11.604518121550082
RY_TO_K = RY_TO_MEV * MEV_TO_K
ROOT = "D:/physics/superconductivity/theory/data/external/cerqueira_raw/extracted/batch-q"


def load(p):
    w, a = [], []
    for line in open(p):
        s = line.split()
        if line.lstrip().startswith('#') or len(s) < 2:
            continue
        try:
            x, y = float(s[0]), float(s[1])
        except ValueError:
            continue
        if x > 0:
            w.append(x); a.append(max(y, 0.0))
    return np.array(w), np.array(a)


def top_eig(w_ry, a2f, T_K, N):
    T_ry = T_K / RY_TO_K
    nu = 2.0 * np.pi * T_ry * np.arange(0, 2 * N + 2)
    integ = 2.0 * a2f[None, :] * w_ry[None, :] / (w_ry[None, :] ** 2 + nu[:, None] ** 2)
    lam_l = np.trapezoid(integ, w_ry, axis=1)
    n = np.arange(N)
    M = lam_l[np.abs(n[:, None] - n[None, :])] + lam_l[n[:, None] + n[None, :] + 1]
    cum = np.concatenate(([0.0], np.cumsum(lam_l[1:N])))
    M[np.diag_indices(N)] -= (2 * n + 1) + lam_l[0] + 2 * cum
    return eigvalsh(M, subset_by_index=[N - 1, N - 1])[0]


def find_tc(w_ry, a2f, N, T_hi):
    """守两端: 低端必须 >0, 高端必须 <0, 否则老实返回 nan (t5 只守了高端)。"""
    T_lo = 1e-14
    if top_eig(w_ry, a2f, T_lo, N) <= 0:
        return float('nan')
    if top_eig(w_ry, a2f, T_hi, N) >= 0:
        return float('inf')
    for _ in range(80):
        m = math.sqrt(T_lo * T_hi)            # 几何二分: Tc 跨数量级
        if top_eig(w_ry, a2f, m, N) > 0:
            T_lo = m
        else:
            T_hi = m
        if T_hi / T_lo < 1 + 1e-12:
            break
    return math.sqrt(T_lo * T_hi)


MATS = [("Al2CoNi_agm003191160", 5), ("Al2CoNi_agm003191160", 6)]
# 再挑两个「成功」侧与「失败」侧紧贴分界的
EXTRA = ["MgCu2Y_agm002225838", "NaS2Rh_agm005113872"]

print("N 收敛检验: Tc(N) 与截断比 nu_N/omega_max  (t5 全程固定 N=400)")
print("判据: Tc(N) 不收敛 ==> 印出的 Tc 是 N 的函数, 不是材料的函数\n")

names = [n for n in os.listdir(ROOT)][:0]
targets = ["Al2CoNi_agm003191160", "NaS2Rh_agm005113872", "AlRu2Ta_agm002635418"]

hdr = f"{'材料':<24}{'lam':>7}{'N':>7}{'Tc(K)':>14}{'nu_N/w_max':>12}{'Tc(N)/Tc(N/2)':>15}"
print(hdr); print("-" * len(hdr))
for name in targets:
    p = os.path.join(ROOT, name, "a2F.dos5")
    if not os.path.isfile(p):
        print(f"{name:<24} [缺文件]"); continue
    w_ry, a2f = load(p)
    lam = np.trapezoid(2 * a2f / w_ry, w_ry)
    w_max = w_ry.max()
    T_hi = 0.3 * w_max * RY_TO_K
    prev = None
    for N in (200, 400, 800, 1600, 3200):
        tc = find_tc(w_ry, a2f, N, T_hi)
        if not np.isfinite(tc):
            print(f"{name if N==200 else '':<24}{lam if N==200 else 0:>7.4f}{N:>7}"
                  f"{'无根(两端同号)':>14}{'':>12}{'':>15}")
            prev = None
            continue
        nu_N = 2 * np.pi * (tc / RY_TO_K) * N / w_max
        rat = f"{tc/prev:>15.4f}" if prev else f"{'-':>15}"
        print(f"{name if N==200 else '':<24}{lam if N==200 else 0:>7.4f}{N:>7}"
              f"{tc:>14.6g}{nu_N:>12.4f}{rat}")
        prev = tc
    print()
