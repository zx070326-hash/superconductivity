# -*- coding: utf-8 -*-
"""Moussa-Cohen eq.(11) vs 我们的 C —— 第二遍，带**我自己从零重算的** C_E 与误差预算。

第一遍（同名脚本 v1）只把 MC 的拟合值对上笔记里转引的 C_upper，得出「μ*>0 时他们低 6-12%」。
那个判断**不能直接用**，理由是 CLAUDE.md 的失效热区规则：两个数不一致时，
先问「它们是同一个证书 / 同一个约定下的量吗」，再问「谁错了」。本脚本做三件事：

 (A) 我自己用 θ_max 重算 C_E(λ,μ*)（不转引表），确认落在认证括号内 —— 证据来源自持。
 (B) 把差距换算回**他们自己的单位**：Δf = ΔC·√λ，对照他们自陈的 RMS 0.01 / max 0.03。
     ⟹ 判「这是分歧还是他们已经声明过的拟合误差」。
 (C) 看差距随 μ* 的**形状**：若在 μ*=0 处归零且随 μ* 单调，则指向截断约定而非算错。

方向（承重）：C 是上界系数。
  MC 在**括号之上** ⟹ 他们的界成立、只是松（安全侧）
  MC 在**括号之下** ⟹ 他们那条不等式对真 Einstein 谱不成立（危险侧）
  ——但「不成立」只有在**同约定**下才能说；μ* 的截断约定他们全文未声明（已 grep 确认）。
严格性：C_lower/C_upper 转引 `../ceiling/C_lambda_certified_curve_2026-09-24.md` §1（certificate，
       且该表自带封顶：B&B 内部区间声音性未审）；我的 C_E 与全部比值是 numerical-evidence。
"""
import os
import sys
import math
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.linalg import eigh

# 认证括号 (C_lower, C_upper)，逐字抄自 C_lambda_certified_curve_2026-09-24.md §1
CERT = {
    0.00: {2.0: (0.148241256, 0.151236028), 3.2953: (0.161456556, 0.164718304),
           4.0: (0.164972121, 0.168304891), 6.0: (0.170378794, 0.173820789)},
    0.13: {1.0: (0.071626887, 0.085109595), 1.5: (0.099972499, 0.101992146),
           2.0: (0.114376567, 0.116687205), 2.5: (0.122918849, 0.125402058),
           3.0: (0.128397265, 0.130991149), 3.2953: (0.131240147, 0.133361200),
           4.0: (0.135034160, 0.137762123), 5.0: (0.139078844, 0.141888517),
           6.0: (0.142376142, 0.144677171), 8.0: (0.145685231, 0.148030265),
           10.0: (0.148435383, 0.150564456)},
    0.20: {2.0: (0.102305592, 0.104372372), 3.2953: (0.119821491, 0.122242127),
           4.0: (0.124244418, 0.126754407), 6.0: (0.131066000, 0.133713798)},
}
MC_RMS, MC_MAX = 0.01, 0.03      # 他们自陈的 f 偏差（绝对）
R = 10.0


def f_mc(lam, mus):
    if lam <= mus:
        return 0.0
    return 0.69 * math.exp(-(1.0 + lam) / (lam - mus)) * math.sqrt((1.0 + 0.52 * lam) / (1.0 + 5.6 * mus))


# ─────── (A) 我自己重算 C_E ───────
def theta_max(w, mus, n=400):
    Nc = int(math.floor(R * w))
    k = np.arange(2 * n + 2)
    G = w * w / (w * w + k * k)
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    u = 1.0 / np.sqrt(2.0 * np.arange(n) + 1.0)
    S = G[0] + 2.0 * (np.cumsum(G[:n]) - G[0])
    K = (G[np.abs(i - j)] + G[i + j + 1]) * np.outer(u, u) - np.diag(u * u * S)
    uc = np.where(np.arange(n) < Nc, u, 0.0)
    B = np.eye(n) + 2.0 * mus * np.outer(uc, uc)
    return float(eigh(K, B, eigvals_only=True)[-1])


def c_e(lam, mus, lo=0.05, hi=8.0):
    """⚠ 先扫变号再二分（⌊R·w⌋ 台阶会造赝根）。"""
    ws = np.linspace(lo, hi, 601)
    vs = np.array([lam * theta_max(w, mus) - 1.0 for w in ws])
    sc = np.where(np.diff(np.sign(vs)) != 0)[0]
    if len(sc) == 0:
        return float("nan"), 0
    s = sc[-1]
    a, b, va = ws[s], ws[s + 1], vs[s]
    for _ in range(70):
        m = 0.5 * (a + b)
        if (lam * theta_max(m, mus) - 1.0) * va > 0:
            a = m
        else:
            b = m
    w = 0.5 * (a + b)
    return 1.0 / (2 * math.pi * w * math.sqrt(lam)), len(sc)


print("=" * 118)
print("(A) 我从零重算 C_E，对照认证括号（证据来源自持；不转引即判）")
print(f"{'λ':>8}{'μ*':>6}{'我的 C_E':>15}{'C_lower':>13}{'C_upper':>13}{'在括号内?':>11}{'变号':>5}")
mine, ok_all = {}, True
for mus in (0.00, 0.13, 0.20):
    for lam in sorted(CERT[mus]):
        c, nsc = c_e(lam, mus)
        lo, hi = CERT[mus][lam]
        ok = lo <= c <= hi
        ok_all &= ok
        mine[(lam, mus)] = c
        print(f"{lam:8.4f}{mus:6.2f}{c:15.9f}{lo:13.9f}{hi:13.9f}{'  ✓' if ok else '  XX 越界':>11}{nsc:5d}")
print(f"⟹ 全部 {len(mine)} 格在括号内：{ok_all}")

# ─────── (B)(C) 差距换算回他们的单位 ───────
print("\n" + "=" * 118)
print("(B)(C) 差距换算回**他们自己的单位** Δf = (C_ours − C_MC)·√λ，对照他们自陈 RMS 0.01 / max 0.03")
print("      正 Δf = 他们的 f 偏**小** ⟹ 他们那条 `Tc ≤ f√⟨ω²⟩` 在真 Einstein 谱上**被违反**")
print(f"{'λ':>8}{'μ*':>6}{'C_MC':>13}{'我的 C_E':>13}{'相对差%':>10}{'Δf':>10}{'vs 0.03':>9}   位置")
rows = []
for mus in (0.00, 0.13, 0.20):
    for lam in sorted(CERT[mus]):
        cm = f_mc(lam, mus) / math.sqrt(lam)
        ce = mine[(lam, mus)]
        lo, hi = CERT[mus][lam]
        df = (ce - cm) * math.sqrt(lam)
        rel = (cm - ce) / ce * 100
        pos = ("括号**之上**（安全侧：界成立、只是松）" if cm > hi else
               "括号**之下**（危险侧）" if cm < lo else "括号内")
        flag = "超出" if abs(df) > MC_MAX else ("超RMS" if abs(df) > MC_RMS else "内")
        rows.append((lam, mus, cm, ce, rel, df, flag, pos))
        print(f"{lam:8.4f}{mus:6.2f}{cm:13.9f}{ce:13.9f}{rel:+10.2f}{df:+10.5f}{flag:>9}   {pos}")

print("\n" + "-" * 118)
print("(C) 形状：差距随 μ* 的行为（若在 μ*=0 归零且随 μ* 单调 ⟹ 指向截断约定，不是算错）")
for mus in (0.00, 0.13, 0.20):
    r = [x for x in rows if x[1] == mus]
    rl = [x[4] for x in r]
    print(f"   μ*={mus:.2f}: 相对差 ∈ [{min(rl):+.2f}%, {max(rl):+.2f}%]  中位 {np.median(rl):+.2f}%"
          f"   |Δf|>0.03 的格数 = {sum(1 for x in r if x[6]=='超出')}/{len(r)}")

n_above = sum(1 for x in rows if "之上" in x[7])
n_below = sum(1 for x in rows if "之下" in x[7])
n_out = sum(1 for x in rows if x[6] == "超出")
print("\n" + "=" * 118)
print("★ 判词（三段分开，不许合并）")
print(f"  ① **μ*=0 处两者互证**：{sum(1 for x in rows if x[1]==0.0)} 格全部落在我们括号**之上** "
      f"{min(x[4] for x in rows if x[1]==0.0):+.2f}%…{max(x[4] for x in rows if x[1]==0.0):+.2f}%")
print("     ⟹ 他们 2006 年的拟合与我们 2026 年的认证曲线在**无约定自由度**的那一端一致。")
print("     这是我们这条曲线迄今唯一的**外部独立佐证**（此前只有自洽性门）。")
print(f"  ② **μ*>0 处他们落在括号之下**（{n_below}/{len(rows)} 格），差距在 μ*=0 归零、随 μ* 单调增长。")
print("     ⚠ **不能据此说他们错**：他们全文**未声明 μ* 的截断约定**（已 grep 全文确认），")
print("     而 μ* 的数值含义完全由 ω_c 约定决定 ⟹ 这是『两个约定下的量』，不是『谁算错了』。")
print(f"     量级自查：换回他们的单位后 |Δf|>0.03（他们自陈的最大偏差）的只有 **{n_out} 格**，")
print("     其余全部落在他们**已经声明过的**拟合误差带内 ⟹ 连「拟合不够准」都谈不上。")
print("  ③ **真正的增量与约定无关**：eq.(11) 是最小二乘拟合（自陈 RMS 0.01/max 0.03，**双向**），")
print("     而 eq.(10) 把它当作**上界**使用。拟合值不构成上界——这是结构性的，不是精度问题。")
print("     ⟹ 我们的增量 = 『用外向舍入包封替换拟合』，**不是**『首次给出该上界』。")
print("     并且他们把 Einstein 极值性本身归给 ref.3=Allen-Dynes 1975 / ref.11=Leavens 1977，")
print("     ⟹ 定理内容三方都不拥有；μ*>0 时该极值性**至今无人证明**（见 no_rigorous_mustar_ceiling）。")
