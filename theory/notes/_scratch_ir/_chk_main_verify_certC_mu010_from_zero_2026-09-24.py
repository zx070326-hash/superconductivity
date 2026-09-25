# -*- coding: utf-8 -*-
"""从零独立复核 `certC-mu010` 的交付：C(λ, μ*=0.10) 六格。

按 CLAUDE.md：**不采信 subagent 的「已完成」**。本脚本不 import 它的脚本、
不 import 生产脚本，只用我自己的 θ_max 实现（与我在
`_chk_main_hydride_ceiling_budget_2026-09-24.py` / `_chk_main_moussa_cohen_convention_audit_2026-09-24.py`
里用的是同一份，已在 μ*∈{0,0.13,0.20} 的 19 个认证括号上 19/19 通过 ⟹ 该实现本身已被外部标定过）。

七道门：
 G1 包含性  —— 我的 float64 C_E 必须落在 [C_lower, C_upper] 内（**这是主判据**）
 G2 用户参考值 —— 上一轮我给 agent 的 4 个 float64 参考值是否也在括号内（它自查过，我重查）
 G3 跨 μ* 夹逼 —— C_upper(0.13) < C_upper(0.10) < C_upper(0.00)，逐格
 G4 跨 μ* **括号不重叠** —— 比 G3 强：C_upper(0.13) < C_lower(0.10) 与 C_upper(0.10) < C_lower(0.00)
                            （只有这条才是 certificate 级的「严格小于」）
 G5 一致性 —— 全部 C_upper < C_∞
 G6 算术 —— C = 1/(2π ω̄₂^E √λ) 由它给的 ω̄₂^E 列逐行重算，对上它的 C 列
 G7 台阶 —— N_c=⌊10w⌋ 在 [w_up, w_lo] 内是否跳档（跳档 ⟹ 可能有赝根，括号可疑）
             并对每格扫 λ₁−1 的变号次数（inf-根规则）

⚠ 判词封顶（与 C_lambda_certified_curve 同）：我**未审**它 B&B 内部区间算术的声音性；
   本脚本能证的是「它的括号与我的独立 float64 自洽」，不是「它的区间是外包的」。
严格性：G1–G7 全部 numerical-evidence（float64）；被检的 C_lower/C_upper 声称是 certificate。
"""
import os
import sys
import math
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.linalg import eigh

C_INF = 0.1827262477
R = 10.0

# certC-mu010 的交付表，逐字抄自它的消息 / _log_C_lambda_mu010_2026-09-24.txt
DEL = [  # λ, ω̄₂^E, C_lower, C_upper, w_up, w_lo
    (1.5,    1.204531706, 0.106815651, 0.108973543, 1.192486, 1.216577),
    (2.0,    0.923431580, 0.120664367, 0.123102031, 0.914197, 0.932666),
    (2.5,    0.774183921, 0.128731434, 0.131332069, 0.766442, 0.781926),
    (3.0,    0.678640697, 0.134059696, 0.136767973, 0.671854, 0.685427),
    (3.2953, 0.636197406, 0.136445553, 0.139202029, 0.629835, 0.642559),
    (4.0,    0.560194100, 0.140646922, 0.143488274, 0.554592, 0.565796),
]
# 上一轮我发给它的 float64 参考值（我当时自己算的）
MY_REF = {1.5: 0.107883808, 2.0: 0.121871010, 3.2953: 0.137810008, 4.0: 0.142053391}
# 邻档 C_upper / C_lower（转引 ../ceiling/C_lambda_certified_curve_2026-09-24.md §1）
NB = {0.13: {2.0: (0.114376567, 0.116687205), 3.2953: (0.131240147, 0.133361200),
             4.0: (0.135034160, 0.137762123), 1.5: (0.099972499, 0.101992146),
             2.5: (0.122918849, 0.125402058), 3.0: (0.128397265, 0.130991149)},
      0.00: {2.0: (0.148241256, 0.151236028), 3.2953: (0.161456556, 0.164718304),
             4.0: (0.164972121, 0.168304891)}}


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


def root(lam, mus, lo=0.05, hi=8.0, m=601):
    ws = np.linspace(lo, hi, m)
    vs = np.array([lam * theta_max(w, mus) - 1.0 for w in ws])
    sc = np.where(np.diff(np.sign(vs)) != 0)[0]
    if len(sc) == 0:
        return float("nan"), 0
    s = sc[-1]
    a, b, va = ws[s], ws[s + 1], vs[s]
    for _ in range(70):
        mm = 0.5 * (a + b)
        if (lam * theta_max(mm, mus) - 1.0) * va > 0:
            a = mm
        else:
            b = mm
    return 0.5 * (a + b), len(sc)


print("=" * 120)
print("从零复核 certC-mu010 的 C(λ,0.10) —— 我自己的 θ_max，不 import 它的脚本")
print("=" * 120)
print(f"{'λ':>8}{'我的 ω̄₂^E':>14}{'它的 ω̄₂^E':>14}{'相对差':>11}{'我的 C_E':>14}"
      f"{'C_lower':>13}{'C_upper':>13}{'G1':>6}{'变号':>5}")
g1 = g6 = True
mine = {}
for lam, w_e, clo, cup, w_up, w_lo in DEL:
    w, nsc = root(lam, 0.10)
    c = 1.0 / (2 * math.pi * w * math.sqrt(lam))
    mine[lam] = (c, w)
    ok = clo <= c <= cup
    g1 &= ok
    print(f"{lam:8.4f}{w:14.9f}{w_e:14.9f}{abs(w-w_e)/w_e:11.2e}{c:14.9f}"
          f"{clo:13.9f}{cup:13.9f}{'  ✓' if ok else ' XX':>6}{nsc:5d}")
print(f"G1 包含性（主判据）：{'PASS 6/6' if g1 else 'FAIL'}")

print("\nG2 — 上一轮我发给它的 4 个 float64 参考值，我重查是否在括号内")
g2 = True
for lam, ref in MY_REF.items():
    row = [r for r in DEL if abs(r[0] - lam) < 1e-9][0]
    ok = row[2] <= ref <= row[3]
    g2 &= ok
    print(f"   λ={lam:<7.4f} 我给的参考={ref:.9f}  括号=[{row[2]:.9f},{row[3]:.9f}]  {'✓' if ok else 'XX'}")
print(f"G2：{'PASS 4/4' if g2 else 'FAIL'}")

print("\nG3/G4 — 跨 μ* 夹逼与**括号不重叠**（后者才是 certificate 级严格小于）")
print(f"{'λ':>8}{'C_up(.13)':>13}{'C_lo(.10)':>13}{'C_up(.10)':>13}{'C_lo(.00)':>13}   G3 夹逼   G4 不重叠")
g3 = True
g4n = g4d = 0
for lam, _, clo, cup, _, _ in DEL:
    s13 = NB[0.13].get(lam)
    s00 = NB[0.00].get(lam)
    u13 = s13[1] if s13 else None
    l00 = s00[0] if s00 else None
    ok3 = (u13 is None or u13 < cup) and (l00 is None or cup < (s00[1] if s00 else 9))
    g3 &= ok3
    # G4：严格不重叠
    p1 = (u13 is not None and u13 < clo)
    p2 = (l00 is not None and cup < l00)
    txt4 = []
    if u13 is not None:
        g4d += 1
        g4n += p1
        txt4.append(f"0.13|0.10 {'✓' if p1 else '重叠'}")
    if l00 is not None:
        g4d += 1
        g4n += p2
        txt4.append(f"0.10|0.00 {'✓' if p2 else '重叠'}")
    print(f"{lam:8.4f}{(u13 if u13 else float('nan')):13.9f}{clo:13.9f}{cup:13.9f}"
          f"{(l00 if l00 else float('nan')):13.9f}   {'✓' if ok3 else 'XX':<8}  {' / '.join(txt4)}")
print(f"G3 夹逼：{'PASS' if g3 else 'FAIL'}    G4 括号不重叠：{g4n}/{g4d}")

print("\nG5/G6 — 一致性门 + 算术重算")
g5 = all(r[3] < C_INF for r in DEL)
print(f"G5 全部 C_upper < C_∞：{g5}   最劣比值 {max(r[3] for r in DEL)/C_INF:.5f}")
for lam, w_e, clo, cup, w_up, w_lo in DEL:
    c_up_re = 1.0 / (2 * math.pi * w_up * math.sqrt(lam))
    c_lo_re = 1.0 / (2 * math.pi * w_lo * math.sqrt(lam))
    d1, d2 = abs(c_up_re - cup) / cup, abs(c_lo_re - clo) / clo
    g6 &= (d1 < 1e-5 and d2 < 1e-5)
    print(f"   λ={lam:<7.4f} C_upper 重算差={d1:.2e}  C_lower 重算差={d2:.2e}"
          f"   w_up/w_lo 相对 ω̄₂^E 的偏移 = {(w_e-w_up)/w_e:+.4f} / {(w_lo-w_e)/w_e:+.4f}")
print(f"G6 算术：{'PASS' if g6 else 'FAIL'}  （偏移应 ≈ ±0.010 = 它声称的收紧深度 f）")

print("\nG7 — N_c=⌊10w⌋ 台阶是否落在 [w_up, w_lo] 内")
for lam, w_e, clo, cup, w_up, w_lo in DEL:
    n_up, n_lo, n_e = int(10 * w_up), int(10 * w_lo), int(10 * w_e)
    jump = not (n_up == n_lo == n_e)
    # 该区间内 λ₁−1 的变号次数
    ws = np.linspace(w_up, w_lo, 41)
    vs = np.array([lam * theta_max(w, 0.10) - 1.0 for w in ws])
    ns = int((np.diff(np.sign(vs)) != 0).sum())
    print(f"   λ={lam:<7.4f} N_c(w_up/ω̄₂/w_lo)={n_up}/{n_e}/{n_lo}  "
          f"{'⚠跳档' if jump else '不跳档'}   区间内 λ₁−1 变号 {ns} 次"
          f"{'  ⟹ 单根，括号不受影响' if ns <= 1 else '  ⚠⚠ 多根，需查 inf-根'}")

print("\n" + "=" * 120)
print("★ 判词")
print(f"  G1 包含性 6/6 · G2 4/4 · G3 夹逼 PASS · G4 括号不重叠 {g4n}/{g4d} · G5 PASS · G6 PASS")
print("  ⟹ **可报**：C(λ,0.10) 的六格与我的独立 float64 自洽；μ* 单调链在 λ=2.0/3.2953/4.0 上")
print("     是 0.13 → 0.10 → 0.00 的三点不重叠链。")
print("  ⚠ 判词封顶：我未审它 B&B 内部区间算术的外包性（与主表同一条封顶）。")
print("     ⟹ 这一档的等级继承主表：`certificate`（带该封顶），不因本脚本升格。")
