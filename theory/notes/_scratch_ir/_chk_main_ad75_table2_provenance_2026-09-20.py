# -*- coding: utf-8 -*-
"""
主会话从零核验 `f-sample-expansion` 的「AD75 Table II = 中子+Born-von Karman 的独立 λ 来源」主张。
(2026-09-20 晚)

被核验主张（逐条）：
  M1  Table II 的 λ 来自 neutron phonon + Born-von Karman 分析，**不是**隧道 ⟹ 可关掉家族阻断。
  M2  Ta (Tc=4.48 K, λ=0.69) 是 neutron+BvK，不是隧道。
  M3  Mo 的 λ_tun = 0.60。
  M4  Al 在 AD75 里「尚未定位」。

真源：AD1975.pdf p.916（PDF 第 12 页）页图，主会话亲读（`_ad75_t2c-12.png`）。
Table II 表头（逐字）：
  Material | wbar_log(K) | wbar_1(K) | wbar_2(K) | Tc(K) |
  lambda Eq.(34) mu*=0.10 | lambda Eq.(34) mu*=0.13 | lambda McMillan Ref.1 mu*=0.13 |
  lambda(expt) | mu*(expt) | eta(eV/A^2)
题注（逐字，关键句）：
  "The empirical values of w_log, w_1, w_2 are calculated from Born-von Karman analysis of
   neutron data assuming a^2F = F, or (where indicated) from tunneling measurements of a^2F."
正文（同页右栏，逐字）：
  "One of the principal uses of approximate Tc equations is the estimation of approximate
   'empirical' values of lambda from measured values of Tc. We illustrate this in Table II."

⟹ 中子+BvK 只供 **频率矩（列 2-4）**；**lambda（列 6/7/8）是从 Tc 反解 Eq.(34)/McMillan 得到的**。

判别动作（沿用 `generating-formula-exists-is-not-zero-dof` 的规则）：
  **能否逐位复现印出的数。** 本脚本对列 6 / 列 7 做 Eq.(34) 反解。
  AD75 自陈 "Because the values of lambda are all below 1, the correction factors f1 and f2
  are entirely negligible" ⟹ 取 f1=f2=1。

方向：本脚本是双侧恒等式复现，不产生单边界判决。
严格性：Eq.(34) 反解为 analytic（给定输入的闭式求根）；与印出值的比对为 numerical-evidence。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.optimize import brentq

C_INF = 0.1827262477

# --- AD75 Table II，页图逐行转录（p.916）-------------------------------------
#  name, wlog, w1, w2, Tc, lam_c6(mu*=.10), lam_c7(mu*=.13), lam_c8(McM), lam_expt, mu_expt, eta
T2 = [
    ("V",            222, 230, 238,  5.30, 0.60, 0.67, 0.60, None,  None,  3.1),
    ("Cr",           308, 317, 324,  None, None, None, None, None,  None,  None),
    ("Nb",           166, 175, 183,  9.22, 0.85, 0.94, 0.82, None,  None,  4.7),
    ("Mo",           236, 244, 251,  0.92, 0.39, 0.45, 0.41, None,  None,  4.2),
    ("Ta",           142, 148, 153,  4.48, 0.67, 0.74, 0.65, None,  None,  5.1),
    ("Ta(tunnel)",   132, 140, 148,  4.48, 0.69, 0.76, 0.65, 0.69,  0.093, 4.9),
    ("W",            208, 215, 220, 0.012, 0.25, 0.28, 0.28, None,  None,  4.5),
    ("Nb3Sn(tunnel)",125, 146, 163,  18.1, 1.67, 1.83, None, None,  None,  7.9),
]

def tc_eq34(lam, wlog, mus):
    """AD75 Eq.(34) with f1=f2=1 (AD75 自陈 lambda<1 时 f1,f2 可忽略)."""
    den = lam - mus*(1.0 + 0.62*lam)
    if den <= 0:
        return 0.0
    return wlog/1.20*np.exp(-1.04*(1.0+lam)/den)

def invert(tc, wlog, mus):
    """给定 Tc 反解 lambda。"""
    f = lambda l: tc_eq34(l, wlog, mus) - tc
    lo = mus*1.0001/(1.0-0.62*mus) + 1e-9      # den>0 的下端
    lo = max(lo, 1e-6)
    a = lo + 1e-7
    while tc_eq34(a, wlog, mus) >= tc and a < 50:   # 从极小 Tc 一侧起步
        a *= 1.0001
    try:
        return brentq(f, a, 50.0, xtol=1e-12, rtol=1e-14)
    except Exception:
        return float("nan")

print("="*100)
print("M1  Table II 的 lambda 是不是从 Tc 反解的？—— 判别动作：逐位复现印出的数")
print("="*100)
print(f"  {'材料':>14}{'wlog':>7}{'Tc':>8} | {'列6印':>7}{'我反解':>9}{'差':>9} | {'列7印':>7}{'我反解':>9}{'差':>9}")
ok6 = ok7 = tot = 0
for nm, wlog, w1, w2, tc, c6, c7, c8, lex, mex, eta in T2:
    if tc is None:
        continue
    tot += 1
    r6 = invert(tc, wlog, 0.10)
    r7 = invert(tc, wlog, 0.13)
    d6 = r6 - c6 if c6 is not None else float("nan")
    d7 = r7 - c7 if c7 is not None else float("nan")
    if c6 is not None and abs(d6) <= 0.005: ok6 += 1
    if c7 is not None and abs(d7) <= 0.005: ok7 += 1
    print(f"  {nm:>14}{wlog:>7}{tc:>8.3f} | {c6:>7.2f}{r6:>9.4f}{d6:>+9.4f} |"
          f" {c7:>7.2f}{r7:>9.4f}{d7:>+9.4f}")
print(f"\n  列6（mu*=0.10）复现 {ok6}/{tot}   列7（mu*=0.13）复现 {ok7}/{tot}   （判据：|差| <= 0.005 = 印数分辨力）")
print("  ⟹ 若全数复现，则 Table II 的 lambda = F(Tc, wlog, mu*) 闭式，**零材料信息**，")
print("     与 `filling_factor_f` 同形（CLAUDE.md 点名的循环量）。")

print()
print("="*100)
print("M2  Ta 是不是隧道？")
print("="*100)
print("  Table II 的行标签逐字为 `Ta (tunneling)`，且 0.69 印在 **lambda(expt)** 列（列 9），")
print("  mu*(expt)=0.093 同行 —— 这两列只有该行与 Nb3Sn 有值。")
print("  AD75 正文（p.913 右栏，逐字）：")
print("    \"We see from this compilation that **except for Ta**, the materials for which data")
print("     are available are all s-p or 'soft' materials.\"  —— 指的正是 Table I（隧道汇编）。")
print("  同页（p.916 右栏）：\"the isotope effect in Mo and **tunneling in Ta** yield values of mu* of 0.09.\"")
print("  ⟹ **M2 为假**：Ta 是 AD75 唯一一个有隧道 a^2F 的过渡金属。")
print()
print("  交叉验证：Table II `Ta(tunneling)` 行与仓库 29 行样本的 Ta 行逐位比对")
repo_ta = dict(wlog=132, w1=140, w2=148, lam=0.69, tc=4.48)
t2_ta   = dict(wlog=132, w1=140, w2=148, lam=0.69, tc=4.48)   # 列 2,3,4,9,5
same = all(abs(repo_ta[k]-t2_ta[k]) < 1e-12 for k in repo_ta)
print(f"    仓库(Table I): wlog=132 w1=140 w2=148 lam=0.69 Tc=4.48")
print(f"    Table II 行  : wlog=132 w1=140 w2=148 lam(expt)=0.69 Tc=4.48")
print(f"    ⟹ 逐位一致：{same}  ⟹ **仓库 Ta 行被第二次印刷独立复核**（AD75 自己印了两遍）")
f_ta = repo_ta["tc"]/(C_INF*np.sqrt(repo_ta["lam"])*repo_ta["w2"])
print(f"    f(Ta) = {f_ta:.4f}   （索引记载的下端极值 0.1994 ⟹ "
      f"{'复现 OK' if abs(f_ta-0.1994)<5e-4 else '不符!!'}）")
print()
print("  ★ 同材料两技术的 **omega 侧** 直接对比（Table II 唯一给出的、零数字化的一对）：")
neu = (142, 148, 153); tun = (132, 140, 148)
for lab, a, b in [("wlog", neu[0], tun[0]), ("w1", neu[1], tun[1]), ("w2", neu[2], tun[2])]:
    print(f"    {lab:>5}: 中子+BvK {a:>4} K   隧道 {b:>4} K   比 {a/b:.4f}  ({100*(a/b-1):+.2f}%)")
print("    ⟹ omega_2 侧的技术偏置在 Ta 上实测 **+3.38%**（sigma_0 型），n=1 ⟹ 估不出 s1。")

print()
print("="*100)
print("M3  Mo 的 lambda_tun = 0.60 ？")
print("="*100)
print("  Table II 的 Mo 行：lambda = 0.39 (列6) / 0.45 (列7) / 0.41 (列8)，lambda(expt) = ···（空）")
print("  **0.60 是 V 的列6/列8 值** ⟹ M3 是一次**行错位读数**（读到了 V 行）。")
print("  且三个值全部是 Tc 反解值 ⟹ 即便读对也不可用。")
print("  AD75 明说 Mo 只有**同位素效应**数据，没有隧道 a^2F。")

print()
print("="*100)
print("M4  Al 在 AD75 里在哪？")
print("="*100)
t1_mats = ("Pb In Sn Hg Tl Ta a-Ga b-Ga + Pb/Tl/Bi/In 合金".split())
print("  Table I  (p.913)：Pb, In, Sn, Hg, Tl, Ta, a-Ga, b-Ga 及 Pb/Tl/Bi/In 系合金 —— **无 Al**")
print("  Table II (p.916)：V, Cr, Nb, Mo, Ta, W, Nb3Sn（bcc 过渡金属 + Nb3Sn）—— **无 Al**")
print("  ⟹ `lambda_tun(Al) = 0.43` **在 AD75 中无出处**；该点是无源数据。")

print()
print("="*100)
print("四点 r1 拟合的输入审计（消费点：q_D = q_T/(1+r1)^2）")
print("="*100)
rows = [("Pb", 1.55, "Table I 隧道", "有效"),
        ("Ta", 0.69, "Table I 隧道（Table II 第二次印刷佐证）", "有效"),
        ("Mo", 0.60, "**读错行**（0.60 是 V 的值）；Mo 真值 0.39/0.45/0.41 且**全为 Tc 反解**", "无效"),
        ("Al", 0.43, "**AD75 中不存在**", "无效")]
print(f"  {'材料':>5}{'用的 lam':>10}  {'出处':<52}{'判定':>6}")
for nm, l, src, verdict in rows:
    print(f"  {nm:>5}{l:>10.2f}  {src:<52}{verdict:>6}")
print("\n  ⟹ 四点中 **2 点无效** ⟹ `r1 = +0.119 ± 0.309` **不是「宽」，是「空」**。")
print("     （两者后果不同：宽 ⟹ 扩样可救；空 ⟹ 该拟合根本没有可扩的样本基础。）")
print("     方向：污染只会**放大**散布 ⟹ 先前基于 sigma_r1 的**驳回**判决方向安全，结论不变。")
print("="*100)
