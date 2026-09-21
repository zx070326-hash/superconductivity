# -*- coding: utf-8 -*-
"""
同位素指数的精确对数可加分解 alpha = alpha_omega + alpha_lambda，
用于 theory/notes/isotope/isotope_pdh_discriminating_test_2026-09-18.md 的二次订正（2026-09-20）。

恒等式（analytic）：任何 Tc = omega_log * F(lambda, mu*) 形式的模型，
  ln Tc = ln omega_log + ln F  =>  alpha := -dlnTc/dlnM = alpha_omega + alpha_lambda，无交叉项。
其中 alpha_omega = -Dln(omega_log)/Dln(M) 不含模型，是纯谱学量。

数据：Errea et al., arXiv:1305.7123, Table I（谐波列 lambda == 1.55 三同位素全同，
作者正文亦明言 "In the harmonic approximation lambda is independent of the mass"）。

跑： OMP_NUM_THREADS=1 python _chk_iso_alpha_decomp_2026-09-20.py
"""
import math
from itertools import product

# CODATA 原子质量 (u)
M = {"H": 1.007825, "D": 2.014102, "T": 3.016049}

# Errea Table I: (lambda, omega_log[K], Tc_ME[K])
HARM = {"H": (1.55, 205.0, 47.0), "D": (1.55, 150.0, 34.0), "T": (1.55, 125.0, 30.0)}
SSCHA = {"H": (0.40, 405.0, 5.0), "D": (0.46, 304.0, 6.5), "T": (0.48, 257.0, 6.9)}

MUSTAR = 0.085  # 本笔记全程采用的值（Errea 亦用 0.085）


def tc_ad(lam, wlog, mus):
    """Allen-Dynes/McMillan 闭式 (K)。"""
    den = lam - mus * (1.0 + 0.62 * lam)
    if den <= 0:
        return float("nan")
    return (wlog / 1.2) * math.exp(-1.04 * (1.0 + lam) / den)


def alpha(x1, x2, m1, m2):
    """alpha = -ln(x2/x1)/ln(m2/m1)，x 为 Tc 或 omega_log。"""
    return -math.log(x2 / x1) / math.log(m2 / m1)


def row(tab, a, b, label):
    lam1, w1, tme1 = tab[a]
    lam2, w2, tme2 = tab[b]
    m1, m2 = M[a], M[b]
    a_w = alpha(w1, w2, m1, m2)                       # 纯谱学，无模型
    F1 = tc_ad(lam1, w1, MUSTAR) / w1                 # F = Tc_AD / omega_log
    F2 = tc_ad(lam2, w2, MUSTAR) / w2
    a_l = alpha(F1, F2, m1, m2)                       # lambda 通道
    a_ad = alpha(tc_ad(lam1, w1, MUSTAR), tc_ad(lam2, w2, MUSTAR), m1, m2)
    a_me = alpha(tme1, tme2, m1, m2)                  # Errea 全 ME 解
    print(f"{label:14s} a_w={a_w:+.4f}  a_lam={a_l:+.4f}  a_AD={a_ad:+.4f} "
          f"(可加性残差 {a_ad - (a_w + a_l):+.2e})  a_ME={a_me:+.4f}  a_ME-a_w={a_me - a_w:+.4f}")
    return a_w, a_l, a_ad, a_me


print("=== 精确分解 alpha = alpha_omega + alpha_lambda ===")
print("(谐波行 lambda 三同位素全同 => alpha_lambda 应恒为 0)")
row(HARM, "H", "D", "谐波 H->D")
hw, hl, had, hme = row(HARM, "D", "T", "谐波 D->T")
row(SSCHA, "H", "D", "SSCHA H->D")
row(SSCHA, "D", "T", "SSCHA D->T")

print("\n=== dln(lambda)/dln(M)（非谐 lambda 增强的质量依赖） ===")
for a, b in (("H", "D"), ("D", "T")):
    d = math.log(SSCHA[b][0] / SSCHA[a][0]) / math.log(M[b] / M[a])
    print(f"  SSCHA {a}->{b}: {d:+.4f}")

print("\n=== 谐波 D->T 残差的舍入极值扫描 ===")
print("Errea 印出 Tc_D=34, Tc_T=30（两位有效数字）=> 真值分别落在 [33.5,34.5) 与 [29.5,30.5)")
best = (None, float("inf"))
worst = (None, 0.0)
N = 201
for i, j in product(range(N), range(N)):
    td = 33.5 + 1.0 * i / (N - 1) * 0.999999
    tt = 29.5 + 1.0 * j / (N - 1) * 0.999999
    dev = alpha(td, tt, M["D"], M["T"]) - hw          # ME 侧 - 纯谱学侧
    if abs(dev) < abs(best[1]):
        best = ((td, tt), dev)
    if abs(dev) > abs(worst[1]):
        worst = ((td, tt), dev)
print(f"  alpha_omega（= 谐波行 AD 的全部）      = {hw:+.4f}")
print(f"  中心值 (34,30) 的 ME 偏离              = {hme - hw:+.4f}")
print(f"  偏离最小（最有利于'全是舍入'的解释）: {best[1]:+.4f}  at Tc_D={best[0][0]:.2f} Tc_T={best[0][1]:.2f}")
print(f"  偏离最大                            : {worst[1]:+.4f}  at Tc_D={worst[0][0]:.2f} Tc_T={worst[0][1]:.2f}")
print("  => 舍入吸收到极限后残差仍 >= 0.0638，且全区间同号（负）。")
print("     方向审查：取的是'最有利于对立解释'的一侧 => 只支持'残差存在'这一 PASS 方向。")
