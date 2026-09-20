# -*- coding: utf-8 -*-
"""_lead_verify_ad_proxy_bias_2026-09-20.py — 主会话从零复核 deficit-real-solve 的三行 AD 代理偏差。

独立于 agent 脚本的部分：
  · 自己的 a2f 解析 + 矩积分（并用文件自带的累计 λ 列做交叉核）
  · 自己的 AD75 Eq.(34)/(35)（f2 无开方）
  · 闭式：μ*=0 时 δ_BRD^AD = 1 − f2(r,λ)/r，r = ω̄₂/ω_log（f1 与指数因子约掉）⟹ 两条路算同一个数
  · 真解 δ_BRD 取 deficit_budget_h3s_2026-09-17.md:36/38/40（已两路径复核），不转引 agent 表
判词由值决定，不预写。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np

ROOT = "theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic"
ROWS = [  # name, path, lam_ref, w2_ref(meV), dBRD_true (deficit_budget_h3s_2026-09-17.md)
    ("H3S-200 anharm", f"{ROOT}/H3S-200GPa/EPW_anharmonic/H3S.a2f", 1.721, 149.26, 0.135),
    ("D3S-200 anharm", f"{ROOT}/D3S-200GPa/EPW_anharmonic/D3S.a2f", 1.784, 105.43, 0.098),
    ("H3S-160 anharm", f"{ROOT}/H3S-160GPa/EPW_anharmonic/H3S.a2f", 1.876, 137.87, 0.100),
]


def read_a2f(path):
    w, a, cum = [], [], []
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 3:
                continue
            try:
                x, y, z = float(p[0]), float(p[1]), float(p[2])
            except ValueError:
                continue
            if x <= 0:
                continue
            w.append(x); a.append(y); cum.append(z)
    return np.array(w), np.array(a), np.array(cum)


def moments(w, a):
    # 网格等距、从 0 起：用梯形，且把 (0,0) 原点补进去（a2F(0)=0）
    w0 = np.concatenate([[0.0], w]); a0 = np.concatenate([[0.0], a])
    g = np.where(w0 > 0, a0 / np.where(w0 > 0, w0, 1.0), 0.0)
    lam = 2 * np.trapezoid(g, w0)
    wlog = np.exp(2 / lam * np.trapezoid(np.where(w0 > 0, g * np.log(np.where(w0 > 0, w0, 1.0)), 0.0), w0))
    w2 = np.sqrt(2 / lam * np.trapezoid(a0 * w0, w0))
    return lam, wlog, w2


def f2_ad(lam, r, mu):
    L2 = 1.82 * (1 + 6.3 * mu) * r
    return 1 + (r - 1) * lam**2 / (lam**2 + L2**2)


def tc_ad(lam, wlog, w2, mu):
    r = w2 / wlog
    L1 = 2.46 * (1 + 3.8 * mu)
    f1 = (1 + (lam / L1) ** 1.5) ** (1 / 3)
    return wlog / 1.2 * f1 * f2_ad(lam, r, mu) * np.exp(-1.04 * (1 + lam) / (lam - mu * (1 + 0.62 * lam)))


print(f"{'row':16} {'lam':>7} {'lam_file':>8} {'w2':>8} {'wlog':>8} {'r':>6} {'dBRD_AD(ratio)':>15} {'dBRD_AD(closed)':>16} {'dBRD_true':>9} {'bias pp':>8}")
bias = []
for name, path, lam_ref, w2_ref, d_true in ROWS:
    w, a, cum = read_a2f(path)
    lam, wlog, w2 = moments(w, a)
    r = w2 / wlog
    d_ratio = 1 - tc_ad(lam, wlog, w2, 0.0) / tc_ad(lam, w2, w2, 0.0)
    d_closed = 1 - f2_ad(lam, r, 0.0) / r
    assert abs(d_ratio - d_closed) < 1e-12, (d_ratio, d_closed)
    assert abs(lam - lam_ref) < 5e-4 * lam_ref and abs(w2 - w2_ref) < 5e-3 * w2_ref, (lam, w2)
    b = (d_closed - d_true) * 100
    bias.append(b)
    print(f"{name:16} {lam:7.4f} {cum[-1]:8.4f} {w2:8.2f} {wlog:8.2f} {r:6.4f} {d_ratio*100:14.2f}% {d_closed*100:15.2f}% {d_true*100:8.1f}% {b:+8.2f}")

bias = np.array(bias)
print(f"\nmean bias = {bias.mean():+.2f} pp   std(pop) = {bias.std():.2f}   std(sample) = {bias.std(ddof=1):.2f}   same sign: {int((bias > 0).sum())}/3")
print(f"AD75 29-row proxy mean 19.7% shifted by mean bias -> {19.7 - bias.mean():.1f}%")
# 真解 δ_BRD 的分辨力：参考表印到 0.1 pp ⟹ 偏差每行 ±0.05 pp 的印数误差，可忽略。
