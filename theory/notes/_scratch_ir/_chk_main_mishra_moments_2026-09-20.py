# -*- coding: utf-8 -*-
"""
主会话从零复核：Mishra α²F 的三个矩（λ, ω_log, ω̄₂），H3S / D3S × 谐波 / 非谐。

被复核的陈述（`materials/h3s_omega_log_discrepancy_resolved_2026-09-20.md`）：
  - 仓库 ω_log(H3S-200-anh) = 109.995646 meV
  - λ(H3S-200-anh) = 1.720958
  - ω̄₂(H3S-200-anh) = 149.26 meV,  ω̄₂(D3S-200-anh) = 105.43 meV
  - 该笔记 §3 行 100 声称这两个 ω̄₂「取自 hydrides_theory.csv:3,9」
    —— 但我实测 CSV 第 9 行（D3S anharmonic）的 omega_2_meV 列是 `not-reported`。
    本脚本一并把 CSV 那两列打出来，判定该出处陈述真伪。

定义（AD75）：
    λ      = 2 ∫ α²F(ω)/ω dω
    ω_log  = exp[ (2/λ) ∫ (α²F/ω) ln ω dω ]
    ω̄₂     = sqrt[ (2/λ) ∫ α²F(ω) ω dω ]        <-- 注意是 α²F·ω，不是 α²F/ω·ω²
                                                     两者相同；写成前者避免 0 处奇点

★ 两套独立求积方案（CLAUDE.md：多精度一致只测舍入，不测离散化）：
    (I)  复合 Simpson（均匀网格，文件自带 500 点）
    (II) 梯形 + Richardson（对半抽稀再外推）
  两者不一致到 1e-6 以上就报警。

严格性：numerical-evidence（float64 求积，两方案交叉）。
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

ROOT = r"D:\physics\superconductivity"
BASE = os.path.join(ROOT, "theory", "data", "external", "mishra_h3s_d3s",
                    "Data-H3S-D3S-PH-EPW-Anharmonic")
RY_TO_MEV = 13605.693122994   # 1 Ry = 13.605693122994 eV


def load_a2f(path):
    """EPW .a2f: 第一列频率(meV)，第二列 a2F，第三列是**运行积分** 2*a2F/w。
    第三列的端点 = 文件自报的 lambda，是一个免费的第三方锚点。
    尾部有 'Integrated el-ph coupling' / 'Summed el-ph coupling' 等非数值行。"""
    rows = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                continue
            if len(vals) >= 2:
                rows.append(vals)
    a = np.array([r[:3] if len(r) >= 3 else r[:2] + [np.nan] for r in rows], dtype=float)
    return a[:, 0], a[:, 1], a[:, 2]


def trap(y, x):
    return np.trapezoid(y, x) if hasattr(np, "trapezoid") else np.trapz(y, x)


def moments(w, a2f, quad):
    """w in meV。剔除 w<=0 与 a2f<0（EPW 尾部会有微小负值）。"""
    m = (w > 1e-9)
    w, a2f = w[m], np.clip(a2f[m], 0.0, None)
    lam = 2.0 * quad(a2f / w, w)
    wlog = np.exp((2.0 / lam) * quad(a2f / w * np.log(w), w))
    w2 = np.sqrt((2.0 / lam) * quad(a2f * w, w))
    return lam, wlog, w2


def richardson(y, x):
    """梯形在 h 与 2h 上的 Richardson 外推（梯形误差 O(h²)）。
    ⚠ 只对**准均匀**网格有效；网格若有跳点，两者都退化为梯形，差值不再是误差估计。"""
    I_h = trap(y, x)
    I_2h = trap(y[::2], x[::2])
    return (4.0 * I_h - I_2h) / 3.0


CASES = [
    ("H3S-200 anh", "H3S-200GPa/EPW_anharmonic/H3S.a2f"),
    ("H3S-200 har", "H3S-200GPa/EPW_harmonic/H3S.a2f"),
    ("D3S-200 anh", "D3S-200GPa/EPW_anharmonic/D3S.a2f"),
    ("D3S-200 har", "D3S-200GPa/EPW_harmonic/D3S.a2f"),
]

print("=" * 88)
print("主会话独立复核：Mishra α²F 的 λ / ω_log / ω̄₂（两套求积交叉）")
print("=" * 88)

res = {}
for tag, rel in CASES:
    p = os.path.join(BASE, rel.replace("/", os.sep))
    if not os.path.exists(p):
        print(f"  {tag:14s}  MISSING {p}")
        continue
    w_raw, a, runlam = load_a2f(p)
    unit = "meV"
    w = w_raw.copy()
    if w.max() < 1.0:
        w = w * RY_TO_MEV
        unit = "Ry->meV"
    # ★ 网格均匀性检查（H3S 谐波末点从 238.45 跳到 250.45 ⟹ 均匀-Simpson 失效）
    dw = np.diff(w)
    nonunif = dw.max() / dw.min() - 1.0
    lam_t, wl_t, w2_t = moments(w, a, trap)
    lam_r, wl_r, w2_r = moments(w, a, richardson)
    res[tag] = (lam_t, wl_t, w2_t)
    print(f"  {tag:14s} npts={len(w):4d} wmax={w.max():8.2f} meV [{unit}] "
          f"网格非均匀度={nonunif:.3f}" + ("  <-- 有跳点" if nonunif > 1e-6 else ""))
    print(f"      梯形       : lam={lam_t:.6f}  wlog={wl_t:.6f} meV  w2={w2_t:.6f} meV")
    print(f"      Richardson : lam={lam_r:.6f}  wlog={wl_r:.6f} meV  w2={w2_r:.6f} meV")
    print(f"      文件自报   : 运行积分端点={runlam[-1]:.7f}   （第三方锚点）")
    d = max(abs(lam_t - lam_r) / lam_t, abs(wl_t - wl_r) / wl_t, abs(w2_t - w2_r) / w2_t)
    e = abs(lam_t - runlam[-1]) / runlam[-1]
    print(f"      两方案最大相对差={d:.2e}   与文件自报 lam 的相对差={e:.2e}"
          + ("   OK" if (d < 1e-4 and e < 1e-3) else "   !! 检查"))

print()
print("-" * 88)
print("A. 对被复核笔记的逐条判决")
tgt = {"H3S-200 anh": (1.720958, 109.995646, 149.262922),
       "D3S-200 anh": (1.783761, 85.535326, 105.430104)}
for tag, (L, W, W2) in tgt.items():
    if tag not in res:
        continue
    l, wl, w2 = res[tag]
    print(f"  {tag}:")
    print(f"     lam   我={l:.6f}   笔记/记忆={L:.6f}   相对差={abs(l-L)/L:.2e}")
    print(f"     wlog  我={wl:.6f}  笔记/记忆={W:.6f}  相对差={abs(wl-W)/W:.2e}")
    print(f"     w2    我={w2:.6f}  笔记/记忆={W2:.6f}  相对差={abs(w2-W2)/W2:.2e}")

print()
print("-" * 88)
print("B. §3 的 1/sqrt(2) 伸缩检验（这是该笔记的物理结论，我自己重算）")
if "H3S-200 anh" in res and "D3S-200 anh" in res:
    _, wlH, w2H = res["H3S-200 anh"]
    _, wlD, w2D = res["D3S-200 anh"]
    inv2 = 1.0 / np.sqrt(2.0)
    print(f"   w2  : D/H = {w2D/w2H:.6f}   1/sqrt2 = {inv2:.6f}   偏离 {(w2D/w2H/inv2-1)*100:+.3f}%")
    print(f"   wlog: D/H = {wlD/wlH:.6f}   1/sqrt2 = {inv2:.6f}   偏离 {(wlD/wlH/inv2-1)*100:+.3f}%")
    print("   （笔记用的是 Errea Table II 的 wlog 比 0.7894=+11.6%；上一行是 Mishra 谱自己的比，")
    print("     两者是不同数据源，不可互相引用。）")
if "H3S-200 har" in res and "D3S-200 har" in res:
    _, wlHh, w2Hh = res["H3S-200 har"]
    _, wlDh, w2Dh = res["D3S-200 har"]
    inv2 = 1.0 / np.sqrt(2.0)
    print(f"   [谐波] w2  : D/H = {w2Dh/w2Hh:.6f}  偏离 {(w2Dh/w2Hh/inv2-1)*100:+.3f}%")
    print(f"   [谐波] wlog: D/H = {wlDh/wlHh:.6f}  偏离 {(wlDh/wlHh/inv2-1)*100:+.3f}%")
    print("   => 谐波行若也偏离，说明 α²F 的**形状**在 H/D 间就不同（不是纯质量伸缩），")
    print("      那么『非谐让 wlog 偏离 1/sqrt2』这个归因就不成立。这是判别性的一行。")

print()
print("-" * 88)
print("C. CSV 出处核验：笔记 §3 行 100 声称 w2 取自 hydrides_theory.csv:3,9")
csv = os.path.join(ROOT, "theory", "data", "hydrides_theory.csv")
with open(csv, "r", errors="replace") as fh:
    lines = fh.read().splitlines()
hdr = lines[0].split(",")
try:
    i_w2 = hdr.index("omega_2_meV")
except ValueError:
    i_w2 = None
print(f"   表头里 omega_2_meV 的列号 = {i_w2}")
for ln in (3, 9):
    if ln - 1 < len(lines):
        f = lines[ln - 1].split(",")
        val = f[i_w2] if (i_w2 is not None and i_w2 < len(f)) else "?"
        print(f"   第 {ln} 行: {f[0]:5s} {f[3] if len(f)>3 else '':12s} omega_2_meV = {val!r}")
print("   （若 D3S 行为 'not-reported'，则笔记 §3 的出处陈述为假：105.43 不在 CSV 里，")
print("     它来自 α²F 积分/主会话消息，出处必须改写。）")

print()
print("-" * 88)
print("D. 判别性检验：λ-累积分位数上的 D/H 频率比（不需要选带阈值）")
print("   机制假说：H3S->D3S 只替换 H，**S 模式完全不伸缩**、H 模式按 1/sqrt2 伸缩。")
print("   若为真 ⟹ 低分位（S 主导）比值≈1，高分位（H 主导）比值≈0.707。")
print("   若整条谱按 1/sqrt2 均匀伸缩（纯质量伸缩假设）⟹ 所有分位比值恒=0.707。")


def quantiles(w, a2f, qs):
    """按 dλ = 2(a2F/w)dω 的累积分布求分位频率。"""
    m = w > 1e-9
    w, a2f = w[m], np.clip(a2f[m], 0.0, None)
    dens = 2.0 * a2f / w
    c = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(w))])
    c = c / c[-1]
    return np.interp(qs, c, w)


QS = np.array([0.10, 0.25, 0.50, 0.75, 0.90])
CD = dict(CASES)
for lab, hk, dk in (("谐波", "H3S-200 har", "D3S-200 har"),
                    ("非谐", "H3S-200 anh", "D3S-200 anh")):
    wh, ah, _ = load_a2f(os.path.join(BASE, CD[hk].replace("/", os.sep)))
    wd, ad, _ = load_a2f(os.path.join(BASE, CD[dk].replace("/", os.sep)))
    qh, qd = quantiles(wh, ah, QS), quantiles(wd, ad, QS)
    print(f"   [{lab}] 比值 : " + "  ".join(f"q{int(q*100):02d}={d/h:.4f}"
                                            for q, h, d in zip(QS, qh, qd)))
    print(f"          H(meV): " + " ".join(f"{v:7.2f}" for v in qh))
    print(f"          D(meV): " + " ".join(f"{v:7.2f}" for v in qd))
print(f"   参照 1/sqrt2 = {1/np.sqrt(2):.4f}")

print("=" * 88)
