# -*- coding: utf-8 -*-
"""从零重算 Mishra a2f 的 λ 与 ω̄₂（D3S 与 H3S，200 GPa，谐波/非谐）。

目的：`omega2-data` 交付了 D3S 的 ω̄₂ = 105.430104（非谐）/ 92.973347（谐波），
并声称「已由主会话独立复核到 1.3e-7」。按 CLAUDE.md，**不采信「已完成」**——
本脚本不 import 任何已有脚本、不读任何笔记里的数，只从 `.a2f` 原始文件重算。

定义（与天花板消费点一致）：
  λ    = 2 ∫ α²F(ω)/ω dω
  ⟨ω²⟩ = (2/λ) ∫ ω α²F(ω) dω          ⟹  ω̄₂ = √⟨ω²⟩
  M2   = λ⟨ω²⟩ = 2 ∫ ω α²F(ω) dω       （= η/M，Hopfield 求和规则；**不经过 λ**）
  ω_log= exp( (2/λ) ∫ ln(ω) α²F(ω)/ω dω )

交叉验证（两套独立方案，不是多精度——多精度只测舍入）：
  (i) 文件自带第 3 列 = 累积 ∫2α²F/ω ⟹ 其末值就是文件作者自己算的 λ；
      再与文件尾部 "Summed el-ph coupling" 对照 ⟹ **三方对撞**。
  (ii) 我自己的梯形积分 vs Simpson 积分（离散化误差，而非舍入）。

⚠ 网格从 ω=0.337 meV 起、等距，末端到 ~170 meV。ω→0 处 α²F/ω 有限（α²F~ω²）。
严格性：numerical-evidence（float64 + 有限网格）。
"""
import os
import sys
import glob
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.integrate import simpson

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                    "data", "external", "mishra_h3s_d3s",
                    "Data-H3S-D3S-PH-EPW-Anharmonic")


def load(path):
    """返回 (w[meV], a2f, cum3, lam_file)。只取前 3 列的数值行。"""
    w, a, c = [], [], []
    lam_file = None
    for ln in open(path, encoding="utf-8", errors="replace"):
        if "Summed el-ph coupling" in ln:
            lam_file = float(ln.split()[-1])
            continue
        p = ln.split()
        if len(p) < 3:
            continue
        try:
            v = [float(x) for x in p[:3]]
        except ValueError:
            continue
        w.append(v[0]); a.append(v[1]); c.append(v[2])
    return np.array(w), np.array(a), np.array(c), lam_file


def moments(w, a):
    """两套积分方案（梯形 / Simpson），返回 dict。"""
    out = {}
    for nm, I in (("trapz", lambda y: np.trapezoid(y, w)), ("simps", lambda y: simpson(y, x=w))):
        lam = 2.0 * I(a / w)
        M2 = 2.0 * I(w * a)                       # = λ⟨ω²⟩ = η/M，**不经过 λ**
        w2 = np.sqrt(M2 / lam)
        wlog = np.exp((2.0 / lam) * I(np.log(w) * a / w))
        out[nm] = dict(lam=lam, M2=M2, w2=w2, wlog=wlog)
    return out


print("=" * 116)
print("从零重算 Mishra .a2f 的 λ / ω̄₂ / M2 / ω_log —— 不 import 任何已有脚本，不读任何笔记里的数")
print("=" * 116)

CASES = [("H3S", "H3S-200GPa/EPW_harmonic/H3S.a2f", "harmonic"),
         ("H3S", "H3S-200GPa/EPW_anharmonic/H3S.a2f", "anharmonic"),
         ("D3S", "D3S-200GPa/EPW_harmonic/D3S.a2f", "harmonic"),
         ("D3S", "D3S-200GPa/EPW_anharmonic/D3S.a2f", "anharmonic")]

# 待核对的外部声明值（来源：CSV 已有列 / omega2-data 的交付消息）
CLAIM = {("H3S", "harmonic"): (2.2929, 107.165, 129.32, "CSV 第2行（已入库）"),
         ("H3S", "anharmonic"): (1.721, 109.995, 149.26, "CSV 第3行（已入库）"),
         ("D3S", "harmonic"): (None, None, 92.973347, "omega2-data 交付（未入库）"),
         ("D3S", "anharmonic"): (None, None, 105.430104, "omega2-data 交付（未入库）")}

res = {}
for mat, rel, treat in CASES:
    p = os.path.join(BASE, rel)
    w, a, c, lam_file = load(p)
    m = moments(w, a)
    res[(mat, treat)] = m
    t, s = m["trapz"], m["simps"]
    print(f"\n── {mat} {treat} @200 GPa  ({len(w)} 个网格点, ω∈[{w[0]:.4f},{w[-1]:.2f}] meV)")
    print(f"   λ      trapz={t['lam']:.7f}  simps={s['lam']:.7f}   离散化差={abs(t['lam']-s['lam']):.2e}")
    print(f"          文件第3列末值={c[-1]:.7f}   文件自陈 'Summed el-ph coupling'={lam_file:.7f}")
    print(f"          ⟹ 三方对撞：|我的 trapz − 文件累积| = {abs(t['lam']-c[-1]):.2e}, "
          f"|我的 − 自陈| = {abs(t['lam']-lam_file):.2e}")
    print(f"   ω̄₂     trapz={t['w2']:.6f}  simps={s['w2']:.6f} meV   离散化差={abs(t['w2']-s['w2']):.2e}")
    print(f"   M2     trapz={t['M2']:.3f}  simps={s['M2']:.3f} meV²   （=η/M，不经过 λ）")
    print(f"   ω_log  trapz={t['wlog']:.6f}  simps={s['wlog']:.6f} meV")
    cl_lam, cl_wlog, cl_w2, src = CLAIM[(mat, treat)]
    print(f"   ↳ 对照 {src}:")
    if cl_lam is not None:
        print(f"       λ 声明={cl_lam:.4f}      我的={t['lam']:.4f}     相对差={abs(t['lam']-cl_lam)/cl_lam*100:+.4f}%")
        print(f"       ω_log 声明={cl_wlog:.3f}  我的={t['wlog']:.3f}   相对差={abs(t['wlog']-cl_wlog)/cl_wlog*100:+.4f}%")
    print(f"       ω̄₂ 声明={cl_w2:.6f}   我的={t['w2']:.6f}   相对差={abs(t['w2']-cl_w2)/cl_w2*100:+.5f}%")

# ── 同位素自查：λ 应同位素无关；ω̄₂ 应 ~1/√2 ──
print("\n" + "=" * 116)
print("同位素自查（这两条是物理必然，任何一条破了就说明我读错了文件）")
for treat in ("harmonic", "anharmonic"):
    h, d = res[("H3S", treat)]["trapz"], res[("D3S", treat)]["trapz"]
    print(f"  {treat:<11} λ_H={h['lam']:.5f}  λ_D={d['lam']:.5f}  相对差={abs(h['lam']-d['lam'])/h['lam']*100:+.3f}%"
          f"   [谐波应≈0（电子量）]")
    print(f"  {'':<11} ω̄₂_H/ω̄₂_D={h['w2']/d['w2']:.6f}   [√2={np.sqrt(2):.6f}，"
          f"偏离 {abs(h['w2']/d['w2']-np.sqrt(2))/np.sqrt(2)*100:+.3f}%]")
    print(f"  {'':<11} M2_H={h['M2']:.1f}  M2_D={d['M2']:.1f}   M2_H/M2_D={h['M2']/d['M2']:.6f}  [应≈2]")

print("\n" + "=" * 116)
print("★ 判词")
print("  (1) D3S 的 ω̄₂ 数值本身：见上面逐行相对差。")
print("  (2) ⚠⚠ **但这两个数不能写进 CSV 第 8/9 行**：那两行的 source_author=`Errea et al.`、")
print("      source_identifier=`arXiv:1502.02832`、λ=2.64/1.87 来自 Errea Table I；")
print("      而这里的 ω̄₂ 是 **Mishra .a2f** 积分值。混填会造出一行「λ 来自 Errea、ω̄₂ 来自 Mishra」")
print("      的**跨构建行** ⟹ 其 M2=λ·ω̄₂² 不是任何一次计算的输出。")
print("      `theory/data/external/MANIFEST.md` 明令 Mishra 与 Errea 两套不得混用。")
print("      正确动作：**新增两行 team-internal D3S 行**（镜像第 2/3 行的 H3S），")
print("      λ 与 ω_log 也用本脚本重算值，整行同一构建。")
