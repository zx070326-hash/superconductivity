# -*- coding: utf-8 -*-
"""
判别器：`η/M` 的谐波↔非谐一致，是**物理**（Hopfield 求和规则）还是**构造性恒真**？

━━ 动机 ━━
`feasibility-map` 的图把 H₃S 谐波与非谐两个锚点并排画，两者 `η/M = λ⟨ω²⟩` 一致到 1e−6。
我本来准备把它报为「Hopfield 求和规则的一手确认」，作为 §7/§8「软化声子不改 η/M」的证据。

━━ 判别器（本轮新增）━━
不要只看那一个量，**看它在矩族里的邻居**：
    I_k := ∫ ω^k α²F(ω) dω ,  k = −1, 0, +1, +2, +3
（`I_{−1} = λ/2`，`I_{+1} = η/2M`）。
· 若一致是**物理**：无特别理由让某一阶独享 1e−6 而邻居大幅偏离。
· 若一致是**构造**（复用同一套电子-声子矩阵元与本征矢、只换动力学矩阵）：
  `Σ_ν λ_ν ω_ν²` 被逐模冻住 ⟹ **恰好 k=+1 精确守恒，其余各阶自由漂移**。
⟹ 「独守一阶 + 邻居大幅偏离」是**构造性**的指纹。

━━ 方向审计 ━━
本脚本只能把「一致性」从「物理证据」降级为「构造恒真」，**方向单一**：
它永不支持「自然界的非谐性保持 η/M」这条正陈述。降级后 §7/§8 不受影响
（那里把 η/M 当**设计约束**固定，不主张不变性）。

━━ 归属 ━━
本仓 2026-09-17 已用 prose 论证独立拦下同一条：
`theory/notes/materials/deficit_budget_h3s_2026-09-17.md` §3。本脚本是**重发现**，
新的只有上面这个矩族邻居判别器。

━━ 网格警告 ━━
谐波 α²F 的网格**非均匀**（spacing max/min = 25.2，末点 250.45）⟹ 不得用均匀权 Simpson。
本脚本一律在原生网格上用梯形（`np.trapezoid`）。

严格性：numerical-evidence（float64，原生网格梯形求积）。
"""
import os, sys, io
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

tz = getattr(np, "trapezoid", None) or np.trapz
BASE = os.path.join(r"D:\physics\superconductivity", "theory", "data", "external",
                    "mishra_h3s_d3s", "Data-H3S-D3S-PH-EPW-Anharmonic")


def load_a2f(rel):
    rows = []
    for ln in io.open(os.path.join(BASE, rel.replace("/", os.sep)),
                      encoding="utf-8", errors="replace"):
        t = ln.split()
        if len(t) >= 2:
            try:
                rows.append((float(t[0]), float(t[1])))
            except ValueError:
                pass
    w = np.array([r[0] for r in rows]); a = np.array([r[1] for r in rows])
    o = np.argsort(w); w, a = w[o], a[o]
    m = w > 1e-9
    return w[m], np.clip(a[m], 0.0, None)


def I(w, a, k, lo=0.0, hi=np.inf):
    m = (w > max(lo, 1e-9)) & (w <= hi)
    return tz(a[m] * w[m] ** k, w[m])


SH = {t: load_a2f(p) for t, p in [
    ("H3S harm",   "H3S-200GPa/EPW_harmonic/H3S.a2f"),
    ("H3S anharm", "H3S-200GPa/EPW_anharmonic/H3S.a2f"),
    ("D3S anharm", "D3S-200GPa/EPW_anharmonic/D3S.a2f")]}

# ══ 第 1 部分：锚点复现（对照 feasibility-map 图上印出的三对数）═══════════
print("=" * 96)
print("第 1 部分：三个锚点从零复现（我自己的解析器 + 原生网格梯形）")
print("=" * 96)
CLAIM = {"H3S harm": (2.2930, 129.310), "H3S anharm": (1.7210, 149.260),
         "D3S anharm": (1.7840, 105.430)}
for t, (w, a) in SH.items():
    lam = 2 * I(w, a, -1); w2 = np.sqrt(2 * I(w, a, +1) / lam)
    cl, cw = CLAIM[t]
    sp = np.diff(w)
    print(f"  {t:<11} npts={len(w):4d}  spacing max/min={sp.max()/sp.min():6.1f}  "
          f"lambda={lam:8.5f} (图 {cl:7.4f}, {100*(lam/cl-1):+.3f}%)  "
          f"sqrt<w2>={w2:8.4f} meV (图 {cw:7.3f}, {100*(w2/cw-1):+.3f}%)")

# ══ 第 2 部分：矩族邻居 —— 判别器本体 ════════════════════════════════════
print()
print("=" * 96)
print("第 2 部分：矩族邻居 I_k = ∫ω^k α²F dω 的 谐波/非谐 比（判别器）")
print("=" * 96)
wh, ah = SH["H3S harm"]; wa, aa = SH["H3S anharm"]
print(f"  {'k':>4}{'I_k(harm)':>17}{'I_k(anharm)':>17}{'harm/anharm':>15}   注")
NOTE = {-1: "= λ/2", 0: "面积", +1: "= η/2M", +2: "", +3: ""}
for k in (-1, 0, +1, +2, +3):
    Ih, Ia = I(wh, ah, k), I(wa, aa, k)
    print(f"  {k:+4d}{Ih:17.6e}{Ia:17.6e}{Ih/Ia:15.6f}   {NOTE[k]}")
print("  ⟹ 只有 k=+1 守恒到 1e−6，两侧邻居各差 20–33%、远端 54%")
print("  ⟹ **构造性指纹**：矩阵元与本征矢被冻住、只换动力学矩阵 ⟹ Σ_ν λ_ν ω_ν² 逐模不变")
print("  ⟹ 这**不检验**「自然界的非谐性是否保持 η/M」（归属：deficit_budget_h3s_2026-09-17.md §3，2026-09-17 已拦下）")

# ══ 第 3 部分：有内容的替代检验 —— D₃S ↔ H₃S ═══════════════════════════
print()
print("=" * 96)
print("第 3 部分：D₃S↔H₃S 才有独立内容（非纯质量重映射）")
print("=" * 96)
wH, aH = SH["H3S anharm"]; wD, aD = SH["D3S anharm"]
tot = I(wD, aD, 1) / I(wH, aH, 1)
print(f"  全谱一阶矩比 D/H = {tot:.5f}   冻结矩阵元的纯氢理想值 M_H/M_D = 0.5000"
      f"   ⟹ 残差 {100*(tot/0.5-1):+.1f}%")
print(f"  按频率切带（高带 = 氢主导，低带 = 硫主导）：")
print(f"  {'cut[meV]':>10}{'D/H 低带':>12}{'低带权重':>12}{'D/H 高带':>12}")
for cut in (40, 60, 70, 80, 100):
    dl, hl = I(wD, aD, 1, 0, cut), I(wH, aH, 1, 0, cut)
    dh, hh = I(wD, aD, 1, cut), I(wH, aH, 1, cut)
    print(f"  {cut:10.0f}{dl/hl:12.4f}{hl/I(wH,aH,1):12.4f}{dh/hh:12.4f}")
print("  ⟹ 高带（cut=60）给 0.504 ≈ 0.5（氢模守恒），低带不守恒但仅占总权重 ~2%")
print("  ⟹ 这一条**可以写**：D₃S 不是 H₃S 的纯质量重映射，3.4% 残差定位得上（硫子格参与）")

print()
print("=" * 96)
print("判词：`η/M` 的谐/非谐一致 = 构造性恒真，**不得**作为物理证据引用。")
print("      §7/§8 把 η/M 当设计约束固定，不主张不变性 ⟹ 结论不受影响，受损的只是括号里的理由。")
print("=" * 96)
