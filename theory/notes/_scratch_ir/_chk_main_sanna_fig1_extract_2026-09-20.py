# -*- coding: utf-8 -*-
"""
Sanna PRB 85,184514 Fig.1 九条 α²F_JJ'(ω) 的矢量提取 + 九个 λ_JJ' 闭合验收门。

诊断动机：f-sample-expansion 报「九条闭合系统性偏低 13%–54%」，并提出「全局 ×2」。
**偏差不是常数 ⟹ 全局缩放必错。** 本脚本从零重做，检验一个具体假说：
  轴标定只有两个刻度 (0, 0.4)，若用 `np.interp` 映射，**它会把范围外的值钳位**，
  而曲线峰值确实超过 0.4（y_pdf 最高 71.3 ⟹ α²F ≈ 0.6）⟹ 峰被削平。
  正确做法是**仿射映射**，不是插值。

主会话 2026-09-20，从零。严格性：numerical-evidence。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, fitz

import os
# 相对脚本自身定位 PDF —— 被笔记引用的脚本必须能从任意 cwd 跑起来
_PDF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_sanna_prb85_184514.pdf")
doc = fitz.open(_PDF); pg = doc[3]

# ---------- 1. 轴标定：用文本 bbox 的**中心**，不是左上角 ----------
def centers(pred):
    out = []
    for w in pg.get_text("words"):
        if pred(w): out.append((w[4], (w[0]+w[2])/2, (w[1]+w[3])/2))
    return out

xt = [(float(t), cx) for t, cx, cy in centers(
        lambda w: abs(w[1]-224.2) < 1.0 and w[4] in ("0","50","100","150","200"))]
xt.sort(key=lambda r: r[0])
A = np.polyfit([c for _, c in xt], [v for v, _ in xt], 1)
print("X 轴标定（meV = a·x_pdf + b）：", [f"{v:g}meV@{c:.1f}" for v, c in xt])
print(f"   a = {A[0]:.6f}  b = {A[1]:.4f}   ⟹ 残差 "
      f"{np.max(np.abs(np.polyval(A,[c for _,c in xt])-[v for v,_ in xt])):.4f} meV")

yt = centers(lambda w: w[4] in ("0","0.4") and 55 < w[0] < 70 and 80 < w[1] < 225)
ylab = sorted([(float(t), cy) for t, cx, cy in yt], key=lambda r: r[1])
panels = []
for i in range(3):
    v4, y4 = ylab[2*i]; v0, y0 = ylab[2*i+1]
    assert (v4, v0) == (0.4, 0.0), (v4, v0)
    panels.append((y0, y4))
    print(f"Y 轴 panel {i+1}: α²F=0 @ y={y0:.2f} , 0.4 @ y={y4:.2f}  (跨度 {y0-y4:.2f} pt)")

# ---------- 2. 取路径点（三次 Bezier 真采样，不拿控制点当数据） ----------
COL = {(0.008,0.545,0.267): 1, (0.929,0.125,0.141): 3, (0.208,0.325,0.643): 2}
def bez(p0,p1,p2,p3,n=12):
    t = np.linspace(0,1,n)[:,None]
    return ((1-t)**3*p0 + 3*(1-t)**2*t*p1 + 3*(1-t)*t**2*p2 + t**3*p3)
buckets = {}
for dr in pg.get_drawings():
    col = dr.get("color")
    if col is None: continue
    J2 = COL.get(tuple(round(v,3) for v in col))
    if J2 is None: continue
    if len(dr["items"]) == 1:          # 图例样线：单段短横线，落在 panel 1 内
        continue                        # 不排除 ⟹ 低频区被污染，1/ω 权重放大（panel 1 +12~30%）
    pts = []
    for it in dr["items"]:
        if it[0] == "l":
            pts += [(it[1].x,it[1].y),(it[2].x,it[2].y)]
        elif it[0] == "c":
            P = bez(*[np.array([it[k].x,it[k].y]) for k in (1,2,3,4)])
            pts += [tuple(q) for q in P]
    for x,y in pts:
        for J,(y0,y4) in enumerate(panels,1):
            if y4-14 <= y <= y0+1.0:
                buckets.setdefault((J,J2),[]).append((x,y)); break

# ---------- 3. 两种 y 映射：仿射 vs np.interp（复现对方的 bug） ----------
PRINT = np.array([[0.301,0.136,0.257],[0.546,0.239,0.479],[0.427,0.198,0.367]])
def lam_of(J,J2,mode):
    P = np.array(sorted(buckets[(J,J2)]))
    w = np.polyval(A,P[:,0]); y0,y4 = panels[J-1]
    if mode=="affine":
        a2f = (y0-P[:,1])*0.4/(y0-y4)
    else:                                   # np.interp 钳位版
        a2f = np.interp(P[:,1], [y4,y0], [0.4,0.0])
    o = np.argsort(w); w,a2f = w[o],a2f[o]
    m = w > 1e-9
    return 2*np.trapezoid(a2f[m]/w[m], w[m]), a2f.max()

print(f"\n{'(J,J\u2032)':>8}{'印出 λ':>9}{'仿射 λ':>9}{'误差':>9}{'  |  '}"
      f"{'interp λ':>9}{'误差':>9}{'仿射 max':>10}{'interp max':>11}")
ea, ei = [], []
for J in (1,2,3):
    for J2 in (1,2,3):
        la,ma = lam_of(J,J2,"affine"); li,mi = lam_of(J,J2,"interp")
        ref = PRINT[J-1,J2-1]
        ea.append(100*(la/ref-1)); ei.append(100*(li/ref-1))
        print(f"{'(%d,%d)'%(J,J2):>8}{ref:>9.3f}{la:>9.3f}{ea[-1]:>8.1f}%{'  |  '}"
              f"{li:>9.3f}{ei[-1]:>8.1f}%{ma:>10.3f}{mi:>11.3f}")
print(f"\n仿射  : 平均 {np.mean(ea):+.1f}%  最劣 {max(np.abs(ea)):.1f}%  RMS {np.sqrt(np.mean(np.square(ea))):.1f}%")
print(f"interp: 平均 {np.mean(ei):+.1f}%  最劣 {max(np.abs(ei)):.1f}%  RMS {np.sqrt(np.mean(np.square(ei))):.1f}%")

# ---------- 4. DOS 加权重建 α²F_iso，算 λ / ω_log / ω̄₂ ----------
N3 = np.array([0.412, 0.104, 0.249])          # Table I 三带 DOS（权重挂**首指标**=panel）
grid = np.linspace(0.05, 200.0, 20000)
def curve(J, J2):
    P = np.array(sorted(buckets[(J,J2)])); w = np.polyval(A, P[:,0])
    y0,y4 = panels[J-1]; a = (y0-P[:,1])*0.4/(y0-y4)
    o = np.argsort(w); w,a = w[o],a[o]
    u,idx = np.unique(w, return_index=True)
    return np.interp(grid, u, a[idx], left=0.0, right=0.0)

iso = sum(N3[J-1]*curve(J,J2) for J in (1,2,3) for J2 in (1,2,3)) / N3.sum()
def moments(f, g=grid):
    lam = 2*np.trapezoid(f/g, g)
    wlog = np.exp(2/lam*np.trapezoid(f*np.log(g)/g, g))
    w2 = np.sqrt(2/lam*np.trapezoid(f*g, g))
    return lam, wlog, w2

lam, wlog, wb2 = moments(iso)
print("\n" + "="*96)
print("4. α²F_iso（按 panel 的 DOS 加权）及其矩")
print("="*96)
print(f"  λ      = {lam:.4f}   （Sanna 印出 0.870；本会话 Table I 重算 0.8685）"
      f"  差 {100*(lam/0.8685-1):+.2f}%")
print(f"  ω_log  = {wlog:.2f} meV        ω̄₂ = {wb2:.2f} meV        ω_log/ω̄₂ = {wlog/wb2:.4f}")

# 尾部敏感度：最高频段面积 ±20%
print("\n  尾部敏感度（ω > 140 meV 的面积 ±20%，δ_ω̄₂）：")
for s in (0.8, 0.9, 1.0, 1.1, 1.2):
    f2 = iso.copy(); f2[grid > 140] *= s
    l2, wl2, w22 = moments(f2)
    print(f"    尾 ×{s:.1f}  ⟹  λ={l2:.4f}  ω_log={wl2:.2f}  ω̄₂={w22:.2f} meV"
          f"  ω_log/ω̄₂={wl2/w22:.4f}")

# 预注册判据
S_LO, S_HI = 0.5446, 0.8919
print("\n" + "="*96)
print("5. 预注册前置分支（§10.3，事后不许改）")
print("="*96)
r = wlog/wb2
print(f"  AD75 的 S 支撑域: ω_log/ω̄₂ ∈ [{S_LO}, {S_HI}]")
print(f"  CaC₆ 实测       : {r:.4f}  ⟹  " +
      ("**(i) 支撑域内 ⟹ 真内插检验**" if S_LO <= r <= S_HI else
       f"**(ii) 外推（{'低于下端' if r < S_LO else '高于上端'}）⟹ 退化为天花板检验**"))
C_INF, KB = 0.1827262477, 0.08617333262
f_obs = 11.5*KB/(C_INF*np.sqrt(lam)*wb2)
print(f"\n  𝔣_obs = Tc·k_B/(C_∞√λ·ω̄₂) = {f_obs:.4f}   "
      f"（天花板 𝔣 ≤ 1 {'**满足**' if f_obs <= 1 else '**违反 —— 报警！**'}）")

# ------------------------------------------------------------------ 6
print()
print("=" * 96)
print("6. ★ 判决对标定误差的免疫性：ω_log/ω̄₂ 在 α²F 的整体缩放下**精确不变**")
print("=" * 96)
print("""  为什么重要：第 3 节量出本流水线有 **+2.4% 的整体标定偏差**（九条闭合的均值）。
  若 ω_log/ω̄₂ 随它漂移，(ii) 支判决就可能是标定伪影。但它不漂——

  [analytic] λ = 2∫f/ω dω，而两个矩都带 2/λ 归一：
      f → c·f  ⟹  λ → cλ  ⟹  (2/(cλ))∫(cf)(·)dω = (2/λ)∫f(·)dω     （c 精确约掉）
  ⟹ ω_log 与 ω̄₂ **各自**不变，比值自然不变。只有**形状**误差进得来（已由 §4 尾扫覆盖）。
  下面是数值确认（不是证明，只作实现的验收）：""")
print(f"\n  {'缩放 c':>10}{'λ':>10}{'ω_log':>11}{'ω̄₂':>11}{'ω_log/ω̄₂':>14}{'𝔣_obs':>10}")
for c in [0.5, 0.976, 1.0, 1.024, 2.0]:
    l_, wl_, w2_ = moments(iso * c)
    print(f"  {c:>10.3f}{l_:>10.4f}{wl_:>11.4f}{w2_:>11.4f}{wl_/w2_:>14.9f}"
          f"{11.5*KB/(C_INF*np.sqrt(l_)*w2_):>10.4f}")
print("""
  ⟹ 比值逐位相同 ⟹ **(ii) 支判决不受 +2.4% 标定偏差影响。**
  ⚠ 反过来要记住：`λ` 与 `𝔣_obs` **不**免疫（λ ∝ c，𝔣_obs ∝ 1/√c）
    ⟹ 若天花板检验曾逼近 1，标定偏差就承重。本例 𝔣=0.0845 离 1 很远 ⟹ 不承重。""")
print("=" * 96)
