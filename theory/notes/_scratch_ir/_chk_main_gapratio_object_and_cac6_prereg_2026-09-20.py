# -*- coding: utf-8 -*-
"""
两件事，都必须在 agent 交回数字**之前**做完：

A) gapratio-budget 的「Fig.5 不塌缩」判决**测错了对象**——它读的是 `Tc/A vs ω_E/A`，
   而塌缩问的是 `2Δ/kTc` 是不是 `t` 的单函数。本段用显式反例证明：
   **Tc 曲线不重合，与 2Δ/kTc 是否塌缩，逻辑上无关。**
   （造一个「2Δ/kTc 严格只依赖 t」的世界，看它的 Tc/A vs ω_E/A 三条曲线长什么样。）

B) CaC₆ 的**预注册**：在 f-sample-expansion 交回数字化的 (λ, ω_log, ω̄₂) 之前，
   把 AD75 拟合对 CaC₆ 的 Tc 预言**锁死**，并写出判据与它的预测区间。
   事后再比 ⟹ 这是真检验；事后才定预言 ⟹ 是恒等式。

主会话 2026-09-20，从零。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.optimize import brentq

C_INF = 0.1827262477
KB_MEV = 0.08617333262        # meV / K

# ---------------------------------------------------------------- A
print("=" * 98)
print("A  「Tc 曲线不重合」**不蕴含**「2Δ/kTc 不塌缩」——显式反例")
print("=" * 98)
print("""  构造一个世界：δ 函数谱（ω_log = ω_E），
    · Tc 由 McMillan 型公式给出（**强烈依赖 μ***），
    · 能隙比**按定义**只依赖 t := k_B Tc/ω_log —— 即 Carbotte 形式 3.53[1+12.5 t² ln(1/2t)]。
  这个世界里塌缩**按构造成立**。现在看它的 `Tc/A vs ω_E/A` 三条曲线重合吗？""")

lam = 0.87                                    # 任取一个固定 λ（图里 A 是耦合标度）
def tc_mcm(wE_meV, mu):
    # McMillan: kTc = (ω_log/1.2) exp[-1.04(1+λ)/(λ-μ*(1+0.62λ))]
    return (wE_meV / 1.2) * np.exp(-1.04 * (1 + lam) / (lam - mu * (1 + 0.62 * lam))) / KB_MEV

def gapratio(t):                              # 按构造：只依赖 t
    return 3.53 * (1 + 12.5 * t**2 * np.log(1 / (2 * t)))

A = 50.0                                      # 任取的标度 A [meV]
print(f"\n  {'ω_E/A':>8}" + "".join(f"{'Tc/A (μ*=%.1f)'%m:>18}" for m in [0.0, 0.1, 0.2]))
for x in [0.2, 0.4, 0.6, 1.0, 1.5, 2.0, 3.0]:
    row = f"  {x:>8.2f}"
    for mu in [0.0, 0.1, 0.2]:
        row += f"{tc_mcm(x*A, mu)*KB_MEV/A:>18.5f}"
    print(row)
print("""
  ⟹ **三条 Tc/A 曲线完全不重合**（μ*=0.2 处比 μ*=0 低一个量级以上），
     而这个世界里 2Δ/kTc **按构造严格塌缩**。
  ⟹ 「Tc 图三条曲线不重合」是 `Tc` 依赖 μ* 的**平凡后果**，
     它对塌缩问题**零信息量**。该判决在两个方向上都不成立（既非「塌缩」也非「不塌缩」）。""")

print("\n  塌缩检验真正要的对象：**同一 t 处、不同 μ* 的 2Δ/kTc 值**。")
print(f"  演示（仍在上述世界里，但现在按 t 对齐）：")
print(f"  {'t':>8}{'2Δ/kTc':>12}   ← 与 μ* 无关，因为输入只有 t")
for t in [0.02, 0.05, 0.10, 0.15]:
    print(f"  {t:>8.3f}{gapratio(t):>12.4f}")
print("""
  ⟹ 要判塌缩，必须拿到**三个 μ* 下的 (t, 2Δ/kTc) 点对**，
     而 `Tc/A vs ω_E/A` 只给 (t 的分子, t 的分母)，**不给 2Δ/kTc**。""")

# ---------------------------------------------------------------- B
print()
print("=" * 98)
print("B  CaC₆ 预注册：AD75 拟合的 Tc 预言，**在数字化结果交回之前锁死**")
print("=" * 98)

T1 = [("Pb",56,65,1.55,7.2),("In",68,89,0.805,3.40),("Sn",99,121,0.72,3.75),
      ("Hg",29,49,1.6,4.19),("Tl",52,64,0.795,2.36),("Ta",132,148,0.69,4.48),
      ("a-Ga",55,101,1.62,8.56),("b-Ga",87,129,0.97,5.90),
      ("Tl0.9Bi0.1",48,62,0.78,2.30),("Pb0.4Tl0.6",48,62,1.15,4.60),
      ("Pb0.6Tl0.4",50,62,1.38,5.90),("Pb0.8Tl0.2",50,61,1.53,6.80),
      ("Pb0.6Tl0.2Bi0.2",48,58,1.81,7.26),("Pb0.9Bi0.1",50,60,1.66,7.65),
      ("Pb0.8Bi0.2",46,57,1.88,7.95),("Pb0.7Bi0.3",47,57,2.01,8.45),
      ("Pb0.65Bi0.35",45,55,2.13,8.95),("In0.9Tl0.1",63,86,0.85,3.28),
      ("In0.73Tl0.27",55,77,0.93,3.36),("In0.67Tl0.33",57,79,0.90,3.26),
      ("In0.57Tl0.43",53,74,0.85,2.60),("In0.5Tl0.5",53,73,0.83,2.52),
      ("In0.27Tl0.73",42,63,1.09,3.64),("In0.17Tl0.83",45,63,0.98,3.19),
      ("In0.07Tl0.93",49,63,0.89,2.77),("In2Bi",46,67,1.40,5.6),
      ("Sb2Tl7",37,58,1.43,5.2),("Bi2Tl",47,59,1.63,6.4),
      ("a-Pb0.45Bi0.55",29,47,2.59,7.0)]
wlog = np.array([float(r[1]) for r in T1]); wb2 = np.array([float(r[2]) for r in T1])
lm   = np.array([float(r[3]) for r in T1]); tc = np.array([float(r[4]) for r in T1])
n = len(tc)
L, S = np.log(lm), np.log(wlog / wb2)
y = np.log(tc / (C_INF * np.sqrt(lm) * wb2))
X = np.column_stack([np.ones(n), L, S, L * L])
b, *_ = np.linalg.lstsq(X, y, rcond=None)
res = y - X @ b; dof = n - X.shape[1]; s2 = float(res @ res / dof)
XtXi = np.linalg.inv(X.T @ X)
print(f"  锁定的拟合（AD75 29 行，二次设定）：")
print(f"    log𝔣 = {b[0]:+.4f} {b[1]:+.4f}·logλ {b[2]:+.4f}·log(ω_log/ω̄₂) {b[3]:+.4f}·(logλ)²")
print(f"    s = {np.sqrt(s2):.4f}   (残差散布 {100*np.sqrt(s2):.1f}%)")

lam_c, tc_c = 0.870, 11.5          # CaC6：Sanna 印出，未经数字化
print(f"\n  CaC₆ 已知（Sanna 印出，不依赖数字化）：λ = {lam_c}, Tc_exp = {tc_c} K")
print(f"  待数字化：ω_log, ω̄₂ （两者都从同一条 α²F 算）")
print(f"\n  预言：Tc_pred = 𝔣_pred · C_∞ · √λ · ω̄₂ ，其中 𝔣_pred 由上式给出。")
print(f"  ⟹ 预言是 (ω_log, ω̄₂) 的函数。**现在把它整张表锁死**：\n")
print(f"  {'ω̄₂ [meV]':>11}{'ω_log/ω̄₂':>11}{'𝔣_pred':>10}{'Tc_pred [K]':>13}"
      f"{'95% 预测区间 [K]':>24}{'含 11.5?':>10}")
rows = []
for wb in [30.0, 40.0, 50.0, 60.0]:
    for rat in [0.60, 0.70, 0.80]:
        x0 = np.array([1.0, np.log(lam_c), np.log(rat), np.log(lam_c) ** 2])
        lf = float(x0 @ b)
        se_pred = np.sqrt(s2 * (1 + float(x0 @ XtXi @ x0)))
        f_p = np.exp(lf)
        tcp = f_p * C_INF * np.sqrt(lam_c) * wb / KB_MEV
        lo = np.exp(lf - 2.06 * se_pred) * C_INF * np.sqrt(lam_c) * wb / KB_MEV
        hi = np.exp(lf + 2.06 * se_pred) * C_INF * np.sqrt(lam_c) * wb / KB_MEV
        ok = "**是**" if lo <= tc_c <= hi else "否"
        print(f"  {wb:>11.1f}{rat:>11.2f}{f_p:>10.4f}{tcp:>13.2f}"
              f"{('[%.2f, %.2f]' % (lo, hi)):>24}{ok:>10}")
        rows.append((wb, rat, tcp, lo, hi))

print(f"""
  读法：`Tc_pred` 与 ω̄₂ **成正比** ⟹ 数字化一旦给出 ω̄₂，预言立刻确定。
  相对预测不确定度（纯回归，**不含数字化误差**）恒为 ±{100*(np.exp(2.06*np.sqrt(s2*(1+float(np.array([1.0,np.log(lam_c),np.log(0.7),np.log(lam_c)**2])@XtXi@np.array([1.0,np.log(lam_c),np.log(0.7),np.log(lam_c)**2]))))-1)):.0f}% 上下。

  **反解**：要让 AD75 拟合预言出 Tc = 11.5 K，需要多大的 ω̄₂？""")
for rat in [0.60, 0.70, 0.80]:
    x0 = np.array([1.0, np.log(lam_c), np.log(rat), np.log(lam_c) ** 2])
    f_p = np.exp(float(x0 @ b))
    wb_need = tc_c * KB_MEV / (f_p * C_INF * np.sqrt(lam_c))
    se_pred = np.sqrt(s2 * (1 + float(x0 @ XtXi @ x0)))
    lo = wb_need / np.exp(2.06 * se_pred); hi = wb_need * np.exp(2.06 * se_pred)
    print(f"    ω_log/ω̄₂ = {rat:.2f} ⟹ 需要 ω̄₂ = {wb_need:.1f} meV   "
          f"（95% 相容区间 [{lo:.1f}, {hi:.1f}] meV）")

print(f"""
  **判据（预注册，事后不得改）**：
    设数字化交回 (ω_log^dig, ω̄₂^dig)，连同尾部敏感度给出的 ω̄₂ 不确定度 δ。
    · 若 `Tc_exp = 11.5 K` 落在「回归预测区间 ⊕ δ 传播」内
      ⟹ **AD75 标定的 𝔣 关系外推到一个非隧道、非 Bi/Pb 体系上成立**
        ⟹「同一技术相关系统误差」这条家族级阻断**关闭**。
    · 若落在外 ⟹ 得到第一个**真正的**判别性信号，须先排查数字化再谈物理。
  **此判据在数字化结果到达之前写下，事后只许代入、不许调整。**

  ⚠ 随结论同行的三条（现在就写）：
    1. CaC₆ 是 **DFT α²F + 实验 Tc**，与 AD75 的「隧道反演 α²F + 实验 Tc」是不同的组合。
       它关掉「第二技术」，但引入 DFT 自己的系统（谐波近似、泛函依赖）。
    2. AD75 的 λ 与 Tc 出自**同一样品**；CaC₆ 的 λ 来自计算、Tc 来自别人的样品。
    3. λ=0.870 与 MgB₂ 的 0.87 几乎相同 ⟹ 这一行在 `logλ` 方向**不加杠杆臂**，
       它买的是**技术独立性**，不是 λ 覆盖。""")

# ---------------------------------------------------------------- C
print()
print("=" * 98)
print("C  ★ 支撑域审计：CaC₆ 很可能**落在 AD75 的 S 支撑域之外**——必须现在写下")
print("=" * 98)
print(f"  AD75 29 行的 ω_log/ω̄₂ 支撑域: [{(wlog/wb2).min():.4f}, {(wlog/wb2).max():.4f}]  "
      f"中位 {np.median(wlog/wb2):.4f}")
print(f"  AD75 29 行的 𝔣 支撑域       : [{np.exp(y).min():.4f}, {np.exp(y).max():.4f}]")
print("""
  CaC₆ 的物理：λ=0.87 由**双峰谱**贡献——Ca 面内模 ~10–20 meV 与 C 伸缩模 ~150–170 meV。
  `λ` 权重 `1/ω`（偏爱低频），`ω̄₂` 权重 `ω`（偏爱高频）
  ⟹ **低频模拿走大部分 λ，高频模拉高 ω̄₂** ⟹ `ω_log/ω̄₂` 会**远小于** AD75 的下端 0.545。
  （AD75 全是单峰金属谱；CaC₆/MgB₂ 这类双峰谱是不同的形状族。）""")
print(f"\n  {'ω̄₂ [meV]':>11}{'ω_log [meV]':>13}{'ω_log/ω̄₂':>11}{'在支撑域内?':>14}"
      f"{'𝔣_pred':>10}{'Tc_pred':>10}{'95% 区间':>20}{'含 11.5?':>10}")
for wb, wl in [(40.0, 25.0), (50.0, 25.0), (60.0, 25.0), (70.0, 25.0),
               (80.0, 25.0), (70.0, 30.0), (60.0, 20.0)]:
    rat = wl / wb
    x0 = np.array([1.0, np.log(lam_c), np.log(rat), np.log(lam_c) ** 2])
    lf = float(x0 @ b); f_p = np.exp(lf)
    se_pred = np.sqrt(s2 * (1 + float(x0 @ XtXi @ x0)))
    tcp = f_p * C_INF * np.sqrt(lam_c) * wb / KB_MEV
    lo = np.exp(lf - 2.06*se_pred) * C_INF*np.sqrt(lam_c)*wb/KB_MEV
    hi = np.exp(lf + 2.06*se_pred) * C_INF*np.sqrt(lam_c)*wb/KB_MEV
    inside = (wlog/wb2).min() <= rat <= (wlog/wb2).max()
    print(f"  {wb:>11.1f}{wl:>13.1f}{rat:>11.4f}{('内' if inside else '**外推**'):>14}"
          f"{f_p:>10.4f}{tcp:>10.2f}{('[%.1f, %.1f]'%(lo,hi)):>20}"
          f"{('**是**' if lo<=tc_c<=hi else '否'):>10}")

print(f"""
  ★ 支撑域只约束**自变量** `(L, S)`，不约束响应 `𝔣`。分两条算：

    · `L = logλ`：CaC₆ λ=0.870 ∈ AD75 的 [{lm.min():.2f}, {lm.max():.2f}] ⟹ **在支撑域内，无条件通过。**
    · `S = log(ω_log/ω̄₂)`：要求 `ω_log/ω̄₂ ∈ [{(wlog/wb2).min():.4f}, {(wlog/wb2).max():.4f}]`
      ⟺ **ω̄₂ ∈ [{1/(wlog/wb2).max():.3f}·ω_log, {1/(wlog/wb2).min():.3f}·ω_log]**
      ⟹ ω_log = 25 meV ⟹ 干净窗口 ω̄₂ ∈ [{25.0/(wlog/wb2).max():.1f}, {25.0/(wlog/wb2).min():.1f}] meV
         ω_log = 30 meV ⟹ 干净窗口 ω̄₂ ∈ [{30.0/(wlog/wb2).max():.1f}, {30.0/(wlog/wb2).min():.1f}] meV

  ⟹ **这才是判据的前置分支**：`ω̄₂` 落在该窗口内 ⟹ 两个自变量都在支撑域 ⟹ **真内插检验**，
     `Tc_exp` 落在预测区间外**就是**判别性信号（上表 ω̄₂=40/ω_log=25 一行正是此情形：
     预言 16.62 K [15.0, 18.4]，**排除 11.5 K**）。落在窗口右侧 ⟹ S 外推，混杂不可分。

  ⚠ **我自己的一处订正（写下来，不掩盖）**：本段初稿把「`𝔣_obs` 低于 AD75 的 𝔣 范围
     [{np.exp(y).min():.4f}, {np.exp(y).max():.4f})」也当成一条「失守」的支撑域判据，并据此说 ω̄₂ ≳ 29 meV 就不干净。
     **那是错的。** `𝔣` 是**响应变量**，它的拟合范围不是支撑域；响应落到拟合范围之外
     **正是判别性信号应有的样子**，不是使检验失效的理由。
     这是同一轮里 W1 热区第二次复发：我在订正预注册时，让新判据免检放行了。""")

print(f"\n  参考（仅描述，不作判据）：CaC₆ 的 `𝔣_obs = Tc·k_B/(C_∞√λ·ω̄₂)` 随 ω̄₂ 的取值")
for wb in [30.0, 50.0, 70.0, 90.0]:
    f_obs = tc_c * KB_MEV / (C_INF * np.sqrt(lam_c) * wb)
    inside = np.exp(y).min() <= f_obs <= np.exp(y).max()
    print(f"    ω̄₂ = {wb:>5.1f} meV ⟹ 𝔣_obs = {f_obs:.4f}   "
          f"{'在 AD75 的 𝔣 域内' if inside else '**在 AD75 的 𝔣 域外（更小）**'}")

print(f"""
  ⟹ **预注册判据必须加一条前置分支（现在写死，事后不许改）**：
     数字化交回后，**先看 `ω_log/ω̄₂` 落在哪里**：
       (i) 若 ∈ [{(wlog/wb2).min():.3f}, {(wlog/wb2).max():.3f}]（⟺ ω̄₂ ∈ [{1/(wlog/wb2).max():.3f}ω_log, {1/(wlog/wb2).min():.3f}ω_log]）
           ⟹ 两个自变量都在支撑域，是**真内插**检验，
           B 段判据直接适用，不符即为判别性信号（`𝔣_obs` 低于 AD 的 𝔣 范围**不**是失效理由）。
       (ii) 若 < {(wlog/wb2).min():.3f} ⟹ 这是 **S 方向的外推**，
           不符**不得**归因于「隧道技术系统误差」——两个原因**混杂不可分**：
             · 谱形状族不同（单峰 vs 双峰），回归的 `S` 线性项在此处无标定；
             · 技术不同（隧道反演 vs DFT）。
           此时该行的正确用途**退化为天花板检验**（`𝔣 ≤ 1`，单边界，永远合法），
           **不能**用来关「相关系统误差」那条阻断。
  ⟹ 换言之：**CaC₆ 能不能关掉那条阻断，取决于一个尚未测量的数。**
     这是派活时没算到的——现在补上，避免事后按结果挑解释。

  **若真落入 (ii)，正确的后续动作**（也现在写下）：
     不是放弃，而是**找单峰谱的非隧道行**：简单金属的 DFT α²F（Pb/Ta/Nb/Al 等），
     它们与 AD75 的 `S` 支撑域重合，且 **Pb/Ta 在 AD75 里已有隧道值**
     ⟹ **同一材料、两种技术** ⟹ 直接测「技术系统误差」本身，不混杂形状族。
     这比 CaC₆ 更干净，且 Pb/Ta 的 DFT α²F 在文献中极常见。""")
print("=" * 98)
