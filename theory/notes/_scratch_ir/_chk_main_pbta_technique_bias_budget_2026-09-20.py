# -*- coding: utf-8 -*-
"""
派活前预算（CLAUDE.md「预算优先」）：`f-sample-expansion` 想数字化 Savrasov Fig.2
把 Pb/Ta 的 DFT `ω̄₂` 抠出来。**先问：它买到什么？能不能威胁到头条？**

头条现在是**曲率**陈述：  log 𝔣 = a + b·logλ + c·logS + q·(logλ)² ，q = −0.525。
家族级阻断是「29 行同出隧道反演一种技术 ⟹ 系统误差可能是人造的」。

本脚本问三件事：
  Q1  [analytic] 若技术偏置是**常数倍率**（λ_D=ρλ_T、ω̄₂_D=σω̄₂_T，ρ、σ 与 λ 无关），
      它动得了 q 吗？—— 符号推演 + 数值确认。
  Q2  用 agent 转述的 Savrasov 数（**未经我核验**）实测 Pb/Ta 的技术偏置有多大，
      它把**两点斜率**推动多少？
  Q3  两个材料能不能约束 q？—— 自由度计数。

方向注记：本脚本不产出单边界判决；Q2 是双侧比较，Q3 是自由度计数（非不等式）。
严格性：Q1 analytic（+数值确认）；Q2/Q3 numerical-evidence，且**输入未经核验**。

主会话 2026-09-20，从零。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

C_INF = 0.1827262477

# ---- AD75 Table I（隧道反演），单位 K，页图转录，𝔣 已独立复现 ----
TUN = {"Pb": dict(lam=1.55, wlog=56.0,  wb2=65.0,  tc=7.20),
       "Ta": dict(lam=0.69, wlog=132.0, wb2=148.0, tc=4.48)}

# ---- Savrasov（DFT）：**agent 转述，我未核验页图** ----
DFT = {"Pb": dict(lam=1.68, wlog=65.0),      # ω̄₂ 未列表 ⟹ 正是要数字化的东西
       "Ta": dict(lam=0.86, wlog=160.0)}

def frak_f(lam, wb2, tc):
    return tc/(C_INF*np.sqrt(lam)*wb2)       # Tc[K] 与 ω̄₂[K] 同单位，无 k_B

print("="*88)
print("Q1 [analytic] 常数倍率技术偏置动得了曲率 q 吗？")
print("="*88)
print("""  设 λ_D = ρ·λ_T ，ω̄₂_D = σ·ω̄₂_T ，ρ、σ 为常数（与 λ 无关）；Tc 两技术相同（同一实验量）。
    log 𝔣_D = log 𝔣_T − ½logρ − logσ        ← 只是**常数平移**
    log λ_D = log λ_T + logρ                 ← 只是**常数平移**
  代入 log𝔣 = a + b·L + q·L²（L:=logλ）：
    (L_T + logρ)² = L_T² + 2logρ·L_T + (logρ)²
  ⟹ L² 的系数**原样不动**；变的只有 b（吸收 2q·logρ）与 a。
  ★ 结论：**q 对常数倍率的技术偏置严格不变** [analytic]。要动 q，偏置必须**随 λ 变且自身有曲率**。""")

# 数值确认 Q1
rng = np.random.default_rng(7)
L = np.log(rng.uniform(0.7, 2.6, 400)); S = np.log(rng.uniform(0.55, 0.89, 400))
y = 0.3 + 1.0*L + 0.7*S - 0.525*L**2 + rng.normal(0, 0.02, 400)
def fit(Lv, yv, Sv):
    A = np.column_stack([np.ones_like(Lv), Lv, Sv, Lv**2])
    return np.linalg.lstsq(A, yv, rcond=None)[0]
b0 = fit(L, y, S)
for rho, sig in [(1.0, 1.0), (1.25, 1.15), (0.8, 1.4)]:
    b1 = fit(L + np.log(rho), y - 0.5*np.log(rho) - np.log(sig), S)
    print(f"  ρ={rho:<5} σ={sig:<5} ⟹ q = {b1[3]:+.6f}  (基准 {b0[3]:+.6f}，"
          f"差 {abs(b1[3]-b0[3]):.2e})   b = {b1[1]:+.4f} (基准 {b0[1]:+.4f})")
print("  ⟹ q 逐位不动、b 动 ⟹ 符号推演确认。")

print()
print("="*88)
print("Q2 实测技术偏置（⚠ Savrasov 数为 agent 转述，未经我核验页图）")
print("="*88)
print(f"  {'':<5}{'λ_tun':>8}{'λ_DFT':>8}{'ρ=λD/λT':>10}"
      f"{'ωlog_tun':>10}{'ωlog_DFT':>10}{'σ_log':>8}")
for m in ("Pb", "Ta"):
    rho = DFT[m]["lam"]/TUN[m]["lam"]; sg = DFT[m]["wlog"]/TUN[m]["wlog"]
    print(f"  {m:<5}{TUN[m]['lam']:>8.2f}{DFT[m]['lam']:>8.2f}{rho:>10.4f}"
          f"{TUN[m]['wlog']:>10.0f}{DFT[m]['wlog']:>10.0f}{sg:>8.4f}")

rho_Pb = DFT["Pb"]["lam"]/TUN["Pb"]["lam"]; rho_Ta = DFT["Ta"]["lam"]/TUN["Ta"]["lam"]
dl = np.log(TUN["Pb"]["lam"]/TUN["Ta"]["lam"])
slope_rho = np.log(rho_Pb/rho_Ta)/dl
sg_Pb = DFT["Pb"]["wlog"]/TUN["Pb"]["wlog"]; sg_Ta = DFT["Ta"]["wlog"]/TUN["Ta"]["wlog"]
slope_sig = np.log(sg_Pb/sg_Ta)/dl
print(f"\n  ρ 不是常数：dlogρ/dlogλ = {slope_rho:+.4f}   （低 λ 端偏置更大）")
print(f"  σ(用 ω_log 代理) ：dlogσ/dlogλ = {slope_sig:+.4f}")

# 两点斜率：隧道 vs DFT（DFT 的 ω̄₂ 用 ω_log 的偏置做代理——正是待测量的假设）
fT = {m: frak_f(TUN[m]["lam"], TUN[m]["wb2"], TUN[m]["tc"]) for m in TUN}
fD = {m: frak_f(DFT[m]["lam"], TUN[m]["wb2"]*DFT[m]["wlog"]/TUN[m]["wlog"], TUN[m]["tc"])
      for m in TUN}
sT = np.log(fT["Pb"]/fT["Ta"])/np.log(TUN["Pb"]["lam"]/TUN["Ta"]["lam"])
sD = np.log(fD["Pb"]/fD["Ta"])/np.log(DFT["Pb"]["lam"]/DFT["Ta"]["lam"])
print(f"\n  𝔣(隧道)  Pb={fT['Pb']:.4f}  Ta={fT['Ta']:.4f}   两点斜率 {sT:+.4f}")
print(f"  𝔣(DFT*)  Pb={fD['Pb']:.4f}  Ta={fD['Ta']:.4f}   两点斜率 {sD:+.4f}")
print(f"  * DFT 的 ω̄₂ 用 ω_log 偏置代理 —— **这正是数字化要替换掉的假设**")
print(f"  ⟹ 两点斜率变 {100*(sD/sT-1):+.1f}%  （对照全样本 b = 0.9972 ± 0.0405）")

print()
print("="*88)
print("Q3 两个材料够不够约束 q？—— 自由度")
print("="*88)
print("""  要在「技术偏置」里看出**曲率**，需要拟合 logρ(λ)、logσ(λ) 的二次项 ⟹ 每条 ≥3 个材料。
  Pb + Ta = **2 点** ⟹ 只能定一条直线，**曲率自由度为 0**。
  ⟹ Pb/Ta **无法**约束 q，无论 ω̄₂ 数字化做得多好。""")

print("""
  ★ 预算判决：
    · 对**曲率头条 q**：预算 = ∞（常数偏置严格不动 q，而 2 点测不出曲率偏置）
      ⟹ 这条活**买不到**头条的安全性，不要以「关掉家族阻断」为名派它。
    · 对**次级陈述 b ≈ 1**：预算**有限且已被咬穿**——两点斜率就动了几十个百分点
      ⟹ 这条活**有信息量**，消费点 = 「b 的技术系统误差」这一行，不是 q。
    · ω̄₂ 数字化的**独有**价值：把「σ 可用 ω_log 代理」这个**假设**换成**测量**。
      上表整列 DFT 𝔣 都建立在该代理上 ⟹ 不数字化就没有一个干净的数。""")
print("="*88)
