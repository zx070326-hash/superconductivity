# -*- coding: utf-8 -*-
"""
核验 f-sample-expansion 报的 MgB₂ 𝔣≈3.85 / CaC₆ 𝔣≈2.25。

𝔣 > 1 意味着天花板 `Tc ≤ C_∞√(λ⟨ω²⟩)` 被违反，而那是**定理** ⟹ 先假定是 bug。

顺带做一件有用的事：ω̄₂ 缺失并**不**妨碍把这两行喂进**天花板检验**——
幂平均不等式给 `ω_log ≤ ω̄₂`（AD75 29/29 实测成立），于是

        𝔣 = Tc/(C_∞√λ·ω̄₂)  ≤  Tc/(C_∞√λ·ω_log) =: 𝔣_up

**方向审计**：`ω̄₂ ≥ ω_log ⟹ 1/ω̄₂ ≤ 1/ω_log ⟹ 𝔣 ≤ 𝔣_up`。
`𝔣_up ≤ 1` ⟹ 天花板 PASS。这是**上**界，按单边界规则**只支持 PASS，不支持 FAIL**：
若 `𝔣_up > 1` 什么也说明不了（必须拿到真 ω̄₂）。

主会话 2026-09-20，从零。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import math

C_INF = 0.1827262477      # 无量纲
KB    = 0.08617333262     # meV / K

print("=" * 92)
print("A. 复现 agent 的数：它到底算了什么？")
print("=" * 92)
# agent: MgB2 Tc=39.4, lam=0.87, "w_log~60 meV" -> f ~ 3.85
#        CaC6 Tc=11.5, lam=0.870,"w_log~30 meV" -> f ~ 2.25
for tag, tc, lam, w_meV, agent_f in [("MgB2", 39.4, 0.87, 60.0, 3.85),
                                     ("CaC6", 11.5, 0.870, 30.0, 2.25)]:
    f_bug = tc / (C_INF * math.sqrt(lam) * w_meV)          # meV 当成 K 用
    w_K   = w_meV / KB
    f_ok  = tc / (C_INF * math.sqrt(lam) * w_K)            # 正确：换算成 K
    print(f"  {tag}: 把 ω 的 meV 数值直接当 K 用 ⟹ 𝔣 = {f_bug:.3f}"
          f"   (agent 报 {agent_f})   差 {abs(f_bug/agent_f-1)*100:.1f}%")
    print(f"        正确换算 ω_log = {w_meV:.1f} meV = {w_K:7.1f} K ⟹ 𝔣_up = {f_ok:.4f}")
    print(f"        两者之比 = 1/k_B = {1/KB:.3f}  ⟹ **单位 bug 定位确凿**")

print()
print("=" * 92)
print("B. 正确的天花板检验（只用 ω_log，ω̄₂ 缺失不妨碍）")
print("=" * 92)
print("   𝔣_up := Tc/(C_∞√λ·ω_log) ≥ 𝔣   （因 ω_log ≤ ω̄₂）")
print()
print(f"{'体系':<8}{'Tc[K]':>8}{'λ':>8}{'ω_log[meV]':>12}{'ω_log[K]':>11}"
      f"{'天花板↓[K]':>12}{'𝔣_up':>9}   判决")
# ω_log 取文献区间，逐点扫（agent 给的是"文献典型值"，未从原文取 ⟹ 我按区间报）
CASES = [("MgB2", 39.4, 0.87,  [50.0, 60.0, 62.0, 70.0]),
         ("CaC6", 11.5, 0.870, [20.0, 25.0, 30.0, 35.0])]
for tag, tc, lam, wl_list in CASES:
    for wl in wl_list:
        w_K  = wl / KB
        ceil_lo = C_INF * math.sqrt(lam) * w_K     # 用 ω_log 代 ω̄₂ ⟹ 天花板的**下**估
        f_up = tc / ceil_lo
        verd = "PASS（天花板成立）" if f_up <= 1.0 else "不判定（需真 ω̄₂）"
        print(f"{tag:<8}{tc:8.1f}{lam:8.3f}{wl:12.1f}{w_K:11.1f}{ceil_lo:12.1f}{f_up:9.4f}   {verd}")

print()
print("   注：天花板用 ω_log 代 ω̄₂ 是**下**估（ceiling_lo ≤ 真天花板）")
print("       ⟹ 『Tc ≤ ceiling_lo』蕴含『Tc ≤ 真天花板』⟹ 只支持 PASS。方向正确。")

print()
print("=" * 92)
print("C. 落在 AD75 29 行的什么位置？（𝔣 ∈ [0.1994, 0.6102]，均值 0.3587）")
print("=" * 92)
for tag, tc, lam, wl in [("MgB2", 39.4, 0.87, 62.0), ("CaC6", 11.5, 0.870, 25.0)]:
    f_up = tc / (C_INF * math.sqrt(lam) * wl / KB)
    print(f"   {tag}: 𝔣 ≤ {f_up:.4f}"
          f"   ⟹ {'落在 AD75 区间内' if 0.1994 <= f_up <= 0.6102 else '区间外'}"
          f"（真值还要更低，因为 ω̄₂ > ω_log）")
print()
print("   ⟹ 两行都**不是**反例，而是与 29 行常规超导体同一量级的普通点。")
print("      𝔣>1 本来就该当场触发警报：天花板是定理，不是拟合。")
print("=" * 92)
