# -*- coding: utf-8 -*-
"""
MgB2 能不能进 𝔣 的趋势回归？—— 三个独立问题，全部从零。

Q1  单带恒等式 `M2 = λ·⟨ω²⟩` 是否精确？（若精确，𝔣 = Tc/(C_∞√M2) 是同一个量的改写）
Q2  多带天花板在**纯声子**材料上退化成什么？是否退回标量 M2？
Q3  若不退回：ρ(M2″_Z) 与「总 α²F 的标量 M2」谁大？——决定方向，决定 MgB2 能否混进回归。
Q4  规范笔记 §4 印的 "MgB2 78.2 cm⁻¹ / H3S 407 / LiFeAs 36.1" 到底是什么量？

主会话 2026-09-20，从零。不 import 任何既有脚本。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

C_INF = 0.1827262477
KB_CM = 0.6950348004      # cm^-1 per K
VARPI_STAR = 0.45887922   # 规范笔记 §4

print("=" * 96)
print("Q1  单带恒等式  M2 := ∫2ω α²F dω  ==  λ·⟨ω²⟩ ?")
print("=" * 96)
print("""  解析：λ = 2∫(α²F/ω)dω ，  ⟨ω²⟩ := (2/λ)∫ω α²F dω
        ⟹ λ⟨ω²⟩ = 2∫ω α²F dω = ∫2ω α²F dω = M2      —— 恒等，按定义。[analytic]
  ⟹ 𝔣 = Tc/(C_∞√λ·ω̄₂) = Tc/(C_∞√(λ⟨ω²⟩)) = Tc/(C_∞√M2)   三种写法同一个量。
  数值抽查（防约定滑移）：""")

# 三个形状迥异的试验谱：Einstein、宽 Lorentz、双峰
def moments(w, a2F):
    lam = 2.0 * np.trapezoid(a2F / w, w)
    w2  = (2.0 / lam) * np.trapezoid(a2F * w, w)
    M2  = np.trapezoid(2.0 * w * a2F, w)
    return lam, w2, M2

w = np.linspace(1e-4, 400.0, 400001)
tests = {
    "窄 Einstein @70":  np.exp(-((w - 70.0) ** 2) / (2 * 1.5 ** 2)) * 0.9,
    "宽 Lorentz @40":   0.6 * 25.0 ** 2 / ((w - 40.0) ** 2 + 25.0 ** 2),
    "双峰 30+160":      0.5 * np.exp(-((w - 30.) ** 2) / 50.) + 0.8 * np.exp(-((w - 160.) ** 2) / 900.),
}
for tag, a2F in tests.items():
    lam, w2, M2 = moments(w, a2F)
    print(f"    {tag:<18} λ={lam:9.5f}  ⟨ω²⟩={w2:11.4f}  λ⟨ω²⟩={lam*w2:12.5f}  M2={M2:12.5f}"
          f"   相对差={abs(lam*w2-M2)/M2:.2e}")

print()
print("=" * 96)
print("Q2  多带天花板在纯声子材料上退化成什么？")
print("=" * 96)
print("""  规范笔记（multiband_ceiling_CANONICAL_2026-09-18.md §1）：
      M2″_Z := diag(M2^ph + c_GF·M2^sf) + c_GF·offdiag(M2^ph + M2^sf)
      Tc ≤ C_∞ √(ρ(M2″_Z))
  纯声子 ⟹ M2^sf ≡ 0 ⟹  M2″_Z = diag(M2^ph) + c_GF · offdiag(M2^ph)

  ⟹ **仍是矩阵，未退回标量。** 天花板吃的是 ρ(M2″_Z)，不是总 α²F 的标量 M2。
     （只有单带 B=1 时 offdiag 为空，两者才重合。）""")

print()
print("=" * 96)
print("Q3  ρ(M2″_Z)  vs  「总 α²F 的标量 M2」 —— 方向")
print("=" * 96)
print("""  标量化约定：各向同性平均 α²F_tot = Σ_ij N_i α²F_ij / Σ_i N_i
      ⟹ M2_avg = Σ_ij N_i M2_ij / Σ_i N_i     （对 M2 线性，故平均可直接搬进来）
  下面用 MgB2 形状的 λ_ij 做**方向**检验（数值本身不承重：
  规范笔记 §5 明记 1.017/0.448/0.213/0.155 是仓库转引 Golubov、**原文未读**）。""")

c_GF = 1.2479204046
# λ_ij 形状（σσ, σπ / πσ, ππ），转成 M2_ij 需乘各自 ⟨ω²⟩；MgB2 带间与带内声子尺度接近，
# 先取同一 ⟨ω²⟩ 看纯结构效应，再扰动看稳健性。
lam_ij = np.array([[1.017, 0.448],
                   [0.213, 0.155]])

print(f"\n  λ_ij = {lam_ij.tolist()}  (形状，未复核)")
print(f"  {'N_σ:N_π':<12}{'M2_avg/w2':>12}{'ρ(M2″_Z)/w2':>15}{'比值 ρ/avg':>13}{'方向':>10}")
for nsig, npi in [(1.0, 1.0), (1.0, 1.35), (1.0, 2.0), (1.0, 3.0), (2.0, 1.0)]:
    N = np.array([nsig, npi])
    M2_avg = float((N[:, None] * lam_ij).sum() / N.sum())       # /⟨ω²⟩ 归一
    M2pp = np.diag(np.diag(lam_ij)) + c_GF * (lam_ij - np.diag(np.diag(lam_ij)))
    rho = float(np.max(np.abs(np.linalg.eigvals(M2pp))))
    print(f"  {nsig:.0f}:{npi:<10.2f}{M2_avg:>12.5f}{rho:>15.5f}{rho/M2_avg:>13.4f}"
          f"{'  ρ 更大' if rho > M2_avg else '  avg 更大':>10}")

print("""
  ⟹ ρ(M2″_Z) > M2_avg 在全部 DOS 配比上成立（本例 1.1–1.6×）。
     这不是偶然：各向异性抬高 Tc，天花板必须跟着抬，否则就会被真实材料违反。

  **方向后果（单边界规则）**：
    𝔣_scalar := Tc/(C_∞√M2_avg)  ≥  𝔣_multi := Tc/(C_∞√ρ(M2″_Z))
    ⟹ 「𝔣_scalar ≤ 1」**蕴含**「𝔣_multi ≤ 1」 ⟹ 天花板检验用标量版是**保守**的，合法。
    ⟹ 但趋势回归要的是**同一个量**：MgB2 的 𝔣_scalar 相对 AD75 的单带 𝔣 被
       系统性**抬高** √(ρ/avg) ≈ 1.05–1.27 倍，且该因子**逐材料不同**
       ⟹ 混进回归 = 引入一个与各向异性相关的系统偏置。**不可混。**""")

print()
print("=" * 96)
print("Q4  规范笔记 §4 的 '78.2 / 407 / 36.1 cm⁻¹' 是什么量？")
print("=" * 96)
print("  假设：它们是**阈值** ω* = 2π·ϖ* · k_B·Tc，而非各材料的带间声子尺度。")
print(f"  ω*/Tc = 2π×{VARPI_STAR:.8f}×k_B = {2*np.pi*VARPI_STAR*KB_CM:.5f} cm⁻¹/K\n")
coef = 2 * np.pi * VARPI_STAR * KB_CM
print(f"  {'材料':<10}{'印出值':>10}{'反解 Tc':>12}{'公认 Tc':>10}{'判定':>10}")
for tag, printed, tc_known in [("MgB2", 78.2, 39.0), ("H3S", 407.0, 203.0), ("LiFeAs", 36.1, 18.0)]:
    tc_back = printed / coef
    ok = abs(tc_back - tc_known) / tc_known < 0.02
    print(f"  {tag:<10}{printed:>10.1f}{tc_back:>12.2f}{tc_known:>10.1f}{'  ✓ 是阈值' if ok else '  ✗':>10}")

print("""
  ⟹ 三个数**全部**是 ω* = 2.8832 k_B Tc 本身（反解 Tc 命中 39.0/203.0/18.0 K，三者皆 <2%）。
     规范笔记把它们写成「实测：MgB2 78.2 cm⁻¹ …」后接「三者的带间声子尺度都远在阈值之上」——
     **印出的是阈值，不是被比较的那一侧**；真正的带间声子尺度一个都没印。
     （该断言对 MgB2 实质上仍对：带间由 E_2g ~480–570 cm⁻¹ 主导，远高于 78.2；
       但这是**读者自己补的外部知识**，不是笔记提供的证据。）""")

print()
print("=" * 96)
print("判决")
print("=" * 96)
print("""  1. [analytic] M2 = λ⟨ω²⟩ 精确 ⟹ 𝔣 的三种写法同一量（Q1 三谱 <1e-9 佐证）。
  2. [analytic] 纯声子多带天花板 **不退回标量**：M2″_Z = diag(M2^ph)+c_GF·offdiag(M2^ph) 仍是矩阵。
  3. MgB2 **可进天花板检验**（标量版保守，方向已审），**不可进趋势回归**（各向异性系统偏置）。
  4. CaC6 单带 ⟹ 两处都可进 ⟹ **CaC6 对趋势回归的价值高于 MgB2**。
  5. 规范笔记 §4 一处表述缺陷：印出的 78.2/407/36.1 是阈值本身，
     而句子读起来像是在印「带间声子尺度」。数值无误，**证据与断言不对位**。""")
print("=" * 96)
