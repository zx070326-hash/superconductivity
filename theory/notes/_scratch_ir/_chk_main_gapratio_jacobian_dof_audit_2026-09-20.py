# -*- coding: utf-8 -*-
"""核实 gapratio-noise 的 Jacobian 比较 —— 先做代数，再做独立仪器交叉。

它的交付（脚本 `gap_ratio_noise_weighted_v2.py`，锚点 λ=1.876, ω_E=106.9 meV, μ*=0.13,
w_cut=500 meV, N=1600, T_low=0.03Tc）：
    J_A = [[+93.854, -385.263], [+0.738, -2.068]]      观测量 {Tc, 2Δ/kTc}
    J_B = [[+93.854, -385.263], [+0.006, -0.154]]      观测量 {Tc, α_isotope}
    误差模型 1（统计）：面积比 A/B = 0.4978，翻转阈 σ(2Δ/kTc)/σ(Tc) = 0.1974
    误差模型 2（总）  ：面积比 A/B = 0.6224，翻转阈 σ(2Δ/kTc)/σ(Tc) = 0.0295

§1  ★ 自由度审计：把「两个阈值」当两个检验之前先数自由度（CLAUDE.md 新规则）。
§2  ★ ∂α/∂λ 是最小、最不可靠的一项，det(J_B) 对它敏感 —— 用 AD 闭式作独立仪器交叉 + 扫描。
§3  v1 脚本「仅 Tc」那一栏是不是奇异矩阵赝物。

方向注记：§1/§3 是恒等式（analytic）；§2 的 AD 交叉是 numerical-evidence，
且**方向已审**：见 §2 末尾。
主会话 2026-09-20，从零。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

J_A = np.array([[+93.854, -385.263], [+0.738, -2.068]])
J_B = np.array([[+93.854, -385.263], [+0.006, -0.154]])
TC0, GR0, AL0 = 188.58, 4.629, 0.4973

print("=" * 94)
print("§1  自由度审计：两个「翻转阈值」其实是同一个数")
print("=" * 94)


def area(J, sig):
    """置信区面积 ∝ 1/√det(JᵀΣ⁻¹J)。对 2×2 满秩 J 恒等于 σ₁σ₂/|det J|。"""
    F = J.T @ np.diag(1.0 / np.asarray(sig) ** 2) @ J
    return 1.0 / np.sqrt(np.linalg.det(F))


dA, dB = np.linalg.det(J_A), np.linalg.det(J_B)
print(f"  det(J_A) = {dA:+.6f}   （逐项：93.854·(−2.068) = {93.854*-2.068:+.4f}，"
      f"−(−385.263)·0.738 = {385.263*0.738:+.4f}）")
print(f"  det(J_B) = {dB:+.6f}   （逐项：93.854·(−0.154) = {93.854*-0.154:+.4f}，"
      f"−(−385.263)·0.006 = {385.263*0.006:+.4f}）")

# 恒等式 area = σ₁σ₂/|det J|，数值验证
for nm, J, d in (("A", J_A, dA), ("B", J_B, dB)):
    s = (np.pi, np.e)                       # 任取两个 σ，验证闭式
    lhs, rhs = area(J, s), s[0] * s[1] / abs(d)
    assert abs(lhs - rhs) / rhs < 1e-10, (nm, lhs, rhs)
print("  恒等式 Area = σ₁σ₂/|det J| 验证通过（1e-10；det 里有 32%/16% 的相消，"
      "吃掉几位有效数字，故非 1e-15）[analytic]")

print()
print("  ⟹ Area_A/Area_B = (σ_Tc·σ_gap/|det J_A|) / (σ_Tc·σ_α/|det J_B|)")
print("                   = (σ_gap/σ_α) · |det J_B|/|det J_A|")
print("     ★ σ_Tc **逐位相消**。判据里只剩一个自由度：")
THR = abs(dA) / abs(dB)
print(f"        A 胜 ⟺ σ_gap/σ_α < |det J_A/det J_B| = {THR:.4f}")

# 数值证明 σ_Tc 相消：扫 σ_Tc 三个数量级，面积比不动
print()
print("  数值佐证（扫 σ_Tc 三个数量级，固定 σ_gap=0.1852, σ_α=0.0500）：")
sg, sa = 0.04 * GR0, 0.05
for st in (0.1, 1.0, 10.0, 100.0):
    r = area(J_A, (st, sg)) / area(J_B, (st, sa))
    print(f"      σ_Tc = {st:>7.2f} K  →  Area_A/Area_B = {r:.10f}")
print("      ⟹ 逐位不变，σ_Tc 确实不进判据。")

print()
print("  复现它印出的两组数，并换算到唯一自由度：")
MODELS = {
    "模型1 统计 (σTc=1%, σgr=4%, σα=0.05 abs)": (0.01 * TC0, 0.04 * GR0, 0.05),
    "模型2 总   (σTc=10%, σgr=7.5%, σα=15% rel)": (0.10 * TC0, 0.075 * GR0, 0.15 * AL0),
}
print(f"  {'误差模型':<44}{'面积比A/B':>11}{'σgap/σα':>11}{'余量×':>9}{'它报的σgap/σTc':>15}")
for nm, (st, sg, sa) in MODELS.items():
    r = area(J_A, (st, sg)) / area(J_B, (st, sa))
    print(f"  {nm:<44}{r:>11.4f}{sg/sa:>11.4f}{THR/(sg/sa):>9.3f}{sg/st:>15.4f}")
print()
print("  ★ 两行的『余量×』= 1/面积比，和它印的 2.01/1.60 一致；")
print("    但它报出的两个「翻转阈值」0.1974 与 0.0295（差 6.7×）**只因 σ_Tc 差 10×**，")
print("    而 σ_Tc 已证明不进判据 ⟹ 那两个数不是两个独立检验，是同一句话的两次改写。")
print(f"    唯一该报的数：σ_gap/σ_α < {THR:.3f}（对两个误差模型**逐位相同**）。")
print("    这正是 CLAUDE.md「把一组数字当作 N 个独立检验之前，先数它的自由度」。")
print("  ⚠ 注意方向：本条**加强**它的结论（判据与误差模型无关），只订正表述。")

print()
print("=" * 94)
print("§2  ∂α/∂λ：最小的一项，det(J_B) 对它敏感 —— 独立仪器交叉 + 扫描")
print("=" * 94)


def alpha_ad(lam, mu):
    """AD/McMillan 闭式同位素指数（ω_E ∝ M^{-1/2}，Tc = (ω_log/1.2)exp(−1.04(1+λ)/(λ−μ*(1+0.62λ)))）。
    α = ½[1 − μ*²·1.04(1+λ)(1+0.62λ)/(λ−μ*(1+0.62λ))²]  —— 解析求导，非差分。"""
    d = lam - mu * (1 + 0.62 * lam)
    return 0.5 * (1.0 - mu ** 2 * 1.04 * (1 + lam) * (1 + 0.62 * lam) / d ** 2)


LAM, MU, DLAM, DMU = 1.876, 0.13, 0.15 * 1.876, 0.05
g_ad = (alpha_ad(LAM + DLAM, MU) - alpha_ad(LAM - DLAM, MU)) / (2 * DLAM)
h_ad = (alpha_ad(LAM, MU + DMU) - alpha_ad(LAM, MU - DMU)) / (2 * DMU)
# 解析导数（复步长），排除 AD 侧自身的差分误差
eps = 1e-7
g_ad_x = (alpha_ad(LAM + eps, MU) - alpha_ad(LAM - eps, MU)) / (2 * eps)
h_ad_x = (alpha_ad(LAM, MU + eps) - alpha_ad(LAM, MU - eps)) / (2 * eps)
print(f"  AD 闭式 α(λ={LAM}, μ*={MU}) = {alpha_ad(LAM, MU):.6f}   （它的 Eliashberg 值 {AL0}）")
print(f"  AD ∂α/∂λ  ：同步长差分 {g_ad:+.6f}  / 解析 {g_ad_x:+.6f}   （它的 {J_B[1,0]:+.6f}）")
print(f"  AD ∂α/∂μ* ：同步长差分 {h_ad:+.6f}  / 解析 {h_ad_x:+.6f}   （它的 {J_B[1,1]:+.6f}）")
print("  ⟹ ∂α/∂λ 两把仪器**同号同量级、相差约 2×**；∂α/∂μ* 相差约 2× 亦同号。")
print("     AD 不是 Eliashberg，2× 不判谁错；它判的是：这一项**小到两把仪器分不清**。")

print()
print("  det(J_B) = 93.854·(∂α/∂μ*) + 385.263·(∂α/∂λ)  ⟹ 扫 ∂α/∂λ 看阈值怎么动：")
print(f"  {'∂α/∂λ':>10}{'det(J_B)':>12}{'阈值 σgap/σα':>15}{'模型1余量':>11}{'模型2余量':>11}   判决")
r1 = (0.04 * GR0) / 0.05
r2 = (0.075 * GR0) / (0.15 * AL0)
for g in (-0.020, -0.010, 0.000, 0.006, 0.0116, 0.020, 0.030):
    d = 93.854 * J_B[1, 1] + 385.263 * g
    t = abs(dA) / abs(d) if d != 0 else np.inf
    v = "A 胜(两模型)" if t > max(r1, r2) else ("A 仅胜模型1" if t > min(r1, r2) else "★ A 败")
    mark = " ←它的值" if g == 0.006 else (" ←AD 交叉" if g == 0.0116 else "")
    print(f"  {g:>10.4f}{d:>12.4f}{t:>15.4f}{t/r1:>11.3f}{t/r2:>11.3f}   {v}{mark}")

print()
print(f"  实测 σ_gap/σ_α：模型1 = {r1:.4f}，模型2 = {r2:.4f}")
print("  ★ 方向审计：det(J_B) 的两项**同号相加**（∂α/∂μ*<0，−(−385)·∂α/∂λ 在 ∂α/∂λ>0 时为正，")
print("    故 |det J_B| 随 ∂α/∂λ 增大而**减小**）⟹ 阈值随 ∂α/∂λ **增大**。")
print("    AD 给的 +0.0116 > 它的 +0.006 ⟹ **朝 A 更胜的方向**。")
print("    ⟹ 「A 胜」这个判决对 ∂α/∂λ 的不确定性是**安全**的；不安全的是阈值的**点值**")
print(f"       （∂α/∂λ ∈ [0, 0.012] 就让阈值在 [{abs(dA)/abs(93.854*J_B[1,1]):.2f}, "
      f"{abs(dA)/abs(93.854*J_B[1,1]+385.263*0.012):.2f}] 之间走）。")
print(f"    翻盘需 ∂α/∂λ ≲ {(abs(dA)/r2*np.sign(-1) - 93.854*J_B[1,1])/385.263:+.4f}"
      f"（模型2）——与两把仪器给的 +0.006/+0.0116 反号 ⟹ 未触及。")

print()
print("=" * 94)
print("§3  v1 脚本「仅 Tc」一栏：1 个观测量定不了 2 个参数")
print("=" * 94)
J1 = J_A[0:1, :]
F1 = J1.T @ np.diag([1.0 / 2.0 ** 2]) @ J1
print(f"  J(1×2) 的 Fisher 阵 2×2，秩 = {np.linalg.matrix_rank(F1)}，"
      f"det = {np.linalg.det(F1):.3e}（理论上恒为 0）")
print("  ⟹ 置信区是**无界条带**，面积 = ∞。v1 `_chk_gapratio_noise_weighted_2026-09-20.py`")
print("     对它算 1/√det 得到有限数，那个数是浮点残差的赝物，不是面积。")
print("     v2 改成 A/B 两个满秩 2×2 对比后此问题消失 ⟹ **只影响 v1，v2 的判决不受影响**。")
print("     （记录在案：若有笔记引用过 v1 的「仅 Tc」数字，须撤。）")
print("=" * 94)
