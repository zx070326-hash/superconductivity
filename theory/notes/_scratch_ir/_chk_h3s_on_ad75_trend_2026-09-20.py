# -*- coding: utf-8 -*-
"""
把 H3S/D3S 放到 AD75 Table I 的 29 材料 (λ, ω_log/ω̄₂) → 𝔣 趋势上。
（结掉 `materials/deficit_factor_ad75_tunneling_sample_2026-09-20.md` §6 的 OPEN。）

主会话 2026-09-20，从零：
  · 29 行 Table I 我本人再次逐行对页图 `_ad75_p9-09.png` 核过（29×9 全中）；
  · H3S/D3S 的 (λ, ω_log, ω̄₂) 由我直接积分 Mishra `.a2f`，不转引 CSV、不转引 agent；
  · 回归我自己重拟（最小二乘的正规方程），不转引既有脚本的系数。

约定对齐（承重）：
  AD75 的 ω̄₂ := [(2/λ)∫α²F ω dω]^{1/2}、ω_log := exp[(2/λ)∫(α²F/ω)ln ω dω]，
  与本仓库 / 本脚本对 .a2f 的定义**逐字相同** ⟹ 可并表。
  单位：AD75 列为 K，Mishra 谱为 meV ⟹ 除以 k_B 换算。

方向注记：本脚本不产生任何单边界判决（回归是双侧的描述性拟合）。
配套笔记：`../materials/h3s_on_ad75_trend_closed_2026-09-20.md`
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

C_INF = 0.1827262477          # 无量纲；P4_g2_analytic.md CERT-B 下端
KB    = 0.08617333262         # meV / K

# ---- AD75 Table I（页图逐行转录；本次由我第二次独立对图核验，29×9 全中）----
# (material, w_log[K], wbar1[K], wbar2[K], lambda, w_ph[K], mu*, Tc[K], eta)
T1 = [
    ("Pb",              56,  60,  65, 1.55, 110, 0.105, 7.2,  2.4),
    ("In",              68,  79,  89, 0.805,179, 0.097, 3.40, 1.3),
    ("Sn",              99, 110, 121, 0.72, 209, 0.092, 3.75, 2.2),
    ("Hg",              29,  38,  49, 1.6,  162, 0.098, 4.19, 1.4),
    ("Tl",              52,  58,  64, 0.795,127, 0.111, 2.36, 1.2),
    ("Ta",             132, 140, 148, 0.69, 228, 0.093, 4.48, 4.9),
    ("a-Ga",            55,  77, 101, 1.62, 291, 0.095, 8.56, 2.1),
    ("b-Ga",            87, 108, 129, 0.97, 285, 0.092, 5.90, 2.0),
    ("Tl0.9Bi0.1",      48,  55,  62, 0.78, 120, 0.099, 2.30, 1.1),
    ("Pb0.4Tl0.6",      48,  56,  62, 1.15, 121, 0.094, 4.60, 1.6),
    ("Pb0.6Tl0.4",      50,  57,  62, 1.38, 119, 0.103, 5.90, 2.0),
    ("Pb0.8Tl0.2",      50,  56,  61, 1.53, 116, 0.101, 6.80, 2.1),
    ("Pb0.6Tl0.2Bi0.2", 48,  53,  58, 1.81, 112, 0.111, 7.26, 2.3),
    ("Pb0.9Bi0.1",      50,  56,  60, 1.66, 108, 0.081, 7.65, 2.2),
    ("Pb0.8Bi0.2",      46,  52,  57, 1.88, 109, 0.093, 7.95, 2.3),
    ("Pb0.7Bi0.3",      47,  52,  57, 2.01, 110, 0.092, 8.45, 2.4),
    ("Pb0.65Bi0.35",    45,  50,  55, 2.13, 110, 0.093, 8.95, 2.4),
    ("In0.9Tl0.1",      63,  75,  86, 0.85, 176, 0.103, 3.28, 1.4),
    ("In0.73Tl0.27",    55,  67,  77, 0.93, 166, 0.110, 3.36, 1.4),
    ("In0.67Tl0.33",    57,  68,  79, 0.90, 167, 0.110, 3.26, 1.4),
    ("In0.57Tl0.43",    53,  64,  74, 0.85, 165, 0.117, 2.60, 1.3),
    ("In0.5Tl0.5",      53,  64,  73, 0.83, 163, 0.110, 2.52, 1.3),
    ("In0.27Tl0.73",    42,  53,  63, 1.09, 151, 0.094, 3.64, 1.4),
    ("In0.17Tl0.83",    45,  55,  63, 0.98, 144, 0.101, 3.19, 1.3),
    ("In0.07Tl0.93",    49,  56,  63, 0.89, 131, 0.107, 2.77, 1.3),
    ("In2Bi",           46,  57,  67, 1.40, 174, 0.096, 5.6,  1.6),
    ("Sb2Tl7",          37,  48,  58, 1.43, 134, 0.102, 5.2,  1.6),
    ("Bi2Tl",           47,  53,  59, 1.63, 120, 0.101, 6.4,  2.1),
    ("a-Pb0.45Bi0.55",  29,  38,  47, 2.59, 128, 0.116, 7.0,  2.1),
]

name = [r[0] for r in T1]
wlog = np.array([r[1] for r in T1], float)      # K
wb2  = np.array([r[3] for r in T1], float)      # K
lam  = np.array([r[4] for r in T1], float)
tc   = np.array([r[7] for r in T1], float)      # K

f_ad   = tc / (C_INF * np.sqrt(lam) * wb2)      # 无量纲
shape  = wlog / wb2

# ---- Mishra α²F：H3S / D3S 的三元组，我自己积分 ----
BASE = os.path.join(r"D:\physics\superconductivity", "theory", "data", "external",
                    "mishra_h3s_d3s", "Data-H3S-D3S-PH-EPW-Anharmonic")
tr = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def moments(rel):
    rows = []
    for line in open(os.path.join(BASE, rel.replace("/", os.sep)), errors="replace"):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        try:
            v = [float(x) for x in s.split()]
        except ValueError:
            continue
        if len(v) >= 2:
            rows.append(v[:2])
    a = np.array(rows)
    w, F = a[:, 0], np.clip(a[:, 1], 0.0, None)
    m = w > 1e-9
    w, F = w[m], F[m]
    L  = 2 * tr(F / w, w)
    wl = np.exp((2 / L) * tr(F / w * np.log(w), w))      # meV
    w2 = np.sqrt((2 / L) * tr(F * w, w))                 # meV
    return L, wl / KB, w2 / KB                           # λ, ω_log[K], ω̄₂[K]

# Tc_exp：压强匹配的干净锚点（Drozdov Nature 2015 Fig.2c annealed @200 GPa）
#   H3S 184.04 K（主会话双路矢量数字化）；D3S 145.5 K（单路，承重前需第二路）
NEW = [("H3S-200anh", "H3S-200GPa/EPW_anharmonic/H3S.a2f", 184.04, "双路"),
       ("D3S-200anh", "D3S-200GPa/EPW_anharmonic/D3S.a2f", 145.5,  "单路"),
       ("H3S-200har", "H3S-200GPa/EPW_harmonic/H3S.a2f",   184.04, "双路"),
       ("D3S-200har", "D3S-200GPa/EPW_harmonic/D3S.a2f",   145.5,  "单路")]

print("=" * 96)
print("A. AD75 29 行：λ 与 shape:=ω_log/ω̄₂ 的覆盖范围（判断 H3S/D3S 是内插还是外推）")
print("=" * 96)
print(f"   λ     ∈ [{lam.min():.3f}, {lam.max():.3f}]")
print(f"   shape ∈ [{shape.min():.4f}, {shape.max():.4f}]   (min={name[int(shape.argmin())]}, "
      f"max={name[int(shape.argmax())]})")
print(f"   𝔣     ∈ [{f_ad.min():.4f}, {f_ad.max():.4f}]   mean {f_ad.mean():.4f}")
print(f"   shape ≤ 1 全成立？ {bool(np.all(shape <= 1))}   （ω_log ≤ ω̄₂，幂平均不等式）")

# ---- 我自己重拟回归（正规方程，不调 polyfit/lstsq 的默认）----
X = np.column_stack([np.ones_like(lam), np.log(lam), np.log(shape)])
y = np.log(f_ad)
XtX = X.T @ X
beta = np.linalg.solve(XtX, X.T @ y)
res  = y - X @ beta
n, p = X.shape
s2   = res @ res / (n - p)
ss_tot = ((y - y.mean()) ** 2).sum()
R2   = 1 - (res @ res) / ss_tot
XtXi = np.linalg.inv(XtX)
se_beta = np.sqrt(np.diag(XtXi) * s2)

print()
print("=" * 96)
print("B. 回归（我从零重拟，正规方程）： log 𝔣 = b0 + b1·log λ + b2·log(ω_log/ω̄₂)")
print("=" * 96)
for k, lbl in enumerate(["b0 (截距)", "b1 (log λ)", "b2 (log shape)"]):
    print(f"   {lbl:16s} = {beta[k]:+.4f}  ± {se_beta[k]:.4f}")
print(f"   R² = {R2:.4f}     残差 RMS(log) = {np.sqrt(res @ res / n):.4f}"
      f"  ⟹ 乘性散布 {(np.exp(np.sqrt(res @ res / n)) - 1) * 100:.2f}%")
print(f"   s (无偏, n−p) = {np.sqrt(s2):.4f}")

print()
print("=" * 96)
print("C. H3S / D3S 落在哪里")
print("=" * 96)
print(f"{'体系':<12}{'λ':>9}{'ω_log[K]':>11}{'ω̄₂[K]':>11}{'shape':>8}{'天花板[K]':>11}"
      f"{'𝔣_obs':>9}{'𝔣_pred':>9}{'偏离':>9}{'z(pred)':>9}")
for tag, rel, tcx, route in NEW:
    L, wl_K, w2_K = moments(rel)
    ceil = C_INF * np.sqrt(L) * w2_K
    sh   = wl_K / w2_K
    fobs = tcx / ceil
    x    = np.array([1.0, np.log(L), np.log(sh)])
    fpr  = np.exp(x @ beta)
    # 预测区间的标准误（含杠杆）：新点在 log 空间
    se_pred = np.sqrt(s2 * (1 + x @ XtXi @ x))
    z = (np.log(fobs) - np.log(fpr)) / se_pred
    inside = "内插" if (lam.min() <= L <= lam.max() and shape.min() <= sh <= shape.max()) else "★外推"
    print(f"{tag:<12}{L:9.4f}{wl_K:11.2f}{w2_K:11.2f}{sh:8.4f}{ceil:11.2f}"
          f"{fobs:9.4f}{fpr:9.4f}{(fobs/fpr-1)*100:+8.2f}%{z:9.2f}   [{inside}, Tc {route}]")

print()
print("   杠杆说明：z = (log𝔣_obs − log𝔣_pred)/se_pred，se_pred 含 (1 + x'(X'X)⁻¹x) 杠杆项。")
print("   |z| < 2 ⟹ 该体系与 29 材料趋势在 95% 预测带内**不可区分**。")

# ---- 敏感性：Tc_exp 锚点换成 legacy 值会怎样 ----
print()
print("=" * 96)
print("D. 对 Tc_exp 锚点的敏感性（legacy 203/150 K vs 压强匹配 184.04/145.5 K）")
print("=" * 96)
for tag, rel, tc_new, tc_old in [("H3S-200anh", "H3S-200GPa/EPW_anharmonic/H3S.a2f", 184.04, 203.0),
                                 ("D3S-200anh", "D3S-200GPa/EPW_anharmonic/D3S.a2f", 145.5, 150.0)]:
    L, wl_K, w2_K = moments(rel)
    ceil = C_INF * np.sqrt(L) * w2_K
    sh = wl_K / w2_K
    x = np.array([1.0, np.log(L), np.log(sh)])
    fpr = np.exp(x @ beta)
    se_pred = np.sqrt(s2 * (1 + x @ XtXi @ x))
    for lbl, t in (("压强匹配", tc_new), ("legacy  ", tc_old)):
        fo = t / ceil
        print(f"   {tag}  {lbl}  Tc={t:7.2f} K  𝔣_obs={fo:.4f}  "
              f"偏离 𝔣_pred={fpr:.4f} 为 {(fo/fpr-1)*100:+7.2f}%   z={(np.log(fo)-np.log(fpr))/se_pred:+.2f}")
print("=" * 96)

# ---------------------------------------------------------------------------
# E. 残差的 μ* 自审：AD75 29 行的 μ* ∈ [0.081,0.117]（均值 ~0.10），而 H3S 文献取 0.13–0.16。
#    𝔣 的**定义式不含 μ***，但**回归是在 μ*≈0.10 的人群上标定的** ⟹ 把它外推到
#    μ* 更大的体系，必然**高估** 𝔣。先问：观测到的 −6~7% 残差，够不够被这件事吃掉？
#    工具：McMillan 指数 −1.04(1+λ)/(λ−μ*(1+0.62λ)) 的 μ* 灵敏度（**AD 族内部**，非独立证据）。
# ---------------------------------------------------------------------------
def mcm_exponent(L, mus):
    return -1.04 * (1 + L) / (L - mus * (1 + 0.62 * L))


print()
print("=" * 96)
print("E. 残差能否由 μ* 偏移解释？（AD75 人群 μ*≈0.10，氢化物文献 0.13–0.16）")
print("=" * 96)
mus_ad = np.array([r[6] for r in T1])
print(f"   AD75 29 行 μ*: min {mus_ad.min():.3f}  mean {mus_ad.mean():.4f}  max {mus_ad.max():.3f}"
      f"  (跨度仅 {mus_ad.max()-mus_ad.min():.3f})")
print()
print(f"{'体系':<12}{'λ':>8}{'残差_obs':>11}   " + "".join(f"μ*={m:<8.2f}" for m in (0.13, 0.16)))
for tag, rel, tcx, _ in NEW[:2]:
    L, wl_K, w2_K = moments(rel)
    sh = wl_K / w2_K
    fobs = tcx / (C_INF * np.sqrt(L) * w2_K)
    fpr = np.exp(np.array([1.0, np.log(L), np.log(sh)]) @ beta)
    line = f"{tag:<12}{L:8.4f}{(fobs/fpr-1)*100:+10.2f}%   "
    for m in (0.13, 0.16):
        d = np.exp(mcm_exponent(L, m) - mcm_exponent(L, mus_ad.mean())) - 1
        line += f"{d*100:+9.2f}%   "
    print(line)
print("   （右侧两列 = 「若真实 μ* 是 0.13 / 0.16，相对 AD75 人群均值 μ*=0.0997 应有的 𝔣 下压」）")
print("   ⚠ 这是 **AD 族内部**的灵敏度估计，不是对残差的独立解释；只回答「量级够不够」。")

# 第三回归量：μ* 在 29 行上本身可辨吗？（若不可辨，上面的外推就没有样本支撑）
X3 = np.column_stack([X, mus_ad])
b3 = np.linalg.solve(X3.T @ X3, X3.T @ y)
r3 = y - X3 @ b3
R2_3 = 1 - (r3 @ r3) / ss_tot
s2_3 = r3 @ r3 / (n - 4)
se3 = np.sqrt(np.diag(np.linalg.inv(X3.T @ X3)) * s2_3)
print()
print("   加入 μ* 作第三回归量（检验：29 行本身能不能看见 μ*）：")
print(f"      系数(μ*) = {b3[3]:+.3f} ± {se3[3]:.3f}   t = {b3[3]/se3[3]:+.2f}   "
      f"R²: {R2:.4f} → {R2_3:.4f}")
dl_dmu = (mcm_exponent(1.5, 0.11) - mcm_exponent(1.5, 0.09)) / 0.02
print(f"      对照：McMillan 在 λ=1.5 处给 dlog Tc/dμ* = {dl_dmu:+.2f}")
print("   ⟹ 若 t 不显著、且 29 行 μ* 跨度只有 0.036，则**该样本对 μ* 基本盲**，")
print("      外推到 μ*=0.13–0.16 属于样本外，必须显式标注（不能只报 |z|<2）。")
print("=" * 96)
