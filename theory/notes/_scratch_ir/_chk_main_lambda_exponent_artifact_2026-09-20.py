# -*- coding: utf-8 -*-
"""
主会话从零复核：`deficit_decomposition_analysis_2026-09-20.md` 任务 A
「λ 指数 ≈1 是 AD 公式伪影」这条判决。

它的证据是：AD 公式的**逐点数值导数** ∂log𝔣/∂logλ 均值 1.092，与回归指数 0.997
「相差 9.5% ⟹ 支持伪影假说」。

本脚本做它**没做**的那个判决性检验：
  把 AD 公式当成「假设 AD 就是真理」的世界，用**同样的 29 组 (λ, ω_log, ω̄₂, μ*)**
  生成 𝔣_AD，跑**同一个回归**，比较 (系数, R², 残差散布)。
  —— 逐点导数的均值不是回归系数的预测；只有重跑回归才是。

数据：AD1975 Table I，29 行，页图逐行转录（见 _chk_ad75_table1_frak_f_2026-09-20.py）。
"""
import numpy as np

C_INF = 0.1827262477

# (material, w_log(K), wbar1, wbar2(K), lambda, w_ph, mu*, Tc(K), eta)
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

wlog = np.array([r[1] for r in T1], float)
wb2  = np.array([r[3] for r in T1], float)
lam  = np.array([r[4] for r in T1], float)
mus  = np.array([r[6] for r in T1], float)
Tc   = np.array([r[7] for r in T1], float)
name = [r[0] for r in T1]


def allen_dynes(lam, wlog, wb2, mus):
    """AD75 Eq.(34)/(35) —— f1 f2 修正的 McMillan 式。w 单位任意（这里 K）。"""
    L1 = 2.46 * (1.0 + 3.8 * mus)
    r = wb2 / wlog
    L2 = 1.82 * (1.0 + 6.3 * mus) * r
    f1 = (1.0 + (lam / L1) ** 1.5) ** (1.0 / 3.0)
    f2 = 1.0 + (r - 1.0) * lam ** 2 / (lam ** 2 + L2 ** 2)
    expo = -1.04 * (1.0 + lam) / (lam - mus * (1.0 + 0.62 * lam))
    return f1 * f2 * (wlog / 1.20) * np.exp(expo)


def frak_f(Tc, lam, wb2):
    return Tc / (C_INF * np.sqrt(lam) * wb2)


def regress(y, X):
    """y = A @ beta, A = [1, X...]. 返回 beta, R2, 残差散布(对数->相对%)。"""
    A = np.column_stack([np.ones_like(y)] + list(X))
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    ss_res = float(resid @ resid)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    R2 = 1.0 - ss_res / ss_tot
    # 对数残差的 rms -> 相对散布
    scat = float(np.sqrt(ss_res / len(y)))
    return beta, R2, np.expm1(scat) * 100.0, resid


x1 = np.log(lam)
x2 = np.log(wlog / wb2)

print("=" * 78)
print("A. 真实数据（实验 Tc）的回归  —— 复现 7.94% 那条")
f_obs = frak_f(Tc, lam, wb2)
b_obs, R2_obs, s_obs, res_obs = regress(np.log(f_obs), [x1, x2])
print(f"   log f = {b_obs[0]:+.4f} {b_obs[1]:+.4f}*log(lam) {b_obs[2]:+.4f}*log(wlog/wb2)")
print(f"   R2 = {R2_obs:.4f}   残差散布 = {s_obs:.2f}%   f in [{f_obs.min():.4f}, {f_obs.max():.4f}]")

print()
print("=" * 78)
print("B. 判决性检验：假设 AD 公式就是真理，用同样 29 组参数生成 f_AD，跑同一个回归")
Tc_AD = allen_dynes(lam, wlog, wb2, mus)
f_AD = frak_f(Tc_AD, lam, wb2)
b_AD, R2_AD, s_AD, res_AD = regress(np.log(f_AD), [x1, x2])
print(f"   log f = {b_AD[0]:+.4f} {b_AD[1]:+.4f}*log(lam) {b_AD[2]:+.4f}*log(wlog/wb2)")
print(f"   R2 = {R2_AD:.4f}   残差散布 = {s_AD:.2f}%   f_AD in [{f_AD.min():.4f}, {f_AD.max():.4f}]")
print()
print(f"   >>> lambda 指数   真实 {b_obs[1]:.4f}  vs  AD 世界 {b_AD[1]:.4f}   "
      f"差 {100*(b_obs[1]-b_AD[1])/b_AD[1]:+.2f}%")
print(f"   >>> shape  指数   真实 {b_obs[2]:.4f}  vs  AD 世界 {b_AD[2]:.4f}   "
      f"差 {100*(b_obs[2]-b_AD[2])/b_AD[2]:+.2f}%")
print(f"   >>> 常数项        真实 {b_obs[0]:.4f}  vs  AD 世界 {b_AD[0]:.4f}")

print()
print("=" * 78)
print("C. AD 公式自己对这 29 行的误差（伪影论的前提）")
rel = (Tc_AD - Tc) / Tc
print(f"   |Tc_AD/Tc_exp - 1|  均值 {np.abs(rel).mean()*100:.2f}%  "
      f"std {np.abs(rel).std()*100:.2f}%  最大 {np.abs(rel).max()*100:.2f}% ({name[int(np.argmax(np.abs(rel)))]})")
# AD 直接当 log f 的预测子（不拟合任何系数）时的 R2 与散布
d = np.log(f_obs) - np.log(f_AD)
ss_res_ad = float(d @ d)
ss_tot = float(((np.log(f_obs) - np.log(f_obs).mean()) ** 2).sum())
print(f"   AD 作为 log f 的**零参数**预测子：R2 = {1 - ss_res_ad/ss_tot:.4f}   "
      f"散布 = {np.expm1(np.sqrt(ss_res_ad/len(d)))*100:.2f}%")
print("   ^^^ 与 A 的 7.94% 对比：回归比 AD 公式本身更紧 ⟹ 回归不可能只是在重新提取 AD")

print()
print("=" * 78)
print("D. 逐点数值导数（agent 用的那个量）——它的散布有多大？")
h = 1e-3
dlogf_dloglam = np.empty(len(lam))
for i in range(len(lam)):
    lp, lm = lam[i] * (1 + h), lam[i] * (1 - h)
    fp = frak_f(allen_dynes(lp, wlog[i], wb2[i], mus[i]), lp, wb2[i])
    fm = frak_f(allen_dynes(lm, wlog[i], wb2[i], mus[i]), lm, wb2[i])
    dlogf_dloglam[i] = (np.log(fp) - np.log(fm)) / (np.log(lp) - np.log(lm))
print(f"   均值 {dlogf_dloglam.mean():.4f}  std {dlogf_dloglam.std():.4f}  "
      f"范围 [{dlogf_dloglam.min():.3f}, {dlogf_dloglam.max():.3f}]  "
      f"跨度 {dlogf_dloglam.max()/dlogf_dloglam.min():.2f}x")
print(f"   与回归指数 {b_AD[1]:.4f}（AD 世界的真回归系数）相比，均值 {dlogf_dloglam.mean():.4f} "
      f"偏 {100*(dlogf_dloglam.mean()-b_AD[1])/b_AD[1]:+.1f}%")
print("   ^^^ 逐点导数的均值 != 回归系数。用前者代替后者本身就是一步未经辩护的替换。")

print()
print("=" * 78)
print("E. 弱耦合闭式给的 lambda 指数（不依赖 AD 的 f1/f2）")
print("   log f = log(shape) - 0.5*log(lam) - 1.04(1+lam)/(lam - mu*(1+0.62lam)) - log(1.20*C_inf)")
print("   d/dlog(lam) 的逐点值：")
def mcm_only(lam_, wlog_, wb2_, mus_):
    expo = -1.04 * (1.0 + lam_) / (lam_ - mus_ * (1.0 + 0.62 * lam_))
    return (wlog_ / 1.20) * np.exp(expo)
dmc = np.empty(len(lam))
for i in range(len(lam)):
    lp, lm = lam[i] * (1 + h), lam[i] * (1 - h)
    fp = frak_f(mcm_only(lp, wlog[i], wb2[i], mus[i]), lp, wb2[i])
    fm = frak_f(mcm_only(lm, wlog[i], wb2[i], mus[i]), lm, wb2[i])
    dmc[i] = (np.log(fp) - np.log(fm)) / (np.log(lp) - np.log(lm))
print(f"   均值 {dmc.mean():.4f}  std {dmc.std():.4f}  范围 [{dmc.min():.3f}, {dmc.max():.3f}]")
f_mc = frak_f(mcm_only(lam, wlog, wb2, mus), lam, wb2)
b_mc, R2_mc, s_mc, _ = regress(np.log(f_mc), [x1, x2])
print(f"   纯 McMillan 世界的回归：lam 指数 {b_mc[1]:.4f}   shape 指数 {b_mc[2]:.4f}   R2 {R2_mc:.4f}")
print("   ^^^ 若纯 McMillan（无 f1/f2，即 AD 的**新**内容全删）也给 ~1，则指数不是 AD 的 f1/f2 带来的")

print()
print("=" * 78)
print("F. mu* 的作用：把 mu* 全部置 0（AD 拟合的那一列被删掉）后重跑")
f_mu0 = frak_f(allen_dynes(lam, wlog, wb2, np.zeros_like(mus)), lam, wb2)
b_m0, R2_m0, s_m0, _ = regress(np.log(f_mu0), [x1, x2])
print(f"   lam 指数 {b_m0[1]:.4f}   shape 指数 {b_m0[2]:.4f}   R2 {R2_m0:.4f}")

print()
print("=" * 78)
print("G. 参数区间的作用：把 lambda 网格换成与样本无关的均匀网格")
lam_g = np.exp(np.linspace(np.log(0.69), np.log(2.59), 29))
sh_g = np.full(29, float(np.mean(wlog / wb2)))
wlog_g = np.full(29, 100.0)
wb2_g = wlog_g / sh_g
f_g = frak_f(allen_dynes(lam_g, wlog_g, wb2_g, np.full(29, 0.10)), lam_g, wb2_g)
A = np.column_stack([np.ones(29), np.log(lam_g)])
bg, *_ = np.linalg.lstsq(A, np.log(f_g), rcond=None)
print(f"   合成网格（shape 固定、mu*=0.10、lam 对数均匀覆盖同一区间）：lam 指数 {bg[1]:.4f}")
print("   ^^^ 若它也 ~1，说明 ~1 是 AD 在**这个 lambda 窗口**上的性质，与样本的具体材料无关")
