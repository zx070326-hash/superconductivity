#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
跨方法核对：Savrasov & Savrasov PRB 54, 16487 (1996)（线性响应 ab initio EPI）
          vs Papaconstantopoulos et al. PRB 15, 4221 (1977)（APW + 刚性 muffin-tin
          近似 / Gaspari-Gyorffy，本仓库 Tier-0 判决的 eta 数据源）。

承重结论（写在 theory/notes/literature/savrasov_vs_pap1977_crossmethod_eta_2026-09-25.md）：

**Savrasov 1996 全文只给出 4 张表**（TABLE I-IV，页图 + pdftotext 双仪器核对，
_savr_text.txt 逐行 grep "TABLE" 只命中 4 处）：
  I   计算 vs 实验声子频率（THz，高对称点）
  II  lambda_calc / lambda_tun / lambda_s-h / N(eps_F) / gamma
  III mu*、omega_log、Delta_0(calc/exp)
  IV  lambda_tr（输运通道）、Theta_tr=sqrt(<omega^2>_tr)、omega_p

**没有 eta（Hopfield 参数）表，也没有超导通道的 <omega^2>（二阶矩）表。**
文中式(2)出现的符号 "eta_{q nu}(R_mu)" 是声子简正模式位移矢量分量，
与 McMillan-Hopfield 的 eta 是**同一个希腊字母、不同的物理量**
（CLAUDE.md：「一个符号在链条里只许有一个意思」——这里恰是反例，需显式分离）。
Theta_tr 是**输运**通道（alpha_tr^2 F，用 lambda_tr 不是 lambda_calc）的二阶矩平方根，
作者自陈 lambda_tr 与超导 lambda_calc 相差可达 20%（正文 p.16497）——
不能拿 Theta_tr 当超导通道 <omega^2> 的替身，否则引入一个不受控的假设。

⟹ **无法从 Savrasov 1996 独立反推 eta_Savrasov**（既没有直接的 eta，
也没有可信的超导通道 <omega^2> 来配 lambda_calc 反推）。

**可达成的、诚实的跨方法检验**改为：直接比较两篇论文各自独立算出的
**lambda（无量纲、两边都已发表，不需要任何单位换算）**：
  lambda_calc(Savrasov，全 BZ 线性响应，1996)
    vs
  lambda_tot(I)/(II)(Pap1977，APW+刚性 muffin-tin 微扰，eta 与 Theta_D 合成，1977)

这弱于「eta 对 eta」，但仍然是**独立方法的交叉验证**，因为
lambda_tot(Pap1977) = eta_tot(Pap1977) / (M * Theta_D^2 / 2)（已证 12/12 一致，
tier0 笔记 §4），比 lambda 就是比 eta 经同一个 (M, Theta_D) 换算后的量；
唯一没被交叉验证的是 Theta_D^2/2 这个 <omega^2> 近似本身。

严格性：numerical-evidence（两边都是他人发表的 DFT/半经验数，非我们能认证的量）。
"""
import math

# ---- Savrasov 1996 Table II, lambda_calc（页图 savr-08.png 核对；无单位，纯数）----
SAVR_LAMBDA = {
    "Al": 0.44, "Pb": 1.68, "V": 1.19, "Nb": 1.26,
    "Ta": 0.86, "Mo": 0.42, "Cu": 0.14, "Pd": 0.35,
}

# ---- Papaconstantopoulos 1977 Table III, lambda_tot(I)/(II)（页图 pap77_p-4.png）----
PAP_LAMBDA_I = {
    "H": 2.102, "Li": 0.363, "Be": 0.034, "Na": 0.033, "Mg": 0.036,
    "Al": 0.135, "K": 0.040, "Ca": 0.128, "Sc": 0.639, "Ti": 0.645,
    "V": 1.055, "Cr": 0.229, "Mn": 0.668, "Fe": 0.520, "Co": 0.376,
    "Ni": 0.167, "Cu": 0.069, "Zn": 0.088, "Ga": 0.094, "Rb": 0.068,
    "Sr": 0.000, "Y": 0.378, "Zr": 0.679, "Nb": 1.222, "Mo": 0.336,
    "Tc": 0.550, "Ru": 0.201, "Rh": 0.224, "Pd": 0.281, "Ag": 0.059,
    "Cd": 0.082, "In": 0.298,
}
PAP_LAMBDA_II = {
    "H": 2.102, "Li": 0.363, "Be": 0.034, "Na": 0.033, "Mg": 0.035,
    "Al": 0.127, "K": 0.039, "Ca": 0.096, "Sc": 0.489, "Ti": 0.426,
    "V": 0.669, "Cr": 0.131, "Mn": 0.364, "Fe": 0.270, "Co": 0.197,
    "Ni": 0.084, "Cu": 0.058, "Zn": 0.086, "Ga": 0.094, "Rb": 0.065,
    "Sr": 0.000, "Y": 0.276, "Zr": 0.514, "Nb": 0.867, "Mo": 0.199,
    "Tc": 0.312, "Ru": 0.109, "Rh": 0.118, "Pd": 0.148, "Ag": 0.051,
    "Cd": 0.077, "In": 0.298,
}
# ---- Papaconstantopoulos 1977 Table III, eta_tot(I)/(II) [eV/A^2]（页图核对）----
PAP_ETA_I = {"Al": 0.590, "Cu": 0.456, "V": 6.894, "Nb": 7.627, "Mo": 5.797, "Pd": 1.992}
PAP_ETA_II = {"Al": 0.558, "Cu": 0.380, "V": 4.370, "Nb": 5.410, "Mo": 3.440, "Pd": 1.050}
# ---- Papaconstantopoulos 1977 Table IV, Theta_D [K]（页图 pap77_p-5.png）----
PAP_THETA_D = {"Al": 428, "Cu": 343, "V": 380, "Nb": 275, "Mo": 450, "Pd": 274}
# 标准原子量（IUPAC，amu）
M_AMU = {"Al": 26.982, "Cu": 63.546, "V": 50.942, "Nb": 92.906, "Mo": 95.95, "Pd": 106.42}

EV = 1.602176634e-19
AMU = 1.66053906660e-27
ANG = 1.0e-10
KB = 1.380649e-23  # not needed here, kept for consistency with tier0 script


def m2_from_eta(eta_evA2, M_amu, theta_D_K):
    """独立复核用：由 eta 与 Theta_D 反推 lambda_tab，核 12/12 一致性（这里对 6 个元素核）"""
    # <omega^2> = (Theta_D * kB/hbar)^2 / 2   [s^-2]   (McMillan Debye 近似, 作者自己的约定)
    HBAR = 1.054571817e-34
    omega_D = KB * theta_D_K / HBAR  # s^-1
    omega2 = omega_D ** 2 / 2.0
    eta_SI = eta_evA2 * EV / ANG ** 2  # J / m^2 ... (Hopfield eta units: energy/length^2 per atom)
    M_SI = M_amu * AMU
    lam = eta_SI / (M_SI * omega2)
    return lam


def poison_test():
    """投毒对照：故意漏 /2 或漏平方，检验单位换算对该扰动敏感"""
    print("=== 投毒对照（复现 tier0 笔记 §4 的 lambda_tab 反推，扩到 6 元素）===")
    ok = 0
    for el in ["Al", "V", "Nb", "Mo", "Cu", "Pd"]:
        lam_correct = m2_from_eta(PAP_ETA_I[el], M_AMU[el], PAP_THETA_D[el])
        lam_tab = PAP_LAMBDA_I[el]
        rel = abs(lam_correct - lam_tab) / lam_tab
        flag = "OK" if rel < 0.02 else "MISMATCH"
        if flag == "OK":
            ok += 1
        print(f"  {el:3s}  lambda_tab={lam_tab:.3f}  lambda_reconstructed={lam_correct:.3f}"
              f"  rel_diff={rel*100:.2f}%  [{flag}]")
    print(f"  -> {ok}/6 落在 2% 容差内（正对照，应几乎全过）")

    print()
    print("  投毒：故意漏掉 '/2'（即用 Theta_D 平方本身当 <omega^2>，多乘 2）")
    n_fail = 0
    for el in ["Al", "V", "Nb", "Mo", "Cu", "Pd"]:
        HBAR = 1.054571817e-34
        omega_D = KB * PAP_THETA_D[el] / HBAR
        omega2_poisoned = omega_D ** 2  # 漏掉 /2
        eta_SI = PAP_ETA_I[el] * EV / ANG ** 2
        M_SI = M_AMU[el] * AMU
        lam_poisoned = eta_SI / (M_SI * omega2_poisoned)
        rel = abs(lam_poisoned - PAP_LAMBDA_I[el]) / PAP_LAMBDA_I[el]
        if rel > 0.3:
            n_fail += 1
        print(f"  {el:3s}  lambda_tab={PAP_LAMBDA_I[el]:.3f}  lambda_poisoned={lam_poisoned:.3f}"
              f"  rel_diff={rel*100:.1f}%")
    print(f"  -> {n_fail}/6 投毒后偏差 >30%（应几乎全部报警，若否则检验对单位常数不敏感）")


def crossmethod_lambda_table():
    print()
    print("=== 跨方法 lambda 比较：Savrasov(1996, 线性响应) vs Pap1977(APW+RMTA) ===")
    print(f"{'El':4s} {'lam_Savr':>10s} {'lam_PapI':>10s} {'lam_PapII':>10s} "
          f"{'ratio_I':>9s} {'ratio_II':>9s}")
    overlap = ["Al", "V", "Nb", "Mo", "Cu", "Pd"]  # Z<=49, present in both papers
    for el in overlap:
        ls = SAVR_LAMBDA[el]
        lpI = PAP_LAMBDA_I[el]
        lpII = PAP_LAMBDA_II[el]
        rI = ls / lpI
        rII = ls / lpII
        print(f"{el:4s} {ls:10.3f} {lpI:10.3f} {lpII:10.3f} {rI:9.3f} {rII:9.3f}")
    print()
    print("Ta, Pb 在 Savrasov 表中但 Z>49，Pap1977 (Z<=49) 未覆盖 -> 未取得，不比较。")
    print("H, Li, ... 等 Pap1977 覆盖但 Savrasov 1996 未算 -> 未取得，不比较。")


if __name__ == "__main__":
    poison_test()
    crossmethod_lambda_table()
