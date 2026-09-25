#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tier-0 天花板检验：kTc <= C_inf * sqrt(eta/M)，**零声子计算**。

数据来源（2026-09-25 主会话亲取、双仪器核对）：
  D. A. Papaconstantopoulos, L. L. Boyer, B. M. Klein, A. R. Williams,
  V. L. Moruzzi, J. F. Janak, PRB 15, 4221 (1977).
  - TABLE III (p.4224): eta_sp, eta_pd, eta_df, eta_tot(I), eta_tot(II)  [eV/A^2]
  - TABLE IV  (p.4225): Theta_D, E_F, n(E_F), mu*, Tc(I), Tc(II), Tc_exp [K]
  取得方式：curl harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.15.4221/fulltext
  两仪器：pdftotext -enc UTF-8（OCR 文本层）+ pdftoppm 页图 Read。
  ⚠ 文本层单独用会错：OCR 把 eta_tot(I) 的列标签印成 "&tot(II)"，
    页图才定死 Nb 的 7.627 属于 (I) 列、(II) 列是 5.410。数值一律以**页图**为准。

严格性：
  - 天花板 kTc <= C_inf*sqrt(M2) 本身：analytic（C_inf = sqrt(g(2))/(2pi)，与 lambda、mu* 无关）
  - 本脚本的判决：numerical-evidence（eta 是作者的 DFT 结果，Tc_exp 是实验值；
    二者都不是我们能认证的量）。这是**外部效度检验**，不是证书。

方向：天花板是**上界**。eta 越大 => 天花板越高 => 越容易 PASS。
  eta_tot(I) >= eta_tot(II) 逐行成立 => 用 (I) 是**宽松**检验，用 (II) 是**严格**检验。
  ⚠ (II) 是作者把 d-f 贡献除以 2 的 ad hoc 订正，且该因子是为改善与 Tc 的一致性而选的
    => (II) 有**循环性风险**（同 hydride-csv-audit-filling-factor-circular 的形状）。
    故主判决用 (I)，(II) 只作敏感性。
"""
import math

C_INF = 0.1827262477          # = sqrt(g(2))/(2*pi)，analytic
HBAR  = 1.054571817e-34       # J s
EV    = 1.602176634e-19       # J
AMU   = 1.66053906660e-27     # kg
KB    = 1.380649e-23          # J/K
ANG   = 1.0e-10               # m

# eta[eV/A^2] / M[amu] -> omega^2 [s^-2] -> hbar*omega [eV]
#   1 eV/A^2 / 1 amu = (EV/ANG^2)/AMU  s^-2
UNIT_S2 = (EV / ANG**2) / AMU
HW_PER_SQRT = HBAR * math.sqrt(UNIT_S2) / EV      # eV per sqrt(eta/M)
K_PER_EV = EV / KB                                 # K per eV

# 元素: (eta_tot_I, eta_tot_II, Tc_exp[K], M[amu])
# eta 与 Tc_exp 逐字抄自页图 TABLE III / TABLE IV；M 为标准原子量（IUPAC）
DATA = [
    ("H",  6.478, 6.478, None,  1.008),   # Tc_exp = "..." (无)
    ("Li", 0.265, 0.265, 0.0,    6.94),
    ("Be", 0.560, 0.560, 0.026,  9.0122),
    ("Na", 0.017, 0.017, 0.0,   22.990),
    ("Mg", 0.123, 0.119, 0.0,   24.305),
    ("Al", 0.590, 0.558, 1.175, 26.982),
    ("K",  0.011, 0.011, 0.0,   39.098),
    ("Ca", 0.239, 0.178, 0.0,   40.078),
    ("Sc", 3.309, 2.530, 0.0,   44.956),
    ("Ti", 4.840, 3.200, 0.39,  47.867),
    ("V",  6.894, 4.370, 5.43,  50.942),
    ("Cr", 4.202, 2.400, 0.0,   51.996),
    ("Mn", 5.482, 2.980, 0.0,   54.938),
    ("Fe", 5.623, 2.920, 0.0,   55.845),
    ("Co", 3.900, 2.050, 0.0,   58.933),
    ("Ni", 1.762, 0.880, 0.0,   58.693),
    ("Cu", 0.456, 0.380, 0.0,   63.546),
    ("Zn", 0.544, 0.530, 0.375, 65.38),
    ("Ga", 0.595, 0.590, 1.080, 69.723),
    ("Rb", 0.016, 0.016, 0.0,   85.468),
    ("Sr", 0.000, 0.000, 0.0,   87.62),
    ("Y",  2.337, 1.700, 0.0,   88.906),
    ("Zr", 4.658, 3.520, 0.53,  91.224),
    ("Nb", 7.627, 5.410, 9.25,  92.906),
    ("Mo", 5.797, 3.440, 0.916, 95.95),
    ("Tc", 8.094, 4.590, 7.73,  98.0),
    ("Ru", 6.510, 3.520, 0.493,101.07),
    ("Rh", 4.722, 2.490, 0.0,  102.91),
    ("Pd", 1.992, 1.050, 0.0,  106.42),
    ("Ag", 0.288, 0.250, 0.0,  107.87),
    ("Cd", 0.359, 0.340, 0.52, 112.41),
    ("In", 0.355, 0.350, 3.40, 114.82),
]


def ceiling_K(eta, M):
    """kTc <= C_inf*sqrt(eta/M)  ->  Tc_ceiling in K"""
    if eta <= 0.0:
        return 0.0
    return C_INF * HW_PER_SQRT * math.sqrt(eta / M) * K_PER_EV


def main():
    print("Tier-0 天花板：kTc <= C_inf*sqrt(eta/M)，C_inf = %.10f" % C_INF)
    print("换算常数：sqrt(1 eV/A^2 / 1 amu) -> hbar*omega = %.8f eV" % HW_PER_SQRT)
    print("        => Tc_ceiling[K] = %.4f * sqrt(eta[eV/A^2]/M[amu])"
          % (C_INF * HW_PER_SQRT * K_PER_EV))
    print()
    hdr = ("%-4s %9s %9s %9s %10s %10s %8s %8s" %
           ("El", "eta(I)", "eta(II)", "Tc_exp", "ceil(I)", "ceil(II)",
            "r(I)", "r(II)"))
    print(hdr)
    print("-" * len(hdr))

    n_I = n_II = viol_I = viol_II = 0
    worst_I = worst_II = (0.0, "")
    rows = []
    for el, eI, eII, tc, M in DATA:
        cI, cII = ceiling_K(eI, M), ceiling_K(eII, M)
        if tc is None:
            print("%-4s %9.3f %9.3f %9s %10.2f %10.2f %8s %8s"
                  % (el, eI, eII, "n/a", cI, cII, "-", "-"))
            continue
        rI = tc / cI if cI > 0 else (float("inf") if tc > 0 else 0.0)
        rII = tc / cII if cII > 0 else (float("inf") if tc > 0 else 0.0)
        n_I += 1; n_II += 1
        if rI > 1.0: viol_I += 1
        if rII > 1.0: viol_II += 1
        if rI > worst_I[0]: worst_I = (rI, el)
        if rII > worst_II[0]: worst_II = (rII, el)
        rows.append((el, tc, cI, cII, rI, rII))
        print("%-4s %9.3f %9.3f %9.3f %10.2f %10.2f %8.4f %8.4f"
              % (el, eI, eII, tc, cI, cII, rI, rII))

    print()
    print("=== 判决（Tc_exp / 天花板；>1 = 违反）===")
    print("  eta_tot(I)  [宽松]: %d/%d 违反；最劣 %s r=%.4f"
          % (viol_I, n_I, worst_I[1], worst_I[0]))
    print("  eta_tot(II) [严格]: %d/%d 违反；最劣 %s r=%.4f"
          % (viol_II, n_II, worst_II[1], worst_II[0]))

    # 方向自查：(II) <= (I) 必须逐行成立，否则"严格/宽松"的说法就反了
    bad = [el for el, eI, eII, _, _ in DATA if eII > eI + 1e-12]
    print()
    print("[方向门] eta_tot(II) <= eta_tot(I) 逐行：%s"
          % ("PASS（%d/%d）" % (len(DATA) - len(bad), len(DATA)) if not bad
             else "FAIL " + str(bad)))

    # 只看真正超导的那些（Tc_exp > 0）——Tc=0 的行对上界是平凡满足
    sc = [r for r in rows if r[1] > 0.0]
    print("[非平凡子集] Tc_exp>0 的有 %d 个：%s"
          % (len(sc), ", ".join(r[0] for r in sc)))
    if sc:
        mI = max(sc, key=lambda r: r[4]); mII = max(sc, key=lambda r: r[5])
        print("  最紧 (I) : %s  Tc=%.3f  ceil=%.2f  r=%.4f" % (mI[0], mI[1], mI[2], mI[4]))
        print("  最紧 (II): %s  Tc=%.3f  ceil=%.2f  r=%.4f" % (mII[0], mII[1], mII[3], mII[5]))
        print("  中位 r(I)=%.4f  r(II)=%.4f"
              % (sorted(r[4] for r in sc)[len(sc)//2],
                 sorted(r[5] for r in sc)[len(sc)//2]))


def unit_crosscheck():
    """
    单位换算是本脚本唯一的承重技术量 => 核它。

    ⚠ 第一版我写的是「Debye 谱下 sqrt(<w^2>)/Theta_D 应 ~0.775」，实测**五个元素全给 0.707**。
      五个独立元素给到三位相同的数**不是一致，是恒等式**（自由度 0，同
      hydride-csv-audit-filling-factor-circular / MgB2 mu* 的形状）。
      追下去发现真因：Table III 的 lambda 列本身就是 (eta/M)/(Theta_D^2/2)
      （McMillan: lambda = eta/(M<w^2>)，作者取 <w^2> = Theta_D^2/2）
      => sqrt(M2/lambda)/Theta_D = 1/sqrt(2) = 0.7071 **恒成立**，与元素无关。

    所以正确的判词有两条，都比原来那条强：
      (a) 单位：我的换算常数**复现了作者自己的算术**，12/12 落在印出精度的舍入容差内。
          若差一个 2pi 或一个量级，lambda_pred 会差同样倍数，12 个全出容差。
      (b) 自由度：Table III 的 lambda 列**没有独立信息**，Table IV 的 Tc(I)/Tc(II) 亦然
          （McMillan 公式从 lambda, mu*, Theta_D 算出）。
          => 真正独立的输入只有 eta、Theta_D、mu*、Tc_exp 四个。
          本脚本用的是 eta、M、**Tc_exp**（实验值）=> **无循环性**。
    """
    print()
    print("=== 单位换算核验（判词：能否复现作者自己的 lambda 列）===")
    probe = [("V", 380, 1.055), ("Nb", 275, 1.222), ("Mo", 450, 0.336),
             ("Zr", 291, 0.679), ("Ti", 420, 0.645), ("Tc", 411, 0.550),
             ("Ru", 600, 0.201), ("Pd", 274, 0.281), ("Al", 428, 0.135),
             ("In", 108, 0.298), ("Sc", 360, 0.639), ("Fe", 467, 0.520)]
    d = {r[0]: (r[1], r[4]) for r in DATA}
    nok = 0
    for el, thD, lam in probe:
        eI, M = d[el]
        M2_K2 = (HW_PER_SQRT * math.sqrt(eI / M) * K_PER_EV) ** 2
        pred = M2_K2 / (thD ** 2 / 2.0)
        diff = abs(pred - lam)
        # 三个印出量各自的半 ulp 传播：eta(3dp)、Theta_D(整数)、lambda(3dp)
        tol = lam * (0.0005 / eI + 2 * 0.5 / thD) + 0.0005
        ok = diff <= tol
        nok += ok
        print("  %-3s lam_tab=%.3f  lam_pred=%.4f  |d|=%.2e  tol=%.2e  %s"
              % (el, lam, pred, diff, tol, "OK" if ok else "X"))
    print("  => %d/%d 落在纯舍入容差内。" % (nok, len(probe)))
    print("     (a) 单位常数与作者一致（差 2pi 则 12/12 全出容差）；")
    print("     (b) lambda 列 = (eta/M)/(Theta_D^2/2)，**零独立信息** =>")
    print("         作者的 Tc(I)/Tc(II) 也是导出量；本检验用 Tc_exp，故**无循环性**。")


if __name__ == "__main__":
    main()
    unit_crosscheck()
