# -*- coding: utf-8 -*-
r"""Brorson PRL 64,2172 (1990) 十行 + Mansart PRB 82,024513 (2010) 一行：超快实测 lam<w^2> 进天花板。

一手来源（主会话双仪器：pdftotext -enc UTF-8 + 页图裁切逐字读）：
  Brorson : harvest.aps.org/.../PhysRevLett.64.2172/fulltext （JBIG2 扫描件，文本层 OCR 退化，
            数值一律以页图 _brorson/crop_p2_table-2.png / crop_p2_tail-2.png 为准）
            Eq.(1)  T_e dT_e/dt = -(3hbar/(pi k_B)) lam<w^2> (T_e - T_l)   （= Allen PRL 59,1460 Eq.(14) 首项）
            <w^2> "defined by McMillan"  ⟹ 与 AD75 / Allen Eq.(11) n=2 同一个量（G3 定义同构门）
  Mansart : harvest.aps.org/.../PhysRevB.82.024513/fulltext （矢量 PDF，文本层 + 页图一致）
            3TM:  2 C_e dT_e/dt = ... - g (T_e - T_1)   ；  g = (6 hbar gamma/(pi k_B)) lam<w^2>  [Ref.28]
            6.8 < g < 8.0 "mJ/K/s/m^3"（单位印错，缺 10 的幂）；gamma = 21 mJ/mol/K^2；lam<w^2> ≅ 64 meV^2；
            Tc = 24 K（x=0.08, "same batches"）

本脚本做五件事（每件都带投毒）：
  (1) Brorson Table I 内部一致性：lam_exp == lam<w^2>/<w^2>，按印出精度传播容差   —— 转录门
  (2) PASS 行：Tc <= K sqrt(X)。X_fit <= X_true（Allen Eq.14 修正项为负；NbN/V3Ga 作者自称下界）
      ⟹ X_fit 是下界 ⟹ PASS 消费点方向安全。报 r 与「翻判决所需 Tc」。
  (3) RMTA 跨方法（eta 级，实验锚定）：kappa := X_RMTA / X_exp，X_RMTA = eta/M 换成 meV^2。
      方向：X_exp 是下界 ⟹ kappa_fit 是 kappa_true 的**上**界（kappa 的高估消费点——r 放大——正需上界）。
  (4) Nb3Ge 头条包封按新 kappa 重算（外推是指示性的，照实标）。
  (5) Mansart：g->X 的换算因子 6 vs 3 判别；天花板判别器预算；印出 ODE 的 2C_e 歧义换算成 r。
"""
import math
import sys

MEV_K = 1.602176634e-22 / 1.380649e-23
C_INF = 0.1827262477
K_MEV = C_INF * MEV_K                      # Tc[K] <= K_MEV * sqrt(X[meV^2])
HBAR = 1.054571817e-34
KB = 1.380649e-23
EV = 1.602176634e-19
AMU = 1.66053906660e-27
NA = 6.02214076e23
# eta[eV/A^2]/M[amu] -> omega^2 -> (hbar omega)^2 in meV^2
ETA_TO_MEV2 = HBAR ** 2 * (EV / 1e-20 / AMU) / (1e-3 * EV) ** 2

FAIL = []


def chk(name, ok, detail=''):
    print('  [%s] %s  %s' % ('OK' if ok else '**FAIL**', name, detail))
    if not ok:
        FAIL.append(name)


def ceil_K(x):
    return K_MEV * math.sqrt(x)


# ---------------------------------------------------------------- Brorson Table I（页图逐字）
# (name, T_e(0) K, X, dX, <w^2>, lam_exp, dlam_exp, lam_lit or None)
BRORSON = [
    ("Cu",   590,  29,  4, 377, 0.08, 0.01, 0.10),
    ("Au",   650,  23,  4, 178, 0.13, 0.02, 0.15),
    ("Cr",   716, 128, 15, 987, 0.13, 0.02, None),
    ("W",   1200, 112, 15, 425, 0.26, 0.04, 0.26),
    ("V",    700, 280, 20, 352, 0.80, 0.06, 0.82),
    ("Nb",   790, 320, 30, 275, 1.16, 0.11, 1.04),
    ("Ti",   820, 350, 30, 601, 0.58, 0.05, 0.54),
    ("Pb",   570,  45,  5,  31, 1.45, 0.16, 1.55),
    ("NbN", 1070, 640, 40, 673, 0.95, 0.06, 1.46),
    ("V3Ga", 1110, 370, 60, 448, 0.83, 0.13, 1.12),
]
AUTHOR_LOWER_BOUND = {"NbN", "V3Ga"}       # p.2174 页图："should cautiously be regarded as lower bounds"
CU_OVERLAYER = {"Nb", "V", "Pb"}           # p.2173 页图："(viz., Nb, V, and Pb) ... overlayer (<= 40 A) of Cu"


def table_consistent(rows):
    bad = []
    for n, _, X, _, w2, lam, _, _ in rows:
        lo = (X - 0.5) / (w2 + 0.5) - 0.005
        hi = (X + 0.5) / (w2 - 0.5) + 0.005
        if not (lo <= lam <= hi):
            bad.append(n)
    return bad


print('== (1) Brorson Table I 内部一致性（转录门）')
bad = table_consistent(BRORSON)
chk('lam_exp = X/<w^2> 在印出精度内 10/10', not bad, 'bad=%s' % bad)
pois = [list(r) for r in BRORSON]
pois[7][2] = 54                             # 投毒：把 Pb 的 45 误读成 54
chk('投毒：Pb X 45->54 必须被抓', table_consistent([tuple(r) for r in pois]) == ['Pb'])

# ---------------------------------------------------------------- (2) PASS 行
# Tc：同样品 Tc 本文未印（Nb/V 带 Cu 覆层的样品 Tc 被邻近效应压到 <1.5 K / <4 K；同炉无覆层样品
# "showed the usual values" 无数字）。取**偏高**的体材料/文献最高值 ⟹ PASS 更难 ⟹ 方向保守。
TC_CONS = {"W": 0.0154, "V": 5.43, "Nb": 9.25, "Ti": 0.39, "Pb": 7.20, "NbN": 17.3, "V3Ga": 16.8}
print('\n== (2) PASS 行：Tc_cons <= K sqrt(X_fit)，X_fit 为下界 ⟹ 方向安全')
print('   %-5s %6s %9s %8s %8s %12s' % ('', 'X', 'ceil(K)', 'Tc_cons', 'r', 'X_fail/X'))
rs = {}
for n, _, X, dX, *_ in BRORSON:
    if n not in TC_CONS:
        continue
    c = ceil_K(X)
    r = TC_CONS[n] / c
    rs[n] = r
    print('   %-5s %6d %9.2f %8.3f %8.4f %12.2f' % (n, X, c, TC_CONS[n], r, 1 / r ** 2))
    # 用 1σ 下端 X 再判一次（更难）
    chk('%s PASS（X 取 1σ 下端 %d）' % (n, X - dX), TC_CONS[n] < ceil_K(X - dX))
worst = max(rs, key=rs.get)
chk('最紧行是 Pb', worst == 'Pb', 'r=%.4f' % rs[worst])
chk('投毒：Pb 取 Tc=20 K 必须 FAIL', not (20.0 < ceil_K(45 - 5)))

# ---------------------------------------------------------------- (3) RMTA vs 实验，eta 级
print('\n== (3) RMTA(Pap1977 / Klein1979) vs Brorson：kappa = X_RMTA / X_exp')
chk('单位常数 = (137.0959/2.120450)^2（Tier-0 脚本常数）',
    abs(ETA_TO_MEV2 - (137.0959 / 2.120450) ** 2) / ETA_TO_MEV2 < 2e-6, '%.3f meV^2' % ETA_TO_MEV2)
# Pap1977 TABLE III（tier0 脚本转录）: eta_tot(I), eta_tot(II) [eV/A^2], M [amu]
PAP = {"Ti": (4.840, 3.200, 47.867), "V": (6.894, 4.370, 50.942), "Cr": (4.202, 2.400, 51.996),
       "Cu": (0.456, 0.380, 63.546), "Nb": (7.627, 5.410, 92.906)}
# Klein1979 V3Ga: eta_A（子晶格合计，见 klein1979 重数约定）, eta_B
X_RMTA = {n: e1 / m * ETA_TO_MEV2 for n, (e1, e2, m) in PAP.items()}
X_RMTA_II = {n: e2 / m * ETA_TO_MEV2 for n, (e1, e2, m) in PAP.items()}
X_RMTA["V3Ga"] = (9.26 / 50.942 + 0.05 / 69.723) * ETA_TO_MEV2
BX = {n: (X, dX) for n, _, X, dX, *_ in BRORSON}
kap = {}
print('   %-5s %9s %9s %8s %8s %8s   %s' % ('', 'X_RMTA', 'X_exp', 'kappa', 'k(-1σ)', 'k(II)', '备注'))
for n in ["Cu", "Ti", "V", "Cr", "Nb", "V3Ga"]:
    X, dX = BX[n]
    k = X_RMTA[n] / X
    kap[n] = k
    k2 = X_RMTA_II[n] / X if n in X_RMTA_II else float('nan')
    note = {'Cr': 'AFM@300K 对象失配', 'V3Ga': '作者自称下界 ⟹ kappa 为上界',
            'V': 'Cu 覆层', 'Nb': 'Cu 覆层'}.get(n, '')
    print('   %-5s %9.1f %9d %8.3f %8.3f %8.3f   %s' % (n, X_RMTA[n], X, k, X_RMTA[n] / (X - dX), k2, note))
chk('V 上 RMTA 高估（kappa>1.5，1σ 内仍 >1.5）', kap['V'] > 1.5 and X_RMTA['V'] / (280 + 20) > 1.5,
    'kappa(V)=%.3f' % kap['V'])
chk('Nb 上 RMTA 与实验 1σ 相符', X_RMTA['Nb'] / (320 + 30) <= 1.0 <= X_RMTA['Nb'] / (320 - 30) or
    abs(kap['Nb'] - 1) < 0.1, 'kappa(Nb)=%.3f' % kap['Nb'])
# 旁证：V 的 X_exp 与作者给的 lam_lit<w^2>_lit（独立于本实验）一致
chk('V: X_exp 与 lam_lit*<w^2> 相差 <5%', abs(280 - 0.82 * 352) / 280 < 0.05, '%.1f' % (0.82 * 352))

# ---------------------------------------------------------------- (4) 头条包封
print('\n== (4) Nb3Ge 头条 r=0.626 在新 kappa 下（元素/化合物 -> Nb3Ge 的外推是指示性的）')
R0 = 0.626
k_flip = (1 / R0) ** 2
print('   翻判决所需 kappa = %.4f' % k_flip)
for lab, k in [('旧包封 Sc(Bose2008)', 1.3969), ('Nb 元素', kap['Nb']),
               ('Nb 元素 1σ', X_RMTA['Nb'] / 290), ('V 元素', kap['V']),
               ('V 元素 1σ', X_RMTA['V'] / 260), ('V3Ga 化合物(上界)', kap['V3Ga']),
               ('V3Ga 化合物 1σ', X_RMTA['V3Ga'] / 310), ('V3Ga 化合物 2σ', X_RMTA['V3Ga'] / 250),
               ('Cr（对象失配）', kap['Cr'])]:
    print('   %-22s kappa=%.4f  r<=%.4f  %s' % (lab, k, R0 * math.sqrt(k), 'PASS' if k < k_flip else '**翻**'))
chk('超导体元素/化合物的 kappa 全部 < 翻判决阈值', max(kap['V'], kap['V3Ga'], kap['Nb'], kap['Ti']) < k_flip)
# 同源（Klein1979 同表）A15 的 kappa 在 1σ 下仍过、2σ 下翻 —— 如实报，余量只剩 ~2%
chk('V3Ga 1σ：Nb3Ge 包封仍 PASS', X_RMTA['V3Ga'] / 310 < k_flip, 'r<=%.4f' % (R0 * math.sqrt(X_RMTA['V3Ga'] / 310)))
chk('V3Ga 2σ：Nb3Ge 包封翻（如实报）', X_RMTA['V3Ga'] / 250 > k_flip, 'r<=%.4f' % (R0 * math.sqrt(X_RMTA['V3Ga'] / 250)))
chk('Cr 的 kappa 超过阈值（如实报，不许吞掉）', kap['Cr'] > k_flip, '%.3f' % kap['Cr'])
# V3Ga：实验分母直接给 r 上界（X_exp 下界 ⟹ r 高估 ⟹ 安全）
r_v3ga = 16.5 / ceil_K(370)
chk('V3Ga 实验分母 r（Tc=16.5，Klein 表同值）== Tier-0 RMTA 版 0.282 × sqrt(kappa)（恒等）',abs(r_v3ga - 0.282 * math.sqrt(kap['V3Ga'])) < 2e-3,
    'r=%.4f（恒等式，非独立证据）' % r_v3ga)

# ---------------------------------------------------------------- (5) Mansart
print('\n== (5) Mansart Ba(Fe0.92Co0.08)2As2')
a, c = 3.96e-10, 12.98e-10                 # 文献常用晶格常数（未一手取；判别只要求误差 << 2x）
Vm = a * a * c / 2 * NA                    # Z=2
gam = 21e-3 / Vm                           # J/m^3/K^2


def X_from_g(g, f):
    return g * math.pi * KB / (f * HBAR * gam) * HBAR ** 2 / (1e-3 * EV) ** 2


def bracket_ok(f, X=64.0):
    """有没有某个 10 的幂 k，使印出尾数 6.8/8.0 在因子 f 下把 X 括住。"""
    for k in range(14, 21):
        lo, hi = X_from_g(6.8 * 10 ** k, f), X_from_g(8.0 * 10 ** k, f)
        if lo <= X <= hi:
            return k, lo, hi
    return None


b6, b3 = bracket_ok(6), bracket_ok(3)
chk('因子 6 下存在 10^k 使 [6.8,8.0] 括住 64', b6 is not None, str(b6))
chk('因子 3 下不存在任何 10^k 括住 64（判别成立）', b3 is None)
chk('投毒：把 64 换成 128，因子 3 必须变成能括住', bracket_ok(3, 128.0) is not None)

Tc_m = 24.0
X_star = (Tc_m / K_MEV) ** 2
print('   判别器阈值 X* = %.2f meV^2；印出 64 ⟹ r=%.4f；所需低估倍数 %.4f' % (X_star, Tc_m / ceil_K(64), X_star / 64))
for Xi in (59.2, 60.9, 66.1, 69.6):
    print('     Table I 行 X=%.1f  r=%.4f' % (Xi, Tc_m / ceil_K(Xi)))
# 印出 ODE 是 2C_e 且印出 g-关系是 6：二者自洽地给出 Allen 弛豫率 3hbar X/(pi kB Te)。
# 若代码实际用的是 C_e（"2" 为排印），则 Allen 约定下 X_A = 2*64：
r_alt = Tc_m / ceil_K(128.0)
print('   若 ODE 的 "2" 是排印（代码用 C_e）⟹ X_A=128 ⟹ r=%.4f' % r_alt)
chk('歧义因子 2 恰好跨过阈值附近（r_alt 与 1 相差 <1%）', abs(r_alt - 1) < 0.01, '%.4f' % r_alt)
chk('⟹ 判别器不得无条件入账：所需低估 %.3fx <= 已点名的未量化机制 2x' % (X_star / 64), X_star / 64 <= 2.0 + 0.01)

print()
if FAIL:
    print('FAILED:', FAIL)
    sys.exit(1)
print('ALL CHECKS PASSED')
