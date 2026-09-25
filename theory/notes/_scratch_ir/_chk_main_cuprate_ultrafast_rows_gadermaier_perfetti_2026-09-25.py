"""
cuprate 超快 lambda<w^2> 行：一手来源 Gadermaier PRL 105,257001 (2010) 与 Perfetti PRL 99,197001 (2007)。
主会话 2026-09-25；页图裁切在 library/pdf/eta_hunt_2026-09-25/{gadermaier,perfetti}_*.png。

消费点：天花板的「声子介导判别器」（ceiling_as_phonon_mediation_discriminator_2026-09-25.md）。
判别器需 X := lambda<w^2> 的**上**界 ⟹ 一律条件式；本脚本把条件换算成「所需低估倍数」并与
**本文自身**已点名的模型因子比（不拍阈值）。

  (1) G1 订正：二手表 YBCO 用 Tc=90/92 K；一手样品是 YBa2Cu3O6.5，Tc=60 K ⟹ r 被抬 1.5x
  (2) Gadermaier：用 NEM Eq.(2) lambda<w^2> = (2pi/3) kB T_l/(hbar tau)；Allen/TTM Eq.(1) = (pi/3) kB T_e/(hbar tau)
      同一 tau 下 X_TTM/X_NEM = T_e/(2 T_l)；T_e 取文中「typical 400–800 K」，T_l = 室温
  (3) Perfetti：三浴 Einstein 模型，lambda Omega0^2 = 360±30（Table I 300–380）；(1-f) 模的耦合被置零
      ⟹ X_fit 只是强耦合子集 ⟹ 低估全谱；判别器失效条件 = 被忽略模的 X 超过子集的 (need-1) 倍
"""
import math

K = 2.120450   # Tc[K] <= K sqrt(X[meV^2])，C_inf * 11.6045181
FAIL = []


def chk(name, ok, info=''):
    print('  [%s] %s  %s' % ('OK' if ok else 'FAIL', name, info))
    if not ok:
        FAIL.append(name)


def r_of(tc, X):
    return tc / (K * math.sqrt(X))


print('== (1) G1：二手 Tc vs 一手同样品 Tc（YBCO）')
r2 = r_of(90, 400); r1 = r_of(60, 400)
print('   二手 (Tc=90, X=400) r=%.4f need=%.3fx ; 一手 (YBa2Cu3O6.5 Tc=60) r=%.4f need=%.4fx' % (r2, r2 ** 2, r1, r1 ** 2))
chk('二手 r 与讨论笔记表 2.122 一致（正对照：我复现的是同一个数）', abs(r2 - 2.122) < 5e-4, '%.4f' % r2)
chk('一手 Tc 把 r 降为 2/3（纯 Tc 比，X 不变）', abs(r1 / r2 - 60 / 90) < 1e-12)
chk('LSCO (Tc=38, X=800) r<1 ⟹ 无判别力', r_of(38, 800) < 1, '%.4f' % r_of(38, 800))

print('\n== (2) Gadermaier YBa2Cu3O6.5：NEM vs TTM(Allen) 换算')
X0, dX = 400.0, 100.0
Tl = 300.0
conv = [Te / (2 * Tl) for Te in (400.0, 800.0)]
print('   X_TTM/X_NEM = T_e/(2T_l) ∈ [%.4f, %.4f]' % tuple(conv))
rows = []
for lab, X in [('中心', X0), ('+1σ', X0 + dX), ('+1σ ×TTM上端', (X0 + dX) * conv[1]),
               ('+2σ ×TTM上端', (X0 + 2 * dX) * conv[1])]:
    rows.append((lab, X, r_of(60, X)))
    print('   %-14s X=%7.2f  r=%.4f' % (lab, X, r_of(60, X)))
chk('中心值所需低估 = 2.0016x（与 Mansart 行数值巧合：60/sqrt(400)=24/8=3）', abs(r1 ** 2 - 2.00157) < 1e-4, '%.5f' % r1 ** 2)
chk('+1σ 叠 TTM 上端仍 r>1', rows[2][2] > 1, '%.4f' % rows[2][2])
chk('+2σ 叠 TTM 上端落在阈值上（|r-1|<1%）⟹ 条件式', abs(rows[3][2] - 1) < 0.01, '%.4f' % rows[3][2])
# 投毒：若把 NEM 系数误抄成 agent 报的 2/3（丢了 pi），换算比会差 pi 倍 —— 本门必须能区分
conv_bad = [Te * (math.pi / 3) / (Tl * 2 / 3) for Te in (400.0, 800.0)]
chk('投毒：丢 pi 的 NEM 系数给出的换算比与正确值差 pi 倍（可区分）',
    abs(conv_bad[1] / conv[1] - math.pi) < 1e-12, '%.4f vs %.4f' % (conv_bad[1], conv[1]))

print('\n== (3) Perfetti Bi2212 (Tc=91 K 同样品)')
for X in (300.0, 330.0, 360.0, 380.0, 390.0):
    r = r_of(91, X)
    print('   X=%5.0f  r=%.4f  need=%.3fx  ⟹ 被忽略 (1-f) 模须 >= %.3f x 子集' % (X, r, r * r, r * r - 1))
rmin = r_of(91, 390.0)
chk('Table I 上端 380 与 abstract +1σ 390 下 r 仍 >2', r_of(91, 380) > 2 and rmin > 2, '%.4f / %.4f' % (r_of(91, 380), rmin))
chk('失效条件：被忽略模的 X >= 3.7x 子集（在 +1σ 390 下）', r_of(91, 390) ** 2 - 1 > 3.7, '%.3f' % (r_of(91, 390) ** 2 - 1))
# 高温极限自检：Perfetti Eq.(1) 在 kT >> Omega0 下 (n_e-n_p) ≈ kB(T_e-T_p)/Omega0 ⟹ 与 Allen 首项同形
kB = 0.0861733  # meV/K
Om, Te, Tp = 50.0, 3000.0, 2900.0
ne = 1 / math.expm1(Om / (kB * Te)); npp = 1 / math.expm1(Om / (kB * Tp))
lhs = Om ** 3 * (ne - npp) / Te; rhs = Om ** 2 * kB * (Te - Tp) / Te
chk('Perfetti Eq.(1) 高温极限 → Allen 首项形（λΩ0^2 即 λ<w^2>）', abs(lhs / rhs - 1) < 0.02, '%.4f' % (lhs / rhs))

print('\nALL CHECKS PASSED' if not FAIL else '\nFAILED: %s' % FAIL)
raise SystemExit(1 if FAIL else 0)
