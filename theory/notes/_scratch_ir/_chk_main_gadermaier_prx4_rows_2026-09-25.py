"""
Gadermaier et al., PRX 4, 011056 (2014) 的 cuprate 行：主会话从零重算（不复用 sub 脚本）。
页图：library/pdf/eta_hunt_2026-09-25/gadermaier_prx4_{tableI-p2-2,Bconst-p4-4,eq1-p3-3}.png

该文不印 lambda<w^2>，只印 Table I 的 (Tc, tau_e-ph)；Eq.(1) 1/tau = B lambda<w^2>/T，图注 B = 3h/2pi k_B。
  (0) h 还是 hbar？——用该文**自己**的两处数字判：p.4 "lambda ~ 0.3-0.5 for cuprates"（sqrt<w^2>=40 meV）
      与同组 PRL105 的 LSCO (tau=45 fs ⟹ 800 meV^2)。hbar 读法两处都复现；字面 h 读法差 2pi。
  (1) 逐行 X = (2pi/3) kB T hbar/tau，T 取正文 300 K（Fig.1 图注 295 K 给更小的 X ⟹ r 更大；取 300 是保守侧）
  (2) 判别器方向：需 X 上界 ⟹ 条件式；该文明写 T_e ~ T_L ⟹ Allen/TTM 读法 X_TTM/X_NEM = T_e/(2T_L) ~ 1/2，
      只会把 r 再抬高 ⟹ 模型换算不再是风险项（PRL105 的 [0.667,1.333] 来自高通量一般陈述，不适用于本文）
  (3) 与二手值（580/450/720）、与 PRL105、与 Perfetti 的交叉
"""
import math

K = 2.120450
KB = 0.08617333262   # meV/K
HBAR = 658.2119569   # meV*fs
FAIL = []


def chk(name, ok, info=''):
    print('  [%s] %s  %s' % ('OK' if ok else 'FAIL', name, info))
    if not ok:
        FAIL.append(name)


def X_nem(tau_fs, T=300.0):
    return (2 * math.pi / 3) * KB * T * HBAR / tau_fs


def r_of(tc, X):
    return tc / (K * math.sqrt(X))


# Table I（页图逐字；Tc 为该文自己印的样品 Tc）
ROWS = [('YBa2Cu3O6.5', 63, 72), ('YBa2Cu3O6.9', 90, 77), ('HgBa2CuO4.1', 98, 62),
        ('Bi2Sr2CaCu2O8.14', 80, 49), ('La1.9Sr0.1CuO4', 30, 42), ('La1.85Sr0.15CuO4', 38, 45)]

print('== (0) h vs hbar')
lam_hbar = [X_nem(t) / 40.0 ** 2 for _, _, t in ROWS]
lam_h = [x / (2 * math.pi) for x in lam_hbar]
print('   lambda (hbar 读法) = %s' % ', '.join('%.3f' % v for v in lam_hbar))
# 容差 = 印刷精度：「0.3-0.5」一位有效数字 ⟹ 舍入区间 [0.25, 0.55)（初版拍 ±0.02 被 0.530 判 FAIL——拍的阈值，已换）
chk('hbar 读法：cuprate lambda 落在该文自述 0.3-0.5 的一位有效数字舍入区间', all(0.25 <= v < 0.55 for v in lam_hbar))
chk('字面 h 读法：lambda 全 <0.1 ⟹ 与该文自述矛盾（投毒应被拒）', all(v < 0.1 for v in lam_h))
# 同上：「800 +- 200」一位有效数字 ⟹ [750, 850)；且 PRL105 自己印的 tau=45+-8 fs
chk('同组 PRL105：LSCO tau=45 fs ⟹ 落在 800 的印刷舍入区间 [750,850)（hbar 读法）', 750 <= X_nem(45) < 850, '%.1f' % X_nem(45))
chk('字面 h 读法下 PRL105 LSCO = %.1f，远在区间外（投毒）' % (X_nem(45) / (2 * math.pi)), not (750 <= X_nem(45) / (2 * math.pi) < 850))
# 注：PRL105 YBCO O6.5 的 tau_a=100+-20 fs 给 356，印 400+-100；在其误差内但不是逐位复现
print('   PRL105 YBCO tau=100 fs -> %.1f（印 400+-100）' % X_nem(100))

print('\n== (1) 逐行')
res = {}
for name, tc, tau in ROWS:
    X = X_nem(tau)
    Xlo, Xhi = X_nem(tau + 0.5), X_nem(tau - 0.5)     # 两位有效数字的舍入区间
    r, r_cons = r_of(tc, X), r_of(tc, Xhi)
    res[name] = (X, r, r_cons)
    print('   %-17s Tc=%3d tau=%3d  X=%6.1f [%6.1f,%6.1f]  r=%.4f  r(舍入上端X)=%.4f  需低估 r^2=%.3f'
          % (name, tc, tau, X, Xlo, Xhi, r, r_cons, r_cons ** 2))
chk('YBCO O6.9 / Hg1201 在舍入上端 X 下 r>1.9', res['YBa2Cu3O6.9'][2] > 1.9 and res['HgBa2CuO4.1'][2] > 1.9)
chk('两个 LSCO r<1 ⟹ 无判别力', res['La1.9Sr0.1CuO4'][1] < 1 and res['La1.85Sr0.15CuO4'][1] < 1)
chk('T=295 K 给更小 X（取 300 K 是保守侧）', X_nem(77, 295) < X_nem(77, 300))

print('\n== (2) 模型换算方向（该文 T_e ~ T_L）')
conv = 300.0 / (2 * 300.0)
print('   X_TTM/X_NEM = T_e/(2T_L) = %.3f ⟹ TTM 读法 r = %.4f（YBCO O6.9）' % (conv, r_of(90, res['YBa2Cu3O6.9'][0] * conv)))
chk('T_e~T_L 下 Allen/TTM 读法只抬高 r', conv < 1)

print('\n== (3) 交叉')
for name, v in [('HgBa2CuO4.1', 580), ('YBa2Cu3O6.9', 450), ('Bi2Sr2CaCu2O8.14', 720)]:
    print('   二手 %-17s %4d vs 反解 %.1f (%+.1f%%)' % (name, v, res[name][0], 100 * (v / res[name][0] - 1)))
chk('三个二手值都在反解值 3% 内', all(abs(v / res[n][0] - 1) < 0.03 for n, v in
                               [('HgBa2CuO4.1', 580), ('YBa2Cu3O6.9', 450), ('Bi2Sr2CaCu2O8.14', 720)]))
Xb = res['Bi2Sr2CaCu2O8.14'][0]
print('   Bi2212：Gadermaier(Tc 80, 全速率) X=%.1f vs Perfetti(Tc 91, 强耦合子集) 360 ⟹ 比 %.3f' % (Xb, Xb / 360))
print('   若把 %.1f 移植到 Perfetti 样品 (Tc 91)：r=%.4f、需低估 %.3f' % (Xb, r_of(91, Xb), r_of(91, Xb) ** 2))
chk('跨样品比 <3.7（Perfetti 失效条件未被本数据触发，但余量缩小）', Xb / 360 - 1 < 3.7)
chk('YBCO O6.5：PRX4(72 fs, Tc 63) 与 PRL105(100+-20 fs, Tc 60) 的 tau 不一致 >1σ', abs(72 - 100) > 20)

print('\nALL CHECKS PASSED' if not FAIL else '\nFAILED: %s' % FAIL)
raise SystemExit(1 if FAIL else 0)
