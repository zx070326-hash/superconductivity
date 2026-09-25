"""
Nb-Ge 隧穿 alpha^2F（Kihlstrom, Mael, Geballe, PRB 29, 150 (1984)）喂天花板 + RMTA 偏差 kappa。
主会话 2026-09-25 从零重算（不复用 agent 数字）。页图：
  library/pdf/eta_hunt_2026-09-25/kihlstrom1984_prb29p150_tableII-5.png   Table II(a) 7 行（16.7–23.9 at.% Ge）
  library/pdf/eta_hunt_2026-09-25/kihlstrom1984_prb29p150_tableIIb-6.png  Table II(b) 6 行（24.3–25.7 at.%，近化学计量）
  library/pdf/eta_hunt_2026-09-25/kihlstrom1984_prb29p150_moments-p3-3.png <omega>、omega_log 的定义式（lambda 加权）

agent 只报了 (a) 7 行，并称 lambda=1.73 的最高 Tc 样品「被作者排除」——假：它就是 (b) 的第 13 列。
消费点：tier0 A15 表 Nb3Ge 头条 r=0.626 用 RMTA 的 X_RMTA=(37.06/K)^2；它在 kappa=X_RMTA/X_true>=2.5518 时翻。

方向类别：隧穿反演（McMillan-Rowell 邻近分析，lambda 与 mu* 一起拟合）**不是单边界** ⟹
这里的 X 是双边估计，不进无条件 PASS 语料；只给 kappa 的指示性读数。

  (1) 转录门：用印出的 lambda, mu*, <omega>, omega_log, <omega^2> 复现该表印出的 McMillan(Dynes) 与 AD Tc，
      输入按印刷精度传播（半个末位；不拍容差）
  (2) G3 尝试：AD 列在 omega2bar := sqrt(<omega^2>) 与投毒 omega2bar := <omega> 下各复现几行
      ——实测判别力约 0（f2 太弱）⟹ 改由 (3b) 预算：只用定义已钉死的 <omega> 走 Jensen，G3 不承重
  (3) 逐行 X、r、kappa
  (4) Freericks Nb3Sn (cond-mat/0201048 Table 1)：agent 的 X=510.3/309.8 = lambda*A^2，即把 A 当 sqrt<omega^2>；
      正文未定义 A；按约定 A = int alpha^2F（面积）则只有 Jensen 下界 X >= 4A^2/lambda。两种读法都算
"""
import itertools
import math

K = 2.120450
MEV_K = 11.6045181
X_RMTA_NB3GE = (37.06 / K) ** 2
X_RMTA_NB3SN = (41.56 / K) ** 2
KAPPA_FLIP_NB3GE = 2.5518
FAIL = []


def chk(name, ok, info=''):
    print('  [%s] %s  %s' % ('OK' if ok else 'FAIL', name, info))
    if not ok:
        FAIL.append(name)


def expo(lam, mu):
    return math.exp(-1.04 * (1 + lam) / (lam - mu * (1 + 0.62 * lam)))


def tc_mcm(lam, mu, w1):
    return w1 * MEV_K / 1.2 * expo(lam, mu)


def tc_ad(lam, mu, wlog, w2bar):
    L1 = 2.46 * (1 + 3.8 * mu)
    L2 = 1.82 * (1 + 6.3 * mu) * (w2bar / wlog)
    f1 = (1 + (lam / L1) ** 1.5) ** (1 / 3)
    f2 = 1 + (w2bar / wlog - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)   # AD75 f2：无开方
    return f1 * f2 * wlog * MEV_K / 1.2 * expo(lam, mu)


def half_ulp(s):
    return 0.5 * 10 ** (-len(s.split('.')[1])) if '.' in s else 0.5


def box(strs):
    """印刷精度盒：每个输入 +-半个末位；返回所有角点"""
    vals = [(float(s), half_ulp(s)) for s in strs]
    return [tuple(v + sg * h for (v, h), sg in zip(vals, sgn)) for sgn in itertools.product((-1, 1), repeat=len(vals))]


# 页图逐字（字符串保留印刷位数）：id, at.%, Tc(实测), lam(+-), mu*, <w>, wlog, <w^2>(+-), Tc_McM, Tc_AD, Tc_a2F
ROWS = [
    ('1', 16.7, '7.0', '0.81', 0.11, '0.13', '14.3', '12.74', '241', 10, '6.4', '5.9', '5.6'),
    ('2', 21.3, '13.2', '0.85', 0.1, '0.05', '13.7', '12.2', '224', 10, '11.3', '10.7', '10.7'),
    ('3', 21.3, '13.3', '0.91', 0.1, '0.06', '13.7', '12.1', '222', 10, '11.4', '10.3', '10.6'),
    ('4', 21.8, '16.8', '1.05', 0.12, '0.05', '13.8', '12.1', '229', 10, '14.6', '14.0', '14.1'),
    ('5', 22.3, '18.2', '1.40', 0.07, '0.09', '12.7', '10.8', '206', 10, '16.1', '15.4', '15.4'),
    ('6', 22.8, '19.5', '1.44', 0.15, '0.07', '12.4', '10.5', '195', 10, '17.2', '16.7', '17.15'),
    ('7', 23.9, '19.8', '1.64', 0.2, '0.12', '12.7', '10.6', '207', 10, '17.3', '16.6', '16.5'),
    ('8', 24.3, '20.1', '1.46', 0.1, '0.06', '12.1', '10.2', '189', 8, '17.6', '17.2', '17.6'),
    ('9', 24.7, '20.3', '1.40', 0.15, '0.04', '12.4', '10.5', '195', 8, '18.05', '17.8', '18.3'),
    ('10', 25.7, '20.1', '1.55', 0.11, '0.07', '12.1', '10.25', '188', 8, '17.85', '17.6', '17.9'),
    ('11', 25.7, '19.8', '1.50', 0.1, '0.05', '11.7', '9.9', '176', 8, '18.0', '17.9', '18.5'),
    ('12', 24.7, '21.2', '1.70', 0.12, '0.11', '12.4', '10.4', '196', 8, '17.85', '17.45', '17.6'),
    ('13', 24.7, '21.2', '1.73', 0.1, '0.12', '12.6', '10.6', '208', 8, '18.1', '17.6', '17.7'),
]

print('== (0) Jensen 链 omega_log <= <omega> <= sqrt<omega^2>（同一 lambda 加权测度的必要条件）')
chk('13/13 行满足', all(float(r[7]) <= float(r[6]) <= math.sqrt(float(r[8])) for r in ROWS))

print('\n== (1)(2) 转录门 + G3：复现印出的 McMillan / AD Tc（输入印刷精度盒 vs 输出半末位）')
n_mcm = n_ad = n_ad_poison = 0
for r in ROWS:
    rid, lam, mu, w1, wl, w2 = r[0], r[3], r[5], r[6], r[7], r[8]
    corners = box([lam, mu, w1, wl, w2])
    m = [tc_mcm(c[0], c[1], c[2]) for c in corners]
    a = [tc_ad(c[0], c[1], c[3], math.sqrt(c[4])) for c in corners]
    p = [tc_ad(c[0], c[1], c[3], c[2]) for c in corners]          # 投毒：omega2bar := <omega>
    pm, pa = float(r[10]), float(r[11])
    hm, ha = half_ulp(r[10]), half_ulp(r[11])
    ok_m = min(m) <= pm + hm and max(m) >= pm - hm
    ok_a = min(a) <= pa + ha and max(a) >= pa - ha
    ok_p = min(p) <= pa + ha and max(p) >= pa - ha
    n_mcm += ok_m; n_ad += ok_a; n_ad_poison += ok_p
    print('   #%-2s McM 印 %5s 算 [%6.2f,%6.2f] %s | AD 印 %5s 算 [%6.2f,%6.2f] %s | 投毒 [%6.2f,%6.2f] %s'
          % (rid, r[10], min(m), max(m), 'ok' if ok_m else '--', r[11], min(a), max(a), 'ok' if ok_a else '--',
             min(p), max(p), 'ok' if ok_p else '--'))
print('   McMillan 复现 %d/13；AD(sqrt<w^2>) 复现 %d/13；AD 投毒(<w>) 复现 %d/13' % (n_mcm, n_ad, n_ad_poison))
chk('转录门：McMillan 列复现 >= 11/13（lambda、mu*、<omega> 转录无系统错）', n_mcm >= 11, '%d/13' % n_mcm)
# 实测 12 vs 11：f2-1 只有几个百分点，被输入印刷精度盒吞掉 ⟹ 该门对 <omega^2> 的定义**判别力约为 0**。
# 不据此宣布 G3 通过；改走 (3b)：只用页 3 印出定义式的 <omega>（lambda 加权），G3 不承重。
chk('G3 复现门判别力近 0（如实记录，不当证据）：|AD - 投毒| <= 2 行', abs(n_ad - n_ad_poison) <= 2, '%d vs %d' % (n_ad, n_ad_poison))

print('\n== (3) 天花板行：X = lambda<omega^2>，r = Tc/(K sqrt X)，kappa = X_RMTA/X')
print('   X_RMTA(Nb3Ge) = %.1f meV^2（tier0 ceil_min 37.06 K 反解）；头条翻转 kappa = %.4f' % (X_RMTA_NB3GE, KAPPA_FLIP_NB3GE))
res = {}
for r in ROWS:
    rid, at, tc, lam, dl, w2, dw2 = r[0], r[1], float(r[2]), float(r[3]), r[4], float(r[8]), r[9]
    X = lam * w2
    Xlo, Xhi = (lam - dl) * (w2 - dw2), (lam + dl) * (w2 + dw2)
    res[rid] = (at, tc, X, Xlo, Xhi)
    print('   #%-2s %4.1f at.%%  Tc=%5.1f  X=%6.1f [%6.1f,%6.1f]  r=%.4f [%.4f,%.4f]  kappa=%.3f [%.3f,%.3f]'
          % (rid, at, tc, X, Xlo, Xhi, tc / (K * math.sqrt(X)), tc / (K * math.sqrt(Xhi)), tc / (K * math.sqrt(Xlo)),
             X_RMTA_NB3GE / X, X_RMTA_NB3GE / Xhi, X_RMTA_NB3GE / Xlo))
near = [k for k, v in res.items() if v[0] >= 24.3]
chk('近化学计量 (b) 6 行在印出误差棒下限 X 处 r 仍 < 1', all(res[k][1] / (K * math.sqrt(res[k][3])) < 1 for k in near))
kap_max = max(X_RMTA_NB3GE / res[k][3] for k in near)
chk('近化学计量 6 行 kappa 上端（X 取误差棒下限）< 翻转值 2.5518', kap_max < KAPPA_FLIP_NB3GE, 'max %.3f' % kap_max)
X13 = res['13'][2]
print('   把 #13 的 X=%.1f 配头条 Tc=23.2（G1 错配，仅示意）：r=%.4f' % (X13, 23.2 / (K * math.sqrt(X13))))
chk('agent 称 lambda=1.73 样品被排除 —— 假：它是 (b) #13', res['13'][2] == 1.73 * 208)

print('\n== (3b) G3 预算：不用 <omega^2>，只用定义已一手钉死的 <omega>=(2/lambda) int alpha^2F')
print('   Jensen：<omega^2>_P >= <omega>_P^2 ⟹ X >= lambda<omega>^2；取 lambda、<omega> 误差棒下端（使 kappa 上界最大）')
kj = {}
for r in ROWS:
    if r[1] < 24.3:
        continue
    lam_lo = float(r[3]) - r[4]
    w1_lo = float(r[6]) - 0.1                 # (b) 表头 <omega> (+-0.1)
    Xj = lam_lo * w1_lo ** 2
    kj[r[0]] = X_RMTA_NB3GE / Xj
    print('   #%-2s X >= %6.1f ⟹ kappa <= %.3f ；r(实测 Tc) <= %.4f' % (r[0], Xj, kj[r[0]], float(r[2]) / (K * math.sqrt(Xj))))
chk('只用 lambda 加权 <omega>：近化学计量 6 行 kappa 上界全 < 2.5518 ⟹ G3 预算为 ∞（不承重）',
    max(kj.values()) < KAPPA_FLIP_NB3GE, 'max %.3f' % max(kj.values()))

print('\n== (4) Freericks Nb3Sn Table 1（lambda, mu*, omega_ln, A, Tc）')
FR = [('nonconstant DOS', 2.738, 13.652), ('constant DOS', 2.501, 11.130)]
print('   X_RMTA(Nb3Sn) = %.1f；头条 r=0.433 翻转需 kappa >= %.3f' % (X_RMTA_NB3SN, 1 / 0.433 ** 2))
for lab, lam, A in FR:
    X_agent = lam * A ** 2
    X_jlo = 4 * A ** 2 / lam            # A = lambda<omega>/2 ⟹ <omega> = 2A/lambda；Jensen <omega^2> >= <omega>^2
    print('   %-16s agent X=lam*A^2=%6.1f（按面积读法 A=int alpha^2F）；Jensen 下界 X >= %6.1f ⟹ kappa <= %.3f'
          % (lab, X_agent, X_jlo, X_RMTA_NB3SN / X_jlo))
chk('agent 的 510.3/309.8 恰为 lam*A^2（确认它把 A 读作 sqrt<w^2>）',
    abs(FR[0][1] * FR[0][2] ** 2 - 510.3) < 0.05 and abs(FR[1][1] * FR[1][2] ** 2 - 309.8) < 0.05)
# ⚠ 该文正文（ar5iv 文本层）**没有给 A 的定义式**；「A = 面积」是 Carbotte 学派的约定读法，不是一手钉死。
# ⟹ 两种读法都算：agent 读法（A ≡ sqrt<w^2>）与面积读法（Jensen 下界）。判词只取两者都成立的部分。
chk('Jensen 下界（面积读法）下两种反演 kappa 上界均 < Nb3Sn 翻转值', all(X_RMTA_NB3SN / (4 * A * A / lam) < 1 / 0.433 ** 2 for _, lam, A in FR))
chk('agent 读法（A = sqrt<w^2>）下 kappa 同样 < 翻转值 ⟹ 判词对 A 的定义不敏感',
    all(X_RMTA_NB3SN / (lam * A * A) < 1 / 0.433 ** 2 for _, lam, A in FR),
    ', '.join('%.3f' % (X_RMTA_NB3SN / (lam * A * A)) for _, lam, A in FR))
chk('面积读法自洽：<omega> = 2A/lambda >= omega_ln', all(2 * A / lam >= wl for (_, lam, A), wl in zip(FR, (7.072, 6.415))))

print('\nALL CHECKS PASSED' if not FAIL else '\nFAILED: %s' % FAIL)
raise SystemExit(1 if FAIL else 0)
