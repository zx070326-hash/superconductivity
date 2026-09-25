# -*- coding: utf-8 -*-
r"""
Wu et al. Nat Commun 15:9683 (2024), PMC11549371 — LaH10±δ @165 GPa 超快泵浦探测
把实测的 lambda<Omega^2> 直接当作 Tier-0 天花板的宗量。

承重恒等式（两端都已取一手页图/仓库一手记录）：
  Allen PRL 59,1460 (1987) Eq.(11):  lambda<w^n> = 2 \int_0^inf dW [a2F(W)/W] W^n
        n=2  =>  lambda<w^2> = 2 \int dW a2F(W) W
  AD75  PRB 12,905 Eq.(5):           <w^n> = (2/lam) \int dW a2F(W) W^{n-1}
        n=2  =>  lam * wbar2^2 = 2 \int dW a2F(W) W
  => lambda<Omega^2>_Wu  ==  lam * wbar2^2  ==  eta/M      [analytic，纯定义链]
  AD75 Eq.(25): Tc = 0.182 (lambda<w^2>)^{1/2}  （lam>~10, mu*=0）——本仓库认证为 C_INF

严格性：本脚本全部为 numerical-evidence（已发表浮点数 + float64），无区间算术。
恒等式识别一节是 analytic。
"""
import math, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

FAIL = []
def chk(name, cond, msg=''):
    print(('  [ OK ] ' if cond else '  [FAIL] ') + name + (('  ' + msg) if msg else ''))
    if not cond:
        FAIL.append(name)

# ---------------- §0 常数 ----------------
MEV_K = 1.602176634e-22 / 1.380649e-23      # K per meV (CODATA 定义值，精确)
C_INF = 0.1827262477                        # 本仓库认证的 lam->inf, mu*=0 天花板常数
NB3GE = 0.6260                              # 目前 0/41 语料里最紧的一行
K_MEV = C_INF * MEV_K                       # Tc[K] <= K_MEV * sqrt(lamW2[meV^2])

print('=' * 74)
print('§0 常数')
print('=' * 74)
print('  meV -> K            : %.7f' % MEV_K)
print('  C_inf               : %.10f' % C_INF)
print('  K_MEV = C_inf*meV_K : %.6f     <= 天花板写成 Tc[K] = K_MEV*sqrt(lamW2[meV^2])' % K_MEV)
# AD75 Eq.(25) 印 0.182。注意 round(C_INF,3) = 0.183 != 0.182 ⟹ **不是四舍五入**。
# 与 floor(截断) 一致——这正是本仓库已记录过的印刷约定（见 quan_pickett 那篇的截断方向讨论）。
# 断言写成截断检验，比放松容差**更锐**：它同时排除了「四舍五入」这一读法。
chk('AD75 印的 0.182 = 截断(C_inf)，且**不等于**四舍五入(C_inf)=0.183',
    math.floor(C_INF * 1000) / 1000 == 0.182 and round(C_INF, 3) == 0.183,
    'floor=%.3f round=%.3f' % (math.floor(C_INF * 1000) / 1000, round(C_INF, 3)))
print('  方向：C_inf 是本仓库**认证的上确界常数**，天花板 Tc <= C_inf*sqrt(lamW2) 用它即定理原文，')
print('        不存在「保守与否」；AD75 的 0.182 只是他们的数值渐近值且**被截断**打印。')
print('        若误用 0.182：分母更小 => r_true 更大 => 会**偏向夺冠**，故不得代换。')

# ---------------- §1 Wu 的实测输入（逐字） ----------------
# 全部取自 EuropePMC 全文 XML（一手，http=200 ssl_verify=0）
LAMW2_T  = 2.45e4   # Fig 3b 温度依赖；式中**只含** lambda<Omega^2>（Allen Eq.14 首项）
LAMW2_F  = 2.30e4   # Fig 4d 通量依赖；式中另含 <Omega> 与 <Omega^3>（需额外矩假设）
LAMW2    = 2.40e4   # 作者取的平均值
LAMW2_SD = 0.10e4   # 注意：这是**两种方法的分歧**，不是统计误差棒
TC_SAME  = 219.0    # 本样品、本实验测得（RT 模型拟合 A_slow/tau_slow）
TC_SAME_E = 4.0
TC_OTHER = [('本样品早期输运（退化前）', 240.0), ('本样品早期输运（退化前）', 250.0),
            ('文献 LaH10 165-190 GPa 下端', 250.0), ('文献 LaH10 上端', 260.0)]
WLOG_163 = 96.4     # Errea @163 的 w_log，Wu 借来做除法
LAM_WU   = 2.58     # Wu 由此报出的 lambda

print()
print('=' * 74)
print('§1 Wu 报的数：那个 lambda 是**除**出来的，不是独立测的')
print('=' * 74)
lam_derived = LAMW2 / WLOG_163 ** 2
sd_derived = LAMW2_SD / WLOG_163 ** 2
print('  lamW2/wlog^2 = %.4e / %.1f^2 = %.4f   (Wu 印 %.2f)' % (LAMW2, WLOG_163, lam_derived, LAM_WU))
print('  sd/wlog^2    = %.4f                      (Wu 印 0.11)' % sd_derived)
chk('逐位复现 Wu 的 lambda = 2.58', abs(lam_derived - LAM_WU) < 0.005,
    'got %.4f' % lam_derived)
chk('逐位复现 Wu 的 +-0.11', abs(sd_derived - 0.11) < 0.005, 'got %.4f' % sd_derived)
print('  => 那个替换是**事后除法**；乘积 lamW2 的拟合本身不经过 w_log ⟹ 我用乘积不继承 @163 的循环性')

# ---------------- §2 天花板：直接用乘积 ----------------
def ceiling_K(lamW2_meV2):
    return K_MEV * math.sqrt(lamW2_meV2)

def r_true(tc, lamW2_meV2):
    return tc / ceiling_K(lamW2_meV2)

print()
print('=' * 74)
print('§2 同一样品配对（唯一合法的那一对）：Tc = 219 +- 4 K')
print('=' * 74)
print('  %-34s %10s %10s' % ('lamW2 [meV^2]', 'ceiling[K]', 'r_true'))
for tag, v in [('Fig3b 温度依赖(矩干净) 2.45e4', LAMW2_T),
               ('作者平均值        2.40e4', LAMW2),
               ('Fig4d 通量依赖(需额外矩) 2.30e4', LAMW2_F)]:
    print('  %-34s %10.2f %10.4f' % (tag, ceiling_K(v), r_true(TC_SAME, v)))

# PASS 判词要 r_true 的**上**界 => max Tc / min 分母
r_up = r_true(TC_SAME + TC_SAME_E, LAMW2 - LAMW2_SD)
# 排名判词要 r_true 的**下**界 => min Tc / max 分母
r_lo = r_true(TC_SAME - TC_SAME_E, LAMW2 + LAMW2_SD)
print()
print('  r_true 区间（两侧角点）: [%.4f, %.4f]' % (r_lo, r_up))
chk('PASS：天花板未被违反（r_true 上界 < 1）', r_up < 1.0, 'r_up = %.4f' % r_up)
chk('r_true 下端 > Nb3Ge 的 0.6260（点上成立，但方向另议）', r_lo > NB3GE,
    'r_lo = %.4f' % r_lo)

# ---------------- §3 方向审计：Allen Eq.(14) 被截断了 ----------------
print()
print('=' * 74)
print('§3 方向审计 —— 这一节决定哪个判词能下、哪个不能')
print('=' * 74)
print('  Allen Eq.(14):  gamma_T = (3 h_bar lam<w2>/pi kB Te) * (1 - d),  d>0')
print('  Wu 的 Fig3b 式 = Eq.(14) 首项（丢掉 d）配 Te = sqrt(TL^2+Theta F)')
print('  给定实测 tau:  X_fit = lam<w2>_true * (1-d)  <  lam<w2>_true')
print('  => 分母被**低**估 => r_true 被**高**估')
print('  ⟹ PASS 判词（r_true<1）保守、合法；排名判词（r_true>0.6260）**偏向夺冠、不得下**')
print()
# 把「不安全」换算成一个可交付的预算：lamW2 要被低估多少，夺冠才失效？
need = (TC_SAME / (K_MEV * NB3GE)) ** 2
print('  夺冠临界: lamW2 需 <= %.4e meV^2 (Tc=219)' % need)
infl = need / LAMW2 - 1.0
print('  即：实测 2.40e4 只要被低估 **%.1f%%**，夺冠即失效' % (infl * 100))
spread = abs(LAMW2_T - LAMW2_F) / LAMW2
print('  论文自身两法分歧 = %.1f%%  （2.45e4 vs 2.30e4）' % (spread * 100))
chk('夺冠余量 > 论文自身两法分歧（否则连条件式都不成立）', infl > spread,
    '%.1f%% vs %.1f%%' % (infl * 100, spread * 100))
print('  ⟹ 判词降级为**条件式**：「若 Wu 的 lamW2 低估不超过 %.1f%%，则该行夺冠」' % (infl * 100))

# ---------------- §4 与 Errea 计算值的对照：一处真实张力 ----------------
print()
print('=' * 74)
print('§4 实验 lamW2 vs 计算 lam*wbar2^2 —— 论文说的 "in line with" 是在比另一个量')
print('=' * 74)
# (P, lam, wlog_meV, r_lo, r_hi)  r = wbar2/wlog，来自我自己的 f2 反解
ERREA = [(129, 3.62,  76.4, 1.3180, 1.3250),
         (163, 2.67,  96.4, 1.1272, 1.1332),   # 已弃用行，仅作诊断
         (214, 2.06, 115.5, 1.2275, 1.2405),
         (264, 1.73, 126.6, 1.2281, 1.2460)]
print('  %5s %7s %9s %12s %12s' % ('P', 'lam', 'wlog', 'lam*wlog^2', 'lam*wbar2^2'))
rows = []
for P, lam, wl, rlo, rhi in ERREA:
    a = lam * wl ** 2
    rm = 0.5 * (rlo + rhi)
    b = lam * (rm * wl) ** 2
    rows.append((P, a, b))
    print('  %5d %7.2f %9.1f %12.4e %12.4e%s' % (P, lam, wl, a, b, '   <- 已弃用' if P == 163 else ''))
a163 = [b for P, a, b in rows if P == 163][0]
print()
print('  Wu 实测 lamW2                 = %.4e' % LAMW2)
print('  Errea @163 的 lam*wlog^2      = %.4e  <- 论文拿来比的就是这个（差 %.1f%%）'
      % (rows[1][1], (rows[1][1] / LAMW2 - 1) * 100))
print('  Errea @163 的 lam*wbar2^2     = %.4e  <- 真正可比的量（差 %.1f%%）'
      % (a163, (a163 / LAMW2 - 1) * 100))
chk('论文的 "in line with" 只在 lam*wlog^2 这一列成立（<10%）',
    abs(rows[1][1] / LAMW2 - 1) < 0.10, '%.1f%%' % ((rows[1][1] / LAMW2 - 1) * 100))
chk('换成真正可比的 lam*wbar2^2 后偏离 >25%',
    abs(a163 / LAMW2 - 1) > 0.25, '%.1f%%' % ((a163 / LAMW2 - 1) * 100))
print('  ⟹ 实验与计算在**同一个量**上相差 30%% 量级，这是一处可被后续实验判死的真张力')

# 自由度审计：@163 在 lam*wbar2^2 上离群，是不是**新**诊断？
seq_w2 = [b for P, a, b in rows]
seq_wl = [a for P, a, b in rows]
mono_wl = all(seq_wl[i] < seq_wl[i + 1] for i in range(3))
mono_w2 = all(seq_w2[i] < seq_w2[i + 1] for i in range(3))
print()
print('  单调性（随 P）: lam*wlog^2 %s ; lam*wbar2^2 %s'
      % ('单调↑' if mono_wl else '非单调', '单调↑' if mono_w2 else '**非单调**'))
chk('lam*wlog^2 对 P 单调', mono_wl)
chk('lam*wbar2^2 因 @163 破坏单调（定位到 r 而非 lam/wlog）', not mono_w2)
print('  ⚠ **自由度**：两列只差一个因子 r^2，且 lam/wlog 那列单调 ⟹ 这不是第四条**独立**诊断，')
print('    而是已知的 "@163 的 r=1.1302 偏低" 换了个量再现一次。按 DOF 规则，不计入独立证据。')

# ---------------- §5 规定性产物：给实验的设计规则 ----------------
print()
print('=' * 74)
print('§5 ★ 规定性产物 —— 零 DFT、零 mu*、单次实验即可判的必要条件')
print('=' * 74)
print('  天花板改写为实验可直接代入的形式：')
print('      Tc[K]  <=  %.6f * sqrt( lambda<Omega^2> [meV^2] )' % K_MEV)
print('      等价：  lambda<Omega^2>  >=  ( Tc / %.6f )^2' % K_MEV)
print()
print('  %10s %22s %18s' % ('目标 Tc[K]', '所需 lamW2 [meV^2]', '所需 sqrt [meV]'))
for tc in [77, 150, 200, 250, 273.15, 300, 350]:
    need_i = (tc / K_MEV) ** 2
    print('  %10.2f %22.4e %18.2f' % (tc, need_i, math.sqrt(need_i)))
print()
need300 = (300.0 / K_MEV) ** 2
print('  室温 300 K 的必要条件: lamW2 >= %.4e meV^2 ; LaH10 实测 %.2e ⟹ %s'
      % (need300, LAMW2, '满足' if LAMW2 >= need300 else '不满足'))
chk('LaH10 实测 lamW2 满足 300 K 的必要条件', LAMW2 >= need300,
    'ratio = %.3f' % (LAMW2 / need300))
print('  ⚠ 这是**必要**条件，不是充分条件（mu*>0 与谱形状只会把 Tc 往下压）')

# ---------------- §6 不作判词的配对（记录） ----------------
print()
print('=' * 74)
print('§6 其它 Tc 配对：**记录，不进任何判词**（跨样品/跨样品态）')
print('=' * 74)
for tag, tc in TC_OTHER:
    print('  Tc=%6.1f K (%s):  r_true = %.4f' % (tc, tag, r_true(tc, LAMW2)))
print('  理由：lamW2 测于**退化后**的样品（本实验测得 Tc=219K，输运初测 240-250K）；')
print('        拿 250/260 K 配退化后的 lamW2 = 混了两个不同对象，且方向**偏向夺冠**。')

# ---------------- §F 投毒负对照 ----------------
print()
print('=' * 74)
print('§F 投毒负对照')
print('=' * 74)
P = []
# P1 漏掉 meV->K
r_p1 = TC_SAME / (C_INF * math.sqrt(LAMW2))
P.append(('P1 漏乘 meV->K', r_p1 > 5, 'r_true = %.3f 荒谬' % r_p1))
# P2 正对照：lam*wlog^2 必须复现 Wu 的 lamW2（证明替换=纯除法）
r_p2 = LAM_WU * WLOG_163 ** 2
P.append(('P2 正对照 lam*wlog^2 复现 lamW2', abs(r_p2 / LAMW2 - 1) < 0.01,
          '%.4e vs %.4e (%.2f%%)' % (r_p2, LAMW2, (r_p2 / LAMW2 - 1) * 100)))
# P3 跨样品配对 => r_true 明显抬高
r_p3 = r_true(260.0, LAMW2)
P.append(('P3 跨样品 Tc=260 抬高 r_true', r_p3 > r_true(TC_SAME, LAMW2) + 0.10,
          '%.4f vs %.4f' % (r_p3, r_true(TC_SAME, LAMW2))))
# P4 改用 Errea 计算的 lam*wbar2^2 => 掉到 0.6260 以下 => 夺冠确实悬于实验值
r_p4 = r_true(TC_SAME, a163)
P.append(('P4 换成计算 lam*wbar2^2 则不夺冠', r_p4 < NB3GE, 'r_true = %.4f' % r_p4))
# P5 方向门：分母增大必须使 r_true 减小
P.append(('P5 方向门 d r_true / d lamW2 < 0',
          r_true(TC_SAME, LAMW2 * 1.10) < r_true(TC_SAME, LAMW2), ''))
# P6 把「点估计」当判词（本篇的方法论正对照）
pt = r_true(TC_SAME, LAMW2)
P.append(('P6 点估计 0.667 > 0.6260 但方向偏向夺冠 ⟹ 不构成排名判词',
          pt > NB3GE and r_lo > NB3GE and infl > 0,
          '点估计 %.4f, 下端 %.4f, 但需低估<%.1f%% 才成立' % (pt, r_lo, infl * 100)))
# P7 单模假设投毒：若 <Omega^2> 被换成 wlog^2（Wu 做的那个替换）当作矩本身
r_p7 = r_true(TC_SAME, LAM_WU * WLOG_163 ** 2)
P.append(('P7 用 lam*wlog^2 冒充 lam<Omega^2> 不改变本行结论（两者数值同源）',
          abs(r_p7 - pt) < 0.01, '%.4f vs %.4f' % (r_p7, pt)))
ok = 0
for n, c, m in P:
    print(('  [ OK ] ' if c else '  [MISS] ') + n + ('  ' + m if m else ''))
    ok += 1 if c else 0
    if not c:
        FAIL.append(n)
print('  投毒 %d/%d' % (ok, len(P)))

print()
print('=' * 74)
if FAIL:
    print('FAILED:', FAIL)
    sys.exit(1)
print('ALL CHECKS PASSED')
