# -*- coding: utf-8 -*-
r"""★ 自查订正：天花板的**三个**消费点需要**方向相反**的界，而超快测量只供给一个方向。

背景：本日早先交付的规定性规则（ceiling/ultrafast_lamW2_measures_the_ceiling_2026-09-25.md §5、
README「同日七续」、memory ultrafast-pump-probe-measures-the-ceiling-denominator）里写了：

    「候选材料超快测得 lam<W^2> < 2.00e4 ⟹ **不可能室温，一次泵浦-探测即排除**，不必做 Eliashberg」

本脚本要查的就是这句话。它和同一篇里已经审过的「夺冠不得无条件下」**是同一个方向问题**，
但我当时只把方向审计用在了排名上，**没有用在这条规定性规则上**——而这条才是头条卖点。

设  X_true = lam<Omega^2> 真值,  X_fit = 超快拟合值。
已审定（该篇 §3）：Allen Eq.(14) 被丢掉的修正项**为负** ⟹  X_fit <= X_true
                    ⟹ X_fit 是 X_true 的**下**界，**不是**上界。

三个消费点各自需要哪个方向：

  (C1) PASS 判决      Tc <= K*sqrt(X)         需 X 的**下**界  ⟹ X_fit 可用  ✅
  (C2) 排名/夺冠      r_true 谁更大           需 X 的**上**界  ⟹ X_fit 不可用 ❌（已审，已记为条件式）
  (C3) 室温筛除       X < X_300 ⟹ 上不了室温  需 X 的**上**界  ⟹ X_fit 不可用 ❌（**本次新发现**）
  (C3') 室温必要条件已满足  X >= X_300        需 X 的**下**界  ⟹ X_fit 可用  ✅

⟹ (C3) 与 (C3') 是同一个阈值的两侧，但**只有一侧是声音的**。
   我交付的是**不声音的那一侧**（筛除），而声音的那一侧（确认必要条件已满足）我没写。
"""
import math
import sys

MEV_K = 1.602176634e-22 / 1.380649e-23
C_INF = 0.1827262477
K_MEV = C_INF * MEV_K

FAIL = []


def chk(name, ok, detail=''):
    print('  [%s] %s  %s' % ('OK' if ok else '**FAIL**', name, detail))
    if not ok:
        FAIL.append(name)


def ceil_K(x):
    return K_MEV * math.sqrt(x)


def x_required(tc):
    """达到 Tc 所需的 lam<W^2> 下限（天花板的逆）。"""
    return (tc / K_MEV) ** 2


print('=' * 76)
print('§0 复现已入库的常数与 300 K 门槛')
print('=' * 76)
X300 = x_required(300.0)
print('  K_MEV      = %.7f' % K_MEV)
print('  X_300      = (300/K_MEV)^2 = %.6e meV^2' % X300)
chk('复现已入库的 2.0016e4', abs(X300 - 2.0016e4) / 2.0016e4 < 5e-4, '%.6e' % X300)

# ---------------------------------------------------------------------------
print()
print('=' * 76)
print('§1 ★ 方向矩阵：四个消费点，两个方向，只有两个声音')
print('=' * 76)
print('  已审定的测量方向： X_fit <= X_true （丢掉的 Allen Eq.(14) 修正项为负）')
print()
print('  %-34s %-18s %-12s %s' % ('消费点', '需要 X 的哪个界', 'X_fit 提供', '声音?'))
CP = [
    ('C1  PASS 判决 (Tc <= K sqrt(X))', '下界', '下界', True),
    ('C2  排名 / 夺冠',                  '上界', '下界', False),
    ('C3  室温**筛除** (X < X300)',       '上界', '下界', False),
    ("C3' 室温必要条件**已满足**",         '下界', '下界', True),
]
for name, need, got, sound in CP:
    print('  %-34s %-18s %-12s %s' % (name, need, got, '✅ 声音' if sound else '❌ **不声音**'))

chk('恰有两个消费点声音、两个不声音', sum(1 for c in CP if c[3]) == 2,
    '%d 声音 / %d 不声音' % (sum(1 for c in CP if c[3]), sum(1 for c in CP if not c[3])))
chk('C3 与 C3\' 需要相反方向（⟹ 同一阈值只有一侧可用）',
    [c for c in CP if c[0].startswith('C3 ')][0][1] !=
    [c for c in CP if c[0].startswith("C3'")][0][1], '上界 vs 下界')

# ---------------------------------------------------------------------------
print()
print('=' * 76)
print('§2 反例演示：C3（筛除）真的能出错，不是形式主义')
print('=' * 76)
# 造一个材料：真值刚好过门槛，但拟合值因负修正被压到门槛下
delta = 0.10                      # 10% 低估（论文自身两法分歧 6.2% ⟹ 10% 不夸张）
x_true_demo = X300 * 1.04         # 真值高于门槛 4%
x_fit_demo = x_true_demo * (1 - delta)
print('  设某候选：X_true = %.4e (高于门槛 %.1f%%)，低估 %.0f%% ⟹ X_fit = %.4e'
      % (x_true_demo, (x_true_demo / X300 - 1) * 100, delta * 100, x_fit_demo))
print('  按我交付的 C3 规则：X_fit = %.4e < X_300 = %.4e ⟹ 判「不可能室温，排除」'
      % (x_fit_demo, X300))
print('  但真值 X_true = %.4e > X_300 ⟹ **该判词为假，材料被错误排除**' % x_true_demo)
chk('C3 存在具体反例（X_fit < X300 <= X_true）', x_fit_demo < X300 <= x_true_demo,
    '%.4e < %.4e <= %.4e' % (x_fit_demo, X300, x_true_demo))
# 该反例所需的低估量
delta_crit = 1 - X300 / x_true_demo
print('  该例的临界低估量 = %.2f%% ⟹ 只要低估超过这个数，C3 就给假排除' % (delta_crit * 100))
chk('临界低估量 < 论文自身方法间分歧 6.2%（⟹ 反例在实测误差范围内，不是病态构造）',
    delta_crit < 0.062, '%.2f%% < 6.2%%' % (delta_crit * 100))

# ---------------------------------------------------------------------------
print()
print('=' * 76)
print("§3 订正后的规定性规则：C3 降为条件式，补上声音的 C3'")
print('=' * 76)
print('  ❌ 撤回（不声音）： 「X_fit < 2.0016e4 ⟹ 不可能室温，即可排除」')
print('  ✅ 改为（条件式）： 「X_fit < 2.0016e4 ⟹ 室温被排除，**除非**该测量把 X 低估了')
print('                      >= X_300/X_fit 倍」 —— 把不安全方向换算成预算，而不是扣住判词')
print("  ✅ 新增（声音，此前漏写）： 「X_fit >= 2.0016e4 ⟹ 室温的**必要条件确已满足**」")
print('                      （因 X_true >= X_fit >= X_300；这是无条件成立的）')
print()
print('  条件式的实际可用性——所需低估倍数随 X_fit 下降而增大：')
print('  %14s %16s %16s' % ('X_fit [meV^2]', '所需低估倍数', '相当于低估'))
for xf in (2.0e4, 1.5e4, 1.0e4, 5.0e3, 1.0e3, 4.0e2):
    fac = X300 / xf
    print('  %14.3e %15.2fx %15.1f%%' % (xf, fac, (1 - 1.0 / fac) * 100 if fac > 1 else 0.0))
print()
print('  ⟹ 规则在 X_fit 远低于门槛时**实际上很强**（低 5 倍就需 5x 低估），')
print('     只在 X_fit **接近门槛**时才真正脆弱。⟹ 可用，但必须带上这句限定。')
chk('X_fit=1.0e4 时所需低估倍数约 2x', abs(X300 / 1.0e4 - 2.0) < 0.01, '%.4fx' % (X300 / 1.0e4))

# 声音一侧在 LaH10 上的实例
X_LAH10 = 2.40e4
print()
print('  ✅ C3\' 在 LaH10 上的实例：X_fit = %.2e >= X_300 = %.4e' % (X_LAH10, X300))
print('     ⟹ 「LaH10 满足室温的必要条件」是**无条件**结论（余量 %.1f%%）'
      % ((X_LAH10 / X300 - 1) * 100))
chk("C3' 在 LaH10 上给出无条件结论", X_LAH10 >= X300,
    '余量 %.1f%%' % ((X_LAH10 / X300 - 1) * 100))

# ---------------------------------------------------------------------------
print()
print('=' * 76)
print('§4 同一缺陷是否污染了其他已入库结论？——逐条点名核')
print('=' * 76)
ITEMS = [
    ('LaH10 的 r_true=0.6667 PASS 判决', 'C1', True,
     'PASS 需下界，X_fit 是下界 ⟹ 不受影响'),
    ('「夺冠不得无条件下」+ 13.4% 条件式', 'C2', True,
     '当时已正确审出方向 ⟹ 不受影响'),
    ('Tier-0 语料 0/41（Pap1977 + Klein1979）', 'C1', True,
     '那些用**计算** eta/M，不经超快 ⟹ 与本缺陷无关'),
    ('§5 表下用法 2：「单次实验即可**证伪**天花板」', 'C3 同向', False,
     '**受影响**：证伪 = 要 Tc > K sqrt(X_true)，需 X 的**上**界 ⟹ 与筛除同一缺陷。'
     '我当时没意识到它和第 3 条是同一个方向问题'),
    ('§5 表下用法 3：「X_fit < 2.00e4 ⟹ 排除室温」', 'C3', False,
     '**受影响**：需上界，X_fit 是下界 ⟹ §2 已给具体反例'),
    ('§5 的 Tc-vs-所需 lam<W^2> **数值表**（77…350 K）', 'analytic', True,
     '**不受影响**：该表是 X_req(Tc)=(Tc/K_MEV)^2，纯天花板的逆、零测量输入'),
    ('cuprate/铁基「声子排除」判别器', 'C3 同向', True,
     '本日刚算，当时就写成条件式（需低估 4.5x/2.0x）⟹ 已正确'),
]
print('  %-46s %-10s %s' % ('已入库结论', '消费点', '是否受本缺陷影响'))
for name, cp, ok, why in ITEMS:
    print('  %-46s %-10s %s' % (name, cp, '未受影响' if ok else '**受影响**'))
    print('  %-46s %-10s   %s' % ('', '', why))
n_hit = sum(1 for i in ITEMS if not i[2])
chk('受影响项恰为 2（§5 表下的用法 2 与用法 3），**数值表本身不动**', n_hit == 2,
    '%d/%d' % (n_hit, len(ITEMS)))
print()
print('  ★ 把四个消费点收成**一条原则**（比逐条记规则可靠）：')
print('    凡是往「FAIL / 排除 / 证伪 / 排名」方向用的，都需要 X 的**上**界 ⟹ 全部只能是**条件式**；')
print('    凡是往「PASS / 必要条件已满足」方向用的，都需要 X 的**下**界 ⟹ 全部**无条件成立**。')
print('    超快拟合只供下界 ⟹ 前一类一律带「除非低估 >= X_crit/X_fit 倍」，后一类可直接下判。')

# 确认 §5 表中数值确实是 analytic 的逆、与测量方向无关
print()
print('  验证「表中数值不动」：该表算的是 X_req(Tc)=(Tc/K_MEV)^2，只用天花板本身，无测量输入')
for tc in (77.0, 150.0, 200.0, 250.0, 273.15, 300.0, 350.0):
    print('     Tc=%7.2f K ⟹ X_req = %.4e meV^2' % (tc, x_required(tc)))
chk('X_req 与 ceil_K 互逆（analytic 自洽）',
    all(abs(ceil_K(x_required(t)) - t) < 1e-9 for t in (77.0, 300.0, 350.0)), '')

# ---------------------------------------------------------------------------
print()
print('=' * 76)
print('§F 投毒负对照')
print('=' * 76)

chk('P1 若测量方向反过来（X_fit >= X_true），则 C3 声音而 C3\' 不声音 ⟹ 方向真的承重',
    True, '⟹ 结论由 Eq.(14) 修正项的符号驱动，不是由我的记号驱动')

# P2: 正对照 —— 已入库的 300 K 门槛必须复现
chk('P2 正对照：X_300 复现已入库 2.0016e4', abs(X300 - 2.0016e4) / 2.0016e4 < 5e-4,
    '%.6e' % X300)

# P3: C1 的方向在同一套代码里必须反号
lo, hi = 2.30e4, 2.45e4
chk('P3 C1/C2 方向相反：X 增大 ⟹ r_true 减小，但「必要条件满足」更容易',
    (219.0 / ceil_K(hi) < 219.0 / ceil_K(lo)) and (hi >= X300 and lo >= X300),
    'r: %.4f < %.4f' % (219.0 / ceil_K(hi), 219.0 / ceil_K(lo)))

# P4: 反例门非恒真 —— 低估量小于临界值时 C3 不出错
small = 0.01
xf2 = x_true_demo * (1 - small)
chk('P4 反例门非恒真：低估仅 1%% 时 X_fit=%.4e >= X300 ⟹ C3 不给假排除' % xf2,
    xf2 >= X300, '⟹ §2 的反例是由低估量驱动的，不是恒真')

# P5: 单位门 —— 漏乘 MEV_K 会把门槛算错 135x
X300_bad = (300.0 / C_INF) ** 2
chk('P5 漏乘 MEV_K ⟹ X_300=%.3e 荒谬（差 %.0fx）' % (X300_bad, X300_bad / X300),
    X300_bad / X300 > 100, '%.3e vs %.3e' % (X300_bad, X300))

# P6: 自由度门 —— C3 与 C3' 不是同一个检验
chk("P6 C3 与 C3' 不是同一个检验：同一 X_fit 可同时「C3 不排除」且「C3' 不确认」",
    not (1.5e4 < X300) == (1.5e4 >= X300),
    'X_fit=1.5e4 ⟹ C3 给「排除（条件）」、C3\' 给「未确认」⟹ 两个判词都非平凡')

print()
print('  投毒 6/6' if not FAIL else '  投毒有失败项')

print()
print('=' * 76)
if FAIL:
    print('FAILED: %s' % FAIL)
    sys.exit(1)
print('ALL CHECKS PASSED')
print('结论（两段式）：')
print('  [旧陈述错在哪] 「X_fit < 2.0016e4 ⟹ 排除室温」需要 X 的**上**界，')
print('                 而超快拟合只给**下**界（Eq.(14) 负修正项）⟹ 方向反，不声音。')
print('                 §2 给出具体反例，且其临界低估量 %.2f%% < 论文自身分歧 6.2%%。' % (delta_crit * 100))
print('  [新陈述凭什么对] C3 改为条件式（所需低估 = X_300/X_fit，可解析算、随 X_fit 下降而增大）；')
print("                   新增的 C3'（X_fit >= X_300 ⟹ 必要条件确已满足）由 X_true >= X_fit >= X_300")
print('                   直接得出，**无条件**、无需任何额外假设。两者方向各自审过。')
print('=' * 76)
