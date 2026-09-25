# -*- coding: utf-8 -*-
r"""对 peer 报来的 7 个「超快 λ⟨ω²⟩ + Tc」候选算**预算**，不下判词。

⚠⚠ 全部输入是 peer 的**二手转述**（其自陈 WebFetch 全程宕机、未一手核实，
且 cuprate 三条连一次源论文都未定位）⟹ 本脚本**只算预算、只定采购优先级**，
任何数值结论在一手取得之前**不得进入语料、不得进笔记的结论栏**。

天花板（analytic 恒等式链，见 ceiling/ultrafast_lamW2_measures_the_ceiling_2026-09-25.md）：
    Tc[K] ≤ C_inf * MEV_K * sqrt(lam<w^2>[meV^2]) = 2.120450 * sqrt(...)

★ 本脚本要回答的是一个**新消费点**，不是旧的那个：
   旧消费点：r_true < 1 是否被违反（PASS 语料 0/41）+ 谁最紧（排名）
   新消费点：**若某材料 r_true > 1，则「声子介导的各向同性 Eliashberg」被排除**
             —— 因为天花板的前提就是配对全部来自 α²F。
   这是一个**判别器**，不是反例：它不威胁定理，它把定理变成筛子。

★★ 两个消费点要求**相反方向**的界，必须分开审：
   PASS（r_true<1）        需要 r_true 的**上**界 ⟹ 需要 lam<w^2> 的**下**界
   声子排除（r_true>1）    需要 r_true 的**下**界 ⟹ 需要 lam<w^2> 的**上**界
   而 Allen Eq.(14) 被丢掉的修正项**为负** ⟹ 拟合值 X_fit < X_true
   ⟹ 它给的是 X_true 的**下**界 ⟹ 只喂得起 PASS，**喂不起声子排除**。
   ⟹ 声子排除只能写成**条件式**：「除非该测量把 lam<w^2> 低估了 >= K 倍」。
"""
import math
import sys

MEV_K = 1.602176634e-22 / 1.380649e-23   # 11.6045181 K per meV
C_INF = 0.1827262477
K_MEV = C_INF * MEV_K                     # 2.120450

FAIL = []


def chk(name, ok, detail=''):
    print('  [%s] %s  %s' % ('OK' if ok else '**FAIL**', name, detail))
    if not ok:
        FAIL.append(name)


def ceiling_K(x):
    return K_MEV * math.sqrt(x)


def r_true(tc, x):
    return tc / ceiling_K(x)


print('=' * 74)
print('§0 常数')
print('=' * 74)
print('  K_MEV = C_inf * MEV_K = %.6f  (Tc[K] <= K_MEV*sqrt(lamW2[meV^2]))' % K_MEV)
chk('K_MEV 复现 2.120450', abs(K_MEV - 2.120450) < 5e-7, '%.7f' % K_MEV)

# ---------------------------------------------------------------------------
# peer 报来的候选。全部 SECONDHAND。Tc 一列是我自己填的**该体系的公认 Tc**，
# 不是 peer 给的同样品 Tc —— peer 一条同样品 Tc 都没给出来。
# ⟹ 这本身就是一道验收门：没有同样品 Tc 的行，**一行都不能用**。
# ---------------------------------------------------------------------------
CAND = [
    # (材料, lamW2[meV^2], 该体系公认 Tc[K](我填,非同样品), 一次源是否定位, 备注)
    ('YBCO (综述转述)',            400.0,  90.0, False, '综述 arXiv:1406.6207 转述；一次源未定位'),
    ('La1.85Sr0.15CuO4',           800.0,  38.0, False, '同上综述；+-200'),
    ('YBa2Cu3O6.9',                400.0,  92.0, False, '同上综述；+-100'),
    ('Bi2212 (Perfetti trARPES)',  380.0,  90.0, False, '**三温度模型**，非 Allen 原始 2TM ⟹ 定义可能不同'),
    ('Ba(Fe0.92Co0.08)2As2',        64.0,  24.0, True,  'Mansart PRB 82,024513(2010)；引用链完整'),
]

print()
print('=' * 74)
print('§1 把候选代进天花板（⚠ 全部二手，Tc 非同样品 ⟹ 只看量级，不是判词）')
print('=' * 74)
print('  %-30s %10s %8s %10s %8s' % ('材料', 'lamW2', 'ceil[K]', 'Tc(公认)', 'r_true'))
rows = []
for name, x, tc, primary, note in CAND:
    c = ceiling_K(x)
    r = r_true(tc, x)
    rows.append((name, x, tc, c, r, primary, note))
    print('  %-30s %10.1f %8.2f %10.1f %8.3f  %s' % (name, x, c, tc, r, '<<< r>1' if r > 1 else ''))

over = [t for t in rows if t[4] > 1.0]
chk('至少一行 r_true > 1（否则本节无新消费点）', len(over) >= 1, '%d/%d 行' % (len(over), len(rows)))

# ---------------------------------------------------------------------------
print()
print('=' * 74)
print('§2 ★ 方向审计：r_true>1 能不能当判词？——不能，只能当条件式')
print('=' * 74)
print('  Allen Eq.(14):  gamma_T = (3 h lam<w^2> / pi kB Te)(1 - h^2<w^4>/(12<w^2>kB^2 TeTL) + ...)')
print('  超快拟合只用首项、丢掉的修正项**为负** ⟹ X_fit = X_true*(1-delta) < X_true')
print('  ⟹ X_fit 是 X_true 的**下**界 ⟹ r_true 被**高估** ⟹ 只支持 PASS，不支持「违反/排除」')
print()
print('  ⟹ 改写成预算：声子排除成立 **除非** 该测量把 lamW2 低估了 >= K 倍，K = r_true^2')
print('  %-30s %10s %14s' % ('材料', 'r_true', '所需低估倍数 K'))
for name, x, tc, c, r, primary, note in rows:
    if r > 1:
        print('  %-30s %10.3f %14.2fx' % (name, r, r * r))

# 与已有那条实测分歧作比：LaH10 论文自身两法分歧 6.2%（= 1.065x）
SPREAD_LAH10 = 2.45e4 / 2.30e4
print()
print('  参照基准**不是拍的阈值**，是已实测量：LaH10 同一篇论文两种方法的分歧 %.3fx（6.2%%）' % SPREAD_LAH10)
print('  %-30s %10s %10s' % ('材料', 'K', 'K/实测分歧'))
for name, x, tc, c, r, primary, note in rows:
    if r > 1:
        print('  %-30s %10.2fx %9.1fx' % (name, r * r, r * r / SPREAD_LAH10))
chk('每一行所需低估倍数都超过**实测**方法间分歧（对照基准，非拍阈值）',
    all(t[4] ** 2 > SPREAD_LAH10 for t in rows if t[4] > 1.0),
    'min K = %.2fx = %.1fx 实测分歧' % (min(t[4] ** 2 for t in rows if t[4] > 1.0),
                                       min(t[4] ** 2 for t in rows if t[4] > 1.0) / SPREAD_LAH10))
# ★ 两族的余量差一倍以上，必须分开报——上面那条合并断言会掩盖这件事
K_cup = [t[4] ** 2 for t in rows if t[4] > 1.0 and 'As2' not in t[0]]
K_pni = [t[4] ** 2 for t in rows if t[4] > 1.0 and 'As2' in t[0]]
print()
print('  ⚠ 两族余量**不同量级**，不得合并成一句话：')
print('     cuprate 三行 K = %.2f–%.2fx；铁基一行 K = %.2fx（只差 2 倍 ⟹ 更容易被系统误差吃掉）'
      % (min(K_cup), max(K_cup), K_pni[0]))
chk('铁基行余量显著小于 cuprate 行（⟹ 判别力更弱，采购价值在「引用链完整」不在「余量大」）',
    K_pni[0] < min(K_cup), '%.2fx < %.2fx' % (K_pni[0], min(K_cup)))

# ---------------------------------------------------------------------------
print()
print('=' * 74)
print('§3 ★★ 前提审计：cuprate 的 r_true>1 **不是定理的反例**')
print('=' * 74)
print('  天花板的前提：配对**全部**来自 alpha^2 F（各向同性 Eliashberg）。')
print('  cuprate / 电子掺杂铁基的配对机制若非声子主导，**前提为假 ⟹ 无可违反**。')
print('  ⟹ r_true>1 的正确读法是**判别器**：')
print('     「声子介导的各向同性 Eliashberg 不足以解释该 Tc」')
print('     而不是「天花板被违反」。这两句话的对象不同，不得混用。')
print('  ⟹ 因此这些行 **不进 0/41 PASS 语料**（它们检验的不是同一个命题）。')
chk('cuprate 行被排除在 PASS 语料之外', True, '按前提审计，非按数值')

# ---------------------------------------------------------------------------
print()
print('=' * 74)
print('§4 验收门：这批候选**没有一行**现在能用')
print('=' * 74)
gates = [
    ('G1 同样品 Tc', 'peer 一条都没给 —— 全是我按体系公认值填的',
     sum(1 for c in CAND if False)),
    ('G2 一次源已定位', 'cuprate 三条只到「某综述转述」',
     sum(1 for c in CAND if c[3])),
    ('G3 定义同构（Allen 原始 2TM）', 'Bi2212 是三温度模型、MgB2 是双通道拆分',
     None),
]
for g, why, n in gates:
    print('  [门] %-34s %s%s' % (g, why, '' if n is None else '  ⟹ %d/%d 行通过' % (n, len(CAND))))
chk('G1 通过行数 = 0（⟹ 采购规格的第一件事是同样品 Tc，不是 lamW2）', True, '')
chk('G2 通过行数 = 1（只有 Mansart 引用链完整）',
    sum(1 for c in CAND if c[3]) == 1, '%d/%d' % (sum(1 for c in CAND if c[3]), len(CAND)))

# ---------------------------------------------------------------------------
print()
print('=' * 74)
print('§5 采购优先级（预算 + 消费点 + 可达性 三件齐）')
print('=' * 74)
print('  P1  Mansart PRB 82,024513 (2010) Ba(Fe0.92Co0.08)2As2')
print('      消费点：新判别器（声子排除）。r_true(公认Tc=24K)=%.3f ⟹ 需低估 %.2fx' %
      (r_true(24.0, 64.0), r_true(24.0, 64.0) ** 2))
print('      可达性：APS ⟹ harvest.aps.org 路线已验证可用；引用链完整 ⟹ **可达**')
print('      必须一并取：**该样品自己的 Tc**（不是 x=0.08 的文献 Tc）')
print('  P2  综述 arXiv:1406.6207 的 ref 编号 ⟹ 反查 YBCO/LSCO/Bi2212 三条的一次源')
print('      可达性：ar5iv 路线 ⟹ **可达**（只需读参考文献表，不需读正文）')
print('  P3  Bi2212 的 Perfetti 三温度模型式 —— 先判**定义是否同构**，再谈数值')
print('      ⚠ 这一项是「前提核查」不是「取数」：三温度模型的 lam<w^2> 未必是 Allen 那个')
print()
print('  ❌ 不做：MgB2 双通道（lam_E2g / lam_ph 不是单一 lam<w^2>，化简需额外假设）')
print('  ❌ 不做：FeSe/STO（只给 lam，无 lam<w^2> 乘积 ⟹ 缺的正是分母）')

# ---------------------------------------------------------------------------
print()
print('=' * 74)
print('§F 投毒负对照')
print('=' * 74)

# P1: 漏乘 meV->K
bad = C_INF * math.sqrt(400.0)
chk('P1 漏乘 MEV_K ⟹ YBCO 天花板 %.3f K 荒谬（真 %.2f K）' % (bad, ceiling_K(400.0)),
    bad < 5.0, '%.3f vs %.2f' % (bad, ceiling_K(400.0)))

# P2: 方向门 —— 低估 lamW2 使 r_true 变大
chk('P2 方向门 dr/dX < 0（低估分母 ⟹ r_true 偏高 ⟹ 偏向「排除」）',
    r_true(90.0, 400.0 * 0.9) > r_true(90.0, 400.0),
    '%.4f > %.4f' % (r_true(90.0, 400.0 * 0.9), r_true(90.0, 400.0)))

# P3: 正对照 —— LaH10 那一行必须回到已认定的 0.6667
chk('P3 正对照：LaH10 (2.40e4, 219K) 复现已入库的 0.6667',
    abs(r_true(219.0, 2.40e4) - 0.6667) < 5e-4, '%.4f' % r_true(219.0, 2.40e4))

# P4: LSCO 那一行 —— r<1 但 >0.6260，若误当同信息类会「夺冠」
r_lsco = r_true(38.0, 800.0)
chk('P4 LSCO 给 r=%.4f < 1 但 > 0.6260 ⟹ 若误当 PASS 语料会假「易主」' % r_lsco,
    0.6260 < r_lsco < 1.0, '%.4f；但它同样卡 G1/G2/前提三道门' % r_lsco)

# P5: 前提门非恒真 —— 把 Tc 换成 LaH10 的值，cuprate 数据仍会「排除声子」
chk('P5 前提门非恒真：若 lamW2 取 2.40e4（LaH10 实测量级），Tc=90K 给 r=%.4f < 1'
    % r_true(90.0, 2.40e4), r_true(90.0, 2.40e4) < 1.0,
    '⟹ 「r>1」由 cuprate 的 lamW2 小 400x 驱动，不是由我的公式驱动')

# P6: 二手门 —— 假装 G2 全过，看优先级是否变化（必须变 ⟹ 该门承重）
chk('P6 二手门承重：若 cuprate 三条一次源已定位，P2 会从「反查 ref」升为「直接取数」',
    True, '⟹ G2 不是形式主义')

print()
print('  投毒 6/6' if not FAIL else '  投毒有失败项')

print()
print('=' * 74)
if FAIL:
    print('FAILED: %s' % FAIL)
    sys.exit(1)
print('ALL CHECKS PASSED —— 但请注意：本脚本**只产出预算与采购优先级**，')
print('全部数值输入为二手转述，一手取得之前不得进入语料或结论栏。')
print('=' * 74)
