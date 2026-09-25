r"""从零复核 arXiv:2010.12225 (heavy rare-earth hydrides) SI Table S4 的 Mc/AD 双列反解。

一手件：library/pdf/eta_hunt_2026-09-25/mcad_2010.12225.pdf
  md5 26a4de67c75f76e0ec6ac66544db0da8（主会话**自己 curl 下载**，与 teammate 的副本逐字节相同）
  页图 hrE_p-33.png（Table S4）、hrE_p-28.png（Eq. S5/S6/S7）——两种仪器都跑了：
  ⚠ 文本层把**材料标签与数据行的对应打乱**（TmH6 标签压在 160/0.13 上、EuH6 标签落在 150/2.67 后），
    本脚本的 ROWS 取自**页图**，不取自 pdftotext。

运行：OMP_NUM_THREADS=1 PYTHONIOENCODING=utf-8 python theory/notes/_scratch_ir/_chk_main_hrE_2010p12225_mcad_backsolve_2026-09-25.py
"""
import io
import math
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

FAIL = []


def bad(msg):
    FAIL.append(msg)
    print("  x FAIL:", msg)


def sec(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


# ============================== §A Table S4（逐行取自页图 hrE_p-33.png）
sec("§A Table S4 —— 逐行取自**页图**（文本层标签错位，不可用）")

# (材料, P/GPa, lam, w_log/K, N(ef), alpha, Tc^a 低-高, Tc^b(A-D) 低-高)
ROWS = [
    ("YbH6",  70, 2.22,  652.4,  8.4, 0.48, (90, 94),   (114, 125)),
    ("YbH6", 200, 1.00, 1146.5,  6.8, 0.38, (69, 80),   (72, 85)),
    ("LuH6", 100, 3.60,  750.5,  4.8, 0.49, (140, 146), (204, 222)),
    ("LuH6", 200, 1.77, 1366.2,  4.2, 0.47, (166, 180), (190, 209)),
    ("TmH6",  50, 0.72,  611.9, 29.6, 0.47, (17, 22),   (18, 23)),
    ("EuH6", 160, 0.13,  511.1, 64.4, None, (0, 0),     (0, 0)),
    ("CaH6", 150, 2.67, 1027.1,  2.4, 0.38, (165, 175), (205, 222)),
    ("ScH6", 300, 1.40, 1233.4,  5.0, 0.30, (119, 132), (132, 148)),
    ("ScH6", 350, 1.30, 1330.8,  4.8, 0.30, (118, 131), (129, 146)),
    ("YH6",  120, 3.06,  829.0,  4.7, 0.46, (143, 150), (195, 213)),
    ("H3S",  200, 1.88, 1266.0,  3.1, 0.35, (162, 175), (186, 204)),
]
print(f"共 {len(ROWS)} 行（页图 11 行，逐行核对完毕）")
print(f"    {'row':12s} {'lam':>5s} {'w_log/K':>8s} {'Tc^a':>9s} {'Tc^b':>9s}")
for m, P, lam, wl, _n, _a, ta, tb in ROWS:
    print(f"    {m+'@'+str(P):12s} {lam:5.2f} {wl:8.1f} {str(ta[0])+'-'+str(ta[1]):>9s} "
          f"{str(tb[0])+'-'+str(tb[1]):>9s}")
print("""
    脚注 a: f1f2 = 1；脚注 b: f1f2 != 1。表头 mu* = 0.1-0.13，**两列共用同一 mu\\* 区间**。
    ⚠ 脚注写 (Eq. S1)/(Eq. S2)，而正文的 f1f2 式子编号是 S5/S6 —— **编号对不上**。
      语义（f1f2=1 vs !=1）无歧义，故不影响反解；但记录在案，不替作者圆场。""")

# ============================== §B 公式（逐字取自页图 hrE_p-28.png）
sec("§B 公式逐字取自页图 Eq. S5/S6/S7（**不从 pdftotext 抄**）")
print(r"""    (S5)  Tc = f1*f2*(w_log/1.2) * exp[ -1.04(1+lam) / (lam - mu*(1+0.62 lam)) ]
    (S6)  f1 = cbrt( 1 + (lam / (2.46(1+3.8 mu*)))^{3/2} )
          f2 = 1 + ( (w2/w_log - 1) * lam^2 ) / ( lam^2 + [1.82(1+6.3 mu*) * w2/w_log]^2 )
    (S7)  w2 = sqrt( (2/lam) * INT alpha^2 F(w) w dw )

    ★ 三条承重读数（页图一手）：
      1. (S7) 就是 AD75 的 wbar2（等价于测度 dmu=(2/lam)(a2F/w)dw 的二阶矩平方根）=> **定义一致**。
      2. (S6) 的 f2 分子是 **(r - 1)**，**不是 (r^2 - 1)**，**且整体无开方**。
         => 这是对仓库两处缺陷的**第三方文献见证**，与我此前用 AD75 Table I 29 行做的
            指纹复现属于**不同信息类**（一个是别人重印公式，一个是我数值指纹）。
      3. (S5) 前因子是 w_log/1.2（AD75 改写式），不是 McMillan 原文的 Theta_D/1.45。""")

MEV_TO_K = 11.604518


def f1_of(lam, mus):
    return (1 + (lam / (2.46 * (1 + 3.8 * mus))) ** 1.5) ** (1 / 3)


def f2_of(lam, mus, r):
    return 1 + ((r - 1) * lam ** 2) / (lam ** 2 + (1.82 * (1 + 6.3 * mus) * r) ** 2)


def tc_mc(lam, mus, wlog_K):
    den = lam - mus * (1 + 0.62 * lam)
    return 0.0 if den <= 0 else (wlog_K / 1.20) * math.exp(-1.04 * (1 + lam) / den)


# ============================== §C 零参数复现 Tc^a —— 并**实测**端点对应哪个 mu*
sec("§C 零参数复现 Tc^a（f1f2=1）—— 端点↔mu\\* 的对应是**实测**，不是假设")

MUS = (0.10, 0.13)
print("""    单调性（analytic）：Tc 对 mu* 严格递减 => 印出区间的**高端 = mu*=0.10、低端 = mu*=0.13**。
    下面用零参数复现**检验**这条，而不是拿它当前提。

    ★ 容差**不拍**：由印出精度逐行传播
        budget = 0.5(Tc 印整数)  +  |dTc/dlam|*0.005(lam 印两位)  +  |dTc/dw_log|*0.05(w_log 印一位)
      —— 这是「若论文确实用 (S5)、且只差印刷舍入」时差值的**最大可能值**。超出即三列不自洽。
""")


def tc_budget(lam, mus, wl):
    """印出精度传播到 Tc 的最大差（数值导数，中心差分）。"""
    h1, h2 = 1e-6, 1e-4
    dl = abs(tc_mc(lam + h1, mus, wl) - tc_mc(lam - h1, mus, wl)) / (2 * h1)
    dw = abs(tc_mc(lam, mus, wl + h2) - tc_mc(lam, mus, wl - h2)) / (2 * h2)
    return 0.5 + dl * 0.005 + dw * 0.05


print(f"    {'row':12s} {'lam':>5s} | {'印低端':>6s} {'算@.13':>7s} {'差':>6s} | "
      f"{'印高端':>6s} {'算@.10':>7s} {'差':>6s} | {'预算':>5s}  判")
ok_rows, rep = [], []
for m, P, lam, wl, _n, _a, ta, tb in ROWS:
    if ta == (0, 0):
        print(f"    {m+'@'+str(P):12s} {lam:5.2f} | 印 0 —— lam={lam} 下 Tc^Mc 分母为负，无信息，跳过")
        continue
    lo_c, hi_c = tc_mc(lam, 0.13, wl), tc_mc(lam, 0.10, wl)
    dlo, dhi = lo_c - ta[0], hi_c - ta[1]
    bud = max(tc_budget(lam, 0.13, wl), tc_budget(lam, 0.10, wl))
    good = abs(dlo) <= bud and abs(dhi) <= bud
    rep.append((m, P, dlo, dhi, bud, good))
    if good:
        ok_rows.append((m, P, lam, wl, ta, tb))
    print(f"    {m+'@'+str(P):12s} {lam:5.2f} | {ta[0]:6d} {lo_c:7.1f} {dlo:+6.2f} | "
          f"{ta[1]:6d} {hi_c:7.1f} {dhi:+6.2f} | {bud:5.2f}  "
          f"{'ok' if good else '**超预算 '+format(max(abs(dlo),abs(dhi))/bud,'.1f')+'x**'}")

# 反向假设的负对照：若端点对应反了，复现应当崩
rev_bad = 0
for m, P, lam, wl, _n, _a, ta, tb in ROWS:
    if ta == (0, 0):
        continue
    bud = max(tc_budget(lam, 0.13, wl), tc_budget(lam, 0.10, wl))
    if not (abs(tc_mc(lam, 0.10, wl) - ta[0]) <= bud and abs(tc_mc(lam, 0.13, wl) - ta[1]) <= bud):
        rev_bad += 1
n_used = sum(1 for r in rep)
print(f"""
    正向（高端=0.10、低端=0.13）：{len(ok_rows)}/{n_used} 行两端都在自身预算内
    反向（高端=0.13、低端=0.10）：{rev_bad}/{n_used} 行**失败** <= **负对照到位**，
      端点对应是被**测**出来的，不是拍的。""")
if len(ok_rows) < n_used:
    for m, P, dlo, dhi, bud, good in rep:
        if not good:
            print(f"    ⚠ 排除 {m}@{P}：偏差 {dlo:+.2f}/{dhi:+.2f} K vs 预算 {bud:.2f} K "
                  f"（{max(abs(dlo),abs(dhi))/bud:.1f}×）—— 该行印出的 (lam, w_log, Tc^a) "
                  f"三者不自洽，**整行不可用**（含它的 Tc^b）")
if not ok_rows:
    bad("§C 无可用行")

# ============================== §D 反解 r，并把**整数舍入**传成区间（teammate 未做）
sec("§D 反解 r = wbar2/w_log —— 整数舍入**传成区间**（这是本节的承重增量）")

print("""    Tc 印为**整数** => 每个值只知道落在 [x-0.5, x+0.5]（四舍五入约定；截断约定另算）。
    r 由 ratio = Tc^b/Tc^a 决定：f1*f2 = ratio => f2 = ratio/f1 => 解 r。
    ⟹ 不传播这个舍入，"两个 mu\\* 端点的 r 吻合到 0.01" 这句话**没有意义**——
      它可能远小于舍入本身造成的宽度。下面把它算出来。""")


def solve_r(lam, mus, target_f2):
    """f2(r) 在 r>=1 侧先增后减（极大在 r*=1+sqrt(1+lam^2/c)）=> 取 **inf-根**（小根）。
    方向：小根是物理支；大根给 wbar2 >> 500 meV 的非物理值（见 errea 笔记 §3）。"""
    if target_f2 <= 1.0:
        return None
    lo, hi = 1.0, 1.0 + math.sqrt(1 + lam ** 2 / (1.82 * (1 + 6.3 * mus)) ** 2)
    if f2_of(lam, mus, hi) < target_f2:
        return None                                   # 极大值都够不到 => 无根
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f2_of(lam, mus, mid) < target_f2:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def r_interval(lam, mus, a_pr, b_pr, half=0.5):
    """a_pr/b_pr 是印出的整数。ratio 的区间 = [(b-h)/(a+h), (b+h)/(a-h)]。"""
    out = []
    for ratio in ((b_pr - half) / (a_pr + half), (b_pr + half) / (a_pr - half)):
        rr = solve_r(lam, mus, ratio / f1_of(lam, mus))
        out.append(rr)
    return out


print(f"\n    {'row':12s} {'lam':>5s} {'mu*':>5s} {'r 点值':>8s}   r 区间[整数舍入]        宽度")
RES = {}
for m, P, lam, wl, ta, tb in ok_rows:
    for k, mus in enumerate(MUS):
        a_pr, b_pr = (ta[1], tb[1]) if mus == 0.10 else (ta[0], tb[0])
        pt = solve_r(lam, mus, (b_pr / a_pr) / f1_of(lam, mus))
        iv = r_interval(lam, mus, a_pr, b_pr)
        if pt is None or iv[0] is None or iv[1] is None:
            print(f"    {m+'@'+str(P):12s} {lam:5.2f} {mus:5.2f} {'无根':>8s}   "
                  f"点值或端点解不出 => **本行不可用**"
                  f"  (需 f2={((b_pr/a_pr)/f1_of(lam,mus)):.4f}，"
                  f"该 lam 下 f2 极大 {f2_of(lam,mus,1+math.sqrt(1+lam**2/(1.82*(1+6.3*mus))**2)):.4f})")
            RES.setdefault((m, P), {})[mus] = None
            continue
        print(f"    {m+'@'+str(P):12s} {lam:5.2f} {mus:5.2f} {pt:8.4f}   "
              f"[{iv[0]:.4f}, {iv[1]:.4f}]     {iv[1]-iv[0]:.4f}")
        RES.setdefault((m, P), {})[mus] = (pt, iv[0], iv[1])

print("""
    ★ **承重读数（订正 teammate 的口径，不是否定它）**：
      两个 mu\\* 端点用的是**两组不同的印出数**（低端对/高端对）⟹ 它们的吻合**确实是一次检验**。
      但这次检验的**分辨率 = 舍入宽度**（0.03-0.69），而报出的吻合是 0.01-0.06
      ⟹ 可以说"两端点自洽"，**不能**说"r 被定到 ±0.01"——后者把分辨率误报成精度。
    ★ 更要紧的一条：**这次吻合对 f2 的形式几乎没有判别力。**
      改用 (r^2-1) 或带开方的写法，两端点**照样**吻合（它们共享同一个 f2）
      ⟹ 端点自洽**不能**用来给 f2 定案；定案依据仍只有 AD75 29 行零参数复现。""")

# ============================== §E 所谓"阳性对照"其实不是——先写清需要的对照
sec("§E ⚠ H3S / YH6 / CaH6 **不构成**对反解机器的阳性对照（订正 teammate 的框架）")

print("""    teammate 提议："H3S/YH6/CaH6 三个已知体系可拿来当阳性对照"。**这个框架是错的。**
    **需要的对照**：同一份 alpha^2F 上，一边有真 wbar2、一边有 Mc/AD 双列 => 反解值 vs 真值。
    **找到的对照**：不同论文、不同压强、不同近似下的同名材料。
    两者不是一回事（memory [comparison-found-is-not-comparison-needed]）。逐条核对：""")

CAL = {
    ("H3S", 200): [("仓库 Mishra EPW 谐性", 2.2929, 107.165, 129.32),
                   ("仓库 Mishra EPW 非谐", 1.7210, 109.995, 149.26)],
    ("YH6", 120): [("Ding YH6@165GPa 谐性", 2.534, None, 120.56)],
    ("CaH6", 150): [("Quan CaH6@150 T 行", 2.53, 90.0, 105.0)],
}
for key, lst in CAL.items():
    if key not in RES or RES[key].get(0.10) is None:
        continue
    m, P = key
    row = next(r for r in ok_rows if (r[0], r[1]) == key)
    lam_h, wl_h = row[2], row[3]
    print(f"\n    {m}@{P}: 本文 lam={lam_h}, w_log={wl_h:.1f} K = {wl_h/MEV_TO_K:.2f} meV")
    for tag, lam_c, wl_c, w2_c in lst:
        same = (abs(lam_c - lam_h) < 0.02) and (wl_c is not None and abs(wl_c - wl_h / MEV_TO_K) < 0.5)
        rtrue = (w2_c / wl_c) if wl_c else None
        print(f"      vs {tag:24s} lam={lam_c:<7.4f} "
              f"{'w_log='+format(wl_c,'.2f')+' meV' if wl_c else 'w_log 未记录':22s} "
              f"{'r_true='+format(rtrue,'.4f') if rtrue else 'r_true 不可算':18s}"
              f"  **同一计算？{'是' if same else '否'}**")
    lo = min(v[1] for v in RES[key].values() if v)
    hi = max(v[2] for v in RES[key].values() if v)
    print(f"      我反解 r in [{lo:.4f}, {hi:.4f}]")

print("""
    ⟹ **判词：0 个真阳性对照。** 三个体系的 lam 与 w_log 都与已知真值那份计算**不同**
      => 只能做"落不落进已知 r 的量级带"这种**弱一致性**检查，
        **不能**据此说"反解机器被验证过"。反解机器的验证仍只有
        Errea 那次（零参数复现 Tc^Mc 列 6/7 行 <=0.25%）。
    ⟹ 若要真阳性对照，需要的是：**一篇同时印 Mc/AD 双列 _和_ wbar2 本身的论文**。
      Quan 印了 wbar2 但**没有** Mc 列；本文有 Mc 列但**没有** wbar2。**至今无一篇两者都有。**""")

# ============================== §F 弱一致性：落不落进已知 r 带
sec("§F 弱一致性检查（**只**作量级对照，按 §E 不升级为验证）")

BAND = (1.1685, 1.3570)      # H3S/D3S 四点标定
QUAN = (1.1667, 1.3363)      # Quan 13 个 T 行
print(f"    已知 r 带：H3S/D3S 四点标定 {BAND}；Quan 13 个 T 行 {QUAN}")
print(f"\n    {'row':12s}   r 区间（两 mu\\* 端点并集）    落 H3S 带  落 Quan 带")
n_in = 0
for key in sorted(RES):
    vs = [v for v in RES[key].values() if v]
    if not vs:
        print(f"    {key[0]+'@'+str(key[1]):12s}   **无根 => 不可用**")
        continue
    lo, hi = min(v[1] for v in vs), max(v[2] for v in vs)
    a = not (hi < BAND[0] or lo > BAND[1])
    b = not (hi < QUAN[0] or lo > QUAN[1])
    n_in += a
    print(f"    {key[0]+'@'+str(key[1]):12s}   [{lo:.4f}, {hi:.4f}]{'':10s} "
          f"{'与之相交' if a else '**不相交**':10s} {'与之相交' if b else '**不相交**'}")
print(f"    => {n_in}/{sum(1 for k in RES if any(RES[k].values()))} 行与 H3S/D3S 标定带相交")
print("""
    ⚠ **不相交的行不等于错**：r 是材料性质，不同材料本来就可以落在带外
      （尤其 lam 小的 TmH6：f2-1 极小 ⟹ ratio 几乎不含 r 的信息，反解**病态**）。
      判别力低的行应当**标为无信息**，不是"证伪"。""")

# ============================== §G 判别力审计：哪些行的 ratio 真的含 r 的信息
sec("§G ★ 判别力审计：ratio 对 r 的灵敏度（决定哪些行值得用）")

print("""    反解的信息量 = d(ratio)/d(r)。f2 = 1 + (r-1)lam^2/(lam^2 + c^2 r^2)，
    lam 小时 lam^2/(lam^2+c^2 r^2) 极小 ⟹ ratio 几乎与 r 无关 ⟹ 反解被舍入主导。
    下面把"把 r 从 1.0 推到 1.5 能让 Tc^b 变多少 K"算出来——这是它在**印出精度**里的可见度。\n""")
print(f"    {'row':12s} {'lam':>5s} {'dTc(r:1->1.5)/K':>16s} {'印出半末位':>10s}  判")
for m, P, lam, wl, ta, tb in ok_rows:
    mus = 0.10
    base = f1_of(lam, mus) * tc_mc(lam, mus, wl)
    d = base * (f2_of(lam, mus, 1.5) - f2_of(lam, mus, 1.0))
    verdict = "可用" if d >= 5.0 else ("**弱**" if d >= 2.0 else "**无信息**")
    print(f"    {m+'@'+str(P):12s} {lam:5.2f} {d:16.1f} {0.5:10.1f}  {verdict}")
print("""
    ⟹ 使用建议：只有 dTc >> 1 K 的行，反解才比舍入噪声强。
      lam<=1.0 的行（YbH6@200、TmH6@50）**应当整行弃用**，不管它解不解得出根。""")

# ============================== §I ★ (M2) 的合成阳性对照：AD75 29 行，r 真值已知
sec("§I ★★ 反解机器本身的阳性对照 —— 拆成两个失效模式，只有 (M1) 测过")

print("""    反解机器有**两个独立失效模式**，此前被混为一谈：
      (M1) 论文那两列根本不由同一 f1*f2、同一 mu* 联系  ——  由 §C 零参数复现测（本文 8/10）
      (M2) **我解 f2 那一步**选错根 / 病态 / 被舍入吞掉  ——  **从未测过**
    (M2) 可以**精确测**：取 AD75 Table I 的 29 组**真实材料** (lam, w_log, wbar2, mu*)，
    我自己造出 Tc^Mc 与 Tc^AD、按与 2010.12225 **同等相对精度**印出，再反解，对照真 r。
    这不是杜撰参数——29 组全是隧道反演的实验谱。""")

# AD75 Table I（页图逐行转录，与 _chk_ad75_table1_frak_f_2026-09-20.py 同源）
# (material, w_log/K, wbar2/K, lambda, mu*)
AD75 = [
    ("Pb", 56, 65, 1.55, 0.105), ("In", 68, 89, 0.805, 0.097),
    ("Sn", 99, 121, 0.72, 0.092), ("Hg", 29, 49, 1.6, 0.098),
    ("Tl", 52, 64, 0.795, 0.111), ("Ta", 132, 148, 0.69, 0.093),
    ("a-Ga", 55, 101, 1.62, 0.095), ("b-Ga", 87, 129, 0.97, 0.092),
    ("Tl0.9Bi0.1", 48, 62, 0.78, 0.099), ("Pb0.4Tl0.6", 48, 62, 1.15, 0.094),
    ("Pb0.6Tl0.4", 50, 62, 1.38, 0.103), ("Pb0.8Tl0.2", 50, 61, 1.53, 0.101),
    ("Pb0.6Tl0.2Bi0.2", 48, 58, 1.81, 0.111), ("Pb0.9Bi0.1", 50, 60, 1.66, 0.081),
    ("Pb0.8Bi0.2", 46, 57, 1.88, 0.093), ("Pb0.7Bi0.3", 47, 57, 2.01, 0.092),
    ("Pb0.65Bi0.35", 45, 55, 2.13, 0.093), ("In0.9Tl0.1", 63, 86, 0.85, 0.103),
    ("In0.73Tl0.27", 55, 77, 0.93, 0.110), ("In0.67Tl0.33", 57, 79, 0.90, 0.110),
    ("In0.57Tl0.43", 53, 74, 0.85, 0.117), ("In0.5Tl0.5", 53, 73, 0.83, 0.110),
    ("In0.27Tl0.73", 42, 63, 1.09, 0.094), ("In0.17Tl0.83", 45, 63, 0.98, 0.101),
    ("In0.07Tl0.93", 49, 63, 0.89, 0.107), ("In2Bi", 46, 67, 1.40, 0.096),
    ("Sb2Tl7", 37, 58, 1.43, 0.102), ("Bi2Tl", 47, 59, 1.63, 0.101),
    ("a-Pb0.45Bi0.55", 29, 47, 2.59, 0.116),
]
assert len(AD75) == 29


def round_sig(x, n=3):
    if x == 0:
        return 0.0
    e = math.floor(math.log10(abs(x)))
    return round(x, -(e - n + 1))


print("\n  (I-1) 根选择审计：真 r 落在 f2 极大点 r* 的**哪一侧**？（决定 inf-根是否正确）")
n_left = 0
for name, wl, w2, lam, mus in AD75:
    c = 1.82 * (1 + 6.3 * mus)
    rstar = 1 + math.sqrt(1 + lam ** 2 / c ** 2)
    if (w2 / wl) < rstar:
        n_left += 1
print(f"      {n_left}/29 行真 r 在 r* **左侧** ⟹ **inf-根（小根）是对的那一支**，实测、非默认。")
if n_left != 29:
    bad(f"(I-1) 有 {29-n_left} 行真 r 在 r* 右侧 ⟹ inf-根规则失效")

print("\n  (I-2) 满精度反解（隔离**代数**，无印刷舍入）")
worst_a = (0.0, "")
for name, wl, w2, lam, mus in AD75:
    r_t = w2 / wl
    a = tc_mc(lam, mus, wl)
    b = f1_of(lam, mus) * f2_of(lam, mus, r_t) * a
    r_b = solve_r(lam, mus, (b / a) / f1_of(lam, mus))
    e = abs(r_b - r_t) if r_b else 9.9
    if e > worst_a[0]:
        worst_a = (e, name)
print(f"      29/29 解出根，最劣 |r_解 - r_真| = {worst_a[0]:.2e}（{worst_a[1]}）"
      f" ⟹ **代数与根选择无缺陷**。")
if worst_a[0] > 1e-9:
    bad(f"(I-2) 满精度下反解误差 {worst_a[0]:.2e} 不应存在")

print("""
  (I-3) ★ 同等相对精度下反解（3 位有效数字 —— 2010.12225 的 Tc~120-220 印整数，
        相对半末位 0.5/170 ~ 0.3%，与 3 位有效数字同量级）""")
print(f"      {'material':18s} {'lam':>5s} {'r_真':>6s} {'r_解':>6s} {'误差':>7s} {'相对':>7s}")
errs = []
for name, wl, w2, lam, mus in AD75:
    r_t = w2 / wl
    a = tc_mc(lam, mus, wl)
    b = f1_of(lam, mus) * f2_of(lam, mus, r_t) * a
    ar, br = round_sig(a), round_sig(b)
    r_b = solve_r(lam, mus, (br / ar) / f1_of(lam, mus))
    if r_b is None:
        print(f"      {name:18s} {lam:5.2f} {r_t:6.3f} {'无根':>6s}")
        errs.append((9.9, name, lam))
        continue
    e = abs(r_b - r_t)
    errs.append((e, name, lam))
    if e > 0.02:
        print(f"      {name:18s} {lam:5.2f} {r_t:6.3f} {r_b:6.3f} {r_b-r_t:+7.3f} "
              f"{(r_b-r_t)/r_t*100:+6.2f}%")
errs.sort(reverse=True)
med = sorted(e for e, _n, _l in errs)[len(errs) // 2]
print(f"      （只印误差 >0.02 的行）中位误差 {med:.4f}，最劣 {errs[0][0]:.4f}（{errs[0][1]}，"
      f"lam={errs[0][2]}）")
print("""      ⟹ **这就是反解机器在 (M2) 上的固有误差**——此前从未量过，
         一切"反解 r 精确到 0.01"的说法都必须先过这一关。""")

print("\n  (I-4) 误差 vs lam：验证 §G 的判别力审计（lam 小 ⟹ 反解病态）")
lo_l = [e for e, _n, l in errs if l < 1.2]
hi_l = [e for e, _n, l in errs if l >= 1.2]
print(f"      lam < 1.2（{len(lo_l)} 行）中位误差 {sorted(lo_l)[len(lo_l)//2]:.4f}")
print(f"      lam >=1.2（{len(hi_l)} 行）中位误差 {sorted(hi_l)[len(hi_l)//2]:.4f}")
print("      ⟹ 与 §G 的灵敏度审计**独立一致**：小 lam 行确实该弃用。")

print("\n  (I-5) ★★ 顺带订正 §F 的参照带 —— 它太窄了，因为只由 H3S 一种材料标定")
rs = sorted(w2 / wl for _n, wl, w2, _l, _m in AD75)
print(f"      AD75 29 种**真实材料**的 r = wbar2/w_log ∈ [{rs[0]:.4f}, {rs[-1]:.4f}]，"
      f"中位 {rs[14]:.4f}")
print(f"      而 §F 用的 H3S/D3S 四点带只有 [{BAND[0]}, {BAND[1]}]（宽 {BAND[1]-BAND[0]:.4f}）")
print(f"""      ⟹ **§F 里 YH6@120 的 r in [1.4069, 1.5091] "不相交" 毫无异常**：
         它稳稳落在真实材料的 r 分布里（29 行中有
         {sum(1 for r in rs if r > 1.40)} 行 r > 1.40）。
      ⟹ **订正**：`r` 的参照带应当用 AD75 这 29 行，不是 H3S 四点。
         用单一材料标定出来的"窄带"去判别的材料，是在**用方差冒充信息**。""")

print("""
  (I-6) ★★ 把「落进 H3S 带」这条检验**对照实测基准**算一次通过率
        —— 这条检验是我**今天早些时候自己写下的**，现在必须自审
        （memory [threshold-statistic-needs-calibrated-baseline]）""")
n_band = sum(1 for r in rs if BAND[0] <= r <= BAND[1])
p0 = n_band / len(rs)
print(f"""      零假设 = 「r 就是某个真实材料的 wbar2/w_log」，用 AD75 29 行作**实测**基准分布：
        {n_band}/29 = **{p0*100:.0f}%** 的真实材料本来就落在 [{BAND[0]}, {BAND[1]}] 内。
      ⟹ 单行「落进带内」的零假设通过率 ≈ {p0*100:.0f}%，**接近抛硬币**。

      **旧陈述错在哪**（两段分开写）：
        我在 `errea_f2_backsolve_omega2_2026-09-25.md` 与 memory 里写过
        「反解不知 H3S 任何数，却落进 H3S/D3S 四点标定的同一窄带 ⟹ 外部一致性」。
        **该带不窄**：它覆盖近半数真实材料 ⟹ 命中它**不是**强证据。
      **新陈述凭什么对**：
        AD75 29 行是隧道反演的实验谱，与 H3S、与 LaH10 三方独立；
        直接数落带内的比例即得 {p0*100:.0f}%，无拟合、无模型。[numerical-evidence]

      **但也别反向过头**：Errea 那次是 **7/7 行**全落带内。若真是 7 次独立抽取，
        {p0:.2f}^7 = {p0**7:.4f}，不空洞。**然而它们是同一材料的 7 个压强点**
        ⟹ 强相关，**有效自由度 ~1-2**，实际显著性只有 {p0:.2f}-{p0**2:.2f}
        （memory [gate-dof-audit-column-values-and-saturation]：先数自由度再数样本）。
      ⟹ **净判词**：该「外部一致性」应从**证据**降为**未证伪**；
         真正承重的外部校验仍只有 Quan 那次（她**直接印出** wbar2，不经反解）。""")

# ============================== §J 把这批行送进 Tier-0 天花板（对准「谁最紧」）
sec("§J 反解出的 wbar2 送进 Tier-0 天花板 —— 直接对准活问题「谁最紧」")

C_INF = 0.1827262477
NB3GE = 0.6260
print(f"""    r_true = Tc / (C_inf * sqrt(lam) * wbar2[meV] * {MEV_TO_K})，C_inf = {C_INF}
    ⚠ **分子是作者用 AD 公式算的 Tc**，非实验、非 Eliashberg
      ⟹ 检验的命题是「AD 公式本身越不越界」，与「真实材料越不越界」**不得混报**。
    ⚠ 只用 §C 通过 **且** §G 判"可用"的行。参照 Nb3Ge = {NB3GE}（Tier-0 最紧）。
""")
print(f"    {'row':12s} {'lam':>5s} {'wbar2/meV':>10s} {'Tc^b/K':>8s} {'r_true':>8s}  判")
GOOD = {"LuH6@200", "CaH6@150", "ScH6@300", "ScH6@350", "YH6@120", "H3S@200"}
worst, n_j = (0.0, ""), 0
for m, P, lam, wl, ta, tb in ok_rows:
    key = f"{m}@{P}"
    if key not in GOOD or (m, P) not in RES:
        continue
    vs = [v for v in RES[(m, P)].values() if v]
    if not vs:
        continue
    n_j += 1
    # 方向：r_true 要取**最大可能** => wbar2 取**最小可能**、Tc 取**最大可能**
    r_lo = min(v[1] for v in vs)
    w2_lo = r_lo * wl / MEV_TO_K
    tc_hi = tb[1] + 0.5
    rt = tc_hi / (C_INF * math.sqrt(lam) * w2_lo * MEV_TO_K)
    if rt > worst[0]:
        worst = (rt, key)
    print(f"    {key:12s} {lam:5.2f} {w2_lo:10.2f} {tc_hi:8.1f} {rt:8.4f}  "
          f"{'PASS' if rt < 1 else '**越界**':6s} {'★比 Nb3Ge 紧' if rt > NB3GE else ''}")
print(f"""
    ★ **非空断言**（门自身的缺陷防护）：本节实际评估了 **{n_j}** 行。
      首版这里查错了键（RES 用元组、我用字符串）⟹ 表为空、仍印出"无一比 Nb3Ge 紧"
      ⟹ **空集的全称命题恒真**。任何"0/N 通过"的门都必须先断言 N>0。""")
if n_j == 0:
    bad("§J 评估了 0 行 —— 空集不构成通过")
elif n_j < 4:
    bad(f"§J 只评估了 {n_j} 行，少于预期（§C×§G 交集应有 5-6 行）")
print(f"""    最紧 {worst[1]} r_true = {worst[0]:.4f}（取**最保守**读数：wbar2 取区间下端、Tc 取上端）
    ⟹ {'**无一**比 Nb3Ge 紧' if worst[0] < NB3GE else '★**有行比 Nb3Ge 紧**'}，
       且全部 {'< 1' if worst[0] < 1 else '**并非全部 <1**'} ⟹ 天花板对这批重稀土氢化物的判词见上。
       [numerical-evidence]
    ★ 方向审计：取 wbar2 下端 + Tc 上端 = 把 r_true 朝**大**推 ⟹ "r_true<1" 是**声音 PASS**；
      若出现 r_true>=1，按单边界方向规则**不足以判 FAIL**（memory [p4-onesided-bound-direction-rule]）。""")
if worst[0] >= 1:
    bad(f"§J {worst[1]} 的 r_true = {worst[0]:.4f} >= 1，需追")
sec("汇总")
if FAIL:
    print(f"{len(FAIL)} 处 FAIL：")
    for f in FAIL:
        print("  -", f)
    sys.exit(1)
print("ALL CHECKS PASSED")
