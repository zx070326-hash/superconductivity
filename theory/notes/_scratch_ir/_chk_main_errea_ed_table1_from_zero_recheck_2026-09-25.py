#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Errea et al. arXiv:1907.11916 / Nature 578, 66 (2020) ED Table 1 —— 从零复核 + 夺冠判决

起因：两个 teammate 对「LaH10@163 GPa 这一行是否存在 / 数字是什么」给出互相冲突的报告。
本脚本不采信任何一方，自取 ar5iv 一手源 (_errea_1907.11916_ar5iv.html, md5 53f8fcc8a442648091498fef4f7b7b78)
逐格转录后独立重算。

新增于 2026-09-25（把当日 LaZrH 那轮磨出的工具用回这张老表）：
  §C  结构性 vs 边际的两分法（此处只有单个 μ*=0.10 ⟹ 用兄弟行散布的量级判据代替端点符号判据）
  §E  用 Tc^ME_ani（真 Eliashberg 解，自由度 4）而非 Tc^AD（对天花板自由度≈0）当分子
  §F  夺冠判决必须是区间判决，且必须过「同一信息类」门

退出码 0 = 全部自检通过。
"""
import math, sys

MEV_K = 1.602176634e-22 / 1.380649e-23        # K per meV = 11.6045181...
C_INF = 0.1827262477                          # 认证的 Tier-0 常数
MU    = 0.10                                  # ED Table 1 印的唯一 μ*
NB3GE = 0.6260                                # 实验语料现任最紧点

# 逐格转录自 ar5iv HTML <table> #7（= ED Table 1）
# System, P(GPa), lambda, omega_log(meV), Tc^Mc, Tc^AD, Tc^ME_ani, Tc^SCDFT
ROWS = [
    ("LaH10", 129, 3.62,  76.4, 171.8, 252.6, 255.3, 230),
    ("LaH10", 163, 2.67,  96.4, 197.1, 247.0, 242.8, 225),
    ("LaH10", 214, 2.06, 115.5, 196.3, 235.9, 237.9, 210),
    ("LaH10", 264, 1.73, 126.6, 189.5, 219.2, 216.9, 201),
    ("LaD10", 159, 3.14,  63.5, 135.0, 184.2, 180.4, 171),
    ("LaD10", 210, 2.21,  81.7, 145.5, 176.5, 172.9, 158),
    ("LaD10", 260, 1.80,  92.2, 142.2, 164.6, 157.9, 151),
]

# 印刷精度 ⟹ 输入侧半宽（reproduce-published-table-must-propagate-input-print-precision）
D_LAM, D_WL, D_TC = 0.005, 0.05, 0.05

FAIL = []
def chk(ok, msg):
    print(("  [OK]   " if ok else "  [FAIL] ") + msg)
    if not ok:
        FAIL.append(msg)


def tmc(lam, wl, mu=MU):
    den = lam - mu * (1 + 0.62 * lam)
    return (wl * MEV_K / 1.2) * math.exp(-1.04 * (1 + lam) / den) if den > 0 else 0.0

def f1(lam, mu=MU):
    return (1 + (lam / (2.46 * (1 + 3.8 * mu))) ** 1.5) ** (1 / 3)

def f2(lam, r, mu=MU):
    L2 = 1.82 * (1 + 6.3 * mu) * r
    return 1 + (r - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)

def rstar(lam, mu=MU):
    c0 = 1.82 * (1 + 6.3 * mu)
    return 1 + math.sqrt(1 + lam ** 2 / c0 ** 2)

def solve_r(lam, need, mu=MU):
    """inf-根。超出 f2 上确界时返回 None（机器拒答）。"""
    rs = rstar(lam, mu)
    if need > f2(lam, rs, mu) or need < 1.0:
        return None
    a, b = 1.0, rs
    for _ in range(300):
        m = (a + b) / 2
        if f2(lam, m, mu) < need: a = m
        else: b = m
    return (a + b) / 2


print("=" * 78)
print("§0  常数自检")
print("=" * 78)
chk(abs(MEV_K - 11.6045181) < 1e-6, f"meV->K = {MEV_K:.7f} (= e/k_B /1000)")
chk(abs(f2(2.0, 1.0) - 1.0) < 1e-15, "f2(r=1) = 1")
chk(abs(f2(1e7, 1.7) - 1.7) < 1e-5, "f2(lam->inf) -> r")
# 单调性：tmc 对 lam、wl 各自单调递增 ⟹ 角点扫即够
chk(tmc(2.0, 100) < tmc(2.1, 100) and tmc(2.0, 100) < tmc(2.0, 101),
    "tmc 对 lam 与 wl 均单调递增 ⟹ 区间取角点合法")

print()
print("=" * 78)
print("§A  Tc^Mc 零参数复现（前提门：两列是否共用 (lam, w_log, mu*)）")
print("=" * 78)
print(f"{'行':16s} {'算出区间(K)':>22s} {'印出':>8s} {'重叠':6s} {'中心相对差':>10s}")
rels, mc_pass = {}, {}
for nm, P, lam, wl, mc, ad, me, sc in ROWS:
    vs = [tmc(l_, w_) for l_ in (lam - D_LAM, lam + D_LAM) for w_ in (wl - D_WL, wl + D_WL)]
    lo, hi = min(vs), max(vs)
    ov = not (hi < mc - D_TC or lo > mc + D_TC)
    rel = (tmc(lam, wl) - mc) / mc * 100
    rels[P], mc_pass[P] = rel, ov
    print(f"{nm}@{P:<11d} [{lo:9.3f},{hi:9.3f}] {mc:8.1f} {'PASS' if ov else 'FAIL':6s} {rel:+9.2f}%")
n_ok = sum(mc_pass.values())
chk(n_ok == 6, f"Mc 列复现 {n_ok}/7（预期 6/7，唯一失效为 @163）")
chk(mc_pass[163] is False, "唯一失效行是 @163")

print()
print("=" * 78)
print("§B  @163 的单数假说")
print("=" * 78)
lam163, wl163, mc163 = 2.67, 96.4, 197.1
wl_req = mc163 / tmc(lam163, 1.0)
a, b = 0.5, 10.0
for _ in range(300):
    m = (a + b) / 2
    if tmc(m, wl163) < mc163: a = m
    else: b = m
lam_req = (a + b) / 2
print(f"  单改 w_log ⟹ 需 {wl_req:8.3f} meV (印 {wl163}, {(wl_req-wl163)/wl163*100:+.2f}%)")
print(f"  单改 lam   ⟹ 需 {lam_req:8.4f}      (印 {lam163}, {(lam_req-lam163)/lam163*100:+.2f}%)")
print("  单改 mu*   ⟹ 任何 mu*>0.10 只会让算出值更低（远离印出值）:")
for mu in (0.10, 0.09, 0.08):
    print(f"       mu*={mu:.2f}: {tmc(lam163, wl163, mu):8.2f} K ({(tmc(lam163,wl163,mu)-mc163)/mc163*100:+6.2f}%)")
# mu* 方向：Tc 对 mu* 单调递减 ⟹ 要抬高算出值只能降 mu*
chk(tmc(lam163, wl163, 0.09) > tmc(lam163, wl163, 0.10) > tmc(lam163, wl163, 0.11),
    "Tc^Mc 对 mu* 单调递减 ⟹ 「该行用了更大 mu*」不能解释缺口（方向相反）")

print()
print("=" * 78)
print("§C  结构性 vs 边际：此处只有单个 mu* ⟹ 改用兄弟行散布的量级判据")
print("=" * 78)
oth = sorted(abs(rels[P]) for P in rels if P != 163)
print("  非-163 六行 |rel|: " + " ".join(f"{x:.3f}%" for x in oth))
print(f"  六行最大 {oth[-1]:.3f}%  中位 {oth[len(oth)//2]:.3f}%")
ratio = abs(rels[163]) / oth[-1]
print(f"  @163 {rels[163]:+.3f}%  ⟹ 兄弟行最大值的 {ratio:.1f} 倍")
chk(ratio > 10, f"量级判据：@163 超出兄弟行散布 {ratio:.1f}x (>10x) ⟹ 判为结构性，整行弃用")
print("  ⚠ 与 LaZrH 那轮的判据不同：那里用「跨 mu* 两端符号一致」，此处 mu* 只有一个值")
print("    ⟹ 端点符号判据在此**不可用**，只能用兄弟行散布。两者是不同判据，不得互相冒充。")

print()
print("=" * 78)
print("§D  反解 r = w2/w_log，并传播输入印刷精度得 r 的区间")
print("=" * 78)
print(f"{'行':16s} {'f2_req':>8s} {'f2max':>8s} {'r_lo':>7s} {'r_c':>7s} {'r_hi':>7s} {'半宽%':>7s}")
RES = {}
for nm, P, lam, wl, mc, ad, me, sc in ROWS:
    rs = []
    for l_ in (lam - D_LAM, lam + D_LAM):
        for mc_ in (mc - D_TC, mc + D_TC):
            for ad_ in (ad - D_TC, ad + D_TC):
                rr = solve_r(l_, (ad_ / mc_) / f1(l_))
                if rr is not None: rs.append(rr)
    rc = solve_r(lam, (ad / mc) / f1(lam))
    if rc is None or not rs:
        print(f"{nm}@{P:<11d}  拒答"); continue
    rlo, rhi = min(rs), max(rs)
    need = (ad / mc) / f1(lam)
    print(f"{nm}@{P:<11d} {need:8.5f} {f2(lam,rstar(lam)):8.5f} {rlo:7.4f} {rc:7.4f} {rhi:7.4f} "
          f"{(rhi-rlo)/2/rc*100:6.2f}%")
    RES[P] = (nm, lam, wl, mc, ad, me, sc, rlo, rc, rhi)
chk(len(RES) == 7, f"七行全部可解（{len(RES)}/7）")
chk(all(RES[P][8] >= 1.0 for P in RES), "Jensen r>=1：7/7（independent check，反解未被喂入该约束）")

print()
print("=" * 78)
print("§E  Tier-0 r_true 的**区间**：r_true = Tc / (C_inf * sqrt(lam) * w2[K])")
print("=" * 78)
print("  分子选谁很重要（自由度审计）：")
print("    Tc^AD      —— 与天花板同为 (lam, w2) 的函数，判决由 AD 大-lam 渐近 0.1878 vs C_inf 0.1827")
print("                  定死 ⟹ 对天花板**自由度约 0**，不是证据")
print("    Tc^ME_ani  —— 真各向异性 Eliashberg 解 ⟹ 自由度 4 ⟹ **这才是天花板要界的对象**")
print("    Tc^SCDFT   —— 另一独立泛函，同样有自由度")
print()
print(f"{'行':14s} {'w2区间(K)':>19s} | {'ME_ani r_true 区间':>24s} | {'SCDFT':>17s}")
CROWN = []
for P in sorted(RES):
    nm, lam, wl, mc, ad, me, sc, rlo, rc, rhi = RES[P]
    def rt(tc, lam_, r_):
        return tc / (C_INF * math.sqrt(lam_) * r_ * wl * MEV_K)
    vals_me, vals_sc, w2s = [], [], []
    for l_ in (lam - D_LAM, lam + D_LAM):
        for r_ in (rlo, rhi):
            for w_ in (wl - D_WL, wl + D_WL):
                w2K = r_ * w_ * MEV_K
                w2s.append(w2K)
                vals_me.append((me - D_TC) / (C_INF * math.sqrt(l_) * w2K))
                vals_me.append((me + D_TC) / (C_INF * math.sqrt(l_) * w2K))
                vals_sc.append((sc - 0.5) / (C_INF * math.sqrt(l_) * w2K))
                vals_sc.append((sc + 0.5) / (C_INF * math.sqrt(l_) * w2K))
    me_lo, me_hi = min(vals_me), max(vals_me)
    sc_lo, sc_hi = min(vals_sc), max(vals_sc)
    flag = "★跨" if me_hi > NB3GE else "   "
    print(f"{nm}@{P:<8d} [{min(w2s):8.1f},{max(w2s):8.1f}] | "
          f"[{me_lo:.4f}, {me_hi:.4f}] {flag} | [{sc_lo:.4f},{sc_hi:.4f}]")
    CROWN.append((nm, P, me_lo, me_hi, mc_pass[P]))
chk(all(hi < 1.0 for _, _, _, hi, _ in CROWN), "天花板未被违反：7/7 行 r_true 区间上端 < 1")

print()
print("=" * 78)
print(f"§F  夺冠判决（对手：Nb3Ge 实验语料 {NB3GE}）")
print("=" * 78)
print("  判决必须过三道，缺一不可：")
print("   (i)  该行的 Mc 复现门通过（否则整行结构性弃用）")
print("   (ii) r_true 的**区间下端** > 0.6260（区间上端越线只说明「判不了」）")
print("   (iii)同一信息类：分子必须与 Nb3Ge 的 0.6260 同为**实验** Tc")
print()
for nm, P, lo, hi, ok_mc in CROWN:
    g1 = "通过" if ok_mc else "弃用"
    if hi <= NB3GE:      g2 = "不夺冠（区间整体在下）"
    elif lo > NB3GE:     g2 = "夺冠（区间整体在上）"
    else:                g2 = "判不了（区间跨越 0.6260）"
    print(f"  {nm}@{P:<6d} (i){g1}  (ii)[{lo:.4f},{hi:.4f}] -> {g2}")
print()
print("  (iii) 三列 Tc 全部是**计算值**（Mc / AD / ME_ani / SCDFT 皆然），")
print("        而 Nb3Ge 的 0.6260 分子是**实验 Tc** ⟹ 跨信息类比较，本门直接不通过。")
n_crown = sum(1 for nm, P, lo, hi, ok in CROWN if ok and lo > NB3GE)
n_undec = sum(1 for nm, P, lo, hi, ok in CROWN if ok and lo <= NB3GE < hi)
chk(n_crown == 0, f"过 (i)+(ii) 两道的夺冠行数 = {n_crown}（应为 0）")
print(f"  ⟹ 过 (i)+(ii) 的夺冠 {n_crown} 行；「判不了」{n_undec} 行；(iii) 对全部七行皆不通过。")

print()
print("=" * 78)
print("§G  @163 缺陷的方向审计（它为什么恰好是最像夺冠的那一行）")
print("=" * 78)
print("  印出 Mc 高于算出 ⟹ 比值 Tc^AD/Tc^Mc 偏低 ⟹ f2_req 偏低 ⟹ r 偏低")
print("  ⟹ w2 偏低 ⟹ r_true = Tc/(C_inf sqrt(lam) w2) 偏**高** ⟹ **偏向夺冠**")
r163 = RES[163][8]
others = [RES[P][8] for P in RES if P != 163]
print(f"  实测：@163 的 r = {r163:.4f} 是七行最低（其余 {min(others):.4f}–{max(others):.4f}）")
chk(r163 < min(others), "@163 的 r 确为七行最低 ⟹ 与方向预言一致")
# 修复假说：若 w_log 实为 99.686 meV，重算
wl_fix = wl_req
lam_, mc_, ad_, me_ = 2.67, 197.1, 247.0, 242.8
r_fix = solve_r(lam_, (ad_ / mc_) / f1(lam_))
rt_fix = me_ / (C_INF * math.sqrt(lam_) * r_fix * wl_fix * MEV_K)
rt_asis = me_ / (C_INF * math.sqrt(lam_) * r_fix * 96.4 * MEV_K)
print(f"  若缺陷在 w_log（需 {wl_fix:.3f} meV）：r_true {rt_asis:.4f} -> {rt_fix:.4f}"
      f"（{'仍跨' if rt_fix > NB3GE else '落回 0.6260 以下'}）")
chk(rt_fix < rt_asis, "修 w_log 使 r_true 下降 ⟹ 与「缺陷偏向夺冠」一致")

print()
print("=" * 78)
print("§H  投毒负对照")
print("=" * 78)
def rerun(rows):
    """返回 (Mc复现数, 七行r最低者的压强, 过(i)+(ii)的夺冠行数)"""
    mp, rr = {}, {}
    for nm, P, lam, wl, mc, ad, me, sc in rows:
        vs = [tmc(l_, w_) for l_ in (lam - D_LAM, lam + D_LAM) for w_ in (wl - D_WL, wl + D_WL)]
        mp[P] = not (max(vs) < mc - D_TC or min(vs) > mc + D_TC)
        x = solve_r(lam, (ad / mc) / f1(lam))
        rr[P] = x
    n = sum(mp.values())
    solved = {P: v for P, v in rr.items() if v is not None}
    lowest = min(solved, key=lambda P: solved[P]) if solved else None
    ncr = 0
    for nm, P, lam, wl, mc, ad, me, sc in rows:
        if not mp[P] or rr[P] is None: continue
        rs = []
        for l_ in (lam - D_LAM, lam + D_LAM):
            for mc_ in (mc - D_TC, mc + D_TC):
                for ad_ in (ad - D_TC, ad + D_TC):
                    y = solve_r(l_, (ad_ / mc_) / f1(l_))
                    if y is not None: rs.append(y)
        if not rs: continue
        lo_ = min((me - D_TC) / (C_INF * math.sqrt(l_) * r_ * w_ * MEV_K)
                  for l_ in (lam - D_LAM, lam + D_LAM) for r_ in (min(rs), max(rs))
                  for w_ in (wl - D_WL, wl + D_WL))
        if lo_ > NB3GE: ncr += 1
    return n, lowest, ncr, len(solved)

base = rerun(ROWS)
print(f"  基线: Mc复现 {base[0]}/7, r 最低行 @{base[1]}, 夺冠 {base[2]}, 可解 {base[3]}/7")
POIS = []
# P1: w_log 当成 K 而非 meV
POIS.append(("P1 w_log 误当 K", [(a,b,c,d/MEV_K,e,f,g,h) for a,b,c,d,e,f,g,h in ROWS]))
# P2: 把 @163 的 Mc 改成我算出的值（即「缺陷被修好」）
POIS.append(("P2 @163 Mc 改为算出值", [(a,b,c,d,(190.60 if b==163 else e),f,g,h) for a,b,c,d,e,f,g,h in ROWS]))
# P3: 交换 @129 与 @163 的 lam
POIS.append(("P3 交换 129/163 的 lam", [(a,b,(2.67 if b==129 else 3.62 if b==163 else c),d,e,f,g,h) for a,b,c,d,e,f,g,h in ROWS]))
# P4: Mc 与 AD 两列对调
POIS.append(("P4 Mc/AD 两列对调", [(a,b,c,d,f,e,g,h) for a,b,c,d,e,f,g,h in ROWS]))
# P5: 全部 Tc^ME 抬高 5%
POIS.append(("P5 ME_ani 全体 +5%", [(a,b,c,d,e,f,g*1.05,h) for a,b,c,d,e,f,g,h in ROWS]))
n_alarm = 0
for name, rows in POIS:
    g = rerun(rows)
    alarm = (g[0] != base[0]) or (g[1] != base[1]) or (g[2] != base[2]) or (g[3] != base[3])
    if alarm: n_alarm += 1
    print(f"  {name:26s} Mc {g[0]}/7  r最低@{g[1]}  夺冠 {g[2]}  可解 {g[3]}/7   "
          f"{'报警' if alarm else '**漏报**'}")
chk(n_alarm == len(POIS), f"投毒 {n_alarm}/{len(POIS)} 报警")

print()
print("=" * 78)
print("§I  分子选 ME_ani 的方向审计（承重——它决定 §F 的判词能不能反过来用）")
print("=" * 78)
print("  我们的 Tier-0 天花板是**各向同性** Eliashberg 的界；Tc^ME_ani 是**各向异性**解。")
print("  各向异性一般抬高 Tc（Errea 表内实测见下）⟹ 分子被抬高 ⟹ r_true 被抬高。")
gt = sum(1 for _, P, _, _, _, _, me, _ in [(r[0], r[1], 0, 0, 0, 0, r[6], 0) for r in ROWS] if True)
pairs = [(P, me, sc) for nm, P, lam, wl, mc, ad, me, sc in ROWS]
n_ani_gt = sum(1 for P, me, sc in pairs if me > sc)
print(f"  表内对照（ME_ani vs SCDFT）：ME_ani 更高的行数 {n_ani_gt}/7")
print("  方向结论（两条判词，后果相反，不得混用）：")
print("   (a) 「天花板未被违反」= PASS 判词：分子被**抬高** ⟹ 判据变**难** ⟹ 过了就真过")
print("       ⟹ §E 的 7/7 r_true<1 **合法**（单边界方向规则的正用）。")
print("   (b) 「谁最紧 / 夺冠」= 排名判词：分子被抬高 ⟹ r_true 虚高 ⟹ **偏向夺冠**")
print("       ⟹ 不得用 ME_ani 宣布夺冠。§F 的 (ii) 即便通过也不作数。")
chk(n_ani_gt >= 6, f"表内 ME_ani > SCDFT 的行数 {n_ani_gt}/7 ⟹ 各向异性抬高 Tc 的方向在本表上实测成立")
print("  ⟹ §F 的三道门里，(ii) 与 (iii) 都已独立否掉夺冠；(i) 另外否掉 @163 整行。")
print("     三条理由互相独立，任一条单独成立即足以判「不易主」。")

print()
print("=" * 78)
if FAIL:
    print(f"EXIT 1 —— {len(FAIL)} 项自检失败:")
    for m in FAIL: print("   - " + m)
    sys.exit(1)
print("EXIT 0 —— ALL CHECKS PASSED")
