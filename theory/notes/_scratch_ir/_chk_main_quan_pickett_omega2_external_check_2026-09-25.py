# -*- coding: utf-8 -*-
r"""Quan-Ghosh-Pickett arXiv:1906.02695 Tables 1-3：一手 omega_2 的外部校验。

目的（四件，互相独立）：
  (1) **钉定义**——不靠 SM。正文写 `Tc^AD(lam, w_log, w_2; mu*)` 且 Eq.(1)
      `lam = int(2/w)a2F dw -> N(0)I^2/(M w_2^2)`，称 w_2 为 "second moment"。
      再用**零参数单位级预测** kappa_j = C*M_j*w_{2,j}^2（C 由 SI 常数算出，不拟合）
      把「w_2 以 meV 印、M 以 amu、kappa 以 eV/A^2」这一整套约定锁死。
  (2) **外部校验我 2026-09-25 的 Errea f2 反解**：Quan 的 LaH10 **直接印** w_2 =>
      形状比 r = w_2/w_log 可与我反解出的 r 带对照（**不可**比绝对 omega_2：
      Errea 是 SSCHA 非谐、Quan 是谐性，压强也不同 => 只有形状量可比）。
  (3) **Tier-0 天花板首次落在「论文印出的 omega_2」的氢化物行**——不经反解、不经 Jensen。
  (4) **eta/M 一行式的第二处文献现身** + 它与 sqrt(lam)*w_2 的恒等性（自由度审计）。

⚠ **mu* 未印**：正文只说 mu* "varies only within the range 0.10-0.15"，从未写用了哪个值。
   => 零参数复现 Tc **不可能**；改为「多行共享 1 个 mu*」的反演。

运行：OMP_NUM_THREADS=1 python theory/notes/_scratch_ir/_chk_main_quan_pickett_omega2_external_check_2026-09-25.py
数据源：本地 theory/notes/_scratch_ir/quan_pickett_1906.02695_ar5iv.html
        （ar5iv.labs.arxiv.org/html/1906.02695，**未用 -k**，TLS 正常校验通过）
"""
import html
import io
import math
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HTMLF = "theory/notes/_scratch_ir/quan_pickett_1906.02695_ar5iv.html"
FAIL = []


def bad(msg):
    FAIL.append(msg)
    print("  x FAIL:", msg)


def sec(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


# ============================================================ §A 解析三张表
sec("§A 从 ar5iv HTML 解析 Table 1/2/3（机器解析，不手抄）")

raw = open(HTMLF, encoding="utf-8").read()
ZW = ["​", "⁢", "⁡", " ", " "]


def strip_tex(cell):
    s = re.sub(r"<[^>]+>", " ", cell)
    s = html.unescape(s)
    for ch in ZW:
        s = s.replace(ch, "")
    return re.sub(r"\s+", " ", s).strip()


def parse_table(tid):
    m = re.search(r'<table id="%s"(.*?)</table>' % re.escape(tid), raw, re.S)
    assert m, tid
    rows = []
    for rm in re.finditer(r"<tr\b(.*?)</tr>", m.group(1), re.S):
        rows.append([strip_tex(c)
                     for c in re.findall(r"<t[hd]\b.*?>(.*?)</t[hd]>", rm.group(1), re.S)])
    return rows


def numify(s):
    s = s.strip()
    if s in ("-", "–", "—", "…", ""):
        return None
    try:
        return float(s)
    except ValueError:
        return None


MATRE = re.compile(r"^([A-Z][a-z]?)H\s?(\d+)\s*[A-Z][a-z]?H_\{\d+\}$")


def harvest(rows):
    """-> [(mat, P, atom, w_log, w_1, w_2, N0, I2, eta, kappa, A, lam, Tc)]"""
    out, mat, P = [], None, None
    for r in rows:
        if not r:
            continue
        c0 = r[0]
        mm = MATRE.match(c0.replace(" ", "").replace("_{", "_{"))
        if mm is None:
            mm = MATRE.match(re.sub(r"\s+", "", c0))
        if mm:                                   # 新材料行：mat, P, atom, 10 个数
            mat = f"{mm.group(1)}H{mm.group(2)}"
            P = int(r[1]); atom = r[2]; vals = r[3:13]
        elif re.fullmatch(r"\d{2,4}", c0):       # 新压强行：P, atom, 10 个数
            P = int(c0); atom = r[1]; vals = r[2:12]
        elif re.fullmatch(r"[A-Z][a-z]?|T", c0):  # 续行：atom, 10 个数
            atom = c0; vals = r[1:11]
        else:
            continue                              # 表头/单位行
        if len(vals) != 10 or mat is None:
            continue
        out.append((mat, P, atom, *[numify(v) for v in vals]))
    return out


DATA = []
for tid in ("A2.T1.2", "A2.T2.2", "A2.T3.2"):
    rows = parse_table(tid)
    got = harvest(rows)
    print(f"  {tid}: {len(rows)} tr -> {len(got)} 数据行")
    DATA += got

print(f"\n共收 {len(DATA)} 行：")
hdr = (f"{'mat':6s} {'P':>4s} {'at':3s} {'w_log':>6s} {'w_1':>6s} {'w_2':>6s} "
       f"{'N0':>5s} {'I2':>5s} {'eta':>6s} {'kap':>6s} {'A':>5s} {'lam':>5s} {'Tc':>5s}")
print(hdr)


def fm(v, w=6, p=2):
    return (" " * (w - 2) + "--") if v is None else f"{v:{w}.{p}f}"


for d in DATA:
    print(f"{d[0]:6s} {d[1]:4d} {d[2]:3s} {fm(d[3],6,0)} {fm(d[4],6,0)} {fm(d[5],6,0)} "
          f"{fm(d[6],5)} {fm(d[7],5,0)} {fm(d[8],6,1)} {fm(d[9],6,1)} {fm(d[10],5,0)} "
          f"{fm(d[11],5)} {fm(d[12],5,0)}")

TROWS = {(d[0], d[1]): d for d in DATA if d[2] == "T"}
# 13 个 (材料,压强) 组 x 3 行（金属/H/T）：SH3 3 + CaH6 4 + MgH6 2 + LaH10 2 + YH10 2
assert len(DATA) == 39, len(DATA)
assert len(TROWS) == 13, len(TROWS)

# --- 对我**先于本脚本、从页面纯文本 dump 逐行手抄**的 T 行做机器比对
#     （解析 bug 与手抄 bug 不能同时藏住；未手抄到的行明确标「未比对」，不冒充覆盖）
HAND = {
    ("SH3", 250): (None, 147, 157, 1.82, 218),   # w_log 当时未落在 dump 窗口内
    ("SH3", 280): (137, 155, 167, 1.61, 199),
    ("CaH6", 200): (107, 120, 128, 1.95, 190),
    ("CaH6", 250): (117, 133, 142, 1.69, 180),
    ("CaH6", 300): (122, 141, 152, 1.59, 175),
    ("MgH6", 300): (124, 138, 149, 2.54, 279),
    ("MgH6", 400): (149, 164, 176, 2.00, 269),
    ("LaH10", 250): (92, 109, 122, 2.46, 206),
    ("LaH10", 300): (115, 133, 145, 1.80, 189),
    ("YH10", 300): (113, 136, 151, 2.42, 249),
    ("YH10", 400): (142, 165, 180, 1.81, 236),
}
print("\n手抄 vs 机器解析（T 行）：")
nchk = 0
for key, hv in HAND.items():
    d = TROWS[key]
    for nm, h, g in zip(("w_log", "w_1", "w_2", "lam", "Tc"),
                        hv, (d[3], d[4], d[5], d[11], d[12])):
        if h is None:
            continue
        nchk += 1
        if abs(h - g) > 1e-9:
            bad(f"{key} {nm}: 手抄 {h} vs 解析 {g}")
nohand = sorted(set(TROWS) - set(HAND))
print(f"  比对 {nchk} 个数字（{len(HAND)}/{len(TROWS)} 个 T 行）"
      f"{'，全部一致' if not FAIL else ''}")
print(f"  **未独立手抄比对**（只有机器解析一路证据）：{nohand}")

# ============================================================ §B 定义钉死门
sec("§B 定义钉死（每条标自由度；identity 只算转抄门，不算物理检验）")

print("""正文（一手，本地 HTML）：
  · "Tc itself can be obtained sufficiently accurately from the Allen-Dynes equation
     Tc^AD(lam, w_log, w_2; mu*) in terms of ... the logarithmic w_log and second w_2
     frequency moments. Specific expressions are provided in the SM."
  · Eq.(1): lam = int (2/w) a2F(w) dw -> N(0)I^2/(M w_2^2) ... "the second moment w_2"
  · "lam_j = eta_j/kappa_j",  "M_j w_{2,j}^2 == kappa_j",  "eta_j = N_j(0) I_j^2"
=> w_2 就是 AD 的 omega_2 槽位；**SM 不是钉定义的必要条件**。""")

print("\n[B1] 幂平均链 w_log <= w_1 <= w_2（dmu=(2/lam)a2F/w dw 下指数 0<1<2 的幂平均）")
n1 = nok1 = 0
for d in DATA:
    wl, w1, w2 = d[3], d[4], d[5]
    if None in (wl, w1, w2):
        continue
    n1 += 1
    if wl <= w1 <= w2:
        nok1 += 1
    else:
        print(f"    违反: {d[0]}@{d[1]} {d[2]}  {wl} {w1} {w2}")
print(f"    {nok1}/{n1} 行满足。自由度：**真**（随机排序下每行先验 1/6）；"
      f"但必要非充分——只排除「w_2 不是 2 阶矩」的一大类误读。")
if nok1 != n1:
    bad(f"B1 幂平均链有 {n1-nok1} 行违反")

print("\n[B2] ★ 零参数单位级预测 kappa_j = C*M_j*w_{2,j}^2，C 由 SI 常数算出（**不拟合**）")
AMU = 1.66053906660e-27
HBAR = 1.054571817e-34
EV = 1.602176634e-19
C_KAPPA = AMU * (1e-3 * EV / HBAR) ** 2 / EV * 1e-20
print(f"    C = m_u*(1e-3 e/hbar)^2/e*1e-20 = {C_KAPPA:.6e}  [eV/A^2 per (amu*meV^2)]")
MASS = {"H": 1.008, "S": 32.06, "Ca": 40.078, "Mg": 24.305, "La": 138.905, "Y": 88.906}
print("""    判据（**由印刷精度导出，不是拍的容差**）：把 kappa 的半末位区间反解成 w_2 区间，
    再问「印出的整数 w_2 与该区间相差几个末位单位」。同时分别记录它与
    「四舍五入」和「截断」两种印刷约定的相容性——约定本身是**实测项**，不是假设。""")
print(f"    {'row':16s} {'w2_pr':>5s} {'w2*_lo':>7s} {'w2*_hi':>7s} {'dist':>5s} round floor")
n2 = nok2 = n2r = n2f = 0
worst2 = 0.0
for d in DATA:
    mat, P, at, w2, kap = d[0], d[1], d[2], d[5], d[9]
    if at == "T" or at not in MASS or None in (w2, kap):
        continue
    n2 += 1
    dk = 0.05 if kap >= 10 else 0.005          # kappa 印到 0.1 / 0.01 的半末位
    lo = math.sqrt((kap - dk) / (C_KAPPA * MASS[at]))
    hi = math.sqrt((kap + dk) / (C_KAPPA * MASS[at]))
    dist = 0.0 if lo <= w2 <= hi else (w2 - hi if w2 > hi else lo - w2)
    worst2 = max(worst2, dist)
    R = not (hi < w2 - 0.5 or lo > w2 + 0.5)   # 与「四舍五入到 w2」相容
    F = not (hi < w2 or lo > w2 + 1.0)         # 与「截断到 w2」相容
    n2r += R
    n2f += F
    ok = dist <= 1.0                           # 末位 1 个单位
    nok2 += ok
    print(f"    {mat+'@'+str(P)+' '+at:16s} {w2:5.0f} {lo:7.2f} {hi:7.2f} {dist:5.2f} "
          f"{'Y' if R else '.':^5s} {'Y' if F else '.':^5s} {'' if ok else ' OUT'}")
print(f"    {nok2}/{n2} 行：印出 w_2 与 kappa 反解区间相差 <= 1 个末位单位；最大 {worst2:.2f} meV")
print(f"    印刷约定（实测）：与四舍五入相容 {n2r}/{n2}；与截断相容 {n2f}/{n2}"
      f" => **截断是更好的假设，但非普适**（{n2-max(n2r,n2f)} 行两者皆不容，kappa 高出约 1.4%）。")
print("""    自由度：**真，且是本节最强的一条**——零可调参数，同时钉死三件事：
      w_2 的单位是 meV、M 是**单个**原子质量(amu)、kappa 列就是 M*w_2^2。
    ★ 顺带钉死**重数约定**：H 行的 kappa 用的是**一个** H 的质量（LaH10 有 10 个 H，
      若 kappa 按化学式单位计会大 10 倍，实测不是）=> **重数全部装在 eta_j = N_j(0)I_j^2
      的 N_j(0) 里，kappa 是逐原子量**。这与 A15 的 eta_A 重数问题同形，
      见 memory [klein1979-a15-multiplicity-convention]。""")
print(f"""    ⚠ **发现（源表内不一致，非我方门写错）**：{n2-nok2} 行的 kappa 隐含 w_2 比印出值高
      1.0-1.1 meV，**两种印刷约定都解释不了**。这是 Quan 表内 kappa 与 w_2 两列
      的一处 ~1.4%(kappa) 不一致，无法从印出精度消化。
      => 本条**不判 FAIL**：门的职责是检验「我的判决是否稳」，不是「源表是否自洽」。
         该不一致是否承重，交 §H 的敏感度门**算**出来，不靠拍阈值。""")
MAX_DIST_W2 = worst2

print("\n[B3] lam_j = eta_j/kappa_j —— Quan 自己的定义式 => **自由度 0**，只是转抄门")
n3 = nok3 = 0
for d in DATA:
    at, eta, kap, lam = d[2], d[8], d[9], d[11]
    if at == "T" or None in (eta, kap, lam):
        continue
    n3 += 1
    pred = eta / kap
    tol = 0.05 / kap + eta * 0.05 / kap ** 2 + 0.005
    if abs(pred - lam) <= tol + 1e-12:
        nok3 += 1
    else:
        print(f"    偏离: {d[0]}@{d[1]} {at}  eta/kap={pred:.4f} vs lam={lam} (tol {tol:.4f})")
print(f"    {nok3}/{n3} 行一致（含印刷舍入传播）")

print("\n[B4] eta_j = N_j(0)*I_j^2 —— 亦为定义式 => **自由度 0**（转抄门）")
print("     ⟹ 它的用途不是验证物理，而是**抓印刷错误**。下面果然抓到一处。")
n4b = nok4b = 0
b4bad = []
for d in DATA:
    at, N0, I2, eta = d[2], d[6], d[7], d[8]
    if at == "T" or None in (N0, I2, eta):
        continue
    n4b += 1
    lo, hi = (N0 - 0.005) * (I2 - 0.5), (N0 + 0.005) * (I2 + 0.5)
    ok = lo - 0.05 <= eta <= hi + 0.05
    nok4b += ok
    if not ok:
        b4bad.append((d, lo, hi))
        print(f"    偏离: {d[0]}@{d[1]} {at}  N0*I2 in [{lo:.2f},{hi:.2f}] vs eta={eta}")
print(f"    {nok4b}/{n4b} 行一致（N0 只印 2 位小数 => 区间被它主导）")

print("\n[B4b] ★ 上面那一行是 Quan Table 2 的**印刷错误**，两条互不相干的诊断都指向它")
print("     诊断①  eta = N0*I2:  0.13*51 = 6.63  vs  印出 eta = 8.8 ；若 N0=0.17 则 8.67 ≈ 8.8")
print(f"     诊断②  N_T = sum_j N_j(0)，在其余 12 组上成立到 ±0.01，只有这一组差 -0.04：")
print(f"       {'group':14s} {'N_X':>5s} {'N_H':>5s} {'sum':>5s} {'N_T':>5s} {'diff':>6s}")
n4c = nok4c = 0
for (mat, P), t in sorted(TROWS.items()):
    ps = [x for x in DATA if x[0] == mat and x[1] == P and x[2] != "T"]
    if len(ps) != 2 or t[6] is None:
        continue
    a, b = ps[0][6], ps[1][6]
    n4c += 1
    diff = a + b - t[6]
    ok = abs(diff) <= 0.015
    nok4c += ok
    print(f"       {mat+'@'+str(P):14s} {a:5.2f} {b:5.2f} {a+b:5.2f} {t[6]:5.2f} {diff:+6.2f}"
          f"{'   <== 唯一超出 ±0.01' if not ok else ''}")
print(f"     {nok4c}/{n4c} 组满足 N_T = sum N_j(0) 到 ±0.01")
print("""     => **判决**：CaH6@300 的 H 行 N↑(0) 印作 0.13，应为 ≈0.17。
        这是从零发现的**已发表表格的一处缺陷**（两条独立诊断，方向一致地指向 0.17）。
     ⚠ **消费点审计**：N↑(0) **不进**我的天花板（天花板只吃 lam、w_2；AD 公式只吃
        lam、w_log、w_2、mu*）=> 该缺陷**不承重**，发现但不消费，也不据此质疑同表的
        lam / w_2 / Tc 三列（那是三个独立命题，各自要证）。""")
assert len(b4bad) == 1 and (b4bad[0][0][0], b4bad[0][0][1], b4bad[0][0][2]) == ("CaH6", 300, "H"), \
    "B4 的偏离行变了——重新诊断，不要直接放宽"

print("\n[B5] lam_T = sum_j lam_j —— 原子分解的**内容**在这一条（自由度：1 条真约束/行）")
n5 = nok5 = 0
for (mat, P), d in sorted(TROWS.items()):
    parts = [x[11] for x in DATA if x[0] == mat and x[1] == P and x[2] != "T"]
    if not parts:
        continue
    n5 += 1
    s = sum(parts)
    ok = abs(s - d[11]) <= 0.005 * (len(parts) + 1) + 1e-12
    nok5 += ok
    print(f"    {mat}@{P}: sum lam_j={s:.3f} vs lam_T={d[11]:.2f}  {'ok' if ok else 'OUT'}")
print(f"    {nok5}/{n5}")
if nok5 < n5:
    bad(f"B5 有 {n5-nok5} 行 lam_T != sum lam_j")

print("\n[B6] ★ A = lam*w_1/2 —— **恒等式** => `A` 列**零新增信息**")
print("""    判据同 B2：由 lam(±0.005)、w_1(±0.5) 的印刷区间推出 A 的区间，
    再问印出的整数距该区间几个末位单位；并分别记录两种约定的相容性。""")
print(f"    {'row':16s} {'A_lo':>6s} {'A_hi':>6s} {'A_pr':>5s} {'dist':>5s} round floor")
n6 = nok6 = n6r = n6f = 0
worst6 = 0.0
for d in DATA:
    lam, w1, A = d[11], d[4], d[10]
    if None in (lam, w1, A) or A == 0:
        continue
    n6 += 1
    lo, hi = (lam - 0.005) * (w1 - 0.5) / 2, (lam + 0.005) * (w1 + 0.5) / 2
    dist = 0.0 if lo <= A <= hi else (A - hi if A > hi else lo - A)
    worst6 = max(worst6, dist)
    R = not (hi < A - 0.5 or lo > A + 0.5)
    F = not (hi < A or lo > A + 1.0)
    n6r += R
    n6f += F
    ok = dist <= 1.0
    nok6 += ok
    if dist > 0.5 or not ok:
        print(f"    {d[0]+'@'+str(d[1])+' '+d[2]:16s} {lo:6.2f} {hi:6.2f} {A:5.0f} {dist:5.2f} "
              f"{'Y' if R else '.':^5s} {'Y' if F else '.':^5s}{'' if ok else ' OUT'}")
print(f"    {nok6}/{n6} 行：印出 A 距区间 <= 1 个末位单位；最大 {worst6:.2f} meV")
print(f"    印刷约定（实测）：与四舍五入相容 {n6r}/{n6}；与截断相容 {n6f}/{n6}"
      f" => **A 列也支持「截断」**（与 B2 同向 ⟹ 是全表统一的印刷约定，不是随机误差）。")
print("""    => **自由度审计结论（承重）**：A := int a2F dw = lam*omega_1/2 是 AD75 定义的直接推论。
    故 Quan 的头条 "maximize the product lam*w_1"（以及 "Tc = S*A + T_o, S=0.150"）
    **不是**独立于 lam、w_1 两列的新信息，而是同一组数的重参数化。
    同形：memory [discriminating-ratio-invariant-under-reparametrization-2026-09-20]。""")
if nok6 < n6:
    bad(f"B6 有 {n6-nok6} 行超出「末位 1 个单位」")

# ============================================================ §C mu* 反演
sec("§C mu* 未印 => 反演（多行共享 1 参数；顺带作 f2 三写法的负对照）")


def Lam1(mus):
    return 2.46 * (1 + 3.8 * mus)


def Lam2(mus, r):
    return 1.82 * (1 + 6.3 * mus) * r


def f1_of(lam, mus):
    return (1 + (lam / Lam1(mus)) ** 1.5) ** (1 / 3)


def f2_of(lam, mus, r, form):
    q = lam ** 2 / (lam ** 2 + Lam2(mus, r) ** 2)
    if form == "r-1":
        return 1 + (r - 1) * q
    if form == "r2-1":
        return 1 + (r ** 2 - 1) * q
    if form == "sqrt":
        return 1 + math.sqrt(abs((r - 1) * q))
    raise ValueError(form)


MEV_TO_K = 11.604518


def tc_mc(lam, mus, wlog_K):
    den = lam - mus * (1 + 0.62 * lam)
    return 0.0 if den <= 0 else (wlog_K / 1.20) * math.exp(-1.04 * (1 + lam) / den)


def tc_ad(lam, mus, wlog, w2, form="r-1"):
    return (f1_of(lam, mus) * f2_of(lam, mus, w2 / wlog, form)
            * tc_mc(lam, mus, wlog * MEV_TO_K))


def solve_mus(lam, wlog, w2, tgt, form="r-1"):
    lo, hi = 1e-9, 0.60
    if tc_ad(lam, lo, wlog, w2, form) < tgt or tc_ad(lam, hi, wlog, w2, form) > tgt:
        return None
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if tc_ad(lam, mid, wlog, w2, form) > tgt:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


ROWS_C = [d for d in DATA if d[12] and d[12] > 0 and None not in (d[3], d[4], d[5])]
print(f"可用行（Tc>0 且三个频率齐全）：{len(ROWS_C)}\n")
res = {}
for form in ("r-1", "r2-1", "sqrt"):
    print(f"--- f2 写法 = {form}")
    print(f"    {'row':18s} {'lam':>5s} {'w_log':>6s} {'w_2':>5s} {'r':>7s} {'Tc':>5s} {'mu*':>8s}")
    vals = []
    for d in ROWS_C:
        lam, wl, w2, tc = d[11], d[3], d[5], d[12]
        mus = solve_mus(lam, wl, w2, tc, form)
        vals.append(mus)
        print(f"    {d[0]+'@'+str(d[1])+' '+d[2]:18s} {lam:5.2f} {wl:6.0f} {w2:5.0f} "
              f"{w2/wl:7.4f} {tc:5.0f} " + (f"{mus:8.4f}" if mus is not None else "     n/a"))
    good = [v for v in vals if v is not None]
    if good:
        lo, hi = min(good), max(good)
        med = sorted(good)[len(good) // 2]
        inr = sum(1 for v in good if 0.10 <= v <= 0.15)
        print(f"    => n={len(good)}/{len(vals)} 可解；mu* in [{lo:.4f}, {hi:.4f}]，"
              f"极差 {hi-lo:.4f}，中位 {med:.4f}；落在论文自述 [0.10,0.15] 的 {inr}/{len(good)}")
        res[form] = (len(good), lo, hi, hi - lo, inr, len(vals))
    print()

print("三写法对照（判别器 = 单一 mu* 能否同时解释所有行 的紧致度 + 落区间比例）：")
print(f"    {'form':8s} {'可解/总':>9s} {'mu*_lo':>8s} {'mu*_hi':>8s} {'极差':>7s} {'在[.10,.15]':>11s}")
for form, (n, lo, hi, sp, ir, tot) in res.items():
    print(f"    {form:8s} {str(n)+'/'+str(tot):>9s} {lo:8.4f} {hi:8.4f} {sp:7.4f} {ir:11d}")
print("""
    ⚠ 方向与限度：本节**不能**升级为「f2 形式由 Quan 定案」——mu* 未印使它成为
      1 参数拟合，三写法都能靠各自的 mu* 吸收一部分差异。f2 的定案依据仍是
      AD75 Table I 的 **29 行零参数**复现（3.85% / 4.84% / 13.15%），见
      ../ceiling/errea_f2_backsolve_omega2_2026-09-25.md §1。本节只是一致性旁证。""")

# ============================================================ §D 外部校验反解
sec("§D ★ 对我 Errea f2 反解的**外部**校验（可比量只有形状比 r）")

ERREA_BS = {
    ("LaH10", 129): (76.4, 100.96, 100.76),
    ("LaH10", 163): (96.4, 108.95, 120.88),
    ("LaH10", 214): (115.5, 142.51, 141.22),
    ("LaH10", 264): (126.6, 156.58, 155.35),
    ("LaD10", 159): (63.5, 79.23, 78.88),
    ("LaD10", 210): (81.7, 97.35, 97.39),
    ("LaD10", 260): (92.2, 109.10, 107.61),
}
print("我的反解（Errea 1907.11916 ED Table I，SSCHA 非谐，mu*=0.1）：")
print(f"    {'row':14s} {'w_log':>7s} {'r(读法一)':>11s} {'r(读法二)':>11s}")
allr, la_r = [], []
for (m, P), (wl, a, b) in sorted(ERREA_BS.items()):
    r1, r2 = a / wl, b / wl
    print(f"    {m+'@'+str(P):14s} {wl:7.1f} {r1:11.4f} {r2:11.4f}")
    allr += [r1, r2]
    if m == "LaH10":
        la_r.append((P, r1, r2))
print(f"    => 全七行 r in [{min(allr):.4f}, {max(allr):.4f}]")
la_ex163 = [v for P, r1, r2 in la_r if P != 163 for v in (r1, r2)]
print(f"    => LaH10 四行中**剔掉 @163** 后 r in [{min(la_ex163):.4f}, {max(la_ex163):.4f}]")
r163 = [v for P, r1, r2 in la_r if P == 163 for v in (r1, r2)]
print(f"    => @163 单独: r = {r163[0]:.4f} / {r163[1]:.4f}")

print("\nQuan（1906.02695 Table 3，谐性，**直接印出** w_2）——LaH10 的 T 行：")
print(f"    {'row':14s} {'w_log':>7s} {'w_2':>6s} {'r':>8s}   r 区间[四舍五入]      r 区间[截断]")
quan_r = []
for key in (("LaH10", 250), ("LaH10", 300)):
    d = TROWS[key]
    wl, w2 = d[3], d[5]
    r = w2 / wl
    rR = ((w2 - 0.5) / (wl + 0.5), (w2 + 0.5) / (wl - 0.5))   # 四舍五入约定
    rF = (w2 / (wl + 1.0), (w2 + 1.0) / wl)                   # 截断约定（B2/B6 实测支持）
    quan_r.append((key, r, min(rR[0], rF[0]), max(rR[1], rF[1])))
    print(f"    {key[0]+'@'+str(key[1]):14s} {wl:7.0f} {w2:6.0f} {r:8.4f}   "
          f"[{rR[0]:.4f}, {rR[1]:.4f}]   [{rF[0]:.4f}, {rF[1]:.4f}]")
qlo = min(x[2] for x in quan_r)
qhi = max(x[3] for x in quan_r)
print(f"    => 两种约定的并集 r in [{qlo:.4f}, {qhi:.4f}]"
      f"  （约定之差只值 ~0.002 于 r ⟹ **不影响本节判读**）")

print("\n全 11 个 T 行的 r（供看整体分布）：")
for (mat, P), d in sorted(TROWS.items()):
    print(f"    {mat+'@'+str(P):14s} r = {d[5]/d[3]:.4f}")
allT_r = [d[5] / d[3] for d in TROWS.values()]
print(f"    => 13 行 r in [{min(allT_r):.4f}, {max(allT_r):.4f}]")

print(f"""
判读（三条，方向都审过）：
  1. **量级与趋势相符。** 我的 LaH10 反解 r（剔掉 @163）= {min(la_ex163):.4f}-{max(la_ex163):.4f}；
     Quan LaH10 = {quan_r[0][1]:.4f}(250 GPa)、{quan_r[1][1]:.4f}(300 GPa)。
     两组都落在同一窄带，且都随压强下降。这是一次**真外部校验**：
     Quan 的 w_2 从未进入我的反演。
  2. ★ **@163 那一行被独立第二次指认为异常。** 它是我七行里 r 最低的一行（{r163[0]:.4f}），
     明显低于 Quan 的 LaH10 区间下端 {qlo:.4f}；而它**同时**是 Tc^Mc 零参数复现
     差 -3.30% 的那一行（其余 6 行 <=0.25%）。两个**互不相干**的诊断指向同一行
     => 「@163 的印出值有问题」从单一诊断升为双诊断。
  3. ⚠ **不可比绝对 omega_2**：Errea 是 SSCHA 非谐 + 不同压强 + mu* 不同；Quan 是谐性。
     绝对值对照无意义，本节只用形状比。也**不得**据本节把反解值升级为 certificate。""")
for key, r, rlo, rhi in quan_r:
    inside = min(allr) <= r <= max(allr)
    extra = "" if inside else f"（高出我的全带上端 {r - max(allr):+.4f}）"
    print(f"    {key[0]}@{key[1]}: r={r:.4f} {'落在' if inside else '**未落在**'}"
          f" 我的反解全带 [{min(allr):.4f}, {max(allr):.4f}] 内{extra}")

# ============================================================ §E Tier-0 天花板
sec("§E ★ Tier-0 天花板首次落在「论文印出的 omega_2」的氢化物行（无反解、无 Jensen）")

C_INF = 0.1827262477
NB3GE = 0.6260
print(f"Tc_ceiling[K] = C_inf*sqrt(lam)*w_2[meV]*{MEV_TO_K}，C_inf = {C_INF}")
print(f"参照：32 元素 + A15 语料最紧点 Nb3Ge r_true = {NB3GE}\n")
print(f"    {'row':18s} {'lam':>5s} {'w_2':>5s} {'Tc^AD':>6s} {'ceiling':>8s} {'r_true':>7s} 判决")
rt = []
for d in DATA:
    lam, w2, tc = d[11], d[5], d[12]
    if None in (lam, w2) or not tc:
        continue
    ceil = C_INF * math.sqrt(lam) * w2 * MEV_TO_K
    r_true = tc / ceil
    rt.append((f"{d[0]}@{d[1]} {d[2]}", r_true))
    if r_true >= 1:
        bad(f"天花板被越过: {d[0]}@{d[1]} {d[2]} r_true={r_true:.4f}")
    print(f"    {d[0]+'@'+str(d[1])+' '+d[2]:18s} {lam:5.2f} {w2:5.0f} {tc:6.0f} "
          f"{ceil:8.1f} {r_true:7.4f} {'PASS' if r_true<1 else '越界'}"
          f"{'  (> Nb3Ge)' if r_true > NB3GE else ''}")
rt.sort(key=lambda x: -x[1])
print(f"\n    {len(rt)}/{len(rt)} 行 PASS。最紧三行：")
for nm, v in rt[:3]:
    print(f"      {nm:18s} r_true = {v:.4f}")
print(f"    比 Nb3Ge 0.6260 更紧的行数 = {sum(1 for _, v in rt if v > NB3GE)}/{len(rt)}")
print("""
    ⚠ **分子是什么，必须写清**：这里的 Tc 是 **Quan 用 AD 公式算出的 Tc**，
      不是实验 Tc、也不是 Eliashberg 数值解。=> 本节检验的命题是
      「**AD 拟合公式本身**是否越过天花板」，与「真实材料是否越过」是两个命题。
      **不得混报。**
    ★ **方向审计（单边界方向规则）**：B2/B6 实测印刷约定偏向**截断** =>
      印出的 w_2 <= 真值 => ceiling = C_inf*sqrt(lam)*w_2 被**低估**
      => r_true = Tc/ceiling 被**高估**。
      判据是 r_true < 1：**高估后仍 < 1 ⟹ PASS 是声音的**（真 r_true 更小）。
      反之若出现 r_true >= 1，**不足以判 FAIL**——那正是单边行只能出「过」的情形。
      本节 13x2+... 行全部 PASS ⟹ 结论方向正确。
    严格性：**numerical-evidence**（float64，输入是论文 1-3 位有效数字的印出值）。""")

# ============================================================ §F eta/M 与 sqrt(lam)w_2
sec("§F ★ eta/M 路线与 w_2 路线是**同一个数**（自由度审计，承重）")

print("""Quan 正文（一手）："The strong coupling limit explored by Allen and Dynes gives
   Tc -> 0.18 sqrt(eta/M) = 0.18 sqrt(N(0)I^2/M), ... also indicating the seeming
   irrelevance of frequencies."
=> 这条 eta/M 一行式**又一次**出现在文献里（此前记于 LAD79）。
   但它写的是 **->（极限）**，不是 **<=（上界）**；对全 lam 的**方向**仍未被声明。

代数事实：在 lam = eta/kappa、kappa = M*w_2^2 之下
      sqrt(lam)*w_2 = sqrt(lam*w_2^2) = sqrt(eta/(M w_2^2) * w_2^2) = sqrt(eta/M)
   —— **恒等式**。=> 二者在**单原子分解**下自由度相同，不是两次检验。""")
kfac = 1.0 / math.sqrt(C_KAPPA)
print(f"\n    {'row':16s} {'sqrt(lam)w_2':>12s} {'sqrt(eta/M)*k':>13s} {'rel%':>7s}")
nF = nokF = 0
for d in DATA:
    at, lam, w2, eta = d[2], d[11], d[5], d[8]
    if at == "T" or at not in MASS or None in (lam, w2, eta):
        continue
    nF += 1
    lhs = math.sqrt(lam) * w2
    rhs = math.sqrt(eta / MASS[at]) * kfac
    rel = (rhs - lhs) / lhs * 100
    ok = abs(rel) < 6.0
    nokF += ok
    print(f"    {d[0]+'@'+str(d[1])+' '+at:16s} {lhs:12.2f} {rhs:13.2f} {rel:+7.2f}"
          f"{'' if ok else '  OUT'}")
print(f"    {nokF}/{nF} 行在 6% 内（eta 只印 2-3 位有效数字 => 主误差来自它的舍入）")
print("""
    => **承重结论**：报「我们同时用 eta 路线与 w_2 路线验证」是**虚假的双重检验**
      （单原子情形下二者恒等）。只有**多原子总量**（lam_T = sum lam_j，无单一 M）
      才使两条路线分离；那时 eta/M 形式需要额外约定（用哪个 M），
      而 sqrt(lam_T)*w_2,T 形式不需要 => **多原子下 w_2 形式才是可用的那一条**。
      同形：memory [gate-dof-audit-column-values-and-saturation]。""")

# ============================================================ §G Leavens 面积
sec("§G Quan 的 Leavens 面积线性拟合：自己重跑最小二乘")


def lsq(pts):
    n = len(pts)
    sx = sum(p[0] for p in pts); sy = sum(p[1] for p in pts)
    sxx = sum(p[0] ** 2 for p in pts); sxy = sum(p[0] * p[1] for p in pts)
    S = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    return n, S, (sy - S * sx) / n


print("论文自述：least-squares fit to Tc = S*A + T_o 给 S = 0.150, T_o = -6 K")
print("（正文说的是 'a linear relation between Tc and H (not total) parameters'）\n")
for tag, sel in (("全部 Tc>0 行", lambda d: True),
                 ("只 T 行", lambda d: d[2] == "T"),
                 ("只 H 行", lambda d: d[2] == "H")):
    pts = [(d[10], d[12]) for d in DATA if d[12] and d[12] > 0 and d[10] and sel(d)]
    n, S, T0 = lsq(pts)
    print(f"    {tag:14s} n={n:2d}  S = {S:.4f}   T_o = {T0:+7.2f} K")
print("""
    => 与 memory [leavens-1975-area-ceiling-and-crossover]、
      [same-information-class-comparison-vs-leavens] 的接口：
      Quan 在**氢化物**上重新标定了 Leavens 的面积-Tc 线性关系。
      但由 §B6（A = lam*w_1/2 恒等式）：这条关系的自变量不是独立观测量
      => 它**不**构成对 Leavens 面积天花板的新证据，只是同一组 (lam, w_1) 的重参数化。
      论文未说用哪些行拟合 => 上表只作量级对照，**不当复现门**。""")

# ============================================================ §H 敏感度门
sec("§H ★ 把 §B2 的源表不一致换算成**下游成本单位**，再判它承不承重")

print(f"""§B2 测到：kappa 列隐含的 w_2 比印出值高至多 {MAX_DIST_W2:.2f} meV（4/26 行 ~1.0-1.1），
两种印刷约定都解释不了。**不拍阈值判它过不过**——按「换了界就要重解」的要求，
把这个 delta 朝**使判据变难**的方向推进 §D/§E，算出各自的**翻转预算**，再比。

方向审计（每条单列）：
  §E 判据 r_true = Tc/(C_inf*sqrt(lam)*w_2*K) < 1 ⟹ 变难 = w_2 **变小**。
     （kappa 的证据其实指向 w_2 **更大** ⟹ 更容易过；这里仍按变难方向做，才是声音的。）
  §D 判据「Quan 的 LaH10 r 区间下端 > 我 @163 的 r={r163[0]:.4f}」⟹ 变难 = Quan 的 w_2 变小。
""")

# --- §E 的翻转预算：每行需要 w_2 下降多少才使 r_true 达到 1
eE = []
for d in DATA:
    lam, w2, tc = d[11], d[5], d[12]
    if None in (lam, w2) or not tc:
        continue
    r_true = tc / (C_INF * math.sqrt(lam) * w2 * MEV_TO_K)
    eE.append((w2 * (1.0 - r_true), f"{d[0]}@{d[1]} {d[2]}", r_true))
eE.sort()
print("§E 翻转预算（w_2 需下降多少 meV 才使该行 r_true 达到 1）——最小三行：")
for dlt, nm, r_true in eE[:3]:
    print(f"    {nm:18s} r_true={r_true:.4f}  需 w_2 下降 {dlt:7.2f} meV")
budE = eE[0][0]

# --- §D 的翻转预算：Quan 的 w_2 需下降多少才让 r 掉到 @163 的 1.1302
eD = []
for key in (("LaH10", 250), ("LaH10", 300)):
    d = TROWS[key]
    wl, w2 = d[3], d[5]
    eD.append((w2 - r163[0] * wl, f"{key[0]}@{key[1]}", w2 / wl))
eD.sort()
print("\n§D 翻转预算（Quan 的 w_2 需下降多少 meV 才使其 r 掉到我 @163 的 r 之下）：")
for dlt, nm, r in eD:
    print(f"    {nm:14s} r={r:.4f}  需 w_2 下降 {dlt:7.2f} meV")
budD = eD[0][0]

bud = min(budE, budD)
print(f"""
判决：
    观测到的不一致          delta = {MAX_DIST_W2:.2f} meV
    §E 最小翻转预算                {budE:.2f} meV   （余量 {budE/MAX_DIST_W2:.0f}x）
    §D 最小翻转预算                {budD:.2f} meV   （余量 {budD/MAX_DIST_W2:.0f}x）
    => 最紧预算 {bud:.2f} meV >> {MAX_DIST_W2:.2f} meV ⟹ **该不一致不承重**，
       §D/§E 的判决在它下面纹丝不动。按「预算优先」：**不派活去追查这 4 行**。
    ⚠ 本条只免除 §D/§E。任何**新**的下游若直接消费 kappa 或 w_2 的第三位有效数字，
      必须**重新**算一次预算，不得引用本条。""")
if bud <= MAX_DIST_W2:
    bad(f"§B2 的不一致 {MAX_DIST_W2:.2f} meV 达到了下游翻转预算 {bud:.2f} meV，承重，必须追查")

# --- 真把 w_2 推到变难方向重跑 §E（不是只比数，按「换了界就要重解」重算一遍）
worstE = max((tc / (C_INF * math.sqrt(lam) * (w2 - MAX_DIST_W2) * MEV_TO_K), f"{m}@{P} {a}")
             for (m, P, a, _wl, _w1, w2, _n, _i, _e, _k, _A, lam, tc) in DATA
             if lam and w2 and tc)
print(f"    重解：w_2 全部下推 {MAX_DIST_W2:.2f} meV 后，最紧行 = {worstE[1]}，"
      f"r_true = {worstE[0]:.4f} < 1 ⟹ 仍 PASS。")
if worstE[0] >= 1:
    bad(f"下推 w_2 后 {worstE[1]} 越界 r_true={worstE[0]:.4f}")

# ============================================================ 汇总
sec("汇总")
if FAIL:
    print(f"{len(FAIL)} 处 FAIL：")
    for f in FAIL:
        print("  -", f)
    sys.exit(1)
print("ALL CHECKS PASSED")
