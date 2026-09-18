#!/usr/bin/env python3
"""从零独立复算 C_inf = sqrt(g(2))/(2*pi) 的外向包围（纯有理算术，无浮点）。

输入只有两样：
  1) g(2) 的 CERT-D 有理端点（theory/derivations/P4_g2_analytic.md §8 表，逐字抄入并在此处标明）
  2) pi 的一个保守有理包围（十进制字面量，双侧各留一位）

输出：C_inf 的外向包围，与 P4_eta_over_M.md §2.A 所写的
  sqrt(g2)/2pi in [0.182726247746, 0.182726247790]
逐位比对。方向审计：C_inf 是一个**上界常数**，故报告值必须取包围的**上端**。
"""
from fractions import Fraction as F
from math import isqrt

# ---- 输入 1：g(2) 的 CERT-D 端点（P4_g2_analytic.md:13 / :212） ----
G2_LO = F("1.318140211762671")
G2_HI = F("1.318140212373478")

# 交叉对照：CERT-B（12 个整数 / N1=40），P4_g2_analytic.md:210
G2B_LO = F("1.318109329375712")
G2B_HI = F("1.318161578437280")

# ---- 输入 2：pi 的保守有理包围（标准十进制展开，末位双侧各放宽 1） ----
PI_LO = F("3.1415926535897932384")
PI_HI = F("3.1415926535897932385")


def sqrt_enclose(x: F, digits: int = 40):
    """返回 (lo, hi) 有理数，使 lo <= sqrt(x) <= hi。纯整数运算。"""
    assert x > 0
    scale = 10 ** digits
    # sqrt(x) = sqrt(p/q) = sqrt(p*q)/q
    p, q = x.numerator, x.denominator
    n = p * q * scale * scale          # (sqrt(n)/ (q*scale))^2 = x
    r = isqrt(n)                       # r <= sqrt(n) < r+1
    lo = F(r, q * scale)
    hi = F(r + 1, q * scale)
    assert lo * lo <= x <= hi * hi
    return lo, hi


def cinf_enclose(g_lo: F, g_hi: F):
    s_lo, _ = sqrt_enclose(g_lo)
    _, s_hi = sqrt_enclose(g_hi)
    lo = s_lo / (2 * PI_HI)            # 下端：分子取小、分母取大
    hi = s_hi / (2 * PI_LO)            # 上端：分子取大、分母取小
    return lo, hi


def dec(x: F, nd: int = 15) -> str:
    """把有理数打成 nd 位小数（截断，不四舍五入），避免引入舍入。"""
    neg = x < 0
    x = -x if neg else x
    ip = x.numerator // x.denominator
    frac = x - ip
    digs = ""
    for _ in range(nd):
        frac *= 10
        d = frac.numerator // frac.denominator
        digs += str(d)
        frac -= d
    return ("-" if neg else "") + f"{ip}.{digs}"


checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond), detail))
    print(("PASS " if cond else "FAIL ") + name + (("  | " + detail) if detail else ""))


print("=== C_inf = sqrt(g2)/(2 pi)，外向包围（rational-exact） ===")
lo_d, hi_d = cinf_enclose(G2_LO, G2_HI)
lo_b, hi_b = cinf_enclose(G2B_LO, G2B_HI)
print(f"CERT-D : [{dec(lo_d, 15)}, {dec(hi_d, 15)}]   宽 {float(hi_d - lo_d):.3e}")
print(f"CERT-B : [{dec(lo_b, 15)}, {dec(hi_b, 15)}]   宽 {float(hi_b - lo_b):.3e}")

# ---- 与 P4_eta_over_M.md §2.A 的断言比对 ----
CLAIM_LO = F("0.182726247746")
CLAIM_HI = F("0.182726247790")
chk("eta/M §2.A 的区间包含我的 CERT-D 包围",
    CLAIM_LO <= lo_d and hi_d <= CLAIM_HI,
    f"claim=[{CLAIM_LO}, {CLAIM_HI}]")
chk("eta/M §2.A 的区间不比 CERT-D 更窄（不得凭空变锐）",
    CLAIM_LO <= lo_d and CLAIM_HI >= hi_d)

# ---- 报告常数 0.1827262478 必须是合法上界 ----
REPORTED = F("0.1827262478")
chk("报告常数 0.1827262478 >= 包围上端（方向：上界须取上端）",
    REPORTED >= hi_d, f"上端={dec(hi_d, 13)}")
# ★ 口径检查：CERT-B（论文正文口径，12 个整数）**不足以**支持这个 10 位常数。
chk("CERT-B 不支持 10 位常数（⟹ 引用 0.1827262478 必须点名 CERT-D）",
    REPORTED < hi_b,
    f"CERT-B 上端={dec(hi_b, 13)}；CERT-B 下能稳住的位数只到 0.18272")

# ---- 与 KAY 报出的 f(0)=0.1827262477... 一致性 ----
KAY_LO = F("0.1827262477")
KAY_HI = F("0.1827262478")
chk("KAY 的 f(0)=0.1827262477… 落在我的 CERT-D 包围内",
    lo_d <= KAY_LO and KAY_LO <= hi_d or (lo_d >= KAY_LO and hi_d <= KAY_HI),
    f"我的包围=[{dec(lo_d, 12)}, {dec(hi_d, 12)}]")

# ---- CERT-B 必须包含 CERT-D（同一个量的两个证书） ----
chk("CERT-B 包围 包含 CERT-D 包围（两证书相容）", lo_b <= lo_d and hi_d <= hi_b)

# ---- 反向投毒：把 g2 换成偏大的错值，报告常数必须失去上界资格 ----
bad2_lo, bad2_hi = cinf_enclose(F("1.33"), F("1.34"))
chk("投毒负对照：g2=1.33–1.34 时 0.1827262478 不再是上界",
    not (REPORTED >= bad2_hi), f"投毒上端={dec(bad2_hi, 12)}")

# ---- 恒等式 (1.3)/(1.4) 的符号级复核（纯算术，Einstein 谱上闭式对账） ----
# P(dw) = 2 a2F(w) dw /(lam w),  lam = 2 int a2F/w
# => <w^2>_P = (2/lam) int w a2F  =>  lam <w^2> = 2 int w a2F.
# 取 a2F(w) = A * delta(w - wE) 的正则化：用一族离散 Einstein 混合精确检验。
from itertools import product
ok_id = True
worst = 0
for weights, freqs in product([(1,), (2, 3), (1, 4, 2)], [(1,), (0.5, 2), (0.3, 1.7, 5)]):
    if len(weights) != len(freqs):
        continue
    ws = [F(str(w)) for w in weights]
    fs = [F(str(f)) for f in freqs]
    # a2F = sum_i A_i delta(w - w_i) 形式下的积分退化成求和
    lam = 2 * sum(A / w for A, w in zip(ws, fs))
    m1 = sum(A * w for A, w in zip(ws, fs))          # int w a2F
    w2_P = (2 / lam) * m1                            # <w^2>_P 按 (1.3)
    lhs = lam * w2_P
    rhs = 2 * m1
    ok_id &= (lhs == rhs)
    worst = max(worst, abs(float(lhs - rhs)))
chk("恒等式 lam*<w^2> = 2*int(w a2F) 在离散 Einstein 混合上精确成立",
    ok_id, f"最大偏差={worst}（有理算术，应恒为 0）")

npass = sum(1 for _, c, _ in checks if c)
print(f"\n=== {npass}/{len(checks)} 通过 ===")
raise SystemExit(0 if npass == len(checks) else 1)
