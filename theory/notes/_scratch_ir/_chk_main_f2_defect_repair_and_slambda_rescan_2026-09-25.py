"""
从零复核：修 `f₂` 的两处仓库缺陷，并**重解**所有建在它上面的下游。

缺陷（两条独立见证）：
  (a) 我自己的 AD75 Table I 29 行零参数复现指纹（_chk_ad75_table1_frak_f_2026-09-20.py）
  (b) arXiv:2010.12225 Eq.(S6) 第三方重印（页图 hrE_p-28.png）
  正确式： f2 = 1 + (r-1)·λ²/(λ² + Λ₂²) ,  Λ₂ = 1.82(1+6.3μ*)·r     [AD75 Eq.(35)]

消费点（本脚本逐个重解，不是只改一行）：
  C1  P1_isotope_master_formula.md §3.1  表「α，Allen–Dynes (r=1.5)」λ=2 五列
  C2  P1_isotope_master_formula.md §3.2a 表「c，AD r=1.5」八列
  C3  P1_isotope_master_formula.md §3.2a 正文「r=1.2→2.0 时 λ=2 处 c 从 −0.26 到 −0.38」
  C4  P1_isotope_master_formula.md §3.2a **承重**符号规则 s_λ ≥ ½ 的 3.7 万格点扫描
  C5  P4_strong_coupling_expansion.md §5 推论 S.3 的 f₂ 闭式 + 小 σ 展开

⚠ C4 的**方向**（先审再跑）：错式 (r²−1) > (r−1) 对 r>1，且 f₂ 随 λ 单调增
   ⟹ 错式**高估** ∂ln f₂/∂lnλ ⟹ **高估** s_λ。
   而判据是 s_λ ≥ ½ 的**下界**。高估一个被要求「≥」的量 = 单边界方向规则里的
   「只支持 PASS」那一侧 ⟹ 旧扫描的「无一低于 ½」**可能是缺陷制造的假 PASS**，必须重扫。

运行：OMP_NUM_THREADS=1 python theory/notes/_scratch_ir/_chk_main_f2_defect_repair_and_slambda_rescan_2026-09-25.py
"""
import math
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

FAILS = []


def bad(msg):
    FAILS.append(msg)
    print("  [FAIL] " + msg)


# ---------------------------------------------------------------- 公式
def f1_of(lam, mus):
    return (1 + (lam / (2.46 * (1 + 3.8 * mus))) ** 1.5) ** (1 / 3)


def f2_right(lam, mus, r):
    """AD75 Eq.(35)：(r-1)，无开方。"""
    L2 = 1.82 * (1 + 6.3 * mus) * r
    return 1 + (r - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)


def f2_wrong_sq(lam, mus, r):
    """p1_ad_sensitivity.py:23 的旧写法 (r**2-1)。"""
    L2 = 1.82 * (1 + 6.3 * mus) * r
    return 1 + (r ** 2 - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)


def f2_wrong_sqrt(lam, mus, r):
    """P4_strong_coupling_expansion.md:123 的旧写法：整块开方。"""
    L2 = 1.82 * (1 + 6.3 * mus) * r
    return math.sqrt((1 + (r - 1) * lam ** 2) / (lam ** 2 + L2 ** 2))


def ln_tc_mc(lam, mus):
    return -1.04 * (1 + lam) / (lam - mus * (1 + 0.62 * lam))


def ln_tc_ad(lam, mus, r, f2fn):
    return math.log(f1_of(lam, mus) * f2fn(lam, mus, r) / 1.2) + ln_tc_mc(lam, mus)


def s_lambda(fn, lam, mus, h=1e-4, **kw):
    return (fn(lam * (1 + h), mus, **kw) - fn(lam * (1 - h), mus, **kw)) / (2 * h)


def alpha_coulomb(fn, lam, mus, **kw):
    h = 1e-5
    d = (fn(lam, mus + h, **kw) - fn(lam, mus - h, **kw)) / (2 * h)
    return 0.5 + 0.5 * mus ** 2 * d


def lam_min(mus):
    """分母 λ − μ*(1+0.62λ) = 0 的根。"""
    return mus / (1 - 0.62 * mus)


# ---------------------------------------------------------------- §0 自检
print("=" * 78)
print("§0 三个写法在 AD75 拟合域中心点上互不相同（确认它们真是三个式子）")
print("=" * 78)
lam0, mu0, r0 = 1.0, 0.13, 1.5
v = [f2_right(lam0, mu0, r0), f2_wrong_sq(lam0, mu0, r0), f2_wrong_sqrt(lam0, mu0, r0)]
print(f"  λ=1.0 μ*=0.13 r=1.5:  正确 {v[0]:.6f}   (r²−1) {v[1]:.6f}   开方 {v[2]:.6f}")
if not (abs(v[0] - v[1]) > 1e-3 and abs(v[0] - v[2]) > 1e-3):
    bad("§0 三写法未区分开")
# 方向断言：r>1 时错式 (r²−1) 必 ≥ 正确式
nsq = 0
for lam in [0.5, 1, 2, 5, 20]:
    for mus in [0.10, 0.13, 0.20, 0.30]:
        for r in [1.05, 1.2, 1.5, 2.0, 2.5]:
            if f2_wrong_sq(lam, mus, r) < f2_right(lam, mus, r) - 1e-15:
                nsq += 1
print(f"  方向检验：r>1 上 (r²−1) ≥ (r−1) 逐点成立 —— 违反 {nsq}/100")
if nsq:
    bad("§0 方向断言被违反")

# ---------------------------------------------------------------- §C1
print()
print("=" * 78)
print("C1  P1 §3.1 表：α（固定裸 μ，ω_c ∝ M^-1/2），λ=2")
print("=" * 78)
print(f"  {'μ*':>5} | {'McMillan':>9} | {'AD r=1.5 旧(r²−1)':>18} | {'AD r=1.5 新(r−1)':>17} | {'Δ':>7}")
C1_NEW = {}
for m in [0.10, 0.13, 0.20, 0.25, 0.30]:
    a_mc = alpha_coulomb(ln_tc_mc, 2.0, m)
    a_old = alpha_coulomb(ln_tc_ad, 2.0, m, r=1.5, f2fn=f2_wrong_sq)
    a_new = alpha_coulomb(ln_tc_ad, 2.0, m, r=1.5, f2fn=f2_right)
    C1_NEW[m] = a_new
    print(f"  {m:5.2f} | {a_mc:9.3f} | {a_old:18.3f} | {a_new:17.3f} | {a_new - a_old:+7.3f}")

# ---------------------------------------------------------------- §C2/C3
print()
print("=" * 78)
print("C2  P1 §3.2a 表：c = ¼ − s_λ/2")
print("=" * 78)
LAMS_C2 = [1.0, 1.5, 1.84, 2.0, 2.64, 3.0, 6.0, 100.0]
mu = 0.13
row_mc, row_old, row_new = [], [], []
for lam in LAMS_C2:
    row_mc.append(0.25 - s_lambda(ln_tc_mc, lam, mu) / 2)
    row_old.append(0.25 - s_lambda(ln_tc_ad, lam, mu, r=1.5, f2fn=f2_wrong_sq) / 2)
    row_new.append(0.25 - s_lambda(ln_tc_ad, lam, mu, r=1.5, f2fn=f2_right) / 2)
print("  λ           " + "".join(f"{x:>9.2f}" for x in LAMS_C2))
print("  c McMillan  " + "".join(f"{x:>9.3f}" for x in row_mc))
print("  c AD 旧     " + "".join(f"{x:>9.3f}" for x in row_old))
print("  c AD 新     " + "".join(f"{x:>9.3f}" for x in row_new))
print(f"\n  「AD 的 c 恒负」在这 8 点上：旧 {sum(x < 0 for x in row_old)}/8，新 {sum(x < 0 for x in row_new)}/8")

print()
print("C3  正文断言「r=1.2→2.0 时 λ=2 处 c 从 −0.26 到 −0.38」")
for r in (1.2, 2.0):
    c_old = 0.25 - s_lambda(ln_tc_ad, 2.0, mu, r=r, f2fn=f2_wrong_sq) / 2
    c_new = 0.25 - s_lambda(ln_tc_ad, 2.0, mu, r=r, f2fn=f2_right) / 2
    print(f"    r={r}:  旧 {c_old:+.4f}   新 {c_new:+.4f}")

# ---------------------------------------------------------------- §C4 承重
print()
print("=" * 78)
print("C4  ★承重★ 符号规则 s_λ ≥ ½ 的重扫（旧扫描用的是错的 f₂）")
print("=" * 78)
print("  网格：μ* 21 点 [0.10,0.30]，r 27 点 [1.2,2.5]，λ 对数 70 点 (λ_min(μ*)·1.001, 5000]")


def scan(f2fn, tag):
    worst = (1e9, None)
    n_below, n_tot = 0, 0
    for i in range(21):
        mus = 0.10 + 0.20 * i / 20
        lo = lam_min(mus) * 1.001
        for j in range(27):
            r = 1.2 + 1.3 * j / 26
            for k in range(70):
                lam = lo * (5000.0 / lo) ** (k / 69)
                s = s_lambda(ln_tc_ad, lam, mus, r=r, f2fn=f2fn)
                n_tot += 1
                if s < 0.5:
                    n_below += 1
                if s < worst[0]:
                    worst = (s, (lam, mus, r))
    print(f"  [{tag}] 格点 {n_tot}，s_λ < ½ 的点 {n_below}，最小 s_λ = {worst[0]:.6f} "
          f"@ λ={worst[1][0]:.4g} μ*={worst[1][1]:.3f} r={worst[1][2]:.3f}")
    return worst, n_below, n_tot


w_old, nb_old, nt = scan(f2_wrong_sq, "旧 (r²−1)")
w_new, nb_new, _ = scan(f2_right, "新 (r−1) ")

print()
if nb_new == 0:
    print("  ⟹ 判决：符号规则 s_λ ≥ ½ 在**正确** f₂ 下**依然成立**（0 个反例）。")
    print("     旧扫描的结论侥幸为真——但它的**证据是无效的**（用了错式，且方向是高估）。")
    print("     两个独立命题：`旧证据无效` 与 `新结论为真`，此处分别成立。")
else:
    print(f"  ⟹ ★判决翻转★：正确 f₂ 下有 {nb_new} 个格点 s_λ < ½ ⟹ 符号规则不再被数值支持。")
print(f"  余量对比：旧最小 {w_old[0]:.6f} → 新最小 {w_new[0]:.6f}"
      f"（新比旧{'低' if w_new[0] < w_old[0] else '高'} {abs(w_new[0]-w_old[0]):.6f}）")
if w_new[0] > w_old[0] + 1e-9:
    bad("C4 方向断言被违反：正确 f₂ 反而给出更大的 s_λ 下界（与 §0 的单调性矛盾）")

# 端点自检：λ→∞ 应趋 ½（f₁ 贡献 ½，f₂ 趋常数）
s_inf_new = s_lambda(ln_tc_ad, 1e7, 0.13, r=1.5, f2fn=f2_right)
s_inf_old = s_lambda(ln_tc_ad, 1e7, 0.13, r=1.5, f2fn=f2_wrong_sq)
print(f"  渐近自检 λ=1e7：新 s_λ={s_inf_new:.6f}，旧 s_λ={s_inf_old:.6f}（两者都应 → 0.5⁺）")
if abs(s_inf_new - 0.5) > 1e-3 or abs(s_inf_old - 0.5) > 1e-3:
    bad("C4 渐近自检未趋 ½")

# ---------------------------------------------------------------- §C5
print()
print("=" * 78)
print("C5  P4 §5 推论 S.3 的开方式：它错的**不只是**一个开方")
print("=" * 78)
print("  旧式  f₂ = {[1 + (r−1)λ²] / [λ² + Λ₂²]}^{1/2}")
print("  正确  f₂ = 1 + (r−1)λ²/(λ² + Λ₂²)")
print("  ⟹ 旧式把 `1 +` 挪进了分子、并让整块开方 ⟹ 量纲与极限全错。检验两个极限：")
for tag, lam in [("λ→0", 1e-6), ("λ→∞", 1e6)]:
    print(f"    {tag}:  正确 {f2_right(lam, 0.13, 1.5):.6f}   旧式 {f2_wrong_sqrt(lam, 0.13, 1.5):.6g}")
print("  正确式满足 f₂(λ→0)=1（无形状修正）与 f₂(λ→∞)=r（饱和于形状比）；")
lim_hi = f2_right(1e9, 0.13, 1.5)
print(f"    实测 f₂(λ=1e9) = {lim_hi:.6f}，应 = r = 1.5")
if abs(lim_hi - 1.5) > 1e-6 or abs(f2_right(1e-6, 0.13, 1.5) - 1.0) > 1e-9:
    bad("C5 正确式的两个极限不符")
print("  旧式两个极限都不是这两个数 ⟹ **不是「少一次代数化简」，是抄错了整条式子**。")

# 小 σ 展开：正确式的 σ² 系数（供 P4 对照表填 AD 侧）
print()
print("  小 σ 展开（正确式，供 P4 对照表的 AD 列）：设 r = 1 + βσ² + O(σ⁴)，则")
print("    f₂ = 1 + (r−1)λ²/(λ²+Λ₂²) = 1 + βσ²·λ²/(λ²+Λ₂(1)²) + O(σ⁴)，Λ₂(1)=1.82(1+6.3μ*)")
print("    ⟹ σ² 系数 = β·λ²/(λ² + [1.82(1+6.3μ*)]²)   —— 闭式，**恒正**，且随 λ 单调增")
for lam in (1.0, 2.0, 5.0):
    L2 = 1.82 * (1 + 6.3 * 0.13)
    print(f"      λ={lam}: 系数/β = {lam**2/(lam**2+L2**2):.6f}")
print("    （旧稿写的 `(βλ² − (1.82…)²β)/2` 由开方式导出 ⟹ 连同母式一并作废）")

# ---------------------------------------------------------------- 汇总
print()
print("=" * 78)
if FAILS:
    print(f"FAILED ({len(FAILS)}):")
    for f in FAILS:
        print("  - " + f)
    sys.exit(1)
print("ALL CHECKS PASSED")
