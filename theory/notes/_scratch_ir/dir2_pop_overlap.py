"""
方向2（μ*>0 Einstein 极值性弱耦合反转）人口重叠计算。
数据源：theory/data/ambient_hydride_candidates.csv (mu*=0.10, AD 公式约定)
        theory/data/hydrides_theory.csv (mustar_used 列逐行给出)
阈值来源：theory/notes/coulomb/mustar_failure_common_mechanism_2026-09-18.md §5.2
  Xi := lambda / mu*(omega_log)
  Xi < 2.6  => 反转区（"omega_log 越大越好"被反转，即天花板方向筛选启发式可疑）
  Xi > 3.1  => 安全
  两阈值本身是在 mu*_rms=0.2023/0.2771 两点测到的【单边下界】（方向规则：
  找到违反的一侧有效，故真实 s_c 只会更高，不会更低）。
  该笔记同时给出一条不同形式的估计：mustar_extremality_failure_map_2026-09-17.md §4
  在物理 mu*=0.13 外推 lambda <~ 0.3（绝对 lambda 阈值，未在算子上直接测量，
  是从弱耦合领头阶判据 lambda < 2 mu*_rms 反推的量级落点，笔记原文写"外推"）。
  本脚本对两种口径都算重叠比例，供人对照；主口径按 Xi（考虑逐材料 mu* 差异）。
"""
import csv, os

REPO = r"D:\physics\superconductivity"
f1 = os.path.join(REPO, "theory", "data", "ambient_hydride_candidates.csv")
f2 = os.path.join(REPO, "theory", "data", "hydrides_theory.csv")

rows = []  # (material, lambda, mustar, source)

with open(f1, encoding="utf-8") as fh:
    r = csv.DictReader(fh)
    for row in r:
        try:
            lam = float(row["lambda"])
        except (ValueError, KeyError):
            continue
        rows.append((row["material"], lam, 0.10, "ambient_hydride_candidates(mu*=0.10 AD convention)"))

with open(f2, encoding="utf-8") as fh:
    r = csv.DictReader(fh)
    for row in r:
        try:
            lam = float(row["lambda"])
            mus = float(row["mustar_used"])
        except (ValueError, KeyError, TypeError):
            continue
        rows.append((row["material"], lam, mus, "hydrides_theory(mustar_used column)"))

print(f"Total rows with usable (lambda, mu*): {len(rows)}")
print(f"  ambient_hydride_candidates.csv usable rows: {sum(1 for r_ in rows if 'ambient' in r_[3])}")
print(f"  hydrides_theory.csv usable rows: {sum(1 for r_ in rows if 'hydrides_theory' in r_[3])}")

# --- histogram of lambda ---
lams = sorted(r_[1] for r_ in rows)
print("\n=== lambda distribution (all rows) ===")
bins = [0,0.3,0.5,0.7,1.0,1.5,2.0,2.5,3.0,10.0]
counts = [0]*(len(bins)-1)
for lam in lams:
    for i in range(len(bins)-1):
        if bins[i] <= lam < bins[i+1]:
            counts[i]+=1
            break
for i,c in enumerate(counts):
    print(f"  [{bins[i]:.2f},{bins[i+1]:.2f}): {c}  {'#'*c}")
print(f"  min={min(lams):.4f} max={max(lams):.4f} n={len(lams)}")

# --- absolute lambda threshold 0.3 (mustar_extremality_failure_map extrapolation) ---
n_below_03 = sum(1 for _,lam,_,_ in rows if lam < 0.3)
n_total = len(rows)
print(f"\n=== 口径 A：绝对阈值 lambda < 0.3 (mu*=0.13 外推口径) ===")
print(f"  below (反转候选侧): {n_below_03}/{n_total} = {100*n_below_03/n_total:.1f}%")
print(f"  above (正常侧): {n_total-n_below_03}/{n_total} = {100*(n_total-n_below_03)/n_total:.1f}%")
below03 = [(m,lam,mu) for m,lam,mu,_ in rows if lam < 0.3]
print("  低于 0.3 的材料:", below03 if below03 else "(none)")

# --- Xi = lambda/mustar threshold ---
print(f"\n=== 口径 B：Xi = lambda/mu* 阈值 (mustar_failure_common_mechanism_2026-09-18 §5.2) ===")
xis = [(m, lam, mu, lam/mu) for m,lam,mu,_ in rows]
xis_sorted = sorted(xis, key=lambda t: t[3])
print("  10 个最小 Xi 材料（离反转区最近）：")
for m,lam,mu,xi in xis_sorted[:10]:
    print(f"    {m:20s}  lambda={lam:.4f}  mu*={mu:.4f}  Xi={xi:.3f}")

for lo,hi,label in [(2.6, None, "Xi<2.6（反转，测到的下界）"),
                     (None, 3.1, "Xi>3.1（安全）"),
                     (2.6, 3.1, "2.6<=Xi<=3.1（灰区，未测）")]:
    if lo is not None and hi is not None:
        n = sum(1 for *_,xi in xis if lo<=xi<=hi)
    elif lo is not None:
        n = sum(1 for *_,xi in xis if xi<lo)
    else:
        n = sum(1 for *_,xi in xis if xi>hi)
    print(f"  {label}: {n}/{n_total} = {100*n/n_total:.1f}%")

# also using the fitted extrapolation from mustar_extremality_failure_map for physical mu*
# s_c_fit(mu*) = 1.7733 + 3.7598*mu*  (numerical-evidence fit, itself flagged as possibly low
#  by the 2026-09-18 note's direction-consistent single-sided finding s_c >= 2.62/3.07 at higher mu*)
def sc_fit(mustar):
    return 1.7733 + 3.7598*mustar

print(f"\n=== 口径 C：逐材料用拟合式 s_c_fit(mu*)=1.7733+3.7598*mu* 做阈值（口径 B 的材料相关版本）===")
n_below_fit = sum(1 for m,lam,mu,xi in xis if xi < sc_fit(mu))
print(f"  Xi < s_c_fit(mu*)（反转候选侧）: {n_below_fit}/{n_total} = {100*n_below_fit/n_total:.1f}%")
below_fit = [(m,lam,mu,xi,sc_fit(mu)) for m,lam,mu,xi in xis if xi < sc_fit(mu)]
if below_fit:
    print("  材料:")
    for m,lam,mu,xi,sc in below_fit:
        print(f"    {m:20s} lambda={lam:.4f} mu*={mu:.4f} Xi={xi:.3f} s_c_fit={sc:.3f}")
else:
    print("  (none)")
