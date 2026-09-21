# -*- coding: utf-8 -*-
"""★ 直接判定 `a` 缺口的机制：**氢化物谱是不是「高 𝔣 形状」？**

上游（`_chk_main_b_gap_bias_inversion` §1）：同阶二次模型下
    数据 a = −1.0003 ± 0.0258 ，理论 a = −0.8707 ± 0.0167 ，差 −0.1296
即**理论（用氢化物谱算）的 𝔣 比 AD75 数据高 12.2%**，而三个形状系数全相合。

上游（`_chk_main_shape_completeness_lamS`）：固定 (λ,S) 上人工形状族给出 a 的极差 0.1069
——量级够（0.83×），但那只说明「存在这么大的自由度」，**没说明氢化物落在哪一端**。

本脚本问的是**可证伪的那一问**：
    在**同一个 (λ, S)** 上，氢化物真实 α²F 的 𝔣，与一个**通用形状**（双 Einstein / 矩形块）的 𝔣，
    谁高？高多少？
预言（若 (ii) 成立）：氢化物应**系统性偏高**，且幅度 ≈ +12%。
若氢化物偏低、或幅度远小于 12%，则 (ii) 被证伪，`a` 缺口另有出处。

这比上游的「极差 vs 缺口」强，因为它**带方向**：极差是双侧的，本测量是单侧定向的。

方向注记：比值双侧点算，不产生单边界判决。
严格性：numerical-evidence（float64 Eliashberg）；「重标度下 (ω_log, ω̄₂) 不变」是 analytic，仅作实现门。
主会话 2026-09-20，从零。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "solvers"))
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K   # noqa: E402

C_INF = 0.1827262477
MU, WRATIO, NMATS = 0.10, 10.0, 400
LAM_GRID = np.exp(np.linspace(np.log(0.69), np.log(2.59), 9))
A_GAP = -0.1296                     # 数据 − 理论（同阶二次），见上游 §1

DATA = ROOT / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"
FILES = {
    "H3S-anh": DATA / "H3S-200GPa/EPW_anharmonic/H3S.a2f",
    "D3S-anh": DATA / "D3S-200GPa/EPW_anharmonic/D3S.a2f",
    "H3S-har": DATA / "H3S-200GPa/EPW_harmonic/H3S.a2f",
    "D3S-har": DATA / "D3S-200GPa/EPW_harmonic/D3S.a2f",
}


def load(path):
    w, a = [], []
    for line in path.read_text().splitlines()[1:]:
        p = line.split()
        if len(p) < 3:
            break
        try:
            x, yv = float(p[0]), float(p[1])
        except ValueError:
            break
        w.append(x); a.append(max(yv, 0.0))
    return np.array(w), np.array(a)


def bisect(f, lo, hi, tol=1e-14, it=200):
    flo, fhi = f(lo), f(hi)
    assert flo * fhi < 0, f"区间不变号：f({lo})={flo:.4g}, f({hi})={fhi:.4g}"
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if hi - lo < tol * max(1.0, abs(mid)):
            break
    return 0.5 * (lo + hi)


def sp_twopeak(lam, S, x, w2_target):
    def lnS(rho):
        return (1 - x) * np.log(rho) - 0.5 * np.log(x + (1 - x) * rho ** 2)
    rho = bisect(lambda r: lnS(r) - np.log(S), 1.0 + 1e-12, 1e6)
    w1 = w2_target / np.sqrt(x + (1 - x) * rho ** 2)
    return Spectrum.two_peak(x * lam, w1, (1 - x) * lam, rho * w1)


def sp_block(lam, S, w2_target):
    def mom(w, d):
        wl = np.exp(np.trapezoid(d / w * np.log(w), w) / np.trapezoid(d / w, w))
        w2 = np.sqrt(np.trapezoid(d * w, w) / np.trapezoid(d / w, w))
        return wl, w2
    R = bisect(lambda r: np.log(np.divide(*mom(np.linspace(1.0, r, 20001),
                                              np.ones(20001)))) - np.log(S),
               1.0 + 1e-6, 1e4)
    w = np.linspace(1.0, R, 20001)
    _, w2 = mom(w, np.ones_like(w))
    w = w * (w2_target / w2)
    sp0 = Spectrum.from_table(w, np.ones_like(w))
    return Spectrum.from_table(w, np.ones_like(w) * (lam / sp0.lambda_total()))


def frak(sp):
    t = Solver(spectrum=sp, mu_star=MU, w_cut=WRATIO * sp.omega_2(),
               n_matsubara=NMATS).tc()
    assert t > 0
    f = t / (C_INF * np.sqrt(sp.lambda_total()) * sp.omega_2() / KB_MEV_PER_K)
    assert f < 1.0, f"𝔣>1 报警：{f:.4f}"
    return f


XS = (0.30, 0.45, 0.60, 0.75)      # 去掉 x=0.15（上游判定为病态：低频权重极小）

print("=" * 100)
print("门：对每条氢化物谱，构造**同 (λ, S)** 的通用形状（4 个双峰 + 1 个矩形块）")
print("=" * 100)
RAW = {k: load(p) for k, p in FILES.items()}
REF = {}
print(f"  {'谱':<10}{'λ_ref':>9}{'S':>9}{'ω̄₂[meV]':>11}  | 通用形状命中 (λ,S) 的最劣相对误差")
for k, (w, a) in RAW.items():
    sp0 = Spectrum.from_table(w, a)
    lam0, w2, wl = sp0.lambda_total(), sp0.omega_2(), sp0.omega_log()
    S = wl / w2
    REF[k] = (lam0, S, w2)
    errs = []
    for mk in [lambda l: sp_twopeak(l, S, x, w2) for x in XS] + [lambda l: sp_block(l, S, w2)]:
        sp = mk(1.30)
        errs.append(max(abs(sp.lambda_total() / 1.30 - 1),
                        abs(sp.omega_log() / sp.omega_2() / S - 1)))
    worst = max(errs)
    print(f"  {k:<10}{lam0:>9.4f}{S:>9.4f}{w2:>11.3f}  | {worst:.2e}")
    assert worst < 5e-6, f"{k} 通用形状未命中 (λ,S)"
print("  ⟹ 每组内部 (λ, S) 逐位相同 ⟹ 𝔣 之差**纯粹**是形状之差。")

print()
print("=" * 100)
print(f"结果：𝔣(氢化物) / 𝔣(通用形状) ——  μ*={MU}, w_cut={WRATIO:.0f}·ω̄₂, N={NMATS}")
print("=" * 100)
print("  比值 > 1 ⟹ 氢化物是**高 𝔣 形状** ⟹ 支持 (ii)；预言幅度 ≈ "
      f"{100*(np.exp(-A_GAP)-1):+.1f}%")
print()
hdr = "  " + f"{'谱':<10}" + "".join(f"{l:>9.2f}" for l in LAM_GRID) + f"{'几何均':>10}"
print(hdr)
ratios_all = []
for k, (w, a) in RAW.items():
    lam0, S, w2 = REF[k]
    gen = [lambda l, x=x: sp_twopeak(l, S, x, w2) for x in XS] + [lambda l: sp_block(l, S, w2)]
    row = []
    for lv in LAM_GRID:
        fh = frak(Spectrum.from_table(w, a * (lv / lam0)))
        fg = np.array([frak(mk(lv)) for mk in gen])
        row.append(fh / float(np.exp(np.mean(np.log(fg)))))
    row = np.array(row)
    ratios_all.append(row)
    print("  " + f"{k:<10}" + "".join(f"{v:>9.4f}" for v in row)
          + f"{float(np.exp(np.mean(np.log(row)))):>10.4f}")

R = np.array(ratios_all)
lr = np.log(R)
gm = float(np.exp(lr.mean()))
print()
print("  " + f"{'四谱几何均':<10}" + "".join(f"{v:>9.4f}" for v in np.exp(lr.mean(axis=0)))
      + f"{gm:>10.4f}")

print()
print("=" * 100)
print("判决")
print("=" * 100)
obs = float(np.exp(-A_GAP))
print(f"  观测要求：理论(氢化物)/数据(AD75) 的 𝔣 比 = exp({-A_GAP:+.4f}) = {obs:.4f}"
      f"（即氢化物需比 AD75 族高 {100*(obs-1):+.1f}%）")
print(f"  实测    ：𝔣(氢化物)/𝔣(同 (λ,S) 通用形状) 的四谱×九 λ 几何均 = {gm:.4f}"
      f"（{100*(gm-1):+.1f}%）")
print(f"  覆盖率  ：log 比 {lr.mean():+.4f} / {-A_GAP:+.4f} = {lr.mean()/(-A_GAP):.2f}×")
print()
if lr.mean() <= 0:
    print("  ★ **判决：(ii) 被证伪，且是反号的。** 氢化物在同 (λ,S) 上是**低** 𝔣 形状，")
    print("     而 `a` 缺口要求它是**高** 𝔣 形状 ⟹ 残余形状自由度把差距**推大**而非解释掉。")
elif lr.mean() / (-A_GAP) >= 0.5:
    print("  ★ **判决：(ii) 方向正确且量级相当**——氢化物确是高 𝔣 形状，")
    print(f"     单这一项就买到缺口的 {100*lr.mean()/(-A_GAP):.0f}%。")
else:
    print(f"  ★ **判决：(ii) 方向正确但量级不足**（只买到 {100*lr.mean()/(-A_GAP):.0f}%）⟹ 主因另有其他。")
print()
print("  ⚠ 承重边界：")
print("    · 「通用形状」是我选的（双峰 x∈{0.30..0.75} + 矩形块），**不是** AD75 那 29 个材料的真实 α²F。")
print("      本测量给的是「氢化物 vs 一族简单形状」的方向与量级，**不是**「氢化物 vs AD75 族」的直接测量。")
print("      要消掉这一条，需要 AD75 材料的真实 α²F（隧道反演谱），那是另一件事、另算预算。")
print("    · 但**方向**这一半是稳的：比值若反号，任何「通用形状」的选择都难以救回 ——")
print("      因为 x 扫描已覆盖双峰权重从 30% 到 75%，且矩形块是完全不同的拓扑。")
print("=" * 100)
