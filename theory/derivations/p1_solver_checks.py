"""
P1 §3.1 / §4.1 的求解器数值：库仑项 δα_μ 的真实 Eliashberg 值，与闭式（McMillan / Allen-Dynes）对照。
运行：python theory/derivations/p1_solver_checks.py   （约 1–3 分钟）
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "solvers"))
from eliashberg_iso import Spectrum, Solver  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p1_ad_sensitivity import alpha_coulomb, ln_tc_allen_dynes, ln_tc_mcmillan  # noqa: E402

WCUT_RATIO = 5.0  # 求解器约定：w_cut = 5·ω_log（implementer 首轮所用）


def alpha_full(spec, mu, fixed):
    solver = Solver(spec, mu_star=mu, w_cut=WCUT_RATIO * spec.omega_log())
    a, _ = solver.isotope_coefficient({"all": spec}, which="all", mu_star_fixed=fixed)
    return a


print("=== A. 全同位素 α，Einstein 谱，λ=1.84（H3S 非谐 λ），ω0=93 meV ===")
sp = Spectrum.einstein(lam=1.84, w0=93.0, name="all")
print(f"{'μ*':>5} | {'α num (μ* 随 ω_c)':>18} | {'α num (μ* 固定)':>16} | {'McMillan':>9} | {'AD r=1':>7}")
for mu in [0.10, 0.13, 0.16, 0.20, 0.25, 0.30]:
    a_var = alpha_full(sp, mu, fixed=False)
    a_fix = alpha_full(sp, mu, fixed=True)
    print(f"{mu:5.2f} | {a_var:18.4f} | {a_fix:16.4f} | {alpha_coulomb(ln_tc_mcmillan, 1.84, mu):9.3f} | {alpha_coulomb(ln_tc_allen_dynes, 1.84, mu, r=1.0):7.3f}")

print("\n=== B. 部分同位素：H3S 型三分量谱（S 模不动，H 弯曲 + H 伸缩缩放）===")
# 近似 Errea 2015 非谐 α²F 的粗结构：S 模 ~45 meV（λ 0.25），H 弯曲 ~60 meV（λ 0.50），H 伸缩 ~150 meV（λ 1.09）
S = Spectrum.einstein(lam=0.25, w0=45.0, name="S")
Hb = Spectrum.einstein(lam=0.50, w0=60.0, name="H_bend")
Hs = Spectrum.einstein(lam=1.09, w0=150.0, name="H_stretch")
full = S.combine(Hb).combine(Hs)
parts = {"S": S, "H_bend": Hb, "H_stretch": Hs}
print(f"λ_total={full.lambda_total():.3f}, ω_log={full.omega_log():.1f} meV, ω2/ω_log={full.omega_2()/full.omega_log():.2f}")
print(f"{'μ*':>5} | {'α_H (μ* 固定)':>13} | {'α_H (μ* 随 ω_c)':>15} | {'δα_μ^(H)':>9} | {'α_all (μ* 固定)':>15} | {'w_H = α_H/α_all':>15}")
for mu in [0.10, 0.16, 0.25]:
    solver = Solver(full, mu_star=mu, w_cut=WCUT_RATIO * full.omega_log())
    aH_fix, _ = solver.isotope_coefficient(parts, which=["H_bend", "H_stretch"], mu_star_fixed=True)
    aH_var, d = solver.isotope_coefficient(parts, which=["H_bend", "H_stretch"], mu_star_fixed=False)
    aA_fix, _ = solver.isotope_coefficient(parts, which="all", mu_star_fixed=True)
    print(f"{mu:5.2f} | {aH_fix:13.4f} | {aH_var:15.4f} | {aH_var - aH_fix:+9.4f} | {aA_fix:15.4f} | {aH_fix / aA_fix:15.3f}")
print(f"  (μ*=0.25 时 ω_c 从 {d['w_cut_minus']:.1f} → {d['w_cut_plus']:.1f} meV，μ* 从 {d['mu_star_minus']:.4f} → {d['mu_star_plus']:.4f})")

print("\n=== C. YH6 型：λ=1.7，ω0=100 meV，Einstein —— α=0.40 需要多大 μ*？===")
sp = Spectrum.einstein(lam=1.7, w0=100.0, name="all")
for mu in [0.20, 0.25, 0.30, 0.35, 0.40]:
    try:
        print(f"  μ*={mu:.2f}: α = {alpha_full(sp, mu, fixed=False):.4f}   (McMillan {alpha_coulomb(ln_tc_mcmillan, 1.7, mu):.3f})")
    except Exception as ex:  # noqa: BLE001
        print(f"  μ*={mu:.2f}: {type(ex).__name__}: {ex}")
