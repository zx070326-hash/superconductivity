"""
P1 v0.2: Allen-Dynes (含 f1 f2) 与纯 McMillan 的 ∂lnTc/∂lnλ 对比，及非谐系数 (1/4 - s/2)。
回应 reviewer F1：纯 McMillan 在强耦合饱和（s→0），Allen-Dynes 才给 s→1/2。
运行：python theory/derivations/p1_ad_sensitivity.py
"""
import sys
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def ln_tc_mcmillan(lam, mu):
    return -1.04 * (1 + lam) / (lam - mu * (1 + 0.62 * lam))


def ln_tc_allen_dynes(lam, mu, r):
    """ln(Tc/omega_log) with f1 f2; r = omega_2/omega_log."""
    L1 = 2.46 * (1 + 3.8 * mu)
    L2 = 1.82 * (1 + 6.3 * mu) * r
    f1 = (1 + (lam / L1) ** 1.5) ** (1 / 3)
    # AD75 Eq.(35) 的形状因子是 (r − 1)，**不是** (r² − 1)。
    # 2026-09-25 订正（原写 (r**2 - 1)）。两条独立见证：
    #   (a) AD75 Table I 29 行零参数复现：(r−1) 3.85% vs (r²−1) 4.84%
    #       —— theory/notes/_scratch_ir/_chk_ad75_table1_frak_f_2026-09-20.py
    #   (b) arXiv:2010.12225 Eq.(S6) 第三方重印（同式，无平方、无开方）
    # 重解记录：theory/notes/_scratch_ir/_chk_main_f2_defect_repair_and_slambda_rescan_2026-09-25.py
    f2 = 1 + (r - 1) * lam**2 / (lam**2 + L2**2)
    return np.log(f1 * f2 / 1.2) + ln_tc_mcmillan(lam, mu)


def s_lambda(fn, lam, mu, h=1e-4, **kw):
    """s = d lnTc / d ln lambda (numerical)."""
    return (fn(lam * (1 + h), mu, **kw) - fn(lam * (1 - h), mu, **kw)) / (2 * h)


def alpha_coulomb(fn, lam, mu, **kw):
    """alpha with fixed bare mu, mu* = mu/(1+mu L), L -> L + (1/2) d ln M (omega_c ∝ M^-1/2).
    d mu*/d ln M = +mu*^2/2  →  alpha = 1/2 - (d lnTc/d mu*) * mu*^2 / 2 ... sign: alpha = -dlnTc/dlnM
    = 1/2 - (∂lnTc/∂mu*)(∂mu*/∂lnM) with ∂mu*/∂lnM = -(mu*^2)(∂L/∂lnM) = -mu*^2 * (1/2)?  L = ln(E_F/omega_c),
    omega_c ∝ M^-1/2 → L increases with M: ∂L/∂lnM = +1/2 → ∂mu*/∂lnM = -mu*^2/2.
    alpha = 1/2 - (∂lnTc/∂mu*)(-mu*^2/2) = 1/2 + (mu*^2/2)(∂lnTc/∂mu*); ∂lnTc/∂mu* < 0 → alpha < 1/2."""
    h = 1e-5
    d = (fn(lam, mu + h, **kw) - fn(lam, mu - h, **kw)) / (2 * h)
    return 0.5 + 0.5 * mu**2 * d


mu = 0.13
print("=== s = ∂lnTc/∂lnλ 与非谐系数 c = 1/4 - s/2 （μ*=0.13）===")
print(f"{'λ':>6} | {'McM s':>7} {'c':>7} | {'AD r=1.2 s':>10} {'c':>7} | {'AD r=1.5 s':>10} {'c':>7} | {'AD r=2.0 s':>10} {'c':>7}")
for lam in [0.5, 1.0, 1.5, 1.84, 2.0, 2.64, 3.0, 4.0, 6.0, 10.0, 100.0]:
    sm = s_lambda(ln_tc_mcmillan, lam, mu)
    row = f"{lam:6.2f} | {sm:7.3f} {0.25 - sm/2:7.3f}"
    for r in [1.2, 1.5, 2.0]:
        sa = s_lambda(ln_tc_allen_dynes, lam, mu, r=r)
        row += f" | {sa:10.3f} {0.25 - sa/2:7.3f}"
    print(row)

print("\n=== 库仑项 α（固定裸 μ，ω_c ∝ M^-1/2），λ=2 ===")
print(f"{'μ*':>5} | {'McMillan':>9} | {'AD r=1.5':>9}")
for m in [0.10, 0.13, 0.20, 0.25, 0.30]:
    print(f"{m:5.2f} | {alpha_coulomb(ln_tc_mcmillan, 2.0, m):9.3f} | {alpha_coulomb(ln_tc_allen_dynes, 2.0, m, r=1.5):9.3f}")

print("\n=== H3S 例子：η≈0.43，λ_anh=1.84 ===")
eta = 1 / 0.70 - 1
for name, fn, kw in [("McMillan", ln_tc_mcmillan, {}), ("AD r=1.5", ln_tc_allen_dynes, {"r": 1.5})]:
    s = s_lambda(fn, 1.84, mu, **kw)
    print(f"{name:9s}: s={s:.3f}, c={0.25 - s/2:+.3f}, δα_anh = c·η = {(0.25 - s/2) * eta:+.3f}")
