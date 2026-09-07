"""
P1 §3.1 / §3.2 的机械推导核对（sympy）。
运行：python theory/derivations/p1_sympy_checks.py
"""
import sympy as sp

print("=== 3.1 库仑项：McMillan 形式下的 alpha ===")
lam, mu, M, L0 = sp.symbols("lambda mu_star M L0", positive=True)
# mu*(M): L = ln(E_F/omega_c), omega_c ∝ M^{-1/2} → L = L0 + (1/2) ln M ; mu* = mu/(1+mu L)
mu_bare = sp.symbols("mu", positive=True)
L = L0 + sp.Rational(1, 2) * sp.log(M)
mu_star_M = mu_bare / (1 + mu_bare * L)
D = lam - mu_star_M * (1 + sp.Rational(62, 100) * lam)
lnTc = -sp.Rational(1, 2) * sp.log(M) + (-sp.Rational(104, 100) * (1 + lam) / D)  # ln omega_log + ln Phi
alpha = -sp.diff(lnTc, M) * M
alpha = sp.simplify(alpha.subs(mu_bare / (1 + mu_bare * L), mu))
# 用 mu* 表达
alpha_expr = sp.simplify(alpha.subs(L0, sp.solve(sp.Eq(mu_bare / (1 + mu_bare * (L0 + sp.log(M) / 2)), mu), L0)[0]))
target = sp.Rational(1, 2) * (1 - sp.Rational(104, 100) * (1 + lam) * (1 + sp.Rational(62, 100) * lam) * mu**2
                              / (lam - mu * (1 + sp.Rational(62, 100) * lam)) ** 2)
diff = sp.simplify(alpha_expr - target)
print("alpha - target =", diff)
for m in [0.10, 0.13, 0.20, 0.25, 0.30]:
    print(f"  lambda=2, mu*={m}: alpha = {float(target.subs({lam: 2, mu: m})):.3f}")

print("\n=== 3.2 非谐项：单模 SCHA 四阶模型 ===")
eta, x = sp.symbols("eta x", positive=True)
# x^2 = 1 + eta/x  →  x^3 - x - eta = 0
sol = sp.series(sp.solve(sp.Eq(x**3 - x - eta, 0), x)[0], eta, 0, 3) if False else None
# 直接用隐函数微分
F = x**3 - x - eta
dx_deta = -sp.diff(F, eta) / sp.diff(F, x)  # = 1/(3x^2-1)
print("dx/deta at x=1:", dx_deta.subs(x, 1))  # 期望 1/2
# 质量依赖：omega0 ∝ M^{-1/2}, eta ∝ M^{-1/2}
Mm = sp.symbols("M", positive=True)
eta_M = sp.symbols("eta0", positive=True) * Mm ** (-sp.Rational(1, 2))
# ln omega = ln omega0 + ln x ; -d ln omega / d ln M
dlnx_dlnM = (dx_deta / x) * sp.diff(eta_M, Mm) * Mm  # (1/x)(dx/deta)(deta/dlnM)
iso_exp = sp.Rational(1, 2) - dlnx_dlnM.subs(x, 1)
print("phonon isotope exponent (leading):", sp.simplify(iso_exp.subs(eta_M, eta)))  # 期望 1/2 + eta/4
# lambda ∝ 1/x^2 → d ln lambda / d ln M = -2 d ln x / d ln M
dlnlam_dlnM = sp.simplify((-2 * dlnx_dlnM).subs(x, 1).subs(eta_M, eta))
print("d ln lambda / d ln M (leading):", dlnlam_dlnM)  # 期望 eta/2
# alpha = (1/2 + eta/4) - (dlnTc/dlnlam) * (eta/2)
s = sp.symbols("s")  # s = ∂lnTc/∂lnλ
delta_alpha = sp.expand(eta / 4 - s * eta / 2)
print("delta alpha_anh =", delta_alpha, " → sign flips at s = 1/2")
# 数值：H3S 估计
eta_h3s = 1 / 0.70 - 1
print(f"H3S: 1/(1+eta)=0.70 → eta≈{eta_h3s:.2f}; s=0.75 → delta alpha ≈ {float(delta_alpha.subs({eta: eta_h3s, s: 0.75})):.3f}")

print("\n=== McMillan 的 ∂lnTc/∂lnλ（λ=2, μ*=0.13）===")
lnPhi = -sp.Rational(104, 100) * (1 + lam) / (lam - mu * (1 + sp.Rational(62, 100) * lam))
dlnTc_dlnlam = sp.simplify(sp.diff(lnPhi, lam) * lam)
print("  =", float(dlnTc_dlnlam.subs({lam: 2, mu: 0.13})))
