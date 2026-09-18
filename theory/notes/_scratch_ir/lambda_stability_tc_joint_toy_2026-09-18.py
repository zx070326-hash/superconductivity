"""
探针（numerical-evidence 级）：把 Semenok–Altshuler–Yuzbashyan 的电声耦合稳定性上界
lambda_*[P] 与 KAY 算子给出的 Tc 严格上界，在同一个谱形测度 P(dω) 上"联合自洽"求解，
看谱展宽 sigma 是否提升还是压低"联合天花板" Tc_max(sigma)。

KAY 侧只用 N=1 截断闭式 h(϶)=϶^2/(϶^2+1)（KAY Thm 4 eq.10），这是已证的 Tc 下界算子
（𝔨^(1) ≤ 𝔨），不是完整算子——本脚本只做 order-of-magnitude 的玩具模型，
不构成对完整 KAY 算子的证书级结论。

Semenok 稳定核 g(x) = 6x^2 + 12x^3 Im[psi'(ix)] + 6x^4 Re[psi''(ix)]
   (P3_tc_max_pressure.md §4, 从 arXiv:2407.12922 Eq. S12 逐字核对)
lambda_*[P] = 1 / max_T int g(ω/2πT) P(dω)

两点对称谱 P = (1/2)delta_{1-s} + (1/2)delta_{1+s}，⟨ω⟩=1 固定，sigma=s（因 <xi^2>=s^2/1^2=s^2），
gamma_3=0（对称）。

对每个 s：
 (1) 数值求 lambda_*(s) = 1/max_T [g((1-s)/2πT)/2 + g((1+s)/2πT)/2]（对 2πT 做一维数值最大化）
 (2) 用 N=1 KAY 核解自洽 Tc 方程： lambda_*(s) * [h((1-s)/2πTc)/2 + h((1+s)/2πTc)/2] = 1
 (3) 记录 Tc(s)/Tc(0)，检查 s 增大时联合天花板升还是降

同时和解析 sigma^2 阶微扰预测（正文推导：dlnTc/dsigma^2 = -dln϶/dsigma^2 = 1.1457，
 使用 P3 的 kappa_S=5.6035、N=1 KAY 的 kappa_1(϶_0)）做交叉核对。
"""
import mpmath as mp

mp.mp.dps = 30


def g_kernel(x):
    x = mp.mpf(x)
    ix = mp.mpc(0, 1) * x
    psi1 = mp.polygamma(1, ix)   # psi'(ix)
    psi2 = mp.polygamma(2, ix)   # psi''(ix)
    return 6*x**2 + 12*x**3*mp.im(psi1) + 6*x**4*mp.re(psi2)


def h1_kernel(x):
    x = mp.mpf(x)
    return x**2/(x**2+1)


def lambda_star_two_point(s, guess_2piT=3.0):
    """max_{2piT>0} avg_ω g(ω/(2piT)) for the two-point spectrum at mean=1, half-width s."""
    w1, w2 = mp.mpf(1)-s, mp.mpf(1)+s

    def negF(u):  # u = 2*pi*T
        if u <= 0:
            return mp.mpf('1e10')
        return -(g_kernel(w1/u) + g_kernel(w2/u))/2

    # golden-section search on u in (1e-3, 20)
    a, b = mp.mpf('1e-3'), mp.mpf(20)
    gr = (mp.sqrt(5)-1)/2
    c = b - gr*(b-a)
    d = a + gr*(b-a)
    for _ in range(200):
        if negF(c) < negF(d):
            b = d
        else:
            a = c
        c = b - gr*(b-a)
        d = a + gr*(b-a)
    u_star = (a+b)/2
    Fmax = -negF(u_star)
    return mp.mpf(1)/Fmax, u_star


def solve_tc_two_point(s, lam, guess_2piT=1.0):
    """solve lam * avg_ω h1(ω/(2piT)) = 1 for u=2piT, two-point spectrum mean=1 half-width s."""
    w1, w2 = mp.mpf(1)-s, mp.mpf(1)+s

    def F(u):
        return (h1_kernel(w1/u) + h1_kernel(w2/u))/2 - 1/lam

    # bisection: F decreasing in u (h1 decreasing in argument as u grows -> arg shrinks -> h1 shrinks)
    lo, hi = mp.mpf('1e-4'), mp.mpf(50)
    for _ in range(200):
        mid = (lo+hi)/2
        if F(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


print("=== s-scan: joint lambda_*(s) + N=1 Tc(s), <omega>=1 fixed, symmetric two-point ===")
results = []
for s in [mp.mpf('0.0001'), mp.mpf('0.02'), mp.mpf('0.05'), mp.mpf('0.10'), mp.mpf('0.15'),
          mp.mpf('0.20'), mp.mpf('0.30')]:
    lam_s, _ = lambda_star_two_point(s)
    u_tc = solve_tc_two_point(s, lam_s)
    results.append((s, lam_s, u_tc))
    print(f"s={float(s):.4f}  lambda_*={float(lam_s):.6f}  2piTc={float(u_tc):.6f}  "
          f"Tc/Tc(s=0)~will normalize below")

u0 = results[0][2]
print()
print("=== normalized (u=2*pi*Tc, so Tc(s)/Tc(0) = u(s)/u(0)) ===")
for s, lam_s, u_tc in results:
    Tc_ratio = u_tc/u0  # Tc = u/(2*pi), so Tc(s)/Tc(0) = u(s)/u(0)
    lam_ratio = lam_s/results[0][1]
    print(f"s={float(s):.4f}  lambda_*/lambda_*(0)={float(lam_ratio):.6f}  "
          f"Tc(s)/Tc(0)={float(Tc_ratio):.6f}")

print()
print("=== sanity: lambda_*(s=0.0001) should match Semenok Einstein value 3.6915 ===")
print(f"lambda_*(s~0) = {float(results[0][1]):.6f}  (expect 3.6915)")

print()
print("=== leading-order sigma^2 perturbative prediction check ===")
# kappa_S from repo P3 (rational-exact, p3_lambda_star_kappa.py): 5.6035
kappa_S = mp.mpf('5.6035')
lam0 = results[0][1]
wbar0 = mp.mpf(1)/u0  # varpi_0 = <omega>/(2 pi Tc0), <omega>=1
h_ = lambda x: x**2/(x**2+1)
hp = lambda x: 2*x/(x**2+1)**2
kappa1 = (2 - 6*wbar0**2)/(wbar0**2+1)**2
print(f"wbar0 = {float(wbar0):.6f}  (predicted ~0.6095)")
print(f"kappa_1(wbar0) = {float(kappa1):.6f}  (predicted ~ -0.1218)")
delta_over_sigma2 = -mp.mpf('0.5')*(kappa_S+kappa1) * h_(wbar0)/hp(wbar0)
print(f"predicted d(wbar)/d(sigma^2) at sigma=0  = {float(delta_over_sigma2):.6f}  (predicted ~ -1.1457)")
print(f"=> predicted d ln Tc / d sigma^2 at sigma=0 = {float(-delta_over_sigma2):.6f}")

# numeric estimate from s=0.05 point (sigma=s since <omega>=1)
s_test = results[2][0]
Tc_ratio_test = results[2][2]/u0
num_dlnTc_dsigma2 = mp.log(Tc_ratio_test) / s_test**2
print(f"numeric  d ln Tc / d sigma^2  (from s={float(s_test)} finite-diff) = {float(num_dlnTc_dsigma2):.6f}")
