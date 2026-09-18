"""
Weak-coupling asymptotics scan for γ-model: explore 𝔤(γ) as γ→∞.

Physical motivation (KAY-I Fig 1): γ=2 gives 𝔤(2)≈1.45, γ→∞ appears to approach
a finite limit 𝔤_∞≈1.3. This script scans γ=1,2,3,4,5,10,20,50,100 at fixed N
to extract:
  1. The asymptotic limit 𝔤_∞
  2. The decay law: 𝔤(γ) ≈ 𝔤_∞ + a/γ^p  (fit p, a)
  3. Physical interpretation of the γ→∞ limit

Method: reuses gamma_rayleigh_cert.py's certificate-level construction,
extracting only mu1_lo (certified lower bound) and mu1_hi (certified upper
bound) for each γ. All bounds are `certificate` level (mpmath.iv outward-
rounded interval arithmetic).

Repo convention: scratch script, stays uncommitted, lives in _scratch_ir/.
"""
import sys
from mpmath import iv, mp
import numpy as np

mp.dps = 30
iv.dps = 30


def lo(x):
    return mp.mpf(x.a) if hasattr(x, 'a') else mp.mpf(x)


def hi(x):
    return mp.mpf(x.b) if hasattr(x, 'b') else mp.mpf(x)


def build_G(gamma_val, N):
    """Build KAY's G(gamma) matrix (N+1 x N+1), all entries as iv.mpf."""
    g = iv.mpf(gamma_val)
    dim = N + 1

    def u(n):
        return 1 / iv.sqrt(iv.mpf(2 * n + 1))

    def powg(k):
        return iv.mpf(k) ** g

    # G1 diagonal: a[n] = (1/(2n+1)) * sum_{k=1}^n 2/k^gamma
    a = [iv.mpf(0)] * dim
    running = iv.mpf(0)
    for n in range(1, dim):
        running = running + 2 / powg(n)
        a[n] = running / iv.mpf(2 * n + 1)

    u_list = [u(n) for n in range(dim)]

    G = [[None] * dim for _ in range(dim)]
    for n in range(dim):
        un = u_list[n]
        for m in range(dim):
            if n == m:
                G[n][n] = -a[n] + un * un / powg(2 * n + 1)
            else:
                um = u_list[m]
                val = un * um * (1 / powg(abs(n - m)) + 1 / powg(n + m + 1))
                G[n][m] = val
    return dict(gamma=g, N=N, dim=dim, G=G, u=u_list)


def to_float_matrix(data):
    G, dim = data['G'], data['dim']
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i][j] = float(G[i][j].mid)
    return A


def rayleigh_lower_iv(G, dim, vfloat):
    """Certified LOWER bound for v^T G v / v^T v."""
    v = [iv.mpf(x) for x in vfloat]
    num = iv.mpf(0)
    for i in range(dim):
        Gv_i = iv.mpf(0)
        Grow = G[i]
        for j in range(dim):
            Gv_i = Gv_i + Grow[j] * v[j]
        num = num + v[i] * Gv_i
    den = iv.mpf(0)
    for i in range(dim):
        den = den + v[i] * v[i]
    return num / den


def gershgorin_mu1_hi(data):
    """Certified upper bound via Gershgorin (all off-diag >=0, no abs needed)."""
    G, dim = data['G'], data['dim']
    max_disc = None
    for n in range(dim):
        s = G[n][n]
        for m in range(dim):
            if m != n:
                s = s + G[n][m]
        disc_hi = hi(s)
        if max_disc is None or disc_hi > max_disc:
            max_disc = disc_hi
    return max_disc


def run_one_gamma(gamma_val, N, verbose=False):
    """Returns (mu1_lo, mu1_hi, lam1_float) for given gamma, N."""
    data = build_G(gamma_val, N)
    G, dim = data['G'], data['dim']

    A = to_float_matrix(data)
    evals, evecs = np.linalg.eigh(A)
    lam1_f = evals[-1]
    v1 = evecs[:, -1]

    mu1_lo_iv = rayleigh_lower_iv(G, dim, v1.tolist())
    mu1_lo = lo(mu1_lo_iv)

    mu1_hi = gershgorin_mu1_hi(data)

    if verbose:
        print(f"  γ={gamma_val:6.1f}: 𝔤 ∈ [{float(mu1_lo):.8f}, {float(mu1_hi):.8f}]  "
              f"width={float(mu1_hi-mu1_lo):.2e}  float64={lam1_f:.8f}")

    return mu1_lo, mu1_hi, lam1_f


def fit_decay(gammas, g_vals, g_inf_guess):
    """
    Fit 𝔤(γ) ≈ 𝔤_∞ + a/γ^p to data.
    Returns (g_inf, a, p) by least-squares on log-log residuals.

    Model: y := 𝔤(γ) - 𝔤_∞,  log(y) ≈ log(a) - p*log(γ)
    Use only γ>=3 for the fit (γ=1,2 may have larger corrections).
    """
    mask = np.array(gammas) >= 3
    x = np.array(gammas)[mask]
    y_raw = np.array(g_vals)[mask]

    # Try a range of g_inf candidates near g_inf_guess
    best_r2 = -np.inf
    best_params = None

    for g_inf_trial in np.linspace(g_inf_guess - 0.05, g_inf_guess + 0.05, 21):
        y = y_raw - g_inf_trial
        if np.any(y <= 0):
            continue
        log_y = np.log(y)
        log_x = np.log(x)
        # Linear fit: log_y = log_a - p*log_x
        A_mat = np.vstack([np.ones_like(log_x), -log_x]).T
        coeffs, residuals, _, _ = np.linalg.lstsq(A_mat, log_y, rcond=None)
        log_a, p = coeffs
        a = np.exp(log_a)

        # R^2
        ss_res = residuals[0] if len(residuals) > 0 else np.sum((log_y - (log_a - p*log_x))**2)
        ss_tot = np.sum((log_y - np.mean(log_y))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        if r2 > best_r2:
            best_r2 = r2
            best_params = (g_inf_trial, a, p, r2)

    return best_params


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    gammas = [1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 50.0, 100.0]

    print(f"γ-model weak-coupling asymptotics scan at N={N}")
    print("=" * 70)

    results = []
    for gamma_val in gammas:
        mu1_lo, mu1_hi, lam1_f = run_one_gamma(gamma_val, N, verbose=True)
        # Use midpoint as best estimate for fitting
        g_mid = float((mu1_lo + mu1_hi) / 2)
        results.append((gamma_val, mu1_lo, mu1_hi, g_mid, lam1_f))

    print()
    print("=" * 70)
    print("FITTING: 𝔤(γ) ≈ 𝔤_∞ + a/γ^p")
    print("=" * 70)

    g_vals = [r[3] for r in results]
    g_inf_guess = g_vals[-1]  # γ=100 value as initial guess

    g_inf, a, p, r2 = fit_decay(gammas, g_vals, g_inf_guess)

    print(f"  𝔤_∞ = {g_inf:.6f}")
    print(f"  a   = {a:.6f}")
    print(f"  p   = {p:.6f}")
    print(f"  R²  = {r2:.6f}")
    print()
    print("Comparison (γ>=3):")
    for gamma_val, mu1_lo, mu1_hi, g_mid, _ in results:
        if gamma_val >= 3:
            g_fit = g_inf + a / gamma_val**p
            error = (g_mid - g_fit) / g_mid * 100
            print(f"  γ={gamma_val:6.1f}: 𝔤_data={g_mid:.6f}  𝔤_fit={g_fit:.6f}  error={error:+.2f}%")

    print()
    print("=" * 70)
    print("CERTIFICATE-LEVEL BOUNDS")
    print("=" * 70)
    for gamma_val, mu1_lo, mu1_hi, g_mid, lam1_f in results:
        width_rel = float((mu1_hi - mu1_lo) / mu1_lo) * 100
        print(f"  γ={gamma_val:6.1f}: 𝔤(γ) ∈ [{float(mu1_lo):.8f}, {float(mu1_hi):.8f}]  "
              f"rel_width={width_rel:.3f}%")
