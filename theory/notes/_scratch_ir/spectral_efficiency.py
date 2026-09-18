#!/usr/bin/env python3
r"""
Spectral efficiency eta(P) and pairing convex optimization computations.

  Part 1: eta_AD (Allen-Dynes level) for hydrides and conventional SC
  Part 2: KAY-level cross-check for 2-atom spectra (small varpi)
  Part 3: Dual-moment (C2, C4) extremal scan with 3-atom spectra (AD level)
  Part 4: mu*>0 concavity check: H_mu''(x) numerical

All results: numerical-evidence level.  Not a gate script.
Convention: scratch, not committed.
"""
import math
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import eigvalsh

MEV_TO_K = 11.604518

# ====================================================================
# KAY operator assembly -- vectorized for speed
# ====================================================================

def build_H_fast(varpi, N):
    """Build N x N Einstein operator H(varpi) = O[bracket(varpi)]."""
    x = varpi ** 2
    u = 1.0 / np.sqrt(2 * np.arange(N) + 1.0)
    ks = np.arange(2 * N + 1, dtype=float)
    fk = x / (ks ** 2 + x)
    fk[0] = 1.0

    nn = np.arange(N)
    mm = np.arange(N)
    diff_nm = np.abs(nn[:, None] - mm[None, :])
    sum_nm = nn[:, None] + mm[None, :] + 1

    L2 = fk[diff_nm] * u[:, None] * u[None, :]
    L3 = fk[sum_nm] * u[:, None] * u[None, :]

    # Diagonal L1
    L1_diag = np.zeros(N)
    for n in range(1, N):
        L1_diag[n] = np.sum(2.0 * fk[1:n + 1]) / (2 * n + 1)

    H = L2 + L3
    np.fill_diagonal(H, 0)
    H += np.diag(fk[2 * nn + 1] * u ** 2 - L1_diag)

    # Fix diagonal: remove self-contribution from L2
    for n in range(N):
        H[n, n] = -L1_diag[n] + fk[2 * n + 1] * u[n] ** 2

    return H


def build_H_fast2(varpi, N):
    """Vectorized build of H(varpi). Correct assembly per P4 notation."""
    x = varpi ** 2
    u = 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)
    ks = np.arange(2 * N + 1, dtype=float)
    fk = x / (ks ** 2 + x)
    fk[0] = 1.0

    nn = np.arange(N)
    mm = np.arange(N)

    # Off-diagonal: H_{nm} = (fk[|n-m|] + fk[n+m+1]) * u_n * u_m  for n != m
    diff_nm = np.abs(nn[:, None] - mm[None, :])
    sum_nm = nn[:, None] + mm[None, :] + 1

    H = (fk[diff_nm] + fk[sum_nm]) * (u[:, None] * u[None, :])

    # Diagonal: H_{nn} = -L1_n + fk[2n+1] * u_n^2
    # L1_n = (1/(2n+1)) * sum_{k=1}^{n} 2*fk[k]
    cumfk = np.cumsum(2.0 * fk[1:N + 1])
    L1 = np.zeros(N)
    L1[1:] = cumfk[:N - 1] / (2 * nn[1:] + 1)

    # The diagonal from the "off-diagonal" formula at n=m gives fk[0]+fk[2n+1] = 1+fk[2n+1],
    # but the correct diagonal is -L1_n + fk[2n+1]*u_n^2. Need to overwrite.
    diag_correct = -L1 + fk[2 * nn + 1] * u ** 2
    np.fill_diagonal(H, diag_correct)

    return H


def h_einstein(varpi, N=None):
    """Top eigenvalue of H(varpi)."""
    if varpi < 1e-10:
        return 0.0
    if N is None:
        N = max(100, int(5 * varpi))
    H = build_H_fast2(varpi, N)
    return float(eigvalsh(H)[-1])


def build_K_fast2(weights_varpi, N):
    """Build K for multi-atom spectrum: K = O[F_T], F_T = sum_i w_i * bracket(varpi_i)."""
    nn = np.arange(N)
    mm = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * nn + 1.0)

    # Aggregate F_T(k) for k = 0..2N
    ks = np.arange(2 * N + 1, dtype=float)
    fk = np.zeros(2 * N + 1)
    for w, vp in weights_varpi:
        x = vp ** 2
        fk += w * x / (ks ** 2 + x)
    fk[0] = sum(w for w, _ in weights_varpi)

    diff_nm = np.abs(nn[:, None] - mm[None, :])
    sum_nm = nn[:, None] + mm[None, :] + 1

    K = (fk[diff_nm] + fk[sum_nm]) * (u[:, None] * u[None, :])

    cumfk = np.cumsum(2.0 * fk[1:N + 1])
    L1 = np.zeros(N)
    L1[1:] = cumfk[:N - 1] / (2 * nn[1:] + 1)
    diag_correct = -L1 + fk[2 * nn + 1] * u ** 2
    np.fill_diagonal(K, diag_correct)
    return K


def Tc_KAY_einstein(lam, omega_E_K, N=None):
    """Tc from KAY for Einstein spectrum."""
    def obj(T_K):
        varpi = omega_E_K / (2 * math.pi * T_K)
        return lam * h_einstein(varpi, N) - 1.0
    T_lo, T_hi = 0.1, omega_E_K * 5
    try:
        if obj(T_lo) < 0: return 0.0
        if obj(T_hi) > 0: T_hi *= 10
        return brentq(obj, T_lo, T_hi, xtol=1e-3)
    except: return 0.0


def Tc_KAY_general(lam, atoms_K, N=None):
    """Tc from KAY for multi-atom spectrum (atoms_K: [(w, omega_K), ...])."""
    def obj(T_K):
        wv = [(w, omega / (2 * math.pi * T_K)) for w, omega in atoms_K]
        max_vp = max(vp for _, vp in wv)
        Nuse = N or max(100, int(5 * max_vp))
        K = build_K_fast2(wv, Nuse)
        return lam * float(eigvalsh(K)[-1]) - 1.0
    max_om = max(om for _, om in atoms_K)
    T_lo, T_hi = 0.1, max_om * 5
    try:
        if obj(T_lo) < 0: return 0.0
        if obj(T_hi) > 0: T_hi *= 10
        return brentq(obj, T_lo, T_hi, xtol=1e-3)
    except: return 0.0


# ====================================================================
# Allen-Dynes formula
# ====================================================================

def f1_AD(lam, mustar):
    return (1 + (lam / (2.46 * (1 + 3.8 * mustar))) ** 1.5) ** (1.0 / 3)

def f2_AD(r, lam, mustar):
    if abs(r - 1) < 1e-15: return 1.0
    return 1 + (r - 1) * lam ** 2 / (lam ** 2 + (1.82 * (1 + 6.3 * mustar) * r) ** 2)

def Tc_AD(lam, omega_log_K, mustar, r=1.0):
    denom = lam - mustar * (1 + 0.62 * lam)
    if denom <= 0: return 0.0
    expo = -1.04 * (1 + lam) / denom
    return (f1_AD(lam, mustar) * f2_AD(r, lam, mustar) * omega_log_K / 1.20) * math.exp(expo)

def Tc_Einstein_AD(lam, omega_E_K, mustar):
    return Tc_AD(lam, omega_E_K, mustar, r=1.0)

def eta_AD(r, lam, mustar):
    """eta_AD = Tc(P; lambda) / Tc^E(lambda, omega_rms).
    At AD level: = f2(r)/r  (the omega_log factor cancels).
    """
    return f2_AD(r, lam, mustar) / r


# ====================================================================
# PART 1
# ====================================================================
print("=" * 95)
print("PART 1: Spectral efficiency eta(P) = Tc(P) / Tc^E(lambda, sqrt<omega^2>)")
print("  At Allen-Dynes level: eta_AD = f2(r,lam,mu*) / r,  r = omega_2/omega_log")
print("  All values: numerical-evidence")
print("=" * 95)

lah10 = [
    dict(P=129, lam=3.62, wlog_K=886.6, Tc_ME=255.3, r=1.194, mustar=0.1),
    dict(P=163, lam=2.67, wlog_K=1118.7, Tc_ME=242.8, r=1.075, mustar=0.1),
    dict(P=214, lam=2.06, wlog_K=1340.3, Tc_ME=237.9, r=1.105, mustar=0.1),
    dict(P=264, lam=1.73, wlog_K=1469.1, Tc_ME=216.9, r=1.095, mustar=0.1),
]

h3s_wlog_K = 92.86 * MEV_TO_K
h3s = [dict(P=200, lam=1.84, wlog_K=h3s_wlog_K, Tc_ME=194.0, r=1.20, mustar=0.16)]

nb_data = [dict(name="Nb (BETE-NET)", lam=1.2325, wlog_K=14.55*MEV_TO_K, r=1.14, mustar=0.13)]
mgb2_data = [dict(name="MgB2 (BETE-NET)", lam=0.6235, wlog_K=65.45*MEV_TO_K, r=1.12, mustar=0.13)]

print()
print("--- LaH10 (Errea 2020, r from P3 closed-form inversion, mu*=0.1) ---")
print(f"{'P(GPa)':>8} {'lam':>6} {'wlog_K':>8} {'r':>6} {'f2':>7} {'eta_AD':>8}"
      f" {'Tc^ME':>7} {'TcE_AD':>8} {'TcP_AD':>8}")
for d in lah10:
    f2v = f2_AD(d['r'], d['lam'], d['mustar'])
    ev = eta_AD(d['r'], d['lam'], d['mustar'])
    om2 = d['wlog_K'] * d['r']
    TcE = Tc_Einstein_AD(d['lam'], om2, d['mustar'])
    TcP = Tc_AD(d['lam'], d['wlog_K'], d['mustar'], d['r'])
    print(f"{d['P']:>8} {d['lam']:>6.2f} {d['wlog_K']:>8.1f} {d['r']:>6.3f} "
          f"{f2v:>7.4f} {ev:>8.5f} {d['Tc_ME']:>7.1f} {TcE:>8.1f} {TcP:>8.1f}")

print()
print("--- H3S (Errea 2015, 200 GPa theory, mu*=0.16; r estimated ~1.20) ---")
for d in h3s:
    f2v = f2_AD(d['r'], d['lam'], d['mustar'])
    ev = eta_AD(d['r'], d['lam'], d['mustar'])
    om2 = d['wlog_K'] * d['r']
    TcE = Tc_Einstein_AD(d['lam'], om2, d['mustar'])
    TcP = Tc_AD(d['lam'], d['wlog_K'], d['mustar'], d['r'])
    print(f"  P={d['P']} GPa  lam={d['lam']:.2f}  wlog={d['wlog_K']:.1f} K  r~{d['r']:.2f}  "
          f"f2={f2v:.4f}  eta_AD={ev:.5f}  TcE={TcE:.1f} K  TcP_AD={TcP:.1f} K")

print()
print("--- Non-hydrides ---")
for d in nb_data + mgb2_data:
    f2v = f2_AD(d['r'], d['lam'], d['mustar'])
    ev = eta_AD(d['r'], d['lam'], d['mustar'])
    om2 = d['wlog_K'] * d['r']
    TcE = Tc_Einstein_AD(d['lam'], om2, d['mustar'])
    TcP = Tc_AD(d['lam'], d['wlog_K'], d['mustar'], d['r'])
    print(f"  {d['name']:20s}  lam={d['lam']:.4f}  wlog={d['wlog_K']:.1f} K  r~{d['r']:.2f}  "
          f"f2={f2v:.4f}  eta_AD={ev:.5f}  TcE={TcE:.1f} K  TcP_AD={TcP:.1f} K")

print()
print("--- eta_AD as function of r (sensitivity) ---")
print(f"{'r':>6} {'lam=2.5,mu*=0.10':>16} {'lam=1.8,mu*=0.16':>17} {'lam=1.0,mu*=0.13':>17}")
for r in [1.00, 1.02, 1.05, 1.08, 1.10, 1.15, 1.20, 1.30, 1.50, 2.00]:
    print(f"{r:>6.2f} {eta_AD(r, 2.5, 0.10):>16.5f} {eta_AD(r, 1.8, 0.16):>17.5f} "
          f"{eta_AD(r, 1.0, 0.13):>17.5f}")


# ====================================================================
# PART 2: KAY cross-check (moderate varpi, small N)
# ====================================================================
print()
print("=" * 95)
print("PART 2: KAY operator cross-check -- 2-atom vs Einstein")
print("=" * 95)

N_KAY = 150  # sufficient for varpi < 10

print()
print("h(varpi) sanity:")
for vp in [0.5, 1.0, 2.0, 3.0, 5.0]:
    print(f"  varpi={vp:.1f}  h={h_einstein(vp, N_KAY):.8f}")

lam_t = 2.0
omega_E = 500.0  # K, small enough to keep varpi < 10
TcE_k = Tc_KAY_einstein(lam_t, omega_E, N_KAY)
TcE_a = Tc_Einstein_AD(lam_t, omega_E, 0.1)
print(f"\nlambda={lam_t}, omega_E={omega_E} K:")
print(f"  Tc^E(KAY) = {TcE_k:.2f} K,  Tc^E(AD,mu*=0.1) = {TcE_a:.2f} K")

print(f"\n{'eps':>6} {'om1':>7} {'om2':>7} {'r':>6} {'TcP_KAY':>9} {'eta_KAY':>9} {'eta_AD':>8}")
for eps, sf in [(0.1, 1.5), (0.2, 1.5), (0.3, 2.0), (0.5, 2.0), (0.1, 3.0)]:
    om2 = sf * omega_E
    om1sq = (omega_E**2 - eps * om2**2) / (1 - eps)
    if om1sq <= 0: continue
    om1 = math.sqrt(om1sq)
    atoms = [(1 - eps, om1), (eps, om2)]
    TcP = Tc_KAY_general(lam_t, atoms, N_KAY)
    eta_k = TcP / TcE_k if TcE_k > 0 else float('nan')
    wlog = math.exp((1-eps)*math.log(om1) + eps*math.log(om2))
    rv = omega_E / wlog
    eta_a = eta_AD(rv, lam_t, 0.1)
    print(f"{eps:>6.2f} {om1:>7.1f} {om2:>7.1f} {rv:>6.3f} {TcP:>9.2f} {eta_k:>9.5f} {eta_a:>8.5f}")


# ====================================================================
# PART 3: Dual-moment extremal scan (AD level)
# ====================================================================
print()
print("=" * 95)
print("PART 3: Dual-moment (C2, C4) extremal -- 3-atom optimization at AD level")
print("=" * 95)

C2_ref = (1.075 * 1118.7) ** 2
omegarms = math.sqrt(C2_ref)
lam_ref = 2.67
mu_ref = 0.1
TcE_ref = Tc_Einstein_AD(lam_ref, omegarms, mu_ref)
print(f"\nRef: LaH10 163 GPa, lambda={lam_ref}, mu*={mu_ref}")
print(f"  omega_rms={omegarms:.1f} K, Tc^E(AD)={TcE_ref:.1f} K")

def solve_w3(om1, om2, om3, C2, C4):
    A = np.array([[1,1,1],[om1**2,om2**2,om3**2],[om1**4,om2**4,om3**4]])
    b = np.array([1, C2, C4])
    try: w = np.linalg.solve(A, b)
    except: return None
    if np.any(w < -1e-10): return None
    w = np.maximum(w, 0); w /= w.sum()
    return w

def Tc_3a_AD(om1, om2, om3, w1, w2, w3, lam, mustar):
    wlog = math.exp(w1*math.log(om1) + w2*math.log(om2) + w3*math.log(om3))
    om2_rms = math.sqrt(w1*om1**2 + w2*om2**2 + w3*om3**2)
    r = om2_rms / wlog
    return Tc_AD(lam, wlog, mustar, r)

print()
print(f"{'C4/C2^2':>8} {'Tc_best':>8} {'Tc^E':>7} {'eta_opt':>8} "
      f"{'om1':>8} {'om2':>8} {'om3':>8} {'w1':>6} {'w2':>6} {'w3':>6}")

for ratio in [1.001, 1.01, 1.05, 1.10, 1.15, 1.20, 1.30, 1.50]:
    C4 = ratio * C2_ref ** 2
    best_Tc = 0.0; best_cfg = None

    for o1f in np.linspace(0.1, 0.95, 25):
        om1 = o1f * omegarms
        for o3f in np.linspace(1.05, 4.0, 25):
            om3 = o3f * omegarms
            for o2f in np.linspace(0.5, 2.0, 12):
                om2 = o2f * omegarms
                if abs(om1-om2)<1 or abs(om2-om3)<1 or abs(om1-om3)<1: continue
                w = solve_w3(om1, om2, om3, C2_ref, C4)
                if w is None: continue
                Tc = Tc_3a_AD(om1, om2, om3, w[0], w[1], w[2], lam_ref, mu_ref)
                if Tc > best_Tc:
                    best_Tc = Tc; best_cfg = (om1, om2, om3, w[0], w[1], w[2])

    if best_cfg:
        o1,o2,o3,w1,w2,w3 = best_cfg
        print(f"{ratio:>8.3f} {best_Tc:>8.1f} {TcE_ref:>7.1f} {best_Tc/TcE_ref:>8.5f} "
              f"{o1:>8.1f} {o2:>8.1f} {o3:>8.1f} {w1:>6.3f} {w2:>6.3f} {w3:>6.3f}")
    else:
        print(f"{ratio:>8.3f}  (no feasible config)")


# ====================================================================
# PART 4: mu*>0 concavity
# ====================================================================
print()
print("=" * 95)
print("PART 4: mu*>0 concavity of H_mu(x) = h_mu(sqrt(x))")
print("  M = H(varpi) - 2*mu_tilde * |u><u|,  h_mu = lambda_1(M)")
print("=" * 95)

N_MU = 200

def h_mu_val(varpi, mu_tilde, N=N_MU):
    if varpi < 1e-10: return 0.0
    Nuse = max(N, int(5 * varpi))
    H = build_H_fast2(varpi, Nuse)
    u = 1.0 / np.sqrt(2.0 * np.arange(Nuse) + 1.0)
    M = H - 2 * mu_tilde * np.outer(u, u)
    return float(eigvalsh(M)[-1])

def H_mu_pp(x, mu_tilde, N=N_MU):
    dx = max(x * 5e-3, 1e-3)
    vals = [h_mu_val(math.sqrt(max(x + i*dx, 1e-8)), mu_tilde, N) for i in range(-2, 3)]
    return (-vals[0] + 16*vals[1] - 30*vals[2] + 16*vals[3] - vals[4]) / (12 * dx**2)

print()
print("(a) lambda_1(M) convexity in P: M = O[F_T] - 2*mu_tilde*|u><u|")
print("  O[F_T] is linear in P, |u><u| is P-independent => M is affine in P")
print("  => lambda_1(M) is convex in P (sup of affine). [analytic]")

print()
print("(b) H_mu''(x) evaluation:")
print(f"{'x':>8} {'mu=0':>12} {'mu=0.03':>12} {'mu=0.05':>12} {'mu=0.10':>12} {'mu=0.13':>12}")

x_test = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0, 100.0]
for x in x_test:
    vals = [H_mu_pp(x, mt, N_MU) for mt in [0.0, 0.03, 0.05, 0.10, 0.13]]
    print(f"{x:>8.1f} " + " ".join(f"{v:>12.4e}" for v in vals))

print()
print("(c) Concavity scan -- any H_mu''(x) > 0?")
any_pos = False
for mt in [0.03, 0.05, 0.10, 0.13, 0.20, 0.26]:
    for x in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0]:
        v = H_mu_pp(x, mt, N_MU)
        if v > 0:
            print(f"  POSITIVE: x={x}, mu_tilde={mt}, H_mu''={v:.4e}")
            any_pos = True
if not any_pos:
    print("  All H_mu''(x) < 0 in scanned range => concavity holds [numerical-evidence]")

print()
print("(d) Direct 2-atom vs Einstein test at mu*>0:")
for mt in [0.05, 0.10, 0.15]:
    vp_rms = 3.0
    h_ein = h_mu_val(vp_rms, mt, N_MU)
    # 2-atom: 0.5*delta(1) + 0.5*delta(sqrt(17))  => <varpi^2>=9
    wv = [(0.5, 1.0), (0.5, math.sqrt(17))]
    Nuse = max(N_MU, int(5 * math.sqrt(17)))
    K = build_K_fast2(wv, Nuse)
    u = 1.0 / np.sqrt(2.0 * np.arange(Nuse) + 1.0)
    M = K - 2 * mt * np.outer(u, u)
    kap = float(eigvalsh(M)[-1])
    print(f"  mu_tilde={mt:.2f}: h_mu(3.0)={h_ein:.6f}  kappa_mu(2-atom)={kap:.6f}  "
          f"ratio={kap/h_ein:.6f}  (<1 = Einstein wins)")

print()
print("=" * 95)
print("DONE. All results: numerical-evidence.")
print("=" * 95)
