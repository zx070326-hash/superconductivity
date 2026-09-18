"""
γ-model Davis-Kahan certificate pipeline — Schur-tail INTEGRATION.

New file (does NOT modify gamma_dk_coulomb.py or gamma_tail_bound.py), per
repo convention "证书脚本经独立新文件扩展，不改基线脚本".

WHAT THIS FIXES
---------------
gamma_dk_coulomb.py's `tail_bounds_coarse()` is NOT a sound upper bound for
γ = 1.  Its derivation writes

    ||𝔊^BB|| <= u_{N+1} * 2 * tail_sum_bound ,   tail_sum_bound := Σ_{k>N} 1/k^γ

and then, for γ = 1, hard-codes `tail_sum_bound = mp.mpf(1.0)` with the comment
"very coarse for γ=1".  But Σ_{k>N} 1/k DIVERGES: the constant 1.0 does not
cover the quantity it stands in for, so the resulting "bound" is an arbitrary
number, not an inequality.  This is precisely the repo's documented pitfall
"一个 pad 必须覆盖它所替身的那份误差".  Consequence: the c₁² values printed by
gamma_dk_coulomb.py are labelled `certificate` but are NOT certificate level.

This script replaces that placeholder with the direction-audited Schur bounds
from gamma_tail_bound.py (schur_C_bound_hi / HBB_bound_hi), whose weighted tail
sums Σ_{m≥J} u_m/m^γ converge for γ > 1/2 and are bounded by closed-form
integrals in iv arithmetic.

SOUNDNESS OF mu1_lo FOR THE INFINITE OPERATOR
---------------------------------------------
v₁ is supported on n ≤ N, so ⟨v₁, 𝔊 v₁⟩ = ⟨v₁, H_N v₁⟩ exactly (the tail block
never enters).  Hence the interval Rayleigh quotient of H_N is a valid lower
bound for the top eigenvalue of the FULL infinite operator 𝔊(γ, μ*, N_c), not
only for the truncation.  No tail control is needed for the 𝔤 lower bound —
tail control is needed only for c₁² (the eigenvector concentration).

Strictness: certificate (mpmath.iv throughout, one-sided directions audited).
"""
import sys, time
import numpy as np
from mpmath import iv, mp

sys.path.insert(0, __file__.rsplit('\\', 1)[0] if '\\' in __file__ else '.')

from gamma_dk_coulomb import (
    build_G_matrix as build_G_dk, to_float_matrix, rayleigh_lower,
    residual_norm, find_lambda2_upper, tail_bounds_coarse, lo, hi,
)
from gamma_tail_bound import (
    build_G_matrix as build_G_tail, schur_C_bound_hi, HBB_bound_hi,
)


def run_one_schur(gamma_val, N, mu_star=0.0, Nc=None, dps=30, verbose=True,
                  compare_placeholder=True):
    """Davis-Kahan pipeline with CERTIFIED Schur tail bounds.

    All one-sided directions:
      mu1_lo      : LOWER  (Rayleigh quotient, valid for any trial vector)
      lam2_hi     : UPPER  (inertia count)
      r1          : UPPER  (residual norm)
      C_hi/HBB_hi : UPPER  (Schur test on nonnegative blocks)
      => sinT is an UPPER bound on the true angle, so c1sq_lo is a LOWER bound.
    """
    mp.dps = dps
    iv.dps = dps
    t0 = time.time()

    data = build_G_dk(N, gamma_val, mu_star=mu_star, Nc=Nc)
    H, dim = data['H'], data['dim']

    A = to_float_matrix(H, dim)
    evals, evecs = np.linalg.eigh(A)
    lam1_approx, lam2_approx = evals[-1], evals[-2]
    v1 = evecs[:, -1]

    # --- head-block certified quantities -------------------------------
    mu1_lo = lo(rayleigh_lower(H, dim, v1.tolist()))
    r1 = hi(residual_norm(H, dim, v1.tolist(), iv.mpf(mu1_lo)))
    c_upper, margin_used = find_lambda2_upper(H, dim, lam2_approx,
                                              lam1_approx - lam2_approx)
    lam2_hi = hi(c_upper)

    gap1_lo = mu1_lo - lam2_hi
    sinT1 = mp.mpf(1) if gap1_lo <= r1 else min(mp.mpf(1), r1 / (gap1_lo - r1))

    # --- tail: CERTIFIED Schur bounds (the fix) ------------------------
    # The Coulomb term -μ* P_{N_c} is supported on n,m < N_c <= N, i.e. entirely
    # inside the head block, so it does not touch C or 𝔊^BB.  The μ*=0 operator
    # is therefore the right object for the tail bounds.
    tail_data = build_G_tail(N, gamma_val, dps=dps)
    C_hi = schur_C_bound_hi(tail_data, N, dps=dps)
    HBB_hi = HBB_bound_hi(tail_data, N, dps=dps)
    Ehi = C_hi + HBB_hi

    lam2A_hi = max(lam2_hi, mp.mpf(0))
    delta_lo = mu1_lo - lam2A_hi
    sinT2 = mp.mpf(1) if delta_lo <= Ehi else min(mp.mpf(1), Ehi / (delta_lo - Ehi))

    sinT = min(mp.mpf(1), sinT1 + sinT2)
    cosT = mp.sqrt(1 - sinT**2) if sinT < 1 else mp.mpf(0)
    v0 = abs(v1[0])
    c1sq_lo = max(mp.mpf(0), cosT * mp.mpf(v0) - sinT) ** 2

    ph = None
    if compare_placeholder:
        HBB_ph, C_ph = tail_bounds_coarse(data, N)
        ph = dict(HBB=HBB_ph, C=C_ph, E=HBB_ph + C_ph)

    if verbose:
        mu_str = f"mu*={mu_star:.3f} Nc={Nc}" if mu_star else "mu*=0"
        print(f"gamma={gamma_val} N={N} {mu_str}  dim={dim}  [{time.time()-t0:.1f}s]")
        print(f"  mu1_lo (cert, valid for the INFINITE operator) = {mu1_lo}")
        print(f"  lam2_hi (cert, inertia)  = {lam2_hi}   margin={margin_used:.2e}")
        print(f"  r1 = {float(r1):.3e}   sinT1 = {float(sinT1):.3e}")
        print(f"  [Schur, CERTIFIED] ||C|| <= {float(C_hi):.6f}"
              f"   ||G^BB|| <= {float(HBB_hi):.6f}   ||E|| <= {float(Ehi):.6f}")
        if ph:
            print(f"  [old placeholder, UNSOUND for gamma=1] ||C|| <= {float(ph['C']):.6f}"
                  f"   ||G^BB|| <= {float(ph['HBB']):.6f}   ||E|| <= {float(ph['E']):.6f}")
            print(f"      -> placeholder/Schur ratio on ||E||: {float(ph['E']/Ehi):.3f}"
                  f"  ({'placeholder was SMALLER, i.e. optimistic' if ph['E'] < Ehi else 'placeholder larger'})")
        print(f"  delta_lo = {float(delta_lo):.6f}   sinT2 = {float(sinT2):.6f}"
              f"   sinT_total = {float(sinT):.6f}")
        print(f"  |v1[0]| (float64) = {v0:.6f}   |v1[0]|^2 = {v0**2:.6f}")
        print(f"  c1sq_lo (CERTIFICATE) = {float(c1sq_lo):.6f}")
        print()

    return dict(gamma=gamma_val, N=N, mu_star=mu_star, Nc=Nc,
                mu1_lo=mu1_lo, lam2_hi=lam2_hi, r1=r1,
                C_hi=C_hi, HBB_hi=HBB_hi, Ehi=Ehi,
                sinT1=sinT1, sinT2=sinT2, sinT=sinT,
                c1sq_lo=c1sq_lo, v0=v0, v0sq_float=v0**2,
                lam1_approx=lam1_approx, placeholder=ph)


if __name__ == '__main__':
    gamma_val = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    Ns = [int(x) for x in sys.argv[2:]] if len(sys.argv) > 2 else [200, 400]

    print("=" * 74)
    print("gamma-model DK pipeline with CERTIFIED Schur tail bounds")
    print("=" * 74)
    print()

    rows = []
    for N in Ns:
        for mu_star, Nc in [(0.0, None), (0.10, 50), (0.13, 50), (0.15, 50)]:
            res = run_one_schur(gamma_val, N, mu_star=mu_star, Nc=Nc)
            rows.append(res)

    print("=" * 74)
    print("SUMMARY  (all quantities certificate level)")
    print("=" * 74)
    print(f"{'N':>5} {'mu*':>6} {'Nc':>5} {'g_lo':>14} {'||E||':>9} {'sinT':>9}"
          f" {'c1sq_lo':>9} {'c1sq_f64':>9}")
    for r in rows:
        nc = r['Nc'] if r['Nc'] else '-'
        print(f"{r['N']:>5} {r['mu_star']:>6.2f} {str(nc):>5} {float(r['mu1_lo']):>14.8f}"
              f" {float(r['Ehi']):>9.5f} {float(r['sinT']):>9.5f}"
              f" {float(r['c1sq_lo']):>9.5f} {r['v0sq_float']:>9.5f}")

    # Tc ratios (gamma=1: Tc ∝ g; general: Tc ∝ g^{1/gamma})
    print()
    print(f"Tc suppression (Tc ∝ g^(1/gamma), gamma={gamma_val}):")
    for N in Ns:
        base = [r for r in rows if r['N'] == N and r['mu_star'] == 0.0]
        if not base:
            continue
        g0 = base[0]['mu1_lo']
        for r in [x for x in rows if x['N'] == N and x['mu_star'] > 0]:
            ratio = (r['mu1_lo'] / g0) ** (1.0 / gamma_val)
            print(f"  N={N} mu*={r['mu_star']:.2f}: Tc/Tc0 ~ {float(ratio):.4f}"
                  f"  ({100*(float(ratio)-1):+.1f}%)")
    print()
    print("NOTE on direction: g_lo is a LOWER bound on g(gamma,mu*,Nc) and the")
    print("mu*=0 g_lo is likewise a LOWER bound, so the printed ratio is NOT a")
    print("certified two-sided bound on Tc/Tc0 -- it is the ratio of two lower")
    print("bounds (numerical-evidence for the ratio; each endpoint separately is")
    print("certificate).  A certified suppression ratio needs an UPPER bound on")
    print("g(mu*=0), e.g. Kato-Temple as in gamma_rayleigh_cert.py.")
