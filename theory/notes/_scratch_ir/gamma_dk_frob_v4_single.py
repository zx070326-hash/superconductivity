"""
gamma_dk_frob_v4_single.py — Test v4 for a single γ value.
"""
import sys, os, time
import numpy as np
from mpmath import iv, mp

os.environ['OMP_NUM_THREADS'] = '1'

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from gamma_dk_coulomb import (
    build_G_matrix, to_float_matrix, rayleigh_lower,
    residual_norm, find_lambda2_upper, lo, hi,
)
from gamma_tail_bound_v2 import C_bound_weighted_hi, lamB_bound_hi
from gamma_frobenius_C import frob_C_bound_hi


def run_one_v4(gamma_val, N, mu_star=0.0, Nc=None, dps=30, verbose=True):
    """v4: use min(Schur, Frobenius) for ||C||."""
    mp.dps = dps
    iv.dps = dps
    t0 = time.time()

    print(f"[run_one_v4] Starting gamma={gamma_val}, N={N}, mu*={mu_star}")
    sys.stdout.flush()

    data = build_G_matrix(N, gamma_val, mu_star=mu_star, Nc=Nc)
    H, dim = data['H'], data['dim']

    A = to_float_matrix(H, dim)
    evals, evecs = np.linalg.eigh(A)
    lam1_approx, lam2_approx = evals[-1], evals[-2]
    v1 = evecs[:, -1]

    print(f"[run_one_v4] Rayleigh...")
    sys.stdout.flush()
    mu1_lo = lo(rayleigh_lower(H, dim, v1.tolist()))
    r1 = hi(residual_norm(H, dim, v1.tolist(), iv.mpf(mu1_lo)))
    c_upper, margin_used = find_lambda2_upper(H, dim, lam2_approx,
                                              lam1_approx - lam2_approx)
    lam2_hi = hi(c_upper)

    gap1_lo = mu1_lo - lam2_hi
    sinT1 = mp.mpf(1) if gap1_lo <= r1 else min(mp.mpf(1), r1 / (gap1_lo - r1))

    # --- tail bounds: weighted Schur
    print(f"[run_one_v4] Weighted Schur...")
    sys.stdout.flush()
    C_schur, A_sup, B_sup = C_bound_weighted_hi(N, gamma_val, dps=dps)

    # --- Frobenius bound
    print(f"[run_one_v4] Frobenius (this will take ~35-70s)...")
    sys.stdout.flush()
    C_frob, finite, tail, M_used = frob_C_bound_hi(N, gamma_val, M=None, dps=dps)

    # --- take minimum (both are upper bounds, so min is sound)
    C_hi = min(C_schur, C_frob)

    # --- block-diagonal DK (v3 structure)
    print(f"[run_one_v4] Block-diagonal DK...")
    sys.stdout.flush()
    lamB_hi, _parts = lamB_bound_hi(N, gamma_val, dps=dps)
    lam2A3_hi = max(lam2_hi, lamB_hi)
    delta3_lo = mu1_lo - lam2A3_hi
    sinT2_3 = (mp.mpf(1) if delta3_lo <= C_hi
               else min(mp.mpf(1), C_hi / (delta3_lo - C_hi)))
    sinT3 = min(mp.mpf(1), sinT1 + sinT2_3)
    cosT3 = mp.sqrt(1 - sinT3**2) if sinT3 < 1 else mp.mpf(0)
    v0 = abs(v1[0])
    c1sq_lo3 = max(mp.mpf(0), cosT3 * mp.mpf(v0) - sinT3) ** 2

    if verbose:
        mu_str = f"mu*={mu_star:.3f} Nc={Nc}" if mu_star else "mu*=0"
        print(f"\ngamma={gamma_val} N={N} {mu_str}  dim={dim}  [{time.time()-t0:.1f}s]")
        print(f"  mu1_lo (cert) = {float(mu1_lo):.8f}")
        print(f"  lam2_hi = {float(lam2_hi):.8f}   r1 = {float(r1):.3e}   sinT1 = {float(sinT1):.3e}")
        print(f"  ||C|| bounds: Schur={float(C_schur):.6f}  Frob={float(C_frob):.6f}"
              f"  min={float(C_hi):.6f}  (M={M_used})")
        print(f"  lam_max(G^BB)<={float(lamB_hi):.6f}  delta={float(delta3_lo):.6f}"
              f"  sinT={float(sinT3):.6f}")
        print(f"  c1sq_lo (CERTIFICATE, v4) = {float(c1sq_lo3):.6f}  |v1[0]|^2={v0**2:.6f}")
        print()

    return dict(gamma=gamma_val, N=N, mu_star=mu_star, Nc=Nc,
                mu1_lo=mu1_lo, lam2_hi=lam2_hi, r1=r1,
                C_schur=C_schur, C_frob=C_frob, C_hi=C_hi, M_used=M_used,
                lamB_hi=lamB_hi, delta3_lo=delta3_lo, sinT3=sinT3, c1sq_lo3=c1sq_lo3,
                v0=v0, v0sq_float=v0**2, lam1_approx=lam1_approx)


if __name__ == '__main__':
    gamma_val = 1.0
    N = 400
    Nc = 50

    print("=" * 78)
    print(f"v4 test: gamma={gamma_val}  N={N}")
    print("=" * 78)
    sys.stdout.flush()

    rows = []
    for mu_star in [0.0, 0.10, 0.13, 0.15]:
        nc = Nc if mu_star > 0 else None
        rows.append(run_one_v4(gamma_val, N, mu_star=mu_star, Nc=nc))

    print("=" * 78)
    print(f"{'N':>5} {'mu*':>6} {'Nc':>5} {'g_lo':>14} {'||C||':>9}"
          f" {'sinT':>9} {'c1sq_lo':>9} {'c1sq_f64':>9}")
    for r in rows:
        nc = r['Nc'] if r['Nc'] else '-'
        print(f"{r['N']:>5} {r['mu_star']:>6.2f} {str(nc):>5}"
              f" {float(r['mu1_lo']):>14.8f} {float(r['C_hi']):>9.5f}"
              f" {float(r['sinT3']):>9.5f} {float(r['c1sq_lo3']):>9.5f}"
              f" {r['v0sq_float']:>9.5f}")

    g0 = rows[0]['mu1_lo']
    print()
    for r in rows[1:]:
        ratio = (r['mu1_lo'] / g0) ** (1.0 / gamma_val)
        print(f"  mu*={r['mu_star']:.2f}: Tc/Tc0 ~ {float(ratio):.4f} "
              f"({100*(float(ratio)-1):+.1f}%)   [numerical-evidence]")
    print()
