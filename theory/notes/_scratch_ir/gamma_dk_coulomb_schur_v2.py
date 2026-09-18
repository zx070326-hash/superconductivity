"""
γ-model Davis-Kahan pipeline, v2 tail bounds (SOUND weighted Schur).

New file; does NOT modify gamma_dk_coulomb_schur.py / gamma_tail_bound.py,
per repo convention.

DIFFERENCE FROM v1 (gamma_dk_coulomb_schur.py)
----------------------------------------------
Tail bounds come from gamma_tail_bound_v2 (weighted Schur, p=q=u) instead of
gamma_tail_bound (unweighted Schur).  See gamma_tail_bound_v2.py for the three
one-sided defects (D1 integral direction, D2 Hankel substitution, D3 worst-row
assumption) that made the v1 derivation not a chain of inequalities.

The v2 bounds are BOTH sound AND numerically smaller than v1's at every γ, N
tested, so the previously reported c₁²_lo values remain valid lower bounds
(sinΘ is increasing in ||E||, c₁²_lo decreasing) — they were underivable as
stated, not false.  This script recomputes them soundly and tighter.

Everything else (Rayleigh lower bound, inertia λ₂ upper bound, residual,
Davis-Kahan two-step angle, Coulomb projector) is unchanged and re-imported.

Strictness: certificate.
"""
import sys, os, time
import numpy as np
from mpmath import iv, mp

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from gamma_dk_coulomb import (
    build_G_matrix as build_G_dk, to_float_matrix, rayleigh_lower,
    residual_norm, find_lambda2_upper, lo, hi,
)
from gamma_tail_bound_v2 import E_bound_hi, lamB_bound_hi
from gamma_tail_bound import (
    build_G_matrix as build_G_tail, schur_C_bound_hi, HBB_bound_hi,
)


# ---------------------------------------------------------------------------
# Davis-Kahan step 2, two reference operators
# ---------------------------------------------------------------------------
# Common lemma (proved here so the direction is explicit).  Let A' be Hermitian
# with λ_max(A') = μ₁ simple, unit eigenvector ê₁, and spec(A')\{μ₁} ⊂ (-∞, μ₁-δ].
# Let 𝔊 = A' + E', ν₁ = λ_max(𝔊) with unit eigenvector w = cosθ ê₁ + sinθ z,
# z ⊥ ê₁.  Project 𝔊w = ν₁w onto ê₁^⊥ with P = I - ê₁ê₁*:
#       (A'|_P - ν₁)(sinθ · z) = -P E' w
#   =>  sinθ · dist(ν₁, spec(A'|_P)) <= ||E'||.
# Weyl gives ν₁ >= μ₁ - ||E'||, so dist(ν₁, spec(A'|_P)) >= δ - ||E'||, hence
#       sinθ <= ||E'|| / (δ - ||E'||)            [requires ||E'|| < δ]
#
#   (v2) A' = 𝔊^(N) ⊕ 0        E' = [[0,C],[Cᵀ, G^BB]]   ||E'|| <= ||C||+||G^BB||
#                              δ  = μ₁ - max(λ₂^head, 0)
#   (v3) A' = 𝔊^(N) ⊕ G^BB     E' = [[0,C],[Cᵀ, 0    ]]   ||E'|| = ||C||
#                              δ  = μ₁ - max(λ₂^head, λ_max(G^BB))
# v3 is far better here because ||G^BB|| is dominated by the NEGATIVE diagonal
# -𝔊₁, which contributes nothing to λ_max(G^BB).  ê₁ = (v₁,0) is an eigenvector
# of A' in BOTH cases (it has zero tail), and μ₁ = λ_max(A') in both (v3 needs
# λ_max(G^BB) < μ₁, true by many orders of magnitude).


def run_one_v2(gamma_val, N, mu_star=0.0, Nc=None, dps=30, verbose=True,
               compare_v1=True):
    mp.dps = dps
    iv.dps = dps
    t0 = time.time()

    data = build_G_dk(N, gamma_val, mu_star=mu_star, Nc=Nc)
    H, dim = data['H'], data['dim']

    A = to_float_matrix(H, dim)
    evals, evecs = np.linalg.eigh(A)
    lam1_approx, lam2_approx = evals[-1], evals[-2]
    v1 = evecs[:, -1]

    mu1_lo = lo(rayleigh_lower(H, dim, v1.tolist()))
    r1 = hi(residual_norm(H, dim, v1.tolist(), iv.mpf(mu1_lo)))
    c_upper, margin_used = find_lambda2_upper(H, dim, lam2_approx,
                                              lam1_approx - lam2_approx)
    lam2_hi = hi(c_upper)

    gap1_lo = mu1_lo - lam2_hi
    sinT1 = mp.mpf(1) if gap1_lo <= r1 else min(mp.mpf(1), r1 / (gap1_lo - r1))

    # --- tail: v2 weighted-Schur bounds --------------------------------
    # The Coulomb term -μ* P_{Nc} lives entirely in the head block (Nc <= N),
    # so the μ*=0 operator is the right object for C and 𝔊^BB.
    Ehi, C_hi, HBB_hi, parts, A_sup, B_sup = E_bound_hi(N, gamma_val, dps=dps)

    lam2A_hi = max(lam2_hi, mp.mpf(0))
    delta_lo = mu1_lo - lam2A_hi
    sinT2 = mp.mpf(1) if delta_lo <= Ehi else min(mp.mpf(1), Ehi / (delta_lo - Ehi))

    sinT = min(mp.mpf(1), sinT1 + sinT2)
    cosT = mp.sqrt(1 - sinT**2) if sinT < 1 else mp.mpf(0)
    v0 = abs(v1[0])
    c1sq_lo = max(mp.mpf(0), cosT * mp.mpf(v0) - sinT) ** 2

    # --- v3: block-diagonal reference head ⊕ G^BB ----------------------
    lamB_hi, _parts = lamB_bound_hi(N, gamma_val, dps=dps)
    lam2A3_hi = max(lam2_hi, lamB_hi)
    delta3_lo = mu1_lo - lam2A3_hi
    sinT2_3 = (mp.mpf(1) if delta3_lo <= C_hi
               else min(mp.mpf(1), C_hi / (delta3_lo - C_hi)))
    sinT3 = min(mp.mpf(1), sinT1 + sinT2_3)
    cosT3 = mp.sqrt(1 - sinT3**2) if sinT3 < 1 else mp.mpf(0)
    c1sq_lo3 = max(mp.mpf(0), cosT3 * mp.mpf(v0) - sinT3) ** 2

    v1res = None
    if compare_v1:
        td = build_G_tail(N, gamma_val, dps=dps)
        C1 = schur_C_bound_hi(td, N, dps=dps)
        H1 = HBB_bound_hi(td, N, dps=dps)
        E1 = C1 + H1
        s2_1 = mp.mpf(1) if delta_lo <= E1 else min(mp.mpf(1), E1 / (delta_lo - E1))
        s_1 = min(mp.mpf(1), sinT1 + s2_1)
        c_1 = mp.sqrt(1 - s_1**2) if s_1 < 1 else mp.mpf(0)
        v1res = dict(C=C1, HBB=H1, E=E1, sinT=s_1,
                     c1sq=max(mp.mpf(0), c_1 * mp.mpf(v0) - s_1) ** 2)

    if verbose:
        mu_str = f"mu*={mu_star:.3f} Nc={Nc}" if mu_star else "mu*=0"
        print(f"gamma={gamma_val} N={N} {mu_str}  dim={dim}  [{time.time()-t0:.1f}s]")
        print(f"  mu1_lo (cert, INFINITE operator) = {mu1_lo}")
        print(f"  lam2_hi = {lam2_hi}   r1 = {float(r1):.3e}   sinT1 = {float(sinT1):.3e}")
        print(f"  [v2 weighted Schur] ||C||<={float(C_hi):.6f}  ||G^BB||<={float(HBB_hi):.6f}"
              f"  ||E||<={float(Ehi):.6f}")
        print("     parts(G^BB): " + ", ".join(f"{k}={float(v):.6f}" for k, v in parts.items()))
        if v1res:
            print(f"  [v1 unweighted   ] ||C||<={float(v1res['C']):.6f}"
                  f"  ||G^BB||<={float(v1res['HBB']):.6f}  ||E||<={float(v1res['E']):.6f}"
                  f"   (v1/v2 = {float(v1res['E']/Ehi):.2f}x)")
        print(f"  delta_lo={float(delta_lo):.6f}  sinT2={float(sinT2):.6f}"
              f"  sinT={float(sinT):.6f}")
        print(f"  |v1[0]|^2 (float64) = {v0**2:.6f}")
        print(f"  c1sq_lo (CERTIFICATE, v2) = {float(c1sq_lo):.6f}"
              + (f"   [v1 gave {float(v1res['c1sq']):.6f}]" if v1res else ""))
        print(f"  [v3 block-diag DK] lam_max(G^BB)<={float(lamB_hi):.6f}  ||E'||=||C||<={float(C_hi):.6f}"
              f"  delta={float(delta3_lo):.6f}  sinT={float(sinT3):.6f}")
        print(f"  c1sq_lo (CERTIFICATE, v3) = {float(c1sq_lo3):.6f}")
        print()

    return dict(gamma=gamma_val, N=N, mu_star=mu_star, Nc=Nc,
                mu1_lo=mu1_lo, lam2_hi=lam2_hi, r1=r1,
                C_hi=C_hi, HBB_hi=HBB_hi, Ehi=Ehi, parts=parts,
                A_sup=A_sup, B_sup=B_sup,
                sinT1=sinT1, sinT2=sinT2, sinT=sinT,
                c1sq_lo=c1sq_lo, v0=v0, v0sq_float=v0**2,
                lamB_hi=lamB_hi, delta3_lo=delta3_lo, sinT3=sinT3, c1sq_lo3=c1sq_lo3,
                lam1_approx=lam1_approx, v1res=v1res)


if __name__ == '__main__':
    gamma_val = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 400

    print("=" * 78)
    print(f"gamma-model DK pipeline -- v2 SOUND weighted-Schur tail bounds  (gamma={gamma_val})")
    print("=" * 78)
    rows = []
    for mu_star, Nc in [(0.0, None), (0.10, 50), (0.13, 50), (0.15, 50)]:
        rows.append(run_one_v2(gamma_val, N, mu_star=mu_star, Nc=Nc,
                               compare_v1=(mu_star == 0.0)))

    print("=" * 78)
    print(f"{'N':>5} {'mu*':>6} {'Nc':>5} {'g_lo':>14} {'||E||v2':>9} {'sinT_v2':>9}"
          f" {'c1sq_v2':>9} {'lamB':>9} {'sinT_v3':>9} {'c1sq_v3':>9} {'c1sq_f64':>9}")
    for r in rows:
        nc = r['Nc'] if r['Nc'] else '-'
        print(f"{r['N']:>5} {r['mu_star']:>6.2f} {str(nc):>5} {float(r['mu1_lo']):>14.8f}"
              f" {float(r['Ehi']):>9.5f} {float(r['sinT']):>9.5f}"
              f" {float(r['c1sq_lo']):>9.5f} {float(r['lamB_hi']):>9.5f}"
              f" {float(r['sinT3']):>9.5f} {float(r['c1sq_lo3']):>9.5f}"
              f" {r['v0sq_float']:>9.5f}")
    g0 = rows[0]['mu1_lo']
    print()
    for r in rows[1:]:
        ratio = (r['mu1_lo'] / g0) ** (1.0 / gamma_val)
        print(f"  mu*={r['mu_star']:.2f}: Tc/Tc0 ~ {float(ratio):.4f} "
              f"({100*(float(ratio)-1):+.1f}%)   [ratio of two LOWER bounds:"
              f" numerical-evidence]")
