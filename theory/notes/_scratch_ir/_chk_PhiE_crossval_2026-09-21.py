"""CROSS-VALIDATION of Phi_E on >=8 grid points with a SECOND INDEPENDENT SCHEME.

Scheme A (production): folded n>=0 SYMMETRIC pencil
    M = H_mat(N, wbar) - mutilde * 2|u_c><u_c|,  lam_1 by np.linalg.eigvalsh,
    Z eliminated analytically inside H_mat (the AD75/KAY-II cancelled form).
    -> _phiE_engine_2026-09-21.py

Scheme B (independent): UNFOLDED n = -N..N-1 Matsubara axis, Z computed by a
    BRUTE-FORCE truncated sum (no analytic cancellation), NON-symmetric operator
    B[n,m] = (lam(n-m) - mu* th_n th_m)/(Z_n |2m+1|), top eigenvalue by POWER
    ITERATION, Toeplitz action by FFT.  Shares no matrix builder, no eigensolver,
    no Z algebra and no index folding with Scheme A.
    -> _phiE_indep_scheme_2026-09-21.py

Also reported: the CONVENTION sensitivity of the Coulomb cutoff.  The literal
Marsiglio eq.(23) / Marsiglio-Carbotte eq.(56) form puts theta on the SOURCE index
m only (a non-self-adjoint operator); the repo symmetrizes to -2 mu* |u_c><u_c|.
Scheme B runs both, so the symmetrization cost is measured here rather than quoted.

LEVEL: numerical-evidence (float64).  Agreement DIGITS are the deliverable.
"""
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib.util as _ilu  # noqa: E402


def _load(stem):
    spec = _ilu.spec_from_file_location(
        stem.replace("-", "_"), os.path.join(HERE, stem + ".py"))
    mod = _ilu.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


A = _load("_phiE_engine_2026-09-21")
B = _load("_phiE_indep_scheme_2026-09-21")

GRID = [(1.0, 0.10, 3.0), (1.0, 0.13, 3.0), (1.0, 0.10, 10.0),
        (1.5, 0.13, 3.0), (2.0, 0.10, 3.0), (2.0, 0.13, 3.0),
        (2.0, 0.10, 10.0), (2.5, 0.13, 5.0), (3.0, 0.10, 3.0),
        (3.0, 0.13, 10.0), (5.0, 0.10, 3.0), (5.0, 0.13, 10.0)]


def digits(a, b):
    if a == b:
        return 99.0
    return -np.log10(abs(a - b) / max(abs(a), abs(b)))


def main():
    t0 = time.time()
    print("=" * 84)
    print("Phi_E CROSS-VALIDATION: scheme A (folded symmetric pencil) vs")
    print("                        scheme B (unfolded brute-force-Z power iteration)")
    print("=" * 84)
    print()
    print("Step 0 -- mu*=0 anchors (both schemes) against the repo's CERTIFIED table")
    print("          P4_global_tc_certified.md sec.4:")
    cert0 = {1.0: (0.1146122, 0.1147841), 2.0: (0.2115556, 0.2118729),
             5.0: (0.3795354, 0.3801047)}
    tc0A, tc0B = {}, {}
    for lam in (1.0, 1.5, 2.0, 2.5, 3.0, 5.0):
        ta, wa, _ = A.tc_over_wE(lam, 0.0, 1.0, N=260)
        tb, wb, _ = B.tc_indep(lam, 0.0, 1.0, N=3000)
        tc0A[lam], tc0B[lam] = ta, tb
        tag = ""
        if lam in cert0:
            lo, hi = cert0[lam]
            tag = f"  cert[{lo:.7f},{hi:.7f}] A_in={lo<=ta<=hi} B_in={lo<=tb<=hi}"
        print(f"  lam={lam:4.1f}  A={ta:.9f}  B={tb:.9f}  digits={digits(ta,tb):5.2f}{tag}")
    print()

    print("Step 1 -- Phi_E on the cross-validation grid (>=8 points required)")
    print(f"{'lam':>5} {'mu*':>5} {'r':>5} | {'Phi_A':>11} {'NcA':>4} | "
          f"{'Phi_B':>11} {'NcB':>4} | {'digits':>7}")
    print("-" * 74)
    dmin = 99.0
    for lam, mus, r in GRID:
        ta, wa, nca = A.tc_over_wE(lam, mus, r, N=260)
        tb, wbv, ncb = B.tc_indep(lam, mus, r, N=3000)
        if ta is None or tb is None:
            print(f"{lam:5.1f} {mus:5.2f} {r:5.1f} |  root not found in one scheme")
            continue
        pa, pb = ta / tc0A[lam], tb / tc0B[lam]
        d = digits(pa, pb)
        dmin = min(dmin, d)
        print(f"{lam:5.1f} {mus:5.2f} {r:5.1f} | {pa:11.8f} {nca:4d} | "
              f"{pb:11.8f} {ncb:4d} | {d:7.2f}")
    print()
    print(f"  WORST agreement over {len(GRID)} grid points: {dmin:.2f} decimal digits")
    print()

    print("Step 2 -- Coulomb-cutoff SYMMETRIZATION convention cost (scheme B, both)")
    print("  theta on both indices (repo) vs theta on the source index only")
    print(f"{'lam':>5} {'mu*':>5} {'r':>5} | {'Tc sym':>11} {'Tc asym':>11} "
          f"{'rel diff':>10}")
    print("-" * 60)
    for lam, mus, r in [(1.0, 0.10, 5.0), (1.0, 0.10, 10.0), (1.0, 0.10, 20.0),
                        (2.0, 0.13, 5.0), (2.0, 0.13, 10.0), (2.0, 0.13, 20.0)]:
        ts, _, _ = B.tc_indep(lam, mus, r, N=3000, sym=True)
        tas, _, _ = B.tc_indep(lam, mus, r, N=3000, sym=False)
        if ts is None or tas is None:
            continue
        print(f"{lam:5.1f} {mus:5.2f} {r:5.1f} | {ts:11.8f} {tas:11.8f} "
              f"{(tas/ts-1)*100:9.3f}%")
    print()
    print(f"wall {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
