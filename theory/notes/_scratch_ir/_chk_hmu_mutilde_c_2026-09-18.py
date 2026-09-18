"""
Round 5: is the weak-coupling extremality boundary a single constant in
mutilde = mu*(omega_c)/lambda, independent of r = omega_c/omega_rms ?

Repo fit (mustar_extremality_failure_map 2026-09-17 section 4.1, r=5 only):
    s_c = lambda_c/mu*_rms = 1.7733 + 3.7598 mu*_rms
Rewritten in mutilde it is nearly CONSTANT (my arithmetic, see note).
Here I test the r-independence directly, with my own operator.

Also: the Tc root is not unique (N_c staircase + the H_mu hump).  Check that
the violation survives at the SMALLEST varpi root (= largest T = the repo's
Tc := sup{T: lambda*lambda_1 >= 1}).

float64 => numerical-evidence.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, importlib.util, pathlib
from scipy.optimize import brentq

_here = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "chk1", _here / "_chk_hmu_concavity_2026-09-18.py")
chk1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(chk1)
OpFFT, f_mix, top_eigs = chk1.OpFFT, chk1.f_mix, chk1.top_eigs


def lam1(weights, w2s, mt, Nc, Nfac=20):
    vmax = np.sqrt(max(w2s))
    N = max(400, 2 * Nc + 10, int(Nfac * vmax))
    kk = np.arange(2 * N, dtype=float)
    op = OpFFT(f_mix(kk, weights, w2s), N, mt=mt, Nc=Nc)
    return top_eigs(op, k=1)[0][0]


def all_tc_roots(lam, mt, r, lo=1.0, hi=3e4, npts=60):
    """all crossings of lambda*h_mu(varpi; N_c=floor(r*varpi)) = 1, ascending."""
    g = np.exp(np.linspace(np.log(lo), np.log(hi), npts))
    F = np.array([lam * lam1([1.0], [v ** 2], mt, max(1, int(np.floor(r * v)))) - 1.0
                  for v in g])
    roots = []
    for i in range(npts - 1):
        if F[i] < 0 <= F[i + 1]:
            def f(lw):
                v = np.exp(lw)
                return lam * lam1([1.0], [v ** 2], mt,
                                  max(1, int(np.floor(r * v)))) - 1.0
            try:
                roots.append(np.exp(brentq(f, np.log(g[i]), np.log(g[i + 1]),
                                           xtol=1e-7)))
            except Exception:
                pass
    return roots


def best_two_atom(vE, Nc, mt, lam):
    """max over a two-atom grid of  lambda*lambda_1(M(P)) - 1 at fixed (T,N_c)."""
    best = (-1e9, None, None)
    for w in (0.97, 0.99, 0.995, 0.998):
        for sp in (5.0, 8.0, 10.0, 15.0, 25.0, 50.0):
            x1 = vE ** 2 / (w + (1 - w) * sp ** 2)
            rel = lam * lam1([w, 1 - w], [x1, x1 * sp ** 2], mt, Nc) - 1.0
            if rel > best[0]:
                best = (rel, w, sp)
    return best


def violates(lam, mt, r, root_choice="first"):
    roots = all_tc_roots(lam, mt, r)
    if not roots:
        return None
    vE = roots[0] if root_choice == "first" else roots[-1]
    Nc = max(1, int(np.floor(r * vE)))
    rel, w, sp = best_two_atom(vE, Nc, mt, lam)
    return rel, vE, Nc, w, sp, len(roots)


def probe_r_independence():
    print("=" * 78)
    print("is mutilde_c independent of r ?   (lambda = 0.50 fixed)")
    print("=" * 78)
    lam = 0.50
    print(f"  {'r':>5} {'mutilde':>9} {'#roots':>7} {'varpi_E':>10} {'N_c':>7} "
          f"{'lam*Q(P)-1':>12} {'w':>7} {'sp':>6}  verdict")
    for r in (3.0, 5.0, 10.0):
        for mt in (0.50, 0.55, 0.58, 0.62, 0.70):
            out = violates(lam, mt, r)
            if out is None:
                print(f"  {r:5.1f} {mt:9.3f}   no Tc root"); continue
            rel, vE, Nc, w, sp, nr = out
            print(f"  {r:5.1f} {mt:9.3f} {nr:7d} {vE:10.3f} {Nc:7d} {rel:+12.3e} "
                  f"{w:7.3f} {sp:6.1f}  {'VIOLATION' if rel>0 else 'ok'}")
        print()


def probe_root_choice():
    print("=" * 78)
    print("does the violation survive at the SMALLEST varpi root (= true Tc)?")
    print("  lambda=0.50, mutilde=0.60, r=5  (the repo's counterexample point)")
    print("=" * 78)
    roots = all_tc_roots(0.50, 0.60, 5.0)
    print(f"  Tc roots (varpi_E): {[f'{x:.3f}' for x in roots]}")
    for k, vE in enumerate(roots):
        Nc = max(1, int(np.floor(5.0 * vE)))
        rel, w, sp = best_two_atom(vE, Nc, 0.50 if False else 0.60, 0.50)
        print(f"  root #{k} varpi_E={vE:10.3f} N_c={Nc:6d} : lam*Q(P)-1 = "
              f"{rel:+.3e}  (w={w}, sp={sp})  "
              f"{'VIOLATION' if rel>0 else 'ok'}")
    print()


if __name__ == "__main__":
    probe_root_choice()
    probe_r_independence()
