"""A SOUND replacement for cert_v2's broken upper bound, tested numerically first.

Temple's inequality: if beta >= lam2(M) and rho > beta, then
      lam1(M) <= rho + R^2/(rho - beta).
Everything here is O(N^2) and interval-able (no eigenvectors needed beyond the trial v,
whose quality affects only WIDTH, not soundness):

  beta := sqrt(||M||_F^2 - rho^2)   is a valid choice, because
      lam2^2 <= sum_{i>=2} lam_i^2 = ||M||_F^2 - lam1^2 <= ||M||_F^2 - rho^2   (rho <= lam1).

Question this script answers: is beta < rho on the actual H3S matrix (Temple's premise),
and how tight is the resulting bound?  [numerical-evidence] float64.
"""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_indep_baseline import read_a2f, P_measure, F_T, O_of_f, A2F

w, a2f, cum = read_a2f(A2F)
lam, p, om2 = P_measure(w, a2f)

print(f"{'N':>5} {'T[meV]':>9} {'lam1':>14} {'rho':>14} {'R':>10} {'||M||_F':>10} "
      f"{'beta':>12} {'beta<rho':>9} {'Temple ub':>14} {'ub-lam1':>11}")
for N in (100, 200, 400):
    for T in (17.146282, 24.0):
        F = F_T(2 * N, w, p, T)
        n, u, L1, A = O_of_f(F, N)
        vp = math.sqrt(om2) / (2.0 * math.pi * T)
        Nc = int(math.floor(10.0 * vp))
        uc = np.where(n < Nc, u, 0.0)
        M = A - 2.0 * (0.13 / lam) * np.outer(uc, uc)

        ev, P_ = np.linalg.eigh(M)
        lam1 = ev[-1]
        v = P_[:, -1]                       # trial vector (float64 quality only)
        Mv = M @ v
        rho = float(v @ Mv)
        R = float(np.linalg.norm(Mv - rho * v))
        fro2 = float(np.sum(M * M))
        beta = math.sqrt(max(fro2 - rho ** 2, 0.0))
        ok = beta < rho
        ub = rho + R ** 2 / (rho - beta) if ok else float("nan")
        print(f"{N:>5} {T:>9.4f} {lam1:>14.10f} {rho:>14.10f} {R:>10.2e} "
              f"{math.sqrt(fro2):>10.5f} {beta:>12.6f} {str(ok):>9} {ub:>14.10f} "
              f"{ub-lam1:>11.2e}", flush=True)

print()
print("If beta < rho holds, Temple gives a SOUND upper bound on lam1(M) using only")
print("rho, R, ||M||_F -- all O(N^2), all outward-roundable in interval arithmetic.")
print("If beta >= rho, Temple is unavailable and a certified lam2 bound is needed")
print("(repo has an analytic route: general rearrangement identity, 2026-09-13).")
