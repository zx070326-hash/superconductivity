"""Corrected beta route (trace route RETRACTED: G is not PSD, beta_trace came out negative).

Sound chain, each step with its hypothesis spelled out:
  (1) lam2(M) <= lam2(A)              M = A - c u_c u_c^T, c >= 0  (Weyl, PSD subtracted)
  (2) lam2(A) <= lam2(G)              A = G - diag(L1), L1 >= 0    (Weyl, PSD subtracted)
  (3) lam2(G) <= sqrt(||G||_F^2 - rho_G^2)   PROVIDED rho_G >= 0, where rho_G is the
      Rayleigh quotient of G at any trial vector.  Proof: rho_G >= 0 and rho_G <= lam1(G)
      give lam1(G)^2 >= rho_G^2, and lam2(G)^2 <= sum_{i>=2} lam_i(G)^2
      = ||G||_F^2 - lam1(G)^2 <= ||G||_F^2 - rho_G^2.   NO positive-definiteness needed.
  (4) Temple: beta >= lam2(M) and rho > beta  =>  lam1(M) <= rho + R^2/(rho - beta).

This script checks the PREMISES (rho_G >= 0, beta >= true lam2(M), beta < rho) across the
whole temperature range a Tc bisection would visit.  [numerical-evidence] float64.
"""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_indep_baseline import read_a2f, P_measure, F_T, O_of_f, A2F, MEV_K

w, a2f, cum = read_a2f(A2F)
lam, p, om2 = P_measure(w, a2f)
N = 400

print(f"N={N}, mu*=0.13, r=10.  Temple premises across the bisection range.")
print(f"{'T[meV]':>9} {'T[K]':>8} {'rho':>12} {'R':>10} {'rho_G':>9} {'beta_fro':>9} "
      f"{'lam2(M)':>10} {'beta>=lam2':>11} {'beta<rho':>9} {'Temple ub - lam1':>17}")
for T in (2.0, 5.0, 10.0, 15.0, 17.146282, 20.0, 24.416658, 30.0, 40.0, 60.0):
    F = F_T(2 * N, w, p, T)
    n, u, L1, A = O_of_f(F, N)
    G = A + np.diag(L1)
    vp = math.sqrt(om2) / (2.0 * math.pi * T)
    Nc = int(math.floor(10.0 * vp))
    uc = np.where(n < Nc, u, 0.0)
    M = A - 2.0 * (0.13 / lam) * np.outer(uc, uc)

    eM, PM = np.linalg.eigh(M)
    v = PM[:, -1]
    Mv = M @ v
    rho = float(v @ Mv)
    R = float(np.linalg.norm(Mv - rho * v))
    vG = np.linalg.eigh(G)[1][:, -1]
    rho_G = float(vG @ (G @ vG))
    beta = math.sqrt(max(float(np.sum(G * G)) - rho_G ** 2, 0.0))

    ok_beta = beta >= eM[-2]
    ok_prem = beta < rho
    ub = rho + R ** 2 / (rho - beta) if ok_prem else float("nan")
    print(f"{T:>9.4f} {T*MEV_K:>8.2f} {rho:>12.8f} {R:>10.2e} {rho_G:>9.5f} {beta:>9.5f} "
          f"{eM[-2]:>10.2e} {str(ok_beta):>11} {str(ok_prem):>9} {ub-eM[-1]:>17.2e}",
          flush=True)
print()
print("Tc criterion crosses 1/lam = %.6f ; the rows near it are the ones that matter."
      % (1.0 / lam))
