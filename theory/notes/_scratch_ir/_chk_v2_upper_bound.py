"""Does cert_v2's claimed UPPER bound  lam1(M) <= rho + R  actually hold?

cert_v2 header (b) asserts: "For symmetric M, Bauer-Fike / Krylov-Bogoliubov residual
theorem guarantees exists eigenvalue mu in [rho-R, rho+R]. Combined with (a) rho <= lam1,
we get lam1 in [rho, rho+R]."

The premise is true; the conclusion does not follow. The residual theorem places SOME
eigenvalue near rho, not the LARGEST one. Test on the actual H3S matrix.
[numerical-evidence] float64.
"""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_indep_baseline import read_a2f, P_measure, F_T, O_of_f, A2F, MEV_K

w, a2f, cum = read_a2f(A2F)
lam, p, om2 = P_measure(w, a2f)

T = 17.146282                      # meV, the Tc(mu*=0.13) found by the baseline
N = 200                            # smaller N: same phenomenon, faster
F = F_T(2 * N, w, p, T)
n, u, L1, A = O_of_f(F, N)
vp = math.sqrt(om2) / (2.0 * math.pi * T)
Nc = int(math.floor(10.0 * vp))
uc = np.where(n < Nc, u, 0.0)
M = A - 2.0 * (0.13 / lam) * np.outer(uc, uc)

ev, P_ = np.linalg.eigh(M)
lam1, lam2 = ev[-1], ev[-2]
print(f"T={T} meV  N={N}  Nc={Nc}")
print(f"lam1(M) = {lam1:.12f}   lam2(M) = {lam2:.12f}   gap = {lam1-lam2:.3e}")
print()


def rho_R(v):
    v = v / np.linalg.norm(v)
    Mv = M @ v
    rho = float(v @ Mv)
    R = float(np.linalg.norm(Mv - rho * v))
    return rho, R


print(f"{'trial vector v':>34} {'rho':>15} {'R':>12} {'rho+R':>15} "
      f"{'rho+R >= lam1?':>15}")
cases = [
    ("exact top eigenvector", P_[:, -1]),
    ("2nd eigenvector", P_[:, -2]),
    ("3rd eigenvector", P_[:, -3]),
    ("top + 0.3*2nd (mixed)", P_[:, -1] + 0.3 * P_[:, -2]),
    ("u_n profile (naive guess)", u.copy()),
    ("e_0 (crudest guess)", np.eye(N)[0]),
]
for label, v in cases:
    rho, R = rho_R(v)
    ok = "PASS" if rho + R >= lam1 else "*** FAILS ***"
    print(f"{label:>34} {rho:>15.9f} {R:>12.3e} {rho+R:>15.9f} {ok:>15}")

print()
print("Meaning of a FAIL row: rho+R is reported as an upper bound on lam1(M) but is")
print("strictly below the true lam1(M). Any Tc UPPER bound derived from it is unsound.")
