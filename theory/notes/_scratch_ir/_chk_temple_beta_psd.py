"""Can a CHEAP certified beta >= lam2(M) with beta < rho be had?  (Temple's premise)

Chain to be tested:
    lam2(M) <= lam2(A)        (M = A - c u_c u_c^T, c>0: negative rank-one interlacing)
    lam2(A) <= lam2(G)        (A = G - diag(L1), diag(L1) >= 0 PSD subtracted)
    G PSD?                    (F(k) = int vp^2/(k^2+vp^2) P is a positive-definite
                               sequence in k; Toeplitz part PSD, Hankel part to be checked)
    then, G PSD  =>  lam2(G) <= min( tr G - rho_G ,  sqrt(||G||_F^2 - rho_G^2) )
                     with rho_G = Rayleigh quotient of G at the trial vector (<= lam1(G)).

Reports the true lam2's alongside, so we can see how much each step gives away.
[numerical-evidence] float64.
"""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_indep_baseline import read_a2f, P_measure, F_T, O_of_f, A2F

w, a2f, cum = read_a2f(A2F)
lam, p, om2 = P_measure(w, a2f)

for N in (200, 400):
    for T in (17.146282, 24.0):
        F = F_T(2 * N, w, p, T)
        n, u, L1, A = O_of_f(F, N)
        G = A + np.diag(L1)                       # strip the -L1 diagonal
        vp = math.sqrt(om2) / (2.0 * math.pi * T)
        Nc = int(math.floor(10.0 * vp))
        uc = np.where(n < Nc, u, 0.0)
        M = A - 2.0 * (0.13 / lam) * np.outer(uc, uc)

        eM = np.linalg.eigvalsh(M); eA = np.linalg.eigvalsh(A); eG = np.linalg.eigvalsh(G)
        vM = np.linalg.eigh(M)[1][:, -1]
        rho = float(vM @ (M @ vM))
        vG = np.linalg.eigh(G)[1][:, -1]
        rho_G = float(vG @ (G @ vG))
        trG = float(np.trace(G))
        froG = float(np.sqrt(np.sum(G * G)))
        b_trace = trG - rho_G
        b_fro = math.sqrt(max(froG ** 2 - rho_G ** 2, 0.0))
        beta = min(b_trace, b_fro)

        print(f"N={N} T={T:.4f}")
        print(f"   lam1(M)={eM[-1]:.8f}  lam2(M)={eM[-2]:.3e}  lam2(A)={eA[-2]:.3e}  "
              f"lam2(G)={eG[-2]:.3e}   min eig G = {eG[0]:.3e} "
              f"({'G PSD' if eG[0] > -1e-12 else '*** G NOT PSD ***'})")
        print(f"   rho={rho:.8f}  tr G={trG:.6f}  ||G||_F={froG:.6f}  rho_G={rho_G:.6f}")
        print(f"   beta_trace={b_trace:.6f}  beta_fro={b_fro:.6f}  ->  beta={beta:.6f}   "
              f"beta<rho? {beta < rho}")
        if beta < rho:
            Mv = M @ vM
            R = float(np.linalg.norm(Mv - rho * vM))
            ub = rho + R ** 2 / (rho - beta)
            print(f"   TEMPLE ub = {ub:.12f}   (lam1 = {eM[-1]:.12f}, slack {ub-eM[-1]:.2e})")
        print(flush=True)
