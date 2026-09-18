"""REVIEWER (audit item iii): does mu* leave kappa(wbar), c3(wbar) unchanged (doc Sec 2.4)?
kappa_X = w^2 <v_X, H'' v_X>/lam1(X) with v_X the top eigvec of X in {H, M=H-mt K_c}.
Cross-check against a direct finite-difference sigma^2 coefficient of
k(P) = lam1( sum_j p_j H(w(1+xi_j)) - mt K_c ), 2pk-sym shape, Richardson in sigma."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import rev_common as rc
N = 400
print(f"{'w':>5} {'mt':>5} {'Nc':>4} {'lam1(M)':>9} {'kappa_M':>9} {'kappa_H':>9} {'ratio':>7} | {'c3_M':>9} {'c3_H':>9} {'ratio':>7} | {'kappa_FD':>9} {'FD/kap_M':>8} {'FD/kap_H':>8}")
for w in (0.7, 1.0, 1.55, 2.0, 5.0):
    H = rc.H_mat(N, w); H2 = rc.H_mat(N, w, 2); H3 = rc.H_mat(N, w, 3)
    eH, VH = np.linalg.eigh(H); vH = VH[:, -1]
    kapH = w**2 * (vH @ H2 @ vH) / eH[-1]; c3H = w**3 * (vH @ H3 @ vH) / (6 * eH[-1])
    for mt, Ncs in ((0.0, (0,)), (0.1, (3, 10, 40)), (0.13, (10,)), (0.26, (10,))):
        for Nc in Ncs:
            M = H - mt * rc.Kc(N, Nc)
            e, V = np.linalg.eigh(M); v = V[:, -1]; lam = e[-1]
            kap = w**2 * (v @ H2 @ v) / lam; c3 = w**3 * (v @ H3 @ v) / (6 * lam)
            def kP(s):
                A = 0.5 * rc.H_mat(N, w * (1 - s)) + 0.5 * rc.H_mat(N, w * (1 + s)) - mt * rc.Kc(N, Nc)
                return np.linalg.eigvalsh(A)[-1]
            s1, s2 = 0.02, 0.01
            r1 = (kP(s1) / lam - 1) / s1**2; r2 = (kP(s2) / lam - 1) / s2**2
            kapFD = 2 * (4 * r2 - r1) / 3
            print(f"{w:5.2f} {mt:5.2f} {Nc:4d} {lam:9.5f} {kap:9.5f} {kapH:9.5f} {kap/kapH:7.4f} | {c3:9.5f} {c3H:9.5f} {c3/c3H:7.4f} | {kapFD:9.5f} {kapFD/kap:8.4f} {kapFD/kapH:8.4f}")
