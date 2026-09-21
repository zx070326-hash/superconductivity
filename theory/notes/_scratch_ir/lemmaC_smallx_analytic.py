"""Fully hand-analytic small-x chain (all constants explicit; series evaluated here only to
read off the reach x0).  ratio(x) = x*kappa4^2 / (gt_lo(x) * F_lo(x)) must be < 1."""
import numpy as np
K = np.arange(1, 200001, dtype=float)
kappa4 = 1.95                      # u-weighted Schur bound on ||O+[1/k^4]|| (row n=1), hand-derivable
g2_lo  = 1.3081                    # CERT-A (paper-and-pencil) lower bound on g(2)
lam2_G2 = 13/25                    # Lemma 3 of P4_g2_analytic
for x in [0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.1]:
    f   = x/(K*K+x); phi = K*K/(K*K+x)**3
    c   = 1/(x*(1+x)**2)
    m1  = np.sum(c*f - phi)                       # = max_n [cT_f - T_phi](n)  (n=0)
    hF  = np.sqrt(np.sum(K*(c*f-phi)**2))         # Frobenius norm of Hankel(cf-phi)
    h_lo = x*(g2_lo - kappa4*x)                   # Rayleigh of CERT-A vector on A(x)
    F_lo = c*h_lo - m1 - hF
    gt_lo = g2_lo - lam2_G2 - 2*kappa4*x          # rescaled gap
    ratio = x*kappa4**2/(gt_lo*F_lo) if F_lo>0 and gt_lo>0 else np.inf
    print(f"x={x:6.3f} w={np.sqrt(x):.3f}  c*h_lo={c*h_lo:.4f} m1={m1:.4f} |H|_F={hF:.4f} F_lo={F_lo:.4f} gt_lo={gt_lo:.4f}  ratio={ratio:.3f}")
