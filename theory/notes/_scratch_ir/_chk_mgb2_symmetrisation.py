import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
# convention (B) as the agent reports it from Mazin Comment / Choi Reply:
#   mu_ij = N_j * K_ij  with K symmetric  =>  mu itself is NOT symmetric.
# Gap equation acts as  Delta_i <- sum_j mu_ij Delta_j , i.e. operator  mu = K diag(N).
# Similarity with psi_i = sqrt(N_i) Delta_i :
#   diag(N)^{1/2} (K diag(N)) diag(N)^{-1/2} = diag(N)^{1/2} K diag(N)^{1/2}
# => the SYMMETRIC object is sqrt(N_i) K_ij sqrt(N_j) = sqrt(N_i/N_j) mu_ij   (NOT sqrt(N_i) mu_ij sqrt(N_j))
sets = {"Mazin contact (3.1:2.6:1.4:1.0)": (3.1, 2.6, 1.4, 1.0),
        "Moon-Kim-Chang   (1.75:2.04:1.61:1.00)": (1.75, 2.04, 1.61, 1.00)}
for tag,(ss, pp, sp, ps) in sets.items():
    r = sp/ps                       # N_pi/N_sigma forced by symmetry of K
    Ns, Np = 1.0, r
    mu = np.array([[ss, sp],[ps, pp]])
    N  = np.array([Ns, Np])
    K  = mu / N[None, :]            # K_ij = mu_ij / N_j
    S  = np.sqrt(N)[:,None]*K*np.sqrt(N)[None,:]        # correct symmetrisation
    W  = np.sqrt(N)[:,None]*mu*np.sqrt(N)[None,:]       # MY formula (wrong for this convention)
    print(f"\n{tag}")
    print(f"  implied N_pi/N_sigma = {r:.4f}")
    print(f"  K symmetric? off-diags {K[0,1]:.6f} vs {K[1,0]:.6f}  -> {np.allclose(K,K.T)}")
    print(f"  K   = {K.ravel().round(5)}  det={np.linalg.det(K):+.5f}  eig={np.linalg.eigvalsh((K+K.T)/2).round(5)}")
    print(f"  S=sqrt(N)K sqrt(N) sym? {np.allclose(S,S.T)}  det={np.linalg.det(S):+.5f}  eig={np.linalg.eigvalsh((S+S.T)/2).round(5)}")
    print(f"  W=sqrt(N)mu sqrt(N) sym? {np.allclose(W,W.T)}  off-diags {W[0,1]:.4f} vs {W[1,0]:.4f}   <-- my formula fails here")
