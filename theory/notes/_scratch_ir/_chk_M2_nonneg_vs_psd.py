import numpy as np
rng = np.random.default_rng(7)

# (G'-2) takes w = Perron vector of M2 and needs w > 0, then
# lam_max(N) <= max_i sum_j N_ij w_j / w_i  = g(2) lam_max(M2)/kappa.
# Question: which property of M2 is load-bearing -- elementwise >= 0, or PSD?

print("=== Q1: for symmetric ELEMENTWISE-NONNEG M2, is lam_max = rho? ===")
bad = 0
for _ in range(50000):
    B = rng.integers(2,6)
    A = rng.uniform(0,3,(B,B)); M = (A+A.T)/2          # elementwise >=0, symmetric
    e = np.linalg.eigvalsh(M)
    if abs(e.max() - abs(e).max()) > 1e-10: bad += 1
print("  elementwise-nonneg: lam_max != rho in", bad, "/50000   (Perron-Frobenius)")
# and are they PSD?
npsd = sum(np.linalg.eigvalsh((lambda A:(A+A.T)/2)(rng.uniform(0,3,(3,3)))).min() < 0
           for _ in range(2000))
print("  ... of which NOT PSD:", npsd, "/2000  -> PSD is NOT implied, and NOT needed")

print("\n=== Q2: with a NEGATIVE off-diagonal entry, does lam_max = rho survive? ===")
bad = 0; ex = None
for _ in range(50000):
    B = rng.integers(2,6)
    A = rng.uniform(-3,3,(B,B)); M = (A+A.T)/2
    e = np.linalg.eigvalsh(M)
    if abs(e.max() - abs(e).max()) > 1e-10:
        bad += 1
        if ex is None and B == 2: ex = (M.copy(), e.copy())
print("  sign-indefinite entries: lam_max != rho in", bad, "/50000")
if ex is not None:
    M,e = ex; print("  example M =", M.round(4), " eig =", e.round(4),
                    " lam_max =", round(e.max(),4), " rho =", round(abs(e).max(),4))

print("\n=== Q3: does the Perron vector still exist (w>0) with a negative entry? ===")
M = np.array([[1.0, -2.0], [-2.0, 1.0]])
e, V = np.linalg.eigh(M)
print("  M =", M.tolist(), "eig =", e, "\n  top eigvec =", V[:,-1].round(4),
      "  (sign-changing -> w_j/w_i step of (G'-2) is undefined)")
print("  lam_max =", e.max(), " rho =", abs(e).max(), " -> ceiling uses sqrt(rho), but chain proves lam_max")

print("\n=== Q4: reducible (decoupled bands) -- w has a zero, division breaks ===")
M = np.array([[2.0,0.0],[0.0,1.0]])
e,V = np.linalg.eigh(M); print("  M block-diagonal, top eigvec =", V[:,-1], "-> w_2 = 0")
print("  patch: apply (G'-2) to M+eps*J (irreducible, w>0), then eps->0 by eigenvalue continuity")
for eps in (1e-2,1e-4,1e-6):
    Me = M + eps*np.ones((2,2)); print(f"    eps={eps:g}  lam_max={np.linalg.eigvalsh(Me).max():.10f}")
