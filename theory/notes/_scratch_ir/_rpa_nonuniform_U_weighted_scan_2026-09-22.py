"""
Scan: does a non-uniform U=diag(U1,U2) let some weighted diagonal-dominance
metric D_w recover D_w[chi^RPA] >= D_w[chi0], where Lemma 2.2 (uniform U)
showed the unweighted/inverse-diagonal weight cannot?

Fixed chi0 = Week2 counterexample [[1.5,0.2],[0.2,0.8]] (kappa=1.875).
Grid over U1,U2 in [0.1,0.8]^2 (covers U1/U2 in [0.5,2.0] as a subset, and
beyond). Filter to Stoner-stable points: spectral radius rho(U@chi0) < 1.

Three weight schemes, each built ONLY from chi0 and U (never from chi^RPA):
  A: w_alpha ~ 1/chi0_{alpha,alpha}   (fixed, independent of U)
  B: w_alpha ~ 1/U_alpha
  C: w_alpha ~ U_alpha

D_w[M] := min_alpha w_alpha * (M_aa - sum_{b!=a} |M_ab|)
Judgment: ratio = D_w[chi^RPA] / D_w[chi0], using the SAME w for both sides.

No conclusion strings baked in below -- print raw values, decide after.
Negative control: uniform U=0.4*I with scheme A must reproduce the known
Lemma 2.2 result (ratio = 0.853, the unweighted 2x2 result), verifying the
harness before trusting the non-uniform grid.
"""
import numpy as np

chi0 = np.array([[1.5, 0.2],
                  [0.2, 0.8]])

def D(M):
    n = M.shape[0]
    return min(M[i, i] - sum(abs(M[i, j]) for j in range(n) if j != i) for i in range(n))

def Dw(M, w):
    n = M.shape[0]
    return min(w[i] * (M[i, i] - sum(abs(M[i, j]) for j in range(n) if j != i)) for i in range(n))

def rpa(U1, U2):
    U = np.diag([U1, U2])
    return chi0 @ np.linalg.inv(np.eye(2) - U @ chi0)

def spectral_radius(U1, U2):
    U = np.diag([U1, U2])
    return max(abs(np.linalg.eigvals(U @ chi0)))

# --- negative control: uniform U=0.4*I, scheme A (inverse diagonal) ---
U1c, U2c = 0.4, 0.4
Mc = rpa(U1c, U2c)
wA_c = np.array([1/chi0[0,0], 1/chi0[1,1]])
ratio_control = Dw(Mc, wA_c) / Dw(chi0, wA_c)
ratio_unweighted_control = D(Mc) / D(chi0)
print("=== negative control (uniform U=0.4I, scheme A) ===")
print("D[chi0]=%.6f D[chi^RPA]=%.6f unweighted ratio=%.6f" % (D(chi0), D(Mc), ratio_unweighted_control))
print("Dw[chi0]=%.6f Dw[chi^RPA]=%.6f weighted(A) ratio=%.6f" % (Dw(chi0, wA_c), Dw(Mc, wA_c), ratio_control))
print()

# --- main grid scan ---
grid = np.linspace(0.1, 0.8, 15)
results = []  # (U1, U2, scheme, ratio, rho)
n_total = 0
n_stable = 0
n_neg_ratio_ge1 = {"A": 0, "B": 0, "C": 0}

for U1 in grid:
    for U2 in grid:
        n_total += 1
        rho = spectral_radius(U1, U2)
        if rho >= 1.0:
            continue
        n_stable += 1
        M = rpa(U1, U2)
        d0 = D(chi0)
        if d0 <= 0:
            continue
        wA = np.array([1/chi0[0,0], 1/chi0[1,1]])
        wB = np.array([1/U1, 1/U2])
        wC = np.array([U1, U2])
        for name, w in [("A", wA), ("B", wB), ("C", wC)]:
            dw0 = Dw(chi0, w)
            dwr = Dw(M, w)
            if dw0 <= 0:
                continue
            ratio = dwr / dw0
            results.append((U1, U2, name, ratio, rho))
            if ratio >= 1.0:
                n_neg_ratio_ge1[name] += 1

print("=== budget scan ===")
print("total grid points:", n_total, "  Stoner-stable (rho<1):", n_stable)
print("evaluated (scheme,ratio) pairs:", len(results))

for name in ["A", "B", "C"]:
    sub = [r for r in results if r[2] == name]
    ratios = [r[3] for r in sub]
    best = max(sub, key=lambda r: r[3])
    print(f"--- scheme {name} ---")
    print(f"  max ratio = {best[3]:.6f} at U1={best[0]:.3f}, U2={best[1]:.3f}, rho={best[4]:.3f}")
    print(f"  #points with ratio>=1: {n_neg_ratio_ge1[name]} / {len(sub)}")
    print(f"  min ratio = {min(ratios):.6f}")

# also report U1/U2 in [0.5,2.0] subset explicitly, as team-lead specified
print()
print("=== restricted to U1/U2 in [0.5,2.0] (team-lead's literal range) ===")
for name in ["A", "B", "C"]:
    sub = [r for r in results if r[2] == name and 0.5 <= r[0]/r[1] <= 2.0]
    if not sub:
        print(f"scheme {name}: no stable points in this sub-range")
        continue
    best = max(sub, key=lambda r: r[3])
    print(f"scheme {name}: max ratio = {best[3]:.6f} at U1={best[0]:.3f}, U2={best[1]:.3f}")
