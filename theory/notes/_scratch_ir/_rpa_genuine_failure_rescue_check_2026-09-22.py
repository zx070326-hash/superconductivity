"""
Final, decisive test: restrict to grid points where the UNWEIGHTED metric
already FAILS (D[chi^RPA]/D[chi0] < 1) -- i.e. exclude points where any
success is already explained without weighting at all (see prior script's
finding that non-uniform U alone often already avoids degradation).
Among those genuine-failure points, does ANY of schemes A/B/C rescue the
weighted ratio to >= 1?
"""
import numpy as np

chi0 = np.array([[1.5, 0.2],
                  [0.2, 0.8]])

def rows(M):
    n = M.shape[0]
    return [M[i, i] - sum(abs(M[i, j]) for j in range(n) if j != i) for i in range(n)]

def D(M):
    return min(rows(M))

def Dw(M, w):
    r = rows(M)
    return min(w[i] * r[i] for i in range(len(r)))

def rpa(U1, U2):
    U = np.diag([U1, U2])
    return chi0 @ np.linalg.inv(np.eye(2) - U @ chi0)

def spectral_radius(U1, U2):
    U = np.diag([U1, U2])
    return max(abs(np.linalg.eigvals(U @ chi0)))

grid = np.linspace(0.05, 0.9, 25)
d0 = D(chi0)

n_genuine_fail = 0
n_rescued = {"A": 0, "B": 0, "C": 0}
worst_after_weight = {"A": 1e18, "B": 1e18, "C": 1e18}

for U1 in grid:
    for U2 in grid:
        rho = spectral_radius(U1, U2)
        if rho >= 1.0:
            continue
        M = rpa(U1, U2)
        unweighted_ratio = D(M) / d0
        if unweighted_ratio >= 1.0:
            continue  # not a genuine-failure point, skip
        n_genuine_fail += 1
        wA = np.array([1/chi0[0,0], 1/chi0[1,1]])
        wB = np.array([1/U1, 1/U2])
        wC = np.array([U1, U2])
        for name, w in [("A", wA), ("B", wB), ("C", wC)]:
            dw0 = Dw(chi0, w)
            if dw0 <= 0:
                continue
            dwr = Dw(M, w)
            wratio = dwr / dw0
            if wratio >= 1.0:
                n_rescued[name] += 1
            worst_after_weight[name] = min(worst_after_weight[name], wratio)

print("genuine-failure grid points (unweighted ratio<1):", n_genuine_fail)
for name in ["A", "B", "C"]:
    print(f"scheme {name}: rescued (weighted ratio>=1) count = {n_rescued[name]} / {n_genuine_fail};  worst weighted ratio among failures = {worst_after_weight[name]:.6f}")
