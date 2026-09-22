"""
Inspect scheme B's 20 rescued points: what do they look like? Any clean
pattern (asymmetry direction, magnitude)? Also check: is rescue happening
via a genuine "different row bottlenecks chi0 vs chi^RPA" effect, or is it
just "row1 was already fine, weight B always favors row1 here"?
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

print(f"{'U1':>6} {'U2':>6} {'U1/U2':>7} {'rho':>6} {'uw_ratio':>9} {'wB_ratio':>9} {'bottleneck0':>4} {'bottleneckM':>4}")
for U1 in grid:
    for U2 in grid:
        rho = spectral_radius(U1, U2)
        if rho >= 1.0:
            continue
        M = rpa(U1, U2)
        unweighted_ratio = D(M) / d0
        if unweighted_ratio >= 1.0:
            continue
        wB = np.array([1/U1, 1/U2])
        r0 = rows(chi0)
        rM = rows(M)
        dw0 = Dw(chi0, wB)
        if dw0 <= 0:
            continue
        dwr = Dw(M, wB)
        wratio = dwr/dw0
        if wratio >= 1.0:
            b0 = int(np.argmin([wB[i]*r0[i] for i in range(2)]))
            bM = int(np.argmin([wB[i]*rM[i] for i in range(2)]))
            print(f"{U1:6.3f} {U2:6.3f} {U1/U2:7.3f} {rho:6.3f} {unweighted_ratio:9.4f} {wratio:9.4f} {b0:>11} {bM:>11}")
