"""
Directional check: scheme B rescue table above always had U1 (paired with
strong diagonal chi0_11=1.5) > U2 (paired with weak diagonal chi0_22=0.8).
Does rescue also occur in the OPPOSITE asymmetry (U2 > U1, i.e. large U on
the weak-diagonal channel)? If rescue only fires in one direction, that's
a real asymmetry to report, not a bug.
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

n_genuine_fail_reverse = 0
n_rescued_reverse = 0
for U2 in grid:
    for U1 in grid:
        if U2 <= U1:
            continue  # only U2 > U1 (weak-diagonal channel has the LARGE U)
        rho = spectral_radius(U1, U2)
        if rho >= 1.0:
            continue
        M = rpa(U1, U2)
        unweighted_ratio = D(M) / d0
        if unweighted_ratio >= 1.0:
            continue
        n_genuine_fail_reverse += 1
        wB = np.array([1/U1, 1/U2])
        dw0 = Dw(chi0, wB)
        if dw0 <= 0:
            continue
        dwr = Dw(M, wB)
        wratio = dwr/dw0
        if wratio >= 1.0:
            n_rescued_reverse += 1

print("reverse-asymmetry (U2>U1) genuine failures:", n_genuine_fail_reverse)
print("reverse-asymmetry rescued by scheme B:", n_rescued_reverse)
