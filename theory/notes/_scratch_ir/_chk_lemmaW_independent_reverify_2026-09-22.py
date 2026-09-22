"""
Independent reverification of team-lead's Lemma (W) numerical claims:
- Stoner-stable / genuine-failure / negative-row / all-positive-row counts
- exact ceiling (sup over positive weights of D_w[rpa]/D_w[chi0]) via the
  kink-point method they describe, cross-checked against brute-force scan
  of w1/w2, on the SAME chi0 = [[1.5,0.2],[0.2,0.8]].
"""
import numpy as np

chi0 = np.array([[1.5, 0.2],[0.2, 0.8]])

def rows(M):
    return [M[i,i]-sum(abs(M[i,j]) for j in range(2) if j!=i) for i in range(2)]

def rpa(U1,U2):
    U = np.diag([U1,U2])
    return chi0 @ np.linalg.inv(np.eye(2)-U@chi0)

def rho(U1,U2):
    U = np.diag([U1,U2])
    return max(abs(np.linalg.eigvals(U@chi0)))

grid = np.linspace(0.02, 0.8, 60)
D0 = rows(chi0)  # fixed
n_stable=0; n_fail=0; n_negrow=0; n_allpos=0
n_ceiling_ge1_among_allpos=0
n_B_rescue=0

def exact_ceiling(R, D):
    # ratio(t) = min(t*R1,R2)/min(t*D1,D2), t=w1/w2>0, D1,D2>0
    cands = [1e-9, 1e9]
    cands.append(D[1]/D[0])
    if R[0] > 0:
        cands.append(R[1]/R[0])
    best = -1e18
    for t in cands:
        if t<=0: continue
        num = min(t*R[0], R[1])
        den = min(t*D[0], D[1])
        if den<=0: continue
        best = max(best, num/den)
    return best

for U1 in grid:
    for U2 in grid:
        if rho(U1,U2) >= 1.0: continue
        n_stable += 1
        M = rpa(U1,U2)
        R = rows(M)
        unweighted_ratio = min(R)/min(D0)
        if unweighted_ratio >= 1.0: continue
        n_fail += 1
        if min(R) < 0:
            n_negrow += 1
            continue
        n_allpos += 1
        ceil = exact_ceiling(R, D0)
        # brute force cross check
        ts = np.logspace(-6,6,2001)
        bf = max(min(t*R[0],R[1])/min(t*D0[0],D0[1]) for t in ts)
        if abs(ceil-bf) > 1e-6*max(1,abs(bf)):
            print("MISMATCH", U1,U2,ceil,bf)
        if ceil >= 1.0:
            n_ceiling_ge1_among_allpos += 1
        wB = np.array([1/U1,1/U2])
        dw0 = min(wB[0]*D0[0], wB[1]*D0[1])
        dwr = min(wB[0]*R[0], wB[1]*R[1])
        if dwr/dw0 >= 1.0:
            n_B_rescue += 1

print("stable:", n_stable, "fail:", n_fail, "negrow:", n_negrow, "allpos:", n_allpos)
print("allpos with exact ceiling>=1:", n_ceiling_ge1_among_allpos, "/", n_allpos)
print("scheme B rescues:", n_B_rescue, "/", n_allpos, "=", n_B_rescue/n_allpos)
