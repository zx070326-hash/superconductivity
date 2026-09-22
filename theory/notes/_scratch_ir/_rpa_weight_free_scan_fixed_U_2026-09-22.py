"""
Sharper test: at the SAME chi0 and the SAME U (U1=U2=0.4, the known Route-1B
counterexample point), can ANY choice of weight (w1,w2), not tied to a
specific formula, push Dw[chi^RPA]/Dw[chi0] >= 1?

Key structural fact (from Lemma 2.2's own argument): Dw's min over rows means
that for any weight, the achieved ratio is EXACTLY the per-row ratio of
whichever row is the bottleneck under that weight for BOTH chi0 and chi^RPA
(when they agree on which row bottlenecks). So the only ratios reachable by
ANY weight are: row1(rpa)/row1(chi0), row2(rpa)/row2(chi0), or something in
between if the bottleneck flips between chi0 and chi^RPA (in which case Dw
ratio can be driven arbitrarily since the two mins are of different rows on
each side -- but then it's an artifact of an unnatural weight, not a genuine
rescue of the SAME structural inequality).

This scans w1/w2 in a wide log-range and reports the max ratio directly,
plus the two "pure row" ratios as reference bounds.
"""
import numpy as np

chi0 = np.array([[1.5, 0.2],
                  [0.2, 0.8]])

def rows(M):
    n = M.shape[0]
    return [M[i, i] - sum(abs(M[i, j]) for j in range(n) if j != i) for i in range(n)]

U1, U2 = 0.4, 0.4
U = np.diag([U1, U2])
M = chi0 @ np.linalg.inv(np.eye(2) - U @ chi0)

r0 = rows(chi0)
rM = rows(M)
print("chi0 rows:", r0, " chi^RPA rows:", rM)
print("pure row-1 ratio:", rM[0]/r0[0])
print("pure row-2 ratio:", rM[1]/r0[1])

best = (-1e18, None)
for logw in np.linspace(-6, 6, 4001):
    w1 = 1.0
    w2 = 10**logw
    dw0 = min(w1*r0[0], w2*r0[1])
    dwr = min(w1*rM[0], w2*rM[1])
    if dw0 <= 0:
        continue
    ratio = dwr/dw0
    if ratio > best[0]:
        best = (ratio, w2)

print("scan w2/w1 in [1e-6,1e6]: best ratio =", best[0], " at w2/w1 =", best[1])
