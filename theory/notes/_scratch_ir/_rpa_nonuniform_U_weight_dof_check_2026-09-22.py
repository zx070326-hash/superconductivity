"""
Follow-up: distinguish "weighting genuinely changes which row is the
bottleneck / the ratio" from the trivial case where weighted ratio equals
the unweighted per-row ratio (which happens whenever the SAME row is the
minimizer for both chi0 and chi^RPA under a given weight -- in that case
Dw ratio = w_a*row_a(rpa) / (w_a*row_a(chi0)) = row_a(rpa)/row_a(chi0),
independent of w entirely, so any weight scheme "succeeds" for a reason
that has nothing to do with the weight).

This directly targets dof-check item 4: is the printed "ratio>=1" set of
points evidence that WEIGHTING helps, or just that the unweighted 2x2
metric itself is non-monotonic in nonuniform U (a fact about chi0/U, zero
extra info from w)?
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

def argminrow(M):
    r = rows(M)
    return int(np.argmin(r))

def rpa(U1, U2):
    U = np.diag([U1, U2])
    return chi0 @ np.linalg.inv(np.eye(2) - U @ chi0)

def spectral_radius(U1, U2):
    U = np.diag([U1, U2])
    return max(abs(np.linalg.eigvals(U @ chi0)))

grid = np.linspace(0.1, 0.8, 15)
d0 = D(chi0)
bottleneck0 = argminrow(chi0)  # which row is chi0's own bottleneck (fixed, weight-independent up to sign of w)

n_stable = 0
n_same_bottleneck_both = {"A": 0, "B": 0, "C": 0}
n_diff_from_unweighted = {"A": 0, "B": 0, "C": 0}
max_genuine_ratio = {"A": (-1e9, None), "B": (-1e9, None), "C": (-1e9, None)}

for U1 in grid:
    for U2 in grid:
        rho = spectral_radius(U1, U2)
        if rho >= 1.0:
            continue
        n_stable += 1
        M = rpa(U1, U2)
        unweighted_ratio = D(M) / d0
        bM = argminrow(M)

        wA = np.array([1/chi0[0,0], 1/chi0[1,1]])
        wB = np.array([1/U1, 1/U2])
        wC = np.array([U1, U2])
        for name, w in [("A", wA), ("B", wB), ("C", wC)]:
            dw0 = Dw(chi0, w)
            dwr = Dw(M, w)
            if dw0 <= 0:
                continue
            wratio = dwr / dw0
            # is the weighted ratio numerically the SAME row-ratio as unweighted?
            # (i.e. weight didn't actually change which row bottlenecks, on either side)
            bottleneck0_w = int(np.argmin([w[i]*rows(chi0)[i] for i in range(2)]))
            bottleneckM_w = int(np.argmin([w[i]*rows(M)[i] for i in range(2)]))
            same_row_both_sides = (bottleneck0_w == bottleneckM_w)
            if same_row_both_sides:
                n_same_bottleneck_both[name] += 1
            if abs(wratio - unweighted_ratio) > 1e-9:
                n_diff_from_unweighted[name] += 1
                if wratio > max_genuine_ratio[name][0]:
                    max_genuine_ratio[name] = (wratio, (U1, U2, unweighted_ratio))

print("stable grid points:", n_stable)
print("chi0's own bottleneck row (weight-independent structure fact):", bottleneck0)
print()
for name in ["A", "B", "C"]:
    print(f"--- scheme {name} ---")
    print(f"  points where SAME row bottlenecks both chi0,chi^RPA under this weight: {n_same_bottleneck_both[name]}/{n_stable}  (=> weighted ratio == raw per-row ratio, weight irrelevant there)")
    print(f"  points where weighted ratio != unweighted-D ratio: {n_diff_from_unweighted[name]}/{n_stable}")
    print(f"  best 'genuinely different' ratio: {max_genuine_ratio[name][0]:.6f}  at (U1,U2,unweighted_ratio)={max_genuine_ratio[name][1]}")
