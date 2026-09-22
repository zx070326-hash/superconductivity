"""
Lemma (W) converse -- ANALYTIC PROOF + from-scratch numerical verification.

Claim (converse half of the dichotomy, previously numerical-evidence 3026/3026):
    chi0 = [[a,b],[b,c]] real symmetric, Delta := ac-b^2 > 0, rows strictly
    diagonally dominant (D1 := a-|b| > 0, D2 := c-|b| > 0);
    U = diag(U1,U2) with U1,U2 >= 0; chi = (chi0^{-1} - U)^{-1} Stoner-stable
    (positive definite).  If BOTH chi-rows are diagonally dominant
    (R1,R2 > 0) then

        ceiling := sup_{t>0} min(t R1, R2) / min(t D1, D2)  >=  1.

PROOF.
    chi0^{-1} = (1/Delta)[[c,-b],[-b,a]].  Write
        p  := c/Delta - U1,   q  := a/Delta - U2,   r := b/Delta,
        p0 := c/Delta,        q0 := a/Delta,        s := p q - r^2 > 0,
        s0 := p0 q0 - r^2 = (ac-b^2)/Delta^2 = 1/Delta.
    Then chi = (1/s)[[q, r],[r, p]], so
        R1 = (q-|r|)/s,  R2 = (p-|r|)/s,
        D1 = (q0-|r|)/s0, D2 = (p0-|r|)/s0   (U=0 specialization).
    The t->0 and t->infinity limits of ratio(t) give
        ceiling >= max(R1/D1, R2/D2),
    so it suffices to show max(R1/D1, R2/D2) >= 1.
    Put rho := |r|, x := q-rho, x0 := q0-rho, y := p-rho, y0 := p0-rho
    (all > 0 by hypothesis), u := x/x0 in (0,1], v := y/y0 in (0,1]
    (<=1 because U1,U2 >= 0).  Since
        s  = pq - rho^2 = xy + rho(x+y),
        s0 = x0 y0 + rho(x0+y0),
    and R1/D1 = u * s0/s,  R2/D2 = v * s0/s,  WLOG u >= v; then
        R1/D1 >= 1
      <=> u (x0 y0 + rho(x0+y0)) >= x y + rho(x+y)
      <=> u x0 y0 (1-v) + rho y0 (u-v) >= 0,
    which holds termwise since v <= 1 and u >= v.  QED
    Equality iff u = v = 1, i.e. U = 0.

So the converse is ANALYTIC for N=2, diagonal U >= 0 -- not numerical-evidence.
Combined with Lemma (W) itself (analytic), the dichotomy
    unsalvageable  <=>  some row deficit of chi^RPA is negative
is now proved in BOTH directions for this class.  Route 2 closes.

This script checks the proof, not the conclusion: it (i) verifies each algebraic
step symbolically-by-random-sampling, (ii) brute-force-scans ceiling on a dense
t-grid and compares with max(R1/D1,R2/D2), (iii) runs a NEGATIVE CONTROL with
U having a negative entry (proof hypothesis U>=0 violated) to confirm the
hypothesis is load-bearing rather than decorative.
"""
import numpy as np

rng = np.random.default_rng(20260922)


def build(a, b, c, U1, U2):
    chi0 = np.array([[a, b], [b, c]], float)
    Delta = a * c - b * b
    chi = np.linalg.inv(np.linalg.inv(chi0) - np.diag([U1, U2]))
    return chi0, chi, Delta


def deficits(M):
    return np.array([M[0, 0] - abs(M[0, 1]), M[1, 1] - abs(M[1, 0])])


def ceiling_brute(R, D, n=200001):
    t = np.logspace(-9, 9, n)
    num = np.minimum(t * R[0], R[1])
    den = np.minimum(t * D[0], D[1])
    return np.max(num / den)


def sample(allow_negative_U=False):
    while True:
        a, c = rng.uniform(0.2, 3.0, 2)
        b = rng.uniform(-1.0, 1.0)
        if a - abs(b) <= 0 or c - abs(b) <= 0 or a * c - b * b <= 0:
            continue
        lo = -0.8 if allow_negative_U else 0.0
        U1, U2 = rng.uniform(lo, 0.9, 2)
        chi0, chi, Delta = build(a, b, c, U1, U2)
        if not np.all(np.linalg.eigvalsh(chi) > 0):
            continue                      # Stoner-unstable / wrong branch
        return a, b, c, U1, U2, chi0, chi, Delta


print("=" * 72)
print("STEP CHECK 1: the p,q,r,s reparametrization reproduces chi exactly")
print("=" * 72)
worst = 0.0
for _ in range(20000):
    a, b, c, U1, U2, chi0, chi, Delta = sample()
    p, q, r = c / Delta - U1, a / Delta - U2, b / Delta
    s = p * q - r * r
    chi_pred = np.array([[q, r], [r, p]]) / s
    worst = max(worst, np.abs(chi - chi_pred).max())
print(f"  max |chi - (1/s)[[q,r],[r,p]]|  = {worst:.3e}   (expect ~1e-14)")

print()
print("=" * 72)
print("STEP CHECK 2: R1/D1 = u*s0/s and R2/D2 = v*s0/s  (the key rewriting)")
print("=" * 72)
worst_u, worst_v = 0.0, 0.0
for _ in range(20000):
    a, b, c, U1, U2, chi0, chi, Delta = sample()
    p, q, r = c / Delta - U1, a / Delta - U2, b / Delta
    p0, q0 = c / Delta, a / Delta
    s, s0 = p * q - r * r, p0 * q0 - r * r
    rho = abs(r)
    x, x0, y, y0 = q - rho, q0 - rho, p - rho, p0 - rho
    R, D = deficits(chi), deficits(chi0)
    if not (R > 0).all():
        continue
    worst_u = max(worst_u, abs(R[0] / D[0] - (x / x0) * s0 / s))
    worst_v = max(worst_v, abs(R[1] / D[1] - (y / y0) * s0 / s))
print(f"  max |R1/D1 - u*s0/s| = {worst_u:.3e}")
print(f"  max |R2/D2 - v*s0/s| = {worst_v:.3e}   (expect ~1e-13)")

print()
print("=" * 72)
print("STEP CHECK 3: the termwise inequality  u*x0*y0*(1-v) + rho*y0*(u-v) >= 0")
print("=" * 72)
worst_term, n3 = np.inf, 0
for _ in range(50000):
    a, b, c, U1, U2, chi0, chi, Delta = sample()
    p, q, r = c / Delta - U1, a / Delta - U2, b / Delta
    p0, q0 = c / Delta, a / Delta
    rho = abs(r)
    x, x0, y, y0 = q - rho, q0 - rho, p - rho, p0 - rho
    if x <= 0 or y <= 0:
        continue
    u, v = x / x0, y / y0
    if u < v:
        u, v, x0, y0 = v, u, y0, x0       # WLOG swap
    worst_term = min(worst_term, u * x0 * y0 * (1 - v) + rho * y0 * (u - v))
    n3 += 1
print(f"  n = {n3};  min of termwise expression = {worst_term:.6e}  (must be >= 0)")

print()
print("=" * 72)
print("STEP CHECK 4: ceiling >= 1 whenever R1,R2 > 0  (the theorem itself)")
print("   and  ceiling >= max(R1/D1, R2/D2) (the t->0 / t->inf limits)")
print("=" * 72)
n_pos, n_fail_ceiling, n_fail_max, worst_ceiling, worst_max = 0, 0, 0, np.inf, np.inf
for _ in range(4000):
    a, b, c, U1, U2, chi0, chi, Delta = sample()
    R, D = deficits(chi), deficits(chi0)
    if not (R > 0).all():
        continue
    n_pos += 1
    mx = max(R[0] / D[0], R[1] / D[1])
    cb = ceiling_brute(R, D)
    worst_max = min(worst_max, mx)
    worst_ceiling = min(worst_ceiling, cb)
    n_fail_max += mx < 1 - 1e-12
    n_fail_ceiling += cb < 1 - 1e-12
print(f"  all-positive-row samples          : {n_pos}")
print(f"  max(R1/D1,R2/D2) < 1 violations   : {n_fail_max}      min value = {worst_max:.6f}")
print(f"  brute-force ceiling < 1 violations: {n_fail_ceiling}      min value = {worst_ceiling:.6f}")

print()
print("=" * 72)
print("NEGATIVE CONTROL: drop the hypothesis U >= 0 (allow a negative U entry).")
print("If the theorem still held, the hypothesis would be decorative.")
print("=" * 72)
n_pos_nc, n_fail_nc, worst_nc = 0, 0, np.inf
for _ in range(8000):
    a, b, c, U1, U2, chi0, chi, Delta = sample(allow_negative_U=True)
    if U1 >= 0 and U2 >= 0:
        continue
    R, D = deficits(chi), deficits(chi0)
    if not (R > 0).all():
        continue
    n_pos_nc += 1
    mx = max(R[0] / D[0], R[1] / D[1])
    worst_nc = min(worst_nc, mx)
    n_fail_nc += mx < 1 - 1e-12
print(f"  all-positive-row samples with some U_a < 0 : {n_pos_nc}")
print(f"  max(R1/D1,R2/D2) < 1 violations            : {n_fail_nc}"
      f"   ({100.0 * n_fail_nc / max(n_pos_nc, 1):.1f}%)   min = {worst_nc:.6f}")
print()
print("Interpretation: violations HERE and none ABOVE  ==>  U >= 0 is load-bearing,")
print("so Step-check 4's 0 violations are not vacuous.")

print()
print("=" * 72)
print("REPRODUCTION of rpa-deep-dive's five families (independent implementation)")
print("=" * 72)
fams = {
    "orig [[1.5,.2],[.2,.8]]": (1.5, 0.2, 0.8),
    "A    [[1.0,.2],[.2,.9]]": (1.0, 0.2, 0.9),
    "B    [[.8,.2],[.2,1.5]]": (0.8, 0.2, 1.5),
    "C    [[1.5,-.2],[-.2,.8]]": (1.5, -0.2, 0.8),
    "D    [[2.0,.2],[.2,.5]]": (2.0, 0.2, 0.5),
    "E    [[1.0,.45],[.45,.5]]": (1.0, 0.45, 0.5),
}
grid = np.linspace(0.02, 0.8, 60)
tot_pos = tot_bad = 0
for name, (a, b, c) in fams.items():
    n_stable = n_fail = n_neg = n_allpos = n_ceil_ge1 = 0
    for U1 in grid:
        for U2 in grid:
            chi0, chi, Delta = build(a, b, c, U1, U2)
            if not np.all(np.linalg.eigvalsh(chi) > 0):
                continue
            n_stable += 1
            R, D = deficits(chi), deficits(chi0)
            if ceiling_brute(R, D, 4001) >= 1 and min(R[0] / D[0], R[1] / D[1]) >= 1:
                pass
            std = min(R) / min(D)
            if std >= 1:
                continue
            n_fail += 1
            if (R < 0).any():
                n_neg += 1
            else:
                n_allpos += 1
                if max(R[0] / D[0], R[1] / D[1]) >= 1:
                    n_ceil_ge1 += 1
    tot_pos += n_allpos
    tot_bad += n_allpos - n_ceil_ge1
    print(f"  {name:26s} stable={n_stable:5d} fail={n_fail:5d} "
          f"negrow={n_neg:5d} allpos={n_allpos:5d} ceil>=1={n_ceil_ge1:5d}")
print(f"\n  TOTAL all-positive-row failure points : {tot_pos}")
print(f"  counterexamples to the converse       : {tot_bad}")
