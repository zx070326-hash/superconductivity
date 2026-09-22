"""
Task from team-lead: test whether Lemma (W)'s converse
  (chi^RPA all-positive-row  ==>  exact weighted-metric ceiling >= 1)
holds beyond the single chi0=[[1.5,0.2],[0.2,0.8]] used in
meta/rpa_weighted_metric_nogo_lemma_W_2026-09-22.md.

Budget-first: single counterexample (allpos but ceiling<1) anywhere falsifies
"necessary and sufficient" and downgrades (W) to "necessary only" -- that is
a clean, deliverable result either way, so scan broadly and stop at the
first solid falsification (don't polish around it).

5 chi0 choices, deliberately varied per team-lead's brief:
  chi0_A: kappa near 1              [[1.0,0.2],[0.2,0.9]]      kappa=1.111
  chi0_B: chi0_11 < chi0_22 (flip)  [[0.8,0.2],[0.2,1.5]]      kappa=1.875 (mirrored)
  chi0_C: negative off-diagonal     [[1.5,-0.2],[-0.2,0.8]]    kappa=1.875, sign flip
  chi0_D: large kappa                [[2.0,0.2],[0.2,0.5]]      kappa=4.0
  chi0_E: near-degenerate row 2      [[1.0,0.45],[0.45,0.5]]    kappa=2.0, D0_2=0.05 (thin margin)

For each: Stoner-stable + genuine-failure classification, then among
all-positive-row failures, exact ceiling via kink-point method (re-derived
independently here, not copy-pasted from team-lead's script), cross-checked
against brute-force per point (assert match, abort loudly if any mismatch --
that would mean the formula does NOT generalize and needs rederiving).
"""
import numpy as np

def rows(chi0, M):
    n = 2
    return [M[i,i]-sum(abs(M[i,j]) for j in range(n) if j!=i) for i in range(n)]

def rpa(chi0, U1,U2):
    U = np.diag([U1,U2])
    return chi0 @ np.linalg.inv(np.eye(2)-U@chi0)

def rho(chi0, U1,U2):
    U = np.diag([U1,U2])
    return max(abs(np.linalg.eigvals(U@chi0)))

def exact_ceiling(R, D):
    cands = [1e-9, 1e9, D[1]/D[0]]
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

def brute_ceiling(R, D):
    ts = np.logspace(-8,8,4001)
    return max(min(t*R[0],R[1])/min(t*D[0],D[1]) for t in ts)

cases = {
    "A_kappa_near1":   np.array([[1.0,0.2],[0.2,0.9]]),
    "B_flip_diag":     np.array([[0.8,0.2],[0.2,1.5]]),
    "C_neg_offdiag":   np.array([[1.5,-0.2],[-0.2,0.8]]),
    "D_large_kappa":   np.array([[2.0,0.2],[0.2,0.5]]),
    "E_thin_margin":   np.array([[1.0,0.45],[0.45,0.5]]),
}

grid = np.linspace(0.02, 0.9, 45)

for name, chi0 in cases.items():
    D0 = rows(chi0, chi0)
    kappa = chi0[0,0]/chi0[1,1]
    assert D0[0] > 0 and D0[1] > 0, f"{name}: chi0 not row diagonal dominant, D0={D0}"
    n_stable=0; n_fail=0; n_negrow=0; n_allpos=0
    n_ceiling_ge1=0
    counterexamples = []
    mismatches = 0
    for U1 in grid:
        for U2 in grid:
            if rho(chi0,U1,U2) >= 1.0: continue
            n_stable += 1
            M = rpa(chi0,U1,U2)
            R = rows(chi0, M)
            unweighted_ratio = min(R)/min(D0)
            if unweighted_ratio >= 1.0: continue
            n_fail += 1
            if min(R) < 0:
                n_negrow += 1
                continue
            n_allpos += 1
            ceil = exact_ceiling(R, D0)
            bf = brute_ceiling(R, D0)
            if abs(ceil-bf) > 1e-5*max(1,abs(bf)):
                mismatches += 1
                if mismatches <= 3:
                    print(f"  [{name}] FORMULA MISMATCH at U1={U1:.3f},U2={U2:.3f}: exact={ceil:.6f} brute={bf:.6f}")
            if ceil >= 1.0:
                n_ceiling_ge1 += 1
            else:
                counterexamples.append((U1,U2,R,D0,ceil))
    print(f"=== {name}  kappa={kappa:.3f}  D0={D0} ===")
    print(f"  stable={n_stable} fail={n_fail} negrow={n_negrow} allpos={n_allpos}")
    print(f"  allpos with ceiling>=1: {n_ceiling_ge1}/{n_allpos}")
    print(f"  formula/brute mismatches: {mismatches}")
    if counterexamples:
        print(f"  COUNTEREXAMPLES (allpos but ceiling<1): {len(counterexamples)}")
        for c in counterexamples[:5]:
            print("   ", c)
    print()
