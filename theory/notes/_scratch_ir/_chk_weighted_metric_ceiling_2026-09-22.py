"""Ceiling for weighted-metric rescue in the 2x2 RPA diagonal-dominance problem.

Decision-relevant question (team-lead, 2026-09-22)
-------------------------------------------------
rpa-deep-dive reports scheme B (w ~ 1/U_alpha) rescues 20 of 210 genuine-failure
grid points.  But their own weight-free scan notes that in 2x2 a weight choice can
essentially only select which row the min() picks.  If so, 210 is the wrong
denominator: the right one is "at how many points does ANY weight rescue?"

Analytic ceiling -- FIRST ATTEMPT WAS WRONG, kept here with its refutation
--------------------------------------------------------------------------
  D_w[M] = min(w1 R1, w2 R2) with R_alpha the row deficit of M.
  With the SAME weights on both sides, put t = w1/w2:
      ratio(t) = min(t R1, R2) / min(t D1, D2)

  WRONG first guess:  sup_t ratio(t) = max(R1/D1, R2/D2).
  It is valid only when BOTH R1 > 0 and R2 > 0.  If the RPA matrix has lost
  diagonal dominance in row 1 (R1 < 0 < R2), then min(t R1, R2) = t R1 for every
  t > 0 -- the numerator is pinned to the negative row and can never select row 2,
  so the value R2/D2 is NOT attainable and the formula over-estimates the ceiling.
  Measured gap on the grid below: up to 8.1e3.  Using that inflated ceiling to say
  "scheme B only reaches 12.9% of it" would be a one-sided over-estimate used to
  support a FAIL verdict -- exactly the direction the one-sided-bound rule forbids.

  CORRECT: ratio(t) is piecewise monotone with breakpoints only where a min()
  switches, i.e. at t = D2/D1 and (when R1 > 0) t = R2/R1.  So the exact sup is the
  max over {t -> 0+, the breakpoints, t -> inf}.  That is what ceiling_exact() does,
  and it is cross-checked against a brute-force sweep.

So the honest score for any weight formula is (rescued) / (ceiling-rescuable),
not (rescued) / (genuine failures) -- with the ceiling computed exactly.

Reproduces rpa-deep-dive's setup: chi0 = [[1.5,0.2],[0.2,0.8]] (kappa=1.875),
U = diag(U1,U2) on a grid, Stoner-stable points only.
"""

import numpy as np

CHI0 = np.array([[1.5, 0.2], [0.2, 0.8]])


def row_deficits(m):
    """R_alpha = M_aa - sum_{b != a} |M_ab| for a 2x2 matrix."""
    return np.array([m[0, 0] - abs(m[0, 1]), m[1, 1] - abs(m[1, 0])])


def rpa(chi0, u_diag):
    u = np.diag(u_diag)
    return chi0 @ np.linalg.inv(np.eye(2) - u @ chi0)


def rho(chi0, u_diag):
    return max(abs(np.linalg.eigvals(chi0 @ np.diag(u_diag))))


def ratio_with_weights(d_bare, d_rpa, w):
    return min(w[0] * d_rpa[0], w[1] * d_rpa[1]) / min(w[0] * d_bare[0], w[1] * d_bare[1])


def ceiling_exact(d_bare, d_rpa):
    """Exact sup over all positive weights of ratio(t), t = w1/w2.

    ratio(t) = min(t R1, R2) / min(t D1, D2) is piecewise monotone; every kink is a
    min() switch, at t = D2/D1 and (only if R1 > 0) t = R2/R1.  Evaluating at the
    kinks plus the two limits therefore gives the exact supremum.
    """
    d1, d2 = d_bare
    r1, r2 = d_rpa
    cands = [d2 / d1]
    if r1 > 0:
        cands.append(r2 / r1)
    vals = [ratio_with_weights(d_bare, d_rpa, np.array([t, 1.0]))
            for t in cands if t > 0]
    # limits t -> 0+ and t -> inf (evaluate just inside, both branches are constant
    # or monotone there, so the endpoint value is the limit)
    tiny = min([t for t in cands if t > 0], default=1.0) * 1e-9
    huge = max([t for t in cands if t > 0], default=1.0) * 1e9
    vals.append(ratio_with_weights(d_bare, d_rpa, np.array([tiny, 1.0])))
    vals.append(ratio_with_weights(d_bare, d_rpa, np.array([huge, 1.0])))
    return max(vals)


def main():
    grid = np.linspace(0.02, 0.8, 60)
    n_stable = n_genuine = 0
    n_ceiling = n_naive_ceiling = n_B = n_A = n_C = 0
    ceiling_but_not_B = []
    naive_inflation = 0.0

    # brute-force t sweep, to confirm the analytic sup claim rather than assume it
    t_sweep = np.concatenate([np.logspace(-6, 6, 4001)])

    max_gap = 0.0
    for u1 in grid:
        for u2 in grid:
            uv = np.array([u1, u2])
            if rho(CHI0, uv) >= 1.0:
                continue
            n_stable += 1
            m = rpa(CHI0, uv)
            d0 = row_deficits(CHI0)
            d1 = row_deficits(m)
            if min(d0) <= 0:
                continue  # chi0 itself not diagonally dominant; criterion undefined
            unweighted = min(d1) / min(d0)
            if unweighted >= 1.0:
                continue  # not a genuine failure
            n_genuine += 1

            per_row = np.array([d1[0] / d0[0], d1[1] / d0[1]])
            naive_sup = per_row.max()          # the refuted formula
            exact_sup = ceiling_exact(d0, d1)  # the correct ceiling
            brute_sup = max(ratio_with_weights(d0, d1, np.array([t, 1.0]))
                            for t in t_sweep)
            max_gap = max(max_gap, abs(exact_sup - brute_sup))
            naive_inflation = max(naive_inflation, naive_sup - exact_sup)

            if exact_sup >= 1.0:
                n_ceiling += 1
            if naive_sup >= 1.0:
                n_naive_ceiling += 1

            wA = 1.0 / np.diag(CHI0)
            wB = 1.0 / uv
            wC = uv.copy()
            rA = ratio_with_weights(d0, d1, wA)
            rB = ratio_with_weights(d0, d1, wB)
            rC = ratio_with_weights(d0, d1, wC)
            n_A += rA >= 1.0
            n_B += rB >= 1.0
            n_C += rC >= 1.0
            if exact_sup >= 1.0 and rB < 1.0:
                ceiling_but_not_B.append((u1, u2, rB, exact_sup))

    print(f"Stoner-stable grid points      : {n_stable}")
    print(f"genuine failures (unweighted<1): {n_genuine}")
    print()
    print(f"max |exact ceiling - brute sweep| = {max_gap:.3e}   (0 => exact)")
    print(f"max over-estimate of the REFUTED formula max(R1/D1,R2/D2)"
          f" = {naive_inflation:.3e}")
    print()
    print(f"CEILING (exact, any weight)    : {n_ceiling}  "
          f"({n_ceiling / n_genuine:.1%} of genuine failures)")
    print(f"  [refuted formula would claim : {n_naive_ceiling}]")
    print(f"  scheme A (1/chi0_aa) rescues : {n_A}")
    print(f"  scheme B (1/U_a)     rescues : {n_B}")
    print(f"  scheme C (U_a)       rescues : {n_C}")
    if n_ceiling:
        print(f"  -> scheme B reaches {n_B / n_ceiling:.1%} of the exact ceiling")
    print()
    if n_ceiling and n_B == n_ceiling:
        print("scheme B attains the ceiling: no weight formula can do better.")
    else:
        print(f"gap: {len(ceiling_but_not_B)} points are rescuable by SOME weight "
              f"but not by scheme B; e.g.")
        for row in ceiling_but_not_B[:5]:
            print(f"    U1={row[0]:.3f} U2={row[1]:.3f}  rB={row[2]:.4f}  "
                  f"ceiling={row[3]:.4f}")


if __name__ == "__main__":
    main()
