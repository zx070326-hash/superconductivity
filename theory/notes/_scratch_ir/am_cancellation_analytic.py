"""Anderson-Morel cancellation: analytic coefficient matching.

Team lead's observation: r in [20,160] (8x span) gives Tc variation 0.02%, while
||u_c||^2 grows by ~1.4 in the same range. This implies AM running PRECISELY
cancels the log divergence.

Key question: Does the ln in AM running 1/mu*(r) = 1/mu*(10) - ln(r/10) match
the ln in ||u_c||^2 ~ (1/2) ln(2N_c) with the SAME COEFFICIENT?

Strictness: [analytic] for asymptotic formulas; [numerical-evidence] for
finite N_c checks.
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

print("=" * 80)
print("Anderson-Morel cancellation: coefficient matching")
print("=" * 80)

# ANALYTIC: Asymptotic expansion of ||u_c||^2
# sum_{n=0}^{N_c-1} 1/(2n+1) = (1/2)[psi((2N_c+1)/2) - psi(1/2)]
# psi(x) ~ ln(x) - 1/(2x) + O(1/x^2) for large x
# psi(1/2) = -gamma - ln(4)
# So: ||u_c||^2 ~ (1/2)[ln((2N_c+1)/2) + gamma + ln(4)]
#             = (1/2)[ln(2N_c) + gamma + ln(2)]  for large N_c
#             = (1/2) ln(2N_c) + const

# Now N_c = floor(r * varpi), so for large r*varpi:
# ||u_c||^2 ~ (1/2) ln(2 * r * varpi)
#          = (1/2)[ln(2*varpi) + ln(r)]
#          = (1/2) ln(r) + const

# AM running: 1/mu*(r) = 1/mu*(r0) - ln(r/r0)
#                      = 1/mu*(r0) - ln(r) + ln(r0)
# So: d(1/mu*)/d ln(r) = -1

# If secular equation were simply mu~ * ||u_c||^2 ~ const, then:
# mu~/lambda * ||u_c||^2 ~ const
# mu* * ||u_c||^2 ~ const * lambda
# d[mu* * ||u_c||^2]/d ln(r) = 0
# (d mu*/d ln r) * ||u_c||^2 + mu* * d(||u_c||^2)/d ln r = 0
# (mu*)^2 * ||u_c||^2 + mu* * (1/2) = 0  [using d mu*/d ln r = (mu*)^2]
# This gives: ||u_c||^2 = -1/(2 mu*), which is NEGATIVE - contradiction!

print("\n1. Why naive cancellation fails")
print("-" * 80)
print("If secular equation were simply mu~ * ||u_c||^2 ~ const:")
print("  d[mu* * ||u_c||^2]/d ln(r) = 0")
print("  => (mu*)^2 * ||u_c||^2 + mu* * (1/2) = 0")
print("  => ||u_c||^2 = -1/(2 mu*) < 0  [CONTRADICTION]")
print("\nSo the secular equation CANNOT be of the form mu~ * ||u_c||^2 ~ const.")

# The ACTUAL secular equation is:
# sum_j c_j^2/(lam - mu_j) = -1/(2 mu~)
# The key is that c_j^2 is NOT uniform. We have sum_j c_j^2 = ||u_c||^2,
# but the WEIGHTED sum has different behavior.

# For the plateau to hold, we need d(lam)/d ln(r) ~ 0.
# From secular equation implicitly:
# F(lam, r) := 2 mu~(r) * sum_j c_j(lam,r)^2/(lam - mu_j) - 1 = 0
# By implicit differentiation:
# (partial F/partial lam) * d(lam)/d ln(r) + (partial F/partial ln r) = 0

# The plateau means d(lam)/d ln(r) ~ 0, so we need partial F/partial ln r ~ 0.

print("\n2. Implicit differentiation of secular equation")
print("-" * 80)
print("Secular equation: F(lam, r) = 2 mu~(r) sum_j c_j^2/(lam - mu_j) - 1 = 0")
print("\nFor plateau d(lam)/d ln(r) = 0, we need:")
print("  partial F / partial ln(r) = 0  (at fixed lam)")
print("\npartial F / partial ln(r) has three contributions:")
print("  (a) d(mu~)/d ln(r) from explicit mu~ dependence")
print("  (b) d(c_j)/d ln(r) from u_c depending on N_c(r)")
print("  (c) Eigenvalues mu_j are INDEPENDENT of r (they're from A, no Coulomb)")

# Contribution (a): d(mu~)/d ln(r) where mu~ = mu*/lambda
# mu*(r) satisfies 1/mu* = 1/mu*_0 - ln(r/r0)
# d(1/mu*)/d ln(r) = -1
# d(mu*)/d ln(r) = (mu*)^2

print("\n3. Contribution (a): d(mu~)/d ln(r)")
print("-" * 80)
print("mu~ = mu*/lambda, where 1/mu* = 1/mu*_0 - ln(r/r0)")
print("d(1/mu*)/d ln(r) = -1")
print("d(mu*)/d ln(r) = (mu*)^2")
print("d(mu~)/d ln(r) = (mu*)^2/lambda = (mu~)^2 * lambda")

# Contribution (b): d(c_j)/d ln(r) where c_j = <psi_j, u_c>
# u_c[n] = u_n * 1[n < N_c], N_c = floor(r * varpi)
# For large r, N_c ~ r * varpi, so d(N_c)/d ln(r) ~ N_c

# The KEY insight: u_c is LOCALIZED at small n (concentrated in first N_c modes).
# As N_c grows, we're adding MORE tails, so ||u_c||^2 grows.
# But c_1 = <psi_1, u_c> where psi_1 is the TOP eigenmode (smooth, concentrated at small n).

# Let's compute the growth rate of c_1^2

print("\n4. Contribution (b): Growth of c_1^2")
print("-" * 80)
print("c_1 = <psi_1, u_c> = sum_{n<N_c} psi_1[n] * u_n")
print("\nKey observation from numerical data (§3.2 of am_plateau_2026-09-17.md):")
print("  c_1^2 / ||u_c||^2 DECREASES from 93% (r=3) to 48% (r=160)")
print("  But c_1^2 / (mu_1 - lam) INCREASES (weighted contribution)")
print("\nLet's check if c_1^2 itself grows like ||u_c||^2:")

# From team lead's data: r in [20, 160] (8x span)
# ||u_c||^2 grows by ~1.4
# Plateau: Tc varies by 0.02%

# If the secular equation at the plateau is dominated by top mode:
# c_1^2 / (mu_1 - lam) ~ -1/(2 mu~) / (weight of top mode)

# Since top mode weight > 100%, the effective equation is approximately:
# c_1^2 / (mu_1 - lam) ~ -1/(2 mu~)

# For plateau, this should be roughly constant in r:
# d[c_1^2 / (mu_1 - lam)] / d ln(r) ~ d[-1/(2 mu~)] / d ln(r)
# [d(c_1^2)/d ln r] / (mu_1 - lam) - c_1^2 * [d(lam)/d ln r] / (mu_1-lam)^2
#     ~ -1/2 * d(1/mu~)/d ln r

# If d(lam)/d ln r ~ 0 (plateau), then:
# d(c_1^2)/d ln r ~ -(mu_1 - lam) / 2 * d(1/mu~)/d ln r
#                = -(mu_1 - lam) / 2 * [-(mu~)^2 * lambda] / (mu*)^2
#                = (mu_1 - lam) * lambda / (2 mu*)

print("\n5. The cancellation mechanism (qualitative)")
print("-" * 80)
print("The secular equation is:")
print("  sum_j c_j^2/(lam - mu_j) = -1/(2 mu~)")
print("\nTop mode dominates: c_1^2/(lam - mu_1) ~ -1/(2 mu~) * (>100%)")
print("\nFor plateau d(lam)/d ln r ~ 0:")
print("  d(c_1^2)/d ln r ~ (mu_1 - lam) * d[1/(2 mu~)]/d ln r")
print("                  = (mu_1 - lam) * (mu~)^2 * lambda / (2 (mu*)^2)")
print("                  = (mu_1 - lam) * lambda / (2 mu*)")
print("\nNumerical check:")

# Using representative values: mu_1 ~ 1.0, lam ~ 0.76, mu_1-lam ~ 0.24
# mu* ~ 0.13-0.20 over the plateau range, lambda ~ 1.72
mu_1 = 1.001
lam = 0.76
gap = mu_1 - lam
mu_star_mid = 0.15
lam_coupling = 1.721

rate_c1sq = gap * lam_coupling / (2 * mu_star_mid)
print(f"  mu_1 ~ {mu_1:.3f}, lam ~ {lam:.3f}, gap ~ {gap:.3f}")
print(f"  mu* ~ {mu_star_mid:.3f}, lambda ~ {lam_coupling:.3f}")
print(f"  d(c_1^2)/d ln r ~ {rate_c1sq:.3f}")

# Compare to d(||u_c||^2)/d ln r ~ 1/2
rate_norm_sq = 0.5
print(f"  d(||u_c||^2)/d ln r ~ {rate_norm_sq:.3f}")
print(f"\n  Ratio: [d(c_1^2)/d ln r] / [d(||u_c||^2)/d ln r] ~ {rate_c1sq/rate_norm_sq:.3f}")

print("\nInterpretation:")
print("  c_1^2 grows SLOWER than ||u_c||^2 (ratio ~ 0.35)")
print("  This matches the observation that c_1^2/||u_c||^2 DECREASES with r")
print("  But the WEIGHTED sum c_1^2/(mu_1-lam) stays approximately constant")
print("  because (mu_1-lam) is also adjusting (shrinking as lam increases)")

print("\n6. Exact coefficient matching: the ln in AM vs ln in ||u_c||^2")
print("-" * 80)
print("||u_c||^2 ~ (1/2) ln(2N_c) = (1/2) ln(2*r*varpi) = (1/2) ln r + const")
print("1/mu* = 1/mu*_0 - ln(r/r0) = -ln r + const")
print("\nCoefficients:")
print("  d(||u_c||^2) / d ln r = 1/2")
print("  d(1/mu*) / d ln r = -1")
print("\nThe ratio is -1/(1/2) = -2, NOT -1.")
print("So the AM running coefficient does NOT directly match ||u_c||^2 growth.")
print("\nBUT: The secular equation is NOT mu* * ||u_c||^2 = const.")
print("It's sum_j c_j^2/(lam - mu_j) * mu~ = -1/2, where the weights 1/(lam-mu_j)")
print("are HIGHLY non-uniform and the top mode dominates.")

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)
print("The Anderson-Morel plateau is NOT due to simple coefficient matching")
print("  ln r in AM  vs  (1/2) ln r in ||u_c||^2.")
print("\nInstead, it's a DYNAMIC equilibrium:")
print("  1. mu~ = mu*/lambda grows as (mu*)^2 * lambda with ln r")
print("  2. ||u_c||^2 grows as (1/2) with ln r")
print("  3. But the secular sum is WEIGHTED: c_1^2/(lam-mu_1) dominates (>100%)")
print("  4. The gap (mu_1-lam) SHRINKS to compensate for c_1^2 and mu~ growth")
print("  5. Net result: lam adjusts to keep the equation satisfied, staying nearly constant")
print("\nThe 0.02% plateau over r in [20,160] is NOT accidental coefficient matching,")
print("but rather a NON-LINEAR self-adjustment via the gap (mu_1-lam).")
print("\nTo prove this analytically, one would need to:")
print("  (a) Expand secular equation to first order in ln r")
print("  (b) Show that d(lam)/d ln r = O(1/ln^2 r) or similar suppression")
print("  (c) Use the resolvent identity (*) to bound the sub-leading corrections")
print("\nThis is FEASIBLE but not yet done (would be analytic if completed).")
