"""
Is "elementwise positivity" basis-dependent, and by how much?

Context: two team deliverables measured the elementwise sign of Im chi0 / Im chi^s
in the ORBITAL basis and drew (hedged) conclusions about Essenberger Eq.(37),
which asserts Im Lambda^SF_{kk'} > 0 with k,k' Fermi-surface (band, momentum)
states.  Lambda^SF = 2 F chi F is a congruence, so the two objects are matrix
elements of related operators in DIFFERENT bases.

This script quantifies how much (if anything) an orbital-basis elementwise sign
census says about the Fermi-surface-basis one, and separately corrects the
team-lead's own earlier power argument.

Q1/Q2 answer the object question in BOTH directions.
Q3/Q4 separate P(data | H_A) from P(data | H_0) -- the distinction the earlier
      "defect B" critique got wrong.
"""
import numpy as np

rng = np.random.default_rng(7)

print("Q1: a PSD matrix that IS elementwise-positive in basis 1 --")
print("    what fraction of its off-diagonals are negative in an unrelated basis?")
N, tot, neg, allpos_rot, trials = 5, 0, 0, 0, 4000
for _ in range(trials):
    # PSD *and* elementwise positive: Gram matrix of nonnegative vectors
    G = rng.random((N, N))
    M = G @ G.T
    assert (M > 0).all() and np.linalg.eigvalsh(M).min() > -1e-12
    Q, _ = np.linalg.qr(rng.standard_normal((N, N)))   # random orthonormal basis
    Mr = Q.T @ M @ Q
    off = Mr[np.triu_indices(N, 1)]
    tot += off.size
    neg += (off < 0).sum()
    allpos_rot += (off > 0).all()
print(f"    negative off-diagonal fraction after rotation : {neg}/{tot} = {100*neg/tot:.1f}%")
print(f"    rotated matrix still elementwise-positive     : {allpos_rot}/{trials} "
      f"= {100*allpos_rot/trials:.2f}%")

print()
print("Q2: reverse -- PSD with NEGATIVE off-diagonals in basis 1;")
print("    can it be elementwise-positive in some other basis?")
hits = 0
for _ in range(4000):
    A = np.array([[1., -.5], [-.5, 1.]])              # PSD, has a negative entry
    th = rng.uniform(0, 2 * np.pi)
    Q = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    hits += (Q.T @ A @ Q)[0, 1] > 0
print(f"    fraction of rotations making it elementwise-positive : "
      f"{hits}/4000 = {100*hits/4000:.1f}%")
print()
print("    => elementwise positivity does NOT transfer between bases in either")
print("       direction, and the no-information baseline is ~50%.")

print()
print("Q3: power against H_A = 'every off-diagonal element > 0' (what Eq.(37) asserts)")
for label, n_neg, n_tot in [("chi0-rebuild 16 groups", 83, 160),
                            ("misc note Im chi^s", 2934, 6400)]:
    print(f"    {label:24s}: {n_neg}/{n_tot} negative  ->  P(data | H_A) = 0 exactly")
print("    A single negative element already falsifies H_A; power = 1, not 0.")
print("    => the earlier 'this design can produce no evidence' verdict was WRONG;")
print("       what had no power was the binary 'was any group all-positive?' summary.")

print()
print("Q4: what the binomial computation actually answered (H_0 = fair coin):")
for label, n_neg, n_tot in [("chi0-rebuild", 83, 160), ("misc note", 2934, 6400)]:
    z = (n_neg - n_tot / 2) / np.sqrt(n_tot * 0.25)
    print(f"    {label:14s}: z = {z:+.3f}  (naive; elements within a q-point are correlated,")
    print(f"    {'':14s}  so the misc note's |z| is inflated and must NOT be read as 6.6 sigma)")
print()
print("    Observed negative fractions: 83/160 = 51.9%, 2934/6400 = 45.8%.")
print("    Both sit on the ~49-50% no-information baseline from Q1/Q2.")
