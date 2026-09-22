"""Real-symmetric (production-relevant) branch of the chi0 gauge test.

Graser's tight-binding H(k) is real symmetric at the k-points used, so numpy's eigh
returns REAL eigenvectors and the residual gauge freedom is only a sign s_m = +/-1 per
band, not a continuous phase.  The question this script answers:

  does the two-factor form's negative-eigenvalue count still move under the residual
  +/-1 gauge?

If yes, the 5-orbital n_neg numbers are an artifact of whatever signs LAPACK happened
to return -- deterministic, reproducible, and physically meaningless.

Complement to _chk_chi0_gauge_invariance_2026-09-22.py (complex/U(1) branch).
"""

import itertools

import numpy as np

rng = np.random.default_rng(920)
NORB = 5


def random_real_symmetric(n, rng):
    a = rng.normal(size=(n, n))
    return (a + a.T) / 2


def chi0_two_factor(u_k, u_kq, bubble):
    chi = np.zeros((NORB, NORB))
    for m in range(NORB):
        for n in range(NORB):
            chi += u_k[:, m][:, None] * bubble[m, n] * u_kq[:, n][None, :]
    return chi


def chi0_four_factor(u_k, u_kq, bubble):
    chi = np.zeros((NORB, NORB))
    for m in range(NORB):
        for n in range(NORB):
            factor1 = u_k[:, m][:, None] * u_k[:, m][None, :]
            factor2 = u_kq[:, n][:, None] * u_kq[:, n][None, :]
            chi += factor1 * factor2 * bubble[m, n]
    return chi


def n_neg(m):
    return int((np.linalg.eigvalsh((m + m.T) / 2) < 0).sum())


def main():
    _, u_k = np.linalg.eigh(random_real_symmetric(NORB, rng))
    _, u_kq = np.linalg.eigh(random_real_symmetric(NORB, rng))
    bubble = rng.normal(size=(NORB, NORB))
    bubble = (bubble + bubble.T) / 2

    base_two = n_neg(chi0_two_factor(u_k, u_kq, bubble))
    base_four = n_neg(chi0_four_factor(u_k, u_kq, bubble))
    print(f"LAPACK-gauge n_neg   two-factor={base_two}   four-factor={base_four}")

    # exhaust all 2^5 x 2^5 = 1024 residual sign gauges
    seen_two, seen_four = {}, {}
    for sk in itertools.product([1, -1], repeat=NORB):
        for skq in itertools.product([1, -1], repeat=NORB):
            a = u_k * np.array(sk)[None, :]
            b = u_kq * np.array(skq)[None, :]
            v2 = n_neg(chi0_two_factor(a, b, bubble))
            v4 = n_neg(chi0_four_factor(a, b, bubble))
            seen_two[v2] = seen_two.get(v2, 0) + 1
            seen_four[v4] = seen_four.get(v4, 0) + 1

    total = 2 ** (2 * NORB)
    print(f"\nover all {total} residual sign gauges (exhaustive):")
    print(f"  two-factor  n_neg distribution : {dict(sorted(seen_two.items()))}")
    print(f"  four-factor n_neg distribution : {dict(sorted(seen_four.items()))}")

    frac = 1 - seen_two.get(base_two, 0) / total
    print(f"\n  two-factor: n_neg differs from the LAPACK value in {frac:.1%} of gauges")
    print(f"  four-factor: single value = gauge invariant"
          if len(seen_four) == 1 else "  four-factor: NOT invariant (unexpected!)")

    # repeat over many random (k, k+q) pairs to show it is not a single lucky draw
    moved = 0
    trials = 200
    for _ in range(trials):
        _, a0 = np.linalg.eigh(random_real_symmetric(NORB, rng))
        _, b0 = np.linalg.eigh(random_real_symmetric(NORB, rng))
        bb = rng.normal(size=(NORB, NORB))
        bb = (bb + bb.T) / 2
        base = n_neg(chi0_two_factor(a0, b0, bb))
        hit = False
        for _ in range(20):
            sk = rng.choice([1, -1], size=NORB)
            skq = rng.choice([1, -1], size=NORB)
            if n_neg(chi0_two_factor(a0 * sk[None, :], b0 * skq[None, :], bb)) != base:
                hit = True
                break
        moved += hit
    print(f"\n  n_neg movable by a sign gauge: {moved}/{trials} random (k,k+q) draws")


if __name__ == "__main__":
    main()
