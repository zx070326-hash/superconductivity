"""
Is H_N(x) = lambda_1(H_N(sqrt x)) concave in x for every finite compression N?
(If yes for all N, Lemma C would follow from finite-dimensional statements +
 h_N -> h pointwise; a counterexample kills the compression route.)
float64, numerical-evidence.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from lemmaC_compact_probe import O_of


def HN(x, N):
    f = lambda k: x / (k * k + x)
    return np.linalg.eigvalsh(O_of(f, N))[-1]


def second_deriv(x, N, rel=1e-3):
    d = rel * x
    return (HN(x + d, N) - 2 * HN(x, N) + HN(x - d, N)) / (d * d)


if __name__ == "__main__":
    xs = np.exp(np.linspace(np.log(1e-3), np.log(1e4), 120))
    for N in [1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100, 200]:
        vals = np.array([second_deriv(x, N) for x in xs])
        bad = xs[vals > 0]
        if len(bad):
            print(f"N={N:4d}: H_N'' > 0 at {len(bad)} pts, x in [{bad.min():.4g}, {bad.max():.4g}]; max x^2 H'' = {np.max(vals*xs*xs):.3e}")
        else:
            print(f"N={N:4d}: concave on grid; max x^2 H'' = {np.max(vals*xs*xs):.3e}")
