r"""
theory/solvers/tests/test_lambda_star.py

Unit tests for theory/solvers/lambda_star.py, the productized
lambda_*[alpha2F] stability-bound functional (P3_tc_max_pressure.md Sec.4).

Run: python -m pytest theory/solvers/tests/test_lambda_star.py -q   (from repo root)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import mpmath as mp
import pytest

import lambda_star as ls


# ---------------------------------------------------------------------------
# Hard anchors (P3_tc_max_pressure.md Sec.4, Semenok 2407.12922 Eq.9)
# ---------------------------------------------------------------------------

def test_einstein_delta_anchor():
    """Einstein (delta-function) spectrum -> lambda_* = 1/g_max = 3.6915 exactly,
    the lower bound of lambda_*."""
    val = float(ls.lambda_star_discrete([1.0], [1.0]))
    assert val == pytest.approx(3.6915, abs=1e-4)
    assert val == pytest.approx(float(ls.LAMBDA0), abs=1e-10)


def test_debye_anchor():
    """Debye spectrum: alpha2F(w) ~ w^2 on (0, w_D] => P(w) = 2*a2F/(lambda*w) ~ w
    (NOT alpha2F itself). lambda_* = 4.7186 (Semenok Eq.9 second anchor)."""
    val = float(ls.lambda_star_continuous(lambda w: w, [0, 1]))
    assert val == pytest.approx(4.7186, abs=1e-3)


def test_kappa_and_c3_constants():
    """P3 Sec.4 closed-form curvature coefficients: kappa=5.6035 (sigma^2
    coefficient is kappa/2), c3=0.8959 (sigma^3 x skewness coefficient)."""
    assert float(ls.KAPPA) == pytest.approx(5.6035, abs=1e-3)
    assert float(ls.C3) == pytest.approx(0.8959, abs=1e-3)
    assert float(ls.LAMBDA0) == pytest.approx(3.6915, abs=1e-4)


# ---------------------------------------------------------------------------
# Cross-checks against theory/derivations/_xcheck_p3_kappa_B.py (independent
# reviewer implementation) -- same model spectra, same exact lambda_* values.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("split,expected", [
    (0.10, 3.79501),
    (0.20, 4.10427),
    (0.30, 4.60028),
])
def test_two_peak_matches_reviewer_xcheck(split, expected):
    ws = [1 - split, 1 + split]
    ps = [1.0, 1.0]
    val = float(ls.lambda_star_discrete(ws, ps))
    assert val == pytest.approx(expected, rel=2e-4)


def test_debye_moments_match_analytic():
    """Debye P(w)~w: analytic sigma=sqrt(1/8), r=e^0.5/sqrt(2) (xcheck_B Part 5)."""
    N = 4000
    ws = [(k + 0.5) / N for k in range(N)]
    ps = list(ws)  # P(w) ~ w
    _, sigma, _, r = ls.moments_discrete(ws, ps)
    assert float(sigma) == pytest.approx(float(mp.sqrt(mp.mpf(1) / 8)), rel=2e-3)
    assert float(r) == pytest.approx(float(mp.e**mp.mpf('0.5') / mp.sqrt(2)), rel=2e-3)


# ---------------------------------------------------------------------------
# sigma^2 + c3*skew*sigma^3 expansion: symmetric spectra (good to <0.4% at
# sigma<~0.2) and the asymmetric skew=+-1.5 case from P3 Sec.4 (reviewer
# claim: adding the sigma^3 term drops the error from ~0.9%/1.7% to
# ~0.07%/0.73%).
# ---------------------------------------------------------------------------

def test_expansion_symmetric_two_peak_sigma_0p1():
    ws, ps = [0.9, 1.1], [1.0, 1.0]
    exact = float(ls.lambda_star_discrete(ws, ps))
    _, sigma, skew, _ = ls.moments_discrete(ws, ps)
    approx = float(ls.lambda_star_expansion(sigma, skew))
    rel_err = abs(approx - exact) / exact
    assert rel_err < 0.004  # doc: symmetric spectra, sigma~0.1-0.2 -> <=0.4%


def test_expansion_symmetric_two_peak_sigma_0p2():
    ws, ps = [0.8, 1.2], [1.0, 1.0]
    exact = float(ls.lambda_star_discrete(ws, ps))
    _, sigma, skew, _ = ls.moments_discrete(ws, ps)
    approx = float(ls.lambda_star_expansion(sigma, skew))
    rel_err = abs(approx - exact) / exact
    assert rel_err < 0.004


def test_skewness_cubic_term_improves_asymmetric_double_peak():
    """P3 Sec.4: asymmetric double-peak with skew=+-1.5, sigma=0.2. Adding the
    c3*skew*sigma^3 term should reduce the 2nd-order-only error (~0.9%/1.7%)
    to ~0.07%/0.73% (per team-lead brief, reviewer-confirmed numbers)."""
    cases = [
        ([1.4, 0.9], [0.2, 0.8], 0.0089),   # skew=+1.5, 2nd-order err ~0.89%
        ([0.6, 1.1], [0.2, 0.8], 0.0172),   # skew=-1.5, 2nd-order err ~1.72%
    ]
    expected_3rd_order_err = [0.0007, 0.0073]  # ~0.07%, ~0.73%
    for (ws, ps, err2_expected), err3_expected in zip(cases, expected_3rd_order_err):
        exact = float(ls.lambda_star_discrete(ws, ps))
        _, sigma, skew, _ = ls.moments_discrete(ws, ps)

        approx2 = float(ls.LAMBDA0) * (1 + float(ls.KAPPA)/2 * float(sigma)**2)
        err2 = abs(approx2 - exact) / exact
        assert err2 == pytest.approx(err2_expected, abs=0.002)

        approx3 = float(ls.lambda_star_expansion(sigma, skew))
        err3 = abs(approx3 - exact) / exact
        assert err3 == pytest.approx(err3_expected, abs=0.002)
        assert err3 < err2  # adding the skewness term must improve accuracy


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
