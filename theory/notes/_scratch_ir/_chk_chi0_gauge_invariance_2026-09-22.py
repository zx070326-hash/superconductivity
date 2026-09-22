"""Gauge-invariance test of the multi-orbital chi0 orbital-basis projection.

Independent check (team-lead, 2026-09-22), written from scratch.

Claim under test
----------------
The two-factor form used in the 5-orbital scripts,

    chi0[i,j] += conj(U_k[i,m]) * bubble(m,n) * U_kq[j,n]

is NOT gauge invariant, hence cannot be a physical response function.

Reason (analytic, no numerics needed):
band eigenvectors are defined only up to a phase, U_k[:,m] -> e^{i theta_{k,m}} U_k[:,m].
Under that, the two-factor term picks up e^{-i theta_{k,m}} e^{+i phi_{k+q,n}}, which does
not cancel.  The correct four-factor (cross-term) form

    chi0[i,j] += [U_k[j,m] conj(U_k[i,m])] * bubble(m,n) * [U_kq[i,n] conj(U_kq[j,n])]

pairs each eigenvector with its own conjugate at the SAME (k, band), so every phase
cancels identically.

A gauge-dependent quantity is not an observable: its eigenvalues (and therefore any
"number of negative eigenvalues" count built on it) are an artifact of whatever phase
convention the diagonalizer happened to return.

This script demonstrates both statements numerically on a random 5-band Hermitian model.
"""

import numpy as np

rng = np.random.default_rng(20260922)
NORB = 5


def random_hermitian(n, rng):
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return (a + a.conj().T) / 2


def eigsys(h):
    w, u = np.linalg.eigh(h)
    return w, u


def chi0_two_factor(u_k, u_kq, bubble):
    """The form appearing in _fese_5orbital_FINAL_2026-09-22.py (lines 83-88)."""
    chi = np.zeros((NORB, NORB), dtype=complex)
    for m in range(NORB):
        for n in range(NORB):
            chi += np.conj(u_k[:, m])[:, None] * bubble[m, n] * u_kq[:, n][None, :]
    return chi


def chi0_four_factor(u_k, u_kq, bubble):
    """Correct cross-term form (two_band_chi0_derivation_2026-09-22.md)."""
    chi = np.zeros((NORB, NORB), dtype=complex)
    for m in range(NORB):
        for n in range(NORB):
            # factor1[i,j] = U_k[j,m] * conj(U_k[i,m])
            factor1 = np.conj(u_k[:, m])[:, None] * u_k[:, m][None, :]
            # factor2[i,j] = U_kq[i,n] * conj(U_kq[j,n])
            factor2 = u_kq[:, n][:, None] * np.conj(u_kq[:, n])[None, :]
            chi += factor1 * factor2 * bubble[m, n]
    return chi


def apply_gauge(u, phases):
    return u * np.exp(1j * phases)[None, :]


def main():
    h_k = random_hermitian(NORB, rng)
    h_kq = random_hermitian(NORB, rng)
    _, u_k = eigsys(h_k)
    _, u_kq = eigsys(h_kq)

    # arbitrary real bubble(m,n); its actual value is irrelevant to the gauge argument
    bubble = rng.normal(size=(NORB, NORB))
    bubble = (bubble + bubble.T) / 2

    two_a = chi0_two_factor(u_k, u_kq, bubble)
    four_a = chi0_four_factor(u_k, u_kq, bubble)

    # regauge: independent random phase on every (k, band)
    ph_k = rng.uniform(0, 2 * np.pi, size=NORB)
    ph_kq = rng.uniform(0, 2 * np.pi, size=NORB)
    two_b = chi0_two_factor(apply_gauge(u_k, ph_k), apply_gauge(u_kq, ph_kq), bubble)
    four_b = chi0_four_factor(apply_gauge(u_k, ph_k), apply_gauge(u_kq, ph_kq), bubble)

    d_two = np.abs(two_a - two_b).max()
    d_four = np.abs(four_a - four_b).max()

    print("gauge sensitivity of chi0 (max |Delta| over matrix elements)")
    print(f"  two-factor  form : {d_two:.6e}   <- code in _fese_5orbital_FINAL")
    print(f"  four-factor form : {d_four:.6e}   <- derivation note")

    e_two_a = np.linalg.eigvalsh((two_a + two_a.conj().T) / 2)
    e_two_b = np.linalg.eigvalsh((two_b + two_b.conj().T) / 2)
    e_four_a = np.linalg.eigvalsh((four_a + four_a.conj().T) / 2)
    e_four_b = np.linalg.eigvalsh((four_b + four_b.conj().T) / 2)

    print("\neigenvalues before / after regauge")
    print(f"  two-factor  before : {np.round(e_two_a, 4)}")
    print(f"  two-factor  after  : {np.round(e_two_b, 4)}")
    print(f"    n_neg: {int((e_two_a < 0).sum())} -> {int((e_two_b < 0).sum())}")
    print(f"  four-factor before : {np.round(e_four_a, 4)}")
    print(f"  four-factor after  : {np.round(e_four_b, 4)}")
    print(f"    n_neg: {int((e_four_a < 0).sum())} -> {int((e_four_b < 0).sum())}")

    # how often does the negative-eigenvalue COUNT itself move under pure regauging?
    flips = 0
    trials = 500
    for _ in range(trials):
        p1 = rng.uniform(0, 2 * np.pi, size=NORB)
        p2 = rng.uniform(0, 2 * np.pi, size=NORB)
        t = chi0_two_factor(apply_gauge(u_k, p1), apply_gauge(u_kq, p2), bubble)
        e = np.linalg.eigvalsh((t + t.conj().T) / 2)
        if int((e < 0).sum()) != int((e_two_a < 0).sum()):
            flips += 1
    print(f"\nn_neg changed by pure regauging (two-factor): {flips}/{trials} trials")
    print("  -> a gauge-dependent n_neg is not a property of the material.")


if __name__ == "__main__":
    main()
