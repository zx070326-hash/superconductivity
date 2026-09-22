"""Is the static Re chi0 non-PSD result a theorem rather than a material property?

Team-lead, 2026-09-22, written from scratch.

Claim (analytic, derived independently of any team script)
----------------------------------------------------------
In the sign convention used throughout this week's scripts,

    chi0_{ll'}(q, w=0) = sum_{m,n} c_{mn} * (A_m o B_n)_{ll'},
        A_m = conj(U_k[:,m]) conj(U_k[:,m])^dag        (rank-1 PSD)
        B_n = U_kq[:,n] U_kq[:,n]^dag                  (rank-1 PSD)
        o   = Hadamard (elementwise) product
        c_mn = (f_m - f_n) / (eps_m - eps_n + i*eta)

1. A_m and B_n are rank-1 PSD, so by the SCHUR PRODUCT THEOREM their Hadamard
   product A_m o B_n is PSD.
2. Re c_mn = (f_m - f_n) * D / (D^2 + eta^2) with D = eps_m - eps_n.  The Fermi
   function is strictly decreasing, so (f_m - f_n) and D always have OPPOSITE
   signs => Re c_mn <= 0 for every (m, n).
3. Therefore Re chi0(q,0) is a NON-POSITIVE combination of PSD matrices
   => Re chi0(q,0) is NEGATIVE SEMIDEFINITE, at EVERY q, for EVERY material,
   with no reference to Fermi surface topology, nesting, or iron pnictides.

Consequence: a full-BZ scan reporting "non-PSD at 256/256 q-points" has ZERO
degrees of freedom.  It measured a theorem, not BaFe2As2.

Exception, and it is a second bug: the code's q=0, m=n degenerate branch is
written  f(1-f)/T  (POSITIVE), whereas the correct limit is
    lim_{n->m} (f_m - f_n)/(eps_m - eps_n) = df/deps = -f(1-f)/T   (NEGATIVE).
So that one branch enters with the WRONG SIGN relative to every other term,
which is why q=0 comes out positive-definite while all q != 0 come out NSD.

This script tests 1-3 numerically on random multi-band models (no repo code
imported), and checks the Fermi-derivative sign claim symbolically-by-limit.
"""

import numpy as np

rng = np.random.default_rng(2209221)


def fermi(e, mu, T):
    return 1.0 / (np.exp((e - mu) / T) + 1.0)


def chi0_four_factor(u_k, e_k, u_kq, e_kq, mu, T, eta, degenerate_sign):
    """Correct 4-U-insertion static chi0.

    degenerate_sign = -1 reproduces the true limit df/deps;
    degenerate_sign = +1 reproduces the sign used in the repo scripts.
    """
    nb = len(e_k)
    chi = np.zeros((nb, nb), dtype=complex)
    f_k, f_kq = fermi(e_k, mu, T), fermi(e_kq, mu, T)
    for m in range(nb):
        for n in range(nb):
            d = e_k[m] - e_kq[n]
            if abs(d) < 1e-12:
                c = degenerate_sign * f_k[m] * (1 - f_k[m]) / T
            else:
                c = (f_k[m] - f_kq[n]) / (d + 1j * eta)
            a = np.conj(u_k[:, m])[:, None] * u_k[:, m][None, :]
            b = u_kq[:, n][:, None] * np.conj(u_kq[:, n])[None, :]
            chi += c * a * b
    return chi


def main():
    print("=" * 68)
    print("TEST 1  Schur product theorem step: is A_m o B_n PSD?")
    worst = 0.0
    for _ in range(2000):
        nb = rng.integers(2, 7)
        u = rng.normal(size=nb) + 1j * rng.normal(size=nb)
        v = rng.normal(size=nb) + 1j * rng.normal(size=nb)
        a = np.conj(u)[:, None] * u[None, :]
        b = v[:, None] * np.conj(v)[None, :]
        h = a * b
        h = (h + h.conj().T) / 2
        worst = min(worst, np.linalg.eigvalsh(h).min())
    print(f"  min eigenvalue over 2000 random rank-1 Hadamard products: {worst:.3e}")
    print("  (>= 0 up to roundoff => Schur step holds)")

    print()
    print("=" * 68)
    print("TEST 2  sign of Re c_mn: is (f_m - f_n)*(eps_m - eps_n) <= 0 always?")
    worst = -np.inf
    for _ in range(200000):
        mu, T = rng.normal(), 10 ** rng.uniform(-3, 0)
        e1, e2 = rng.normal(scale=3, size=2)
        prod = (fermi(e1, mu, T) - fermi(e2, mu, T)) * (e1 - e2)
        worst = max(worst, prod)
    print(f"  max of (f_m-f_n)*(eps_m-eps_n) over 2e5 draws: {worst:.3e}")
    print("  (<= 0 => Re c_mn <= 0, since Re c = that product /(D^2+eta^2))")

    print()
    print("=" * 68)
    print("TEST 3  full claim: is static Re chi0 NSD at every q?")
    print("  random multi-band H(k), H(k+q) -- no repo code, no lattice model")
    for nb in (2, 3, 5):
        viol = 0
        worst_pos = -np.inf
        trials = 400
        for _ in range(trials):
            def rh():
                a = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb))
                return (a + a.conj().T) / 2
            e_k, u_k = np.linalg.eigh(rh())
            e_kq, u_kq = np.linalg.eigh(rh())
            mu, T, eta = rng.normal(), 10 ** rng.uniform(-2, -0.5), 1e-3
            c = chi0_four_factor(u_k, e_k, u_kq, e_kq, mu, T, eta, -1)
            re = (c.real + c.real.T) / 2
            ev = np.linalg.eigvalsh(re)
            worst_pos = max(worst_pos, ev.max())
            if ev.max() > 1e-9:
                viol += 1
        print(f"  nb={nb}: PSD-violating-the-NSD-claim {viol}/{trials}, "
              f"largest eigenvalue seen {worst_pos:+.3e}")

    print()
    print("=" * 68)
    print("TEST 4  the degenerate-branch sign bug (q=0, m=n)")
    T = 0.02
    mu = 0.0
    for e in (-0.05, -0.005, 0.0, 0.005, 0.05):
        h = 1e-7
        num = (fermi(e + h, mu, T) - fermi(e - h, mu, T)) / (2 * h)
        repo = +fermi(e, mu, T) * (1 - fermi(e, mu, T)) / T
        true = -fermi(e, mu, T) * (1 - fermi(e, mu, T)) / T
        print(f"  eps={e:+.3f}  d f/d eps (numeric) = {num:+.4f}   "
              f"repo branch = {repo:+.4f}   correct = {true:+.4f}")
    print("  => repo's degenerate branch has the OPPOSITE sign to df/deps")


if __name__ == "__main__":
    main()
