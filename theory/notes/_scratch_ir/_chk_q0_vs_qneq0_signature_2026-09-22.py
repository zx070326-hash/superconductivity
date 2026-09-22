"""Does the theorem PREDICT the reported q=0 / q=(pi,0) eigenvalue signature?

Team-lead, 2026-09-22.  Companion to _chk_static_chi0_nsd_is_a_theorem_2026-09-22.py.

The team reported, for corrected 5-orbital Graser chi0 (static):
    q = (0,0)    eigenvalues [ 0.005, 0.010, 0.092, 0.703, 2.828]   -> 0 negative
    q = (pi,0)   eigenvalues [-1.796,-1.425,-0.610,-0.073, 0.000]   -> 4 negative (+1 zero)

My theorem says this signature is FORCED and carries no material information:
  * q != 0: all bands non-degenerate => every Re c_mn <= 0 => Re chi0 is NSD.
  * q =  0: the m == n terms hit the code's degenerate branch, which is written
    +f(1-f)/T instead of the correct -f(1-f)/T.  Those terms are O(1/T) and sit on
    the Fermi surface, so they dominate and flip the whole matrix to PD.

Prediction, therefore: with COMPLETELY RANDOM 5-band parameters (nothing to do with
BaFe2As2), the repo branch must give ~0 negative at q=0 and all-non-positive at
q != 0.  If that reproduces, the "256/256 non-PSD" scan measured the convention and
the sign bug, not the material.

Control: with the CORRECT degenerate sign, q = 0 must be NSD as well.

Note on the zero eigenvalue at q=(pi,0): at any q the 5x5 Re chi0 is a sum of
Hadamard products of rank-1 PSD matrices; exact zero modes are generic when the
number of (m,n) pairs carrying weight is limited by Fermi occupation, so "4 negative
+ 1 zero" is NSD with a kernel, not "4 of 5 directions unstable".
"""

import numpy as np

rng = np.random.default_rng(20260922)

NORB = 5


def fermi(e, mu, T):
    return 1.0 / (1.0 + np.exp(np.clip((e - mu) / T, -500, 500)))


def random_hopping():
    """Random real-symmetric 5-orbital tight-binding: H(k) = sum_d t_d cos/sin(k.d)."""
    dvecs = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 0), (0, 2)]
    ts = []
    for _ in dvecs:
        a = rng.normal(scale=0.3, size=(NORB, NORB))
        ts.append((a + a.T) / 2)
    onsite = np.diag(rng.normal(scale=0.5, size=NORB))
    return dvecs, ts, onsite


def hk(kx, ky, dvecs, ts, onsite):
    h = onsite.copy()
    for (dx, dy), t in zip(dvecs, ts):
        h = h + 2 * t * np.cos(kx * dx + ky * dy)
    return h


def chi0_static(q, dvecs, ts, onsite, mu, T, nk, degenerate_sign, eta=1e-3):
    qx, qy = q
    chi = np.zeros((NORB, NORB), dtype=complex)
    ks = 2 * np.pi * np.arange(nk) / nk
    for kx in ks:
        for ky in ks:
            e_k, u_k = np.linalg.eigh(hk(kx, ky, dvecs, ts, onsite))
            e_q, u_q = np.linalg.eigh(hk(kx + qx, ky + qy, dvecs, ts, onsite))
            f_k, f_q = fermi(e_k, mu, T), fermi(e_q, mu, T)
            for m in range(NORB):
                a = np.conj(u_k[:, m])[:, None] * u_k[:, m][None, :]
                for n in range(NORB):
                    d = e_k[m] - e_q[n]
                    if abs(d) < 1e-10:
                        c = degenerate_sign * f_k[m] * (1 - f_k[m]) / T
                    else:
                        c = (f_k[m] - f_q[n]) / (d + 1j * eta)
                    b = u_q[:, n][:, None] * np.conj(u_q[:, n])[None, :]
                    chi += c * a * b
    return chi / nk ** 2


def n_neg(m):
    h = (m + m.conj().T) / 2
    ev = np.linalg.eigvalsh(h.real if np.isrealobj(h) else h)
    return int((ev < -1e-9).sum()), ev


def main():
    T, nk = 0.02, 12
    qs = {"(0,0)": (0.0, 0.0),
          "(pi,0)": (np.pi, 0.0),
          "(pi,pi)": (np.pi, np.pi),
          "(0.37pi,0.81pi)": (0.37 * np.pi, 0.81 * np.pi)}

    print("=" * 74)
    print("RANDOM 5-orbital models -- no BaFe2As2, no Graser parameters at all")
    print("=" * 74)
    hdr = f"{'model':>6} {'q':>17} {'repo (+f(1-f)/T)':>22} {'correct (-f(1-f)/T)':>22}"
    print(hdr)
    summary = {}
    for trial in range(4):
        dvecs, ts, onsite = random_hopping()
        e_all = np.concatenate([np.linalg.eigvalsh(hk(kx, ky, dvecs, ts, onsite))
                                for kx in 2 * np.pi * np.arange(nk) / nk
                                for ky in 2 * np.pi * np.arange(nk) / nk])
        mu = np.median(e_all)  # put the chemical potential in the band manifold
        for label, q in qs.items():
            nb_repo, ev_repo = n_neg(chi0_static(q, dvecs, ts, onsite, mu, T, nk, +1))
            nb_corr, ev_corr = n_neg(chi0_static(q, dvecs, ts, onsite, mu, T, nk, -1))
            print(f"{trial:>6} {label:>17} {nb_repo:>13} neg / 5 {nb_corr:>14} neg / 5")
            summary.setdefault(label, []).append((nb_repo, nb_corr))
            if trial == 0:
                print(f"{'':>6} {'  repo ev':>17} "
                      + " ".join(f"{v:+8.3f}" for v in ev_repo))
                print(f"{'':>6} {'  correct ev':>17} "
                      + " ".join(f"{v:+8.3f}" for v in ev_corr))

    print()
    print("=" * 74)
    print("VERDICT")
    for label, rows in summary.items():
        repo = sorted({r for r, _ in rows})
        corr = sorted({c for _, c in rows})
        print(f"  q={label:>17}: repo n_neg in {repo},  correct n_neg in {corr}")
    print()
    print("  If repo gives 0 negative at q=(0,0) and ~5 negative elsewhere on RANDOM")
    print("  parameters, the reported Graser signature is parameter-free: it is the")
    print("  sign convention plus the degenerate-branch sign bug, not the material.")


if __name__ == "__main__":
    main()
