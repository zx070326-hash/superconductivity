"""Are "256/256 non-PSD" and "RPA improve/worsen 0/0" BOTH forced theorems?

Team-lead, 2026-09-22.  Third in the series; companions:
  _chk_static_chi0_nsd_is_a_theorem_2026-09-22.py     (static Re chi0 is NSD)
  _chk_q0_vs_qneq0_signature_2026-09-22.py            (q=0 is the sign-bug point)

Chain of claims, each tested below.

(C1) Static Re chi0 (4-U-insertion form, code's sign convention) is negative
     definite at every q != 0.                       [established in companion 1]

(C2) If chi0 is negative definite and U is PSD, then
         chi^s = chi0 (1 - U chi0)^{-1} = (chi0^{-1} - U)^{-1}
     is also negative definite.
     Proof: chi0^{-1} is ND; x*(chi0^{-1} - U)x = x*chi0^{-1}x - x*Ux < 0 since
     x*Ux >= 0.  So chi0^{-1} - U is ND, hence so is its inverse.  ANALYTIC.
     => n_neg(chi^s) = n_neg(chi0) = N at EVERY q, for EVERY U PSD.
     => a scan reporting "RPA improved 0, worsened 0 out of 256" measured this
        identity, not the interplay of interaction and nesting.

(C3) Worse: eigenvalues of U chi0 equal those of U^{1/2} chi0 U^{1/2}, which is NSD,
     so every eigenvalue of (1 - U chi0) is >= 1.  The RPA denominator can NEVER
     vanish.  The Stoner instability is UNREACHABLE in this sign convention, for any
     PSD U and any interaction strength.
     => the physical susceptibility is chi0^phys = -chi0^code; the "non-PSD finding"
        is a missing overall minus sign, and the object that should have been tested
        for PSD-ness is the negated one (which is then PSD by C1 -- also vacuous).

(C4) Is Graser's U matrix actually PSD?  If not, C2/C3 need the extra hypothesis
     stated, so this is checked against the standard Hubbard-Kanamori orbital
     vertex rather than assumed.
"""

import numpy as np

rng = np.random.default_rng(922003)


def nd_check(m):
    ev = np.linalg.eigvalsh((m + m.conj().T) / 2)
    return ev.max()


def main():
    print("=" * 72)
    print("(C2) chi0 ND + U PSD  =>  chi^s ND.  Random test, N = 2..6")
    worst = -np.inf
    n_bad = 0
    trials = 4000
    for _ in range(trials):
        n = rng.integers(2, 7)
        a = rng.normal(size=(n, n))
        chi0 = -(a @ a.T) - 1e-3 * np.eye(n)          # negative definite
        b = rng.normal(size=(n, n)) * rng.uniform(0.1, 5)
        u = b @ b.T                                    # PSD
        chis = chi0 @ np.linalg.inv(np.eye(n) - u @ chi0)
        top = nd_check(chis)
        worst = max(worst, top)
        if top > -1e-12:
            n_bad += 1
    print(f"  chi^s largest eigenvalue over {trials} trials: {worst:+.3e}")
    print(f"  trials where chi^s was NOT negative definite: {n_bad}/{trials}")

    print()
    print("=" * 72)
    print("(C3) can the RPA denominator (1 - U chi0) ever be singular?")
    worst_min = np.inf
    for _ in range(4000):
        n = rng.integers(2, 7)
        a = rng.normal(size=(n, n))
        chi0 = -(a @ a.T) - 1e-3 * np.eye(n)
        b = rng.normal(size=(n, n)) * rng.uniform(0.1, 50)   # up to huge U
        u = b @ b.T
        ev = np.linalg.eigvals(np.eye(n) - u @ chi0)
        worst_min = min(worst_min, ev.real.min())
    print(f"  min over 4000 trials of min Re eig(1 - U chi0): {worst_min:.6f}")
    print("  (>= 1 => denominator never vanishes => no Stoner instability, ever)")

    print()
    print("=" * 72)
    print("(C4) is the Hubbard-Kanamori orbital interaction matrix PSD?")
    print("  spin channel, orbital-diagonal basis: U_aa = U, U_ab = J (a != b)")
    print("  (the form used for chi^s = chi0 (1 - U chi0)^-1 in 5-orbital RPA)")
    for U, J in [(1.0, 0.25), (1.2, 0.30), (0.8, 0.0), (1.0, 0.5), (1.0, -0.2)]:
        for n in (2, 5):
            m = np.full((n, n), J)
            np.fill_diagonal(m, U)
            ev = np.linalg.eigvalsh(m)
            tag = "PSD" if ev.min() >= -1e-12 else "NOT PSD"
            print(f"    N={n}  U={U:+.2f} J={J:+.2f}: eig min={ev.min():+.4f} "
                  f"max={ev.max():+.4f}  -> {tag}")
    print("  (U > 0, 0 <= J <= U  =>  M = (U-J)I + J*ones  =>  eigenvalues")
    print("   U-J (multiplicity N-1) and U+(N-1)J, both >= 0  => PSD.  ANALYTIC.)")

    print()
    print("=" * 72)
    print("(C2') negative control: does the claim need U PSD?  try indefinite U")
    n_flip = 0
    for _ in range(4000):
        n = rng.integers(2, 7)
        a = rng.normal(size=(n, n))
        chi0 = -(a @ a.T) - 1e-3 * np.eye(n)
        b = rng.normal(size=(n, n))
        u = (b + b.T) / 2                              # indefinite
        d = np.eye(n) - u @ chi0
        if abs(np.linalg.det(d)) < 1e-9:
            continue
        chis = chi0 @ np.linalg.inv(d)
        if nd_check(chis) > 1e-12:
            n_flip += 1
    print(f"  chi^s acquired a positive eigenvalue in {n_flip}/4000 trials")
    print("  (> 0 => the PSD hypothesis on U is load-bearing, not decoration)")


if __name__ == "__main__":
    main()
