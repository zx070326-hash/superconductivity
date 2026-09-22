"""
team-lead, 2026-09-22.  Does the (ML3-abs) cell REALLY need  Lambda^sf_ij >= 0
elementwise (= Essenberger Eq.(37)), or does a strictly weaker pair of
conditions already deliver  L(k) >= 0 ?

Context:  ceiling/sf_psd_route_fails_signless_laplacian_2026-09-22.md (Team B)
established that replacing "elementwise nonneg" by "band-space PSD" ALONE fails
(46.6% of PSD inputs give a non-PSD signless Laplacian).  That is correct.
But its section 2 SPLITS

      L := diag(R^Z) - Lam^Delta  =  [diag(R^ph) - Lam^ph]  +  [diag(R^sf) + Lam^sf]
                                  =:        L^ph            +          M

and then demands  M >= 0  on its own.  The split DISCARDS any surplus in L^ph.
It is justified only in the pure-SF limit (Lam^ph = 0), which is one entry of
the coverage census -- not the FeSe verdict, where Lam^ph != 0.

UNSPLIT route tested here:
    Lam^Delta <= Lam^Z  (PSD order)  =>  L >= diag(R^Z) - Lam^Z =: L^Z
    and   x'L^Z x = 1/2 sum_{i!=j} Lam^Z_ij ( sqrt(w_j/w_i) x_i - sqrt(w_i/w_j) x_j )^2
    => L^Z >= 0  whenever the OFF-DIAGONAL entries of Lam^Z are >= 0
       (diagonal entries drop out of the identity entirely).

With Lam^Z = Lam^ph + Lam^sf, Lam^Delta = Lam^ph - Lam^sf:
    (i)  Lam^Z - Lam^Delta = 2 Lam^sf >= 0   (PSD)      <- FREE: congruence /
         Sylvester inertia, since Im Lam^SF = 2 f^xc (Im chi) f^xc with Im chi >= 0
    (ii) (Lam^ph + Lam^sf)_ij >= 0  for i != j          <- weaker than Lam^sf_ij >= 0

OLD condition (what Eq.(37) is used for):  Lam^sf_ij >= 0 for all i,j.
(i)+(ii) is implied by OLD (given Lam^ph >= 0) but not conversely.

Checks below:
  A. the Laplacian identity, numerically
  B. (i)+(ii) => L >= 0 : violation count must be 0
  C. negative control 1: drop (i)            -> violations must appear
  D. negative control 2: drop (ii)           -> must reproduce Team B's ~46.6%
  E. GAIN: fraction of samples satisfying NEW but violating OLD (if ~0 the
     relaxation is empty in practice and buys nothing)
  F. Team B's hand counterexample, and how much phonon weight rescues it

Run:  OMP_NUM_THREADS=1 python _chk_unsplit_LZ_route_relaxes_eq37_2026-09-22.py
"""
import numpy as np

rng_global = np.random.default_rng(20260922)


def sym(A):
    return (A + A.T) / 2.0


def R_of(Lam, w):
    """R_i = (Lam w)_i / w_i."""
    return (Lam @ w) / w


def L_of(LamZ, LamD, w):
    """L := diag(R^Z) - Lam^Delta."""
    return np.diag(R_of(LamZ, w)) - LamD


def min_eig(A):
    return float(np.linalg.eigvalsh(sym(A)).min())


def offdiag_min(A):
    B = A.copy()
    np.fill_diagonal(B, np.inf)
    return float(B.min())


# ---------------------------------------------------------------- A. identity
def check_identity(n_trials=4000):
    print("=== A. Laplacian identity for L^Z = diag(R^Z) - Lam^Z ===")
    worst = 0.0
    for _ in range(n_trials):
        B = rng_global.integers(2, 7)
        LamZ = sym(rng_global.normal(size=(B, B)))      # ARBITRARY sign: identity
        w = rng_global.uniform(0.2, 3.0, size=B)        # is an identity, not a bound
        x = rng_global.normal(size=B)
        lhs = float(x @ (np.diag(R_of(LamZ, w)) - LamZ) @ x)
        rhs = 0.0
        for i in range(B):
            for j in range(B):
                if i == j:
                    continue
                rhs += 0.5 * LamZ[i, j] * (np.sqrt(w[j] / w[i]) * x[i]
                                           - np.sqrt(w[i] / w[j]) * x[j]) ** 2
        worst = max(worst, abs(lhs - rhs) / max(1.0, abs(lhs)))
    print(f"  max relative residual over {n_trials} random (B,Lam,w,x): {worst:.3e}")
    print("  -> diagonal entries of Lam^Z drop out; only i!=j entries enter.\n")


# -------------------------------------------------------- sampler for B--E
def sample(mode, B, rng):
    """Return (Lam_ph, Lam_sf, w).  Lam_ph is always elementwise >= 0 (phonon
    spectral function), which is the uncontroversial input."""
    w = rng.uniform(0.3, 2.5, size=B)
    # phonon: elementwise nonneg, symmetric
    scale_ph = rng.uniform(0.0, 3.0)
    Lam_ph = sym(rng.uniform(0.0, 1.0, size=(B, B))) * scale_ph

    if mode in ("new", "drop_ii"):
        # Lam^sf PSD, off-diagonals free in sign
        S = rng.normal(size=(B, B))
        Lam_sf = S @ S.T
        Lam_sf *= rng.uniform(0.1, 2.0) / max(1e-12, np.abs(Lam_sf).max())
    elif mode == "drop_i":
        # Lam^sf indefinite (NOT PSD), but we will filter on (ii) only
        Lam_sf = sym(rng.normal(size=(B, B)))
        Lam_sf *= rng.uniform(0.1, 2.0) / max(1e-12, np.abs(Lam_sf).max())
    else:
        raise ValueError(mode)
    return Lam_ph, Lam_sf, w


def cond_i(Lam_sf):
    return min_eig(Lam_sf) >= -1e-12


def cond_ii(Lam_ph, Lam_sf):
    return offdiag_min(Lam_ph + Lam_sf) >= -1e-12


def cond_old(Lam_sf):
    return float(Lam_sf.min()) >= -1e-12          # Eq.(37): ALL entries >= 0


def run(mode, target, n_keep, seed):
    """Draw until n_keep samples satisfy `target`, report L>=0 violation rate."""
    rng = np.random.default_rng(seed)
    kept = viol = 0
    worst = np.inf
    worst_norm = np.inf
    drawn = 0
    while kept < n_keep and drawn < 4_000_000:
        drawn += 1
        B = int(rng.integers(2, 7))
        Lam_ph, Lam_sf, w = sample(mode, B, rng)
        if not target(Lam_ph, Lam_sf):
            continue
        kept += 1
        LamZ, LamD = Lam_ph + Lam_sf, Lam_ph - Lam_sf
        m = min_eig(L_of(LamZ, LamD, w))
        if m < -1e-9:
            viol += 1
            worst = min(worst, m)
            scale = max(1e-12, np.abs(LamZ).max())
            worst_norm = min(worst_norm, m / scale)
    return kept, viol, worst, worst_norm, drawn


def main():
    check_identity()

    print("=== B. NEW conditions (i) Lam^sf PSD  +  (ii) offdiag(Lam^ph+Lam^sf) >= 0 ===")
    k, v, wst, wn, d = run("new", lambda p, s: cond_i(s) and cond_ii(p, s), 20000, 11)
    print(f"  kept {k} (from {d} draws);  L not PSD: {v}  ({100.0*v/k:.3f}%)")
    print(f"  worst min-eig: {wst if v else 0.0}\n")

    print("=== C. negative control 1: drop (i) -- Lam^sf NOT PSD, (ii) still held ===")
    k, v, wst, wn, d = run("drop_i", lambda p, s: (not cond_i(s)) and cond_ii(p, s), 20000, 12)
    print(f"  kept {k} (from {d} draws);  L not PSD: {v}  ({100.0*v/k:.3f}%)")
    print(f"  worst min-eig {wst:.4g}, worst normalised {wn:.4g}\n")

    print("=== D. negative control 2: drop (ii) -- Lam^sf PSD, some offdiag(Lam^Z) < 0 ===")
    k, v, wst, wn, d = run("drop_ii", lambda p, s: cond_i(s) and (not cond_ii(p, s)), 20000, 13)
    print(f"  kept {k} (from {d} draws);  L not PSD: {v}  ({100.0*v/k:.3f}%)")
    print(f"  worst min-eig {wst:.4g}, worst normalised {wn:.4g}")
    print("  (Team B reported 46.6% / -2736 for the Lam^ph=0 slice of this.)\n")

    print("=== E. GAIN: does NEW ever hold where OLD (Eq.37) fails? ===")
    rng = np.random.default_rng(14)
    n_new = n_old = n_new_not_old = n_old_not_new = 0
    N = 200000
    for _ in range(N):
        B = int(rng.integers(2, 7))
        Lam_ph, Lam_sf, w = sample("new", B, rng)
        new = cond_i(Lam_sf) and cond_ii(Lam_ph, Lam_sf)
        old = cond_old(Lam_sf)
        n_new += new
        n_old += old
        n_new_not_old += (new and not old)
        n_old_not_new += (old and not new)
    print(f"  over {N} draws (Lam^sf PSD by construction):")
    print(f"    OLD  (Lam^sf_ij >= 0 for all i,j) holds : {n_old:7d}  ({100.0*n_old/N:.3f}%)")
    print(f"    NEW  ((i) and (ii))              holds : {n_new:7d}  ({100.0*n_new/N:.3f}%)")
    print(f"    NEW and not OLD  (the actual gain)     : {n_new_not_old:7d}  ({100.0*n_new_not_old/N:.3f}%)")
    print(f"    OLD and not NEW  (must be 0: OLD=>NEW) : {n_old_not_new:7d}\n")

    print("=== F. Team B's hand counterexample, and the phonon weight that rescues it ===")
    Lam_sf = np.array([[4.0, -1.9], [-1.9, 1.0]])
    w = np.array([1.0, 1.0])
    print(f"  Lam^sf eig = {np.linalg.eigvalsh(Lam_sf)}  (PSD, off-diagonal -1.9 < 0)")
    M = np.diag(R_of(Lam_sf, w)) + Lam_sf
    print(f"  Team B's split object M = diag(R^sf)+Lam^sf : eig = {np.linalg.eigvalsh(M)}"
          f"  det = {np.linalg.det(M):+.4f}   -> NOT PSD, reproduced")
    print("  unsplit route, adding a uniform off-diagonal phonon coupling p:")
    for p in (0.0, 1.0, 1.8, 1.9, 1.95, 2.5, 4.0):
        Lam_ph = np.array([[0.0, p], [p, 0.0]])
        LamZ, LamD = Lam_ph + Lam_sf, Lam_ph - Lam_sf
        ok_ii = offdiag_min(LamZ) >= -1e-12
        print(f"    p={p:4.2f}:  offdiag(Lam^Z) = {LamZ[0,1]:+.3f}  (ii){'OK ' if ok_ii else 'X  '}"
              f"  min eig L = {min_eig(L_of(LamZ, LamD, w)):+.4f}"
              f"   min eig L^Z = {min_eig(np.diag(R_of(LamZ, w)) - LamZ):+.4f}")
    print("  -> (ii) flips exactly at p = 1.9, and L >= 0 from there on;")
    print("     Team B's split demands M>=0, which p can never fix (M is p-free).")


if __name__ == "__main__":
    main()
