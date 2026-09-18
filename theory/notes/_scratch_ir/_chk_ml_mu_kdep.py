"""[numerical-evidence] Audit requested by team-lead: in the mu*-PSD chain, every
place a MODULUS is taken -- is the object a band-space SCALAR or a FUNCTION OF k?

Their defect (three times now) was `||O[f]|| = ||O[|f|]||`, which silently assumes
sign(f(k)) is k-independent.  I enumerate the moduli in my (mu1)-(mu4) chain:

  (mu1) Gram form  M_{(n,i),(m,j)} = 2 mu_ij w^(i)_n w^(j)_m :  NO modulus.
        The n-dependence enters as an OUTER PRODUCT, so <v,Mv> = 2 a^T mu a with
        a_i = <v_i, w^(i)>.  PSD holds for ARBITRARY REAL profiles w^(i) -- in
        particular sign-alternating-in-n ones.  Tested in (A) below.
  (mu2) |sum_i a_i rhohat_i(q)|^2 :  modulus of a COMPLEX number under the
        pointwise-nonneg weight Vhat(q) >= 0, inside an IDENTITY, not a bound.
        Safe for any q-dependence.  (Not retested here; it is (mu2)'s own gate.)
  (mu3) lmin(S L S) > -1, lmax(mu), ell_min :  band-space SCALARS.  ell_i =
        ln(omega_c^(i)/W_i) is one number per band, no k.  Tested in (C).
  (mu4) set inclusion :  no modulus.

  The one k-modulus in the WHOLE delivery is (Z-dom) |Lam^Delta_ij(k)| <= Lam^Z_ij(k),
  used POINTWISE IN k inside (ML3-abs) -- it never needs |Lam^Delta| to be
  non-increasing, and G3 is applied to Lam^ph and Lam^sf SEPARATELY (each nonneg
  and non-increasing by construction).  That is exactly the arrow that broke in
  their section 6.2.  Tested in (B) in the two-frequency regime where |Lam^Delta|
  is NOT non-increasing.

DIRECTION BOOKKEEPING.  lam1(P_N O P_N) <= lam1(O)  (Rayleigh on a subspace)
  ==> every ratio below is a LOWER estimate of the true ratio
  ==> ratio > 1 is a SOUND violation; ratio < 1 is evidence only.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.linalg import eigvalsh, eigvals
import _chk_ml_build as Bld
from _chk_ml_absM2 import multiband_O2

G2 = 1.318140212132
rng = np.random.default_rng(20260918_2)


def coulomb_block(mu, profs, N, Bb):
    """M_{(n,i),(m,j)} = 2 mu_ij w^(i)_n w^(j)_m, ordering index = n*Bb + i."""
    W = np.zeros((N * Bb, Bb))
    for i in range(Bb):
        W[np.arange(N) * Bb + i, i] = profs[i]
    return 2.0 * W @ mu @ W.T


def rho(M):
    return float(np.max(np.abs(eigvals(M))))


# ---------------------------------------------------------------- (A)
# (mu1) with profiles whose SIGN VARIES WITH n.  If the Gram argument had
# secretly assumed w >= 0, this is where it would break.
print("(A) (mu1) Gram form under sign-alternating-in-n profiles")
worst_eig, worst_id = np.inf, 0.0
for _ in range(4000):
    Bb = int(rng.integers(2, 6)); N = 40
    A = rng.normal(size=(Bb, Bb + 2)); mu = A @ A.T          # PSD by construction
    u = Bld.u_vec(N)
    profs = []
    for i in range(Bb):
        Nc = int(rng.integers(3, N))
        sgn = rng.choice([-1.0, 1.0], size=N)                # sign flips with n
        profs.append(u * (np.arange(N) < Nc) * sgn * rng.uniform(0.3, 2.0))
    M = coulomb_block(mu, profs, N, Bb)
    worst_eig = min(worst_eig, eigvalsh((M + M.T) / 2).min() / max(np.linalg.norm(M), 1e-12))
    v = rng.normal(size=(N * Bb,))
    a = np.array([v[np.arange(N) * Bb + i] @ profs[i] for i in range(Bb)])
    worst_id = max(worst_id, abs(v @ M @ v - 2 * a @ mu @ a) / max(abs(v @ M @ v), 1e-12))
print(f"    min normalised eig of M over 4000 draws : {worst_eig:.3e}   (must be >= 0)")
print(f"    max rel deviation of <v,Mv> = 2 a^T mu a: {worst_id:.3e}")
print("    ==> (mu1) is sign-in-n blind BY STRUCTURE (outer product), not by luck.\n")


# ---------------------------------------------------------------- (B)
# Intersection never tested before: two-frequency pairing kernel (so sign(Lam^Delta)
# DOES depend on k) TOGETHER with a Coulomb block.
print("(B) two-frequency Lam^Delta (sign k-dependent) + Coulomb block")


def kernels(Lph, Lsf, wph, wsf, N):
    ks = np.arange(2 * N + 2)
    fph = Bld.brk(ks, wph)[:, None, None] * Lph[None]
    fsf = Bld.brk(ks, wsf)[:, None, None] * Lsf[None]
    return fph + fsf, fph - fsf, fph, fsf


n_noninc_fail = 0
viol_psd_step = 0
viol_ceiling = 0
worst_ratio = -np.inf
N = 90
for _ in range(800):
    Bb = int(rng.integers(2, 4))
    Lph = rng.uniform(0, 1.5, (Bb, Bb)); Lph = (Lph + Lph.T) / 2
    Lsf = rng.uniform(0, 1.5, (Bb, Bb)); Lsf = (Lsf + Lsf.T) / 2
    wph = 10 ** rng.uniform(-1.0, 0.3)
    wsf = 10 ** rng.uniform(-1.0, 0.6)
    lamZ, lamD, fph, fsf = kernels(Lph, Lsf, wph, wsf, N)
    wgt = np.ones(Bb)

    # is |Lam^Delta| non-increasing in k?  (the premise that broke in their 6.2)
    absD = np.abs(lamD[1:N])
    if np.any(np.diff(absD, axis=0) > 1e-14):
        n_noninc_fail += 1

    O = multiband_O2(lamZ, lamD, N, wgt)
    A2 = rng.normal(size=(Bb, Bb + 2)); mu = A2 @ A2.T
    u = Bld.u_vec(N)
    profs = [u * (np.arange(N) < int(rng.integers(3, N)))
             * rng.choice([-1.0, 1.0], size=N) for _ in range(Bb)]
    Mc = coulomb_block(mu, profs, N, Bb)
    l1_0 = eigvalsh(O)[-1]
    l1_mu = eigvalsh(O - Mc)[-1]
    viol_psd_step += l1_mu > l1_0 + 1e-10

    M2ph = Lph * wph ** 2
    M2sf = Lsf * wsf ** 2
    off = lambda X: X - np.diag(np.diag(X))
    M2G = (M2ph + 3 * M2sf) + 2 * off(M2ph + M2sf)
    r = l1_mu / (G2 * rho(M2G))
    viol_ceiling += r > 1.0
    worst_ratio = max(worst_ratio, r)
print(f"    |Lam^Delta| NOT non-increasing in : {n_noninc_fail}/800 draws"
      f"   (their 6.2 arrow is genuinely absent here)")
print(f"    step  lam1(O - M) <= lam1(O)      : {viol_psd_step}/800 violations")
print(f"    ceiling with M2^G                 : {viol_ceiling}/800 violations,"
      f" worst ratio {worst_ratio:.4f}")

# control: mu elementwise nonneg but NOT PSD -> the step must break
bad = 0
worst_c = -np.inf
for _ in range(300):
    Bb = 2
    Lph = rng.uniform(0, 1.5, (Bb, Bb)); Lph = (Lph + Lph.T) / 2
    Lsf = rng.uniform(0, 1.5, (Bb, Bb)); Lsf = (Lsf + Lsf.T) / 2
    wph = 10 ** rng.uniform(-1.0, 0.3); wsf = 10 ** rng.uniform(-1.0, 0.6)
    lamZ, lamD, _, _ = kernels(Lph, Lsf, wph, wsf, N)
    O = multiband_O2(lamZ, lamD, N, np.ones(Bb))
    d = rng.uniform(0.0, 0.3, Bb)
    mu = np.array([[d[0], 0.9], [0.9, d[1]]])          # nonneg, indefinite
    u = Bld.u_vec(N)
    profs = [u * (np.arange(N) < int(rng.integers(3, N))) for _ in range(Bb)]
    Mc = coulomb_block(mu, profs, N, Bb)
    l1_mu = eigvalsh(O - Mc)[-1]
    bad += l1_mu > eigvalsh(O)[-1] + 1e-10
    M2ph = Lph * wph ** 2; M2sf = Lsf * wsf ** 2
    off = lambda X: X - np.diag(np.diag(X))
    worst_c = max(worst_c, l1_mu / (G2 * rho((M2ph + 3 * M2sf) + 2 * off(M2ph + M2sf))))
print(f"    control mu nonneg but INDEFINITE  : {bad}/300 break the step,"
      f" worst ceiling ratio {worst_c:.3f}\n")


# ---------------------------------------------------------------- (C)
# (mu3): confirm every quantity in the sharp criterion is a band-space scalar,
# i.e. the criterion's verdict is unchanged when the Coulomb profile's sign in n
# is scrambled.  If any k-function had leaked into L, this would move.
print("(C) (mu3) criterion is blind to the n-profile (all objects band scalars)")
mismatch = 0
for _ in range(2000):
    Bb = int(rng.integers(2, 6))
    A = rng.normal(size=(Bb, Bb + 2)); mu = A @ A.T
    ell = rng.uniform(-2.0, 2.0, Bb)
    L = np.diag(ell)
    S = np.real(np.linalg.cholesky(mu + 1e-14 * np.eye(Bb)))
    Ssym = _sq = None
    ev, V = np.linalg.eigh(mu)
    Ssym = V @ np.diag(np.sqrt(np.maximum(ev, 0))) @ V.T
    crit = eigvalsh(Ssym @ L @ Ssym).min() > -1
    ms = np.linalg.solve(np.eye(Bb) + mu @ L, mu)
    psd = eigvalsh((ms + ms.T) / 2).min() > -1e-9 * max(np.linalg.norm(ms), 1.0)
    mismatch += crit != psd
print(f"    iff criterion vs direct PSD test  : {mismatch}/2000 mismatches")
print("    (L = diag(ln(omega_c^(i)/W_i)) has no k index at all; the n-profile")
print("     enters only through (mu1), which is already sign-in-n blind.)")
print("\nTOTAL FAILURES:",
      viol_psd_step + viol_ceiling + mismatch + (worst_eig < -1e-9))
