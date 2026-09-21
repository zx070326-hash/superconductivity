"""Independent re-derivation of the mu* extension of the multiband ceiling.

Level: the lemmas below are analytic; everything this script prints is a CONSISTENCY
CHECK (float64, finite N), never load-bearing.

Re-derived from zero (NOT copied from mustar_psd_multiband_ceiling_2026-09-18.md):

 (G)  Gram block.  The multiband Coulomb block in the symmetrized basis is
        M_{(n,i),(m,j)} = 2 mu_ij w^{(i)}_n w^{(j)}_m
      with a per-band window profile w^{(i)}.  For any v,
        <v, M v> = 2 sum_ij mu_ij a_i a_j = 2 a^T mu a,  a_i := <v_i, w^{(i)}>.
      So mu PSD  ==>  M PSD.  No Kronecker structure used.

 (P)  Parseval, not Bochner.  mu_ij = <rho_i, V rho_j> with V a convolution kernel:
        a^T mu a = int dq Vhat(q) |sum_i a_i rhohat_i(q)|^2 >= 0  whenever Vhat >= 0.
      No square root, no sign condition on rho_i.

 (A)  Anderson-Morel with a MATRIX L.  mu* = (I + mu L)^{-1} mu.  Then
        mu* = S (I + S L S)^{-1} S,   S := mu^{1/2},
      because S(I + S L S) = (I + S S L) S.  Hence L symmetric PSD (NOT necessarily
      scalar, NOT necessarily commuting with mu) already gives mu* symmetric PSD.
      This is strictly weaker than "L is a scalar": band-dependent cutoffs give
      L = diag(ln(E_F/omega_c^(i))) > 0, which is PSD.  Counterexample hunt below
      shows the hypothesis is not vacuous: L indefinite CAN break PSD.
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.linalg import eigvalsh, inv, norm

import _chk_ml_build as B

rng = np.random.default_rng(20260918)
FAIL = 0


def report(tag, ok, detail=""):
    global FAIL
    if not ok:
        FAIL += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {tag}  {detail}")


def coul_block(mu, profiles, N):
    """M_{(n,i),(m,j)} = 2 mu_ij w^(i)_n w^(j)_m, ordering index = n*Bb + i."""
    Bb = mu.shape[0]
    W = np.asarray(profiles)  # (Bb, N)
    M = np.zeros((N, Bb, N, Bb))
    for i in range(Bb):
        for j in range(Bb):
            M[:, i, :, j] = 2.0 * mu[i, j] * np.outer(W[i], W[j])
    return M.reshape(N * Bb, N * Bb)


def window(N, Nc):
    u = B.u_vec(N)
    return np.where(np.arange(N) < Nc, u, 0.0)


print("=" * 72)
print("(G) Gram block: PSD of mu  ==>  PSD of the Coulomb block, per-band windows")
print("=" * 72)
N = 60
worst_eig = np.inf
worst_id = 0.0
for trial in range(4000):
    Bb = int(rng.integers(2, 6))
    A = rng.normal(size=(Bb, Bb))
    mu = A @ A.T                                   # PSD by construction
    Ncs = rng.integers(1, N, size=Bb)              # DELIBERATELY unequal per band
    profs = [window(N, Ncs[i]) for i in range(Bb)]
    if rng.random() < 0.3:                         # also non-window profiles, signed
        profs = [rng.normal(size=N) for _ in range(Bb)]
    M = coul_block(mu, profs, N)
    worst_eig = min(worst_eig, eigvalsh(M)[0])
    # identity check
    v = rng.normal(size=N * Bb)
    a = np.array([v.reshape(N, Bb)[:, i] @ profs[i] for i in range(Bb)])
    worst_id = max(worst_id, abs(v @ M @ v - 2.0 * a @ mu @ a) / (abs(v @ M @ v) + 1e-12))
report("min eig of Coulomb block over 4000 draws", worst_eig > -1e-11, f"min={worst_eig:.3e}")
report("identity <v,Mv> = 2 a^T mu a", worst_id < 1e-10, f"max rel dev={worst_id:.3e}")

# control: mu elementwise nonneg but NOT PSD
bad = 0.0
for trial in range(2000):
    Bb = 2
    m11, m22 = rng.uniform(0.02, 0.4, size=2)
    m12 = rng.uniform(np.sqrt(m11 * m22), 3.0)     # violates m11 m22 >= m12^2
    mu = np.array([[m11, m12], [m12, m22]])
    Ncs = rng.integers(1, N, size=2)
    M = coul_block(mu, [window(N, Ncs[0]), window(N, Ncs[1])], N)
    bad = min(bad, eigvalsh(M)[0])
report("control: elementwise-nonneg non-PSD mu gives an indefinite block",
       bad < -1e-6, f"min eig={bad:.4f}")

print()
print("=" * 72)
print("(G') ceiling with the Coulomb block, per-band windows:  lam1(OO - M) <= h*lmax")
print("=" * 72)
Nc_test = 80
worst_psd = -np.inf
worst_nonpsd = -np.inf
for trial in range(400):
    Bb = 2
    a_, b_ = rng.uniform(0.05, 3.0, size=2)
    c_ = rng.uniform(0.0, np.sqrt(a_ * b_))
    Lt = np.array([[a_, c_], [c_, b_]])
    varpi = float(rng.uniform(0.1, 8.0))
    w = np.array([1.0, float(rng.uniform(0.2, 5.0))])
    lam = B.lam_factorized(Lt, varpi, Nc_test)
    OO = B.multiband_O(lam, Nc_test, w=w)
    denom = B.h_of(varpi, Nc_test) * eigvalsh(Lt)[-1]
    Ncs = rng.integers(1, Nc_test, size=2)
    profs = [window(Nc_test, Ncs[0]), window(Nc_test, Ncs[1])]
    # PSD mu
    Am = rng.normal(size=(2, 2)) * 0.3
    mu = Am @ Am.T
    worst_psd = max(worst_psd, B.lam1(OO - coul_block(mu, profs, Nc_test)) / denom)
    # elementwise nonneg, non-PSD mu
    m11, m22 = rng.uniform(0.02, 0.3, size=2)
    m12 = rng.uniform(np.sqrt(m11 * m22), 2.0)
    mu2 = np.array([[m11, m12], [m12, m22]])
    worst_nonpsd = max(worst_nonpsd, B.lam1(OO - coul_block(mu2, profs, Nc_test)) / denom)
report("PSD mu, UNEQUAL per-band cutoffs: max ratio <= 1", worst_psd <= 1.0 + 1e-9,
       f"max={worst_psd:.9f}")
report("control: non-PSD mu breaks the ceiling", worst_nonpsd > 1.0,
       f"max={worst_nonpsd:.3f}")

print()
print("=" * 72)
print("(P) Parseval: mu_ij from a convolution kernel with Vhat >= 0 is PSD")
print("=" * 72)
ng = 64
q = np.fft.fftfreq(ng) * ng
QX, QY = np.meshgrid(q, q, indexing="ij")
q2 = QX ** 2 + QY ** 2
worst = np.inf
worst_par = 0.0
for trial in range(500):
    ks2 = float(rng.uniform(0.5, 40.0))
    Vhat = 1.0 / (q2 + ks2)                        # screened Coulomb, strictly > 0
    Bb = int(rng.integers(2, 6))
    rho = rng.normal(size=(Bb, ng, ng))            # SIGNED on purpose
    rh = np.fft.fft2(rho, axes=(1, 2))
    mu = np.real(np.einsum("iqr,qr,jqr->ij", rh.conj(), Vhat, rh)) / ng ** 2
    mu = 0.5 * (mu + mu.T)
    e = eigvalsh(mu)
    worst = min(worst, e[0] / max(e[-1], 1e-300))
    aa = rng.normal(size=Bb)
    lhs = aa @ mu @ aa
    rhs = np.sum(Vhat * np.abs(np.einsum("i,iqr->qr", aa, rh)) ** 2) / ng ** 2
    worst_par = max(worst_par, abs(lhs - rhs) / (abs(lhs) + 1e-12))
report("Vhat>0, signed rho: min normalized eig >= 0", worst > -1e-12, f"min={worst:.3e}")
report("Parseval identity a^T mu a = int Vhat |sum a_i rhohat_i|^2", worst_par < 1e-10,
       f"max rel dev={worst_par:.3e}")

# control: overscreening (Vhat < 0 somewhere) can break PSD
neg = 0.0
for trial in range(500):
    Vhat = 1.0 / (q2 + 2.0) - 0.9 * np.exp(-(q2 - 9.0) ** 2 / 4.0)   # dips negative
    rho = rng.normal(size=(3, ng, ng))
    rh = np.fft.fft2(rho, axes=(1, 2))
    mu = np.real(np.einsum("iqr,qr,jqr->ij", rh.conj(), Vhat, rh)) / ng ** 2
    mu = 0.5 * (mu + mu.T)
    e = eigvalsh(mu)
    neg = min(neg, e[0] / max(abs(e[-1]), 1e-300))
report("control: Vhat with a negative window can give non-PSD mu", neg < -1e-6,
       f"min normalized eig={neg:.4f}  (Vhat_min={Vhat.min():.3f})")

print()
print("=" * 72)
print("(A) Anderson-Morel with MATRIX L:  L PSD suffices (scalar L not needed)")
print("=" * 72)


def am(mu, L):
    return inv(np.eye(mu.shape[0]) + mu @ L) @ mu


worst_sym = 0.0
worst_eig = np.inf
worst_comm = 0.0
worst_fact = 0.0
for trial in range(4000):
    Bb = int(rng.integers(2, 6))
    A = rng.normal(size=(Bb, Bb))
    mu = A @ A.T
    # L symmetric PSD but NOT scalar and NOT commuting with mu
    if rng.random() < 0.5:
        L = np.diag(rng.uniform(0.05, 8.0, size=Bb))          # band-dependent cutoffs
    else:
        C = rng.normal(size=(Bb, Bb))
        L = C @ C.T
    ms = am(mu, L)
    worst_sym = max(worst_sym, norm(ms - ms.T) / max(norm(ms), 1e-300))
    e = eigvalsh(0.5 * (ms + ms.T))
    worst_eig = min(worst_eig, e[0] / max(e[-1], 1e-300))
    worst_comm = max(worst_comm, norm(mu @ L - L @ mu) / max(norm(mu) * norm(L), 1e-300))
    # factorization mu* = S (I + S L S)^{-1} S
    ev, U = np.linalg.eigh(mu)
    S = U @ np.diag(np.sqrt(np.maximum(ev, 0.0))) @ U.T
    alt = S @ inv(np.eye(Bb) + S @ L @ S) @ S
    worst_fact = max(worst_fact, norm(alt - ms) / max(norm(ms), 1e-300))
report("mu* symmetric for symmetric L", worst_sym < 1e-9, f"max asym={worst_sym:.3e}")
report("mu* PSD for L PSD (non-scalar, non-commuting)", worst_eig > -1e-11,
       f"min normalized eig={worst_eig:.3e}; max ||[mu,L]||={worst_comm:.3f}")
report("factorization mu* = S(I+SLS)^{-1}S", worst_fact < 1e-9,
       f"max rel dev={worst_fact:.3e}")

# is the hypothesis L >= 0 vacuous?  hunt for L indefinite breaking PSD
hit = 0
worst_bad = 0.0
for trial in range(20000):
    Bb = 2
    A = rng.normal(size=(Bb, Bb))
    mu = A @ A.T
    C = rng.normal(size=(Bb, Bb))
    L = 0.5 * (C + C.T)                                # symmetric, INDEFINITE allowed
    if eigvalsh(L)[0] >= 0:
        continue
    Mm = np.eye(Bb) + mu @ L
    if abs(np.linalg.det(Mm)) < 1e-6:
        continue
    ms = am(mu, L)
    e = eigvalsh(0.5 * (ms + ms.T))
    r = e[0] / max(abs(e[-1]), 1e-300)
    if r < -1e-8:
        hit += 1
        worst_bad = min(worst_bad, r)
report("hypothesis L PSD is NOT vacuous: indefinite L breaks PSD", hit > 0,
       f"{hit} hits, worst normalized eig={worst_bad:.4f}")

print()
print("=" * 72)
print(f"TOTAL FAILURES: {FAIL}")
print("=" * 72)
