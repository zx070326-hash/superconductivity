"""Independent check of universality-lift's congruence route for Anderson-Morel PSD.

Claim under test (their Lemma mu-3 v2):
    mu* := (I + mu L)^{-1} mu  ==  S (I + S L S)^{-1} S,   S := mu^{1/2}
and hence L symmetric PSD (NOT necessarily scalar, NOT necessarily commuting
with mu) already suffices for mu* symmetric PSD.

My own refinement under test:
    the sharp condition is  I + S L S  positive definite,
    which is strictly weaker than L PSD -- some ell_i < 0 is tolerated.

Nothing here is load-bearing; the algebra is analytic.  This is a consistency
check plus a count of how often the refinement actually buys something.
"""
import numpy as np

rng = np.random.default_rng(20260918)


def rand_psd(B, rng, eps=0.0):
    A = rng.normal(size=(B, B))
    return A @ A.T + eps * np.eye(B)


def sqrtm_sym(M):
    w, V = np.linalg.eigh(M)
    w = np.clip(w, 0.0, None)
    return (V * np.sqrt(w)) @ V.T


def mustar_direct(mu, L):
    return np.linalg.solve(np.eye(len(mu)) + mu @ L, mu)


def mustar_congruence(mu, L):
    S = sqrtm_sym(mu)
    return S @ np.linalg.solve(np.eye(len(mu)) + S @ L @ S, S)


def norm_min_eig(M):
    s = np.linalg.norm(M)
    return np.linalg.eigvalsh((M + M.T) / 2).min() / (s if s > 0 else 1.0)


# ---- 1. identity + PSD for L = diag(ell_i) > 0 (per-band cutoffs) ----
bad_id = bad_sym = bad_psd = 0
maxcomm = 0.0
for _ in range(4000):
    B = rng.integers(2, 6)
    mu = rand_psd(B, rng, eps=1e-3)
    ell = rng.uniform(0.05, 8.0, size=B)          # ln(W_i / omega_c^(i)) > 0
    L = np.diag(ell)
    A, C = mustar_direct(mu, L), mustar_congruence(mu, L)
    bad_id += np.max(np.abs(A - C)) > 1e-8 * max(1.0, np.linalg.norm(A))
    bad_sym += np.max(np.abs(C - C.T)) > 1e-9 * max(1.0, np.linalg.norm(C))
    bad_psd += norm_min_eig(C) < -1e-9
    maxcomm = max(maxcomm, np.linalg.norm(mu @ L - L @ mu))
print(f"L=diag(ell>0), 4000 draws: identity fails {bad_id}, "
      f"asym {bad_sym}, non-PSD {bad_psd}; max||[mu,L]|| = {maxcomm:.4f}")

# ---- 2. non-diagonal, non-commuting symmetric PSD L ----
bad_id = bad_psd = 0
for _ in range(4000):
    B = rng.integers(2, 6)
    mu, L = rand_psd(B, rng, eps=1e-3), rand_psd(B, rng)
    A, C = mustar_direct(mu, L), mustar_congruence(mu, L)
    bad_id += np.max(np.abs(A - C)) > 1e-8 * max(1.0, np.linalg.norm(A))
    bad_psd += norm_min_eig(C) < -1e-9
print(f"L sym PSD, non-commuting, 4000 draws: identity fails {bad_id}, non-PSD {bad_psd}")

# ---- 3. sharp criterion: I + S L S > 0, with L INDEFINITE ----
n_neg_ell = n_broken = n_saved = 0
n_crit_agree = 0
for _ in range(20000):
    B = rng.integers(2, 6)
    mu = rand_psd(B, rng, eps=1e-3)
    ell = rng.uniform(-3.0, 3.0, size=B)          # some omega_c^(i) > W_i
    L = np.diag(ell)
    if (ell >= 0).all():
        continue
    n_neg_ell += 1
    S = sqrtm_sym(mu)
    crit = np.linalg.eigvalsh(np.eye(B) + S @ L @ S).min() > 0
    try:
        C = mustar_direct(mu, L)
    except np.linalg.LinAlgError:
        n_broken += 1
        continue
    psd = norm_min_eig(C) > -1e-9
    n_crit_agree += (crit == psd)
    if psd:
        n_saved += 1
    else:
        n_broken += 1
print(f"L=diag with at least one ell_i<0: {n_neg_ell} draws -> "
      f"mu* still PSD in {n_saved}, broken in {n_broken}")
print(f"  sharp criterion (I+SLS>0) agrees with PSD verdict: "
      f"{n_crit_agree}/{n_neg_ell}")

# ---- 4. a fully indefinite L, to confirm the premise is non-empty ----
n_broken = 0
for _ in range(20000):
    B = rng.integers(2, 6)
    mu = rand_psd(B, rng, eps=1e-3)
    A = rng.normal(size=(B, B)) * 3.0
    L = (A + A.T) / 2
    if np.linalg.eigvalsh(L).min() >= 0:
        continue
    try:
        C = mustar_direct(mu, L)
    except np.linalg.LinAlgError:
        n_broken += 1
        continue
    if norm_min_eig(C) < -1e-9:
        n_broken += 1
print(f"general indefinite sym L: mu* non-PSD in {n_broken} / 20000")
