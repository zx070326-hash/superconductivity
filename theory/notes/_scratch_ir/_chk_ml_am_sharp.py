"""[numerical-evidence, consistency only]  Independent re-derivation of team-lead's
refinement of (mu3) + my own explicit sufficient condition + Metzler check on O_pair.

(a) Sharp criterion.  mu* = S (I+SLS)^{-1} S, S = mu^{1/2}.  <v,mu* v> = <Sv,(I+SLS)^{-1} Sv>.
    mu > 0  ==>  S invertible  ==>  Sv sweeps R^B  ==>  mu* >= 0  <==>  (I+SLS)^{-1} >= 0.
    (I+SLS) symmetric, invertible (det(I+SLS)=det(I+mu L) by Sylvester)  ==>  inverse PSD
    <==>  every eigenvalue of I+SLS is > 0  <==>  lambda_min(SLS) > -1.   (iff, mu > 0)
    mu singular: only "<==" (Sv sweeps range(S)); mu=0 gives mu*=0 PSD for any L.

(b) Explicit sufficient condition (mine, one line).  L = diag(ell), ell_min := min ell_i < 0:
    L >= ell_min I  ==>  SLS >= ell_min S^2 = ell_min mu >= ell_min ||mu|| I   (ell_min<0 flips)
    ==>  lambda_min(SLS) >= ell_min lambda_max(mu).
    Sufficient:  lambda_max(mu) * max_i ln(omega_c^(i)/W_i) < 1.
    Scalar case reproduces  1 + mu*ell > 0  <==>  omega_c < W e^{1/mu}  exactly.
    DIRECTION: (b) is a lower estimate of lambda_min(SLS); passing (b) implies passing (a);
    failing (b) implies nothing.

(c) O_pair Metzler?  In the Matsubara basis the same-n, i!=j entry is
    u_n^2 [ Lam(2n+1) - 2 sum_{k=1}^{n} Lam(k) ]_ij ,
    negative for every n>=1 whenever Lam(k)_ij is positive and non-increasing in k
    (Lam(2n+1) <= Lam(1) <= 2 sum_{k<=n} Lam(k)).  So O_pair is NOT Metzler for any
    B>=2 with a nonzero off-diagonal band element and N>=2 -- analytic, agrees with
    team-lead's staircase-basis conclusion.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.linalg import eigvalsh, solve, norm

rng = np.random.default_rng(918)
fails = 0


def sqrtm(M):
    w, V = np.linalg.eigh(M)
    return (V * np.sqrt(np.clip(w, 0, None))) @ V.T


def psd(M, tol=1e-9):
    s = norm(M)
    return eigvalsh((M + M.T) / 2).min() / (s if s > 0 else 1) > -tol


# ---------------- (a) iff, mu > 0 ----------------
n_a = agree = 0
n_singular_iff_fail = 0
for _ in range(20000):
    B = rng.integers(2, 6)
    A = rng.normal(size=(B, B)); mu = A @ A.T + 1e-2 * np.eye(B)   # mu > 0
    ell = rng.uniform(-3, 3, size=B); L = np.diag(ell)
    S = sqrtm(mu)
    if abs(np.linalg.det(np.eye(B) + mu @ L)) < 1e-8:
        continue
    crit = eigvalsh(S @ L @ S).min() > -1
    ms = solve(np.eye(B) + mu @ L, mu)
    n_a += 1; agree += (crit == psd(ms))
print(f"(a) mu>0, L=diag indefinite: criterion lam_min(SLS)>-1 vs PSD agree {agree}/{n_a}")
fails += (agree != n_a)

# singular mu: only <== ; look for PSD-but-criterion-false (allowed) and criterion-true-but-nonPSD (forbidden)
forb = allowed = n_s = 0
for _ in range(20000):
    B = rng.integers(2, 6)
    A = rng.normal(size=(B, B - 1)); mu = A @ A.T        # rank B-1
    ell = rng.uniform(-3, 3, size=B); L = np.diag(ell)
    S = sqrtm(mu)
    if abs(np.linalg.det(np.eye(B) + mu @ L)) < 1e-8:
        continue
    crit = eigvalsh(S @ L @ S).min() > -1
    p = psd(solve(np.eye(B) + mu @ L, mu))
    n_s += 1
    forb += (crit and not p); allowed += (p and not crit)
print(f"    mu singular ({n_s}): crit&!PSD (forbidden) {forb}; PSD&!crit (allowed, only '<==') {allowed}")
fails += (forb != 0)

# ---------------- (b) explicit sufficient condition ----------------
n_b = n_b_pass = n_b_pass_psd = n_a_pass = 0
for _ in range(20000):
    B = rng.integers(2, 6)
    A = rng.normal(size=(B, B)); mu = A @ A.T + 1e-2 * np.eye(B)
    ell = rng.uniform(-3, 3, size=B); L = np.diag(ell)
    if (ell >= 0).all():
        continue
    n_b += 1
    S = sqrtm(mu)
    lmin = eigvalsh(S @ L @ S).min()
    lower = ell.min() * eigvalsh(mu).max()
    if lower > lmin + 1e-12:
        fails += 1; print("  (b) lower estimate ABOVE true lam_min -- BUG")
    b_pass = lower > -1
    a_pass = lmin > -1
    n_b_pass += b_pass; n_a_pass += a_pass
    if b_pass and not a_pass:
        fails += 1; print("  (b) passes but (a) fails -- BUG")
print(f"(b) at least one ell_i<0 ({n_b}): sharp (a) PSD {n_a_pass}, explicit (b) PSD {n_b_pass}, "
      f"(b)=>(a) violations 0 expected")

# scalar sanity: mu=0.3 -> omega_c/W < e^{1/0.3}
for mu0 in (0.1, 0.2, 0.3, 0.5):
    print(f"    scalar mu={mu0}: mu* stays PSD up to omega_c/W = e^(1/mu) = {np.exp(1/mu0):.1f}x")

# ---------------- (c) O_pair Metzler ----------------
import importlib.util, sys
spec = importlib.util.spec_from_file_location(
    "B", os.path.join(os.path.dirname(os.path.abspath(__file__)), "_chk_ml_build.py"))
Bm = importlib.util.module_from_spec(spec); spec.loader.exec_module(Bm)
N = 60
n = np.arange(N); u = 1 / np.sqrt(2 * n + 1.0)
for w0 in (0.5, 1.0, 3.0):
    lam = lambda k: w0 ** 2 / (k ** 2 + w0 ** 2)
    # off-diag band entry lam_12 = c>0, same n: u_n^2 [lam(2n+1) - 2 sum_{k=1}^n lam(k)]
    ks = np.arange(0, 2 * N + 2, dtype=float); lv = lam(ks)
    same_n = u ** 2 * (lv[2 * n + 1] - 2 * np.concatenate(([0.], np.cumsum(lv[1:N]))))
    neg = int((same_n[1:] < 0).sum())
    # cross-check against the scalar builder's diagonal (O_pair block (1,2) = c * O[lam])
    Osc = Bm.scalar_O(lv, N)
    d = None if Osc is None else np.max(np.abs(np.diag(Osc) - same_n))
    print(f"(c) w={w0}: same-n off-band entry of O_pair negative for {neg}/{N-1} n>=1; "
          f"n=0 entry {same_n[0]:+.4f}; builder-diag maxdiff {d}")
    fails += (neg != N - 1)

print(f"\nTOTAL FAILURES: {fails}")
