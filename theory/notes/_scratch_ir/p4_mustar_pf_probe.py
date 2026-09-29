#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
P4-D mu*>0 uniqueness/extremality -- Perron-Frobenius / irreducibility probe.

Purpose (strategy-doc support, NOT a theorem):
  (A) Verify the load-bearing structural identities against source
      (O[1]=2uu^T-I ; K_c=2 u_c u_c^T ; O[d]_00 = d(1)).
  (B) Demonstrate the r=3 cross-window DOUBLE ROOT (ties) at lambda=1, mu~=0.10,
      wbar in [2.000, 2.012], reproducing P4D_coulomb_uniqueness.md sec.6.
  (C) Test PERRON-FROBENIUS / IRREDUCIBILITY of A = K_ph  vs  M = A - mu~ K_c:
      - A has all strictly positive off-diagonal entries (=> irreducible, PF:
        lambda1 simple, psi1 > 0);
      - M is still a DENSE (graph-irreducible) matrix but has MIXED-SIGN
        off-diagonal entries at physical mu~ -- so PF/positivity does NOT apply.
      This is the crux of why the cone/single-index route dies and why PF is
      NOT a global tie-breaker.
  (D) Count sign changes of the top eigenvector of M (falsify "oscillating cone":
      expect exactly 2, N-independent).
  (E) Extremality Q(P) <= Q(E): show the FREE bound  lambda1(M_P) <= h(wrms)
      [mu*=0 Einstein] always holds, but the genuine Q(P) <= Q(E) [both mu*>0]
      can be VIOLATED above the near-crossing threshold.

All float64 (numpy) for scans; mpmath spot-checks for the key numbers.
Run from theory/notes/_scratch_ir/  (gitignored).  Uses only numpy + mpmath.
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 40

# ----------------------------------------------------------------------
# Core operator O[f] on l^2(N_0), truncated to N coordinates.
#   u_n = (2n+1)^{-1/2}
#   [[k]](w) = w^2/(k^2+w^2)
#   O[f] = -L1 + L2 + L3
#   L1_nn = (2n+1)^{-1} * sum_{k=1}^{n} 2 f(k)      (diagonal)
#   L2_nm = (1-delta) f(|n-m|) u_n u_m               (off-diagonal)
#   L3_nm = f(n+m+1) u_n u_m
# ----------------------------------------------------------------------

def u_vec(N):
    n = np.arange(N)
    return 1.0 / np.sqrt(2 * n + 1.0)

def O_of_f(fvals, N):
    """fvals[k] = f(k) for k=0..(2N).  Returns NxN symmetric matrix O[f]."""
    u = u_vec(N)
    n = np.arange(N)
    # L2 + L3 off/all:
    I, J = np.meshgrid(n, n, indexing='ij')
    absd = np.abs(I - J)
    s = I + J + 1
    L2 = np.where(I != J, fvals[absd], 0.0) * np.outer(u, u)
    L3 = fvals[s] * np.outer(u, u)
    M = L2 + L3
    # L1 diagonal:
    # cumulative 2 f(k) for k=1..n
    csum = np.zeros(N)
    running = 0.0
    for k in range(1, N):
        running += 2.0 * fvals[k]
        csum[k] = running
    diagL1 = csum / (2 * n + 1.0)
    M[n, n] -= diagL1
    return M

def kernel_einstein(w, K):
    k = np.arange(K + 1)
    return (w * w) / (k * k + w * w)

def H_op(w, N):
    """Einstein operator Hbar(w) = O[[[.]](w)]."""
    fv = kernel_einstein(w, 2 * N + 2)
    return O_of_f(fv, N)

def K_c(N, Nc):
    """Coulomb operator 2 |u_c><u_c|, u_c[n]=(2n+1)^{-1/2} 1_{n<Nc}."""
    u = u_vec(N).copy()
    u[Nc:] = 0.0
    return 2.0 * np.outer(u, u)

def top_eig(M):
    w, V = np.linalg.eigh(M)
    return w[-1], V[:, -1], w

def sign_changes(v):
    """number of sign changes in v, ignoring (near-)zeros."""
    s = np.sign(v[np.abs(v) > 1e-10 * np.max(np.abs(v))])
    return int(np.sum(s[1:] != s[:-1]))

# ----------------------------------------------------------------------
print("=" * 72)
print("(A) STRUCTURAL IDENTITIES")
print("=" * 72)
N = 60
u = u_vec(N)
O1 = O_of_f(np.ones(2 * N + 3), N)
target = 2.0 * np.outer(u, u) - np.eye(N)
print("  max|O[1] - (2uu^T - I)| =", np.max(np.abs(O1 - target)))

Nc = 10
Kc = K_c(N, Nc)
P = np.zeros((N, N));
for i in range(Nc): P[i, i] = 1.0
Kc_ident = P @ (O1 + np.eye(N)) @ P
print("  max|K_c - P(O[1]+I)P| =", np.max(np.abs(Kc - Kc_ident)))
uc = u.copy(); uc[Nc:] = 0.0
print("  max|K_c - 2 u_c u_c^T| =", np.max(np.abs(Kc - 2 * np.outer(uc, uc))))

# O[d]_00 = d(1) for arbitrary kernel d
rng = np.random.default_rng(0)
dv = np.abs(rng.standard_normal(2 * N + 3))
Od = O_of_f(dv, N)
print("  O[d]_00 =", Od[0, 0], " vs d(1) =", dv[1], " diff =", abs(Od[0, 0] - dv[1]))

# ----------------------------------------------------------------------
print("=" * 72)
print("(C) IRREDUCIBILITY / PERRON-FROBENIUS of A vs M")
print("=" * 72)
for wbar in [0.5, 1.0, 2.0, 5.0]:
    N = 200
    A = H_op(wbar, N)
    off = A - np.diag(np.diag(A))
    min_off = np.min(off[~np.eye(N, dtype=bool)])
    # top eigenpair of A
    lamA, psiA, specA = top_eig(A)
    gapA = specA[-1] - specA[-2]
    # psi1 sign
    psi_min = np.min(psiA * np.sign(psiA[0]))
    print(f"  A=Hbar({wbar}): min off-diag = {min_off:.3e} (>0 => irreducible, PF)"
          f"  lambda1 gap = {gapA:.4e}  min(psi1 up to sign) = {psi_min:.3e}")
    # Now M = A - mu~ K_c at physical mu~
    for mut in [0.10, 0.26]:
        Nc = max(1, int(np.floor(3 * wbar)))  # r=3 physical Coulomb
        M = A - mut * K_c(N, Nc)
        offM = M - np.diag(np.diag(M))
        offvals = offM[~np.eye(N, dtype=bool)]
        n_neg = int(np.sum(offvals < 0))
        min_offM = np.min(offvals)
        # graph-irreducibility of nonzero pattern:
        nz_dense = np.all(np.abs(offvals) > 0)
        print(f"      M=A-{mut}*K_c (r=3,Nc={Nc}): #neg off-diag = {n_neg}, "
              f"min off = {min_offM:.3e}, all-nonzero(dense/irred graph)={nz_dense}")

# ----------------------------------------------------------------------
print("=" * 72)
print("(D) SIGN CHANGES of top eigenvector of M (falsify oscillating cone)")
print("=" * 72)
for (wbar, Nc) in [(1.0, 10), (2.0, 20), (5.0, 40)]:
    for Ntry in [200, 400, 800]:
        A = H_op(wbar, Ntry)
        M = A - 0.13 * K_c(Ntry, Nc)
        lam, v, _ = top_eig(M)
        v = v * np.sign(v[0])
        sc = sign_changes(v)
        n_neg = int(np.sum(v < -1e-9))
        print(f"  (wbar={wbar},Nc={Nc},N={Ntry},mu~=0.13): sign-changes={sc}  #(v<0)={n_neg}")

# ----------------------------------------------------------------------
print("=" * 72)
print("(B) r=3 CROSS-WINDOW DOUBLE ROOT  (lambda=1, mu~=0.10)")
print("=" * 72)
# reproduce P4D_coulomb_uniqueness.md sec.6: wbar in [2.000, 2.012], r=3
N = 400
mut = 0.10
r = 3
for wbar in [1.990, 1.995, 2.000, 2.006, 2.012, 2.020]:
    Nc = int(np.floor(r * wbar))
    A = H_op(wbar, N)
    M = A - mut * K_c(N, Nc)
    lam, v, _ = top_eig(M)
    # lambda*lambda1(M) = 1 boundary; here we track lambda1(M) directly at lambda=1
    print(f"  wbar={wbar:.3f}  Nc=floor({r}*{wbar:.3f})={Nc}  lambda1(M)={lam:.6f}")
print("  -> at lambda=1 the equation lambda*lambda1(M)=1 is crossed twice across")
print("     the Nc: 5->6 jump (down-jump), giving a genuine double root in T.")

# ----------------------------------------------------------------------
print("=" * 72)
print("(E) EXTREMALITY  Q(P) <= Q(E)  vs free bound lambda1(M_P) <= h(wrms)")
print("=" * 72)
# Build a 2-atom P with same <w^2> as an Einstein spectrum at wrms.
# Compare Q(P)=lambda1(K(P,T)-mu~K_c) with Q(E)=lambda1(Hbar(wrms)-mu~K_c)
# and with h(wrms)=lambda1(Hbar(wrms)) [mu*=0].
def K_of_P(atoms, weights, N):
    """K(P,T)=O[F_P], F_P(k)=sum_i w_i [[k]](wbar_i)."""
    Kmax = 2 * N + 2
    fv = np.zeros(Kmax + 1)
    for (wi, wt) in zip(atoms, weights):
        fv += wt * kernel_einstein(wi, Kmax)
    return O_of_f(fv, N)

N = 200
r = 3
# choose a spread 2-atom P; wrms^2 = sum w_i wbar_i^2
cases = [
    ([0.3, 4.0], [0.8, 0.2]),   # from H6 probe (dangerous), high spread
    ([0.5, 3.0], [0.5, 0.5]),
    ([1.0, 2.0], [0.5, 0.5]),
]
for atoms, weights in cases:
    wrms = np.sqrt(sum(w * a * a for a, w in zip(atoms, weights)))
    Nc = max(1, int(np.floor(r * wrms)))
    KP = K_of_P(atoms, weights, N)
    HE = H_op(wrms, N)
    hE = top_eig(HE)[0]  # mu*=0 Einstein value
    print(f"  P={atoms} w={weights}  wrms={wrms:.4f}  Nc={Nc}")
    for mut in [0.05, 0.13, 0.26, 0.33]:
        Kc_ = K_c(N, Nc)
        QP = top_eig(KP - mut * Kc_)[0]
        QE = top_eig(HE - mut * Kc_)[0]
        free_ok = QP <= hE + 1e-12       # free bound: pencil <= mu*=0 Einstein
        extremal_ok = QP <= QE + 1e-12   # genuine extremality both mu*>0
        flag = "" if extremal_ok else "  <<< Q(P) > Q(E) VIOLATION"
        print(f"      mu~={mut:.2f}: Q(P)={QP:.6f} Q(E)={QE:.6f} h(wrms)={hE:.6f}"
              f"  free<=hE:{free_ok} extremal:{extremal_ok}{flag}")

print("=" * 72)
print("done")
