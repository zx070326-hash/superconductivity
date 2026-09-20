#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_zdom_budget_2026-09-21.py  --  2026-09-21, ceiling-fe agent, FROM ZERO.

TASK D (budget) of the "multiband ceiling closure on iron-based s+-" job.

QUESTION: (Z-dom) is the last surviving physical premise of the multiband
ceiling.  If it FAILS on real s+- materials, what is left?  Per repo discipline
(CLAUDE.md "predicate-first budgeting"), degrade (Z-dom) in the direction that
makes the criterion HARDER and find where the verdict flips.

DEGRADATION PARAMETER:
    (Z-dom, kappa_Z)   |Lam^Delta_ij(k)| <= kappa_Z * Lam^Z_ij(k)  elementwise,
                       kappa_Z >= 1.   kappa_Z = 1 is the repo's (Z-dom).

Everything here is built from the DEFINITIONS, not from any repo script.

Scalar convention (repo P4_main_theorems.md Sec.0):
    u_n = (2n+1)^{-1/2}
    O[f]_{nm} = u_n u_m [ f(|n-m|) 1_{n!=m} + f(n+m+1) ]
                - delta_nm * 2 u_n^2 sum_{k=1}^{n} f(k)
    g(2) := lam_1( O[1/k^2] ) = 1.318140212132...

Multiband linearised Eliashberg on the Matsubara axis (omega_n = piT(2n+1)):
    Lam^X_ij(k) = int 2w a2F^X_ij(w) / (w^2 + (2 pi T k)^2) dw       (X = Z, Delta)
    Z_i(n)       = 1 + u_n^2 ( R^Z_i(0) + 2 sum_{k=1}^{n} R^Z_i(k) ),
                   R^Z_i(k) = sum_j Lam^Z_ij(k)
    Z_i(n)Delta_i(n) = sum_{j,m} [Lam^D_ij(|n-m|) + Lam^D_ij(n+m+1)] Delta_j(m)/(2m+1)
  psi_i(n) = u_n Delta_i(n)  =>  eigenvalue-1 criterion  O psi = psi  with
    O_{(n,i),(m,j)} = u_n u_m [Lam^D_ij(|n-m|) 1_{n!=m} + Lam^D_ij(n+m+1)]
                      + delta_nm [ u_n^2 Lam^D_ij(0)
                                   - delta_ij u_n^2 (R^Z_i(0)+2 sum_{k=1}^n R^Z_i(k)) ]
  (symmetrised: detailed balance N_i lam_ij = N_j lam_ji, w_i = sqrt(N_i),
   Lam~_ij = sqrt(N_i N_j) V_ij symmetric; below we work directly with
   symmetric Lam~ so the conjugation is already done.)

DIRECTION AUDIT (single-sided-bound rule, CLAUDE.md):
   lam_1(P_N O P_N) <= lam_1(O)  (Rayleigh quotient restricted to a subspace)
   => a truncated numerator is a LOWER estimate of lam_1
   => ratio = lam_1(P_N O P_N) * kappa / ( g(2) rho(M2'') ) > 1 is a TRUE
      violation (sound FAIL);  ratio <= 1 is NOT by itself a proof of PASS.
   All FAIL verdicts below are therefore sound; all "0 violations" rows are
   reported as evidence only.

OMP_NUM_THREADS=1 is set below (repo mechanical convention).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

# ----------------------------------------------------------------- scalar side
def u_vec(N):
    n = np.arange(N)
    return 1.0 / np.sqrt(2.0 * n + 1.0)


def scalar_O(f, N):
    """f: callable on integers k>=0.  Returns the NxN matrix O[f]."""
    u = u_vec(N)
    n = np.arange(N)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    fD = np.vectorize(f)(D).astype(float)
    fS = np.vectorize(f)(S).astype(float)
    np.fill_diagonal(fD, 0.0)                     # the 1_{n!=m}
    M = np.outer(u, u) * (fD + fS)
    tail = np.array([2.0 * u[i] ** 2 * sum(f(k) for k in range(1, i + 1))
                     for i in range(N)])
    M -= np.diag(tail)
    return M


def gate_g2(N=600):
    O = scalar_O(lambda k: 0.0 if k == 0 else 1.0 / k ** 2, N)
    return np.linalg.eigvalsh(0.5 * (O + O.T))[-1]


# -------------------------------------------------------------- multiband side
def multiband_O(LamD, LamZ, N):
    """LamD, LamZ : callables k(int) -> (B,B) symmetric numpy arrays.
       Builds O in the basis ordered (n, i), n=0..N-1, i=0..B-1."""
    B = LamD(0).shape[0]
    u = u_vec(N)
    dim = N * B
    O = np.zeros((dim, dim))
    # cache
    maxk = 2 * N
    LD = [LamD(k) for k in range(maxk + 2)]
    LZ = [LamZ(k) for k in range(maxk + 2)]
    for n in range(N):
        for m in range(N):
            k1 = abs(n - m)
            k2 = n + m + 1
            blk = LD[k2].copy()
            if n != m:
                blk = blk + LD[k1]
            blk = blk * (u[n] * u[m])
            if n == m:
                blk = blk + u[n] ** 2 * LD[0]          # the k=0 pairing term
                RZ = np.array([LZ[0][i].sum() + 2.0 * sum(LZ[k][i].sum()
                               for k in range(1, n + 1)) for i in range(B)])
                blk = blk - np.diag(u[n] ** 2 * RZ)
            O[n * B:(n + 1) * B, m * B:(m + 1) * B] = blk
    return 0.5 * (O + O.T)


def O_pair(LamD, N):
    """The relaxed operator: scalar O acting elementwise on the matrix kernel."""
    B = LamD(0).shape[0]
    dim = N * B
    P = np.zeros((dim, dim))
    for i in range(B):
        for j in range(B):
            Oij = scalar_O(lambda k: LamD(k)[i, j], N)
            P[i::B, j::B] = Oij
    # re-order: above used stride B with offset i => index = n*B + i . consistent
    return 0.5 * (P + P.T)


def D_extra(LamZ, N, B):
    """Diagonal operator  D = sum_k Z_k (x) diag(R^Z(k)),  D_{(n,i)} = Z_i(n)-1."""
    u = u_vec(N)
    LZ = [LamZ(k) for k in range(N + 2)]
    d = np.zeros(N * B)
    for n in range(N):
        for i in range(B):
            d[n * B + i] = u[n] ** 2 * (LZ[0][i].sum()
                                        + 2.0 * sum(LZ[k][i].sum()
                                                    for k in range(1, n + 1)))
    return d


# ------------------------------------------------------------- Einstein kernels
def einstein(lam, wbar):
    """Lam(k) = lam * wbar^2/(wbar^2 + k^2)  (wbar = Omega/(2 pi T), dimensionless).
       M2/kappa = lam * wbar^2 ,  kappa = (2 pi T)^2."""
    def f(k):
        return lam * wbar ** 2 / (wbar ** 2 + k ** 2)
    return f


def mat_kernel(lams, wbars):
    """lams, wbars: (B,B) arrays.  Returns callable k -> (B,B)."""
    lams = np.asarray(lams, float)
    wbars = np.asarray(wbars, float)

    def f(k):
        return lams * wbars ** 2 / (wbars ** 2 + k ** 2)
    return f


C_GF = None  # filled in main


def m2_matrix(lams, wbars):
    return np.asarray(lams, float) * np.asarray(wbars, float) ** 2


# ------------------------------------------------------------------------ main
def main():
    np.set_printoptions(precision=6, suppress=False)
    print("=" * 78)
    print("GATE 0  -- scalar reproduction")
    print("=" * 78)
    for N in (200, 400, 800):
        print("   g(2) at N=%4d : %.12f" % (N, gate_g2(N)))
    g2 = gate_g2(1200)
    print("   g(2) at N=1200: %.12f    (repo rational-exact CERT-D bracket"
          " [1.318140211762671, 1.318140212373478])" % g2)
    g2_repo = 1.318140212132
    print("   |mine - repo|  : %.3e" % abs(g2 - g2_repo))
    cGF = (np.pi ** 2 / 6.0) / g2
    Cinf = np.sqrt(g2) / (2 * np.pi)
    print("   c_GF = (pi^2/6)/g(2) = %.10f ;  C_inf = sqrt(g2)/2pi = %.10f"
          % (cGF, Cinf))

    print()
    print("=" * 78)
    print("GATE 1  -- B=1 multiband operator == scalar operator, elementwise")
    print("=" * 78)
    for (lam, wb) in ((0.7, 0.5), (2.5, 1.7), (1.0, 0.13)):
        N = 60
        f = einstein(lam, wb)
        Om = multiband_O(lambda k: np.array([[f(k)]]),
                         lambda k: np.array([[f(k)]]), N)
        Os = scalar_O(f, N)
        print("   lam=%.2f wbar=%.2f : maxdiff = %.3e" % (lam, wb,
              np.abs(Om - Os).max()))

    print()
    print("=" * 78)
    print("GATE 2  -- (ML1) identity:  O = O_pair[Lam^D] - sum_k Z_k (x) L(k),")
    print("           L(k) = diag(R^Z(k)) - Lam~^D(k).   Check on SIGNED Lam^D.")
    print("=" * 78)
    rng = np.random.default_rng(20260921)
    for trial in range(4):
        B = rng.integers(2, 5)
        N = 30
        lamD = rng.normal(size=(B, B)); lamD = lamD + lamD.T
        lamZ = np.abs(rng.normal(size=(B, B))) + 0.2; lamZ = lamZ + lamZ.T
        wb = np.abs(rng.normal(size=(B, B))) + 0.3; wb = wb + wb.T
        KD = mat_kernel(lamD, wb)
        KZ = mat_kernel(lamZ, wb)
        O = multiband_O(KD, KZ, N)
        # rebuild via (ML1)
        P = np.zeros_like(O)
        for i in range(B):
            for j in range(B):
                P[i::B, j::B] = scalar_O(lambda k: KD(k)[i, j], N)
        P = 0.5 * (P + P.T)
        u = u_vec(N)
        Corr = np.zeros_like(O)
        for n in range(N):
            # sum_{k>=0} Z_k(n) L(k) ; Z_0 = u_n^2, Z_k = 2u_n^2 1_{k<=n}
            acc = np.zeros((B, B))
            for k in range(0, n + 1):
                Lk = np.diag(KZ(k).sum(axis=1)) - KD(k)
                acc += (u[n] ** 2 if k == 0 else 2 * u[n] ** 2) * Lk
            Corr[n * B:(n + 1) * B, n * B:(n + 1) * B] = acc
        print("   B=%d  maxdiff |O - (O_pair - Corr)| = %.3e"
              % (B, np.abs(O - (P - Corr)).max()))

    print()
    print("=" * 78)
    print("TASK D -- BUDGET.  Analytic claim derived by hand, verified here:")
    print("""
   Relax (Z-dom) to  |Lam^D_ij(k)| <= kappa_Z Lam^Z_ij(k)  elementwise.
   Then  diag(kappa_Z R^Z(k)) - Lam~^D(k) >= 0  (same AM-GM proof), so

     O = O_pair[Lam^D] - sum_k Z_k(x)(diag(R^Z) - Lam^D)
       = O_pair[Lam^D] - sum_k Z_k(x)(diag(kappa_Z R^Z) - Lam^D)
                       + (kappa_Z - 1) sum_k Z_k (x) diag(R^Z)
       <= O_pair[Lam^D] + (kappa_Z - 1) * D ,     D = diag(Z_i(n) - 1) >= 0.

   D's largest entry is at n=0:  D_{(0,i)} = R^Z_i(0) = lam^Z_i := sum_j lam^Z_ij
   and for n>=1, D_{(n,i)} = (lam^Z_i + 2 sum_{k=1}^n R^Z_i(k))/(2n+1)
                           <= (lam^Z_i + 2n lam^Z_i)/(2n+1) = lam^Z_i
   because R^Z_i(k) is non-increasing in k.  So  max D = lam^Z_max.   [analytic]

   At T = Tc, lam_1(O) = 1, hence

     1 <= g(2) rho(M2''_Z)/(2 pi Tc)^2 + (kappa_Z - 1) lam^Z_max
   => Tc <= C_inf sqrt( rho(M2''_Z) / (1 - (kappa_Z-1) lam^Z_max) )
                                      valid iff (kappa_Z-1) lam^Z_max < 1
   => NO BOUND AT ALL once  kappa_Z >= 1 + 1/lam^Z_max.          [analytic]
""")
    print("=" * 78)
    print("   D0 check: max_n,i D_{(n,i)} == lam^Z_max ?")
    for trial in range(4):
        B = int(rng.integers(2, 5)); N = 80
        lamZ = np.abs(rng.normal(size=(B, B))) + 0.05; lamZ = lamZ + lamZ.T
        wb = np.abs(rng.normal(size=(B, B))) + 0.2; wb = wb + wb.T
        KZ = mat_kernel(lamZ, wb)
        d = D_extra(KZ, N, B)
        print("      B=%d  max D = %.10f   lam^Z_max = %.10f   diff %.2e"
              % (B, d.max(), lamZ.sum(axis=1).max(),
                 d.max() - lamZ.sum(axis=1).max()))

    print()
    print("   Budget table: 1 + 1/lam^Z_max  for real materials (lam^Z from lit.)")
    for name, lam in (("Nb", 1.0), ("MgB2 (sigma)", 1.0), ("LiFeAs (tot)", 1.5),
                      ("Ba-122 s+- (Ummarino lam31=3.866,nu1=1)", 4.34),
                      ("H3S", 2.3), ("SmFeAsO (Ummarino lam31=14.52)", 15.4)):
        print("      %-42s lam^Z_max = %5.2f  ->  kappa_Z budget = %.4f"
              % (name, lam, 1.0 + 1.0 / lam))

    # ---- numerical: does the ACTUAL ceiling get violated as kappa_Z grows? ----
    print()
    print("=" * 78)
    print("   NUMERICAL: max observed ratio  lam_1(P_N O P_N)*kappa /"
          " (g2 * rho(M2''_Z))")
    print("   over random (Z-dom,kappa_Z)-admissible configs.  ratio>1 = TRUE")
    print("   violation (truncation under-estimates the numerator).")
    print("=" * 78)
    N = 90
    for kz in (1.0, 1.1, 1.25, 1.5, 2.0, 3.0):
        worst = 0.0
        worst_cfg = None
        nviol = 0
        ntr = 240
        for _ in range(ntr):
            B = int(rng.integers(2, 4))
            # build a physical (K1)-style split then break it by kappa_Z
            lph = np.abs(rng.normal(size=(B, B))) * rng.uniform(0.2, 3.0)
            lsf = np.abs(rng.normal(size=(B, B))) * rng.uniform(0.2, 3.0)
            lph = lph + lph.T; lsf = lsf + lsf.T
            wph = np.abs(rng.normal(size=(B, B))) + 0.1; wph = wph + wph.T
            wsf = np.abs(rng.normal(size=(B, B))) + 0.1; wsf = wsf + wsf.T
            KZ = lambda k, a=lph, b=lsf, wa=wph, wb_=wsf: (
                a * wa ** 2 / (wa ** 2 + k ** 2) + b * wb_ ** 2 / (wb_ ** 2 + k ** 2))
            # Delta channel: phonon - sf, then AMPLIFIED by kappa_Z (breaks Z-dom)
            KD = lambda k, a=lph, b=lsf, wa=wph, wb_=wsf, s=kz: s * (
                a * wa ** 2 / (wa ** 2 + k ** 2) - b * wb_ ** 2 / (wb_ ** 2 + k ** 2))
            O = multiband_O(KD, KZ, N)
            lam1 = np.linalg.eigvalsh(O)[-1]
            if lam1 <= 0:
                continue
            M2ph = lph * wph ** 2
            M2sf = lsf * wsf ** 2
            M2pp = np.diag(np.diag(M2ph + cGF * M2sf)) + cGF * (
                (M2ph + M2sf) - np.diag(np.diag(M2ph + M2sf)))
            rho = np.abs(np.linalg.eigvals(M2pp)).max()
            ratio = lam1 / (g2 * rho)      # kappa already absorbed (wbar scaled)
            if ratio > 1:
                nviol += 1
            if ratio > worst:
                worst = ratio
                worst_cfg = (B, lph.max(), lsf.max())
        print("   kappa_Z = %.2f : %3d/%d violations, worst ratio = %.4f  %s"
              % (kz, nviol, ntr, worst, worst_cfg))

    print()
    print("   NOTE: this random family samples the compensation ratio")
    print("   r = M2^sf/M2^ph from a continuous distribution, so it gives")
    print("   ZERO mass to the degenerate face r=1 (continuous-type degeneracy,")
    print("   CLAUDE.md).  Directed scan at r -> 1 follows.")
    print()
    for kz in (1.0, 1.15, 1.3):
        for r in (1.0, 0.99, 0.9):
            worst = 0.0; nviol = 0; ntr = 60
            for _ in range(ntr):
                B = 2
                lph = np.abs(rng.normal(size=(B, B))) + 0.1; lph = lph + lph.T
                wph = np.abs(rng.normal(size=(B, B))) + 0.15; wph = wph + wph.T
                wsf = np.abs(rng.normal(size=(B, B))) + 0.15; wsf = wsf + wsf.T
                M2ph = lph * wph ** 2
                M2sf = r * M2ph
                lsf = M2sf / wsf ** 2
                KZ = lambda k, a=lph, b=lsf, wa=wph, wb_=wsf: (
                    a * wa**2/(wa**2+k**2) + b * wb_**2/(wb_**2+k**2))
                KD = lambda k, a=lph, b=lsf, wa=wph, wb_=wsf, s=kz: s * (
                    a * wa**2/(wa**2+k**2) - b * wb_**2/(wb_**2+k**2))
                O = multiband_O(KD, KZ, N)
                lam1 = np.linalg.eigvalsh(O)[-1]
                if lam1 <= 0:
                    continue
                M2pp = np.diag(np.diag(M2ph + cGF*M2sf)) + cGF*(
                    (M2ph+M2sf) - np.diag(np.diag(M2ph+M2sf)))
                rho = np.abs(np.linalg.eigvals(M2pp)).max()
                ratio = lam1/(g2*rho)
                nviol += (ratio > 1)
                worst = max(worst, ratio)
            print("   kappa_Z=%.2f  r=%.2f : %2d/%d violations, worst %.4f"
                  % (kz, r, nviol, ntr, worst))


if __name__ == "__main__":
    # c_GF needs g2 first; recompute cheaply
    _g2 = gate_g2(1200)
    cGF = (np.pi ** 2 / 6.0) / _g2
    g2 = _g2
    main()
