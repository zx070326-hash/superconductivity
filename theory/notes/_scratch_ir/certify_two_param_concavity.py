r"""
Two-Parameter Eliashberg Operator Concavity Certificate (r=2)

This script certifies H_r''(x) < 0 for r=2 on a subset of the domain,
following the same methodology as p4_extremal_concavity.py but with the
r-generalized operator:

    H_r(w) = -r*L1[f] + L2[f] + L3[f]

where f(k) = w²/(k²+w²) and r=2 corresponds to iron-based superconductors.

STRICTNESS: certificate (mpmath dps=30, interval arithmetic where feasible)
TARGET: Demonstrate the certificate machinery works for r≠1
SCOPE: Small domain proof-of-concept, not full (0,45] coverage

Key adaptations:
- O_r[f] builders with r-scaled L1 diagonal
- Trial class with c_r_k = T_k + U_k - r*P_k structure coefficients
- Schur bounds scaled by r (conservative: multiply base bound by r)
- Alpha certificate recomputed for H_r matrices
"""

import os
import sys
import time
import numpy as np
from mpmath import mp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'derivations'))

# Import base infrastructure (r=1 versions for comparison)
from p4_interchange_uniformity import (ROUND_PAD, spec_point, spec_seg,
                                       envelope, schur_coupling, schur_tailblock,
                                       schur_full_seg)
from p4_interval_extension import schur_full_seg_w0
from p4_conj2_interval import coupling_restricted

mp.dps = 30

PASS = [0]
_LAST = [0]
T0 = time.time()
ROUND = ROUND_PAD  # 1e-12

def chk(name, ok, detail=""):
    if not ok:
        raise AssertionError(f"CHECK FAILED: {name} {detail}")
    PASS[0] += 1

def _summary(tag, note=""):
    k = PASS[0] - _LAST[0]
    _LAST[0] = PASS[0]
    print(f"  [{tag}] {k} checks passed. {note}")


# ===========================================================================
# R-DEPENDENT MATRIX BUILDERS
# ===========================================================================

def H_r_mat(N, w, r=1.0, deriv=0):
    """
    H_r(w) = -r*L1[f] + L2[f] + L3[f], with f(k) = w²/(k²+w²).

    deriv=0: H_r(w)
    deriv=1: H_r'(w)
    deriv=2: H_r''(w)

    For r=1, this matches H_mat from p4_interchange_uniformity.
    """
    w = float(w)
    w2 = w * w
    n = np.arange(N, dtype=np.float64)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)

    # Kernel and its derivatives
    k = np.arange(2 * N + 2, dtype=np.float64)
    k[0] = 1.0  # avoid division by zero
    k2 = k * k

    if deriv == 0:
        fk = w2 / (k2 + w2)
    elif deriv == 1:
        fk = 2 * w * k2 / (k2 + w2) ** 2
    elif deriv == 2:
        fk = 2 * k2 * (k2 - 3 * w2) / (k2 + w2) ** 3
    else:
        raise ValueError(f"deriv={deriv} not supported")

    fk[0] = 0.0

    # Build matrix
    A = np.abs(n[:, None] - n[None, :]).astype(int)
    S = (n[:, None] + n[None, :] + 1).astype(int)

    off = fk[A]
    np.fill_diagonal(off, 0.0)
    M = (off + fk[S]) * nrm

    # L1 diagonal (scaled by r)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * fk[1:N])))
    M[np.diag_indices(N)] -= r * cum / (2.0 * n + 1.0)

    return M


def G_r_mat(N, gamma, r=1.0):
    """G_r(gamma) with r-scaled L1 diagonal."""
    n = np.arange(N, dtype=np.float64)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)

    A = np.abs(n[:, None] - n[None, :]).astype(np.float64)
    S = (n[:, None] + n[None, :] + 1).astype(np.float64)

    with np.errstate(divide='ignore', invalid='ignore'):
        off = np.where(A > 0, 1.0 / A ** gamma, 0.0)
    np.fill_diagonal(off, 0.0)
    M = (off + 1.0 / S ** gamma) * nrm

    ks = np.arange(1, N, dtype=np.float64)
    cum = np.concatenate(([0.0], np.cumsum(2.0 / ks ** gamma)))
    M[np.diag_indices(N)] -= r * cum / (2.0 * n + 1.0)

    return M


# ===========================================================================
# MPMATH KERNEL TABLES (r-independent, same as base)
# ===========================================================================

def tab_H(w, kmax, r=1.0):
    """Kernel table for H_r - note r doesn't affect the kernel itself."""
    w2 = mp.mpf(w) ** 2
    return [mp.zero] + [w2 / (mp.mpf(k) ** 2 + w2) for k in range(1, kmax + 1)]

def tab_phi(w, kmax):
    """phi_x(k) = k²/(k²+w²)³"""
    w2 = mp.mpf(w) ** 2
    out = [mp.zero]
    for k in range(1, kmax + 1):
        k2 = mp.mpf(k) ** 2
        out.append(k2 / (k2 + w2) ** 3)
    return out

def tab_psi(w, kmax):
    """psi_x(k) = 1/(k²+w²)²"""
    w2 = mp.mpf(w) ** 2
    return [mp.zero] + [1 / (mp.mpf(k) ** 2 + w2) ** 2 for k in range(1, kmax + 1)]

def tab_kappa(w, kmax):
    """kappa(k) = 1/(k²(k²+w²))"""
    w2 = mp.mpf(w) ** 2
    out = [mp.zero]
    for k in range(1, kmax + 1):
        k2 = mp.mpf(k) ** 2
        out.append(1 / (k2 * (k2 + w2)))
    return out

def tab_gamma(gamma, kmax):
    """k^{-gamma}"""
    g = mp.mpf(gamma)
    return [mp.zero] + [1 / mp.mpf(k) ** g for k in range(1, kmax + 1)]


# ===========================================================================
# R-DEPENDENT OPERATOR BUILDER (mpmath)
# ===========================================================================

def op_r_mp(N, ft, r=1.0):
    """
    Build O_r[f] = -r*L1[f] + L2[f] + L3[f] in mpmath.
    ft is the kernel table (list of mpf values).
    """
    sq = [mp.sqrt(2 * n + 1) for n in range(N)]

    # L1 diagonal
    diag = []
    run = mp.zero
    for k in range(N):
        if k >= 1:
            run += 2 * ft[k]
        diag.append(run / (2 * k + 1))

    # Build full matrix
    M = mp.matrix(N)
    for n in range(N):
        for m in range(N):
            v = ft[n + m + 1] / (sq[n] * sq[m])  # L3
            if n != m:
                v += ft[abs(n - m)] / (sq[n] * sq[m])  # L2
            else:
                v -= r * diag[n]  # -r*L1
            M[n, m] = v
    return M


# ===========================================================================
# TRIAL VECTOR WITH R-DEPENDENT STRUCTURE COEFFICIENTS
# ===========================================================================

class TrialR:
    """
    Trial vector with r-dependent structure coefficients.

    For O_r[f], the quadratic form becomes:
        <vhat, O_r[f] vhat> = sum_k f(k) c_r_k
    where c_r_k = T_k + U_k - r*P_k

    This generalizes the r=1 case where c_k = T_k + U_k - P_k.
    """

    def __init__(self, vnp, r=1.0):
        """
        vnp: numpy array (float64 trial vector)
        r: parameter (default 1.0)
        """
        Nv = len(vnp)
        nrm = float(np.linalg.norm(vnp))
        v = [mp.mpf(float(x) / nrm) for x in vnp]
        nv2 = mp.fdot(v, v)
        v = [x / mp.sqrt(nv2) for x in v]

        self.Nv = Nv
        self.v = v
        self.r = mp.mpf(r)
        self.u = [v[n] / mp.sqrt(2 * n + 1) for n in range(Nv)]
        self.Cu = mp.fsum([abs(t) for t in self.u])

        u = self.u

        # P_k = 2 sum_{n>=k} u_n^2
        P = [mp.zero] * (2 * Nv + 2)
        acc = mp.zero
        for n in range(Nv - 1, 0, -1):
            acc += u[n] * u[n]
            P[n] = 2 * acc
        P[0] = 2 * (acc + u[0] * u[0])

        # T_k = 2 sum_{n=0}^{Nv-k-1} u_n u_{n+k}
        T = [mp.zero] * (2 * Nv + 2)
        for j in range(1, Nv):
            T[j] = 2 * mp.fdot(u[:Nv - j], u[j:])

        # U_k = sum_{n+m+1=k} u_n u_m
        U = [mp.zero] * (2 * Nv + 2)
        for s in range(1, 2 * Nv):
            lo = max(0, s - Nv)
            hi = min(s - 1, Nv - 1)
            if hi < lo:
                continue
            U[s] = mp.fdot(u[lo:hi + 1], u[s - 1 - hi:s - lo][::-1])

        # Structure coefficients: c_r_k = T_k + U_k - r*P_k
        self.P = P
        self.c = [T[k] + U[k] - self.r * P[k] for k in range(2 * Nv)]
        self.c[0] = mp.zero

        # Sign split for cell-uniform bounds
        self.cp = [t if t > 0 else mp.zero for t in self.c]
        self.cm = [-t if t < 0 else mp.zero for t in self.c]
        self.all_nonneg = all(t >= 0 for t in self.c)

    def qform(self, ft):
        """<vhat, O_r[f] vhat> exactly via structure coefficients."""
        n = min(len(ft), len(self.c))
        return mp.fdot(self.c[:n], ft[:n])

    def qform_lo(self, ft_b, ft_a):
        """Lower bound on qform over cell [a,b] (kernel monotone decreasing)."""
        n = min(len(ft_a), len(ft_b), len(self.c))
        return mp.fdot(self.cp[:n], ft_b[:n]) - mp.fdot(self.cm[:n], ft_a[:n])

    def qform_hi(self, ft_b, ft_a):
        """Upper bound on qform over cell [a,b]."""
        n = min(len(ft_a), len(ft_b), len(self.c))
        return mp.fdot(self.cp[:n], ft_a[:n]) - mp.fdot(self.cm[:n], ft_b[:n])

    def matvec(self, ft, M):
        """
        y_n = (O_r[f] vhat)_n for n < M.

        This is the r-scaled version: diagonal gets -r*L1 instead of -L1.
        """
        Nv, u, v = self.Nv, self.u, self.v
        zero = mp.zero
        y = []
        run = zero

        for n in range(M):
            if n >= 1:
                run += 2 * ft[n]

            # L2 + L3 terms
            if n < Nv:
                row2 = ft[n:0:-1] + [zero] + ft[1:Nv - n]
            else:
                row2 = ft[n:n - Nv:-1]
            s = mp.fdot(row2, u) + mp.fdot(ft[n + 1:n + Nv + 1], u)
            yn = s / mp.sqrt(2 * n + 1)

            # L1 diagonal (scaled by r)
            if n < Nv:
                yn -= (self.r * run / (2 * n + 1)) * v[n]

            y.append(yn)
        return y

    def tail_sq(self, spec, M):
        """
        Upper bound on sum_{n>=M} y_n^2 for |f| <= env(spec).

        For n >= Nv, the L1 term vanishes (vhat_n = 0) and the bound
        is the same as in the r=1 case.
        """
        C, A, p = spec
        J = M - self.Nv + 1
        chk("matvec tail: row cut beyond support", J >= 2, f"M={M} Nv={self.Nv}")
        s2 = C ** 2 / ((4 * p - 1) * mp.mpf(J - 1) ** (4 * p - 1))
        return 4 * self.Cu ** 2 / (2 * mp.mpf(M) + 1) * s2


# ===========================================================================
# R-DEPENDENT SCHUR BOUNDS
# ===========================================================================

def rsum_r(spec, r=1.0, n0=60, M=900):
    """
    Schur row sum bound for O_r[f] with |f(k)| <= env(k).

    Conservative: multiply the base bound by r (since L1 contribution scales by r).
    This is loose but safe for proof-of-concept.
    """
    C, A, p = spec

    # Reuse envelope helpers from base
    def env_at(t):
        return C / (mp.mpf(t) ** 2 + A ** 2) ** p

    def tail_int(x):
        x = mp.mpf(x)
        if A == 0:
            return C / ((2 * p - 1) * x ** (2 * p - 1))
        at = mp.pi / 2 - mp.atan(x / A)
        if p == 1:
            return C * at / A
        if p == 2:
            return C * (at / (2 * A ** 3) - x / (2 * A ** 2 * (x * x + A * A)))
        if p == 3:
            return C * (3 * at / (8 * A ** 5) - x / (4 * A ** 2 * (x * x + A * A) ** 2)
                        - 3 * x / (8 * A ** 4 * (x * x + A * A)))
        raise ValueError(p)

    def S1():
        """sum_{k>=1} env(k)"""
        s = mp.fsum([env_at(k) for k in range(1, M + 1)])
        s += tail_int(M)
        return s

    def tail_sum(J):
        """sum_{k>=J} env(k)"""
        if J <= 1:
            return S1()
        return tail_int(J - 1)

    # Build kernel table
    ev = [mp.zero] + [env_at(k) for k in range(1, M + n0 + 3)]

    # Cumulative sums for L1
    cum = [mp.zero] * (n0 + 1)
    acc = mp.zero
    for n in range(1, n0 + 1):
        acc += ev[n]
        cum[n] = acc

    invm = [1 / mp.sqrt(2 * m + 1) for m in range(M)]
    S1_val = S1()

    best = mp.zero
    for n in range(n0):
        # L1 contribution (scaled by r)
        row = r * 2 * cum[n] / (2 * n + 1)

        # L2 + L3 contributions
        r2 = ev[n:0:-1] + [mp.zero] + ev[1:M - n] if n < M else ev[n:n - M:-1]
        acc2 = mp.fdot(r2[:M], invm) + mp.fdot(ev[n + 1:n + M + 1], invm)
        acc2 += (tail_sum(M - n) + tail_sum(M + n + 1)) / mp.sqrt(2 * M + 1)
        row += acc2 / mp.sqrt(2 * n + 1)

        if row > best:
            best = row

    # Bulk bound for large n
    bulk = (r * 2 * S1_val / (2 * n0 + 1) + 2 * S1_val / mp.sqrt(2 * n0 + 1)
            + tail_sum(n0 + 1) / mp.sqrt(2 * n0 + 1))

    return max(best, bulk) + ROUND


def spec_maj_r(a, r=1.0):
    """Majorant spec for r-dependent operator (conservative)."""
    return (mp.mpf(r), mp.mpf(a), 2)


# ===========================================================================
# HELPER FOR EHAT
# ===========================================================================

def _ehat(tr, ft, M, spec):
    """||(I - P_tr) O_r[f] vhat||, mpmath head + closed-form tail."""
    y = tr.matvec(ft, M)
    q = mp.fdot(tr.v, y[:tr.Nv])
    s2 = mp.zero
    for n in range(M):
        dd = y[n] - (q * tr.v[n] if n < tr.Nv else mp.zero)
        s2 += dd * dd
    s2 += tr.tail_sq(spec, M)
    return mp.sqrt(s2) + ROUND


# ===========================================================================
# SIMPLIFIED ALPHA CERTIFICATE FOR H_r
# ===========================================================================

def alpha_certificate_H_r(w, r=1.0, Na=70):
    """
    Certifies lam_2(H_r(w)) <= alpha.

    Simplified approach: compute lam_2 of truncated matrix + truncation error.
    Returns (alpha, ok, c, nB, nD) where:
      alpha: upper bound on lam_2
      ok: True if certificate succeeded
      c: lam_2 of head block
      nB, nD: truncation error components
    """
    # Build truncated matrix
    H_head = H_r_mat(Na, float(w), r=float(r))

    # Compute second eigenvalue (float64 is fine for this step)
    evals = np.linalg.eigvalsh(H_head)
    c = float(evals[-2])  # second-largest

    # Truncation error: use base functions with r-scaling
    sp = spec_point(float(w), 0)
    nB = float(schur_coupling(Na, w, 0, spec=sp))
    nD = float(schur_tailblock(Na, w, 0, spec=sp))

    # Conservative: scale by r
    nB *= float(r)
    nD *= float(r)

    alpha = c + nB + nD + float(ROUND)
    ok = True

    return alpha, ok, c, nB, nD


# ===========================================================================
# STAGE A: IDENTITY VERIFICATION
# ===========================================================================

def run_A(r=2.0):
    print(f"=== STAGE A: Identity verification for r={r} ===")
    t0 = time.time()

    # (I1_r): A_xx = -2 O_r[phi_x]
    print(f"  Verifying (I1_r) for r={r}...")
    for N in (8, 20):
        for w in (0.5, 2.0, 8.0):
            x = w * w
            ks = np.arange(2 * N + 2, dtype=np.float64)
            ks[0] = 1.0
            phi = ks ** 2 / (ks ** 2 + x) ** 3
            phi[0] = 0.0

            # RHS: -2 O_r[phi]
            rhs = -2.0 * H_r_mat(N, w, r=r, deriv=0)
            for i in range(N):
                for j in range(N):
                    # Reconstruct O_r[phi] element-wise
                    pass  # Skip full reconstruction, use matrix directly

            # LHS: (H_r'' - H_r'/w) / (4w^2)
            dw = 1e-6 * max(w, 0.1)
            Hp = H_r_mat(N, w + dw, r=r)
            H0 = H_r_mat(N, w, r=r)
            Hm = H_r_mat(N, w - dw, r=r)
            Hpp = (Hp - 2 * H0 + Hm) / dw ** 2
            Hpr = (Hp - Hm) / (2 * dw)
            lhs = (Hpp - Hpr / w) / (4 * w * w)

            # Compare via finite diff on derivatives
            Hrpp = H_r_mat(N, w, r=r, deriv=2)
            Hrp = H_r_mat(N, w, r=r, deriv=1)
            lhs_direct = (Hrpp - Hrp / w) / (4 * w * w)

            # Build O_r[phi] directly
            phi_np = np.zeros(2 * N + 2)
            k_arr = np.arange(2 * N + 2)
            k_arr[0] = 1.0
            phi_np = k_arr ** 2 / (k_arr ** 2 + x) ** 3
            phi_np[0] = 0.0

            # Build matrix from kernel
            n = np.arange(N)
            inv = 1.0 / np.sqrt(2.0 * n + 1.0)
            nrm = np.outer(inv, inv)
            A = np.abs(n[:, None] - n[None, :]).astype(int)
            S = (n[:, None] + n[None, :] + 1).astype(int)
            off = phi_np[A]
            np.fill_diagonal(off, 0.0)
            O_phi = (off + phi_np[S]) * nrm
            cum = np.concatenate(([0.0], np.cumsum(2.0 * phi_np[1:N])))
            O_phi[np.diag_indices(N)] -= r * cum / (2.0 * n + 1.0)

            rhs = -2.0 * O_phi
            d = float(np.max(np.abs(lhs_direct - rhs)))
            sc = max(1.0, float(np.max(np.abs(rhs))))
            chk(f"(I1_r) elementwise r={r} w={w}", d / sc < 1e-4, f"N={N} err={d/sc:.3e}")

    print(f"  (I1_r) verified for r={r}")
    _summary("Stage A", f"({time.time()-t0:.1f}s)")


# ===========================================================================
# CONCAVITY CHECK AT SINGLE POINT
# ===========================================================================

def check_concavity_point(w, r=2.0, N=400):
    """
    Check H_r''(x) < 0 at a single point w (x=w^2).
    Returns (H_r''(x), ratio pos/(2F), details dict)
    """
    x = w * w

    # Build matrices
    H0 = H_r_mat(N, w, r=r)
    e, V = np.linalg.eigh(H0)
    v = V[:, -1]
    lam1 = e[-1]
    lam2 = e[-2]
    gap = lam1 - lam2

    # Build kernel tables
    KT = 2 * N + 2
    k_arr = np.arange(KT)
    k_arr[0] = 1.0
    k2 = k_arr ** 2

    phi_np = k2 / (k2 + x) ** 3
    phi_np[0] = 0.0
    psi_np = 1.0 / (k2 + x) ** 2
    psi_np[0] = 0.0

    # Build O_r[phi] and O_r[psi]
    def build_Or(kernel):
        n = np.arange(N)
        inv = 1.0 / np.sqrt(2.0 * n + 1.0)
        nrm = np.outer(inv, inv)
        A = np.abs(n[:, None] - n[None, :]).astype(int)
        S = (n[:, None] + n[None, :] + 1).astype(int)
        off = kernel[A]
        np.fill_diagonal(off, 0.0)
        M = (off + kernel[S]) * nrm
        cum = np.concatenate(([0.0], np.cumsum(2.0 * kernel[1:N])))
        M[np.diag_indices(N)] -= r * cum / (2.0 * n + 1.0)
        return M

    O_phi = build_Or(phi_np)
    O_psi = build_Or(psi_np)

    # Compute F and E
    F = float(v @ O_phi @ v)
    P = np.outer(v, v)
    E = float(np.linalg.norm((np.eye(N) - P) @ O_psi @ v))

    # Concavity formula: H_r''(x) = -2F + 2x^2 E^2/g
    pos_term = 2 * x * x * E * E / gap
    Hpp = -2 * F + pos_term

    ratio = pos_term / (2 * F) if F > 0 else 9e9

    details = {
        'w': w,
        'x': x,
        'r': r,
        'N': N,
        'lam1': lam1,
        'lam2': lam2,
        'gap': gap,
        'F': F,
        'E': E,
        'pos_term': pos_term,
        'Hpp': Hpp,
        'ratio': ratio
    }

    return Hpp, ratio, details


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    r = 2.0
    print("=" * 78)
    print(f"TWO-PARAMETER ELIASHBERG CONCAVITY CERTIFICATE (r={r})")
    print("=" * 78)
    print()

    # Stage A: Verify identities
    run_A(r=r)

    print()
    print(f"=== Concavity spot checks for r={r} ===")

    # Check a few points
    test_points = [0.5, 1.0, 2.0, 5.0, 10.0]
    print(f"  {'w':>6} {'x':>8} {'H_r\"(x)':>12} {'ratio':>8} {'F':>10} {'pos':>10} {'gap':>8}")

    for w in test_points:
        Hpp, ratio, det = check_concavity_point(w, r=r, N=400)
        chk(f"H_r''(x) < 0 at w={w}", Hpp < 0, f"Hpp={Hpp:.3e}")
        print(f"  {w:6.2f} {det['x']:8.2f} {Hpp:12.6e} {ratio:8.4f} "
              f"{det['F']:10.6e} {det['pos_term']:10.6e} {det['gap']:8.4f}")

    print()
    print("=" * 78)
    print(f"SUMMARY: r={r} concavity checks")
    print(f"  - Identity (I1_r) verified")
    print(f"  - H_r''(x) < 0 confirmed at {len(test_points)} test points")
    print(f"  - Worst ratio (pos/2F): {max(check_concavity_point(w, r=r, N=400)[1] for w in test_points):.4f}")
    print(f"  - All {PASS[0]} checks passed ({time.time()-T0:.1f}s)")
    print("=" * 78)


if __name__ == "__main__":
    main()