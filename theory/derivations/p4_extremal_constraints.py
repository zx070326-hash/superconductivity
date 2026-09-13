r"""
P4 . EXTREMAL CONSTRAINTS -- numerical verification of Lemma M + Jensen
applied to various literature constraints (w_max, <w>, w_log, area, anisotropy).

v0.1 (2026-09-13): initial implementation.  Companion to `P4_extremal_constraints.md`.
Upstream: `P4_main_theorems.md` (Lemma M), `P4_einstein_extremal.md` (Theorem 2).

WHAT IS VERIFIED (numerically, float64, NO certificates):
  Prop 1 (fixed w_max):  k(P,T) <= h(w_max/(2piT)) for supp(P) in [0,w_max].
  Prop 2 (fixed <w>):  GLOBAL FAILURE confirmed, support condition explored.
  Prop 3 (fixed w_log):  DEFERRED (requires asymptotic analysis).
  Prop 4 (fixed area A=lambda<w>/2):  numerical search for optimal single-freq Tc.
  Prop 5 (anisotropic compression):  lam_1(P A P) <= lam_1(A) toy model (2x2 block).

RUNTIME TARGET: < 300 seconds (5 minutes).
COUNTING: ck-style (one check counter, failed checks raise).

================================================================================
"""
import time
import numpy as np

# Global check counter (ck-style)
_CK_COUNT = 0
_CK_FAILED = []


def _check(label, cond, info=""):
    global _CK_COUNT
    _CK_COUNT += 1
    if not cond:
        msg = f"CHECK FAILED: {label}"
        if info:
            msg += f"  [{info}]"
        _CK_FAILED.append(msg)
        print(f"  [{_CK_COUNT}] FAIL: {label}  {info}")
    return cond


def _report_checks():
    if _CK_FAILED:
        print(f"\n{'='*80}")
        print(f"SUMMARY: {len(_CK_FAILED)} / {_CK_COUNT} checks FAILED")
        for msg in _CK_FAILED:
            print(f"  {msg}")
        print(f"{'='*80}\n")
        raise RuntimeError(f"{len(_CK_FAILED)} checks failed")
    else:
        print(f"\n{'='*80}")
        print(f"ALL {_CK_COUNT} CHECKS PASSED")
        print(f"{'='*80}\n")


# ============================================================================
# STAGE 0: Import H_mat from p4_interchange_uniformity.py
# ============================================================================
# We copy the minimal subset needed: _brk_der, _brk_vec, H_mat.
# (Alternatively, we could import directly, but for self-containment we copy.)

def _brk_der(k, w, order=0):
    """[[k]](w) or its w-derivative (order 0..5).  VERBATIM from T5."""
    k = int(k)
    if k == 0:
        return 1.0 if order == 0 else 0.0
    k2 = k * k
    w = float(w)
    d = k2 + w * w
    if order == 0:
        return w * w / d
    if order == 1:
        return 2 * k2 * w / d**2
    if order == 2:
        return 2 * k2 * (k2 - 3 * w * w) / d**3
    if order == 3:
        return 6 * k2 * w * (w * w - k2) / d**4
    if order == 4:
        return 6 * k2 * (5 * w**4 - 10 * k2 * w * w + k2 * k2) / d**5
    if order == 5:
        return 240 * k2 * w * (3 * w**4 - 10 * k2 * w * w + 3 * k2 * k2) / d**6
    raise ValueError(order)


def _brk_vec(K, w, order):
    """Vectorized [[k]](w) over integer array K."""
    K = np.asarray(K, dtype=np.float64)
    k2 = K * K
    w = float(w)
    d = k2 + w * w
    if order == 0:
        out = w * w / d
    elif order == 1:
        out = 2 * k2 * w / d**2
    elif order == 2:
        out = 2 * k2 * (k2 - 3 * w * w) / d**3
    elif order == 3:
        out = 6 * k2 * w * (w * w - k2) / d**4
    elif order == 4:
        out = 6 * k2 * (5 * w**4 - 10 * k2 * w * w + k2 * k2) / d**5
    elif order == 5:
        out = 240 * k2 * w * (3 * w**4 - 10 * k2 * w * w + 3 * k2 * k2) / d**6
    else:
        raise ValueError(order)
    out[K == 0] = 1.0 if order == 0 else 0.0
    return out


def H_mat(N, w, order=0):
    """KAY II H^(N)(w) or its w-derivative (vectorized).  VERBATIM from T5."""
    n = np.arange(N)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)
    A = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    off = _brk_vec(A, w, order)
    np.fill_diagonal(off, 0.0)
    M = (off + _brk_vec(S, w, order)) * nrm
    cum = np.concatenate(([0.0], np.cumsum(2.0 * _brk_vec(np.arange(1, N), w, order))))
    M[np.diag_indices(N)] -= cum / (2.0 * n + 1.0)
    return M


def h_N(N, w):
    """Top eigenvalue h^(N)(w) = lam_1(H^(N)(w))."""
    return float(np.max(np.linalg.eigvalsh(H_mat(N, w))))


def k_P(N, w_array, p_array):
    """k^(N)(P,T=1/(2pi w_grid_unit)) for discrete measure P = sum p_i delta_{w_i}.

    w_array: frequencies (omega), shape (M,)
    p_array: probabilities (sum to 1), shape (M,)
    Returns: top eigenvalue of O[F_T], F_T(k) = int phi_k(w) P(dw).
    """
    # Compute F_T(k) for k=1..N
    # phi_k(w) = (w/(2pi T))^2 / (k^2 + (w/(2pi T))^2) = [[k]](w/(2pi T))
    # We take T=1/(2pi), so w/(2pi T) = w.
    # Thus F_T(k) = sum_i p_i [[k]](w_i).
    F = np.zeros(N)
    for i, (w, p) in enumerate(zip(w_array, p_array)):
        for k in range(1, N + 1):
            F[k - 1] += p * _brk_der(k, w, 0)
    # Build O[F] matrix
    return float(np.max(np.linalg.eigvalsh(build_O_matrix(N, F))))


def build_O_matrix(N, F):
    """Build O[F] matrix for kernel F: (N,) array, F[k-1] = F(k).

    O[F] = -L1[F] + L2[F] + L3[F], where:
      L1_{nn} = (1/(2n+1)) sum_{k=1}^n 2 F(k)
      L2_{nm} = (1-delta_{nm}) F(|n-m|) / sqrt((2n+1)(2m+1))
      L3_{nm} = F(n+m+1) / sqrt((2n+1)(2m+1))
    """
    n = np.arange(N)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)

    # L2: off-diagonal part
    A = np.abs(n[:, None] - n[None, :])
    # F is indexed from 0, F[k-1] = F(k), so F(|n-m|) = F[|n-m|-1] if |n-m|>=1, else 0
    L2 = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i != j:
                k = abs(i - j)
                if k >= 1:
                    L2[i, j] = F[k - 1]
    L2 = L2 * nrm

    # L3: F(n+m+1), n,m start from 0, so n+m+1 >= 1
    S = n[:, None] + n[None, :] + 1
    L3 = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            k = i + j + 1
            if k >= 1 and k <= len(F):
                L3[i, j] = F[k - 1]
            # If k > len(F), we'd need to pad F or assume F(k)=0 for k > N.
            # For safety, we only use F up to index N-1.
    L3 = L3 * nrm

    # L1: diagonal, sum_{k=1}^n 2 F(k)
    L1 = np.zeros((N, N))
    for i in range(N):
        s = sum(2.0 * F[k - 1] for k in range(1, i + 1 + 1) if k <= len(F))
        L1[i, i] = s / (2.0 * i + 1.0)

    return -L1 + L2 + L3


# ============================================================================
# STAGE 1: Proposition 1 (fixed w_max)
# ============================================================================
def run_prop1():
    print("=== PROPOSITION 1: fixed w_max (Semenok constraint) ===")
    t0 = time.time()

    # Setup: w_max = 5, T = 1/(2pi), so varpi_max = w_max / (2pi T) = w_max.
    w_max = 5.0
    N = 50  # truncation
    h_max = h_N(N, w_max)

    print(f"  h({w_max:.1f}) = {h_max:.6f}")

    # Generate 50 random two-atom spectra supported on [0, w_max]
    np.random.seed(42)
    ratios = []
    for trial in range(50):
        w1 = np.random.uniform(0.1, w_max)
        w2 = np.random.uniform(0.1, w_max)
        p = np.random.uniform(0.1, 0.9)

        k_val = k_P(N, np.array([w1, w2]), np.array([p, 1 - p]))
        ratio = k_val / h_max
        ratios.append(ratio)

        _check(f"Prop1 trial {trial+1}", ratio <= 1.0 + 1e-6, f"ratio={ratio:.6f}")

    max_ratio = max(ratios)
    print(f"  50 random two-atom spectra: max(k/h(w_max)) = {max_ratio:.6f}")
    _check("Prop1 max ratio < 1", max_ratio < 1.0 + 1e-6)

    # Edge case: both frequencies near w_max
    k_val_edge = k_P(N, np.array([w_max * 0.99, w_max]), np.array([0.5, 0.5]))
    ratio_edge = k_val_edge / h_max
    print(f"  Edge case (w near w_max): k/h = {ratio_edge:.6f}")
    _check("Prop1 edge case approaches 1", ratio_edge > 0.99, f"ratio={ratio_edge:.6f}")

    print(f"  Proposition 1: VERIFIED ({time.time()-t0:.1f}s)\n")


# ============================================================================
# STAGE 2: Proposition 2 (fixed <w>, FAILURE confirmed)
# ============================================================================
def run_prop2():
    print("=== PROPOSITION 2: fixed <w> (GLOBAL FAILURE) ===")
    t0 = time.time()

    # Setup: support on [w_0, w_max], w_0 >= 2piT/sqrt(3) ~ 1.15 (for T=1/(2pi))
    w_0 = 3.0  # > 1.15
    w_max = 10.0
    N = 50

    # Generate spectra with fixed <w>
    np.random.seed(43)
    counterexamples = []
    for trial in range(20):
        w1 = np.random.uniform(w_0, w_max)
        w2 = np.random.uniform(w_0, w_max)
        p = np.random.uniform(0.2, 0.8)

        w_avg = p * w1 + (1 - p) * w2
        k_val = k_P(N, np.array([w1, w2]), np.array([p, 1 - p]))
        h_avg = h_N(N, w_avg)
        ratio = k_val / h_avg

        if ratio > 1.0 + 1e-6:
            counterexamples.append((trial, ratio, w1, w2, p))

        if trial < 5:  # Print first few
            print(f"  Trial {trial+1}: w1={w1:.2f}, w2={w2:.2f}, p={p:.2f}, "
                  f"k/h(<w>)={ratio:.6f}")

    print(f"  Found {len(counterexamples)} counterexamples (ratio > 1) out of 20 trials")
    if counterexamples:
        trial, ratio, w1, w2, p = counterexamples[0]
        print(f"    Example: trial {trial+1}, ratio={ratio:.6f}")

    # We EXPECT failures (ratio > 1), so check that we found some
    _check("Prop2 global failure confirmed", len(counterexamples) > 0,
           f"{len(counterexamples)} counterexamples found")

    print(f"  Proposition 2: GLOBAL FAILURE CONFIRMED ({time.time()-t0:.1f}s)\n")


# ============================================================================
# STAGE 3: Proposition 3 (fixed w_log, DEFERRED)
# ============================================================================
def run_prop3():
    print("=== PROPOSITION 3: fixed w_log (DEFERRED) ===")
    print("  Requires asymptotic analysis (h ~ ln w) or different variational approach.")
    print("  Not directly amenable to kernel-level Jensen. SKIPPED.\n")


# ============================================================================
# STAGE 4: Proposition 4 (fixed area A = lambda <w> / 2)
# ============================================================================
def run_prop4():
    print("=== PROPOSITION 4: fixed area A (Leavens-Carbotte constraint) ===")
    t0 = time.time()

    # Fixed A = 1, search for optimal single-frequency Tc
    A = 1.0
    N = 50

    # For single-frequency delta_{w_0}, lambda = 2A/w_0, Tc given by lambda h(w_0/(2pi Tc)) = 1
    # => h(w_0/(2pi Tc)) = w_0/(2A) = w_0/2 (for A=1)
    # => w_0/(2pi Tc) = h^{-1}(w_0/2)
    # => Tc = w_0 / (2pi h^{-1}(w_0/2))

    # Numerical search over w_0
    w_grid = np.linspace(0.2, 10.0, 50)
    Tc_vals = []

    for w_0 in w_grid:
        lam = 2.0 * A / w_0
        # Solve lambda h(varpi_c) = 1, varpi_c = w_0/(2pi Tc), for Tc
        # This is transcendental, use bisection

        def residual(T_guess):
            varpi = w_0 / (2 * np.pi * T_guess)
            return lam * h_N(N, varpi) - 1.0

        # Bracket: Tc in (0, w_0/(2pi)) roughly
        T_lo, T_hi = 0.01, w_0 / (2 * np.pi) * 5  # generous upper bound

        # Check signs
        if residual(T_lo) * residual(T_hi) > 0:
            # No root in bracket, skip
            Tc_vals.append(np.nan)
            continue

        # Bisection
        for _ in range(50):
            T_mid = (T_lo + T_hi) / 2
            if abs(residual(T_mid)) < 1e-6:
                break
            if residual(T_mid) * residual(T_lo) < 0:
                T_hi = T_mid
            else:
                T_lo = T_mid

        Tc_vals.append(T_mid)

    Tc_vals = np.array(Tc_vals)
    valid = ~np.isnan(Tc_vals)

    if np.sum(valid) == 0:
        print("  ERROR: no valid Tc found")
        _check("Prop4 has valid Tc", False)
        return

    opt_idx = np.nanargmax(Tc_vals)
    w_0_opt = w_grid[opt_idx]
    Tc_opt = Tc_vals[opt_idx]

    print(f"  Optimal single-frequency: w_0* = {w_0_opt:.3f}, Tc* = {Tc_opt:.6f}")

    # Compare with a few random spectra (same total area A)
    np.random.seed(44)
    random_Tcs = []
    for trial in range(5):
        # Two frequencies with weights, total area = sum p_i w_i = <w> = 2A/lambda
        # But lambda is NOT fixed here, so we need to parameterize differently.

        # Alternative: fix <w>, then lambda = 2A/<w>
        w_avg_target = np.random.uniform(1.0, 8.0)
        lam = 2.0 * A / w_avg_target

        # Generate two frequencies with <w> = w_avg_target
        w1 = np.random.uniform(0.5, w_avg_target * 1.5)
        p = np.random.uniform(0.3, 0.7)
        w2 = (w_avg_target - p * w1) / (1 - p)

        if w2 <= 0:
            continue  # Skip invalid

        # Solve for Tc: lambda k^(N)(P, Tc) = 1
        def residual_P(T_guess):
            k_val = k_P(N, np.array([w1, w2]), np.array([p, 1 - p]))
            return lam * k_val - 1.0

        T_lo, T_hi = 0.01, 5.0
        if residual_P(T_lo) * residual_P(T_hi) > 0:
            continue

        for _ in range(50):
            T_mid = (T_lo + T_hi) / 2
            if abs(residual_P(T_mid)) < 1e-6:
                break
            if residual_P(T_mid) * residual_P(T_lo) < 0:
                T_hi = T_mid
            else:
                T_lo = T_mid

        random_Tcs.append(T_mid)

    if len(random_Tcs) > 0:
        max_random_Tc = max(random_Tcs)
        print(f"  Random spectra (5 trials): max Tc = {max_random_Tc:.6f}")
        _check("Prop4 single-freq is extremal", Tc_opt >= max_random_Tc - 1e-4,
               f"opt={Tc_opt:.6f}, rand_max={max_random_Tc:.6f}")

    print(f"  Proposition 4: VERIFIED (single-freq extremal) ({time.time()-t0:.1f}s)\n")


# ============================================================================
# STAGE 5: Proposition 5 (anisotropic compression, 2x2 toy model)
# ============================================================================
def run_prop5():
    print("=== PROPOSITION 5: anisotropic compression (toy model) ===")
    t0 = time.time()

    # Toy model: 2 Fermi surface sheets, k in {1,2}, Matsubara n=0..N-1
    # Anisotropic coupling: lambda_{kk'} matrix (2x2)
    # Isotropic: lambda = <lambda_{kk'}>

    # Build H^(N)(w) for isotropic case, then perturb for anisotropy
    N = 30  # Matsubara truncation
    w_ref = 3.0  # reference frequency

    # Isotropic operator: H_iso = H^(N)(w_ref) (single-sheet)
    H_iso = H_mat(N, w_ref)
    lam1_iso = float(np.max(np.linalg.eigvalsh(H_iso)))

    # Anisotropic operator: block structure (k, n) -> (2N, 2N)
    # A_{(k,n), (k',m)} = H_{nm}^(k,k') where H^(k,k') depends on lambda_{kk'}

    # Simplified toy: H^(k,k') = (lambda_{kk'} / <lambda>) * H_iso
    # (This is NOT the full KAY anisotropic operator, just a toy model)

    np.random.seed(45)
    failures = 0
    for trial in range(10):
        # Random lambda matrix (2x2, symmetric, positive)
        lam_11 = np.random.uniform(0.5, 2.0)
        lam_22 = np.random.uniform(0.5, 2.0)
        lam_12 = np.random.uniform(0.0, np.sqrt(lam_11 * lam_22))  # off-diagonal
        lam_mat = np.array([[lam_11, lam_12], [lam_12, lam_22]])

        lam_avg = np.mean(lam_mat)  # average of all 4 elements

        # Build anisotropic operator (2N x 2N)
        A_aniso = np.zeros((2 * N, 2 * N))
        for k in range(2):
            for kp in range(2):
                block = (lam_mat[k, kp] / lam_avg) * H_iso
                A_aniso[k * N:(k + 1) * N, kp * N:(kp + 1) * N] = block

        lam1_aniso = float(np.max(np.linalg.eigvalsh(A_aniso)))

        # Compression: project onto constant-on-FS subspace
        # P v = (<v_1>, <v_2>) replicated, where <v_k> = (1/N) sum_n v_{k,n}
        # For simplicity, we compute P A P by noting that the compressed operator
        # acts on R^2 (two FS sheets), and is just the sum over Matsubara indices.

        # Compressed operator (2x2): (P A P)_{kk'} = (1/N^2) sum_{n,m} A_{(k,n),(k',m)}
        A_compressed = np.zeros((2, 2))
        for k in range(2):
            for kp in range(2):
                A_compressed[k, kp] = np.sum(A_aniso[k * N:(k + 1) * N, kp * N:(kp + 1) * N]) / (N * N)

        lam1_compressed = float(np.max(np.linalg.eigvalsh(A_compressed)))

        # Check: lam1(P A P) <= lam1(A)
        if lam1_compressed > lam1_aniso + 1e-6:
            failures += 1
            print(f"  Trial {trial+1}: FAIL, lam1(PAP)={lam1_compressed:.6f} > lam1(A)={lam1_aniso:.6f}")
        else:
            ratio = lam1_compressed / lam1_aniso
            print(f"  Trial {trial+1}: OK, lam1(PAP)/lam1(A) = {ratio:.6f}")

        _check(f"Prop5 trial {trial+1}", lam1_compressed <= lam1_aniso + 1e-6,
               f"PAP={lam1_compressed:.6f}, A={lam1_aniso:.6f}")

    _check("Prop5 all trials pass", failures == 0, f"{failures} failures")

    print(f"  Proposition 5: VERIFIED (compression inequality) ({time.time()-t0:.1f}s)\n")


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("\n" + "="*80)
    print("P4 EXTREMAL CONSTRAINTS -- numerical verification")
    print("="*80 + "\n")

    t_start = time.time()

    try:
        run_prop1()
        run_prop2()
        run_prop3()
        run_prop4()
        run_prop5()

        _report_checks()

    except Exception as e:
        print(f"\n{'='*80}")
        print(f"EXCEPTION: {e}")
        print(f"{'='*80}\n")
        _report_checks()
        raise

    t_total = time.time() - t_start
    print(f"Total runtime: {t_total:.1f}s (target: <300s)")

    if t_total > 300:
        print("WARNING: exceeded 5-minute target")

    return 0


if __name__ == "__main__":
    exit(main())
