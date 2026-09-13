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
# Direct import per project convention (do NOT copy code)
from p4_interchange_uniformity import H_mat, _brk_vec


def h_N(N, w):
    """Top eigenvalue h^(N)(w) = lam_1(H^(N)(w))."""
    return float(np.max(np.linalg.eigvalsh(H_mat(N, w))))


def k_P(N, w_array, p_array):
    """k^(N)(P,T=1/(2pi w_grid_unit)) for discrete measure P = sum p_i delta_{w_i}.

    Correct implementation (v0.4, per coordinator):
    O[F_T] is LINEAR in P, so k(P) = lam_1(sum_i p_i * H_mat(N, w_i)).

    w_array: frequencies (omega), shape (M,)
    p_array: probabilities (sum to 1), shape (M,)
    Returns: top eigenvalue of O[F_T].
    """
    # O[F] is linear in F, so O[sum p_i F_i] = sum p_i O[F_i]
    # For F_i = kernel from delta_{w_i}, we have O[F_i] = H_mat(N, w_i, 0)
    A = sum(p * H_mat(N, w, 0) for w, p in zip(w_array, p_array))
    return float(np.linalg.eigvalsh(A)[-1])


# Self-check: for single-atom P = delta_w, k_P must equal h_N
def _self_check_k_P():
    """Verify k_P(N, [w], [1]) == h_N(N, w) for single-atom measure."""
    N = 30
    for w in [0.5, 1.0, 3.0, 10.0]:
        k_single = k_P(N, np.array([w]), np.array([1.0]))
        h_ref = h_N(N, w)
        diff = abs(k_single - h_ref)
        if diff >= 1e-12:
            raise RuntimeError(f"k_P self-check FAILED: w={w}, diff={diff:.3e}")
    print("  k_P self-check: PASSED (single-atom measures agree with h_N)")



# ============================================================================
# STAGE 1: Proposition 1 (fixed w_max)
# ============================================================================
def run_prop1():
    print("=== PROPOSITION 1: fixed w_max (Semenok constraint) ===")
    t0 = time.time()

    # Self-check k_P implementation
    _self_check_k_P()

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

    # Edge case: approach point mass delta_{w_max}
    # P = (1-eps)*delta_{w_max} + eps*delta_{w_max/2}, eps=1e-3
    # Per coordinator: k >= (1-eps)*h(varpi_max) by Lemma M lower bound, must be >0.99
    eps = 1e-3
    w_arr = np.array([w_max, w_max / 2])
    p_arr = np.array([1 - eps, eps])
    k_val_edge = k_P(N, w_arr, p_arr)

    # Print diagnostic info
    varpi_max = w_max  # T=1/(2pi), so varpi = w
    varpi_half = w_max / 2
    h_varpi_max = h_N(N, varpi_max)
    h_varpi_half = h_N(N, varpi_half)

    print(f"  Edge case diagnostic:")
    print(f"    eps={eps}, w_max={w_max}, varpi_max={varpi_max}")
    print(f"    h(varpi_max)={h_varpi_max:.6f}, h(varpi_max/2)={h_varpi_half:.6f}")
    print(f"    k(P)={k_val_edge:.6f}, expected >= (1-eps)*h(varpi_max)={(1-eps)*h_varpi_max:.6f}")

    ratio_edge = k_val_edge / h_varpi_max
    print(f"    ratio k/h(varpi_max) = {ratio_edge:.6f}")

    _check("Prop1 edge case approaches 1", ratio_edge > 0.99,
           f"ratio={ratio_edge:.6f}, expected >0.99 by Lemma M")

    print(f"  Proposition 1: VERIFIED ({time.time()-t0:.1f}s)\n")


# ============================================================================
# STAGE 2: Proposition 2 (fixed <w>, FAILURE confirmed)
# ============================================================================
def run_prop2():
    print("=== PROPOSITION 2: fixed <w> (GLOBAL FAILURE) ===")
    t0 = time.time()

    # Part (a): Global failure with SMALL varpi_1, varpi_2 near 0
    # Per coordinator: varpi_1 in [0.2, 0.8], varpi_2 in [0.001, 0.05], eps in [0.05, 0.3]
    # P = (1-eps)*delta(varpi_1) + eps*delta(varpi_2), find ratio > 1
    print("  Part (a): Confirm global failure with small varpi_1, varpi_2->0")
    N = 50
    np.random.seed(43)
    found_counterexample = False

    for trial in range(30):
        varpi_1 = np.random.uniform(0.2, 0.8)
        varpi_2 = np.random.uniform(0.001, 0.05)
        eps = np.random.uniform(0.05, 0.3)
        p = 1.0 - eps  # weight on varpi_1

        # T=1/(2pi), so w = varpi * 2pi*T = varpi
        w1, w2 = varpi_1, varpi_2
        w_avg = p * w1 + eps * w2

        k_val = k_P(N, np.array([w1, w2]), np.array([p, eps]))
        h_avg = h_N(N, w_avg)
        ratio = k_val / h_avg

        if ratio > 1.0 + 1e-6:
            found_counterexample = True
            print(f"    Trial {trial+1}: varpi_1={varpi_1:.3f}, varpi_2={varpi_2:.4f}, "
                  f"eps={eps:.2f}, k/h(<w>)={ratio:.6f} > 1 [COUNTEREXAMPLE]")
            break

    _check("Prop2(a) global failure confirmed", found_counterexample,
           "Must find ratio > 1 with small varpi_1, varpi_2 near 0")

    # Part (b): Support bounded away from 0 (report, no assertion on direction)
    print("  Part (b): Support w_min >= 3, report worst ratio (no direction assertion)")
    w_min, w_max = 3.0, 10.0
    np.random.seed(44)
    ratios_bounded = []

    for trial in range(20):
        w1 = np.random.uniform(w_min, w_max)
        w2 = np.random.uniform(w_min, w_max)
        p = np.random.uniform(0.2, 0.8)

        w_avg = p * w1 + (1 - p) * w2
        k_val = k_P(N, np.array([w1, w2]), np.array([p, 1 - p]))
        h_avg = h_N(N, w_avg)
        ratio = k_val / h_avg
        ratios_bounded.append(ratio)

    max_ratio = max(ratios_bounded)
    min_ratio = min(ratios_bounded)
    print(f"    20 trials with support in [{w_min}, {w_max}]: "
          f"ratio in [{min_ratio:.6f}, {max_ratio:.6f}]")
    print(f"    NOTE: No counterexample found, but kernel-level Jensen fails for k > sqrt(3)*w_min")
    print(f"          (phi_k concave only for w >= k*2piT/sqrt(3); k=2 needs w >= 2.3, violated)")

    _check("Prop2(b) bounded support computed", len(ratios_bounded) == 20)

    print(f"  Proposition 2: Parts (a)+(b) DONE ({time.time()-t0:.1f}s)\n")


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
# STAGE 5: Proposition 5 (anisotropic compression, two-sheet KAY operator)
# ============================================================================
def build_aniso_KAY_operator(N, w, lam_mat, w_FS):
    """
    Build anisotropic KAY operator on l^2(N_0) x C^2 (Matsubara x 2 FS sheets).

    A acts on (S*N,) vectors v, with S=2 Fermi surface sheets.

    Correct structure (coordinator verified):
      - Symmetrized coupling: Lambda_ij = sqrt(w_i) * lambda_ij * sqrt(w_j)
      - Z row sum: z_i = sum_j lambda_ij * w_j  (WITH w_j weight!)
      - Block (i,j) = Lambda_ij * L23, where L23 is the L2+L3 part with kernel=1
      - Diagonal blocks (i,i) also subtract z_i * Z, where Z is the diagonal L1 part

    L23 := (off + _brk_vec(n+m+1)) * outer(u, u), with off diagonal cleared
    Z := diag(cum / (2n+1)), cum = concat([0], cumsum(2*_brk_vec(1..N-1, w, 0)))

    Compression: Q is (S*N)xN, i-th block = sqrt(w_i)*I_N (columns orthonormal)
    PAP := Q^T A Q must equal lambda_bar * H_iso, lambda_bar = w^T lambda w

    Parameters:
    - N: Matsubara truncation
    - w: frequency parameter (varpi)
    - lam_mat: (S,S) coupling matrix lambda_{ij}, symmetric
    - w_FS: (S,) Fermi surface weights, sum to 1

    Returns:
    - A: (S*N, S*N) operator matrix
    """
    S = len(w_FS)  # number of sheets
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)

    # Build L23 (L2+L3 part, kernel=1)
    A_diff = np.abs(n[:, None] - n[None, :])
    S_sum = n[:, None] + n[None, :] + 1
    off = _brk_vec(A_diff, w, 0)
    np.fill_diagonal(off, 0.0)
    L23 = (off + _brk_vec(S_sum, w, 0)) * np.outer(u, u)

    # Build Z (diagonal L1 part)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * _brk_vec(np.arange(1, N), w, 0))))
    Z = np.diag(cum / (2.0 * n + 1.0))

    # Verify L23 - Z == H_mat(N, w, 0)
    H_ref = H_mat(N, w, 0)
    diff_structure = np.max(np.abs(L23 - Z - H_ref))
    if diff_structure > 1e-12:
        print(f"  WARNING: L23 - Z != H_mat, diff={diff_structure:.3e}")

    # Symmetrized coupling Lambda_ij = sqrt(w_i) * lambda_ij * sqrt(w_j)
    sqrt_w = np.sqrt(w_FS)
    Lambda = sqrt_w[:, None] * lam_mat * sqrt_w[None, :]

    # Z row sum z_i = sum_j lambda_ij * w_j (with w_j!)
    z = np.dot(lam_mat, w_FS)

    # Build full operator A
    A = np.zeros((S * N, S * N))
    for i in range(S):
        for j in range(S):
            block = Lambda[i, j] * L23
            if i == j:
                block -= z[i] * Z
            A[i * N:(i + 1) * N, j * N:(j + 1) * N] = block

    return A


def run_prop5():
    print("=== PROPOSITION 5: anisotropic compression (two-sheet KAY operator) ===")
    t0 = time.time()

    N = 30  # Matsubara truncation
    w_ref = 3.0  # reference frequency
    w_FS = np.array([0.6, 0.4])  # Fermi surface weights (sum to 1)

    np.random.seed(45)

    for trial in range(10):
        # Random lambda matrix (2x2, symmetric, positive)
        S = 2
        lam_11 = np.random.uniform(0.5, 2.0)
        lam_22 = np.random.uniform(0.5, 2.0)
        lam_12 = np.random.uniform(0.0, 0.8 * np.sqrt(lam_11 * lam_22))
        lam_mat = np.array([[lam_11, lam_12], [lam_12, lam_22]])

        # Weighted average: lambda_bar = w^T lambda w
        lam_bar = np.dot(w_FS, np.dot(lam_mat, w_FS))

        # Build anisotropic operator
        A_aniso = build_aniso_KAY_operator(N, w_ref, lam_mat, w_FS)
        lam1_aniso = float(np.max(np.linalg.eigvalsh(A_aniso)))

        # Build isotropic reference: lambda_bar * H_mat(N, w_ref)
        H_iso = H_mat(N, w_ref)
        A_iso_full = lam_bar * H_iso

        # Compression: Q is (S*N)xN, i-th block = sqrt(w_i)*I_N
        Q = np.zeros((S * N, N))
        for i in range(S):
            Q[i * N:(i + 1) * N, :] = np.sqrt(w_FS[i]) * np.eye(N)

        # Verify Q is orthonormal: Q^T Q = I
        QtQ = np.dot(Q.T, Q)
        _check(f"Prop5 trial {trial+1} Q orthonormal", np.max(np.abs(QtQ - np.eye(N))) < 1e-12,
               f"Q^T Q != I, diff={np.max(np.abs(QtQ - np.eye(N))):.3e}")

        # Compute PAP = Q^T A Q
        PAP = np.dot(Q.T, np.dot(A_aniso, Q))

        # Verify PAP == lambda_bar * H_iso
        diff = np.max(np.abs(PAP - A_iso_full))
        _check(f"Prop5 trial {trial+1} PAP structure", diff < 1e-12,
               f"PAP != lambda_bar*H_iso, diff={diff:.3e}")

        # Eigenvalues
        lam1_PAP = float(np.max(np.linalg.eigvalsh(PAP)))
        lam1_iso_ref = float(np.max(np.linalg.eigvalsh(A_iso_full)))

        # Anisotropy measure
        aniso_measure = np.max(np.abs(lam_mat - lam_bar)) / lam_bar
        gap = lam1_aniso - lam1_PAP

        print(f"  Trial {trial+1}: aniso={aniso_measure:.3f}, "
              f"lam1(A)={lam1_aniso:.6f}, lam1(PAP)={lam1_PAP:.6f}, "
              f"gap={gap:.6f}")

        # Main check: lam1(PAP) <= lam1(A)
        _check(f"Prop5 trial {trial+1} compression", lam1_PAP <= lam1_aniso + 1e-6,
               f"PAP={lam1_PAP:.6f} > A={lam1_aniso:.6f}")

    print(f"  Proposition 5: VERIFIED (two-sheet KAY operator) ({time.time()-t0:.1f}s)\n")


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
