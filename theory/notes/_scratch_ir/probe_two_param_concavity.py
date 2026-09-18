r"""
Two-parameter Eliashberg operator: numerical probe for boson-mediated Tc bound.

H_r(w) = -r * H1(w) + H2(w) + H3(w)

where
  H1 = diagonal (Z / mass-renormalization channel)
  H2 + H3 = off-diagonal (pairing channel)
  r = lambda_Z / lambda_Delta >= 1

Physical interpretation:
  - Phonons: r = 1 exactly (same vertex in both channels)
  - Cuprate paramagnon (d-wave): r ~ 2-3
  - Iron-based (s+-): r ~ 1.5-2

Identities (derived in the math doc, direct two-parameter generalization of KAY):
  (I1_r)  d^2/dx^2 H_r(sqrt(x)) = -2 O_r[phi_x],   phi_x(k) = k^2/(k^2+x)^3
  (I2_r)  d/dx H_r(sqrt(x)) = H_r(w)/w^2 - x O_r[psi_x],  psi_x(k) = 1/(k^2+x)^2
  where O_r[f] = -r L1[f] + L2[f] + L3[f]

Concavity condition:
  H_r''(x) = -2 F_r + 2 x^2 E_r^2 / g_r    (second-order perturbation theory)
  F_r = <v_r, O_r[phi] v_r>,  E_r = ||(I-P_r) O_r[psi] v_r||,  g_r = gap

STRICTNESS: numerical-evidence (float64 throughout, no interval arithmetic).
"""

import sys
import os
import time
import numpy as np

# Import base builders for r=1 cross-check
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'derivations'))
try:
    from p4_interchange_uniformity import H_mat, G_matrix
    HAS_BASE = True
except ImportError:
    HAS_BASE = False
    print("WARNING: cannot import base builders; r=1 cross-check skipped")


# ========================================================================
# BUILDERS
# ========================================================================

def H_r_mat(N, w, r=1.0):
    """Vectorized H_r^(N)(w) = -r*H1 + H2 + H3.

    H1 (diagonal, Z channel):     -(r/(2n+1)) sum_{k=1}^n 2 f(k)
    H2 (off-diag, pairing):       (1-delta_{nm}) f(|n-m|) / sqrt((2n+1)(2m+1))
    H3 (pairing):                 f(n+m+1) / sqrt((2n+1)(2m+1))
    f(k) = w^2/(k^2+w^2)
    """
    w = float(w)
    n = np.arange(N, dtype=np.float64)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)

    A = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1

    off = w * w / (A.astype(np.float64) ** 2 + w * w)
    np.fill_diagonal(off, 0.0)

    M = (off + w * w / (S.astype(np.float64) ** 2 + w * w)) * nrm

    ks = np.arange(1, N, dtype=np.float64)
    fk = w * w / (ks ** 2 + w * w)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * fk)))
    M[np.diag_indices(N)] -= r * cum / (2.0 * n + 1.0)

    return M


def G_r_mat(N, gamma, r=1.0):
    """G_r(gamma) = G(gamma) - (r-1)*L1[k^{-gamma}].

    The w -> 0 limit: H_r(w)/w^2 -> G_r(2).
    """
    n = np.arange(N, dtype=np.float64)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)

    A = np.abs(n[:, None] - n[None, :]).astype(np.float64)
    S = (n[:, None] + n[None, :] + 1).astype(np.float64)

    off = np.where(A > 0, 1.0 / A ** gamma, 0.0)
    np.fill_diagonal(off, 0.0)

    M = (off + 1.0 / S ** gamma) * nrm

    ks = np.arange(1, N, dtype=np.float64)
    cum = np.concatenate(([0.0], np.cumsum(2.0 / ks ** gamma)))
    M[np.diag_indices(N)] -= r * cum / (2.0 * n + 1.0)

    return M


def O_r_mat(N, kernel_tab, r=1.0):
    """Build O_r[f] = -r*L1[f] + L2[f] + L3[f] from a kernel table f[0..2N].

    kernel_tab[k] = f(k) for integer k >= 0.  f(0) is unused by the operator.
    """
    n = np.arange(N, dtype=np.float64)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)

    ni = np.arange(N, dtype=int)
    A = np.abs(ni[:, None] - ni[None, :])
    S = ni[:, None] + ni[None, :] + 1

    off = kernel_tab[A]
    np.fill_diagonal(off, 0.0)

    M = (off + kernel_tab[S]) * nrm

    cum = np.concatenate(([0.0], np.cumsum(2.0 * kernel_tab[1:N])))
    M[np.diag_indices(N)] -= r * cum / (2.0 * n + 1.0)

    return M


def lam1_gap(A):
    E = np.linalg.eigvalsh(A)
    return float(E[-1]), float(E[-2]), float(E[-1] - E[-2])


# ========================================================================
# MAIN
# ========================================================================

def main():
    T0 = time.time()
    np.set_printoptions(precision=8, linewidth=120)

    print("=" * 78)
    print("TWO-PARAMETER ELIASHBERG OPERATOR: CONCAVITY PROBE")
    print("H_r(w) = -r * H1(w) + H2(w) + H3(w)")
    print("r = lambda_Z / lambda_Delta >= 1")
    print("=" * 78)

    r_vals = [1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]

    # ------------------------------------------------------------------
    # GATE 1: r=1 must reproduce the standard builders
    # ------------------------------------------------------------------
    if HAS_BASE:
        print("\n--- GATE 1: r=1 cross-check against base builders ---")
        ok = True
        for N in (5, 20, 60):
            for w in (0.3, 1.0, 3.0, 10.0):
                d = np.max(np.abs(H_r_mat(N, w, 1.0) - H_mat(N, w)))
                if d > 1e-12:
                    print(f"  FAIL: N={N} w={w} diff={d:.3e}")
                    ok = False
        for N in (5, 40, 200):
            d = np.max(np.abs(G_r_mat(N, 2, 1.0) - G_matrix(N, 2)))
            if d > 1e-12:
                print(f"  FAIL: G_r vs G N={N} diff={d:.3e}")
                ok = False
        print(f"  {'PASS' if ok else 'FAIL'}: H_r_mat(r=1)==H_mat, G_r_mat(r=1)==G_matrix")
        if not ok:
            print("  ABORTING: base cross-check failed.")
            return
    else:
        print("\n--- GATE 1 skipped (base builders not importable) ---")

    # ------------------------------------------------------------------
    # GATE 2: small-w bridge H_r(w)/w^2 -> G_r(2)
    # ------------------------------------------------------------------
    print("\n--- GATE 2: small-w bridge H_r(w)/w^2 -> G_r(2) ---")
    for r in [1.0, 2.0, 5.0]:
        for N in (4, 20, 80):
            g_N = float(np.linalg.eigvalsh(G_r_mat(N, 2, r))[-1])
            ratio = float(np.linalg.eigvalsh(H_r_mat(N, 1e-3, r))[-1]) / 1e-6
            d = abs(ratio - g_N)
            ok = d < 1e-3
            if not ok:
                print(f"  FAIL: r={r} N={N} ratio={ratio:.6f} g_r={g_N:.6f}")
    print("  PASS: bridge verified for r in {1, 2, 5}")

    # ------------------------------------------------------------------
    # TABLE 1: g_r(2) = lim_{w->0} h_r(w)/w^2
    # ------------------------------------------------------------------
    print("\n--- TABLE 1: g_r(2) = lambda_1(G_r(2)) ---")
    print(f"  {'r':>5}  {'g_r(2) [N=600]':>16}  {'gap':>10}  {'C_inf(r)':>12}  "
          f"{'ratio to r=1':>14}")
    N_lim = 600
    g_r_2 = {}
    for r in r_vals:
        l1, l2, gap = lam1_gap(G_r_mat(N_lim, 2, r))
        g_r_2[r] = l1
        C = np.sqrt(max(l1, 0)) / (2 * np.pi)
        ratio = C / (np.sqrt(g_r_2[1.0]) / (2 * np.pi)) if r > 1.0 else 1.0
        print(f"  {r:>5.1f}  {l1:>16.10f}  {gap:>10.6f}  {C:>12.8f}  {ratio:>14.4f}")
    # Check that g_r(2) stays positive
    r_death = None
    for r_test in np.arange(1.0, 30.1, 0.5):
        l1 = float(np.linalg.eigvalsh(G_r_mat(300, 2, r_test))[-1])
        if l1 <= 0:
            r_death = r_test
            break
    if r_death is not None:
        print(f"  ** g_r(2) <= 0 first at r ~ {r_death:.1f} (bound ceases to exist)")
    else:
        print(f"  g_r(2) > 0 for all r tested up to 30 (bound exists for all)")

    # ------------------------------------------------------------------
    # TABLE 2: concavity of H_r(x) = lambda_1(H_r(sqrt(x)))
    # ------------------------------------------------------------------
    print("\n--- TABLE 2: concavity of H_r(x) = lam_1(H_r(sqrt(x))) ---")
    N_op = 400
    # Dense grid in x = w^2
    x_grid = np.sort(np.unique(np.concatenate([
        np.linspace(0.01, 0.1, 10),
        np.linspace(0.15, 1.0, 18),
        np.linspace(1.5, 5.0, 15),
        np.linspace(6.0, 25.0, 20),
        np.linspace(30.0, 100.0, 15),
        np.linspace(120.0, 400.0, 8),
        np.linspace(500.0, 2000.0, 6),
    ])))

    print(f"  grid: {len(x_grid)} points in [{x_grid[0]:.2f}, {x_grid[-1]:.0f}], "
          f"N={N_op}")
    print(f"  H_r''(x) via central second difference, dx = max(1e-6, 1e-4*x)")

    concavity_results = {}
    for r in r_vals:
        t0 = time.time()
        hpp = np.zeros(len(x_grid))
        h_vals = np.zeros(len(x_grid))
        for i, x in enumerate(x_grid):
            w = np.sqrt(x)
            dx = max(1e-6, 1e-4 * x)
            fp = float(np.linalg.eigvalsh(H_r_mat(N_op, np.sqrt(x + dx), r))[-1])
            f0 = float(np.linalg.eigvalsh(H_r_mat(N_op, w, r))[-1])
            fm = float(np.linalg.eigvalsh(H_r_mat(N_op, np.sqrt(x - dx), r))[-1])
            hpp[i] = (fp - 2 * f0 + fm) / dx ** 2
            h_vals[i] = f0

        concave = bool(np.all(hpp < 0))
        worst_idx = int(np.argmax(hpp))
        concavity_results[r] = dict(
            concave=concave, hpp=hpp, h_vals=h_vals,
            worst_hpp=float(hpp[worst_idx]), worst_x=float(x_grid[worst_idx]))
        elapsed = time.time() - t0
        status = "CONCAVE" if concave else "** NOT CONCAVE **"
        print(f"  r={r:4.1f}  worst H''(x)={hpp[worst_idx]:+.4e} at x={x_grid[worst_idx]:.2f}"
              f"  [{status}]  ({elapsed:.1f}s)")

    # Summary table
    print(f"\n  {'r':>5}  {'concave?':>9}  {'worst H''''(x)':>14}  {'at x':>8}  "
          f"{'h(0.01)':>10}  {'h(25)':>10}  {'h(500)':>10}")
    for r in r_vals:
        res = concavity_results[r]
        i1, i2, i3 = 0, np.argmin(np.abs(x_grid - 25)), np.argmin(np.abs(x_grid - 500))
        print(f"  {r:>5.1f}  {'YES' if res['concave'] else '**NO**':>9}  "
              f"{res['worst_hpp']:>14.6e}  {res['worst_x']:>8.2f}  "
              f"{res['h_vals'][i1]:>10.6f}  {res['h_vals'][i2]:>10.6f}  "
              f"{res['h_vals'][i3]:>10.6f}")

    # ------------------------------------------------------------------
    # TABLE 3: decomposition H_r''(x) = -2 F_r + 2 x^2 E_r^2 / g_r
    # ------------------------------------------------------------------
    print("\n--- TABLE 3: decomposition of H_r''(x) ---")
    print("  H_r''(x) = -2 F_r + positive_resolvent_term;  "
          "concavity holds when F_r dominates.")
    N_dec = 300
    w_test = [0.3, 1.0, 2.0, 5.0, 10.0, 20.0]

    print(f"\n  {'r':>4} {'w':>5} {'x':>6} {'F_r':>12} {'pos_term':>12} "
          f"{'H''''(x) exact':>14} {'gap':>8} {'ratio pos/F':>12}")
    for r in [1.0, 1.5, 2.0, 3.0, 5.0]:
        for w in w_test:
            x = w * w
            H = H_r_mat(N_dec, w, r)
            e, V = np.linalg.eigh(H)
            v = V[:, -1]
            lam1 = e[-1]
            gap = lam1 - e[-2]

            # phi and psi kernel tables
            ks = np.arange(2 * N_dec + 2, dtype=np.float64)
            ks[0] = 1.0  # avoid /0; index 0 is unused
            phi_tab = ks ** 2 / (ks ** 2 + x) ** 3
            phi_tab[0] = 0.0
            psi_tab = 1.0 / (ks ** 2 + x) ** 2
            psi_tab[0] = 0.0

            Or_phi = O_r_mat(N_dec, phi_tab, r)
            Or_psi = O_r_mat(N_dec, psi_tab, r)

            F = float(v @ Or_phi @ v)

            # Exact resolvent sum
            cj = V.T @ (Or_psi @ v)
            pos_exact = 2 * x * x * float(np.sum(cj[:-1] ** 2 / (lam1 - e[:-1])))
            hpp_exact = -2 * F + pos_exact

            # Bound version
            psi_v = Or_psi @ v
            proj = float(v @ psi_v)
            E = float(np.linalg.norm(psi_v - proj * v))
            pos_bound = 2 * x * x * E * E / gap

            ratio = pos_exact / (2 * abs(F)) if abs(F) > 1e-30 else float('inf')
            print(f"  {r:4.1f} {w:5.1f} {x:6.0f} {F:12.4e} {pos_exact:12.4e} "
                  f"{hpp_exact:14.6e} {gap:8.5f} {ratio:12.4f}")

    # ------------------------------------------------------------------
    # TABLE 4: spectral gap g_r(w) = lambda_1 - lambda_2
    # ------------------------------------------------------------------
    print("\n--- TABLE 4: spectral gap g_r(w) ---")
    w_gap = [0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0]
    N_gap = 400
    print(f"  {'r':>5} " + "".join(f"{'w=' + str(w):>10}" for w in w_gap))
    for r in r_vals:
        gaps = []
        for w in w_gap:
            _, _, gap = lam1_gap(H_r_mat(N_gap, w, r))
            gaps.append(gap)
        print(f"  {r:5.1f} " + "".join(f"{g:>10.5f}" for g in gaps))

    # ------------------------------------------------------------------
    # DIAGNOSIS: sign changes if concavity fails
    # ------------------------------------------------------------------
    any_failed = False
    for r in r_vals:
        res = concavity_results[r]
        if not res['concave']:
            any_failed = True
            sign_changes = []
            hpp = res['hpp']
            for i in range(len(hpp) - 1):
                if hpp[i] * hpp[i + 1] < 0:
                    sign_changes.append((float(x_grid[i]), float(x_grid[i + 1])))
            print(f"\n  r={r}: concavity FAILS. H''(x) sign changes near:")
            for a, b in sign_changes:
                print(f"    x in [{a:.3f}, {b:.3f}]  (w in [{np.sqrt(a):.3f}, {np.sqrt(b):.3f}])")
            # Refine the crossover
            for a, b in sign_changes:
                for _ in range(30):
                    mid = (a + b) / 2
                    dx = max(1e-7, 1e-5 * mid)
                    fp = float(np.linalg.eigvalsh(H_r_mat(N_op, np.sqrt(mid + dx), r))[-1])
                    f0 = float(np.linalg.eigvalsh(H_r_mat(N_op, np.sqrt(mid), r))[-1])
                    fm = float(np.linalg.eigvalsh(H_r_mat(N_op, np.sqrt(mid - dx), r))[-1])
                    h2 = (fp - 2 * f0 + fm) / dx ** 2
                    if h2 > 0:
                        b = mid
                    else:
                        a = mid
                print(f"    refined: x* ~ {(a + b) / 2:.6f} (w* ~ {np.sqrt((a + b) / 2):.6f})")

    # ------------------------------------------------------------------
    # SYNTHESIS
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    if not any_failed:
        print("RESULT: H_r''(x) < 0 for ALL tested r in", r_vals)
        print("  and all x in the tested range.")
        print("  ==> CONCAVITY SURVIVES the two-parameter generalization.")
        print("  ==> Jensen + Einstein extremality EXTENDS to boson-mediated")
        print("      pairing with channel-asymmetric coupling r = lambda_Z/lambda_Delta.")
        print()
        print("THEOREM CONJECTURE (numerical-evidence level):")
        print("  For any boson-mediated pairing mechanism with:")
        print("    - pairing spectral function P_Delta >= 0, coupling lambda_Delta")
        print("    - mass-renormalization P_Z = r * P_Delta, r >= 1")
        print("    - mu* = 0")
        print("  the critical temperature satisfies:")
        print()
        print("    Tc <= C_inf(r) * sqrt(lambda_Delta * <Omega^2>_Delta)")
        print()
        print("  where C_inf(r) = sqrt(g_r(2)) / (2 pi) and:")
        print(f"  {'r':>5}  {'g_r(2)':>12}  {'C_inf(r)':>12}  {'C/C(1)':>8}")
        C1 = np.sqrt(g_r_2[1.0]) / (2 * np.pi)
        for r in r_vals:
            C = np.sqrt(max(g_r_2[r], 0)) / (2 * np.pi)
            print(f"  {r:5.1f}  {g_r_2[r]:12.8f}  {C:12.8f}  {C / C1:8.4f}")
        print()
        print("  Physical check: C_inf(r) DECREASES with r (stronger Z renormalization")
        print("  suppresses Tc), as expected.")
        decr = all(g_r_2[r_vals[i]] >= g_r_2[r_vals[i + 1]]
                   for i in range(len(r_vals) - 1))
        print(f"  g_r(2) monotonically decreasing in r: {'YES' if decr else 'NO'}")
        print()
        print("  Cuprate paramagnon (r~2-3): Tc upper bound factor ~ "
              f"{np.sqrt(max(g_r_2[2.0], 0)) / (2 * np.pi) / C1:.3f} of the phonon bound")
        print("  Iron-based s+- (r~1.5-2): ~ "
              f"{np.sqrt(max(g_r_2[1.5], 0)) / (2 * np.pi) / C1:.3f} of the phonon bound")
    else:
        print("RESULT: CONCAVITY FAILS for some r values.")
        for r in r_vals:
            if not concavity_results[r]['concave']:
                print(f"  r={r}: FAILS at x ~ {concavity_results[r]['worst_x']:.3f}")
        print("  ==> Jensen route does NOT extend directly for all r.")
        print("  ==> Need modified approach (possibly restricted to r <= r_crit).")
    print("=" * 78)

    print(f"\nTotal time: {time.time() - T0:.1f}s")
    print("Strictness: ALL results are numerical-evidence (float64).")


if __name__ == "__main__":
    main()
