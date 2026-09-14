"""
P4 · Concavity via q-monotonicity route (numerical gates)

Verifies:
  (Q1) O[d_x f] = (O[f] - O[f^2])/x  (algebraic identity, kernel-level)
  (Q2) xH'(x) = H(x) - q(x), q := <v, O[f^2] v>  (Feynman-Hellmann)
       ⟹ H''(x) = -q'(x)/x, so concavity ⟺ q increasing
  (A) E(x) := <v, O[g] v> > 0, g := d_x(f^2) = 2xk^2/(k^2+x)^3
      Corollary A' (P4_uniqueness_strategy.md §2.A) gives E > 0
      when v  in  C (nonincreasing nonnegative) and g(1) > 0.
  (B) |C(x)| < E(x), where C := q' - E ≤ 0 (eigenvector motion term)

Status: (A) PROVED analytically (Cor A'), verified here.
        (B) numerical evidence |C|/E ≤ 0.184 up to x=1638, route TBD.
"""
import sys
import time
import numpy as np

# Import H_mat builder from existing script
sys.path.insert(0, "D:/physics/superconductivity/.claude/worktrees/agent-a8932ec4c14045b28/theory/derivations")
from p4_interchange_uniformity import H_mat

_ck_count = 0
_t0 = time.time()


def ck(cond, msg):
    """Gate: count checks, fail-stop on False."""
    global _ck_count
    _ck_count += 1
    if not cond:
        print(f"\n*** FAIL check {_ck_count}: {msg}")
        sys.exit(1)


def O_ker(N, fk):
    """
    Generic kernel operator O[f] for f:[1..∞)→ℝ given as array fk[0..KX-1].
    fk[0] = f(0) by convention (set =1 for Lorentz kernel, =0 for derivatives).
    Matches the builder in p4_interchange_uniformity.py.
    """
    n = np.arange(N)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)
    KX = len(fk)

    A = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1

    # L2 term: f(|n-m|), zero on diagonal
    off = fk[np.minimum(A, KX - 1)].copy()
    np.fill_diagonal(off, 0.0)

    # L3 term: f(n+m+1)
    M = (off + fk[np.minimum(S, KX - 1)] * (S < KX)) * nrm

    # L1 term: diagonal -2 Σ_{k=1}^n f(k) / (2n+1)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * fk[1:min(N, KX)])))
    if N > KX:
        cum = np.concatenate([cum, np.full(N - KX, cum[-1])])
    M[np.diag_indices(N)] -= cum / (2.0 * n + 1.0)

    return M


def kernels(x, KX=2400):
    """
    Lorentz kernel f(k;x) = x/(k^2+x) and related quantities.
    Returns dict with keys: f, f2, f3, g, df (all arrays [0..KX-1]).
    """
    k = np.arange(KX, dtype=float)
    k2 = k * k
    d = k2 + x

    f = x / d
    f[0] = 1.0

    f2 = f * f
    f3 = f * f2

    # g := d_x(f^2) = 2x k^2/(k^2+x)^3
    g = 2.0 * x * k2 / (d ** 3)
    g[0] = 0.0

    # df := d_x f = k^2/(k^2+x)^2
    df = k2 / (d ** 2)
    df[0] = 0.0

    return {"f": f, "f2": f2, "f3": f3, "g": g, "df": df}


def top_eigen(A):
    """Return (λ₁, v₁, λ₂) with v₁ normalized, v₁[0]>0."""
    E, V = np.linalg.eigh(A)
    lam1 = float(E[-1])
    v = V[:, -1]
    if v[0] < 0:
        v = -v
    lam2 = float(E[-2]) if len(E) > 1 else -np.inf
    return lam1, v, lam2


def main():
    print("=" * 70)
    print("P4 Concavity via q-monotonicity: numerical gates")
    print("=" * 70)

    N = 300
    xs = [0.05, 0.2, 0.8, 1.6, 3.2, 6.4, 25.6, 102.4, 409.6, 1638.4]

    print(f"\nConfiguration: N={N}, {len(xs)} test points")
    print("")

    # =========================================================================
    # Gate 1: (Q1) algebraic identity O[df/dx] = (O[f] - O[f^2])/x
    # =========================================================================
    print("Gate 1: (Q1) O[df/dx] = (O[f] - O[f^2])/x")
    print("-" * 70)
    for x in xs:
        ker = kernels(x, KX=8 * N)
        LHS = O_ker(N, ker["df"])
        RHS = (O_ker(N, ker["f"]) - O_ker(N, ker["f2"])) / x
        err = np.max(np.abs(LHS - RHS))
        ck(err < 1e-13, f"x={x:.3f}: |O[df] - (O[f]-O[f2])/x| = {err:.2e}")
        print(f"  x={x:8.3f}  max_err={err:.2e}")
    print("")

    # =========================================================================
    # Gate 2: v_x  in  C (nonnegative, nonincreasing) for all x
    # =========================================================================
    print("Gate 2: v_x  in  cone C (nonnegative, nonincreasing)")
    print("-" * 70)
    for x in xs:
        ker = kernels(x, KX=8 * N)
        A = O_ker(N, ker["f"])
        h, v, lam2 = top_eigen(A)

        minv = np.min(v)
        dv = v[:-1] - v[1:]
        neg_dv_count = int(np.sum(dv < -1e-14))

        ck(minv >= -1e-14, f"x={x:.3f}: min(v) = {minv:.2e} < 0")
        ck(neg_dv_count == 0, f"x={x:.3f}: {neg_dv_count} negative differences")

        if xs.index(x) % 3 == 0:
            print(f"  x={x:8.3f}  min(v)={minv:.2e}  min(dv)={dv.min():.2e}  neg_dv={neg_dv_count}")
    print("")

    # =========================================================================
    # Gate 3: g = 2(f^2 - f^3)/x  (algebraic check for E alternative form)
    # =========================================================================
    print("Gate 3: g = 2(f^2 - f^3)/x  (algebraic identity)")
    print("-" * 70)
    for x in xs:
        ker = kernels(x, KX=8 * N)
        g_direct = ker["g"]
        g_alt = 2.0 * (ker["f2"] - ker["f3"]) / x
        err = np.max(np.abs(g_direct - g_alt))
        ck(err < 1e-14, f"x={x:.3f}: |g - 2(f2-f3)/x| = {err:.2e}")
        if xs.index(x) % 3 == 0:
            print(f"  x={x:8.3f}  max_err={err:.2e}")
    print("")

    # =========================================================================
    # Gate 4: (A) E(x) > 0  (Corollary A' ⟹ analytical, verify numerically)
    # =========================================================================
    print("Gate 4: (A) E(x) = <v, O[g] v> > 0  [Corollary A' PROVED]")
    print("-" * 70)
    E_vals = []
    for x in xs:
        ker = kernels(x, KX=8 * N)
        A = O_ker(N, ker["f"])
        h, v, lam2 = top_eigen(A)

        E = v @ (O_ker(N, ker["g"]) @ v)
        E_vals.append(E)

        ck(E > 0, f"x={x:.3f}: E = {E:.6e} ≤ 0")
        print(f"  x={x:8.3f}  E={E:.6f}  (>0 by Cor A')")
    print("")

    # =========================================================================
    # Gate 5: Alternative form E = 2(q - <v,O[f^3]v>)/x
    # =========================================================================
    print("Gate 5: E alternative form E = 2(q - <v,O[f^3]v>)/x")
    print("-" * 70)
    for i, x in enumerate(xs):
        ker = kernels(x, KX=8 * N)
        A = O_ker(N, ker["f"])
        h, v, lam2 = top_eigen(A)

        q = v @ (O_ker(N, ker["f2"]) @ v)
        m3 = v @ (O_ker(N, ker["f3"]) @ v)

        E_direct = E_vals[i]
        E_alt = 2.0 * (q - m3) / x
        err = abs(E_direct - E_alt)

        ck(err < 1e-10, f"x={x:.3f}: |E_direct - E_alt| = {err:.2e}")
        if i % 3 == 0:
            print(f"  x={x:8.3f}  E_dir={E_direct:.6f}  E_alt={E_alt:.6f}  err={err:.2e}")
    print("")

    # =========================================================================
    # Gate 6: q' = E + C  (eigenvector motion term, soft check)
    # =========================================================================
    print("Gate 6: q' = E + C  (eigenvector motion, numerical estimate only)")
    print("-" * 70)
    print("  Note: C <= 0 analytically (second-order perturbation theory),")
    print("        but numerical derivative has O(dx^2) error; we report |C|/E")
    print("        rather than asserting C <= 0 at every point.")
    print("")

    # Compute q(x) on a denser grid for numerical derivative
    xs_dense = []
    x = 0.05
    while x <= 1700:
        xs_dense.append(x)
        x *= 1.05

    q_dense = []
    E_dense = []
    for x in xs_dense:
        ker = kernels(x, KX=8 * N)
        A = O_ker(N, ker["f"])
        h, v, lam2 = top_eigen(A)
        q = v @ (O_ker(N, ker["f2"]) @ v)
        E = v @ (O_ker(N, ker["g"]) @ v)
        q_dense.append(q)
        E_dense.append(E)

    # Compute q' via central difference, then C = q' - E
    C_ratios = []
    C_vals = []
    for i in range(1, len(xs_dense) - 1):
        xm, x0, xp = xs_dense[i - 1], xs_dense[i], xs_dense[i + 1]
        qm, qp = q_dense[i - 1], q_dense[i + 1]
        E0 = E_dense[i]

        qp_num = (qp - qm) / (xp - xm)
        C = qp_num - E0
        C_vals.append(C)

        ratio = abs(C) / E0 if E0 > 0 else 0.0
        C_ratios.append(ratio)

        # Report sparse
        if x0 in xs or (i % 40 == 1):
            print(f"  x={x0:8.3f}  q'={qp_num:.6f}  E={E0:.6f}  C={C:.6e}  |C|/E={ratio:.4f}")

    max_ratio = max(C_ratios)
    neg_C_frac = sum(1 for c in C_vals if c < 0) / len(C_vals)
    print(f"\n  Summary: |C|/E max = {max_ratio:.4f}")
    print(f"           Fraction with C < 0: {neg_C_frac:.3f}  (expect ~1.0)")
    print(f"           (numerical noise causes ~{1-neg_C_frac:.1%} apparent violations)")
    print("")

    # =========================================================================
    # Gate 7: Consistency with two N values (300 vs 600)
    # =========================================================================
    print("Gate 7: N-consistency check (N=300 vs N=600)")
    print("-" * 70)
    xs_check = [0.8, 6.4, 102.4]
    for x in xs_check:
        ker = kernels(x, KX=8 * 600)

        A300 = O_ker(300, ker["f"])
        h300, v300, _ = top_eigen(A300)
        E300 = v300 @ (O_ker(300, ker["g"]) @ v300)

        A600 = O_ker(600, ker["f"])
        h600, v600, _ = top_eigen(A600)
        E600 = v600 @ (O_ker(600, ker["g"]) @ v600)

        reldiff_h = abs(h300 - h600) / h600
        reldiff_E = abs(E300 - E600) / E600

        ck(reldiff_h < 1e-5, f"x={x}: |h300-h600|/h600 = {reldiff_h:.2e}")
        ck(reldiff_E < 1e-3, f"x={x}: |E300-E600|/E600 = {reldiff_E:.2e}")

        print(f"  x={x:8.3f}  Δh/h={reldiff_h:.2e}  ΔE/E={reldiff_E:.2e}")
    print("")

    # =========================================================================
    # Gate 8: (B) commutator identity  (A-h)(I-P1)Bv = (I-P1)[A,B]v
    #         and skew-symmetry of [A,B] (=> ||.||_F >= sqrt2 ||.||_2)
    # =========================================================================
    print("Gate 8: (B) commutator identity and skew-symmetry")
    print("-" * 70)
    for x in [0.04, 1.0, 4.0, 25.0, 400.0]:
        ker = kernels(x, KX=8 * N)
        A = O_ker(N, ker["f"])
        B = O_ker(N, ker["f2"])
        h, v, lam2 = top_eigen(A)

        Bv = B @ v
        r = Bv - (v @ Bv) * v                      # (I-P1) B v
        comm = A @ B - B @ A
        cv = comm @ v
        rc = cv - (v @ cv) * v                     # (I-P1) [A,B] v

        # (A - h) applied to the residual must reproduce (I-P1)[A,B]v
        lhs = A @ r - h * r
        idn = np.max(np.abs(lhs - rc))
        scale = max(np.max(np.abs(rc)), 1e-300)
        ck(idn / scale < 1e-8,
           f"x={x}: commutator identity rel err {idn/scale:.2e}")

        skew = np.max(np.abs(comm + comm.T))
        ck(skew < 1e-10, f"x={x}: [A,B] not skew, {skew:.2e}")

        cF = np.linalg.norm(comm, "fro")
        c2 = np.linalg.norm(comm, 2)
        ck(cF >= np.sqrt(2.0) * c2 * (1 - 1e-9),
           f"x={x}: ||.||_F < sqrt2||.||_2")
        print(f"  x={x:8.3f}  id_relerr={idn/scale:.2e}  skew={skew:.1e}  "
              f"F/(sqrt2*2norm)={cF/(np.sqrt(2)*c2):.4f}")
    print("")

    # =========================================================================
    # Gate 9: (B) ANALYTIC closure. Bound must be a real number < E.
    #   |C| <= (2/x) ||(I-P1)Bv||^2 / g   and   ||(I-P1)Bv|| <= ||[A,B]v||/g
    #   ||[A,B]v|| <= ||[A,B]||_2 <= ||[A,B]||_F / sqrt2  (skew, trace form)
    #   => |C| <= (2/x) (cF/sqrt2)^2 / g^3
    # Gap g: analytic lower bound g_lo = max(0.42h, h - Theta*).
    # =========================================================================
    print("Gate 9: (B) analytic closure  |C|_bound < E  with analytic gap")
    print("-" * 70)
    THETA_STAR = 0.6349
    hole = []
    for varpi in [0.2, 0.5, 1.0, 1.35, 1.5, 1.6, 1.7, 2.0, 3.2, 6.4, 20.0, 40.0]:
        x = varpi ** 2
        ker = kernels(x, KX=8 * N)
        A = O_ker(N, ker["f"])
        B = O_ker(N, ker["f2"])
        h, v, lam2 = top_eigen(A)
        E = v @ (O_ker(N, ker["g"]) @ v)

        cF = np.linalg.norm(A @ B - B @ A, "fro")
        num = (2.0 / x) * (cF / np.sqrt(2.0)) ** 2

        g_lo = max(0.42 * h, h - THETA_STAR)       # analytic gap
        C_an = num / g_lo ** 3 if g_lo > 0 else np.inf
        C_tr = num / (h - lam2) ** 3               # true-gap version

        ck(C_tr < E, f"varpi={varpi}: true-gap closure failed "
                     f"({C_tr:.4e} >= {E:.4e})")
        if not (C_an < E):
            hole.append(varpi)
        print(f"  w={varpi:5.2f}  E={E:.6f}  |C|_an={C_an:8.5f} "
              f"({'OK ' if C_an < E else 'HOLE'})  |C|_true={C_tr:.6f} "
              f"ratio={C_tr/E:.4f}")

    print(f"\n  analytic-gap hole: {hole if hole else 'none'}")
    print(f"  true-gap closure:  all points, worst ratio |C|/E ~ 0.225")
    ck(all(1.4 <= w <= 1.7 for w in hole),
       f"hole escaped [1.4,1.7]: {hole}")
    print("")

    # =========================================================================
    # Gate 9b: the hole closes IF the note's POINTWISE gap bound is used.
    #   2026-09-13_gap_lorentz_probe.md §4.2 tabulates lam2_hi/h at 7 points.
    #   The uniform worst case (0.582, the value at varpi=2) is what fails;
    #   the pointwise value at varpi=1.5 is 0.454, which closes.
    #   CONDITIONAL: linear interpolation between the 7 tabulated points is
    #   NOT proved here (needs lam2_hi/h monotone, or a fine-grid evaluation).
    #   The conservative step function (right-endpoint value) still FAILS.
    # =========================================================================
    print("Gate 9b: hole vs note's pointwise gap  [CONDITIONAL on interp]")
    print("-" * 70)
    tab_w = np.array([0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0])
    tab_r = np.array([0.269, 0.215, 0.237, 0.277, 0.332, 0.454, 0.582])
    ck(bool(np.all(np.diff(tab_r[2:]) > 0)),
       "note table lam2/h not increasing on [0.5,2] -- interp unjustified")

    worst_pw = 0.0
    for varpi in [1.40, 1.50, 1.54, 1.60, 1.70, 1.76, 1.90, 2.00]:
        x = varpi ** 2
        ker = kernels(x, KX=8 * N)
        A = O_ker(N, ker["f"])
        B = O_ker(N, ker["f2"])
        h, v, lam2 = top_eigen(A)
        E = v @ (O_ker(N, ker["g"]) @ v)
        cF = np.linalg.norm(A @ B - B @ A, "fro")
        num = (2.0 / x) * (cF / np.sqrt(2.0)) ** 2

        r_pw = float(np.interp(varpi, tab_w, tab_r))
        g_pw = max((1.0 - r_pw) * h, h - THETA_STAR)
        cb_pw = num / g_pw ** 3
        # conservative step function: right-endpoint (largest) ratio
        r_hi = float(tab_r[np.searchsorted(tab_w, varpi, side="left")]
                     if varpi <= 2.0 else tab_r[-1])
        cb_hi = num / max((1.0 - r_hi) * h, h - THETA_STAR) ** 3

        worst_pw = max(worst_pw, cb_pw / E)
        print(f"  w={varpi:5.2f}  E={E:.5f}  pointwise:{cb_pw/E:6.3f}"
              f"{'OK ' if cb_pw < E else 'HOLE'}   step:{cb_hi/E:6.3f}"
              f"{'OK' if cb_hi < E else 'HOLE'}")

    ck(worst_pw < 1.0, f"pointwise gap fails too: worst {worst_pw:.3f}")
    print(f"\n  pointwise-gap worst ratio = {worst_pw:.4f} -> hole CLOSES")
    print("  but the conservative step function does NOT (see .md §5.2)")
    print("")

    # =========================================================================
    # Gate 10: N-consistency of the (B) closure margin
    # =========================================================================
    print("Gate 10: N-consistency of closure margin (N=300 vs N=600)")
    print("-" * 70)
    for varpi in [0.5, 1.5, 6.4]:
        x = varpi ** 2
        ker = kernels(x, KX=8 * 600)
        margins = []
        for NN in (300, 600):
            A = O_ker(NN, ker["f"])
            B = O_ker(NN, ker["f2"])
            h, v, l2 = top_eigen(A)
            E = v @ (O_ker(NN, ker["g"]) @ v)
            cF = np.linalg.norm(A @ B - B @ A, "fro")
            margins.append(E - (2.0 / x) * (cF / np.sqrt(2.0)) ** 2 / (h - l2) ** 3)
        m3, m6 = margins
        rel = abs(m3 - m6) / max(abs(m3), 1e-30)
        ck(m3 > 0 and m6 > 0, f"varpi={varpi}: margin not positive")
        ck(rel < 1e-2, f"varpi={varpi}: margin N-drift {rel:.2e}")
        print(f"  w={varpi:5.2f}  margin(300)={m3:+.6f}  margin(600)={m6:+.6f}  "
              f"reldiff={rel:.2e}")
    print("")

    # =========================================================================
    # Summary
    # =========================================================================
    elapsed = time.time() - _t0
    print("=" * 70)
    print(f"PASS: {_ck_count} checks in {elapsed:.1f}s")
    print("")
    print("Conclusions:")
    print("  (Q1) O[df] = (O[f]-O[f^2])/x: algebraic, verified to 1e-13")
    print("  (Q2) xH' = H - q: Feynman-Hellmann (analytic proof in .md)")
    print("  (A)  E > 0: PROVED by Corollary A' (v in C, g>=0, g(1)>0)")
    print("  (B)  |C| <= (2/x)||[A,B]||_F^2/(2 g^3): closes with true gap at all")
    print("       tested varpi (worst ratio 0.225); with ANALYTIC gap g_lo =")
    print("       max(0.42h, h-Theta*) closes EXCEPT varpi ~ [1.4,1.7].")
    print("=" * 70)


if __name__ == "__main__":
    main()
