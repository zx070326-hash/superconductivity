r"""
Independent P_Z / P_Delta probe (2026-09-16), for team-lead's exploration task
"探索1.1：独立双测度问题".

Question: push-universal-bound proved concavity survives for H_r(w) = -r*H1(w)
+ H2(w) + H3(w), i.e. Z-channel and Delta-channel SHARE one frequency shape w
(same measure P, only overall coupling ratio r = lambda_Z/lambda_Delta differs).
Here we DECOUPLE the shapes: Z-channel diagonal built from its own w_Z,
Delta-channel off-diagonal built from its own w_Delta, independently.

  H_two(x_Delta, x_Z; r) = -r*H1(sqrt(x_Z)) + H2(sqrt(x_Delta)) + H3(sqrt(x_Delta))
  h(x_Delta, x_Z; r)     = lambda_1(H_two)

Double-Einstein extremality (max over independent P_Delta, P_Z with fixed
first/second moments) reduces, via the SAME two Jensen steps as the single-
measure Theorem E, to:  is h JOINTLY CONCAVE in (x_Delta, x_Z) on R_+^2?
(Step 1, lambda_1 convex in the affine pair (P_Delta,P_Z), is automatic and
direction-safe regardless of what we find here; only Step 2 needs checking.)

STRICTNESS: numerical-evidence only (float64, ARPACK/eigvalsh, finite-difference
Hessian). No interval arithmetic in this probe.
"""
import sys, os, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 '..', '..', 'derivations'))
from p4_interchange_uniformity import H_mat  # r=1, w_Delta=w_Z cross-check


def H_two_mat(N, w_delta, w_z, r=1.0):
    """H_two(w_delta, w_z; r) = -r*H1(w_z) + H2(w_delta) + H3(w_delta)."""
    w_d = float(w_delta); w_z = float(w_z)
    n = np.arange(N, dtype=np.float64)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)

    A = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1

    off = w_d * w_d / (A.astype(np.float64) ** 2 + w_d * w_d)
    np.fill_diagonal(off, 0.0)
    M = (off + w_d * w_d / (S.astype(np.float64) ** 2 + w_d * w_d)) * nrm

    ks = np.arange(1, N, dtype=np.float64)
    fk = w_z * w_z / (ks ** 2 + w_z * w_z)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * fk)))
    M[np.diag_indices(N)] -= r * cum / (2.0 * n + 1.0)
    return M


def h_val(N, x_delta, x_z, r=1.0):
    A = H_two_mat(N, np.sqrt(x_delta), np.sqrt(x_z), r)
    return float(np.linalg.eigvalsh(A)[-1])


def gate1(N=60):
    """r=1, w_delta=w_z=w must reproduce H_mat exactly."""
    ok = True
    for w in (0.3, 1.0, 3.0, 10.0):
        d = np.max(np.abs(H_two_mat(N, w, w, 1.0) - H_mat(N, w)))
        if d > 1e-12:
            print(f"  FAIL gate1: w={w} diff={d:.3e}")
            ok = False
    print(f"GATE 1 (r=1, w_delta=w_z=w == H_mat): {'PASS' if ok else 'FAIL'}")
    return ok


def hessian_fd(N, x_delta, x_z, r, h=1e-3):
    """Central finite-difference Hessian of h(x_delta,x_z;r) at a point."""
    hd, hz = h * max(x_delta, 1.0), h * max(x_z, 1.0)

    def f(xd, xz):
        return h_val(N, xd, xz, r)

    f00 = f(x_delta, x_z)
    fpp = f(x_delta + hd, x_z + hz)
    fpm = f(x_delta + hd, x_z - hz)
    fmp = f(x_delta - hd, x_z + hz)
    fmm = f(x_delta - hd, x_z - hz)
    fp0 = f(x_delta + hd, x_z)
    fm0 = f(x_delta - hd, x_z)
    f0p = f(x_delta, x_z + hz)
    f0m = f(x_delta, x_z - hz)

    Hdd = (fp0 - 2 * f00 + fm0) / hd ** 2
    Hzz = (f0p - 2 * f00 + f0m) / hz ** 2
    Hdz = (fpp - fpm - fmp + fmm) / (4 * hd * hz)
    return Hdd, Hzz, Hdz


def main():
    T0 = time.time()
    gate1()

    print("\n--- JOINT HESSIAN OF h(x_delta, x_z; r) ---")
    print("  concave point requires: Hdd<0, Hzz<0, Hdd*Hzz - Hdz^2 >= 0")
    N = 400
    r_vals = [1.0, 1.5, 2.0, 3.0]
    # grid of (ratio x_z/x_delta, x_delta scale)
    ratios = [0.1, 0.3, 1.0, 3.0, 10.0]
    scales = [0.5, 2.0, 8.0, 30.0]

    n_fail = 0
    n_total = 0
    worst = None
    for r in r_vals:
        for xd in scales:
            for ratio in ratios:
                xz = xd * ratio
                Hdd, Hzz, Hdz = hessian_fd(N, xd, xz, r)
                det = Hdd * Hzz - Hdz ** 2
                n_total += 1
                joint_concave = (Hdd < 0) and (Hzz < 0) and (det >= -1e-8 * abs(Hdd * Hzz))
                if not joint_concave:
                    n_fail += 1
                    print(f"  NOT jointly concave: r={r} x_d={xd} x_z={xz:.3f} "
                          f"Hdd={Hdd:.4e} Hzz={Hzz:.4e} Hdz={Hdz:.4e} det={det:.4e}")
                if worst is None or det < worst[0]:
                    worst = (det, r, xd, xz, Hdd, Hzz, Hdz)
    print(f"\n  {n_total - n_fail}/{n_total} grid points jointly concave "
          f"(det=HddHzz-Hdz^2 >= 0 & both diag <0)")
    print(f"  worst det: {worst}")

    print(f"\n[{time.time()-T0:.1f}s] done")


if __name__ == "__main__":
    main()


# ========================================================================
# EXPLICIT TWO-POINT MIXTURE TEST (does a non-Einstein P_Z beat double-
# Einstein at fixed lambda_Z, <omega^2>_Z, lambda_Delta, <omega^2>_Delta ?)
# ========================================================================

def H1_only(N, w, r=1.0):
    """-r*H1(w), the pure Z-channel (diagonal) piece."""
    n = np.arange(N, dtype=np.float64)
    M = np.zeros((N, N))
    ks = np.arange(1, N, dtype=np.float64)
    fk = w * w / (ks ** 2 + w * w)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * fk)))
    M[np.diag_indices(N)] = -r * cum / (2.0 * n + 1.0)
    return M


def H23_only(N, w):
    """H2(w)+H3(w), the pure Delta-channel (off-diagonal + Hankel) piece."""
    n = np.arange(N, dtype=np.float64)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)
    A = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    off = w * w / (A.astype(np.float64) ** 2 + w * w)
    np.fill_diagonal(off, 0.0)
    return (off + w * w / (S.astype(np.float64) ** 2 + w * w)) * nrm


def mixture_probe(N=400, x_delta=2.0, x_z_mean=2.0, r=1.0, spread=0.9):
    """P_Delta = Einstein at x_delta (fixed).  P_Z = two-point mixture on
    x_z with mean x_z_mean, vs double-Einstein at x_z_mean.  spread in (0,1)
    controls how far the two-point support is pulled apart (as a fraction
    of the way to the boundary {0, ...}); p is solved so the mean matches
    exactly for a symmetric two-point spread first, then an asymmetric one.
    """
    w_delta = np.sqrt(x_delta)
    D = H23_only(N, w_delta)

    # double-Einstein baseline
    Z0 = H1_only(N, np.sqrt(x_z_mean), r)
    kappa_dE = float(np.linalg.eigvalsh(Z0 + D)[-1])

    results = []
    # symmetric two-point mixture: x_z in {x_z_mean*(1-spread), x_z_mean*(1+spread)}, p=1/2
    a = x_z_mean * (1 - spread)
    b = x_z_mean * (1 + spread)
    Za = H1_only(N, np.sqrt(max(a, 1e-8)), r)
    Zb = H1_only(N, np.sqrt(b), r)
    Zmix = 0.5 * Za + 0.5 * Zb
    kappa_mix = float(np.linalg.eigvalsh(Zmix + D)[-1])
    results.append(("symmetric p=1/2", a, b, 0.5, kappa_mix))

    # asymmetric: p at 0+ (mass near omega->0, x_z~eps), (1-p) at c so that
    # p*eps + (1-p)*c = x_z_mean  (extreme two-point, mimics a low-frequency
    # peak + a high-frequency peak -- the kind of shape a real phonon DOS has)
    for p in (0.1, 0.3, 0.5, 0.7, 0.9):
        eps = 1e-6
        c = (x_z_mean - p * eps) / (1 - p)
        Ze = H1_only(N, np.sqrt(eps), r)
        Zc = H1_only(N, np.sqrt(c), r)
        Zmix2 = p * Ze + (1 - p) * Zc
        kappa_mix2 = float(np.linalg.eigvalsh(Zmix2 + D)[-1])
        results.append((f"eps+high p={p}", eps, c, p, kappa_mix2))

    print(f"\n--- MIXTURE PROBE: x_delta={x_delta} (Einstein, fixed), "
          f"x_z mean={x_z_mean}, r={r} ---")
    print(f"  double-Einstein kappa = {kappa_dE:.10f}")
    for name, a_, b_, p_, k_ in results:
        beats = "BEATS double-Einstein" if k_ > kappa_dE + 1e-10 else "does not beat"
        print(f"  {name:>18s}: a={a_:.4g} b={b_:.4g} p={p_:.2f}  "
              f"kappa={k_:.10f}  delta={k_-kappa_dE:+.3e}  [{beats}]")


if __name__ == "__main__":
    for r in (1.0, 2.0):
        for xd in (0.5, 2.0, 8.0):
            mixture_probe(N=400, x_delta=xd, x_z_mean=xd, r=r, spread=0.9)
