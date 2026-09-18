"""[numerical-evidence] Structural audit of the LINEARIZED vertex kernel K^V_Delta,
against the four "incompatibility" claims of vertex_sign_theorem_recon_2026-09-17.md
section 3.

The recon report says the vertex kernel cannot enter P4's operator class O[f] because
  (1) it is a double sum, not a single integral over P(domega);
  (2) "even linearized at Delta->0 the kernel is still a FOUR-INDEX tensor
      lambda_{m-m'} lambda_{m'-m''}, not a two-index f(|m-m'|) kernel";
  (3) the Z-channel and Delta-channel kernels differ, while O[f] uses one f;
  (4) "the vertex kernel is no longer a Pick-function integral, because
      gamma = Delta/sqrt(omega^2+Delta^2) is NOT concave in Delta -- its second
      derivative ~ -3 Delta/omega^2 is negative for Delta>0 but its overall sign
      depends on omega", so kernel-level Jensen does not apply.

Claims (2) and (4) are checked here directly.  (4) is pure calculus (see the sympy
one-liner in the session log): d2/dDelta2 [Delta/sqrt(w^2+D^2)] = -3 D w^2/(w^2+D^2)^{5/2},
which is < 0 for EVERY w != 0 and every D > 0 -- gamma^Delta is strictly concave on
Delta >= 0, uniformly in omega.  The function whose convexity DOES flip with Delta is
gamma^omega = omega/sqrt(w^2+D^2), with d2 = w(2D^2-w^2)/(w^2+D^2)^{5/2}, sign change
at D = w/sqrt(2).  So claim (4) names the wrong function; and P4's Jensen is over the
PHONON measure P(domega) in F(k) = int varpi^2/(k^2+varpi^2) P(domega), i.e. in the
variable x = varpi^2, not in Delta at all -- so it also names the wrong variable.

Claim (2) is checked here numerically: K^V_Delta[m,k] = -lam_{m-k}(C(m,k)+C(k,m)) +
D(m,k) is a genuine MATRIX (the diagram's free index is summed inside C and D).  What
we test is which of P4's structural hypotheses it actually violates:
  H1  Toeplitz:      K[m,k] depends only on m-k ?
  H2  symmetry:      K = K^T ?
  H3  elementwise sign-definiteness (O[f] has L2/L3 >= 0 offdiagonal) ?
  H4  monotone decay in |m-k| along a row ?

Run: OMP_NUM_THREADS=1 python _chk_vertex_kernel_structure.py
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "derivations")))
import p2_schrodi_fsr as sf


def toeplitz_defect(K):
    """max deviation of K from the best Toeplitz matrix with the same diagonals."""
    m = K.shape[0]
    worst, worst_d = 0.0, None
    for d in range(-(m - 1), m):
        diag = np.diag(K, d)
        if diag.size < 2:
            continue
        sp = float(diag.max() - diag.min())
        sc = float(np.abs(diag).max())
        rel = sp / sc if sc > 0 else 0.0
        if rel > worst:
            worst, worst_d = rel, d
    return worst, worst_d


if __name__ == "__main__":
    # Einstein spectrum at the H3S omega_log used in P2 5.6; any lambda_pos works,
    # the structural questions do not depend on the spectrum's details.
    for n_pos, lam_val, Omega_mev in ((24, 1.72, 110.0), (40, 1.72, 110.0)):
        m = 2 * n_pos
        T_mev = 190.0 * sf.KB_MEV_PER_K
        nu = 2 * np.pi * T_mev * np.arange(m)
        lam_pos = lam_val * Omega_mev ** 2 / (nu ** 2 + Omega_mev ** 2)

        lam, C, D, E = sf.vertex_kernels(lam_pos, n_pos)
        KV = sf.kv_delta(lam, C, D)

        print(f"=== n_pos={n_pos}  (m={m})  Einstein Omega={Omega_mev} meV, "
              f"lambda={lam_val} ===")
        print(f"  Migdal kernel Lambda: Toeplitz defect = "
              f"{toeplitz_defect(lam)[0]:.2e}   symmetric = "
              f"{np.allclose(lam, lam.T)}   min elem = {lam.min():.4f}")

        td, dd = toeplitz_defect(KV)
        print(f"  H1 Toeplitz?   max rel spread within a diagonal = {td:.3f}"
              f" (worst offset d={dd})  -> {'TOEPLITZ' if td < 1e-10 else 'NOT Toeplitz'}")
        asym = np.abs(KV - KV.T).max() / np.abs(KV).max()
        print(f"  H2 symmetric?  ||K-K^T||_max/||K||_max = {asym:.2e}"
              f"  -> {'SYMMETRIC' if asym < 1e-10 else 'NOT symmetric'}")
        off = KV - np.diag(np.diag(KV))
        npos_frac = float((off > 0).sum()) / off.size
        print(f"  H3 sign-definite offdiagonal?  frac(>0) = {npos_frac:.3f},"
              f"  min = {off.min():.3f}, max = {off.max():.3f}"
              f"  -> {'SIGN-DEFINITE' if npos_frac in (0.0, 1.0) else 'INDEFINITE'}")
        row = KV[n_pos, :]
        d = np.abs(np.arange(m) - n_pos)
        order = np.argsort(d, kind="stable")
        rv = row[order]
        mono = bool(np.all(np.diff(np.abs(rv)) <= 1e-12))
        print(f"  H4 monotone in |m-k| along the central row? {mono}"
              f"   (row max {row.max():.3f}, min {row.min():.3f})")
        # is the vertex kernel even a perturbation of the same sign everywhere?
        print(f"  extra: is K^V_Delta <= 0 elementwise (a pure weakening)?  "
              f"{bool((KV <= 0).all())}   frac(<0) = {float((KV<0).sum())/KV.size:.3f}")
        print()

    print("=== claim (4) re-check, exact ===")
    import sympy as sp
    Dl, w = sp.symbols("Delta omega", positive=True)
    print("   d2/dD2 gamma^Delta =", sp.simplify(sp.diff(Dl / sp.sqrt(w**2 + Dl**2), Dl, 2)),
          "  -> < 0 for ALL omega>0, Delta>0  => gamma^Delta IS concave")
    print("   d2/dD2 gamma^omega =", sp.simplify(sp.diff(w / sp.sqrt(w**2 + Dl**2), Dl, 2)),
          "  -> sign flips at Delta = omega/sqrt(2)")
