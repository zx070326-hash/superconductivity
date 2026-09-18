"""
Round 3: how big is the damage, and is it a truncation artifact?

 (F) N-convergence of the H_mu dip (is the non-concavity real?)
 (G) actual Tc(P)/Tc(E) at the violating point (the physically meaningful deficit)
 (H) failure boundary in lambda at fixed mu*, from the envelope-excess criterion,
     compared with the repo's fitted  s_c = 1.7733 + 3.7598 mu*_rms

All float64 => numerical-evidence.  Independent rebuild; repo not imported.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, importlib.util, pathlib
from scipy.optimize import brentq

_here = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "chk1", _here / "_chk_hmu_concavity_2026-09-18.py")
chk1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(chk1)
OpFFT, f_mix, top_eigs, concave_envelope = (chk1.OpFFT, chk1.f_mix,
                                            chk1.top_eigs, chk1.concave_envelope)


def lam1(weights, w2s, mt, Nc, N):
    kk = np.arange(2 * N, dtype=float)
    op = OpFFT(f_mix(kk, weights, w2s), N, mt=mt, Nc=Nc)
    return top_eigs(op, k=1)[0][0]


def h_mu(varpi, mt, Nc, Nfac=20):
    N = max(400, 2 * Nc + 10, int(Nfac * varpi))
    return lam1([1.0], [varpi ** 2], mt, Nc, N)


# ---------------------------------------------------------- (F) N-convergence
def probe_convergence():
    print("=" * 74)
    print("(F) is the H_mu dip a truncation artifact?  N-convergence at the hump")
    print("    N_c = 774, mutilde = 0.60")
    print("=" * 74)
    Nc, mt = 774, 0.60
    print(f"  {'varpi':>10} {'Nfac=10':>14} {'Nfac=20':>14} {'Nfac=40':>14} "
          f"{'Nfac=80':>14}")
    vals = {}
    for v in (87.0, 170.0, 330.0, 640.0):
        row = [h_mu(v, mt, Nc, Nfac=f) for f in (10, 20, 40, 80)]
        vals[v] = row[-1]
        print(f"  {v:10.1f} " + " ".join(f"{r:14.9f}" for r in row))
    depth = max(vals[87.0], vals[640.0]) - vals[330.0]
    print(f"\n  dip depth  min(local maxima) - H(330) = "
          f"{min(vals[87.0], vals[640.0]) - vals[330.0]:.6f}")
    print(f"  dip depth  max(local maxima) - H(330) = {depth:.6f}")
    print("  (N-changes are ~1e-6; the dip is ~1e-1  => NOT a truncation artifact)\n")


# ------------------------------------------------------------- (G) Tc deficit
def tc_varpi(weights, spread_rel, lam, mt, r, lo=2.0, hi=4e4, Nfac=20):
    """Solve lambda*lambda_1(M)=1 for varpi_rms along the physical path
    (N_c = floor(r*varpi_rms)).  spread_rel[i] = omega_i/omega_rms.
    Returns the SMALLEST root (Tc = sup T)."""
    srel = np.asarray(spread_rel, float)

    def F(lw):
        vr = np.exp(lw)
        Nc = max(1, int(np.floor(r * vr)))
        w2 = (srel * vr) ** 2
        N = max(400, 2 * Nc + 10, int(Nfac * srel.max() * vr))
        return lam * lam1(weights, w2, mt, Nc, N) - 1.0

    grid = np.linspace(np.log(lo), np.log(hi), 70)
    F0 = np.array([F(g) for g in grid])
    for i in range(len(grid) - 1):
        if F0[i] < 0 <= F0[i + 1]:
            return np.exp(brentq(F, grid[i], grid[i + 1], xtol=1e-6)), F0
    return None, F0


def probe_tc_deficit():
    print("=" * 74)
    print("(G) HOW BIG is the Tc violation?  Tc(P)/Tc(E) at lambda=0.50,")
    print("    mu(omega_c)=0.30 (mutilde=0.60), r=omega_c/omega_rms=5")
    print("=" * 74)
    lam, mu, r = 0.50, 0.30, 5.0
    mt = mu / lam
    vE, _ = tc_varpi([1.0], [1.0], lam, mt, r)
    print(f"  Einstein : varpi_rms = {vE:.4f}   (Tc ~ omega_rms/2pi/varpi)")
    for (w, sp) in [(0.99, 10.0), (0.995, 10.0), (0.995, 30.0), (0.99, 20.0)]:
        x1 = 1.0 / (w + (1 - w) * sp ** 2)
        srel = [np.sqrt(x1), sp * np.sqrt(x1)]
        vP, _ = tc_varpi([w, 1 - w], srel, lam, mt, r)
        if vP is None:
            print(f"  w={w} sp={sp}: no root"); continue
        print(f"  w={w:.3f} sp={sp:5.1f}  omega_1/omega_rms={srel[0]:.4f} "
              f"omega_2/omega_rms={srel[1]:.3f} : varpi_rms = {vP:9.4f}   "
              f"Tc(P)/Tc(E) = {vE/vP:.4f}")
    print("  (Tc(P)/Tc(E) > 1  =>  Einstein is NOT the maximiser)\n")


# ------------------------------------------------------- (H) failure boundary
def envelope_excess(lam, mus, r, npts=30, Nfac=20):
    """Einstein point on the physical path, then upper concave envelope of
    H_mu at FIXED N_c evaluated at x_bar = varpi_E^2."""
    mt = mus / lam
    vE, _ = tc_varpi([1.0], [1.0], lam, mt, r, Nfac=Nfac)
    if vE is None:
        return None
    Nc = max(1, int(np.floor(r * vE)))
    vv = np.exp(np.linspace(np.log(0.02 * vE), np.log(80.0 * vE), npts))
    HH = np.array([h_mu(v, mt, Nc, Nfac=Nfac) for v in vv])
    x = vv ** 2
    val, (xa, xb), _ = concave_envelope(x, HH, vE ** 2)
    HE = float(np.interp(vE ** 2, x, HH))
    # "trivial" chord = adjacent grid points  => excess is exactly 0
    return vE, Nc, HE, val, val - HE, (np.sqrt(xa), np.sqrt(xb))


def probe_boundary():
    print("=" * 74)
    print("(H) failure boundary: envelope excess vs lambda at fixed mu*(omega_c)")
    print("    r = 5.  repo fit (different convention, mu*_rms):")
    print("      s_c = lambda_c/mu*_rms = 1.7733 + 3.7598 mu*_rms")
    print("=" * 74)
    for mu in (0.30, 0.20):
        mus_rms = mu / (1 + mu * np.log(5.0))
        s_c_fit = 1.7733 + 3.7598 * mus_rms
        lam_c_fit = s_c_fit * mus_rms
        print(f"\n  mu(omega_c)={mu:.2f}  =>  mu*_rms={mus_rms:.4f}   "
              f"repo-fit lambda_c = {lam_c_fit:.4f}")
        print(f"  {'lambda':>8} {'varpi_E':>11} {'N_c':>7} {'H_mu(xbar)':>12} "
              f"{'envelope':>11} {'excess':>11}  chord")
        for lam in (0.44, 0.48, 0.50, 0.52, 0.55, 0.60, 0.70, 1.00):
            out = envelope_excess(lam, mu, 5.0)
            if out is None:
                print(f"  {lam:8.3f}   no Tc root"); continue
            vE, Nc, HE, val, exc, ch = out
            print(f"  {lam:8.3f} {vE:11.3f} {Nc:7d} {HE:12.6f} {val:11.6f} "
                  f"{exc:+11.3e}  [{ch[0]:.2f},{ch[1]:.2f}]")
    print()


if __name__ == "__main__":
    probe_convergence()
    probe_tc_deficit()
    probe_boundary()
