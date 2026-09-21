"""
Round 2 probes (all float64 => numerical-evidence):

 (A) exact Jensen reduction survives mu*>0 :  lambda_1(M(P)) <= int h_mu dP
 (B) exact 2x2 secular equation for h_mu
 (C) where does the non-concavity of H_mu sit?  scaling with N_c
 (D) top-eigenvector sign pattern (cone exit) location vs N_c
 (E) strong-coupling / physical-hydride check: is H_mu concave there?

Independent rebuild; imports nothing from the repository.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.fft import rfft, irfft
from scipy.sparse.linalg import LinearOperator, eigsh
import importlib.util, pathlib

_here = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "chk1", _here / "_chk_hmu_concavity_2026-09-18.py")
chk1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(chk1)          # reuse MY OWN round-1 builder only

OpFFT, O_dense, f_mix, top_eigs = chk1.OpFFT, chk1.O_dense, chk1.f_mix, chk1.top_eigs


def lam1(weights, w2s, mt=0.0, Nc=0, N=None, vref=None, k=1):
    if N is None:
        N = max(400, 2 * Nc + 10, int(20 * vref))
    kk = np.arange(2 * N, dtype=float)
    fv = f_mix(kk, weights, w2s)
    op = OpFFT(fv, N, mt=mt, Nc=Nc)
    vals, vecs = top_eigs(op, k=k)
    return vals, vecs, N


def h_mu(varpi, mt, Nc, N=None):
    v, _, _ = lam1([1.0], [varpi ** 2], mt=mt, Nc=Nc, N=N, vref=varpi)
    return v[0]


# ------------------------------------------------------- (A) Jensen reduction
def probe_jensen():
    print("=" * 74)
    print("(A) does the layer-one Jensen step survive mu*>0 ?")
    print("    claim [analytic]: M(P) = int [ H(varpi_w) - 2 mt u_c u_c^T ] P(dw)")
    print("    because K_c is P-independent and P is a PROBABILITY measure")
    print("    => lambda_1(M(P)) <= int h_mu(varpi) P(dw)   (convexity of lambda_1)")
    print("=" * 74)
    rng = np.random.default_rng(7)
    Nc, mt, N = 40, 0.30, 1200
    print(f"  N_c={Nc}  mutilde={mt}  N={N}")
    print(f"  {'#atoms':>7} {'lam1(M(P))':>13} {'int h_mu dP':>13} {'slack':>12}  ok?")
    worst = -1e9
    for trial in range(12):
        na = rng.integers(2, 5)
        w = rng.dirichlet(np.ones(na))
        vs = np.exp(rng.uniform(np.log(0.3), np.log(60.0), na))
        v1, _, _ = lam1(w, vs ** 2, mt=mt, Nc=Nc, N=N)
        rhs = sum(wi * h_mu(vi, mt, Nc, N=N) for wi, vi in zip(w, vs))
        slack = rhs - v1[0]
        worst = max(worst, -slack)
        print(f"  {na:7d} {v1[0]:13.8f} {rhs:13.8f} {slack:12.3e}  "
              f"{'OK' if slack >= -1e-9 else 'VIOLATED'}")
    print(f"  worst violation of the inequality: {worst:.3e}  "
          f"(negative => inequality holds)\n")


# ------------------------------------------------------- (B) secular equation
def probe_secular():
    print("=" * 74)
    print("(B) exact 2x2 secular equation  [analytic, verified numerically]")
    print("    H(x) = 2 u u^T - I - W(x),  B := I + W(x) (so H = 2uu^T - B)")
    print("    M    = -B + [u,u_c] diag(2,-2mt) [u,u_c]^T")
    print("    alpha=<u,(h+B)^-1 u>, beta=<u,(h+B)^-1 u_c>, gamma=<u_c,(h+B)^-1 u_c>")
    print("    h solves   2 alpha(h) = 1 + 4 mt beta(h)^2 / (1 + 2 mt gamma(h))")
    print("=" * 74)
    N = 500
    for (varpi, Nc, mt) in [(5.0, 10, 0.20), (12.0, 30, 0.50), (3.0, 60, 0.13)]:
        kk = np.arange(2 * N, dtype=float)
        fv = f_mix(kk, [1.0], [varpi ** 2])
        H = O_dense(fv, N)
        n = np.arange(N); u = 1.0 / np.sqrt(2 * n + 1.0)
        uc = np.zeros(N); uc[:min(Nc, N)] = u[:min(Nc, N)]
        B = 2 * np.outer(u, u) - H          # = I + W(x)
        M = H - 2 * mt * np.outer(uc, uc)
        h = np.linalg.eigvalsh(M)[-1]
        R = np.linalg.inv(h * np.eye(N) + B)
        a, b, g = u @ R @ u, u @ R @ uc, uc @ R @ uc
        lhs, rhs = 2 * a, 1 + 4 * mt * b ** 2 / (1 + 2 * mt * g)
        print(f"  varpi={varpi:5.1f} N_c={Nc:3d} mt={mt:.2f}: h={h:.10f}  "
              f"LHS-RHS = {lhs - rhs:+.3e}")
    print("  (mt=0 reduces to the repo's rank-one secular 2<u,(h+B)^-1 u> = 1)\n")


# --------------------------------------------- (C) location of non-concavity
def probe_location(mt=0.60):
    print("=" * 74)
    print("(C) WHERE is H_mu non-concave?  scaling of the hump with N_c")
    print(f"    mutilde = {mt}")
    print("=" * 74)
    print(f"  {'N_c':>6} {'varpi@localmax':>15} {'H_max':>9} {'varpi@localmin':>15} "
          f"{'H_min':>9} {'vmax/N_c':>9} {'vmin/N_c':>9} {'1/H_max':>8}")
    for Nc in (50, 100, 200, 400, 774, 1500):
        vs = np.exp(np.linspace(np.log(0.01 * Nc), np.log(20.0 * Nc), 42))
        H = np.array([h_mu(v, mt, Nc) for v in vs])
        # local max = first interior maximum; local min = the following minimum
        imax = None
        for i in range(1, len(vs) - 1):
            if H[i] > H[i - 1] and H[i] > H[i + 1]:
                imax = i; break
        imin = None
        if imax is not None:
            for i in range(imax + 1, len(vs) - 1):
                if H[i] < H[i - 1] and H[i] < H[i + 1]:
                    imin = i; break
        if imax is None or imin is None:
            print(f"  {Nc:6d}   no interior hump found on grid")
            continue
        print(f"  {Nc:6d} {vs[imax]:15.3f} {H[imax]:9.4f} {vs[imin]:15.3f} "
              f"{H[imin]:9.4f} {vs[imax]/Nc:9.4f} {vs[imin]/Nc:9.4f} "
              f"{1/H[imax]:8.4f}")
    print("  => hump sits at a fixed FRACTION of N_c: the Tolmachev crossover")
    print("     omega ~ omega_c, a scale set by the COULOMB cutoff, not the phonons\n")


# ------------------------------------------------- (D) eigenvector sign pattern
def probe_cone():
    print("=" * 74)
    print("(D) top eigenvector: where does it leave the nonnegative cone?")
    print("=" * 74)
    for (varpi, Nc, mt) in [(154.8, 774, 0.60), (50.0, 250, 0.60), (12.0, 60, 0.30)]:
        vals, vecs, N = lam1([1.0], [varpi ** 2], mt=mt, Nc=Nc, vref=varpi, k=2)
        v = vecs[:, 0]
        if v[np.argmax(np.abs(v))] < 0:
            v = -v
        sgn = np.sign(v)
        flips = np.where(np.diff(sgn) != 0)[0]
        nneg = int((v < 0).sum())
        # mass-weighted: where is the first sign change carrying real weight?
        first = flips[0] + 1 if len(flips) else None
        print(f"  varpi={varpi:7.1f} N_c={Nc:5d} mt={mt:.2f} N={N}:  lam1={vals[0]:.5f} "
              f"lam2={vals[1]:.5f}  gap={vals[0]-vals[1]:.5f}")
        print(f"      #sign changes={len(flips)}  first flip at n={first}  "
              f"#neg={nneg}   n/N_c of first flip = "
              f"{first/Nc if first else float('nan'):.3f}")
        # mu*=0 control
        vals0, vecs0, _ = lam1([1.0], [varpi ** 2], mt=0.0, Nc=0, N=N, k=2)
        v0 = vecs0[:, 0]
        if v0[np.argmax(np.abs(v0))] < 0:
            v0 = -v0
        print(f"      mu*=0 control: #neg={int((v0<0).sum())}  gap="
              f"{vals0[0]-vals0[1]:.5f}")
    print()


# --------------------------------------- (E) strong coupling / physical values
def probe_strong():
    print("=" * 74)
    print("(E) PHYSICAL / STRONG-COUPLING SIDE: is H_mu concave there?")
    print("    hydride-like: lambda=2, mu*=0.13 => mutilde=0.065;  r = omega_c/omega = 10")
    print("=" * 74)
    for (lam, mus, r) in [(2.0, 0.13, 10.0), (1.0, 0.13, 10.0), (0.7, 0.13, 10.0),
                          (2.0, 0.20, 5.0), (0.35, 0.13, 10.0)]:
        mt = mus / lam
        # Einstein Tc on the physical path (N_c = floor(r*varpi)): find smallest root
        vs = np.exp(np.linspace(np.log(0.05), np.log(4000.0), 90))
        Hp = np.array([h_mu(v, mt, max(1, int(np.floor(r * v)))) for v in vs])
        idx = np.where(np.diff(np.sign(lam * Hp - 1.0)) > 0)[0]
        if len(idx) == 0:
            print(f"  lam={lam} mu*={mus} r={r}: no Tc root on grid"); continue
        i = idx[0]
        vE = float(np.exp(np.interp(1.0 / lam, [Hp[i], Hp[i + 1]],
                                    [np.log(vs[i]), np.log(vs[i + 1])])))
        Nc = max(1, int(np.floor(r * vE)))
        # now the Jensen object at FIXED N_c
        vv = np.exp(np.linspace(np.log(0.02 * vE), np.log(60.0 * vE), 34))
        HH = np.array([h_mu(v, mt, Nc) for v in vv])
        x = vv ** 2
        d2 = []
        for j in range(1, len(vv) - 1):
            h1, h2 = x[j] - x[j - 1], x[j + 1] - x[j]
            d2.append(2 * ((HH[j + 1] - HH[j]) / h2 - (HH[j] - HH[j - 1]) / h1)
                      / (h1 + h2))
        d2 = np.array(d2)
        bad = np.where(d2 > 0)[0]
        val, (xa, xb), _ = chk1.concave_envelope(x, HH, vE ** 2)
        HE = float(np.interp(vE ** 2, x, HH))
        print(f"  lam={lam:4.2f} mu*={mus:.2f} r={r:4.1f} mt={mt:.4f}: "
              f"varpi_E={vE:9.3f} N_c={Nc:6d}")
        print(f"      #convex grid pts = {len(bad)}/{len(d2)}"
              + (f"  varpi in [{vv[bad[0]+1]:.3f},{vv[bad[-1]+1]:.3f}]"
                 f" = [{vv[bad[0]+1]/Nc:.3f},{vv[bad[-1]+1]/Nc:.3f}]*N_c" if len(bad) else ""))
        print(f"      Einstein H_mu={HE:.6f}  envelope={val:.6f}  "
              f"excess={val-HE:+.3e}  (chord {np.sqrt(xa):.3f}..{np.sqrt(xb):.3f})")
    print()


if __name__ == "__main__":
    probe_jensen()
    probe_secular()
    probe_location()
    probe_cone()
    probe_strong()
