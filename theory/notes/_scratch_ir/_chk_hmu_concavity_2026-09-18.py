"""
Independent from-scratch rebuild of the KAY kernel-generic operator O[f], the
Coulomb pencil M = O[f_P] - 2*mutilde*u_c u_c^T, and probes for:

  (1) self-checks  : O[1] = 2 u u^T - I ; h(10), h(20) vs repo table (B1)
  (2) reproduction : weak-coupling mu*>0 extremality counterexample
  (3) gap/cone     : is the top spectral gap collapsed at the violating point?
  (4) H_mu(x)      : concavity scan at fixed (N_c, mutilde) over the relevant x
  (5) envelope     : upper concave envelope of H_mu and the optimal 2-atom P

NOTHING is imported from the repository.  Operator matrix elements are written
directly from the definition quoted in theory/derivations/P4_conjecture2.md:190

   O[f] := -L1[f] + L2[f] + L3[f]
   L1_nn = (2/(2n+1)) * sum_{k=1..n} f(k)
   L2_nm = (1-delta_nm) f(|n-m|) / sqrt((2n+1)(2m+1))
   L3_nm =              f(n+m+1) / sqrt((2n+1)(2m+1))

   H(varpi) = O[ k -> varpi^2/(k^2+varpi^2) ]
   K_c      = 2 |u_c><u_c| ,  (u_c)_n = 1_{n<N_c} / sqrt(2n+1)
   M        = O[F_P] - mutilde * K_c ,   mutilde = mu*/lambda
   Tc condition: lambda * lambda_1(M) = 1

Rigor level of everything produced here: numerical-evidence (float64).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.fft import rfft, irfft
from scipy.sparse.linalg import LinearOperator, eigsh
from scipy.optimize import brentq

# ----------------------------------------------------------------- kernels
def f_lorentz(k, w2):
    """single Einstein mode at varpi^2 = w2 :  f(k) = w2/(k^2+w2)"""
    return w2 / (k * k + w2)

def f_mix(k, weights, w2s):
    out = np.zeros_like(k, dtype=float)
    for wt, w2 in zip(weights, w2s):
        out += wt * (w2 / (k * k + w2))
    return out

# ----------------------------------------------------------------- dense build
def O_dense(fvals, N):
    """fvals[j] = f(j), j = 0 .. 2N-1.  Returns dense N x N O[f]."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    I, J = np.meshgrid(n, n, indexing="ij")
    T = fvals[np.abs(I - J)]
    np.fill_diagonal(T, 0.0)
    K = fvals[I + J + 1]
    A = (T + K) * np.outer(u, u)
    # L1 diagonal: (2/(2n+1)) sum_{k=1..n} f(k)
    csum = np.concatenate(([0.0], np.cumsum(fvals[1:N])))  # csum[n] = sum_{k=1..n} f(k)
    A[n, n] -= 2.0 * u * u * csum
    return A

# ----------------------------------------------------------------- FFT matvec
class OpFFT(LinearOperator):
    """O[f] as a LinearOperator, O(N log N) matvec, plus optional rank-one
    Coulomb term  -2*mt*u_c u_c^T  (u_c supported on n < Nc)."""

    def __init__(self, fvals, N, mt=0.0, Nc=0):
        super().__init__(dtype=np.float64, shape=(N, N))
        self.N = N
        n = np.arange(N)
        self.u = 1.0 / np.sqrt(2 * n + 1.0)
        self.f0 = fvals[0]
        csum = np.concatenate(([0.0], np.cumsum(fvals[1:N])))
        self.diagL1 = 2.0 * self.u * self.u * csum
        L = 1
        while L < 2 * N + 2:
            L *= 2
        self.L = L
        # Toeplitz circulant first column: [f0, f1, ..., f_{N-1}, 0, f_{N-1}, ..., f1]
        col = np.zeros(L)
        col[:N] = fvals[:N]
        col[L - N + 1:] = fvals[1:N][::-1]
        self.Fc = rfft(col)
        # Hankel via convolution with g(j) = f(j+1), j = 0 .. 2N-2
        g = np.zeros(L)
        g[:2 * N - 1] = fvals[1:2 * N]
        self.Fg = rfft(g)
        self.mt = mt
        self.uc = np.zeros(N)
        if Nc > 0:
            m = min(Nc, N)
            self.uc[:m] = self.u[:m]

    def _matvec(self, x):
        x = np.asarray(x).ravel()
        N, L = self.N, self.L
        z = self.u * x
        pz = np.zeros(L); pz[:N] = z
        toep = irfft(self.Fc * rfft(pz), L)[:N] - self.f0 * z
        zr = np.zeros(L); zr[:N] = z[::-1]
        hank = irfft(self.Fg * rfft(zr), L)[N - 1:2 * N - 1]
        y = self.u * (toep + hank) - self.diagL1 * x
        if self.mt != 0.0:
            y -= (2.0 * self.mt * (self.uc @ x)) * self.uc
        return y


def top_eigs(op, k=1, v0=None, tol=0.0):
    vals, vecs = eigsh(op, k=k, which="LA", v0=v0, tol=tol, maxiter=20000)
    idx = np.argsort(vals)[::-1]
    return vals[idx], vecs[:, idx]


def lam1_mix(weights, w2s, varpi_ref, mt=0.0, Nc=0, Nfac=40, Nmin=400, k=1):
    """lambda_1 (and lambda_2 if k=2) of O[F_P] - 2 mt u_c u_c^T."""
    N = max(Nmin, int(Nfac * varpi_ref))
    kk = np.arange(2 * N, dtype=float)
    fv = f_mix(kk, weights, w2s)
    op = OpFFT(fv, N, mt=mt, Nc=Nc)
    vals, vecs = top_eigs(op, k=k)
    return vals, vecs, N


def h_mu(varpi, mt=0.0, Nc=0, Nfac=40, Nmin=400):
    v, _, _ = lam1_mix([1.0], [varpi ** 2], varpi, mt=mt, Nc=Nc, Nfac=Nfac, Nmin=Nmin)
    return v[0]


# ================================================================= (1) checks
def check_self():
    print("=" * 74)
    print("(1) SELF-CHECKS  (independent construction, nothing imported)")
    print("=" * 74)
    N = 60
    fv = np.ones(2 * N)
    A = O_dense(fv, N)
    n = np.arange(N); u = 1.0 / np.sqrt(2 * n + 1.0)
    ref = 2 * np.outer(u, u) - np.eye(N)
    print(f"  O[1] == 2uu^T - I :  max|diff| = {np.abs(A - ref).max():.3e}   (analytic identity)")

    # dense vs FFT
    N = 300; w2 = 7.0 ** 2
    kk = np.arange(2 * N, dtype=float)
    fv = f_lorentz(kk, w2)
    Ad = O_dense(fv, N)
    op = OpFFT(fv, N)
    X = np.random.default_rng(0).standard_normal((N, 4))
    err = max(np.abs(op @ X[:, j] - Ad @ X[:, j]).max() for j in range(4))
    print(f"  FFT matvec vs dense build : max|diff| = {err:.3e}")

    # repo table (B1) anchors, computed by my own code
    for varpi, ref_h in [(10.0, 2.873537558690), (20.0, 3.549054516265),
                         (50.0, 4.448928310050)]:
        # Richardson in N exactly as the repo describes (h - h_N ~ C (varpi/N)^4)
        hN = h_mu(varpi, Nfac=80, Nmin=800)
        h2N = h_mu(varpi, Nfac=160, Nmin=1600)
        hstar = h2N + (h2N - hN) / 15.0
        print(f"  h({varpi:5.0f}) = {hstar:.12f}   repo (B1) {ref_h:.12f}   "
              f"diff {hstar - ref_h:+.2e}")
    print()


# ============================================ (2) weak-coupling counterexample
def weak_coupling_case(lam, mu_bare, r, Nfac=20, verbose=True):
    """Solve Einstein Tc (i.e. varpi_E with lambda*h_mu = 1, N_c = floor(r*varpi)),
    then test two-atom spectra at the same T (same varpi_E-defined T, same N_c)."""
    mt = mu_bare / lam
    mustar_rms = mu_bare / (1.0 + mu_bare * np.log(r))

    def resid(lw):
        v = np.exp(lw)
        Nc = int(np.floor(r * v))
        return lam * h_mu(v, mt=mt, Nc=Nc, Nfac=Nfac, Nmin=400) - 1.0

    lo, hi = np.log(2.0), np.log(3.0e4)
    lwE = brentq(resid, lo, hi, xtol=1e-10, rtol=1e-12)
    vE = np.exp(lwE)
    Nc = int(np.floor(r * vE))
    if verbose:
        print(f"  lambda={lam:.3f}  mu(omega_c)={mu_bare:.3f}  mutilde={mt:.4f}  "
              f"r={r:g}  mu*_rms={mustar_rms:.4f}  2mu*/lambda={2*mustar_rms/lam:.3f}")
        print(f"  varpi_E = {vE:.4f}   N_c = {Nc}")
    return vE, Nc, mt, mustar_rms


def two_atom(vE, w, sp):
    """w at varpi_1, (1-w) at varpi_2 = sp*varpi_1, with w*x1+(1-w)*x2 = vE^2."""
    x1 = vE ** 2 / (w + (1 - w) * sp ** 2)
    x2 = x1 * sp ** 2
    return [w, 1 - w], [x1, x2]


def scan_counterexample(lam=0.50, mu_bare=0.30, r=5.0, Nfac=20):
    print("=" * 74)
    print("(2) WEAK-COUPLING COUNTEREXAMPLE  (independent recompute)")
    print("=" * 74)
    vE, Nc, mt, mustar = weak_coupling_case(lam, mu_bare, r, Nfac=Nfac)
    qE, _, Nused = lam1_mix([1.0], [vE ** 2], vE, mt=mt, Nc=Nc, Nfac=Nfac)
    print(f"  N used = {Nused};  lambda*Q(E)-1 = {lam*qE[0]-1:+.3e}  (brentq residual)")
    print(f"  {'w':>6} {'sp':>9} {'varpi_1':>10} {'varpi_2':>11} "
          f"{'lam*Q(P)-1':>13}  verdict")
    best = None
    for w in (0.90, 0.95, 0.97, 0.99, 0.995):
        for sp in (3.0, 10.0, 30.0, 100.0, 300.0):
            ws, x2s = two_atom(vE, w, sp)
            vref = max(np.sqrt(x2s[1]), vE)
            q, _, _ = lam1_mix(ws, x2s, vref, mt=mt, Nc=Nc, Nfac=Nfac)
            rel = lam * q[0] - 1.0
            flag = "VIOLATION" if rel > 0 else "ok"
            print(f"  {w:6.3f} {sp:9.1f} {np.sqrt(x2s[0]):10.3f} "
                  f"{np.sqrt(x2s[1]):11.2f} {rel:+13.3e}  {flag}")
            if best is None or rel > best[0]:
                best = (rel, w, sp, ws, x2s)
    print(f"\n  best: rel = {best[0]:+.4e} at w={best[1]}, sp={best[2]}")
    return vE, Nc, mt, best


# =============================================== (3) gap + cone at violation
def gap_and_cone(vE, Nc, mt, ws, x2s, Nfac=20):
    print("=" * 74)
    print("(3) IS THE TOP SPECTRAL GAP COLLAPSED AT THE VIOLATING POINT?")
    print("=" * 74)
    for label, wts, xs, vref in [("Einstein", [1.0], [vE ** 2], vE),
                                 ("two-atom", ws, x2s, max(np.sqrt(x2s[1]), vE))]:
        vals, vecs, N = lam1_mix(wts, xs, vref, mt=mt, Nc=Nc, Nfac=Nfac, k=2)
        v = vecs[:, 0]
        if v[np.argmax(np.abs(v))] < 0:
            v = -v
        nneg = int(np.sum(v < 0))
        # cone C = {v >= 0, nonincreasing}
        noninc = bool(np.all(np.diff(v) <= 1e-14))
        print(f"  {label}:  lam1={vals[0]:.6f}  lam2={vals[1]:.6f}  "
              f"gap={vals[0]-vals[1]:.6f}  gap/lam1={(vals[0]-vals[1])/vals[0]:.4f}")
        print(f"            N={N}  #negative components of v = {nneg}  "
              f"nonincreasing = {noninc}")
    # compare with mu*=0
    vals0, _, _ = lam1_mix([1.0], [vE ** 2], vE, mt=0.0, Nc=0, Nfac=Nfac, k=2)
    print(f"  Einstein, mu*=0 : lam1={vals0[0]:.6f} lam2={vals0[1]:.6f} "
          f"gap={vals0[0]-vals0[1]:.6f}")
    print()


# ================================================== (4) H_mu concavity scan
def concavity_scan(vE, Nc, mt, npts=25, span=(0.02, 60.0), Nfac=20):
    """H_mu(x) = lambda_1( H(sqrt(x)) - 2 mt u_c u_c^T ), N_c FIXED.
    Scan x = (s*vE)^2 on a log grid; report discrete second derivative in x."""
    print("=" * 74)
    print("(4) CONCAVITY OF H_mu(x) AT FIXED N_c   (x = varpi^2)")
    print("=" * 74)
    s = np.exp(np.linspace(np.log(span[0]), np.log(span[1]), npts))
    vs = s * vE
    H = np.array([h_mu(v, mt=mt, Nc=Nc, Nfac=Nfac, Nmin=400) for v in vs])
    x = vs ** 2
    print(f"  N_c = {Nc} (fixed),  mutilde = {mt:.4f},  x_bar = varpi_E^2 = {vE**2:.4g}")
    print(f"  {'varpi':>11} {'x':>13} {'H_mu':>12} {'2nd-diff sign':>15}")
    d2 = np.full(npts, np.nan)
    for i in range(1, npts - 1):
        # divided second difference on a nonuniform grid
        h1, h2 = x[i] - x[i - 1], x[i + 1] - x[i]
        d2[i] = 2 * ((H[i + 1] - H[i]) / h2 - (H[i] - H[i - 1]) / h1) / (h1 + h2)
    for i in range(npts):
        tag = "" if np.isnan(d2[i]) else ("CONVEX(+)" if d2[i] > 0 else "concave(-)")
        print(f"  {vs[i]:11.4f} {x[i]:13.5g} {H[i]:12.6f} {tag:>15}")
    bad = [(vs[i], d2[i]) for i in range(npts) if not np.isnan(d2[i]) and d2[i] > 0]
    print(f"\n  grid points with H_mu'' > 0 : {len(bad)} of {npts-2}")
    if bad:
        print(f"    varpi range of convexity: [{min(b[0] for b in bad):.3f}, "
              f"{max(b[0] for b in bad):.3f}]")
    return vs, x, H


# ============================================ (5) concave envelope of H_mu
def concave_envelope(x, H, xbar):
    """Upper concave envelope of the sampled (x,H) on the log-spaced grid,
    evaluated at xbar; returns (value, contact points)."""
    # upper convex hull of the point set
    pts = sorted(zip(x, H))
    hull = []
    for p in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            # remove hull[-1] if it lies below the chord hull[-2]->p
            if (y2 - y1) * (p[0] - x1) <= (p[1] - y1) * (x2 - x1):
                hull.pop()
            else:
                break
        hull.append(p)
    hx = np.array([p[0] for p in hull]); hy = np.array([p[1] for p in hull])
    val = np.interp(xbar, hx, hy)
    j = int(np.searchsorted(hx, xbar))
    j = max(1, min(j, len(hx) - 1))
    return val, (hx[j - 1], hx[j]), (hy[j - 1], hy[j])


def envelope_report(x, H, vE, lam):
    print("=" * 74)
    print("(5) UPPER CONCAVE ENVELOPE of H_mu  -> optimal two-atom / Tc deficit")
    print("=" * 74)
    xbar = vE ** 2
    HE = float(np.interp(xbar, x, H))
    val, (xa, xb), (ya, yb) = concave_envelope(x, H, xbar)
    print(f"  H_mu(x_bar)          = {HE:.6f}   (Einstein)")
    print(f"  envelope Hhat(x_bar) = {val:.6f}")
    print(f"  excess               = {val - HE:+.3e}   "
          f"(>0  <=>  Einstein NOT extremal)")
    if xb > xa:
        wstar = (xb - xbar) / (xb - xa)
        print(f"  supporting chord touches varpi = {np.sqrt(xa):.4f} and "
              f"{np.sqrt(xb):.4f};  optimal weight w = {wstar:.4f}")
    print(f"  lambda*Hhat - 1      = {lam*val - 1:+.3e}")
    print()


if __name__ == "__main__":
    check_self()
    vE, Nc, mt, best = scan_counterexample()
    gap_and_cone(vE, Nc, mt, best[3], best[4])
    vs, x, H = concavity_scan(vE, Nc, mt)
    envelope_report(x, H, vE, 0.50)
