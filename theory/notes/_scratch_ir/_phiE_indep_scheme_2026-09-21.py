"""INDEPENDENT SECOND SCHEME for T_c of the Einstein spectrum with mu* > 0.

Deliberately shares NOTHING with the repo pencil machinery except the physics:
  - full integer-Z Matsubara index set n = -N..N-1 (NO folding to n>=0),
  - Z(n) computed by BRUTE-FORCE truncated sum (NO analytic
    "lam(0) + 2 sum_{j=1..n} lam(j)" cancellation -- that cancellation is exactly
    the algebra that H_mat encodes, so re-using it would not be an independent test),
  - NON-symmetric operator  B[n,m] = (lam(n-m) - mu* th_n th_m) / (Z_n |2m+1|),
    top eigenvalue by ARPACK / power iteration (NOT eigvalsh of a symmetric matrix),
  - lam NOT factored out of the pencil; T_c found by bisecting wbar for rho(wbar)=1,
  - Toeplitz matvec via FFT (O(N log N)), so N can be pushed to test truncation.

Model (matched to the repo convention so the comparison is apples-to-apples):
    lam(j)  = lam * wbar^2/(wbar^2 + j^2),      wbar = omega_E/(2 pi T)
    th_n    = 1 for -Nc <= n <= Nc-1,           Nc = floor(r*wbar), r = omega_c/omega_E
    Z_n     = 1 + (1/(2n+1)) sum_m lam(n-m) sgn(2m+1)
    Z_n D_n = sum_m [lam(n-m) - mu* th_n th_m] D_m / |2m+1|
    T_c: rho(wbar) = 1;  T_c/omega_E = 1/(2 pi wbar).

An ASYMMETRIC-Coulomb variant (th_m only, the literal Marsiglio eq.23 /
Marsiglio-Carbotte eq.56 form) is also provided to quantify the symmetrization cost
that P4D_coulomb_operator.md sec.1 declares benign.

LEVEL: numerical-evidence (float64, finite N; N-convergence reported).
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

TWO_PI = 2.0 * np.pi


def _toeplitz_matvec(lam_full, x, N):
    """(L x)_n = sum_{m=0}^{2N-1} lam_full[n-m+(2N-1)] x_m, n=0..2N-1, via FFT.

    lam_full has length 4N-1 and index offset 2N-1 (i.e. lam_full[j+2N-1] = lam(j)).
    """
    M = 2 * N
    L = 1
    while L < (len(lam_full) + M):
        L *= 2
    F = np.fft.rfft(lam_full, L) * np.fft.rfft(x, L)
    conv = np.fft.irfft(F, L)
    # conv[k] = sum_m lam_full[k-m] x[m]; want k = n + (2N-1)
    return conv[(M - 1):(M - 1) + M]


class IndepEliashberg:
    def __init__(self, lam, mu_star, wbar, r, N=2000, symmetric_coulomb=True):
        self.lam, self.mu, self.w, self.r, self.N = lam, mu_star, wbar, r, N
        self.sym = symmetric_coulomb
        M = 2 * N
        self.n = np.arange(-N, N)              # n = -N .. N-1
        self.d = np.abs(2.0 * self.n + 1.0)    # |2n+1|
        j = np.arange(-(M - 1), M)             # -(2N-1) .. 2N-1
        self.lam_full = lam * wbar**2 / (wbar**2 + j.astype(float)**2)
        self.Nc = int(np.floor(r * wbar + 1e-9))
        th = np.zeros(M)
        th[(self.n >= -self.Nc) & (self.n <= self.Nc - 1)] = 1.0
        self.th = th
        # Z by brute-force truncated sum (independent of any analytic cancellation)
        s = np.sign(2.0 * self.n + 1.0)
        self.Z = 1.0 + _toeplitz_matvec(self.lam_full, s, N) / (2.0 * self.n + 1.0)

    def matvec(self, x):
        y = x / self.d
        out = _toeplitz_matvec(self.lam_full, y, self.N)
        if self.mu != 0.0:
            c = float(np.dot(self.th, y))
            if self.sym:
                out = out - self.mu * self.th * c
            else:
                out = out - self.mu * c      # th on m only -> constant in n
        return out / self.Z

    def rho(self, tol=1e-13, itmax=20000):
        """Largest eigenvalue by power iteration with Rayleigh-quotient readout."""
        x = np.exp(-np.abs(2.0 * self.n + 1.0) / (4.0 * max(self.w, 0.3) + 4.0))
        x /= np.linalg.norm(x)
        lamprev = 0.0
        for it in range(itmax):
            y = self.matvec(x)
            nrm = np.linalg.norm(y)
            if nrm == 0.0:
                return 0.0
            x = y / nrm
            if it % 5 == 4:
                r_ = float(np.dot(x, self.matvec(x)))
                if abs(r_ - lamprev) < tol * max(1.0, abs(r_)):
                    return r_
                lamprev = r_
        return float(np.dot(x, self.matvec(x)))


def rho_of(lam, mu_star, wbar, r, N=2000, sym=True):
    return IndepEliashberg(lam, mu_star, wbar, r, N, sym).rho()


def tc_indep(lam, mu_star, r, N=2000, w_max=40.0, sub=10, sym=True, eps=1e-9):
    """T_c/omega_E via the SUP rule, window-aware, in the independent scheme.

    Within a window Nc is frozen: we bisect on the CONTINUOUS branch
    rho(wbar; Nc fixed) = 1.
    """
    def rho_fixed(w, Nc):
        obj = IndepEliashberg(lam, mu_star, w, r, N, sym)
        obj.Nc = Nc
        th = np.zeros(2 * N)
        th[(obj.n >= -Nc) & (obj.n <= Nc - 1)] = 1.0
        obj.th = th
        return obj.rho()

    kmax = int(np.ceil(r * w_max))
    wfloor = 0.02
    for k in range(0, kmax + 1):
        a = max(k / r, wfloor)
        b = min((k + 1) / r, w_max)
        if b <= a:
            continue
        ws = np.linspace(a + eps, b - eps, sub)
        gs = [rho_fixed(w, k) for w in ws]
        if gs[0] >= 1.0:
            return 1.0 / (TWO_PI * ws[0]), float(ws[0]), k
        for i in range(1, sub):
            if gs[i - 1] < 1.0 <= gs[i]:
                lo, hi = ws[i - 1], ws[i]
                for _ in range(45):
                    m = 0.5 * (lo + hi)
                    if rho_fixed(m, k) < 1.0:
                        lo = m
                    else:
                        hi = m
                wstar = 0.5 * (lo + hi)
                return 1.0 / (TWO_PI * wstar), float(wstar), k
    return None, None, None


if __name__ == "__main__":
    print("N-convergence of the independent scheme (rho at a fixed point):")
    for N in (500, 1000, 2000, 4000, 8000):
        print(f"  N={N:5d}  rho(lam=2,mu*=0.10,wbar=0.9,r=3) = "
              f"{rho_of(2.0, 0.10, 0.9, 3.0, N):.12f}")
    print()
    print("mu*=0 anchor vs the certified table:")
    for lam, lo, hi in ((1.0, 0.1146122, 0.1147841), (2.0, 0.2115556, 0.2118729)):
        t, w, _ = tc_indep(lam, 0.0, 1.0, N=2000)
        print(f"  lam={lam}  Tc/wE={t:.8f}   cert[{lo:.7f},{hi:.7f}]  in={lo<=t<=hi}")
