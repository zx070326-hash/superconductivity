"""Robustness of the nesting-function local-maxima locations to the Lorentzian
broadening eta (a free numerical parameter, not physical -- must check the
qualitative peak-location claim isn't an eta artifact before reporting it)."""
import importlib.util, os
import numpy as np

_here = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    "graser1003_hamiltonian", os.path.join(_here, "_graser1003_hamiltonian_2026-09-22.py"))
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
build_H_graser = _mod.build_H_graser

N = 128
NBAND = 5
kx = np.linspace(-np.pi, np.pi, N, endpoint=False)
ky = kx.copy()
bands = np.empty((N, N, NBAND))
for i, kxi in enumerate(kx):
    for j, kyj in enumerate(ky):
        bands[i, j, :] = np.linalg.eigvalsh(build_H_graser(kxi, kyj, 0.0))
flat = bands.reshape(-1, NBAND)
lo, hi = flat.min() - 1, flat.max() + 1
for _ in range(60):
    mu = 0.5 * (lo + hi)
    n = 2 * np.mean(flat < mu) * NBAND
    lo, hi = (mu, hi) if n < 6.0 else (lo, mu)
mu = 0.5 * (lo + hi)
qx = np.fft.fftshift(np.fft.fftfreq(N, d=1.0/N)) * (2*np.pi/N)
qy = qx.copy()
ci, cj = N//2, N//2

for ETA in (0.01, 0.02, 0.03, 0.05, 0.08, 0.15):
    A = np.sum(ETA / np.pi / ((bands - mu) ** 2 + ETA ** 2), axis=2)
    F = np.fft.fft2(A)
    xi = np.fft.fftshift(np.fft.ifft2(F * np.conj(F)).real)
    avg = xi.mean()
    def xi_at(qx0, qy0):
        i = np.argmin(np.abs(qx - qx0)); j = np.argmin(np.abs(qy - qy0))
        return xi[i, j]
    xg = xi_at(0, 0); xx = xi_at(np.pi, 0); xm = xi_at(np.pi, np.pi)
    # search diagonal 0.3pi-0.7pi for the true diagonal peak
    best = (-1, 0, 0)
    for qd in np.linspace(0.3*np.pi, 0.7*np.pi, 41):
        v = xi_at(qd, qd)
        if v > best[0]:
            best = (v, qd, qd)
    print(f"eta={ETA:.2f}: xi(G)={xg:9.1f} xi(X)={xx:9.1f}({xx/avg:.3f}x) "
          f"xi(M)={xm:9.1f}({xm/avg:.3f}x) diag-peak q={best[1]/np.pi:.3f}pi "
          f"xi={best[0]:9.1f}({best[0]/avg:.3f}x)")
