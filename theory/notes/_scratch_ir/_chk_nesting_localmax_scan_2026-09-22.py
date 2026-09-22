"""Is the reported 'global argmax at q=(-0.0625pi,0)' a genuine second nesting
feature, or just the smooth shoulder/tail of the broad q=0 self-correlation peak
(xi(0)=85442.57, ~sum of DOS^2-like term that is large and smooth near q=0, not a
sharp nesting signature)? The original script only masked a radius-3 disk (~0.047pi)
around Gamma before taking argmax -- too small if the self-peak has a broad smooth
tail extending well past that.

Method: (1) rerun with progressively larger exclusion radii around Gamma and see if
the argmax point keeps retreating outward (=> it IS the Gamma tail, no real
competing peak) or stabilizes at a fixed q (=> genuine separate local maximum).
(2) independently, find ALL local maxima of xi(q) on the grid (points exceeding all
8 nearest neighbors), away from Gamma, and report the strongest ones with location.
"""
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
ETA = 0.03

kx = np.linspace(-np.pi, np.pi, N, endpoint=False)
ky = np.linspace(-np.pi, np.pi, N, endpoint=False)
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

A = np.sum(ETA / np.pi / ((bands - mu) ** 2 + ETA ** 2), axis=2)
F = np.fft.fft2(A)
xi = np.fft.ifft2(F * np.conj(F)).real
xi = np.fft.fftshift(xi)
qx = np.fft.fftshift(np.fft.fftfreq(N, d=1.0/N)) * (2*np.pi/N)
qy = qx.copy()
ci, cj = N//2, N//2

print(f"mu = {mu:.4f} eV, xi(Gamma)={xi[ci,cj]:.2f}, BZ avg={xi.mean():.2f}")
print()
print("(1) argmax(excluding q=0) vs exclusion-disk radius (grid points):")
for rmask in (1, 2, 3, 5, 8, 12, 16, 20, 24, 28, 32):
    m = xi.copy()
    for di in range(-rmask, rmask+1):
        for dj in range(-rmask, rmask+1):
            if di*di+dj*dj <= rmask*rmask:
                m[(ci+di) % N, (cj+dj) % N] = -np.inf
    im, jm = np.unravel_index(np.argmax(m), m.shape)
    print(f"  rmask={rmask:3d} ({rmask*2*np.pi/N/np.pi:.3f}pi radius): "
          f"argmax q=({qx[im]/np.pi:+.4f},{qy[jm]/np.pi:+.4f})pi  xi={xi[im,jm]:10.2f}")

print()
print("(2) local maxima strictly away from Gamma (exceed all 8 neighbors, |q|>0.3pi):")
maxima = []
for i in range(N):
    for j in range(N):
        r = np.hypot(qx[i], qy[j])
        if r < 0.3 * np.pi:
            continue
        v = xi[i, j]
        neigh = [xi[(i+di) % N, (j+dj) % N] for di in (-1,0,1) for dj in (-1,0,1) if (di,dj) != (0,0)]
        if v > max(neigh):
            maxima.append((v, qx[i], qy[j]))
maxima.sort(reverse=True)
for v, qxv, qyv in maxima[:12]:
    print(f"  xi={v:10.2f}  q=({qxv/np.pi:+.4f},{qyv/np.pi:+.4f})pi  enhancement={v/xi.mean():.3f}x")
