"""Acceptance-gate check (team-lead requirement, 2026-09-22): does the corrected
Graser H(k) reproduce the qualitative FS topology (Gamma-centered HOLE pockets,
M-centered ELECTRON pockets) claimed in the paper's Fig.1-2, not just "some band
crosses mu near this point" (which the earlier _graser_find_mu_fermi_surface script
tested and is too weak to distinguish hole vs electron character)?

Method: at each high-symmetry point (Gamma, X, M), find the band whose energy is
nearest mu, then compute the 2D curvature (trace of the finite-difference Hessian
along kx,ky) of that band there.
  - band value ABOVE mu at the point + NEGATIVE curvature (concave, local max)
    => hole pocket (occupied only in a small disk around this point, since the band
       dips below mu moving outward... wait: need band value < mu growing outward)
  - more robust discriminator used below: sign of (E(k0)-mu) combined with sign of
    curvature. Hole pocket <=> band has a LOCAL MAXIMUM at k0 that pokes above mu
    (curvature < 0, E(k0) can be above OR just below mu -- what matters for "hole
    pocket AT k0" is that occupied states (E<mu) form a small pocket CENTERED at k0,
    i.e. E(k0) < mu and E increases away from k0 => curvature > 0 is WRONG framing;
    correct: occupied disk around k0 <=> E(k) - mu < 0 near k0 and > 0 further out
    <=> E(k0) < mu (interior of pocket) with band curving UP away from k0 (curvature>0)
    ... but that describes an electron pocket by the usual convention (small pocket of
    conduction states). Use the standard convention instead: HOLE pocket = pocket of
    band maxima that lie just above mu, i.e. unoccupied-outside/occupied-inside is
    backwards -- resolve by explicit reasoning documented in comments printed below.)
"""
import importlib.util, os
import numpy as np

_here = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    "graser1003_hamiltonian", os.path.join(_here, "_graser1003_hamiltonian_2026-09-22.py"))
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
build_H_graser = _mod.build_H_graser

NBAND = 5

def bands(kx, ky):
    return np.linalg.eigvalsh(build_H_graser(kx, ky, 0.0))

# recompute mu fresh (post-correction H(k), N=128, matches nesting-scan script)
N = 128
ks = np.linspace(-np.pi, np.pi, N, endpoint=False)
allb = np.empty((N, N, NBAND))
for i, kx in enumerate(ks):
    for j, ky in enumerate(ks):
        allb[i, j, :] = bands(kx, ky)
flat = allb.reshape(-1, NBAND)
lo, hi = flat.min() - 1, flat.max() + 1
for _ in range(60):
    mu = 0.5 * (lo + hi)
    n = 2 * np.mean(flat < mu) * NBAND
    if n < 6.0:
        lo = mu
    else:
        hi = mu
mu = 0.5 * (lo + hi)
print(f"mu (6 e/Fe, N={N}) = {mu:.4f} eV")

pts = {"Gamma": (0.0, 0.0), "X=(pi,0)": (np.pi, 0.0), "Y=(0,pi)": (0.0, np.pi),
       "M=(pi,pi)": (np.pi, np.pi)}

h = 0.02
for name, (kx0, ky0) in pts.items():
    e0 = bands(kx0, ky0)
    iband = np.argmin(np.abs(e0 - mu))
    e_val = e0[iband]

    def eb(kx, ky, n=iband):
        return bands(kx, ky)[n]

    d2x = (eb(kx0 + h, ky0) - 2 * e_val + eb(kx0 - h, ky0)) / h**2
    d2y = (eb(kx0, ky0 + h) - 2 * e_val + eb(kx0, ky0 - h)) / h**2
    curv = d2x + d2y  # trace of Hessian (isotropic proxy)
    # standard convention: near a band EXTREMUM, E(k) ~ E0 + (curv/2)|k-k0|^2.
    # A pocket of OCCUPIED states (E<mu) surrounding k0 requires E0<mu and curv>0
    # (band curves UP away from k0, so E<mu only in the interior disk) -- but that is
    # by definition a region of extra occupied states relative to the "sea" outside,
    # i.e. an ELECTRON-like Fermi surface only if k0 sits in an otherwise-empty band
    # (conduction band minimum poking below mu). Symmetrically, a HOLE pocket is a
    # small disk of EMPTY states (E>mu) surrounding a band MAXIMUM that pokes above
    # mu: E0>mu, curv<0, and E>mu only in the interior disk (states there are holes
    # relative to an otherwise-filled band).
    if e_val < mu and curv > 0:
        kind = "electron-like (band min pokes below mu, curv>0)"
    elif e_val > mu and curv < 0:
        kind = "hole-like (band max pokes above mu, curv<0)"
    else:
        kind = f"AMBIGUOUS/saddle (E-mu={e_val-mu:+.4f}, curv={curv:+.3f}) -- not a clean pocket extremum at this exact point"
    print(f"{name:10s} band#{iband} E-mu={e_val-mu:+.4f}  curv(trace Hess)={curv:+.4f}  -> {kind}")
