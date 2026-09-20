"""[numerical-evidence] What does SOUNDNESS actually cost, measured in the unit
the baseline script was built to measure: number of CERT-B cells?

Established:
  1. _chk_nphi_convergence.py : the N=600 eigvalsh values ARE the true norms to
     13-15 digits => the baseline number is right at numerical-evidence tier;
     what is missing is a PROOF that it is an upper bound.
  2. _chk_nphi_budget.py      : at the per-w optimal (alpha,beta), the ratio is
     0.28-0.31 at kappa=1 and only 0.42-0.45 at kappa=3.24 -- nowhere near
     TARGET=0.90.  So the binding constraint is CELL WIDTH (one (alpha,beta) per
     cell must cover an interval), not certificate slack.

So the honest question is not "do 3 cells fail" but "how many more cells does a
sound certificate cost".  Re-run the SAME greedy with nphi,npsi -> kappa*norm.

Direction: ratio is increasing in kappa (nphi enters a lower bound negatively,
npsi enters an upper bound positively) => larger kappa = harder => a PASS at
kappa is inherited by every enclosure with constant <= kappa.
"""
import os, io, contextlib, importlib, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
with contextlib.redirect_stdout(io.StringIO()):
    gc = importlib.import_module("_chk_grid_coarsen")

BETAS, ALPHAS, NGRID, TARGET = gc.BETAS, gc.ALPHAS, gc.NGRID, gc.TARGET
THETA_STAR = gc.THETA_STAR
WS = np.geomspace(0.30, 20.0, NGRID)


def base(w, N=600):
    x = w * w
    A = gc.O(lambda k: x / (k * k + x), N)
    Ophi = gc.O(lambda k: k * k / (k * k + x) ** 3, N)
    Opsi = gc.O(lambda k: 1.0 / (k * k + x) ** 2, N)
    return dict(x=x, N=N, A=A, Ophi=Ophi, Opsi=Opsi,
                h=np.linalg.eigvalsh(A)[-1],
                n0phi=np.abs(np.linalg.eigvalsh(Ophi)).max(),
                n0psi=np.abs(np.linalg.eigvalsh(Opsi)).max())


def ratio(bs, al, l2b, kap):
    g = bs["h"] - l2b
    if g <= 0:
        return np.inf
    nphi, npsi = kap * bs["n0phi"], kap * bs["n0psi"]
    N = bs["N"]; n = np.arange(N); a = al * np.sqrt(bs["x"])
    t = (1 / np.sqrt(2 * n + 1.)) * a * a / (a * a + (n + 0.5) ** 2)
    t /= np.linalg.norm(t)
    At = bs["A"] @ t
    res = np.linalg.norm(At - (t @ At) * t)
    Pt = bs["Ophi"] @ t; Ft = t @ Pt; e_phi = np.linalg.norm(Pt - Ft * t)
    Kt = bs["Opsi"] @ t; qpsi = t @ Kt; e_psi = np.linalg.norm(Kt - qpsi * t)
    s = min(res / g, 1.0)
    F_lo = (1 - s * s) * Ft - 2 * s * e_phi - s * s * nphi
    E_hi = e_psi + s * (abs(qpsi) + npsi)
    return bs["x"] ** 2 * E_hi ** 2 / (g * F_lo) if F_lo > 0 else np.inf


def greedy(R):
    cells, i = [], 0
    while i < NGRID - 1:
        j, best = i, None
        while j + 1 < NGRID:
            m = R[i:j + 2].max(axis=0)
            if m.min() > TARGET:
                break
            best = (j + 1, float(m.min())); j += 1
        if best is None:
            cells.append((WS[i], WS[i + 1], np.inf)); i += 1; continue
        cells.append((WS[i], WS[best[0]], best[1])); i = best[0]
    return cells


print(f"precomputing {NGRID} grid points ...", flush=True)
BS = [base(w) for w in WS]
PB = np.array([gc.peel_all_betas(w) for w in WS])
print("done.\n")

print(f"cells needed to cover [0.30, 20.0] at ratio <= {TARGET}")
print("kappa = certified_bound / true_norm   (kappa=1 is the unsound baseline)\n")
print("  kappa   cells   worst-cell ratio   infeasible cells")
print("  " + "-" * 52)
for kap in (1.0, 1.5, 2.0, 2.5, 3.24, 5.0, 10.0):
    R = np.empty((NGRID, len(ALPHAS), len(BETAS)))
    for i, bs in enumerate(BS):
        for ai, al in enumerate(ALPHAS):
            for bi in range(len(BETAS)):
                R[i, ai, bi] = ratio(bs, al,
                                     min(PB[i, bi], THETA_STAR, 0.6 * bs["h"]),
                                     kap)
    cells = greedy(R)
    inf_n = sum(1 for c in cells if not np.isfinite(c[2]))
    fin = [c[2] for c in cells if np.isfinite(c[2])]
    print(f"  {kap:5.2f}   {len(cells):5d}   {max(fin) if fin else float('nan'):16.4f}"
          f"   {inf_n:16d}")
