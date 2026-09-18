"""[numerical-evidence] How much can the Lemma-C grid be COARSENED with the new
lambda_2 bound?

CERT-B cost ~ (number of cells) x (rationals per cell).  The peeled bound adds
one rational per cell (beta) but buys a 3.5x headroom in the DK ratio.  This
script asks the only question that matters for cost: with g_an = h - min(P_beta,
THETA*, 0.6h), how few cells cover [0.30, 20.0] at ratio <= TARGET?

Greedy from the left: extend the right edge over a fixed geometric varpi grid as
far as the best (alpha, beta) pair keeps every grid point in the cell below
TARGET.  Compared against the same greedy run with the OLD gap.

Run:  OMP_NUM_THREADS=1 python _chk_grid_coarsen.py
"""
import os, io, contextlib, importlib
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

MT, NROW = 20000, 300  # NROW=300: the analytic tail bound V_beta <= C0*w^2/(2n+1),
# C0=16.510148, only dominates the peak beyond n=300 (at 150 it does not; worst
# bd/peak = 1.37).  See theta_tilde_keeping_D_2026-09-17.md section 12.4.
BETAS = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0])
ALPHAS = np.concatenate([np.geomspace(1.10, 6.0, 22)])
NGRID = 49
TARGET = 0.90


def O(phi, N):
    n = np.arange(N); I = n[:, None]; J = n[None, :]
    pv = phi(np.arange(0, 2 * N + 2, dtype=float))
    T = pv[np.abs(I - J)].copy(); np.fill_diagonal(T, 0.0)
    u = 1 / np.sqrt(2 * n + 1.)
    M = (T + pv[I + J + 1]) * u[:, None] * u[None, :]
    M[n, n] -= 2.0 * np.concatenate(([0.], np.cumsum(pv[1:N]))) / (2 * n + 1.)
    return M


with contextlib.redirect_stdout(io.StringIO()):
    THETA_STAR = importlib.import_module("lemmaC_compact_probe").THETA_STAR


def peel_all_betas(w):
    x = w * w
    f = lambda t: x / (t * t + x)
    m = np.arange(1, MT, dtype=float)
    um = 1.0 / np.sqrt(2.0 * m + 1.0)
    pm = (2.0 * m + 1.0) ** (-BETAS[:, None] / 2.0)
    kk = np.arange(0, NROW + 1, dtype=float); fk = f(kk); fk[0] = 0.0
    cum = np.cumsum(fk)
    best = np.full(len(BETAS), -np.inf)
    for lo in range(1, NROW, 25):
        hi = min(lo + 25, NROW)
        ii = np.arange(lo, hi, dtype=float)[:, None]
        t = f(np.abs(ii - m[None, :])) + f(ii + m[None, :] + 1.0)
        for r, i in enumerate(range(lo, hi)):
            t[r, i - 1] = f(float(2 * i + 1))
        un = 1.0 / np.sqrt(2.0 * ii[:, 0] + 1.0)
        idx = np.arange(lo, hi)
        core = t * um[None, :]
        for bi in range(len(BETAS)):
            pn = (2.0 * ii[:, 0] + 1.0) ** (-BETAS[bi] / 2.0)
            rs = un * (core @ pm[bi]) / pn
            pad = un * 2.0 * x / (np.sqrt(2.0 * MT + 1.0) * (MT - ii[:, 0] - 1.0))
            best[bi] = max(best[bi],
                           (rs + pad - 2.0 * cum[idx] / (2.0 * idx + 1.0)).max())
    return best


def base(w, N=600):
    x = w * w
    A = O(lambda k: x / (k * k + x), N)
    Ophi = O(lambda k: k * k / (k * k + x) ** 3, N)
    Opsi = O(lambda k: 1.0 / (k * k + x) ** 2, N)
    return dict(x=x, N=N, A=A, Ophi=Ophi, Opsi=Opsi,
                h=np.linalg.eigvalsh(A)[-1],
                nphi=np.abs(np.linalg.eigvalsh(Ophi)).max(),
                npsi=np.abs(np.linalg.eigvalsh(Opsi)).max())


def ratio(bs, al, l2b):
    g = bs["h"] - l2b
    if g <= 0:
        return np.inf
    N = bs["N"]; n = np.arange(N); a = al * np.sqrt(bs["x"])
    t = (1 / np.sqrt(2 * n + 1.)) * a * a / (a * a + (n + 0.5) ** 2)
    t /= np.linalg.norm(t)
    At = bs["A"] @ t
    res = np.linalg.norm(At - (t @ At) * t)
    Pt = bs["Ophi"] @ t; Ft = t @ Pt; e_phi = np.linalg.norm(Pt - Ft * t)
    Kt = bs["Opsi"] @ t; qpsi = t @ Kt; e_psi = np.linalg.norm(Kt - qpsi * t)
    s = min(res / g, 1.0)
    F_lo = (1 - s * s) * Ft - 2 * s * e_phi - s * s * bs["nphi"]
    E_hi = e_psi + s * (abs(qpsi) + bs["npsi"])
    return bs["x"] ** 2 * E_hi ** 2 / (g * F_lo) if F_lo > 0 else np.inf


if __name__ == "__main__":
    WS = np.geomspace(0.30, 20.0, NGRID)
    print(f"precomputing {NGRID} grid points ...", flush=True)
    BS = [base(w) for w in WS]
    PB = np.array([peel_all_betas(w) for w in WS])          # (NGRID, nbeta)
    print("done.\n")

    # ratio tables: R[mode][i, ai, bi]
    Rnew = np.empty((NGRID, len(ALPHAS), len(BETAS)))
    Rold = np.empty((NGRID, len(ALPHAS)))
    for i, (w, bs) in enumerate(zip(WS, BS)):
        l2_old = 0.6 * bs["h"] if w <= 2 else THETA_STAR
        for ai, al in enumerate(ALPHAS):
            Rold[i, ai] = ratio(bs, al, l2_old)
            for bi in range(len(BETAS)):
                Rnew[i, ai, bi] = ratio(bs, al,
                                        min(PB[i, bi], THETA_STAR, 0.6 * bs["h"]))

    def greedy(R, has_beta):
        cells, i = [], 0
        while i < NGRID - 1:
            j = i
            best = None
            while j + 1 < NGRID:
                seg = R[i:j + 2]
                if has_beta:
                    m = seg.max(axis=0)                 # (nalpha, nbeta)
                else:
                    m = seg.max(axis=0)                 # (nalpha,)
                if m.min() > TARGET:
                    break
                best = (j + 1, np.unravel_index(int(np.argmin(m)), m.shape),
                        float(m.min()))
                j += 1
            if best is None:
                cells.append((WS[i], WS[i + 1], None, np.inf)); i += 1; continue
            cells.append((WS[i], WS[best[0]], best[1], best[2]))
            i = best[0]
        return cells

    for tag, R, hb in (("OLD gap", Rold, False), ("NEW gap", Rnew, True)):
        cells = greedy(R, hb)
        print(f"=== {tag}:  {len(cells)} cells at target ratio <= {TARGET} ===")
        for lo, hi, idx, r in cells:
            if idx is None:
                print(f"  [{lo:6.3f},{hi:6.3f}]  INFEASIBLE")
            elif hb:
                print(f"  [{lo:6.3f},{hi:6.3f}]  a={ALPHAS[idx[0]]:5.2f} "
                      f"b={BETAS[idx[1]]:4.2f}  worst={r:.4f}")
            else:
                print(f"  [{lo:6.3f},{hi:6.3f}]  a={ALPHAS[idx[0]]:5.2f} "
                      f"           worst={r:.4f}")
        print()
