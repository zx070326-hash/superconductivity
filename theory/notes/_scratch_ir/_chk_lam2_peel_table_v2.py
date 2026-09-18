"""[numerical-evidence] 12-cell DK table, PEELED bound with a per-cell weight.

lambda_2(A) <= P_beta(varpi) := sup_{n>=1} [ sum_{m>=1} B_nm p_m/p_n - Delta_nn ],
p_n = (2n+1)^{-beta/2}.  Analytic gap g_an = h - min(P_beta, THETA_STAR, 0.6h).
One beta per cell (like the existing per-cell alpha), chosen by the cell's worst
point.

Direction discipline: P is an UPPER bound on lambda_2 => LOWER bound on the gap
=> the DK ratio can only come out WORSE than truth.  Tail of the m-sum is padded
UPWARD, never dropped; the sup over n is truncated at NROW and the argmax is
printed so the truncation can be seen to be non-binding.

Run:  OMP_NUM_THREADS=1 python _chk_lam2_peel_table_v2.py
"""
import os, io, contextlib, importlib
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

MT = 20000
NROW = 150
BETAS = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0])


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

CELLS = [(0.30, 0.42, 5.0), (0.42, 0.60, 3.6), (0.60, 0.85, 2.7), (0.85, 1.20, 2.1),
         (1.20, 1.70, 1.8), (1.70, 2.40, 1.6), (2.40, 3.40, 1.5), (3.40, 4.80, 1.42),
         (4.80, 6.80, 1.35), (6.80, 9.60, 1.30), (9.60, 13.6, 1.27), (13.6, 20.0, 1.24)]


def peel_all_betas(w):
    """Returns (P[beta] for beta in BETAS, argmax n for each beta)."""
    x = w * w
    f = lambda t: x / (t * t + x)
    m = np.arange(1, MT, dtype=float)
    um = 1.0 / np.sqrt(2.0 * m + 1.0)
    pm = (2.0 * m + 1.0) ** (-BETAS[:, None] / 2.0)          # (nb, MT-1)
    kk = np.arange(0, NROW + 1, dtype=float); fk = f(kk); fk[0] = 0.0
    cum = np.cumsum(fk)
    best = np.full(len(BETAS), -np.inf)
    arg = np.zeros(len(BETAS), dtype=int)
    for lo in range(1, NROW, 25):
        hi = min(lo + 25, NROW)
        ii = np.arange(lo, hi, dtype=float)[:, None]
        t = f(np.abs(ii - m[None, :])) + f(ii + m[None, :] + 1.0)
        for r, i in enumerate(range(lo, hi)):
            t[r, i - 1] = f(float(2 * i + 1))
        un = 1.0 / np.sqrt(2.0 * ii[:, 0] + 1.0)
        idx = np.arange(lo, hi)
        core = (t * um[None, :])                              # (nr, MT-1)
        for bi in range(len(BETAS)):
            pn = (2.0 * ii[:, 0] + 1.0) ** (-BETAS[bi] / 2.0)
            rs = un * (core @ pm[bi]) / pn
            pad = un * 2.0 * x / (np.sqrt(2.0 * MT + 1.0) * (MT - ii[:, 0] - 1.0))
            v = rs + pad - 2.0 * cum[idx] / (2.0 * idx + 1.0)
            j = int(np.argmax(v))
            if v[j] > best[bi]:
                best[bi] = v[j]; arg[bi] = idx[j]
    return best, arg


def pieces(w, al, N=600):
    x = w * w
    A = O(lambda k: x / (k * k + x), N)
    Ophi = O(lambda k: k * k / (k * k + x) ** 3, N)
    Opsi = O(lambda k: 1.0 / (k * k + x) ** 2, N)
    h = np.linalg.eigvalsh(A)[-1]
    n = np.arange(N); a = al * w
    t = (1 / np.sqrt(2 * n + 1.)) * a * a / (a * a + (n + 0.5) ** 2)
    t /= np.linalg.norm(t)
    res = np.linalg.norm(A @ t - (t @ (A @ t)) * t)
    Ft = t @ (Ophi @ t); e_phi = np.linalg.norm(Ophi @ t - Ft * t)
    Kt = Opsi @ t; qpsi = t @ Kt; e_psi = np.linalg.norm(Kt - qpsi * t)
    nphi = np.abs(np.linalg.eigvalsh(Ophi)).max()
    npsi = np.abs(np.linalg.eigvalsh(Opsi)).max()
    return dict(x=x, h=h, res=res, Ft=Ft, e_phi=e_phi, qpsi=qpsi,
                e_psi=e_psi, nphi=nphi, npsi=npsi)


def ratio(pc, l2b):
    g = pc["h"] - l2b
    if g <= 0:
        return np.inf
    s = min(pc["res"] / g, 1.0)
    F_lo = (1 - s * s) * pc["Ft"] - 2 * s * pc["e_phi"] - s * s * pc["nphi"]
    E_hi = pc["e_psi"] + s * (abs(pc["qpsi"]) + pc["npsi"])
    return pc["x"] ** 2 * E_hi ** 2 / (g * F_lo) if F_lo > 0 else np.inf


if __name__ == "__main__":
    print("PEELED lambda_2 bound with per-cell weight beta, fed into the "
          "12-cell DK criterion   [numerical-evidence]")
    print(f"MT={MT}  NROW={NROW}  betas={list(BETAS)}")
    print()
    hdr = (f"{'cell':>15} {'a':>5} {'b*':>5} | {'old':>8} {'new':>8} {'gain':>7} "
           f"| {'P worst':>8} {'Theta*':>7} {'argmax':>7}")
    print(hdr); print("-" * len(hdr))
    wo = wn = 0.0
    for wa, wb, al in CELLS:
        pts = np.geomspace(wa, wb, 11)
        pcs = [pieces(w, al) for w in pts]
        Ps = []; args = []
        for w in pts:
            b, a_ = peel_all_betas(w); Ps.append(b); args.append(a_)
        Ps = np.array(Ps); args = np.array(args)
        co = max(ratio(pc, 0.6 * pc["h"] if w_ <= 2 else THETA_STAR)
                 for pc, w_ in zip(pcs, pts))
        bestc, bi = np.inf, 0
        for j in range(len(BETAS)):
            c = max(ratio(pc, min(Ps[i, j], THETA_STAR, 0.6 * pc["h"]))
                    for i, pc in enumerate(pcs))
            if c < bestc:
                bestc, bi = c, j
        wo = max(wo, co); wn = max(wn, bestc)
        print(f"[{wa:5.2f},{wb:5.2f}] {al:5.2f} {BETAS[bi]:5.2f} | {co:8.4f} "
              f"{bestc:8.4f} {1 - bestc / co:7.1%} | {Ps[:, bi].max():8.4f} "
              f"{THETA_STAR:7.4f} {args[:, bi].max():7d}")
    print("-" * len(hdr))
    print(f"{'12-cell worst':>15} {'':5} {'':5} | {wo:8.4f} {wn:8.4f} "
          f"{1 - wn / wo:7.1%}")
