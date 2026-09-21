"""[numerical-evidence] team-lead 从零独立重建：M2^Delta=0 反例的 lam1(OO)。

争议：线 C agent 报 lam1 = 0.123471，仓库三处笔记写 0.115511（差 6.9%）。
本脚本不 import 任何仓库脚本，逐字按 P4_g2_analytic.md §1 的算子定义重写，
并分别构造两个算子：
  (I)  单核 O[f]，f = Lam^ph - Lam^sf        <- 线 C agent 用的
  (II) 双核 OO(Lam^Z=ph+sf 驱动对角, Lam^D=ph-sf 驱动配对)  <- 仓库反例用的
校准：O[1/k^2] 的 lam1 必须复现 g2 = 1.3181402121。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np


def u_vec(N):
    return 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)


def O_single(f, N):
    """O[f]_{nm} = -d_{nm} 2 u_n^2 sum_{k=1..n} f(k) + (1-d_{nm}) f(|n-m|) u_n u_m
                   + f(n+m+1) u_n u_m.   f is an array indexed by k, len >= 2N."""
    u = u_vec(N)
    n = np.arange(N)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    uu = np.outer(u, u)
    M = uu * (np.where(D == 0, 0.0, f[D]) + f[S])
    cum = np.concatenate([[0.0], np.cumsum(f[1:N])])   # cum[n] = sum_{k=1..n} f(k)
    M[n, n] -= 2.0 * u ** 2 * cum
    return M


def O_two(fZ, fD, N):
    """Z channel drives the diagonal row-sum subtraction, D channel drives pairing."""
    u = u_vec(N)
    n = np.arange(N)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    uu = np.outer(u, u)
    M = uu * (np.where(D == 0, 0.0, fD[D]) + fD[S])
    M[n, n] += u ** 2 * fD[0]
    cum = np.concatenate([[0.0], np.cumsum(fZ[1:N])])
    M[n, n] -= u ** 2 * (fZ[0] + 2.0 * cum)
    return M


def lam1(M):
    return float(np.linalg.eigvalsh(M)[-1])


def brk(k, varpi):
    return varpi ** 2 / (k ** 2 + varpi ** 2)


# ---- calibration: g2 = lam1(O[1/k^2]) -------------------------------------
print("calibration  lam1(O[1/k^2])  (repo g2 = 1.3181402121)")
for N in (400, 800, 1600, 3200):
    k = np.arange(2 * N + 2, dtype=float)
    f = np.zeros_like(k)
    f[1:] = 1.0 / k[1:] ** 2
    print(f"  N={N:5d}  lam1={lam1(O_single(f, N)):.10f}")

# ---- the counterexample ---------------------------------------------------
w_ph, w_sf = 0.10, 5.0
lam_ph = 10.0
M2_ph = lam_ph * w_ph ** 2
lam_sf = M2_ph / w_sf ** 2
M2_sf = lam_sf * w_sf ** 2
print(f"\nM2_ph={M2_ph:.6f}  M2_sf={M2_sf:.6f}  M2^D={M2_ph - M2_sf:.3e}")

print("\n   N   lam1 single-kernel O[Lam^D]   lam1 two-kernel OO(Z,D)")
for N in (100, 200, 400, 800, 1600, 3200):
    k = np.arange(2 * N + 2, dtype=float)
    lph = lam_ph * brk(k, w_ph)
    lsf = lam_sf * brk(k, w_sf)
    fD = lph - lsf
    fZ = lph + lsf
    print(f"{N:5d}   {lam1(O_single(fD, N)):.10f}                {lam1(O_two(fZ, fD, N)):.10f}")
