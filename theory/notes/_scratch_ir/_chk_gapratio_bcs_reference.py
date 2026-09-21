"""独立核 gapratio 的 "BCS limit check fails" 断言。
它测 Einstein lambda=0.5, w=50meV, mu*=0，得 2D/kTc=4.05，与 3.53 比、判 +15% BLOCKING。

问题：3.53 是 lambda->0 的 BCS 极限，不是 lambda=0.5 的参考值。
lambda=0.5 已有可观的强耦合修正。正确参考是 Carbotte (RMP 62, 1027 (1990)) 的
    2D0/kTc = 3.53 [1 + 12.5 (Tc/wlog)^2 ln(wlog/(2Tc))]
这里先用自建 KAY 算子解出 Tc，再代入。算子与 _chk_prop1_independent.py 同一份代码。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.optimize import brentq


def O(phi_vals, N):
    n = np.arange(N)
    I = n[:, None]; J = n[None, :]
    T = phi_vals[np.abs(I - J)].copy()
    np.fill_diagonal(T, 0.0)
    Hk = phi_vals[I + J + 1]
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    M = (T + Hk) * u[:, None] * u[None, :]
    cum = np.concatenate(([0.0], np.cumsum(phi_vals[1:N])))
    M[n, n] -= 2.0 * cum / (2.0 * n + 1.0)
    return M


def lam1(varpi, N=800):
    k = np.arange(0, 2 * N, dtype=float)
    f = varpi**2 / (k * k + varpi**2)       # f(0)=1 自动
    return np.linalg.eigvalsh(O(f, N))[-1]


for lam in (0.5, 0.3, 0.2, 0.1):
    target = 1.0 / lam
    vp = brentq(lambda v: lam1(v) - target, 0.05, 400.0, xtol=1e-10)
    wE = 50.0                               # meV, Einstein => wlog = wE
    Tc = wE / (2 * np.pi * vp)              # meV
    t = Tc / wE
    ratio = 3.53 * (1 + 12.5 * t**2 * np.log(1.0 / (2 * t)))
    print(f"lambda={lam:<5} varpi_c={vp:9.5f}  Tc={Tc:8.5f} meV = {Tc*11.6045:7.3f} K"
          f"   Tc/wlog={t:.5f}   Carbotte 2D/kTc={ratio:.4f}")
