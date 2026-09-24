# -*- coding: utf-8 -*-
"""同信息类比较对库仑截断约定 r 的稳健性。（v2：降规模，加 N 收敛检查）

两个天花板都只吃 (λ, M2)，都给出「η/M=9.172、T=300K 下所需的 λ 下限」：
  ·我们：λ*(r) = Einstein 族上 λ₁=1 的根（②(c) 已证族 sup 由 Einstein 取到）
  ·Leavens：kTc ≤ c(μ*,r)·A + Cauchy–Schwarz A ≤ √(λM2)/2 ⟹ λ_area = 1/(π²c²M2)
    无量纲化：A/(2πTc)=λω̄₂/2 ⟹ kTc/A = 1/(π λ ω̄₂_root(λ;r))，c = max_λ 该式。
方向：两者都是「真所需 λ」的**下**估 ⟹ 数越大排除越强，两侧同向可比。
⚠ 台阶：N_c=⌊r·ω̄₂⌋ 不连续 ⟹ λ₁=1 可能无根；检出即标记，不用二分端点冒充根。
严格性：numerical-evidence（float64）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import numpy as np

_d = os.path.dirname(os.path.abspath(__file__))
def _load(n):
    s = importlib.util.spec_from_file_location(n, os.path.join(_d, n+".py"))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
KT = _load("_chk_main_kato_temple_two_atom_prototype_2026-09-23")
M2, mus = 1.453151466, 0.13


def l1(L, r_, N):
    wl = float(np.sqrt(M2/L))
    return KT.lam_k([L], [wl], mus, int(np.floor(r_*wl)), N), wl


def root(r_, N):
    lo, hi = 1.5, 8.0
    for _ in range(50):
        mid = 0.5*(lo+hi)
        if l1(mid, r_, N)[0] < 1.0: lo = mid
        else: hi = mid
    return 0.5*(lo+hi), l1(lo, r_, N)[0], l1(hi, r_, N)[0]


def cmax(r_, N, grid):
    best = (-9.0, None)
    for L in grid:
        try: wr = KT.einstein_root(L, mus, r_, N)
        except Exception: continue
        v = 1.0/(np.pi*L*wr)
        if v > best[0]: best = (v, L)
    if best[1] is None: return best
    lo, hi = best[1]-0.15, best[1]+0.15          # 细化
    for L in np.linspace(lo, hi, 31):
        try: wr = KT.einstein_root(L, mus, r_, N)
        except Exception: continue
        v = 1.0/(np.pi*L*wr)
        if v > best[0]: best = (v, L)
    return best


N = 400
G = np.linspace(1.2, 6.0, 25)
print(f"N 收敛检查（r=10）：", flush=True)
for NN in (300, 400, 600):
    c_, La_ = cmax(10, NN, G); Ls_, _, _ = root(10, NN)
    print(f"  N={NN}: c={c_:.8f} (argmax λ={La_:.3f})  λ*={Ls_:.6f}", flush=True)

print(f"\n{'r':>5}{'λ*(我们)':>12}{'台阶?':>8}{'c(μ*,r)':>12}{'argmax λ':>10}"
      f"{'λ_area':>10}{'λ*/λ_area':>12}", flush=True)
for r_ in [4, 5, 6, 7, 8, 10, 12, 15, 20, 30]:
    Ls, v_lo, v_hi = root(r_, N)
    gap = "跳过" if min(abs(v_lo-1), abs(v_hi-1)) > 1e-6 else "—"
    c, La = cmax(r_, N, G)
    la = 1.0/(np.pi**2*c*c*M2)
    print(f"{r_:>5}{Ls:>12.6f}{gap:>8}{c:>12.8f}{La:>10.3f}{la:>10.5f}"
          f"{Ls/la:>12.4f}", flush=True)
