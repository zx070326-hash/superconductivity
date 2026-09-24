# -*- coding: utf-8 -*-
"""预算优先（CLAUDE.md §1）：把 (Phi-A) 的抑制量 Delta 朝**判据变难**的方向退化
(Delta -> kappa*Delta, kappa<=1 ⟹ ceiling 变大 ⟹ Tc 上界变大 ⟹ lam 阈值变小)，
扫出 lam_ceil 在哪里跌破面积天花板阈值 lam_area=2.5382。

消费点：判词「(Phi-A) 的 Einstein 阈值是否高于 Leavens 面积阈值」。
若 kappa=0（**完全丢掉** Delta，即只用 mu*=0 的 AD75 解析极值性、不给库仑任何信用）
仍然 > 2.5382 ⟹ Delta 不承重，不要为「把 Delta 推广到一般谱」派活。
严格性：numerical-evidence (float64, N=400)。"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util, numpy as np
from scipy.linalg import eigh
spec = importlib.util.spec_from_file_location(
    "m", "theory/notes/_scratch_ir/_chk_main_leavens_c_and_phiA_rescale_2026-09-23.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

N, mus, r = 400, 0.13, 10.0
eta_over_M = 9.172
pref = m.HB * np.sqrt(eta_over_M)
LAM_AREA = 2.5382

def ceil_kappa(wbar, lam, kappa):
    """lam-标度下的 (Phi-A) ceiling，Delta 乘 kappa<=1（退化方向：ceiling 变大）。"""
    Nc = int(np.floor(r * wbar))
    Khat, u = m.pieces(wbar, N)
    ev, V = eigh(Khat)
    m1h, m2h = ev[-1], ev[-2]
    uc = np.where(np.arange(N) < Nc, u, 0.0)
    c1s = float(V[:, -1] @ uc) ** 2
    un = float(uc @ uc)
    den = lam * (m1h - m2h) - 2.0 * mus * un          # 前提 (PA)，正确标度
    if den <= 0.0:
        return None
    sp = 2.0 * mus * (un - c1s) / den
    return lam * m1h - kappa * 2.0 * mus * c1s / (1.0 + sp)

def Tc_K(lam, kappa):
    _, w = m.roots_in_wbar(lambda w: ceil_kappa(w, lam, kappa), 1.0, 0.05, 12.0)
    if w is None: return None
    return (pref / np.sqrt(lam)) / (2.0 * np.pi * w) / m.KB

def lam_ceil(kappa):
    f = lambda l: (Tc_K(l, kappa) or -1e9) - 300.0
    lo, hi = 1.05, 12.0
    if not (f(lo) < 0 < f(hi)): return None
    for _ in range(42):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

print("=" * 78)
print("(Phi-A) Delta 的预算：Delta -> kappa*Delta（kappa=0 即完全丢掉库仑抑制）")
print(f"面积天花板（一般谱）阈值 lam_area = {LAM_AREA}")
print(f"{'kappa':>7}{'lam_ceil':>12}{'vs lam_area':>14}")
for kap in (0.0, 0.25, 0.5, 0.75, 1.0):
    lc = lam_ceil(kap)
    tag = "胜" if (lc and lc > LAM_AREA) else "负"
    print(f"{kap:7.2f}{(('%.4f'%lc) if lc else 'None'):>12}{tag:>14}")

# 精确 Einstein Tc（不用 ceiling）作端到端方向门：ceiling 阈值必须 <= 精确阈值
def Tc_exact_K(lam):
    _, w = m.roots_in_wbar(lambda w: m.lam1_exact(w, lam, mus, r, N), 1.0, 0.05, 12.0)
    if w is None: return None
    return (pref / np.sqrt(lam)) / (2.0 * np.pi * w) / m.KB
f = lambda l: (Tc_exact_K(l) or -1e9) - 300.0
lo, hi = 1.05, 12.0
for _ in range(42):
    mid = 0.5 * (lo + hi)
    if f(mid) < 0: lo = mid
    else: hi = mid
lam_exact = 0.5 * (lo + hi)
print(f"\n端到端方向门：Einstein **精确** 300K 阈值 lam_exact = {lam_exact:.4f}")
print(f"  (Phi-A) kappa=1 阈值 = {lam_ceil(1.0):.4f} ；必须 <= lam_exact："
      f"{'OK' if lam_ceil(1.0) <= lam_exact + 1e-6 else '**方向违反**'}")
