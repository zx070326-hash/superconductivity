# -*- coding: utf-8 -*-
"""两个修补：
(A) 上一版的 `max−E` 列自由度为 0 —— argmax 落在 q=σ 对角线，而**整条对角线就是 Einstein**
    （QS 坐标：σ=q ⟹ ω₁=ω₂=ω̄₂）。改报**非对角**成员的 max 减 Einstein，并印 argmax。
(B) 判决对库仑截断的敏感性：阈值 λ*(r) —— 固定 η/M（即 M2）与 μ*，
    沿 Einstein 族解 λ₁=1 的 λ，扫 r。N_c=⌊r·ω̄₂⌋ 是**台阶**函数 ⟹ λ*(r) 分段常数式跳变。
⚠ λ₁ 对 N_c **非单调**（u_c 随 N_c 换方向），不得用单调性省格。
严格性：全部 numerical-evidence（float64，N=400/600，有限抽样）。
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
QS = _load("_chk_main_qsigma_coords_and_deficit_identity_2026-09-23")

M2, mus, lam = 1.453151466, 0.13, 3.29
w = float(np.sqrt(M2/lam)); N = 400
rng = np.random.default_rng(11)
pts = [(rng.uniform(1e-4, 1-1e-4), rng.uniform(1e-4, 1-1e-4)) for _ in range(300)]
pts += [(q, s) for q in np.linspace(.02, .98, 12) for s in np.linspace(.02, .98, 12)
        if abs(q-s) > 1e-12]
print(f"(A) 非对角成员 vs Einstein，λ={lam}，{len(pts)} 个**非对角**点（|q−σ|>0）", flush=True)
print(f"{'N_c':>5}{'λ₁(Einstein)':>16}{'非对角 max':>16}{'argmax (q,σ)':>24}{'非对角−E':>13}", flush=True)
for Nc in [0, 3, 4, 5, 6, 7, 8, 10, 15, 25]:
    lE = KT.lam_k([lam], [w], mus, Nc, N)
    best = (-9.0, None)
    for q, s in pts:
        L_, W_ = QS.mem_qs(lam, w, q, s)
        v = KT.lam_k(L_, W_, mus, Nc, N)
        if v > best[0]:
            best = (v, (q, s))
    print(f"{Nc:>5}{lE:>16.9f}{best[0]:>16.9f}   ({best[1][0]:.6f},{best[1][1]:.6f})"
          f"{best[0]-lE:>13.3e}{'  <-- 超过 E' if best[0] > lE else ''}", flush=True)

print(f"\n(B) 阈值 λ*(r)：固定 M2={M2}, μ*={mus}，沿 Einstein 解 λ₁=1", flush=True)
def lam1_E(L, r_):
    wl = float(np.sqrt(M2/L)); return KT.lam_k([L], [wl], mus, int(np.floor(r_*wl)), 600), wl
print(f"{'r':>6}{'λ*':>12}{'ω̄₂(λ*)':>12}{'N_c':>5}{'λ₁ 检验':>14}", flush=True)
for r_ in [4, 5, 6, 8, 10, 12, 15, 20, 30]:
    lo, hi = 1.5, 8.0
    for _ in range(60):
        mid = 0.5*(lo+hi)
        if lam1_E(mid, r_)[0] < 1.0: lo = mid
        else: hi = mid
    L = 0.5*(lo+hi); v, wl = lam1_E(L, r_)
    print(f"{r_:>6}{L:>12.6f}{wl:>12.6f}{int(np.floor(r_*wl)):>5}{v:>14.9f}", flush=True)
