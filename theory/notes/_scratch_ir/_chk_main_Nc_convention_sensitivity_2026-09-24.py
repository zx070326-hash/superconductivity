# -*- coding: utf-8 -*-
"""N_c 约定敏感性：②(c) 的 certificate 是在**固定** N_c=⌊r·ω̄₂(Einstein)⌋=6 下取得的。
本脚本测两件事：
  (1) 固定 N_c 时，Einstein 是否仍是族内极值成员；
  (2) 族 sup λ₁ 随 N_c 怎么变 —— 即「λ=3.29 被排除」这一判决对 N_c 取值有多敏感。
⚠ 前提订正：λ₁ 对 N_c **非单调**（u_c 随 N_c 换方向，两个秩一之差不定号；
   顶本征向量出节点后加正分量反而减小 ⟨u_c,v⟩）。故不得用单调性省掉任何一格。
严格性：全部 numerical-evidence（float64，N=400 截断，有限抽样）。
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
print(f"λ={lam}  Einstein ω̄₂={w:.9f}   r=10 ⟹ N_c=6   截断 N={N}", flush=True)
rng = np.random.default_rng(11)
pts = [(rng.uniform(1e-4, 1-1e-4), rng.uniform(1e-4, 1-1e-4)) for _ in range(300)]
pts += [(q, s) for q in np.linspace(.02, .98, 12) for s in np.linspace(.02, .98, 12)]
print(f"抽样 {len(pts)} 点/每个 N_c", flush=True)
print(f"\n{'N_c':>5}{'λ₁(Einstein)':>16}{'族内 max λ₁':>16}{'argmax (q,σ)':>24}{'max−E':>13}", flush=True)
for Nc in [0, 3, 4, 5, 6, 7, 8, 10, 15, 25]:
    lE = KT.lam_k([lam], [w], mus, Nc, N)
    best = (-9.0, None)
    for q, s in pts:
        L_, W_ = QS.mem_qs(lam, w, q, s)
        v = KT.lam_k(L_, W_, mus, Nc, N)
        if v > best[0]:
            best = (v, (q, s))
    flag = "  <-- 族 sup > 1" if best[0] > 1 else ""
    print(f"{Nc:>5}{lE:>16.9f}{best[0]:>16.9f}   ({best[1][0]:.6f},{best[1][1]:.6f})"
          f"{best[0]-lE:>13.3e}{flag}", flush=True)
