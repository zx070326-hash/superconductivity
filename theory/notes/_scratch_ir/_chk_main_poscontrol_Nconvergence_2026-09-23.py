# -*- coding: utf-8 -*-
"""
主会话从零：**正对照的 N 收敛复核**（回应 teammate `extremality-hunt` 的质疑）。

━━ 被质疑的对象（承重）━━
笔记 `ceiling/mustar_extremality_two_routes_closed_2026-09-23.md` §6 写：
  「正对照必须亮：μ*=0.30、r=5、λ=0.540 → **12/323 违反**，max(λ₁−1)=+2.385e−04，
    最优 (w₁,w₂)/ω̄₂ = (0.95, 100)」
这条 PASS 是 §4「门 4 对失效是盲的」与 §6「判据区 0/3172 有意义」的**唯一**正对照。
若它是松原截断赝物，则：§6 的 0/3172 退化为 VOID（检验可能本来就盲），
§4 的「局部极大对失效盲」这一否定性结论**也失去支撑**（它要求正对照在门 3 亮、门 4 不亮）。

teammate 的质疑（值得当真）：最优点 w₂ = 100·ω̄₂ = 15000 meV，而 2πT ~ 0.3 meV
 ⟹ 需 N ~ w₂/(2πT) ~ 5e4 才能被松原网格分辨；且他实测 Tc_E 在 N=300→600 之间跳 3×。

━━ 本脚本判词设计（方向审计）━━
· 唯一判词：**λ₁(split) − 1 的符号是否随 N 稳定**。
  它若在 N 增大时翻负 ⟹ 正对照是赝物 ⟹ 上述两条结论必须撤回。
  它若保号且量级稳定 ⟹ 正对照成立，teammate 的 Tc_E 跳变是**另一个量**（解 Tc 的二分）
  的不收敛，不自动传染到**钉温比较**（两谱在同一 T、同一 N 下比，共模误差部分抵消）。
· ⚠ 这两件事必须分开测，不许用其一替另一：本脚本**同时**报
  (A) 钉温筛 λ₁(split)−1 @ T=Tc_E(N)   —— §6 实际用的量
  (B) 直接解 Tc_split 与 Tc_E 之比       —— teammate 实际测的量
· N_c = ⌊ω_c/2πT⌋ 逐 N 印出；N_c > N 的行按仓库既有门标 ⚠（库仑项覆盖全矩阵）。
· 大 N 用 Lanczos(eigsh)；**先在 N=400 上对 eigh 逐位核对**，不核对不许用。

严格性：numerical-evidence。
"""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh

KB = 8.617333262e-2
WB = 150.0                       # ω̄₂ [meV]，标度无关


def build(T, G, N, mus, wc):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    L1 = np.zeros(N)
    L1[1:] = 2.0 * np.cumsum(G[1:N]) / (2.0 * n[1:] + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    M = (G[np.abs(i - j)] + G[i + j + 1]) * u[i] * u[j] - np.diag(L1)
    M[n, n] -= G[0] * u * u
    if mus > 0.0:
        Nc = int(np.floor(wc / (2.0 * np.pi * T)))
        if Nc > 0:
            uc = np.where(n < Nc, u, 0.0)
            M -= 2.0 * mus * np.outer(uc, uc)
    return M, u


def lam1(T, G, N, mus, wc, use_sparse=None):
    M, _ = build(T, G, N, mus, wc)
    if use_sparse is None:
        use_sparse = N > 900
    if not use_sparse:
        return float(eigh(M, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])
    # Lanczos：M 谱在 [-c, ~1]，'LA' 取代数最大；加 shift 提升收敛
    v0 = np.ones(N) / np.sqrt(N)
    w = eigsh(M, k=1, which="LA", v0=v0, tol=1e-12, maxiter=20000,
              return_eigenvectors=False)
    return float(w[0])


def Gk(kmax, T, aw):
    x2 = (2.0 * np.pi * T * np.arange(kmax, dtype=float)) ** 2
    G = np.zeros(kmax)
    for a, w in aw:
        G += 2.0 * a * w / (w * w + x2)
    return G


def Tc(aw, N, mus, wc, lo=0.005, hi=1200.0, nbis=50):
    f = lambda T: lam1(T, Gk(2 * N + 2, T, aw), N, mus, wc) - 1.0
    lo, hi = lo * KB, hi * KB
    if not (f(lo) > 0 > f(hi)):
        return np.nan
    for _ in range(nbis):
        m = np.sqrt(lo * hi)
        if f(m) > 0: lo = m
        else: hi = m
    return np.sqrt(lo * hi) / KB


def two_atom(lam, wbar, w1, wB):
    p, q = lam / 2.0, lam * wbar ** 2 / 2.0
    d = wB ** 2 - w1 ** 2
    if d <= 1e-12: return None
    a1 = w1 * (p * wB ** 2 - q) / d
    a2 = wB * (q - p * w1 ** 2) / d
    if a1 < -1e-14 or a2 < -1e-14: return None
    return [(max(a1, 0.0), w1), (max(a2, 0.0), wB)]


MUS, R, LAM = 0.30, 5.0, 0.540
WC = R * WB
EIN = [(LAM * WB / 2.0, WB)]

# ── 门 0：eigsh 对 eigh 逐位核对（不核对不许用大 N）────────────────────
print("=" * 100)
print("门 0：Lanczos(eigsh) vs 稠密 eigh —— 不一致就停")
print("=" * 100)
_N = 400
_T = 5.0 * KB
_G = Gk(2 * _N + 2, _T, EIN)
_d = lam1(_T, _G, _N, MUS, WC, use_sparse=False)
_s = lam1(_T, _G, _N, MUS, WC, use_sparse=True)
print(f"  N={_N}  eigh = {_d:.14f}   eigsh = {_s:.14f}   |差| = {abs(_d-_s):.2e}")
assert abs(_d - _s) < 1e-10, "eigsh 与 eigh 不符，停"
print("  ✓ 通过\n")

# ── 主表：正对照点的 N 收敛 ───────────────────────────────────────────
print("=" * 100)
print(f"正对照 N 收敛：μ*={MUS}, r={R:.0f} (ω_c={WC:.0f} meV), λ={LAM}, ω̄₂={WB:.0f} meV")
print("  (A) 钉温筛 λ₁(split)−1 @ T=Tc_E(N)  —— §6 实际用的量，>0 即「违反」")
print("  (B) Tc_split/Tc_E                   —— teammate 实际测的量，>1 即「违反」")
print("=" * 100)
SPLITS = [(0.95, 100.0), (0.95, 20.0), (0.95, 10.0), (0.8, 5.0), (0.5, 3.0)]
print(f"  {'N':>6}{'Tc_E[K]':>11}{'N_c':>8}{'N_c/N':>7}   "
      + "".join(f"{'(%.2f,%g)' % s:>16}" for s in SPLITS))
rows = {}
for N in (300, 400, 800, 1600, 3200):
    t0 = time.time()
    tE = Tc(EIN, N, MUS, WC)
    if not np.isfinite(tE):
        print(f"  {N:6d}   Einstein 本身不超导（截断地板之下）"); continue
    T = tE * KB
    Nc = int(np.floor(WC / (2.0 * np.pi * T)))
    GE = Gk(2 * N + 2, T, EIN)
    lE = lam1(T, GE, N, MUS, WC)
    cells = []
    for (w1f, w2f) in SPLITS:
        a = two_atom(LAM, WB, w1f * WB, w2f * WB)
        ls = lam1(T, Gk(2 * N + 2, T, a), N, MUS, WC) if a is not None else np.nan
        cells.append(ls - 1.0)
    rows[N] = (tE, Nc, lE, cells)
    print(f"  {N:6d}{tE:11.4f}{Nc:8d}{Nc/N:7.2f}   "
          + "".join(f"{c:+16.3e}" for c in cells)
          + ("  ⚠N_c>N" if Nc > N else "") + f"   [{time.time()-t0:.0f}s]", flush=True)
print(f"\n  （λ₁(Einstein)−1 逐行：" + ", ".join(f"N={N}:{v[2]-1:+.1e}" for N, v in rows.items())
      + "  —— 应 ≈0，二分收敛自检）")

# ── (B) 直接解 Tc 之比（teammate 的量）─────────────────────────────────
print()
print("=" * 100)
print("(B) 直接解 Tc_split / Tc_E（teammate 报「N=300→600 跳 3×」的那个量）")
print("=" * 100)
print(f"  {'N':>6}{'Tc_E[K]':>11}" + "".join(f"{'(%.2f,%g)' % s:>16}" for s in SPLITS))
for N in (300, 400, 800, 1600):
    tE = Tc(EIN, N, MUS, WC)
    if not np.isfinite(tE):
        print(f"  {N:6d}   Einstein 不超导"); continue
    cs = []
    for (w1f, w2f) in SPLITS:
        a = two_atom(LAM, WB, w1f * WB, w2f * WB)
        ts = Tc(a, N, MUS, WC) if a is not None else np.nan
        cs.append(ts / tE if np.isfinite(ts) else np.nan)
    print(f"  {N:6d}{tE:11.4f}" + "".join(f"{c:16.4f}" for c in cs), flush=True)

print()
print("=" * 100)
print("判读规则（落笔前对照）：")
print("  · (A) 列符号随 N 翻负 ⟹ §6 正对照是截断赝物 ⟹ 撤回 §4/§6 的两条结论。")
print("  · (A) 保号且量级稳定 ⟹ 正对照成立；teammate 的 Tc_E 跳变属于 (B) 的二分不收敛，")
print("        **不自动传染到 (A)**（钉温比较里 Einstein 与 split 共用同一 T、同一 N）。")
print("  · (A)(B) 判词相反 ⟹ 两个量测的不是同一件事，必须逐一说明，不许择优报告。")
print("=" * 100)
