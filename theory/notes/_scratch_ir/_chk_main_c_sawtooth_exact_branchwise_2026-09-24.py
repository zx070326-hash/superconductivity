# -*- coding: utf-8 -*-
"""
主会话从零第三次独立重算 `c(mu*,r) = max_lam kTc/A`（2026-09-24）。

背景：同一个量现在有三个互不相同的值：
  我（黄金分割）      0.16562437  @ lam=1.6979   —— 已自认前提为假（单峰）
  我（网格+细化）      0.16573686  @ lam~1.70    —— 写进了笔记 §2/§5 与 memory
  leavens-bridge      0.16574437  @ lam=1.600    —— 本轮来报，说是更高的齿
三者都是**下估**（都是「扫到的最大值」），真 sup 只会更大。

━━━ 本脚本的新东西：换参数化，把锯齿变成显式分支，并把内层二分整个消掉 ━━━

判据 `lam1(lam*Khat(w) - 2mu* P(w)) = 1`，其中 P = u_c u_c^T、`Nc = floor(r*w)`。
**引理（analytic，本脚本唯一的解析内容）**：固定 w，记 B := I + 2mu* P（正定），则
    lam1(lam*Khat - 2mu* P) >= 1
  <=> exists v: lam<Khat v,v> - 2mu*<Pv,v> >= <v,v>
  <=> exists v: lam<Khat v,v> >= <Bv,v>                       （移项，B 的定义）
  <=> exists v（含 <Khat v,v> > 0）: lam >= <Bv,v>/<Khat v,v>
  <=> lam >= 1/theta_max,   theta_max := max_v <Khat v,v>/<Bv,v>.
所以**最小根** lam*(w) = 1/theta_max(w)，theta_max 是广义对称本征问题
`Khat v = theta B v` 的最大广义本征值 —— **一次 eigh 精确给出，不需要任何二分**。
（这同时证明了 lam1 关于 lam 的「首次穿越」是良定义的，不会有赝根。）

于是根轨迹用 w 参数化：(lam*(w), w)，而 `Nc=floor(r*w)` 在 w 上是**显式**台阶，
只在 w = k/r 处跳。⟹ 在每个开区间 (k/r,(k+1)/r) 内 c(w) := 1/(pi*lam*(w)*w) 是光滑的，
逐分支取 max 即得每个「齿」的齿顶，不存在「优化器落进哪个吸引域」的运气成分。

这是与 `_chk_main_leavens_c_and_phiA_rescale_2026-09-23.py` 的
（lam 扫 + w 二分）路线**完全不同的第二套算法**，用来交叉验证而非复用。

━━━ 方向审计（承重）━━━
`lam_area = 1/(pi^2 c^2 M2)` ⟹ `d lam_area/dc < 0`。
c 是「扫到的最大值」⟹ **下估** ⟹ lam_area 被**高**估 ⟹ Leavens 的排除力被**高**估
⟹ 比值 lam*/lam_area 被**低**估 ⟹ 对「我们更强」这个判词是**保守**方向。
（leavens-bridge 报的是反过来的方向；本脚本的 §3 把它数值化钉死。）

严格性：引理 analytic；一切数值 numerical-evidence（float64, eigh）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh

M2_BASE = 1.453151466        # eta/M = 9.172, T = 300 K（松原单位）


def pieces(wbar, N):
    """lam=1 归一的 Khat（无库仑）与 u。与 LV.pieces 同式，但本脚本自持一份。"""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    k = np.arange(2 * N + 2)
    G = wbar * wbar / (wbar * wbar + k * k)
    Kp, Km = G[np.abs(i - j)], G[i + j + 1]
    S = (Kp - Km).sum(axis=1)
    return (Kp + Km) * np.outer(u, u) - np.diag(u * u * S), u


def lam_root_exact(wbar, mus, r, N):
    """引理给出的**最小** lam 使 lam1 = 1：lam = 1/theta_max，一次广义 eigh。"""
    Khat, u = pieces(wbar, N)
    B = np.eye(N)
    Nc = int(np.floor(r * wbar))
    if mus > 0.0 and Nc > 0:
        uc = np.where(np.arange(N) < Nc, u, 0.0)
        B = B + 2.0 * mus * np.outer(uc, uc)
    th = float(eigh(Khat, B, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])
    return (1.0 / th) if th > 0 else np.inf, Nc


def c_of_w(wbar, mus, r, N):
    lam, Nc = lam_root_exact(wbar, mus, r, N)
    if not np.isfinite(lam) or lam <= 0:
        return None, None, Nc
    return 1.0 / (np.pi * lam * wbar), lam, Nc


def branch_max(k, mus, r, N, npts, eps=1e-9):
    """第 k 齿：w in (k/r, (k+1)/r) 内取 max c。返回 (c, w, lam, Nc)。"""
    lo, hi = k / r + eps, (k + 1) / r - eps
    ws = np.linspace(lo, hi, npts)
    best = (-1.0, None, None, None)
    for w in ws:
        c, lam, Nc = c_of_w(w, mus, r, N)
        if c is not None and c > best[0]:
            best = (c, w, lam, Nc)
    # 分支内再细化一轮（分支内光滑，可用局部三点加密）
    if best[1] is not None:
        h = (hi - lo) / (npts - 1)
        a, b = max(lo, best[1] - h), min(hi, best[1] + h)
        for w in np.linspace(a, b, 41):
            c, lam, Nc = c_of_w(w, mus, r, N)
            if c is not None and c > best[0]:
                best = (c, w, lam, Nc)
    return best


def main():
    mus, r, N = 0.13, 10.0, 300
    print("=" * 78)
    print("GATE C1  交叉验证：新算法（广义 eigh，零二分） vs 旧算法（w 二分）")
    print("=" * 78)
    import importlib.util as ilu
    spec = ilu.spec_from_file_location(
        "LV", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "_chk_main_leavens_c_and_phiA_rescale_2026-09-23.py"))
    LV = ilu.module_from_spec(spec); spec.loader.exec_module(LV)
    print("  lam        w*(旧:二分)     w*(新:1/theta 反解)   |差|")
    for lam in (1.0, 1.6, 1.69786562, 2.5, 3.2953):
        _, w_old = LV.roots_in_wbar(lambda w: LV.lam1_exact(w, lam, mus, r, N), 1.0, 0.05, 12.0)
        # 新路线：解 lam_root_exact(w)=lam，对 w 二分（w 上 lam*(w) 单调）
        a, b = 0.05, 12.0
        for _ in range(60):
            m = 0.5 * (a + b)
            if lam_root_exact(m, mus, r, N)[0] > lam:
                a = m
            else:
                b = m
        w_new = 0.5 * (a + b)
        print(f"  {lam:9.6f}  {w_old:14.9f}  {w_new:18.9f}   {abs(w_old-w_new):.3e}")

    print()
    print("=" * 78)
    print(f"GATE C2  逐齿 max c（mu*={mus}, r={r}, N={N}）—— 锯齿的显式分支")
    print("=" * 78)
    print("   k   Nc     齿顶 w         齿顶 lam        c(齿顶)")
    rows = []
    for k in range(4, 31):
        c, w, lam, Nc = branch_max(k, mus, r, N, npts=60)
        if c is None:
            continue
        rows.append((c, k, w, lam, Nc))
        print(f"  {k:3d}  {Nc:3d}  {w:12.9f}  {lam:12.8f}  {c:14.10f}")
    rows.sort(reverse=True)
    print()
    print("  ── 最高 5 齿 ──")
    for c, k, w, lam, Nc in rows[:5]:
        print(f"  k={k:3d}  lam={lam:10.7f}  w={w:11.8f}  c={c:.10f}")
    cmax, kmax, wmax, lmax, ncmax = rows[0]

    print()
    print("=" * 78)
    print("GATE C3  三个既有报数的复核（各自 lam 处的 c，用新算法）")
    print("=" * 78)
    for tag, lam in (("我(黄金分割)", 1.69786562), ("leavens-bridge", 1.600),
                     ("我(网格+细化)", 1.70)):
        a, b = 0.05, 12.0
        for _ in range(60):
            m = 0.5 * (a + b)
            if lam_root_exact(m, mus, r, N)[0] > lam:
                a = m
            else:
                b = m
        w = 0.5 * (a + b)
        print(f"  {tag:16s} lam={lam:10.6f}  w*={w:11.8f}  "
              f"Nc={int(np.floor(r*w)):3d}  c={1.0/(np.pi*lam*w):.10f}")
    print(f"  {'本脚本逐齿 max':16s} lam={lmax:10.6f}  w*={wmax:11.8f}  "
          f"Nc={ncmax:3d}  c={cmax:.10f}")

    print()
    print("=" * 78)
    print("GATE C4  N 收敛（齿顶处）")
    print("=" * 78)
    for NN in (200, 300, 400, 600, 800):
        c, w, lam, Nc = branch_max(kmax, mus, r, NN, npts=40)
        print(f"  N={NN:4d}   c={c:.10f}   lam={lam:11.7f}   w={w:11.8f}")

    print()
    print("=" * 78)
    print("GATE C5  方向审计：c -> lam_area 的单调性与三个报数的下游")
    print("=" * 78)
    print(f"  lam_area = 1/(pi^2 c^2 M2),  M2 = {M2_BASE}")
    print("   c 报数          lam_area     比值 3.29/lam_area   来源")
    for tag, cc in (("0.16562437", 0.16562437), ("0.16573686", 0.16573686),
                    ("0.16574437", 0.16574437), (f"{cmax:.8f}", cmax)):
        la = 1.0 / (np.pi ** 2 * cc * cc * M2_BASE)
        print(f"  {tag:14s} {la:11.6f}   {3.29/la:14.4f}      "
              f"{'本脚本' if abs(cc-cmax)<1e-12 else ''}")
    print()
    print("  ⟹ d(lam_area)/dc < 0：c 越大 ⟹ lam_area 越小 ⟹ Leavens 排除力越弱")
    print("  ⟹ 比值 lam*/lam_area 越大 ⟹ **对「我们更强」是有利方向**。")
    print("  ⟹ 所有「扫到的 max」都是 c 的下估 ⟹ 比值被低估 ⟹ 判词保守，方向安全。")
    print("     （leavens-bridge 报「更高的齿 ⟹ lam_area 继续变大 ⟹ 对我们不利」= 方向反了。）")


if __name__ == "__main__":
    main()
