"""
独立核验 team-lead 指出的 §3.2 缺陷：PF 论证的前提列表漏了「不可约性」，
只写「带内非负+带间非负」两条不足以保证 band-sign-changing。

不抄 team-lead 给的具体参数，自己另选一组数值重构两个反例：

反例 A（Lam^sf ≡ 0，满足"两条前提"，但 O 块对角，顶向量完全局域化在一条带上，不是变号）：
  两带故意不对称：带0 Einstein(w=0.8,amp=2.2)，带1 Einstein(w=1.7,amp=0.9)。
  预期：lambda1 != lambda2（无对称保护的简并），且顶向量另一带分量恰为 0，不可约=False。

反例 B（Lam^ph ≡ 0，Lam^sf 只在单个 k 值非零，Metzler 仍成立但图不连通）：
  Lam^sf(k) = amp0 if k==7 else 0（换一个 k 值，不用 team-lead 的 k=5）。
  预期：DOD 非对角最小元 >= 0（Metzler 仍成立），但图按 |n-m|=7 或 n+m+1=7 分裂成若干弱连通片，
  不可约=False（用 scipy 或手写 BFS 连通分量检查 > 1 个分量）。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, os.path.dirname(__file__))
import importlib
import numpy as np
from collections import deque

mod = importlib.import_module('_chk_spm_within_P_construction')
build_O = mod.build_O
einstein = mod.einstein


def connected_components(adj_nonzero, dim):
    """adj_nonzero: dim x dim bool 邻接（对称，非零即视为边）。返回连通分量个数。"""
    visited = np.zeros(dim, dtype=bool)
    ncomp = 0
    for s in range(dim):
        if visited[s]:
            continue
        ncomp += 1
        q = deque([s])
        visited[s] = True
        while q:
            u = q.popleft()
            nbrs = np.nonzero(adj_nonzero[u])[0]
            for v in nbrs:
                if not visited[v]:
                    visited[v] = True
                    q.append(v)
    return ncomp


def counterexample_A(N=60):
    print("=== 反例 A: Lam^sf ≡ 0（满足两条前提），块对角，顶向量局域化在一条带 ===")

    def Lam_ph(i, j, k):
        if i != j:
            return 0.0
        if i == 0:
            return einstein(k, w=0.8, amp=2.2)
        else:
            return einstein(k, w=1.7, amp=0.9)

    def Lam_sf_zero(i, j, k):
        return 0.0

    O = build_O(N, Lam_ph, Lam_sf_zero, 2)
    assert np.max(np.abs(O - O.T)) < 1e-10

    off = O - np.diag(np.diag(O))
    adj = np.abs(off) > 1e-12
    ncomp = connected_components(adj, 2 * N)

    w, v = np.linalg.eigh(O)
    lam1, lam2 = w[-1], w[-2]
    gap = lam1 - lam2
    degenerate = gap < 1e-8
    top = v[:, -1].reshape(N, 2)
    band0_norm = np.linalg.norm(top[:, 0])
    band1_norm = np.linalg.norm(top[:, 1])
    print(f"N={N}: lambda1={lam1:.6f}  lambda2={lam2:.6f}  gap={gap:.6e}  简并?={degenerate}")
    print(f"顶向量在带0上的范数={band0_norm:.6e}  带1上的范数={band1_norm:.6e}")
    print(f"连通分量数={ncomp}（>1 ⟹ 不可约=False）")
    localized = (band0_norm < 1e-8) or (band1_norm < 1e-8)
    print(f"顶向量完全局域化在一条带（另一带恰为0）？{localized}")
    if degenerate:
        print("⚠ gap≈0：顶本征值简并，任何打印出的 band_signs 都是 LAPACK 在简并子空间里的")
        print("  规范选取，不携带信息——判据只能用 gap 和 ncomp，不能读符号。")
        print(f"⟹ 满足两条前提但结论「顶本征值单重」直接为假（gap={gap:.2e}），"
              f"「顶本征向量恒变号」在简并子空间上无定义：True\n")
    else:
        print(f"⟹ 满足两条前提但结论（band-sign-changing）不成立（单带局域化，非简并）："
              f"{localized and ncomp > 1}\n")
    return localized, degenerate, ncomp


def counterexample_B(N=30, kstar=7):
    print(f"=== 反例 B: Lam^ph ≡ 0, Lam^sf(k)=amp·1_{{k={kstar}}}（Metzler仍成立但图不连通）===")

    def Lam_ph_zero(i, j, k):
        return 0.0

    def Lam_sf(i, j, k):
        if i == j:
            return 0.0
        return 3.5 if k == kstar else 0.0

    O = build_O(N, Lam_ph_zero, Lam_sf, 2)
    assert np.max(np.abs(O - O.T)) < 1e-10

    D = np.ones(2 * N)
    for n in range(N):
        D[n * 2 + 1] = -1.0
    Dm = np.diag(D)
    DOD = Dm @ O @ Dm
    off = DOD - np.diag(np.diag(DOD))
    metzler_min = np.min(off)

    adj = np.abs(O - np.diag(np.diag(O))) > 1e-12
    ncomp = connected_components(adj, 2 * N)

    w, v = np.linalg.eigh(O)
    lam1, lam2 = w[-1], w[-2]
    gap = lam1 - lam2
    print(f"N={N}, k*={kstar}: DOD 非对角最小元={metzler_min:+.6e}（Metzler 仍成立？{metzler_min > -1e-10}）")
    print(f"  lambda1={lam1:.6f}  lambda2={lam2:.6f}  gap={gap:.6e}  简并?={gap<1e-8}")
    if gap < 1e-8:
        print("  警告：gap≈0，顶本征值简并；脚本若打出「反号=True」同样是 LAPACK 规范假象，"
              "不能当作「结论仍成立」的证据——本反例的作用只是证明「Lam^sf≢0」单独不足以给出不可约。")
    print(f"连通分量数={ncomp}（>1 ⟹ 不可约=False，Metzler 不足以保证结论）")
    return metzler_min, ncomp


def check_H3prime(N=50):
    """
    独立核验 team-lead 订正后的充分条件 (H3'):
        Lam^ph_ii(k) > 0 对所有 k（两带各自，Einstein/任何 P 元素自动满足）
        且 Lam^sf_ij 不恒为零（支撑可以是单个 k，振幅可以任意小）
    ⟹ 带内块处处正 => 每条带自身是完全图；一条非零带间边把两个完全图接通 => 全图连通 => 不可约。
    自己另选参数（与本文件其它函数的参数不同，且与 team-lead 消息里的数字不同）：
    带0 Einstein(w=1.2,amp=1.8)，带1 Einstein(w=0.5,amp=3.0)；Lam^sf 只在 k=9 非零，
    振幅分别测 3.0 / 1e-6 / 1e-9 三档，验证"振幅多小都照样强制不可约"。
    """
    print("=== (H3') 核验：Lam^ph 两带处处>0 + Lam^sf 支撑单点(振幅任意小) ⟹ 不可约 ===")

    def Lam_ph(i, j, k):
        if i != j:
            return 0.0
        w_, a_ = (1.2, 1.8) if i == 0 else (0.5, 3.0)
        return einstein(k, w=w_, amp=a_)  # > 0 对所有 k（Einstein 核处处正）

    for amp in (3.0, 1e-6, 1e-9):
        def Lam_sf(i, j, k, amp=amp):
            if i == j:
                return 0.0
            return amp if k == 9 else 0.0

        O = build_O(N, Lam_ph, Lam_sf, 2)
        assert np.max(np.abs(O - O.T)) < 1e-10

        adj = np.abs(O - np.diag(np.diag(O))) > 1e-14
        ncomp = connected_components(adj, 2 * N)

        D = np.ones(2 * N)
        for n in range(N):
            D[n * 2 + 1] = -1.0
        DOD = np.diag(D) @ O @ np.diag(D)
        metzler_min = np.min(DOD - np.diag(np.diag(DOD)))

        w, v = np.linalg.eigh(O)
        lam1, lam2 = w[-1], w[-2]
        gap = lam1 - lam2
        top = v[:, -1].reshape(N, 2)
        j0 = np.argmax(np.abs(top[:, 0]))
        j1 = np.argmax(np.abs(top[:, 1]))
        anti_sign = np.sign(top[j0, 0]) != np.sign(top[j1, 1])
        print(f"  Lam^sf 振幅={amp:.0e}: 连通分量数={ncomp}  Metzler={metzler_min > -1e-10}  "
              f"lambda1={lam1:.6f}  gap(lam1-lam2)={gap:.6e}  单重?={gap>1e-8}  "
              f"顶向量两带反号?={anti_sign}")
        assert ncomp == 1 and gap > 1e-8 and anti_sign, "H3' 应强制不可约+单重+反号"
    print("  ⟹ (H3') 三档振幅全部：不可约、顶本征值单重、顶本征向量严格反号。\n")


if __name__ == "__main__":
    counterexample_A()
    counterexample_B()

    print("=== 对照：把两个反例的谱改成处处>0（比如 Einstein，属𝒫）时不可约自动恢复 ===")
    N = 40

    def Lam_ph(i, j, k):
        if i != j:
            return 0.0
        w_, a_ = (0.8, 2.2) if i == 0 else (1.7, 0.9)
        return einstein(k, w=w_, amp=a_)

    def Lam_sf_P(i, j, k):
        if i == j:
            return 0.0
        return einstein(k, w=1.1, amp=1.4)  # 处处 >0，属𝒫

    O = build_O(N, Lam_ph, Lam_sf_P, 2)
    off = O - np.diag(np.diag(O))
    adj = np.abs(off) > 1e-12
    ncomp = connected_components(adj, 2 * N)
    w, v = np.linalg.eigh(O)
    top = v[:, -1].reshape(N, 2)
    signs = [np.sign(top[np.argmax(np.abs(top[:, i])), i]) for i in range(2)]
    print(f"Lam^sf∈𝒫（处处>0）: 连通分量数={ncomp}（应为1，不可约恢复）  band_signs={signs}  "
          f"lambda1={w[-1]:.6f}")

    check_H3prime()

