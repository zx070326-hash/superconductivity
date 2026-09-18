"""
独立从零推导/核验 Perron-Frobenius 版 s± 论证（team-lead 提出的升级：去掉两带对称假设、
去掉对定理 I 引理 1+2 的依赖）。不 import team-lead 的任何脚本，只复用本目录下已经修好的
build_O（来自 _chk_spm_within_P_construction.py，这是我自己写并修过 bug 的脚本，不是转抄）。

论证要点（自己重新推一遍，不抄团队消息里的措辞）：

  设 B=2。构造要求只有：
    - Lambda^ph_ij(k) 只在 i=j 时非零，且 >=0（带内声子，任意强度，两带不必相同）；
    - Lambda^sf_ij(k) 只在 i!=j 时非零，且 >=0（带间自旋涨落，两带不必对称）。
  则 Lambda^Delta_ij = Lambda^ph_ij - Lambda^sf_ij：
    - i=j 块：等于 Lambda^ph_ii >= 0（带内，Lambda^sf 贡献为 0）；
    - i!=j 块：等于 -Lambda^sf_ij <= 0（带间，Lambda^ph 贡献为 0）。

  算子公式：
    O_{(n,i),(m,j)} = u_n u_m[Lam^Delta_ij(|n-m|) 1_{n!=m} + Lam^Delta_ij(n+m+1)]
                      + delta_nm[u_n^2 Lam^Delta_ij(0) - delta_ij u_n^2 R_i(n)]
  R_i(n) 只影响对角元 (n,i)=(n,i)，不影响非对角结构。

  ⟹ O 的「组合指标非对角元」(n,i)!=(m,j)（含 n=m,i!=j）符号只由 i,j 是否相等决定：
     - i=j, n!=m: 系数 u_n u_m Lambda^ph_ii(...) >= 0
     - i!=j (任意 n,m，含 n=m): 系数 u_n u_m * (-Lambda^sf_ij(...)) <= 0

  取 D = diag(+1 on band1, -1 on band2)（正交对角阵，D^2=I）。共轭 D O D：
     - i=j 块不变（D_i D_j = 1）；
     - i!=j 块乘 D_1 D_2 = -1，符号翻转，变成 >= 0。
  ⟹ DOD 的所有非对角元（组合指标意义下）>= 0，即 Metzler；对角元任意（Metzler 定义不管对角）。
  取 c 充分大使 DOD + cI 的对角也 >= 0（只需 c >= -min_diag）⟹ DOD+cI 逐元非负。

  不可约性：只要 Einstein 谱在所用 k 上处处 >0（本构造如此），则任意 (n,i) 与 (m,j) 之间
  都有非零边（同带靠 Lambda^ph_ii(|n-m|)>0，跨带靠 Lambda^sf_ij(|n-m|)>0）⟹ 对应有向图完全连通
  ⟹ 不可约（其实还更强：矩阵处处严格正，是 primitive 的，不需要不可约这么弱的条件）。

  两步得到「顶本征向量 band-sign-changing 且顶本征值单重」：
  (i) lambda_max(DOD+cI) = rho(DOD+cI)（谱半径）。这一步不用完整 PF 定理，只用一个初等事实：
      对任意逐元非负对称阵 M，任意实向量 x 满足 x^T M x <= |x|^T M |x|（因为每一项 M_ij x_i x_j
      <= M_ij |x_i x_j|，用了 M_ij>=0），所以把 Rayleigh 商的最大化限制在非负象限不减小最大值
      ⟹ lambda_max(M) 有一个非负的最大化向量 ⟹ lambda_max(M) = max_{x>=0,|x|=1} x^T M x。
      而对逐元非负矩阵，谱半径 rho(M) 也等于这个非负象限上的最大值（标准 Perron-Frobenius/
      Collatz-Wielandt），所以 lambda_max(M) = rho(M)。
  (ii) M=DOD+cI 不可约（甚至处处正）非负 ⟹ 标准 Perron 定理：rho(M) 是单重本征值，
       对应本征向量可取**严格正**。
  ⟹ lambda_max(O+cI) = lambda_max(DOD+cI)（D 正交相似不变谱），本征向量 w=D v，
     v 是 DOD+cI 的顶本征向量（严格正）⟹ O 的顶本征向量 = D w，第一带分量 = +w>0，
     第二带分量 = -w<0 ⟹ **band-sign-changing，且单重**。

  以下代码：(A) 随机不对称扫描验证结论；(B) 独立验证 Metzler/非负化步骤；
  (C) 验证 lambda_max(O)=lambda_max(DOD+cI)-c 精确成立。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, os.path.dirname(__file__))
import importlib
import numpy as np

mod = importlib.import_module('_chk_spm_within_P_construction')
build_O = mod.build_O
einstein = mod.einstein


def random_two_band_case(rng, N):
    # 两带完全不对称：带内声子(w,amp)各自独立随机，带间自旋涨落一个(w,amp)（对称核 Lambda_ij=Lambda_ji）
    w1, a1 = rng.uniform(0.3, 2.0), rng.uniform(0.2, 3.0)
    w2, a2 = rng.uniform(0.3, 2.0), rng.uniform(0.2, 3.0)
    wsf, asf = rng.uniform(0.2, 2.0), rng.uniform(0.2, 3.0)

    def Lam_ph(i, j, k):
        if i != j:
            return 0.0
        w, a = (w1, a1) if i == 0 else (w2, a2)
        return einstein(k, w=w, amp=a)

    def Lam_sf(i, j, k):
        if i == j:
            return 0.0
        return einstein(k, w=wsf, amp=asf)

    return Lam_ph, Lam_sf, (w1, a1, w2, a2, wsf, asf)


def check_case(N, Lam_ph, Lam_sf, B=2):
    O = build_O(N, Lam_ph, Lam_sf, B)
    assert np.max(np.abs(O - O.T)) < 1e-9

    D = np.ones(N * B)
    for n in range(N):
        D[n * B + 1] = -1.0  # band index 1 (second band) -> -1
    D = np.diag(D)
    DOD = D @ O @ D
    assert np.max(np.abs(DOD - DOD.T)) < 1e-9

    # Metzler check: off-(combined)-diagonal entries of DOD >= -tol
    off = DOD - np.diag(np.diag(DOD))
    metzler_min = np.min(off)

    c = -np.min(np.diag(DOD)) + 1.0
    M = DOD + c * np.eye(N * B)
    entrywise_min = np.min(M)

    w_O, v_O = np.linalg.eigh(O)
    w_M, v_M = np.linalg.eigh(M)
    lam_O = w_O[-1]
    lam_M = w_M[-1]
    match = abs(lam_O - (lam_M - c)) < 1e-8

    top_O = v_O[:, -1].reshape(N, B)
    signs = [np.sign(top_O[np.argmax(np.abs(top_O[:, i])), i]) for i in range(B)]
    anti = (signs[0] != signs[1])

    top_M = v_M[:, -1].reshape(N, B)
    m_all_pos = np.min(top_M) > -1e-9 or np.max(top_M) < 1e-9  # sign-definite (strict Perron vector)

    return dict(metzler_min=metzler_min, entrywise_min=entrywise_min,
                lam_O=lam_O, lam_M_shifted=lam_M - c, match=match,
                signs=signs, anti=anti, M_positive_vec=m_all_pos)


if __name__ == "__main__":
    print("=== (A) 随机不对称两带扫描（Lam_ph 两带独立随机、Lam_sf 单独随机）===")
    rng = np.random.default_rng(20260918)
    N = 100
    n_trials = 12
    all_anti = True
    all_metzler = True
    all_match = True
    for t in range(n_trials):
        Lam_ph, Lam_sf, params = random_two_band_case(rng, N)
        r = check_case(N, Lam_ph, Lam_sf, B=2)
        all_anti = all_anti and r['anti']
        all_metzler = all_metzler and (r['metzler_min'] > -1e-10)
        all_match = all_match and r['match']
        print(f"trial {t:2d} params={tuple(round(p,3) for p in params)} "
              f"metzler_min={r['metzler_min']:+.3e} lam_O={r['lam_O']:.6f} "
              f"lam_M-c={r['lam_M_shifted']:.6f} match={r['match']} signs={r['signs']} anti={r['anti']} "
              f"PF正向量={r['M_positive_vec']}")
    print(f"\n汇总: 全部 Metzler 成立={all_metzler}  全部 lambda_max(O)==lambda_max(DOD+cI)-c={all_match}  "
          f"全部反号={all_anti}  ({n_trials}/{n_trials})")

    print("\n=== (B) N 收敛（单组随机参数）===")
    rng2 = np.random.default_rng(7)
    Lam_ph, Lam_sf, params = random_two_band_case(rng2, 400)
    for N in (40, 100, 200, 400):
        r = check_case(N, Lam_ph, Lam_sf, B=2)
        print(f"N={N:4d}: lam_O={r['lam_O']:.8f}  signs={r['signs']}  anti={r['anti']}")
