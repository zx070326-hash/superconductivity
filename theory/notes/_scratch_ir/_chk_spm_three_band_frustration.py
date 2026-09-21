"""
独立核验「B>=3 全带间排斥时符号阻挫」——PF 论证要求存在一个符号阵 D=diag(+-1,...,+-1)使
D O D 的非对角元全部 >=0。对固定的 3 带、三条带间边全为排斥（Lam^sf_12,Lam^sf_13,Lam^sf_23
全部严格正、Lam^ph 只在带内）：

  Lambda^Delta_ij (i!=j) = -Lambda^sf_ij < 0 对三条边都成立。
  D 共轭后 (i,j) 块乘 D_i D_j；要让三条边同时变号（负变正），需要
     D_1 D_2 > 0,  D_1 D_3 > 0,  D_2 D_3 > 0
  三式相乘：(D_1 D_2 D_3)^2 > 0，这总成立，不排除任何东西——真正的约束是要三个都同号，
  D_i in {+1,-1}，D_1 D_2>0 意味着 D_1=D_2；同理 D_1=D_3, D_2=D_3；三者只能全同号，
  这样一来 D_1 D_2 D_3 全同号仍然满足两两同号 -> 矛盾在哪？直接穷举 8 种 D 更清楚，不空谈。

  这与「符号图是否平衡」的标准结果一致：三角形三条边全负 = 不平衡（奇负环）；
  B=2 只有一条边，环长 1（不存在非平凡环）恒平衡。

  代码里对全部 2^3=8 个 D in {+1,-1}^3，直接算 DOD 的三个带间 2x2 块符号，判断是否非负。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, os.path.dirname(__file__))
import importlib
import numpy as np
import itertools

mod = importlib.import_module('_chk_spm_within_P_construction')
build_O = mod.build_O
einstein = mod.einstein


def three_band_frustrated(N, rng):
    # 三带：Lam^ph 只带内(各自随机)，Lam^sf 三条带间边(12,13,23)全排斥、各自随机正强度
    w = {i: (rng.uniform(0.4, 1.8), rng.uniform(0.3, 2.5)) for i in range(3)}
    wsf = {}
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        wsf[(i, j)] = (rng.uniform(0.3, 1.8), rng.uniform(0.3, 2.5))

    def Lam_ph(i, j, k):
        if i != j:
            return 0.0
        wi, ai = w[i]
        return einstein(k, w=wi, amp=ai)

    def Lam_sf(i, j, k):
        if i == j:
            return 0.0
        key = (min(i, j), max(i, j))
        wk, ak = wsf[key]
        return einstein(k, w=wk, amp=ak)

    return Lam_ph, Lam_sf


def check_all_signings(N, Lam_ph, Lam_sf, B=3, tol=1e-10):
    O = build_O(N, Lam_ph, Lam_sf, B)
    results = []
    for signs in itertools.product([1, -1], repeat=B):
        D = np.diag(np.repeat(signs, N)[np.argsort(np.tile(np.arange(N), B))]) \
            if False else None
        # 构造 D 作用在组合索引 (n,i) 上：band i 对应符号 signs[i]
        Dvec = np.zeros(N * B)
        for n in range(N):
            for i in range(B):
                Dvec[n * B + i] = signs[i]
        Dm = np.diag(Dvec)
        DOD = Dm @ O @ Dm
        # 检查带间块 (i!=j) 的所有元素是否 >=0（Metzler 的跨带部分）
        bad = False
        worst = 0.0
        for n in range(N):
            for m in range(N):
                for i in range(B):
                    for j in range(B):
                        if i == j:
                            continue
                        val = DOD[n * B + i, m * B + j]
                        if val < worst:
                            worst = val
                        if val < -tol:
                            bad = True
        results.append((signs, not bad, worst))
    return results


if __name__ == "__main__":
    rng = np.random.default_rng(31415)
    N = 12  # 小 N 足够判定符号结构，穷举 8 种签名 * 全部 (n,m,i,j) 组合
    Lam_ph, Lam_sf = three_band_frustrated(N, rng)
    results = check_all_signings(N, Lam_ph, Lam_sf, B=3)
    print("=== 三带全排斥：穷举全部 2^3=8 种符号阵 D，检查 D O D 带间块是否处处 >=0 ===")
    any_ok = False
    for signs, ok, worst in results:
        print(f"D={signs}: 带间块 Metzler(>=0)？ {ok}   最劣值={worst:+.6f}")
        any_ok = any_ok or ok
    print(f"\n存在使 Metzler 成立的签名？ {any_ok}  (预期 False，即三带全排斥无法被任何符号阵 Metzler 化)")

    print("\n=== 数值顶本征向量的实际带符号（numerical-evidence，非 PF 强制）===")
    for N in (40, 100, 200):
        O = build_O(N, Lam_ph, Lam_sf, 3)
        w, v = np.linalg.eigh(O)
        top = v[:, -1].reshape(N, 3)
        signs = [np.sign(top[np.argmax(np.abs(top[:, i])), i]) for i in range(3)]
        print(f"N={N}: lambda1={w[-1]:.6f}  band_signs={signs}")

    print("\n=== 对照：三带，只有 12/13 排斥、23 吸引（无阻挫环）应能被 Metzler 化 ===")
    def Lam_sf_unfrustrated(i, j, k):
        if i == j:
            return 0.0
        # (0,1) 与 (0,2) 用负号(排斥, 走 Lam^sf 通道), (1,2) 也走 Lam^sf 但强制符号一致不产生阻挫：
        # 简化: 直接构造成 "星形"耦合: 带0是中心, 与带1/带2都排斥, 带1-带2无耦合(0耦合)
        if (i, j) in [(1, 2), (2, 1)]:
            return 0.0
        return einstein(k, w=1.0, amp=2.0)
    O2 = build_O(60, Lam_ph, Lam_sf_unfrustrated, 3)
    w2, v2 = np.linalg.eigh(O2)
    top2 = v2[:, -1].reshape(60, 3)
    signs2 = [np.sign(top2[np.argmax(np.abs(top2[:, i])), i]) for i in range(3)]
    print(f"星形(无阻挫环) N=60: lambda1={w2[-1]:.6f}  band_signs={signs2}  (带1,带2之间无边，可自由取相同符号)")
