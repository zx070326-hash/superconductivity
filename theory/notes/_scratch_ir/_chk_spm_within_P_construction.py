"""
直接构造：两片各自逐元非负（"广义 P 类"）的声子/自旋涨落谱，在 K1 约定
(Lambda^Z = Lambda^ph+Lambda^sf, Lambda^Delta = Lambda^ph-Lambda^sf) 下，
顶本征向量在带指标上确实反号（s±）。

用途：team-lead 派活「造出变号解才算证伪推论」。这里独立从零构造（不 import 仓库模块）,
验证 s± 型带间符号结构可以由纯正的(在带指标各自意义上属于 P 类)谱函数产生,
不需要任何"P 之外"的物理输入——变号来自 Z/Delta 两通道对同一段正谱权取不同线性组合(加号/减号)
这个结构性事实，而不是谱本身带符号。

算子形式取自 theory/notes/ceiling/multiband_ceiling_lift_2026-09-18.md §1（我已独立复核并从零重写）：

  O_{(n,i),(m,j)} = u_n u_m [ Lam_ij(|n-m|) 1_{n!=m} + Lam_ij(n+m+1) ]
                    + delta_{nm} [ u_n^2 Lam_ij(0) - delta_ij u_n^2 ( R_i(0) + 2 sum_{k=1}^n R_i(k) ) ]

其中 Delta 通道用 Lam=Lambda^Delta, Z 通道的行和用 R_i = sum_j (Lambda^Z)_ij(k)。
"""
import numpy as np

def u(n):
    return (2*n+1)**-0.5

def einstein(k, w, amp):
    # KAY 型 Einstein 谱核: f(k) = amp * w^2/(w^2+k^2), k>=0 (含 k=0)
    return amp * w*w/(w*w+k*k)

def build_O(N, Lam_ph, Lam_sf, B):
    """Lam_ph, Lam_sf: 函数 (i,j,k)->float, 各自逐元非负（每个 (i,j) 对独立 Einstein 谱）"""
    def LamZ(i,j,k):
        return Lam_ph(i,j,k) + Lam_sf(i,j,k)
    def LamD(i,j,k):
        return Lam_ph(i,j,k) - Lam_sf(i,j,k)

    # R_i(k) = sum_j LamZ(i,j,k)
    def R(i,k):
        return sum(LamZ(i,j,k) for j in range(B))

    dim = N*B
    O = np.zeros((dim,dim))
    idx = lambda n,i: n*B+i
    for n in range(N):
        for m in range(N):
            un, um = u(n), u(m)
            for i in range(B):
                for j in range(B):
                    val = 0.0
                    if n != m:
                        val += un*um*LamD(i,j,abs(n-m))
                    val += un*um*LamD(i,j,n+m+1)
                    if n==m:
                        # BUG FIX (team-lead caught 2026-09-18): this term must apply for
                        # ALL band pairs i,j at n=m, not just i==j -- it comes from the
                        # 1_{n!=m} guard on the FIRST bracket only, per the formula in
                        # multiband_ceiling_lift_2026-09-18.md sec 1 and this note's sec 3:
                        #   O = u_n u_m[Lam(|n-m|)1_{n!=m} + Lam(n+m+1)]
                        #       + delta_nm[u_n^2 Lam_ij(0) - delta_ij u_n^2(R_i(0)+2 sum R_i(k))]
                        # Original code wrongly gated u_n^2*Lam_ij(0) on i==j too, dropping
                        # the off-diagonal-band n=m contribution -- this changed lambda1
                        # from the correct 2.081 down to 0.1485 for the s+- construction below.
                        val += un*un*LamD(i,j,0)
                        if i==j:
                            s = R(i,0) + 2*sum(R(i,k) for k in range(1,n+1))
                            val -= un*un*s
                    O[idx(n,i), idx(m,j)] = val
    return O

def run_case(N, B, Lam_ph, Lam_sf, label, elemwise_check):
    O = build_O(N, Lam_ph, Lam_sf, B)
    assert np.max(np.abs(O-O.T)) < 1e-10, "算子应对称（细致平衡=1 简化，w_i=1）"
    w, v = np.linalg.eigh(O)
    top = v[:, -1]
    # reshape 成 (N, B)
    top = top.reshape(N, B)
    # 每条带上的符号：用该带上 |分量| 最大的几个 n 的符号（去掉小分量噪声）
    band_sign = []
    for i in range(B):
        col = top[:, i]
        j = np.argmax(np.abs(col))
        band_sign.append(np.sign(col[j]))
    print(f"[{label}] N={N} B={B} lambda1={w[-1]:.6f} lambda2={w[-2]:.6f} band_signs={band_sign}")
    # 逐元非负检查（elemwise_check 抽样若干 (i,j,k)）
    ok = True
    for (i,j,k) in elemwise_check:
        if Lam_ph(i,j,k) < -1e-12 or Lam_sf(i,j,k) < -1e-12:
            ok = False
    print(f"   逐片逐元非负抽样检查: {'PASS' if ok else 'FAIL'}")
    return w[-1], w[-2], band_sign

if __name__ == "__main__":
    N = 60
    # 带内声子: 只有对角非零, 中等吸引
    def Lam_ph(i,j,k):
        if i != j:
            return 0.0
        return einstein(k, w=1.0, amp=1.6)
    # 带间自旋涨落交换: 只有非对角非零(标准 s± 机制: Cooper 对通过带间自旋涨落交换散射)
    def Lam_sf(i,j,k):
        if i == j:
            return 0.0
        return einstein(k, w=0.6, amp=2.6)

    checks = [(i,j,k) for i in range(2) for j in range(2) for k in [0,1,2,5]]
    l1, l2, signs = run_case(N, 2, Lam_ph, Lam_sf, "s±构造(K1约定)", checks)

    assert signs[0] != 0 and signs[1] != 0
    if signs[0] != signs[1]:
        print("\n=> 顶本征向量在两带上符号相反：s± 型带间反号能隙，在 Lambda^ph, Lambda^sf 均逐元非负（各自属 P 类）下产生。")
        print("=> 变号完全来自 K1 约定的 Z=加号/Delta=减号 结构性线性组合，不需要任何 P 之外的谱密度。")
    else:
        print("\n=> 本组参数未产生带间反号（需要调大 Lam_sf/Lam_ph 比例）")

    # 对照：把 Lam_sf 关掉 -> 纯声子, 应回到同号(s++,单带极值型)
    def Lam_sf0(i,j,k):
        return 0.0
    l1b, l2b, signs0 = run_case(N, 2, Lam_ph, Lam_sf0, "对照:纯声子(sf=0)", checks)
    print(f"对照顶带符号: {signs0} (期望同号)")

    # N 收敛检查（承重结论需要）
    for Ncheck in [40, 60, 90, 120]:
        l1c, l2c, sc = run_case(Ncheck, 2, Lam_ph, Lam_sf, "N收敛", checks)
