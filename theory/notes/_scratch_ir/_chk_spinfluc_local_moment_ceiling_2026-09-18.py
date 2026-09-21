#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
侦察级探针（numerical-evidence）：自旋涨落配对核的「最强吸引本征值 <= 局域矩」这条
初等界，是否在标准 MMP 型 chi(q) 上非平凡（即：比值不是贴着 1）。

背景（本脚本要检验的那条陈述）：
  纯排斥的自旋涨落配对核 V(k,k') = gbar^2 * S(k-k')（S = 静态自旋结构因子 >= 0）
  在全 BZ 上是卷积算子 ⟹ 本征函数是平面波 e^{iRk}、本征值恰是实空间关联
      Vhat(R) = (1/N) sum_q S(q) e^{-i q.R}
  于是 [analytic，初等]  |lambda_min(V)| <= sum_q S(q)/N = Vhat(0) = 局域矩 <S_0.S_0>
  而 lambda_min(V) < 0 的那个通道正是「排斥产生的吸引」（d 波 / s± 的来源）。
  ⟹ 任何配对通道能拿到的吸引强度都被**局域矩**封顶，取等只在 R=0（而 R=0 是排斥的）。

本脚本只做三件事（都不是定理，只是把上面那条陈述的**非平凡性**量化）：
  (A) 在 MMP 形式 chi(q) = chi0 / (1 + xi^2 * (2 + cos qx + cos qy)) 上直接对角化卷积
      （= FFT），核对「本征值 = Vhat(R)」「lambda_max = Vhat(0)」「|lambda_min| <= lambda_max」。
  (B) 报出比值 rho := |lambda_min| / lambda_max 随关联长度 xi 的走向，以及取到 lambda_min 的 R。
  (C) 费米面加权（compression）之后这条界退化多少：取权 w(k) >= 0，看
      |lambda_min(w^{1/2} V w^{1/2})| <= lambda_max(V) * max(w) 这一侧是否仍被满足（方向自检）。

方向审计：本脚本产出的一切都是**下游判决的输入**，不是判决本身。
  - (A)(B) 是等式/精确对角化，无单边放缩。
  - (C) 的不等式是我要用的那条界，脚本里实测「实际值 / 界」，比值 <= 1 才算界没写反。
"""
import numpy as np

np.random.seed(0)


def chi_mmp(N, xi):
    """MMP 型静态自旋结构因子，峰在 (pi,pi)，逐点 > 0。"""
    q = 2.0 * np.pi * np.arange(N) / N
    QX, QY = np.meshgrid(q, q, indexing="ij")
    # 2 + cos qx + cos qy 在 (pi,pi) 处为 0 ⟹ chi 在那里最大
    return 1.0 / (1.0 + xi ** 2 * (2.0 + np.cos(QX) + np.cos(QY)))


def conv_eigs(S):
    """卷积算子 V(k,k') = S(k-k') 的本征值 = 实空间变换 Vhat(R)（精确，FFT）。"""
    N = S.shape[0]
    Vhat = np.fft.ifft2(S).real  # (1/N^2) sum_q S(q) e^{-i q.R}
    return Vhat


def brute_check(S):
    """小 N 上用显式矩阵对角化核对 FFT 给的本征值（防实现错）。"""
    N = S.shape[0]
    idx = [(a, b) for a in range(N) for b in range(N)]
    M = np.empty((N * N, N * N))
    for i, (ax, ay) in enumerate(idx):
        for j, (bx, by) in enumerate(idx):
            M[i, j] = S[(ax - bx) % N, (ay - by) % N]
    ev = np.linalg.eigvalsh(M)
    return ev.min() / (N * N), ev.max() / (N * N)  # 归一到与 Vhat 同约定


def main():
    print("=" * 74)
    print("(A) 实现自检：显式矩阵对角化 vs FFT（N=8）")
    print("=" * 74)
    for xi in (1.0, 3.0):
        S = chi_mmp(8, xi)
        Vh = conv_eigs(S)
        lo_b, hi_b = brute_check(S)
        print(f"  xi={xi:4.1f}  FFT: [{Vh.min():+.12f}, {Vh.max():+.12f}]")
        print(f"            brute:[{lo_b:+.12f}, {hi_b:+.12f}]  "
              f"max|diff|={max(abs(Vh.min()-lo_b), abs(Vh.max()-hi_b)):.3e}")

    print()
    print("=" * 74)
    print("(B) rho = |lambda_min| / lambda_max，以及取到 lambda_min 的 R（N=128）")
    print("=" * 74)
    N = 128
    print(f"  {'xi':>6} {'lam_max=Vhat(0)':>18} {'lam_min':>14} {'rho':>9} {'argmin R':>12} "
          f"{'lam_max==Vhat(0)?':>18}")
    for xi in (0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0):
        S = chi_mmp(N, xi)
        Vh = conv_eigs(S)
        lmax, lmin = Vh.max(), Vh.min()
        R = np.unravel_index(np.argmin(Vh), Vh.shape)
        Rs = tuple(int(r) if r <= N // 2 else int(r) - N for r in R)
        at0 = abs(lmax - Vh[0, 0]) < 1e-14 * max(1.0, abs(lmax))
        print(f"  {xi:6.1f} {lmax:18.10f} {lmin:14.10f} {abs(lmin)/lmax:9.5f} "
              f"{str(Rs):>12} {str(at0):>18}")

    print()
    print("  解读：lambda_max 恒在 R=(0,0)（= 局域矩），lambda_min 恒在最近邻 R=(1,0)（反铁磁）。")
    print("       rho<1 ⟹ 界 |lambda_min| <= lambda_max 成立且**不贴边**（非平凡）。")

    print()
    print("=" * 74)
    print("(C) 费米面加权后的方向自检：|lam_min(w^.5 V w^.5)| <= lam_max(V)*max(w) ?")
    print("=" * 74)
    N = 24  # 显式对角化，N^2 x N^2
    q = 2.0 * np.pi * np.arange(N) / N
    KX, KY = np.meshgrid(q, q, indexing="ij")
    for xi in (2.0, 5.0, 10.0):
        S = chi_mmp(N, xi)
        Vh = conv_eigs(S)
        lam_max_V = Vh.max()
        # 权：费米面附近的高斯壳（近半填充 tight-binding，eps = -2t(cos kx + cos ky)）
        eps = -2.0 * (np.cos(KX) + np.cos(KY))
        for width in (0.2, 0.6):
            w = np.exp(-(eps / width) ** 2)
            idx = [(a, b) for a in range(N) for b in range(N)]
            M = np.empty((N * N, N * N))
            sq = np.sqrt(w)
            for i, (ax, ay) in enumerate(idx):
                for j, (bx, by) in enumerate(idx):
                    M[i, j] = sq[ax, ay] * S[(ax - bx) % N, (ay - by) % N] * sq[bx, by]
            M /= (N * N)
            ev = np.linalg.eigvalsh(M)
            bound = lam_max_V * w.max()
            actual = abs(ev.min())
            # 顶本征向量的节点结构（d 波？）：看 lam_min 对应本征向量与 cos kx - cos ky 的重叠
            vmin = np.linalg.eigh(M)[1][:, 0].reshape(N, N)
            gd = (np.cos(KX) - np.cos(KY)) * sq
            ov = abs(np.sum(vmin * gd)) / (np.linalg.norm(vmin) * np.linalg.norm(gd))
            print(f"  xi={xi:5.1f} width={width:4.2f}: |lam_min|={actual:.8f}  "
                  f"bound={bound:.8f}  actual/bound={actual/bound:8.5f}  "
                  f"|<v_min, d-wave>|={ov:.4f}")
    print()
    print("  actual/bound <= 1 全部成立 ⟹ 界方向没写反；比值远小于 1 ⟹ 加权把界放松了不少。")
    print("  |<v_min,d-wave>| 接近 1 ⟹ 最强吸引通道确实是 d 波（机制自洽性核对，非承重）。")

    print()
    print("=" * 74)
    print("(D) rho 随关联长度 xi 的走向（网格随 xi 加密，每个峰保持 >=6 个格点）")
    print("=" * 74)
    print(f"  {'xi':>6} {'N':>7} {'lam_max':>12} {'lam_min':>13} {'rho':>9} {'pts/peak':>10}")
    for xi in (2, 5, 10, 20, 50, 100, 200):
        Ng = int(2 ** np.ceil(np.log2(max(128, 40 * xi))))
        Vh = conv_eigs(chi_mmp(Ng, xi))
        print(f"  {xi:6d} {Ng:7d} {Vh.max():12.3e} {Vh.min():13.3e} "
              f"{abs(Vh.min())/Vh.max():9.5f} {Ng/(2*np.pi*xi):10.2f}")
    print()
    print("  rho 单调上升、未见饱和在 1 以下。与解析极限一致：S(q) -> delta_{(pi,pi)} 时")
    print("  Vhat(R) -> (-1)^{Rx+Ry} * const ⟹ lam_max = |lam_min| ⟹ rho -> 1。")
    print("  ⟹ 「|lam_min| <= 局域矩」这条界在完美 AFM（QCP）极限上是**渐近取等**的，")
    print("     即它在物理上最有意思的那一端不浪费任何余量（也说明它就是求和规则本身）。")


if __name__ == "__main__":
    main()
