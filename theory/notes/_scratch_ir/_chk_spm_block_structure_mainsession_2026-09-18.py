#!/usr/bin/env python3
"""主会话从零复核 pclass-nogo 的 s± 构造，并把它从 numerical-evidence 升级到 analytic。

我不 import 也不改 pclass-nogo 的脚本。三件事：

(1) 逐字按 spm_gap_needs_no_exotic_physics_2026-09-18.md §3 写下来的公式重建算子，
    与该笔记脚本里**实际写的** if 分支对比 —— 两者不一致（脚本把 u_n²Λ^Δ_ij(0)
    也guard 在 i==j 下，公式里那一项只 guard 在 n==m 下）。

(2) 块结构：Λ^ph 纯带内 + Λ^sf 纯带间 ⟹ 在带空间里
        O = I₂⊗A + σ_x⊗(−B),   A = T[P] − L1[R],  B = T[S]
    ⟹ 本征向量必为 (1,1)⊗v（对称，本征值 λ(A−B)）或 (1,−1)⊗v（反对称，λ(A+B)）。
    若 B ⪰ 0 则 λ₁(A+B) ≥ λ₁(A−B) ⟹ **顶本征向量必定反号**。这是 analytic 的，
    不需要任何数值实验。数值只用来确认 B ⪰ 0 与严格不等号。

(3) B ⪰ 0 是否成立：T[f]_{nm} = u_n u_m[f(|n−m|) + f(n+m+1)] 对 f∈𝒫 是反射延拓
    引理（定理 I 引理 2）。脚本版本的 B 少了对角的 u_n²f(0) ⟹ B_code = T[S] − diag(u_n²S(0))，
    **未必 PSD** ⟹ 脚本恰好把这条 analytic 论证破坏掉了。两个版本都算。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

W_PH, A_PH = 1.0, 1.6      # 带内声子 Einstein 谱
W_SF, A_SF = 0.6, 2.6      # 带间自旋涨落 Einstein 谱


def ein(k, w, amp):
    k = np.asarray(k, dtype=float)
    return amp * w * w / (w * w + k * k)


def P(k):  return ein(k, W_PH, A_PH)
def S(k):  return ein(k, W_SF, A_SF)


def T(fun, N, include_f0_diag=True):
    """T[f]_{nm} = u_n u_m [f(|n−m|) + f(n+m+1)]，n=m 处 f(|n-m|)=f(0)。"""
    n = np.arange(N)
    u = (2.0 * n + 1.0) ** -0.5
    d = np.abs(n[:, None] - n[None, :])
    s = n[:, None] + n[None, :] + 1
    M = np.outer(u, u) * (fun(d) + fun(s))
    if not include_f0_diag:                       # 复现脚本里 i≠j 分支的遗漏
        M = M - np.diag(u * u * float(fun(0)))
    return M


def L1(N):
    """带内的 Z 行和项 u_n²(R(0) + 2Σ_{k=1}^n R(k))，R = Λ^ph + Λ^sf（两带行和相同）。"""
    n = np.arange(N)
    u = (2.0 * n + 1.0) ** -0.5
    R = lambda k: P(k) + S(k)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * R(np.arange(1, N + 1)))))
    return np.diag(u * u * (float(R(0)) + cum[:N]))


def report(N, faithful):
    A = T(P, N) - L1(N)                       # 带内块（含 u_n²P(0)，i=j 分支两版本一致）
    B = T(S, N, include_f0_diag=faithful)     # 带间块（faithful=公式版；False=脚本版）
    lam_anti = np.linalg.eigvalsh(A + B)[-1]  # (1,−1) 扇区
    lam_sym = np.linalg.eigvalsh(A - B)[-1]   # (1,+1) 扇区
    bmin = np.linalg.eigvalsh(B)[0]
    return lam_anti, lam_sym, bmin


print("=== 块分解：λ₁ = max(λ₁(A+B) [反号], λ₁(A−B) [同号]) ===")
print(f"{'N':>5} | {'版本':<10}{'λ₁(A+B)反号':>14}{'λ₁(A−B)同号':>14}{'λmin(B)':>11}  顶向量")
for N in (40, 60, 90, 120, 200, 400):
    for faithful, tag in ((True, "笔记公式"), (False, "脚本实现")):
        la, ls, bmin = report(N, faithful)
        top = "反号(s±)" if la > ls else "同号(s++)"
        print(f"{N:>5} | {tag:<10}{la:14.6f}{ls:14.6f}{bmin:11.3e}  {top}")

print("\n=== 直接建全算子做交叉验证（不用块分解，2N×2N，笔记公式版） ===")
def build_full(N, faithful):
    n = np.arange(N)
    u = (2.0 * n + 1.0) ** -0.5
    dim = 2 * N
    O = np.zeros((dim, dim))
    LamD = lambda i, j, k: (P(k) if i == j else -S(k))
    R = lambda k: P(k) + S(k)
    for a in range(N):
        for b in range(N):
            for i in range(2):
                for j in range(2):
                    v = 0.0
                    if a != b:
                        v += u[a] * u[b] * LamD(i, j, abs(a - b))
                    v += u[a] * u[b] * LamD(i, j, a + b + 1)
                    if a == b:
                        if faithful or i == j:            # ← 争议点就在这一行
                            v += u[a] ** 2 * LamD(i, j, 0)
                        if i == j:
                            s = float(R(0)) + 2.0 * float(np.sum(R(np.arange(1, a + 1))))
                            v -= u[a] ** 2 * s
                    O[a * 2 + i, b * 2 + j] = v
    return O

for faithful, tag in ((True, "笔记公式"), (False, "脚本实现")):
    N = 120
    O = build_full(N, faithful)
    assert np.max(np.abs(O - O.T)) < 1e-12
    w, v = np.linalg.eigh(O)
    top = v[:, -1].reshape(N, 2)
    sg = [np.sign(top[np.argmax(np.abs(top[:, i])), i]) for i in range(2)]
    la, ls, _ = report(N, faithful)
    print(f"{tag}: λ₁={w[-1]:.6f}  块分解 max={max(la, ls):.6f}  "
          f"符合={abs(w[-1]-max(la,ls))<1e-10}  两带符号={sg}")

print("\n=== 引理 2（反射延拓）：T[f] ⪰ 0 对 Einstein f，以及脚本版是否破坏它 ===")
for N in (60, 200, 400):
    print(f"  N={N}: λmin(T[S])={np.linalg.eigvalsh(T(S,N))[0]:+.4e}   "
          f"λmin(T[S]−diag u²S(0))={np.linalg.eigvalsh(T(S,N,False))[0]:+.4e}")

print("""
判读：
  · T[S] ⪰ 0（反射延拓引理，f∈𝒫）⟹ λ₁(A+B) ≥ λ₁(A−B) ⟹ 顶本征向量**必定**两带反号。
    这是 analytic 的一步，不依赖任何数值实验 —— 构造不是"造出来一个例子"，
    是"在这套约定下反号是强制的"。
  · 脚本版把 B 减掉一个正对角，λmin(B)<0 ⟹ 上面那步不再成立，
    结论退回成"这组参数碰巧反号"（numerical-evidence）。
  · 两版的 λ₁ 数值不同 ⟹ 笔记里印的 0.148520 是脚本版的值，需订正。
""")
