"""
任务：team-lead 派活 2026-09-22 —— 能否把 Lambda^SF 走 PSD 减法而非逐元折进 (K1)，
从而绕开 Eq.(37) 的逐元正断言。

核心恒等式（来自 multiband_ceiling_lift_2026-09-18.md (ML1)(ML2)(ML3-abs)）：
    L := diag(R^Z) - Lambda_tilde^Delta ,  R^Z_i = (Lambda_tilde^Z w)_i / w_i
天花板证明需要 L >= 0 (PSD)。
Lambda^Z = Lambda^ph + Lambda^sf, Lambda^Delta = Lambda^ph - Lambda^sf,
Lambda^ph >= 0 逐元（声子谱函数非负，无争议）。

拆分：
    L = L^ph + M ,   L^ph := diag(R^ph) - Lambda_tilde^ph  (PSD 无条件，Lambda^ph>=0 逐元的标准图 Laplacian)
                      M    := diag(R^sf) + Lambda_tilde^sf  ("signless" 型)

问题：若只知道 Lambda_tilde^sf 是 PSD（带空间 Gram 矩阵，允许负非对角元），
M 是否必然 PSD？

若"是"：SF 片可以只用 PSD（涨落-耗散定理，Team B 已建立的 B.1 结论）代入，不需要 Eq.(37) 的逐元正。
若"否"：给出显式反例，说明 PSD 路由在这一步真正断开。
"""
import numpy as np

def L_signless(A, w):
    """M = diag(R) + A,  R_i = (A w)_i / w_i"""
    B = len(w)
    R = (A @ w) / w
    M = np.diag(R) + A
    return M

def min_eig(M):
    return np.linalg.eigvalsh(M).min()

# ---- 手推反例：B=2, w=(1,1) ----
w2 = np.array([1.0, 1.0])
A_ex = np.array([[4.0, -1.9], [-1.9, 1.0]])
# 确认 A_ex 本身 PSD
eig_A = np.linalg.eigvalsh(A_ex)
print("A_ex (Lambda_tilde^sf) eigenvalues:", eig_A, "-> PSD" if eig_A.min() >= -1e-12 else "-> NOT PSD")

M_ex = L_signless(A_ex, w2)
print("M = diag(R^sf) + Lambda_tilde^sf =\n", M_ex)
eig_M = np.linalg.eigvalsh(M_ex)
print("M eigenvalues:", eig_M, "-> PSD" if eig_M.min() >= -1e-12 else "-> NOT PSD (COUNTEREXAMPLE)")

print()
print("=" * 70)
print("随机扫描：B=2..6, 随机 PSD Lambda_tilde^sf (含负非对角), 随机正 w")
print("=" * 70)

rng = np.random.default_rng(20260922)
n_trials = 20000
n_fail = 0
n_psd_input = 0
worst_min_eig = 0.0
worst_case = None

for trial in range(n_trials):
    B = rng.integers(2, 7)
    w = rng.uniform(0.1, 5.0, size=B)
    # 造一个保证 PSD 的对称矩阵：A = S S^T，S 随机（允许负号 => 非对角可正可负）
    S = rng.normal(size=(B, B)) * rng.uniform(0.2, 3.0)
    A = S @ S.T
    # 确认 A 本身 PSD（构造保证，只做 sanity check）
    eA = np.linalg.eigvalsh(A)
    if eA.min() < -1e-8:
        continue  # 数值噪声，跳过
    n_psd_input += 1

    M = L_signless(A, w)
    eM = np.linalg.eigvalsh(M)
    m = eM.min()
    if m < -1e-8:
        n_fail += 1
        if m < worst_min_eig:
            worst_min_eig = m
            worst_case = (B, w.copy(), A.copy())

print(f"PSD 输入样本数: {n_psd_input}/{n_trials}")
print(f"M 非 PSD 的样本数: {n_fail} ({100*n_fail/max(n_psd_input,1):.2f}%)")
print(f"最劣 min-eig(M): {worst_min_eig:.6f}")
if worst_case is not None:
    B, w, A = worst_case
    print("最劣构型 B =", B)
    print("w =", w)
    print("Lambda_tilde^sf =\n", A)

print()
print("=" * 70)
print("负对照：Lambda_tilde^sf 逐元非负（而不仅仅 PSD）时，M 是否总 PSD？")
print("=" * 70)
n_trials2 = 5000
n_fail2 = 0
for trial in range(n_trials2):
    B = rng.integers(2, 7)
    w = rng.uniform(0.1, 5.0, size=B)
    A = rng.uniform(0, 3.0, size=(B, B))
    A = (A + A.T) / 2  # 对称化，仍逐元非负
    M = L_signless(A, w)
    eM = np.linalg.eigvalsh(M)
    if eM.min() < -1e-8:
        n_fail2 += 1
print(f"逐元非负输入下 M 非 PSD 样本数: {n_fail2}/{n_trials2}  (应为 0，验证 signless Laplacian 标准结论)")
