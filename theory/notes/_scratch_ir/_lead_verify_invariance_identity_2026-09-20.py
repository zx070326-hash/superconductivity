"""[analytic 验证 + numerical-evidence 演示] 不变性恒等式及其边界。

主会话命题：把观测量从「O 的值」改报成「O 自己零点的位置 M*」，判别性比值逐字不变。
  dM*/dp = -(dO/dp)/(dO/dM)  =>  |dM*/dp Dp| / |dM*/dq Dq| = |dO/dp Dp| / |dO/dq Dq|

本脚本：
 (1) 用一个与 alpha(M) 无关的随机非线性 O(M;p,q) 数值确认恒等式（不依赖仓库任何模型）；
 (2) 对主会话报的 PdH 数字做算术自洽核验；
 (3) ★ 检验我怀疑的边界：判别性分母里若有一项**不是** |dO/dnu|*Dnu 形式
     （典型：直接attach在观测量上的实验误差 sigma），该项不随同一 Jacobian 变换 =>
     不变性对它失效。这决定了「重参数化候选可否合并成一个等价类」。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import math

# ---------- (1) 与仓库模型无关的数值确认 ---------------------------------
# 随便一个有单根的非线性 O
def O(M, p, q):
    return math.sin(0.7 * M) - 0.3 * p * M**1.5 + 0.11 * q * math.exp(-0.4 * M) - 0.25


def bisect(f, lo, hi, tol=1e-14):
    flo = f(lo)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (f(mid) > 0) == (flo > 0):
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


p0, q0 = 0.20, 0.50
h = 1e-6
Mstar = bisect(lambda M: O(M, p0, q0), 0.3, 3.0)
dO_dM = (O(Mstar + h, p0, q0) - O(Mstar - h, p0, q0)) / (2 * h)
dO_dp = (O(Mstar, p0 + h, q0) - O(Mstar, p0 - h, q0)) / (2 * h)
dO_dq = (O(Mstar, p0, q0 + h) - O(Mstar, p0, q0 - h)) / (2 * h)

# M* 的导数：两条独立路径（隐函数式 vs 直接重解根）
dMs_dp_impl = -dO_dp / dO_dM
dMs_dq_impl = -dO_dq / dO_dM
dMs_dp_num = (bisect(lambda M: O(M, p0 + h, q0), 0.3, 3.0)
              - bisect(lambda M: O(M, p0 - h, q0), 0.3, 3.0)) / (2 * h)
dMs_dq_num = (bisect(lambda M: O(M, p0, q0 + h), 0.3, 3.0)
              - bisect(lambda M: O(M, p0, q0 - h), 0.3, 3.0)) / (2 * h)

Dp, Dq = 0.09, 0.0638
print("(1) 与仓库无关的独立确认")
print(f"  M*={Mstar:.12f}  dO/dM={dO_dM:.10f}")
print(f"  dM*/dp: implicit={dMs_dp_impl:.10f}  re-solved={dMs_dp_num:.10f}"
      f"  reldiff={abs(dMs_dp_impl-dMs_dp_num)/abs(dMs_dp_impl):.2e}")
print(f"  dM*/dq: implicit={dMs_dq_impl:.10f}  re-solved={dMs_dq_num:.10f}"
      f"  reldiff={abs(dMs_dq_impl-dMs_dq_num)/abs(dMs_dq_impl):.2e}")
r_O = abs(dO_dq * Dq) / abs(dO_dp * Dp)
r_M = abs(dMs_dq_impl * Dq) / abs(dMs_dp_impl * Dp)
print(f"  判别性比值  报 O 的值 = {r_O:.12f}")
print(f"              报 O 的零点 = {r_M:.12f}   商 = {r_M/r_O:.12f}")

# ---------- (2) 主会话 PdH 数字的算术自洽 --------------------------------
print("\n(2) 主会话 PdH 数字的算术自洽核验（我未重建 alpha 模型，只核算术）")
da_dmu, da_dd, da_dM = -3.1286, -1.0000, +0.3079
print(f"  dM*/dmu* = {-da_dmu/da_dM:+.4f} u   (主会话报 +10.1605)")
print(f"  dM*/ddelta = {-da_dd/da_dM:+.4f} u   (主会话报 +3.2477)")
rr_a = abs(da_dd * 0.0638) / abs(da_dmu * 0.09)
rr_M = abs((-da_dd/da_dM) * 0.0638) / abs((-da_dmu/da_dM) * 0.09)
print(f"  比值 报 alpha = {rr_a:.6f} ; 报 M* = {rr_M:.6f} ; 商 = {rr_M/rr_a:.8f}")

# ---------- (3) ★ 边界：分母里非 |dO/dnu|*Dnu 形式的项 -------------------
# 实验误差 sigma 直接 attach 在观测量上，不经 dO/dnu 进入。
# 情形 A：M* 只能由「测 O 再解根」得到 => sigma_M* = sigma_O/|dO/dM| => 仍然相消。
# 情形 B：M* 有独立测量协议（直接扫组分找零点），sigma_M* 是独立输入 => 不相消。
sigma_O = 0.02
print("\n(3) 边界检验：分母含实验误差项时")
denom_O = abs(dO_dp * Dp) + sigma_O
num_O = abs(dO_dq * Dq)
print(f"  报 O:  判别性 = {num_O/denom_O:.6f}   (sigma_O={sigma_O})")

sigma_M_propagated = sigma_O / abs(dO_dM)          # 情形 A
denom_MA = abs(dMs_dp_impl * Dp) + sigma_M_propagated
num_M = abs(dMs_dq_impl * Dq)
print(f"  报 M* 情形A(sigma 由 O 传播) = {num_M/denom_MA:.6f}"
      f"   商 = {(num_M/denom_MA)/(num_O/denom_O):.8f}  <- 仍然不变")

for sigma_M_indep in (0.002, 0.05):                 # 情形 B
    denom_MB = abs(dMs_dp_impl * Dp) + sigma_M_indep
    print(f"  报 M* 情形B(独立协议 sigma_M*={sigma_M_indep}) = {num_M/denom_MB:.6f}"
          f"   商 = {(num_M/denom_MB)/(num_O/denom_O):.4f}  <- 不变性失效")
