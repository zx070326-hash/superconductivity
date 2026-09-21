#!/usr/bin/env python
"""t5 的 `solve_einstein_tc` 审计: 它解的根本不是 AD75 方程, 且 wbar2 全程不进来。

t5 笔记 §4.1 把 13 个 f/f_E>1 判为「N=400 截断赝物」(断言, 从未实测)。
在测 N 收敛之前, 先审被比较的两个对象本身 —— CLAUDE.md「判词对象错」。

发现 (待本脚本证实):
  (A) 实谱走 AD75 完整对称矩阵 (top eigenvalue = 0);
      Einstein 走 max_n[(pi/T) lam_n/(2 w_n)] = 1 的**对角单模**判据。两个不同的方程。
  (B) 该对角判据把 T 放在**分母** (正确的松原核是 pi*T*lam_n/(2 w_n), T 在分子)
      ==> 判据退化为 lam_0(T)/(2 T^2) = 1, 低温极限 T = sqrt(lam/2) [meV]
      ==> **Einstein 频率 wbar2 完全不进来**。
  (C) ==> f/f_E = Tc_real/(8.2057 sqrt(lam)) [K], 是 (Tc_real, lam) 的函数, 零谱信息。
      「Einstein 极值性检验」没有测到 Einstein 谱。第五起「先数自由度」同形事故。

验收门 (t5 自己的附录, 不是我编的): Al2CoNi lam=0.2198 -> t5 印 Tc(Einstein)=3.84 K。
等级 numerical-evidence (纯属复现对方代码的行为)。
"""
import os, sys, math, importlib.util
import numpy as np

# 被审对象在 t5 的分支上 (UNMERGED)。合并后走仓库内相对路径; 未合并时回落到 worktree 绝对路径。
REL = "theory/notes/_scratch_ir/_cerqueira_eliashberg_solver_2026-09-21.py"
WT = ("D:/physics/superconductivity/.claude/worktrees/agent-a251c94abbc4f0ae2/" + REL)
SOLVER = REL if os.path.isfile(REL) else WT
if not os.path.isfile(SOLVER):
    sys.exit("找不到被审求解器: t5 分支未合并且 worktree 路径不存在 -> " + WT)
spec = importlib.util.spec_from_file_location("t5solver", SOLVER)
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)

MEV_TO_K = 11.604518121550082
RY_TO_MEV = 13605.693122994
RY_TO_K = RY_TO_MEV * MEV_TO_K

print("=" * 78)
print("(B)+(C) Einstein 解的 wbar2 不变性:  T = sqrt(lam/2) meV = 8.2057 sqrt(lam) K ?")
print("=" * 78)
print(f"{'lam':>8} " + "".join(f"{'w2=%g' % w:>12}" for w in (5, 25, 100, 400))
      + f"{'闭式 8.2057Vlam':>18}")
for lam in (0.029261, 0.084, 0.2198, 0.486):
    row = f"{lam:>8.4f} "
    for w2 in (5.0, 25.0, 100.0, 400.0):
        row += f"{S.solve_einstein_tc(lam, w2, N=400):>12.6f}"
    row += f"{MEV_TO_K * math.sqrt(lam / 2.0):>18.6f}"
    print(row)

print()
print("=" * 78)
print("验收门: t5 附录印 Al2CoNi lam=0.2198 -> Tc(Einstein)=3.84 K")
print("=" * 78)
print(f"  闭式 8.2057*sqrt(0.2198) = {MEV_TO_K*math.sqrt(0.2198/2):.4f} K"
      f"   t5 印 3.84 K")

# ---------------------------------------------------------------- 实谱侧
ROOT = "D:/physics/superconductivity/theory/data/external/cerqueira_raw/extracted/batch-q"
MATS = ["NaS2Rh_agm005113872", "Ca2AuBi_agm002259704", "Li2CaPb_agm002693901",
        "AlTaOs2_agm002632444", "AlRu2Ta_agm002635418", "Al2CoNi_agm003191160"]
NOTE = {"NaS2Rh_agm005113872": 1.001246, "Ca2AuBi_agm002259704": 1.000722,
        "Li2CaPb_agm002693901": 1.000188, "AlTaOs2_agm002632444": 1.000190,
        "AlRu2Ta_agm002635418": 1.000169, "Al2CoNi_agm003191160": 0.9868}

print()
print("=" * 78)
print("(A) 直接跑 t5 自己的两个求解器; 并查实谱侧二分是否**根本没有括到根**")
print("    (t5 的 find_tc_bisection 只守 lam_high>0, 不守 lam_low<=0)")
print("=" * 78)
hdr = (f"{'材料':<22}{'lam':>8}{'w2(meV)':>9}{'eig(T_lo)':>11}{'eig(T_hi)':>11}"
       f"{'Tc_real':>11}{'Tc_E':>9}{'ratio':>10}{'笔记':>10}")
print(hdr); print("-" * len(hdr))
for name in MATS:
    for sm in (5, 6):
        p = os.path.join(ROOT, name, f"a2F.dos{sm}")
        if os.path.isfile(p):
            break
    else:
        print(f"{name:<22} [缺 a2F.dos5/6]"); continue
    w_ry, a2f, _, _ = S.load_a2f_dos(p)
    lam, w2_mev, area_ry, M1_ry = S.compute_moments(w_ry, a2f)
    w_mev = w_ry * RY_TO_MEV

    T_lo = 1e-12
    T_hi = 0.3 * np.max(w_ry) * RY_TO_K
    e_lo = S.solve_eliashberg_isotropic(w_mev, a2f, T_lo, 400)
    e_hi = S.solve_eliashberg_isotropic(w_mev, a2f, T_hi, 400)

    tc = S.find_tc_bisection(w_mev, a2f, N=400)
    te = S.solve_einstein_tc(lam, w2_mev, N=400)
    r = tc / te if (tc and te) else float('nan')
    flag = "  <== 两端同号, 二分无根!" if (e_lo <= 0 and e_hi <= 0) else ""
    print(f"{name:<22}{lam:>8.4f}{w2_mev:>9.3f}{e_lo:>11.4f}{e_hi:>11.4f}"
          f"{tc:>11.3e}{te:>9.3f}{r:>10.6f}{NOTE.get(name, float('nan')):>10.6f}{flag}")

print()
print("dos5 存在? ", os.path.isfile(os.path.join(ROOT, MATS[0], "a2F.dos5")),
      " dos6 存在? ", os.path.isfile(os.path.join(ROOT, MATS[0], "a2F.dos6")))
