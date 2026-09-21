# -*- coding: utf-8 -*-
"""
主会话从零复核：gapratio-budget agent 的「预算 1.0× => 实轴转换 ESSENTIAL」判决。

被复核的对象：theory/notes/_scratch_ir/_chk_realaxis_budget_2026-09-20.py
它的链条是：
    O_real = c * O_mat,   c := (Carbotte 实轴公式) / (松原轴 O)
    x := (dc/dlam)*Dlam/c,  y := (dc/dmu*)*Dmu*/c
    R_real = |dO_lam + O*x| / |dO_mu + O*y|
    「预算」= 使分子归零所需的 x = -dO_lam/O = -0.04178
    「经验估计」= Carbotte 给的 x = +0.04349      => |比| = 1.04 => ESSENTIAL

我要查三件事（每件都可能推翻那个判决）：
  D1 符号：预算是 -0.0418，估计是 +0.0435。正的 x **放大**分子、把 R 推高，
     它永远不可能「吃掉 lam 轴信号」。取绝对值比较 = 跨符号比较。
  D2 分母冻结：agent 算 c 时把 O_mat 固定成 O_base（脚本 207-208 行
     `O_mat_plus_approx = O_base`），于是 Carbotte 的**全部**强耦合增量都被记到 c 头上。
     但 O_mat 自己也随 lam 涨 —— dO_lam/O 就是那个涨幅，agent 手上已经有这个数。
     c 是**比值**，两边同涨要相消。
  D3 O_base 取值：agent 用 O_base=4.97。那是我 09-20 已**订正掉**的 AD-闭式支的量级
     （自洽支 O ∈ [4.41,4.82]）。用被撤回的支给分母 = 预算本身错。

  D4 一致性：若真的把 c := O_carb/O_mat 代回 O_real = c*O_mat，则 O_real ≡ O_carb。
     => 「用 Carbotte 当实轴估计」时 R_real **就是** R_Carbotte，不需要任何 x/y 中间量。
     我 09-20 的表里已经有 R_Carb = 1.3777。本脚本从零重算它。

严格性：D1/D2/D4 是 analytic（纯代数）；下面的数值是 numerical-evidence。
"""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, r"D:\physics\superconductivity")
import numpy as np
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

LAM0, W0, WCUT, MU0, N = 1.876, 106.9, 500.0, 0.13, 1600
DMU, DLAM = 0.05, 0.15
KB_EV = 8.617333262e-5


def carbotte(tc_K, wlog_meV):
    """Carbotte RMP 1990: 2D/kTc = 3.53[1 + 12.5 t^2 ln(1/2t)], t = kB Tc / wlog。"""
    t = KB_EV * tc_K / (wlog_meV * 1e-3)
    return 3.53 * (1.0 + 12.5 * t * t * np.log(1.0 / (2.0 * t))), t


def obs(lam, mu):
    spec = Spectrum.einstein(lam=lam, w0=W0)
    s = Solver(spectrum=spec, mu_star=mu, w_cut=WCUT, n_matsubara=N)
    tc = s.tc()
    wn, d, z, c = s.solve_gap(0.03 * tc, mixing=0.25, max_iter=3000)
    return tc, 2.0 * d[0] / (KB_MEV_PER_K * tc)


print("=" * 74, flush=True)
print("主会话独立复核：实轴修正预算（自洽 solver.tc() 支，N=%d）" % N, flush=True)
print("=" * 74, flush=True)

t0 = time.time()
pts = {}
for tag, lam, mu in (("0  ", LAM0, MU0),
                     ("mu+", LAM0, MU0 + DMU), ("mu-", LAM0, MU0 - DMU),
                     ("lam+", LAM0 * (1 + DLAM), MU0), ("lam-", LAM0 * (1 - DLAM), MU0)):
    tc, o = obs(lam, mu)
    oc, tpar = carbotte(tc, W0)
    pts[tag.strip()] = (tc, o, oc, tpar)
    print(f"  {tag:>4}  lam={lam:6.4f} mu*={mu:5.3f}  Tc={tc:8.3f}K  "
          f"O_mat={o:.4f}  t={tpar:.4f}  O_carb={oc:.4f}   [{time.time()-t0:.0f}s]", flush=True)

tc0, O0, oc0, _ = pts["0"]
dO_mu = 0.5 * (pts["mu+"][1] - pts["mu-"][1])
dO_lam = 0.5 * (pts["lam+"][1] - pts["lam-"][1])
dC_mu = 0.5 * (pts["mu+"][2] - pts["mu-"][2])
dC_lam = 0.5 * (pts["lam+"][2] - pts["lam-"][2])

print()
print("-" * 74)
print("A. 两条轴上的响应（半中心差；Dmu*=0.05, Dlam/lam=0.15）")
print(f"   松原轴   dO_mu  = {dO_mu:+.5f}   dO_lam  = {dO_lam:+.5f}   R_mat  = {abs(dO_lam/dO_mu):.4f}")
print(f"   Carbotte dC_mu  = {dC_mu:+.5f}   dC_lam  = {dC_lam:+.5f}   R_carb = {abs(dC_lam/dC_mu):.4f}")
print(f"   O_mat(基) = {O0:.4f}   <-- agent 用的是 4.97（已撤回的 AD 支）  [D3]")

print()
print("-" * 74)
print("B. [D4, analytic] c := O_carb/O_mat  =>  O_real = c*O_mat === O_carb")
print("   => 「以 Carbotte 为实轴估计」时，R_real 就是 R_carb，不经过 x/y。")
print(f"   R_real(Carbotte 估计) = {abs(dC_lam/dC_mu):.4f}")
print(f"   相对翻转阈 R=1 的余量 = {abs(dC_lam/dC_mu)-1:+.4f}")
print(f"   已消耗的余量比例 = (R_mat - R_carb)/(R_mat - 1) = "
      f"{(abs(dO_lam/dO_mu)-abs(dC_lam/dC_mu))/(abs(dO_lam/dO_mu)-1)*100:.1f}%")

# x / y 的正确算法：c 是比值，两边都动
x_correct = dC_lam / oc0 - dO_lam / O0          # = Dln c 沿 lam
y_correct = dC_mu / oc0 - dO_mu / O0            # = Dln c 沿 mu*
budget = -dO_lam / O0                            # 使分子归零所需的 x

print()
print("-" * 74)
print("C. [D1/D2] agent 的 x 与正确的 x")
print(f"   预算（使 lam 信号归零所需）  x_budget = -dO_lam/O = {budget:+.6f}")
print(f"   agent 的估计（冻结分母）     x_agent  = Dln O_carb   = {dC_lam/oc0:+.6f}")
print(f"   正确的估计（c 是比值）       x_correct= Dln O_carb - Dln O_mat = {x_correct:+.6f}")
print(f"       其中 Dln O_carb = {dC_lam/oc0:+.6f},  Dln O_mat = {dO_lam/O0:+.6f}")
print()
print(f"   [D1 符号] x_budget<0 而 x_agent>0：agent 把一个**同向放大** lam 信号的量")
print(f"             拿去和一个**抵消** lam 信号的预算比绝对值 —— 跨符号比较，")
print(f"             它那个 1.04x 在方向上不成立。")
print(f"             (自我订正：我落笔前以为 x_correct 也 >0；实测 x_correct={x_correct:+.6f}")
print(f"              **是负的** —— 方向确实在抵消侧，只是量级不够。D1 只管 agent 的估计。)")
print(f"   [D2 冻结] agent 把 Dln O_mat 记成 0；它手上已有该数 = {dO_lam/O0:+.6f}，")
print(f"             与 Dln O_carb={dC_lam/oc0:+.6f} 相消 "
      f"{(1-abs(x_correct/(dC_lam/oc0)))*100:.1f}%。这是主缺陷。")
print(f"   => 正确比值 |x_budget/x_correct| = {abs(budget/x_correct):.1f}x  "
      f"（agent 报 1.0x）=> 按它自己的分档规则属于「预算 >5x，不承重」")

print()
print("-" * 74)
print("D. 用 agent 自己的 R_real 公式代入正确的 (x,y)，检查是否回到 R_carb")
num = abs(dO_lam + O0 * x_correct)
den = abs(dO_mu + O0 * y_correct)
print(f"   x_correct={x_correct:+.6f}  y_correct={y_correct:+.6f}")
print(f"   R_real(公式) = {num/den:.4f}   vs   R_carb(直接) = {abs(dC_lam/dC_mu):.4f}"
      f"   差 {abs(num/den - abs(dC_lam/dC_mu)):.1e}")
print("   （两者应一致到 O_mat/O_carb 的一阶；不一致说明公式或 x/y 定义有错。）")

print()
print("-" * 74)
print("E. 余量到底被哪条轴吃掉？（agent 的预算整篇只写在 lam 轴上）")
print(f"   lam 轴：分子 {dO_lam:+.5f} -> {dO_lam + O0*x_correct:+.5f}   "
      f"变化 {(abs(dO_lam + O0*x_correct)/abs(dO_lam)-1)*100:+.1f}%")
print(f"   mu* 轴：分母 {dO_mu:+.5f} -> {dO_mu + O0*y_correct:+.5f}   "
      f"变化 {(abs(dO_mu + O0*y_correct)/abs(dO_mu)-1)*100:+.1f}%")
print("   => R 从 2.008 掉到 1.378 主要来自**分母被放大**（mu* 轴），不是分子被抵消。")
print("      agent 的预算写在 lam 轴上 => 它连自己要算的那条轴都选错了。")

print()
print("=" * 74)
Rc = abs(dC_lam / dC_mu)
print(f"判决：R_real（唯一可用的实轴估计）= {Rc:.4f} > 1 => 判决**不翻转**。")
print(f"      实轴修正确实承重（把 R 从 {abs(dO_lam/dO_mu):.3f} 压到 {Rc:.3f}，吃掉 "
      f"{(abs(dO_lam/dO_mu)-Rc)/(abs(dO_lam/dO_mu)-1)*100:.0f}% 的余量），")
print(f"      但 agent 的『1.0× 刀口、必须先造实轴求解器』不成立。")
print("=" * 74)
