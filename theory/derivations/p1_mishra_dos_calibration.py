"""
用 Mishra 数据集的真实 H3S/D3S α²F + 真实 QE DOS（band/H3S.dos，E_F=17.619 eV）跑本求解器的能量依赖 DOS 路径，
与数据集里 EPW 的 fbw_adia（全带宽、绝热）Tc 对照，标定 DOS 实现。
运行：python theory/derivations/p1_mishra_dos_calibration.py   （约 10–20 分钟，含 fix_mu="conserve" 两列）

fix_mu="conserve" 两列是为了对照 EPW 的 muchem（电子数守恒的化学势修正，见
theory/solvers/README.md「化学势守恒」一节）：
  - "QE-table+cons"      用与 "QE-table" 完全相同的 DOS 表和 dos_energy_cutoff_ratio=2.0；
  - "QE-table+cons(±1eV)" 额外做的探索性对照：把 DOS 表先截到 EPW fsthick=1.0 eV 的窗口再求
    mu_shift，检验「mu_shift 对能量窗口很敏感」这个已知局限是否是 gap 没有如预期收窄的原因。

脚本末尾另有一段独立的探索：把 QE DOS 表整段（不截窗口，两端 N̄ 自然 -> 0）喂给
fix_mu="conserve"，检验「用真实能带的自然边界代替人为窗口，能否消除已知局限里
"贯穿整个积分域的不对称"这个 mu_shift 不收敛机制」。结论：这个机制确实被消除了，
但换来一个不同性质的问题——标准 dos_n_eps=200-800 在覆盖整段表格（~67eV 宽）时
求积混叠、mu_shift 混沌跳动，不满足「dos_n_eps=200/400/800 之间平滑收敛」这个
门槛，因此没有被采纳为推荐用法，详见该小节的打印与 theory/solvers/README.md
「化学势守恒」一节末尾的记录。
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "solvers"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from eliashberg_iso import Spectrum, Solver, DOS  # noqa: E402
from p1_mishra_a2f import load_epw_a2f, FILES  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

DATA = ROOT / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"
E_F_EV = 17.619
WINDOW_EV = 3.0  # 只取 E_F ± 3 eV 的 DOS 表（远处对 Tc 无影响，且避免 tan 代换网格浪费）

# --- 读 QE DOS，归一化到 N̄(E_F)=1，能量以 meV 相对 E_F ---
raw = np.loadtxt(DATA / "H3S-200GPa/band/H3S.dos", comments="#")
E = raw[:, 0]; N = raw[:, 1]
mask = np.abs(E - E_F_EV) < WINDOW_EV
eps_meV = (E[mask] - E_F_EV) * 1e3
nbar = N[mask] / np.interp(E_F_EV, E, N)
dos_tab = DOS.from_table(eps_meV, nbar)
print(f"QE DOS: N(E_F)={np.interp(E_F_EV, E, N):.4f} states/eV; window ±{WINDOW_EV} eV, {mask.sum()} points")
for e in [-500, -200, -100, -50, 0, 50, 100, 200, 500]:
    print(f"  N̄({e:+5d} meV) = {np.interp(e, eps_meV, nbar):.3f}")

# 探索性对照表：截到 EPW fsthick=1.0 eV 的窗口（只给 fix_mu="conserve" 用，
# 检验 mu_shift 对窗口大小的敏感性，见模块 docstring）
mask_1ev = np.abs(E - E_F_EV) < 1.0
eps_meV_1ev = (E[mask_1ev] - E_F_EV) * 1e3
nbar_1ev = N[mask_1ev] / np.interp(E_F_EV, E, N)
dos_tab_1ev = DOS.from_table(eps_meV_1ev, nbar_1ev)

# 全带 QE DOS 表（不截窗口）：E 覆盖整段计算区间，ε = E - E_F ∈ [-25239, +41641] meV，
# 两端 N̄ 自然 -> 0（真实能带边界，见下方打印），不像 ±3eV 窗口版本那样在窗口边界出现两侧
# 不对称的 N̄(edge) 突变（0.799 vs 0.455）。DOS.from_table 边界外按端点值延拓，窗口版本的
# 不对称延拓会在 mu_shift 求根中引入远场 ε→-ε 不对称的二阶残差（不收敛的根源之一，见
# README「已知局限」）；全带表两端对称 ->0，理论上能消除这个机制，见下方"探索"小节的检验。
eps_meV_full = (E - E_F_EV) * 1e3
nbar_full = N / np.interp(E_F_EV, E, N)
dos_tab_full = DOS.from_table(eps_meV_full, nbar_full)
print(f"\nQE DOS 全表（不截窗口）：ε ∈ [{eps_meV_full[0]:.0f}, {eps_meV_full[-1]:.0f}] meV，"
      f"{len(eps_meV_full)} 点，两端 N̄ = {nbar_full[0]:.3g} / {nbar_full[-1]:.3g}（应 ≈0）")

# --- 谱 ---
sp = {k: load_epw_a2f(p) for k, p in FILES.items()}

# EPW 参考值（p1_mishra_tc_table.py，200 GPa，μ*=0.16，wscut=2.5 eV）
ref = {("har", "fsr"): (245, 180), ("har", "fbw"): (234, 173), ("anh", "fsr"): (192, 152), ("anh", "fbw"): (189, 149)}

print("\n=== 本求解器 vs EPW（μ*=0.16 固定；本求解器 ω_c = 2.5 eV 与 EPW 的 wscut 一致） ===")
print(f"{'kind':4s} {'DOS':22s} {'Tc_H':>7} {'Tc_D':>7} {'α_H':>6} | {'EPW Tc_H':>8} {'EPW Tc_D':>8} {'EPW α_H':>7}")
WCUT = 2500.0
# fix_mu="conserve" 的数值求根成本远高于 "none"（每个 T 上多一层二分求 mu_shift，
# 见 README「化学势守恒」一节），这里放松 dos_tol/mu_shift_tol/tc(xtol=...) 换取可行
# 的运行时间（单次 tc() 从默认容差下的十分钟量级降到 ~25–35s）。只影响新增的两个
# "+cons" 列，不改变 "const"/"QE-table" 两列的 dos_n_eps=400、默认容差配置。
FAST_CONSERVE_KW = dict(dos_tol=1e-6, mu_shift_tol=1e-2, mu_shift_max_expand=20, mu_shift_max_iter=30)

mu_shift_log = []  # (kind, dos_name, mat, mu_shift_meV, converged)，跑完后打印一份小结
for kind in ["har", "anh"]:
    for dos_name, dos, fix_mu, dos_n_eps in [
        ("const", DOS.constant(), "none", 400),
        ("QE-table", dos_tab, "none", 400),
        ("QE-table+cons", dos_tab, "conserve", 200),
        ("QE-table+cons(±1eV)", dos_tab_1ev, "conserve", 200),
    ]:
        extra_kw = FAST_CONSERVE_KW if fix_mu == "conserve" else {}
        tc_kw = dict(xtol=1e-4) if fix_mu == "conserve" else {}
        tcs = []
        for mat in ["H3S", "D3S"]:
            s = Solver(sp[(mat, kind)], mu_star=0.16, w_cut=WCUT, dos=dos, dos_n_eps=dos_n_eps,
                       dos_energy_cutoff_ratio=2.0, fix_mu=fix_mu, **extra_kw)
            tc = s.tc(**tc_kw)
            tcs.append(tc)
            if fix_mu == "conserve":
                mu, conv = s.mu_shift_at_T(tc)
                mu_shift_log.append((kind, dos_name, mat, mu, conv))
        a = -np.log(tcs[1] / tcs[0]) / np.log(2)
        key = (kind, "fsr" if dos_name == "const" else "fbw")
        rh, rd = ref[key]
        print(f"{kind:4s} {dos_name:22s} {tcs[0]:7.1f} {tcs[1]:7.1f} {a:6.3f} | {rh:8.0f} {rd:8.0f} {-np.log(rd/rh)/np.log(2):7.3f}")

print("\n=== fix_mu=\"conserve\" 的 mu_shift（在各自 Tc 处求得，meV，负号=化学势下移） ===")
for kind, dos_name, mat, mu, conv in mu_shift_log:
    print(f"{kind:4s} {dos_name:22s} {mat:4s} mu_shift={mu:9.1f} meV  converged={conv}")

# --- 探索：QE-table(full) + conserve 是否可行？----------------------------------
# 团队负责人的假设：±3eV 窗口版本的 mu_shift 不收敛（随 dos_energy_cutoff_ratio 增大而单调
# 增长，见 README「已知局限」）是因为窗口边界两侧 N̄ 不对称（0.799 vs 0.455），DOS.from_table
# 在窗口外按端点值延拓，破坏了远场 ε→-ε 对称性；真实能带两端 N̄→0，天然对称，喂全表应该能收敛。
#
# 检验：dos_energy_cutoff_ratio=17.0 让 Solver._dos_e_max() 的
# max(ratio*w_cut, |eps_tab[0]|, |eps_tab[-1]|) 覆盖到整段表格边界（|eps|~25-42 eV），
# 在各 kind/mat 组合的 EPW fbw 参考 Tc 处，扫 dos_n_eps=200/400/800 下 mu_shift 是否收敛——
# 这是本次任务要求先过的门槛，通过了才值得再跑完整 tc() 标定。
print("\n=== 探索：QE-table(full) + conserve，mu_shift 对 dos_n_eps 的收敛性 ===")
FULL_RATIO = 17.0  # max(|eps_full[0]|,|eps_full[-1]|)/WCUT ≈ 42500/2500 ≈ 17.0，覆盖整段表格边界
full_conv_log = []
for kind in ["har", "anh"]:
    for mat in ["H3S", "D3S"]:
        T = float(ref[(kind, "fbw")][0 if mat == "H3S" else 1])  # EPW fbw 参考 Tc，仅作探测点
        row = []
        for n_eps in (200, 400, 800):
            s = Solver(sp[(mat, kind)], mu_star=0.16, w_cut=WCUT, dos=dos_tab_full, dos_n_eps=n_eps,
                       dos_energy_cutoff_ratio=FULL_RATIO, fix_mu="conserve", **FAST_CONSERVE_KW)
            mu, conv = s.mu_shift_at_T(T)
            row.append(mu)
        full_conv_log.append((kind, mat, T, row))
        print(f"{kind:4s} {mat:4s} T={T:6.1f}K  mu_shift(n_eps=200/400/800) = "
              f"{row[0]:9.1f} / {row[1]:9.1f} / {row[2]:9.1f} meV")

# 判据：200/400/800 三点之间相邻差应 < 100 meV（平滑收敛到 O(10-100meV) 级别）；
# 观察到的实际情况是几十到近五千 meV 之间无单调趋势的跳动——不是"慢收敛"，是求积混叠。
spread_ok = all(max(abs(r[1] - r[0]), abs(r[2] - r[1])) < 100.0 for _k, _m, _T, r in full_conv_log)
print(f"\n结论：dos_n_eps=200/400/800 之间 mu_shift {'平滑收敛' if spread_ok else '不收敛（混沌跳动，非单调）'}。")
if not spread_ok:
    print("  诊断：tan 代换 + 固定节点 Gauss-Legendre 求积在覆盖整段表格的 ~67 eV 宽域上（相对"
          "\n  ±3eV 窗口宽 ~20 倍），用本代码库其它地方常用的 dos_n_eps=200-800 会严重欠采样"
          "\n  H3S.dos 的精细结构，产生求积混叠，不是物理不收敛——这与「边界不对称延拓」是两个"
          "\n  不同的机制，全带表确实消除了后者，但暴露出了前者。"
          "\n  人工扫描 dos_n_eps 到 1600 以上才重新平滑收敛（H3S/har 单点核实：n_eps=1600，"
          "\n  Tc=229.05K、mu_shift=-250meV，vs EPW fbw 234K 差 4.95K；单次 tc() 耗时 ~407s，"
          "\n  比本脚本其它列贵一个数量级，且尚未确认对其余三个组合同样成立）。"
          "\n  按任务约定的停止条件：dos_n_eps=200/400/800 下不收敛，此处如实记录并停止，不再"
          "\n  往下算全带表的完整 tc()/α_H（那些数字建立在不可靠的 mu_shift 上会误导人）。"
          "\n  结论见 README「化学势守恒」节：DOS 路径保留用于 α 差分，不推荐用于绝对 Tc 标定。")

