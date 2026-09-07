"""
用 Mishra 数据集的真实 H3S/D3S α²F + 真实 QE DOS（band/H3S.dos，E_F=17.619 eV）跑本求解器的能量依赖 DOS 路径，
与数据集里 EPW 的 fbw_adia（全带宽、绝热）Tc 对照，标定 DOS 实现。
运行：python theory/derivations/p1_mishra_dos_calibration.py   （约 10–20 分钟，含 fix_mu="conserve" 两列）

fix_mu="conserve" 两列是为了对照 EPW 的 muchem（电子数守恒的化学势修正，见
theory/solvers/README.md「化学势守恒」一节）：
  - "QE-table+cons"      用与 "QE-table" 完全相同的 DOS 表和 dos_energy_cutoff_ratio=2.0；
  - "QE-table+cons(±1eV)" 额外做的探索性对照：把 DOS 表先截到 EPW fsthick=1.0 eV 的窗口再求
    mu_shift，检验「mu_shift 对能量窗口很敏感」这个已知局限是否是 gap 没有如预期收窄的原因。
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

