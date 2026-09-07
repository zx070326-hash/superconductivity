"""
用 Mishra 数据集的真实 H3S/D3S α²F + 真实 QE DOS（band/H3S.dos，E_F=17.619 eV）跑本求解器的能量依赖 DOS 路径，
与数据集里 EPW 的 fbw_adia（全带宽、绝热）Tc 对照，标定 DOS 实现。
运行：python theory/derivations/p1_mishra_dos_calibration.py   （约 5–15 分钟）
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

# --- 谱 ---
sp = {k: load_epw_a2f(p) for k, p in FILES.items()}

# EPW 参考值（p1_mishra_tc_table.py，200 GPa，μ*=0.16，wscut=2.5 eV）
ref = {("har", "fsr"): (245, 180), ("har", "fbw"): (234, 173), ("anh", "fsr"): (192, 152), ("anh", "fbw"): (189, 149)}

print("\n=== 本求解器 vs EPW（μ*=0.16 固定；本求解器 ω_c = 2.5 eV 与 EPW 的 wscut 一致） ===")
print(f"{'kind':4s} {'DOS':8s} {'Tc_H':>7} {'Tc_D':>7} {'α_H':>6} | {'EPW Tc_H':>8} {'EPW Tc_D':>8} {'EPW α_H':>7}")
WCUT = 2500.0
for kind in ["har", "anh"]:
    for dos_name, dos in [("const", DOS.constant()), ("QE-table", dos_tab)]:
        tcs = []
        for mat in ["H3S", "D3S"]:
            s = Solver(sp[(mat, kind)], mu_star=0.16, w_cut=WCUT, dos=dos, dos_n_eps=400, dos_energy_cutoff_ratio=2.0)
            tcs.append(s.tc())
        a = -np.log(tcs[1] / tcs[0]) / np.log(2)
        key = (kind, "fsr" if dos_name == "const" else "fbw")
        rh, rd = ref[key]
        print(f"{kind:4s} {dos_name:8s} {tcs[0]:7.1f} {tcs[1]:7.1f} {a:6.3f} | {rh:8.0f} {rd:8.0f} {-np.log(rd/rh)/np.log(2):7.3f}")
