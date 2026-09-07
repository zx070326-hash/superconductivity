"""
用 Mishra et al. Figshare 数据集（CC BY 4.0, doi:10.6084/m9.figshare.31136635.v1）里的真实 α²F(ω)
跑各向同性 Eliashberg：H3S / D3S，谐波 / 非谐，200 GPa。直接得到部分同位素系数 α_H。
运行：python theory/derivations/p1_mishra_a2f.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "solvers"))
from eliashberg_iso import Spectrum, Solver, mu_star_at_cutoff  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

DATA = ROOT / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"
FILES = {
    ("H3S", "har"): DATA / "H3S-200GPa/EPW_harmonic/H3S.a2f",
    ("D3S", "har"): DATA / "D3S-200GPa/EPW_harmonic/D3S.a2f",
    ("H3S", "anh"): DATA / "H3S-200GPa/EPW_anharmonic/H3S.a2f",
    ("D3S", "anh"): DATA / "D3S-200GPa/EPW_anharmonic/D3S.a2f",
}
WCUT_RATIO = 5.0


def load_epw_a2f(path: Path) -> Spectrum:
    """EPW .a2f：第一行表头，之后 'w[meV] a2f int' 直到 'Integrated el-ph coupling' 行。"""
    w, a = [], []
    for line in path.read_text().splitlines()[1:]:
        parts = line.split()
        if len(parts) < 3:
            break
        try:
            x, y = float(parts[0]), float(parts[1])
        except ValueError:
            break
        w.append(x)
        a.append(max(y, 0.0))
    w = np.array(w)
    a = np.array(a)
    return Spectrum.from_table(w, a, name=path.stem)


spectra = {k: load_epw_a2f(p) for k, p in FILES.items()}

print("=== 谱的矩（本求解器积分） ===")
print(f"{'case':10s} {'λ':>6} {'ω_log(meV)':>11} {'ω2/ω_log':>9} {'ω_max':>7}")
for (mat, kind), sp in spectra.items():
    print(f"{mat}-{kind:4s} {sp.lambda_total():6.3f} {sp.omega_log():11.2f} {sp.omega_2()/sp.omega_log():9.3f} {sp.w_max:7.1f}")

print("\n=== Tc 与部分同位素系数 α_H = -ln(Tc_D/Tc_H)/ln2，ω_c = 5·ω_log(H) ===")
print("约定 A：μ* 固定（Errea 2015 / Mishra 的做法）；约定 B：固定裸 μ，μ*(ω_c) 随 D 谱的 ω_c 一起重整化")
for kind in ["har", "anh"]:
    spH, spD = spectra[("H3S", kind)], spectra[("D3S", kind)]
    wcH = WCUT_RATIO * spH.omega_log()
    wcD = WCUT_RATIO * spD.omega_log()
    print(f"\n--- {kind} ---")
    print(f"{'μ*':>5} | {'Tc_H':>7} {'Tc_D(A)':>8} {'α_H(A: μ*固定)':>15} | {'μ*_D(B)':>8} {'Tc_D(B)':>8} {'α_H(B: μ*随ω_c)':>16} | {'AD Tc_H':>8} {'AD Tc_D':>8} {'α_AD':>6}")
    for mu in [0.10, 0.13, 0.16, 0.20]:
        tcH = Solver(spH, mu_star=mu, w_cut=wcH).tc()
        tcD_A = Solver(spD, mu_star=mu, w_cut=wcD).tc()
        muD = mu_star_at_cutoff(mu, wcH, wcD)
        tcD_B = Solver(spD, mu_star=muD, w_cut=wcD).tc()
        aA = -np.log(tcD_A / tcH) / np.log(2)
        aB = -np.log(tcD_B / tcH) / np.log(2)
        adH, adD = spH.allen_dynes(mu), spD.allen_dynes(mu)
        print(f"{mu:5.2f} | {tcH:7.1f} {tcD_A:8.1f} {aA:15.3f} | {muD:8.4f} {tcD_B:8.1f} {aB:16.3f} | {adH:8.1f} {adD:8.1f} {-np.log(adD/adH)/np.log(2):6.3f}")

print("\n=== 两通道分解（μ*=0.16 固定，本求解器 Tc） ===")
for kind in ["har", "anh"]:
    spH, spD = spectra[("H3S", kind)], spectra[("D3S", kind)]
    tcH = Solver(spH, 0.16, WCUT_RATIO * spH.omega_log()).tc()
    tcD = Solver(spD, 0.16, WCUT_RATIO * spD.omega_log()).tc()
    a = -np.log(tcD / tcH) / np.log(2)
    ch_w = -np.log(spD.omega_log() / spH.omega_log()) / np.log(2)
    dlnlam = np.log(spD.lambda_total() / spH.lambda_total()) / np.log(2)
    print(f"{kind}: α_H={a:.3f}  ω_log 通道={ch_w:.3f}  dlnλ/dlnM={dlnlam:+.3f}  λ 通道(余项)={a-ch_w:+.3f}")
