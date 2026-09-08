"""
P2 §5：静态极限的 DOS 缺损泛函
  P_s/λ = (ω₀/2)∫dε N̄(ε)/(|ε|+ω₀)²  −  1
由 GPS99 Eq.(8) 第一项（P_pol）在 T=0、ω_k=0、q=0、ω_q→0 取得，加上 P_ex 的静态极限 −λ（Eq. 20，与 DOS 形状无关）。
用 Mishra 数据集里的三种 DOS 文件算 H3S/D3S 在 160/200 GPa 的值，并与第一性原理顶角修正 δ=ln(Tc_na/Tc_adia) 对比。
运行：python theory/derivations/p2_dos_deficit.py
"""
import sys
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parents[1] / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"


def load_qe_dos(path, ef_ev):
    raw = np.loadtxt(path, comments="#")
    eps = (raw[:, 0] - ef_ev) * 1e3
    nbar = raw[:, 1] / np.interp(ef_ev, raw[:, 0], raw[:, 1])
    return eps, nbar


def load_epw_dos(path):
    ef = None
    rows = []
    for line in open(path):
        if line.startswith("#"):
            if "EFermi[eV]" in line and ef is None:
                ef = float(line.split()[2])
            continue
        p = line.split()
        if len(p) >= 2:
            try:
                rows.append((float(p[0]), float(p[1])))
            except ValueError:
                pass
    a = np.array(rows)
    eps = (a[:, 0] - ef) * 1e3
    nbar = a[:, 1] / np.interp(0.0, eps, a[:, 1])
    return eps, nbar


def deficit_integral(w0, eps, nbar):
    """(ω₀/2)∫dε N̄(ε)/(|ε|+ω₀)²；N̄≡1、无限带时 → 1；有限带 [−E_F, E_F] → 1 − ω₀/(E_F+ω₀)。"""
    return 0.5 * w0 * np.trapezoid(nbar / (np.abs(eps) + w0) ** 2, eps)


def tail_const(w0, e_lo, e_hi, nbar_tail=1.0):
    """把窗口外按常数 N̄=nbar_tail 延拓到 ±∞ 的补偿项。"""
    return 0.5 * w0 * nbar_tail * (1.0 / (abs(e_lo) + w0) + 1.0 / (abs(e_hi) + w0))


W0 = {"H3S": 110.0, "D3S": 77.8}
# 第一性原理（Mishra 数据集，FSR，μ*=0.16，非谐）：δ = ln(Tc_na/Tc_adia)
DATA = {("H3S", 200): -0.208, ("D3S", 200): -0.149, ("H3S", 160): -0.076, ("D3S", 160): -0.045}
# 每个 (体系, 压力) 自己的 λ 与 ω_log（本求解器对 Figshare 非谐 .a2f 的积分；reviewer M2：不能把 200 GPa 的值复用到 160 GPa）
LAM = {("H3S", 200): 1.721, ("D3S", 200): 1.784, ("H3S", 160): 1.876, ("D3S", 160): 1.941}
WLOG = {("H3S", 200): 110.0, ("D3S", 200): 85.5, ("H3S", 160): 106.9, ("D3S", 160): 82.0}
S_LAMBDA = 1.27  # 本求解器对 200 GPa 非谐谱的 ∂lnTc/∂lnλ（P2 §2）；160 GPa 未单独算，沿用

print("静态极限 DOS 缺损泛函  P_s/λ = I − 1,  I ≡ (ω₀/2)∫dε N̄(ε)/(|ε|+ω₀)²")
print("数据侧对照：P̄ ≡ δ/(s_λ·λ)，λ=1.72（200 GPa 非谐），s_λ=1.27\n")

# 1) QE 全带表（只有 200 GPa）
eps_q, nb_q = load_qe_dos(ROOT / "H3S-200GPa/band/H3S.dos", 17.619)
print("A. QE 全带 DOS，200 GPa（ε ∈ [−25.2, +41.6] eV）")
for mat in ["H3S", "D3S"]:
    I = deficit_integral(W0[mat], eps_q, nb_q)
    print(f"  {mat}: I = {I:.4f} → P_s/λ = {I-1:+.4f};  E_F^eff ≡ ω₀(1/|P_s/λ| − 1) = {W0[mat]*(1/abs(I-1)-1)/1e3:.2f} eV")
print(f"  同位素比 P_s(D)/P_s(H) = {(deficit_integral(77.8,eps_q,nb_q)-1)/(deficit_integral(110,eps_q,nb_q)-1):.3f}  (M^-1/2 = 0.707)")
print(f"  常数 DOS 对照：E_F=24.1 eV → {-(110/(24135+110)):+.4f}；E_F=1.0 eV → {-(110/(1000+110)):+.4f}")

# 2) EPW 窗口 DOS（fsthick=1 eV，就是 EPW fbw 计算实际用到的态）
print("\nB. EPW 窗口 DOS（±0.95 eV，这是 EPW fbw 实际包含的电子态；两种处理：截断 / 以 N̄(边界) 常数延拓）")
files = {("H3S", 200): "H3S-200GPa/EPW_harmonic/fsr_na/H3S.dos", ("H3S", 160): "H3S-160GPa/EPW_anharmonic/H3S.dos"}
res = {}
for (mat, P), f in files.items():
    eps, nb = load_epw_dos(ROOT / f)
    for iso in ["H3S", "D3S"]:
        w0 = W0[iso]
        I_win = deficit_integral(w0, eps, nb)
        I_ext = I_win + tail_const(w0, eps.min(), eps.max(), nbar_tail=0.5 * (nb[0] + nb[-1]))
        res[(iso, P)] = (I_win - 1, I_ext - 1)
        print(f"  {iso} DOS@{P} GPa (ω₀={w0:.0f} meV): 截断 P_s/λ = {I_win-1:+.4f} | 延拓 P_s/λ = {I_ext-1:+.4f}   N̄(−100/+100)={np.interp(-100,eps,nb):.2f}/{np.interp(100,eps,nb):.2f}")

print("\nC. 与第一性原理顶角修正对照（P̄_data = δ/(s_λ λ)，线性化；各 (体系,压力) 用各自 λ）")
for (mat, P), d in DATA.items():
    pbar = d / (S_LAMBDA * LAM[(mat, P)])
    w, e = res[(mat, P)]
    print(f"  {mat} {P} GPa (λ={LAM[(mat,P)]:.2f}, ω_log={WLOG[(mat,P)]:.0f}): P̄_data = {pbar:+.4f}   静态泛函（截断/延拓） = {w:+.4f} / {e:+.4f}")

print("\nD. 窗口敏感性（reviewer M1）：QE 全带 DOS 按不同窗口截断，200 GPa H3S")
for win in [500, 950, 2000, 5000, None]:
    if win is None:
        I = deficit_integral(110.0, eps_q, nb_q); label = "全带"
    else:
        m = np.abs(eps_q) < win; I = deficit_integral(110.0, eps_q[m], nb_q[m]); label = f"±{win/1000:.2f} eV"
    print(f"  {label:9s}: P_s/λ = {I-1:+.4f}")
print(f"\n  压力比 160/200（H3S）：数据 {DATA[('H3S',160)]/DATA[('H3S',200)]:.2f}；泛函截断 {res[('H3S',160)][0]/res[('H3S',200)][0]:.2f}；延拓 {res[('H3S',160)][1]/res[('H3S',200)][1]:.2f}")
print(f"  同位素比 D/H（200 GPa）：数据 {DATA[('D3S',200)]/DATA[('H3S',200)]:.2f}；泛函截断 {res[('D3S',200)][0]/res[('H3S',200)][0]:.2f}；延拓 {res[('D3S',200)][1]/res[('H3S',200)][1]:.2f}")
