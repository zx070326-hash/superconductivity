# -*- coding: utf-8 -*-
"""承重值 q = −0.5565 会不会被库仑截断的**阶梯赝物**污染？

背景：`_chk_main_alpha_smoothness_dense_scan` 查实 Solver 第 565–568 行的库仑斜坡
`clip((w_cut−ω_n)/(2πT)+0.5,0,1)` 使 Tc 在 (λ,μ*) 上 C⁰ 而**非** C¹：
每当 `w_cut/(2πTc)` 穿过整数，α（本身是导数）出现扭结，|Δ²α| 最大/中位达 108×–541×。

q 也是一个「曲率型」量（logλ 的二次系数）⟹ 必须问同一个问题。
但两者的**结构不同**，先写下我的先验，再去测（免得测完才编理由）：
  · α 是**单锚点的逐点导数**，锚点距整数边界仅 0.10 格 ⟹ 扭结直接进入判决；
  · q 是 λ∈[0.69,2.59] 上 **9 点拟合**的系数，比值扫过很宽一段 ⟹ 多次穿越，
    阶梯应当落进**拟合残差**（rms）而不是系数——除非它恰好与 L² 对齐。
先验≠核实。本脚本测三件事：

§1  36 个工作点的 `w_cut/(2πTc)` 与它们到整数边界的距离（看有无贴边的点）
§2  ★ 判别实验：把 w_cut/ω̄₂ 在 [9.6, 10.4] 上**密扫 17 档**重拟合 q。
    阶梯赝物是**非光滑**的（扭结），ω_c 的物理依赖是**光滑**的
    （仓库已知 Tc(μ*>0) 对 ω_c 非单调，故 q 随 ω_c 动是预期内的**物理**）
    ⟹ 判据只能是 q(ratio) 的**光滑性**，不能是「q 动没动」。
    用与 α 同一把尺子：|Δ²q| 最大/中位 < 5 判光滑。
§3  阶梯若真进了链条，它进的是残差还是系数——把 rms 一并印出来对比。

方向注记：全为双侧点算与最小二乘，无单边界放缩。
严格性：numerical-evidence。
主会话 2026-09-20，从零。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, time
from pathlib import Path
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "solvers"))
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K   # noqa: E402

C_INF = 0.1827262477
MU = 0.10
NM = 400
LAM_GRID = np.exp(np.linspace(np.log(0.69), np.log(2.59), 9))
DATA = ROOT / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"
FILES = {"H3S-anh": DATA / "H3S-200GPa/EPW_anharmonic/H3S.a2f",
         "D3S-anh": DATA / "D3S-200GPa/EPW_anharmonic/D3S.a2f",
         "H3S-har": DATA / "H3S-200GPa/EPW_harmonic/H3S.a2f",
         "D3S-har": DATA / "D3S-200GPa/EPW_harmonic/D3S.a2f"}


def load(p):
    w, a = [], []
    for line in p.read_text().splitlines()[1:]:
        q = line.split()
        if len(q) < 3:
            break
        try:
            x, yv = float(q[0]), float(q[1])
        except ValueError:
            break
        w.append(x); a.append(max(yv, 0.0))
    return np.array(w), np.array(a)


RAW = {k: load(p) for k, p in FILES.items()}


def sweep(ratio):
    """四谱 × 9 个 λ，返回 (L, Σ, lnf, 比值 w_cut/(2πTc))。"""
    L, S, Y, NC = [], [], [], []
    for k, (w, a) in RAW.items():
        sp0 = Spectrum.from_table(w, a)
        l0, w2, wl = sp0.lambda_total(), sp0.omega_2(), sp0.omega_log()
        wc = ratio * w2
        for lv in LAM_GRID:
            sp = Spectrum.from_table(w, a * (lv / l0))
            t = Solver(spectrum=sp, mu_star=MU, w_cut=wc, n_matsubara=NM).tc()
            f = t / (C_INF * np.sqrt(sp.lambda_total()) * w2 / KB_MEV_PER_K)
            assert 0.0 < f < 1.0, f"𝔣 越界报警：{f:.4f}"     # 天花板当场报警
            L.append(np.log(sp.lambda_total())); S.append(np.log(wl / w2))
            Y.append(np.log(f)); NC.append(wc / (2 * np.pi * t * KB_MEV_PER_K))
    return map(np.array, (L, S, Y, NC))


def fit(L, S, Y):
    X = np.column_stack([np.ones(len(L)), L, S, L * L])
    c, *_ = np.linalg.lstsq(X, Y, rcond=None)
    r = Y - X @ c
    return c, float(np.sqrt(r @ r / len(Y)))


t0 = time.time()
print("=" * 96)
print(f"§1  生产工作点（w_cut = 10·ω̄₂, μ*={MU}, N={NM}）的 w_cut/(2πTc) 与到整数的距离")
print("=" * 96)
L, S, Y, NC = sweep(10.0)
c10, rms10 = fit(L, S, Y)
keys = list(RAW)
print(f"  {'谱':>9}" + "".join(f"{f'λ={v:.2f}':>9}" for v in LAM_GRID))
worst = (9.9, None)
for i, k in enumerate(keys):
    seg = NC[i * 9:(i + 1) * 9]
    print(f"  {k:>9}" + "".join(f"{v:>9.3f}" for v in seg))
    d = np.minimum(seg % 1, 1 - seg % 1)
    print(f"  {'距整数':>9}" + "".join(f"{v:>9.3f}" for v in d))
    if d.min() < worst[0]:
        worst = (d.min(), f"{k} @ λ={LAM_GRID[int(np.argmin(d))]:.2f}")
xs = sum(int(np.count_nonzero(np.diff(np.floor(NC[i * 9:(i + 1) * 9])))) for i in range(4))
print(f"  ★ 最贴边的点：{worst[1]}，距整数 {worst[0]:.4f} 格"
      f"（α 那次锚点是 0.1032 格 ⟹ {'同样危险' if worst[0] < 0.15 else '比 α 那次安全'}）")
print(f"  ★ 四条 λ 扫描合计穿越整数边界 {xs} 次 ⟹ 阶梯在拟合里被**平均**，"
      "与 α 的单点取值不是同一处境")
print(f"  生产拟合：a={c10[0]:.6f}  b={c10[1]:.6f}  c={c10[2]:.6f}  q={c10[3]:.6f}  rms={rms10:.3e}")

print()
print("=" * 96)
print("§2  ★ 判别实验：w_cut/ω̄₂ 密扫 [9.6, 10.4] 共 17 档，看 q(ratio) 光不光滑")
print("=" * 96)
print("  判据是**光滑性**，不是「q 动没动」——q 随 ω_c 平滑移动是已知物理（Tc 对 ω_c 非单调）。")
print(f"  {'ratio':>8}{'q':>12}{'c':>12}{'b':>12}{'rms':>11}{'最小距整数':>12}")
RAT = np.linspace(9.6, 10.4, 17)
QS, CS, RM = [], [], []
for r in RAT:
    L, S, Y, NC = sweep(r)
    cc, rr = fit(L, S, Y)
    QS.append(cc[3]); CS.append(cc[2]); RM.append(rr)
    dmin = float(np.min(np.minimum(NC % 1, 1 - NC % 1)))
    print(f"  {r:>8.3f}{cc[3]:>12.6f}{cc[2]:>12.6f}{cc[1]:>12.6f}{rr:>11.3e}{dmin:>12.4f}")
QS, CS, RM = np.array(QS), np.array(CS), np.array(RM)
PP_Q, PP_C = float(np.ptp(QS)), float(np.ptp(CS))   # numpy 2.x 已删除 ndarray.ptp 方法

d2 = np.abs(np.diff(QS, 2))
med = float(np.median(d2))
print()
print(f"  q 的跨度：[{QS.min():.6f}, {QS.max():.6f}]，峰峰 {PP_Q:.2e}"
      f"（生产值 {c10[3]:.6f}）")
print(f"  ★ |Δ²q|：中位 {med:.3e}，最大 {d2.max():.3e}，最大/中位 = {d2.max()/med:.2f}×"
      f"  @ ratio={RAT[int(np.argmax(d2))+1]:.3f}")
smooth = d2.max() / med < 5.0
print(f"     判决：{'光滑 ⟹ 阶梯**未**进入 q' if smooth else '★ 不光滑 ⟹ q 被阶梯污染，须换约定'}"
      f"（同一把尺子：α 那次 μ* 向 108.5×、λ 向 541.1×）")

print()
print("=" * 96)
print("§3  阶梯落在哪里：系数 vs 残差")
print("=" * 96)
print(f"  q 在整段 ratio 上的峰峰          = {PP_Q:.3e}")
print(f"  c 在整段 ratio 上的峰峰          = {PP_C:.3e}")
print(f"  拟合 rms（模型失配）中位         = {np.median(RM):.3e}")
print(f"  ⟹ q 的全部 ω_c 敏感度 {PP_Q:.1e} 是 rms {np.median(RM):.1e} 的 "
      f"{PP_Q/np.median(RM):.2f}×")
print()
print("  预算换算（CLAUDE.md 预算优先）：")
print(f"    · 仓库记录 dq/dμ* ≈ −2.3 ⟹ q 动 {PP_Q:.3e} 等价于 μ* 动 "
      f"{PP_Q/2.3:.2e}")
print(f"    · 理论侧 q 的误差棒 ±0.026、与数据 −0.525 差 1.2σ")
print(f"    ⟹ 本项对「q 对上数据」这条结论的预算 = "
      f"{'∞（不承重，不派活）' if PP_Q < 0.026 else '★ 有限，须收口'}")
print()
print("=" * 96)
print(f"完成，用时 {time.time()-t0:.0f} 秒")
print("=" * 96)
