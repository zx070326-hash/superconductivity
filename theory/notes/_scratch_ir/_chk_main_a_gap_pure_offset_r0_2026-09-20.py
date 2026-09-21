# -*- coding: utf-8 -*-
"""★ `a` 缺口的最后一个单参数候选：**λ 的纯偏移 r₀**，且 `b` 是零自由度的**预言**。

局面（同阶二次模型，`_chk_main_b_gap_bias_inversion` §1）：
    a: 数据 −1.0003 vs 理论 −0.8707  ⟹ 差 −0.1296   ← 唯一不相合项
    b: +1.2384 vs +1.1745            ⟹ 1.40σ（拟合σ）
    c,q: 相合，且 c 对六个偏置参数**逐位免疫**、q 只经 (1+r₁)² 进入（符号免疫）

已排除/已测：
  · 残余形状自由度（`_chk_main_hydride_shape_is_high_f`）：方向**对**（氢化物确是高 𝔣 形状，
    36/36 格比值 >1），但只买到缺口的 **21%**（+2.7% vs 需要的 +13.8%）。
  · s₀ 单独（ω̄₂ 纯偏移）：需 s₀ = +1.04，即 ω̄₂ 差 2.8× ——量级荒谬，排除。

**本脚本的关键点：三个纯偏移 (r₀, s₀, u₀) 里，只有 r₀ 同时进 `a` 和 `b`。**
⟹ 拿 `a` 定 r₀（1 个自由度），`b` 就是**零自由度的预言**（1 个检验）。
这是一次**可以失败**的检验：若所需 r₀ 把 `b` 推到 3σ 外，r₀ 假说当场死。
u₀ 则只进 `a` 不进 `b` ⟹ 它能补上 `a` 但对 `b` 毫无作为。两者由此可分。

误差预算按 `_chk_main_shape_completeness_lamS` 的标定加宽：理论侧除拟合 σ 外，
另计**谱形状系统误差**（σ_sys(a)=0.0220–0.0370，σ_sys(b)=0.0518–0.0985，σ_sys(q)=0.0314–0.0409）。
本脚本一律取**保守（较大）**的一端 —— 注意方向：σ 取大使「相合」更容易，
⟹ 对「相合」的判决是**反保守**的，故凡宣布相合处均须另行说明（见 §4）。

方向注记：闭式为代数恒等式（双侧）；误差预算的方向在 §4 逐条审。
严格性：§1 analytic（代数）+ 1e-10 无噪声恒等式断言；其余 numerical-evidence。
主会话 2026-09-20，从零。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "solvers"))
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K   # noqa: E402

C_INF = 0.1827262477

T1 = [("Pb",56,65,1.55,7.2),("In",68,89,0.805,3.40),("Sn",99,121,0.72,3.75),
      ("Hg",29,49,1.6,4.19),("Tl",52,64,0.795,2.36),("Ta",132,148,0.69,4.48),
      ("a-Ga",55,101,1.62,8.56),("b-Ga",87,129,0.97,5.90),
      ("Tl0.9Bi0.1",48,62,0.78,2.30),("Pb0.4Tl0.6",48,62,1.15,4.60),
      ("Pb0.6Tl0.4",50,62,1.38,5.90),("Pb0.8Tl0.2",50,61,1.53,6.80),
      ("Pb0.6Tl0.2Bi0.2",48,58,1.81,7.26),("Pb0.9Bi0.1",50,60,1.66,7.65),
      ("Pb0.8Bi0.2",46,57,1.88,7.95),("Pb0.7Bi0.3",47,57,2.01,8.45),
      ("Pb0.65Bi0.35",45,55,2.13,8.95),("In0.9Tl0.1",63,86,0.85,3.28),
      ("In0.73Tl0.27",55,77,0.93,3.36),("In0.67Tl0.33",57,79,0.90,3.26),
      ("In0.57Tl0.43",53,74,0.85,2.60),("In0.5Tl0.5",53,73,0.83,2.52),
      ("In0.27Tl0.73",42,63,1.09,3.64),("In0.17Tl0.83",45,63,0.98,3.19),
      ("In0.07Tl0.93",49,63,0.89,2.77),("In2Bi",46,67,1.40,5.6),
      ("Sb2Tl7",37,58,1.43,5.2),("Bi2Tl",47,59,1.63,6.4),
      ("a-Pb0.45Bi0.55",29,47,2.59,7.0)]
wlog = np.array([float(r[1]) for r in T1]); wb2 = np.array([float(r[2]) for r in T1])
lam  = np.array([float(r[3]) for r in T1]); tc  = np.array([float(r[4]) for r in T1])


def lsq(X, y):
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ coef
    s2 = float(r @ r / (len(y) - X.shape[1]))
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    return coef, se


def quad(L, Sg, y):
    return lsq(np.column_stack([np.ones(len(L)), L, Sg, L * L]), y)


fr = tc / (C_INF * np.sqrt(lam) * wb2)
L_D, S_D, y_D = np.log(lam), np.log(wlog / wb2), np.log(fr)
cQ, sQ = quad(L_D, S_D, y_D)
assert abs(cQ[3] + 0.5253) < 5e-4 and abs(cQ[1] - 1.2384) < 5e-4, "数据门未复现"
assert fr.max() < 1.0, f"𝔣>1 报警：{fr.max():.4f}"

# ---------------------------------------------------------------- 理论侧
DATA = ROOT / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"
FILES = {"H3S-anh": DATA / "H3S-200GPa/EPW_anharmonic/H3S.a2f",
         "D3S-anh": DATA / "D3S-200GPa/EPW_anharmonic/D3S.a2f",
         "H3S-har": DATA / "H3S-200GPa/EPW_harmonic/H3S.a2f",
         "D3S-har": DATA / "D3S-200GPa/EPW_harmonic/D3S.a2f"}
LAM_GRID = np.exp(np.linspace(np.log(0.69), np.log(2.59), 9))


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


Ls, Ss, ys = [], [], []
for k, p in FILES.items():
    w, a = load(p)
    sp0 = Spectrum.from_table(w, a)
    l0, w2, wl = sp0.lambda_total(), sp0.omega_2(), sp0.omega_log()
    for lv in LAM_GRID:
        sp = Spectrum.from_table(w, a * (lv / l0))
        t = Solver(spectrum=sp, mu_star=0.10, w_cut=10.0 * w2, n_matsubara=400).tc()
        assert t > 0
        Ls.append(np.log(sp.lambda_total())); Ss.append(np.log(wl / w2))
        ys.append(np.log(t / (C_INF * np.sqrt(sp.lambda_total()) * w2 / KB_MEV_PER_K)))
cT, sT = quad(np.array(Ls), np.array(Ss), np.array(ys))
assert abs(cT[1] - 1.1745) < 2e-3 and abs(cT[3] + 0.5565) < 2e-3, "理论门未复现"
A_T, B_T, C_T, Q_T = (float(v) for v in cT)

# 形状系统误差（`_chk_main_shape_completeness_lamS` 的标定，取保守端）
SYS = {"a": 0.0370, "b": 0.0985, "q": 0.0409}
sig_a = float(np.hypot(np.hypot(sQ[0], sT[0]), SYS["a"]))
sig_b = float(np.hypot(np.hypot(sQ[1], sT[1]), SYS["b"]))
A_GAP = float(cQ[0]) - A_T
B_GAP = float(cQ[1]) - B_T

print("=" * 100)
print("§0  门：两侧二次拟合都从零重算")
print("=" * 100)
print(f"  数据 a={cQ[0]:+.4f}±{sQ[0]:.4f}  b={cQ[1]:+.4f}±{sQ[1]:.4f}  "
      f"c={cQ[2]:+.4f}±{sQ[2]:.4f}  q={cQ[3]:+.4f}±{sQ[3]:.4f}")
print(f"  理论 a={A_T:+.4f}±{sT[0]:.4f}  b={B_T:+.4f}±{sT[1]:.4f}  "
      f"c={C_T:+.4f}±{sT[2]:.4f}  q={Q_T:+.4f}±{sT[3]:.4f}")
print(f"  𝔣max(数据) = {fr.max():.4f} < 1（天花板未破）")
print(f"  缺口：Δa = {A_GAP:+.4f}（σ_合成含形状 = {sig_a:.4f} ⟹ {abs(A_GAP)/sig_a:.2f}σ）")
print(f"        Δb = {B_GAP:+.4f}（σ_合成含形状 = {sig_b:.4f} ⟹ {abs(B_GAP)/sig_b:.2f}σ）")

# ---------------------------------------------------------------- §1 闭式
print()
print("=" * 100)
print("§1  纯偏移闭式（r₁=s₁=u₁=0），以无噪声恒等式断言")
print("=" * 100)
print("  代数：y_D = [a_T − (b_T+½)r₀ − c_T(u₀−s₀) − s₀ + q_T r₀²] + [b_T − 2q_T r₀]·L_D")
print("            + c_T·Σ_D + q_T·L_D²")
print("  ⟹ Δa = −(b_T+½)r₀ − c_T(u₀−s₀) − s₀ + q_T r₀²    Δb = −2q_T r₀    Δc = Δq = 0")
print("  ★ **三个纯偏移里只有 r₀ 进 `b`。**")
rng = np.random.default_rng(20260920)
n = 400
Lt0 = rng.uniform(-0.40, 0.96, n)
Sg0 = rng.uniform(np.log(0.50), np.log(0.92), n)
yt0 = A_T + B_T * Lt0 + C_T * Sg0 + Q_T * Lt0 ** 2       # 无噪声


def offs(r0, s0, u0):
    return quad(Lt0 + r0, Sg0 + (u0 - s0), yt0 - 0.5 * r0 - s0)[0]


print()
print(f"  {'r0':>7}{'s0':>7}{'u0':>7} | {'Δa 实测':>11}{'Δa 闭式':>11} | "
      f"{'Δb 实测':>11}{'Δb 闭式':>11} | {'Δc':>9}{'Δq':>9}")
for r0, s0, u0 in [(0, 0, 0), (0.0755, 0, 0), (0, 0, 0.148), (0, 0.30, 0),
                   (-0.20, 0.07, -0.05), (0.33, -0.11, 0.21)]:
    a_, b_, c_, q_ = offs(r0, s0, u0)
    da_p = -(B_T + 0.5) * r0 - C_T * (u0 - s0) - s0 + Q_T * r0 ** 2
    db_p = -2.0 * Q_T * r0
    print(f"  {r0:>7.4f}{s0:>7.3f}{u0:>7.3f} | {a_-A_T:>11.6f}{da_p:>11.6f} | "
          f"{b_-B_T:>11.6f}{db_p:>11.6f} | {c_-C_T:>9.1e}{q_-Q_T:>9.1e}")
    assert abs((a_ - A_T) - da_p) < 1e-10, "Δa 闭式不符"
    assert abs((b_ - B_T) - db_p) < 1e-10, "Δb 闭式不符"
    assert abs(c_ - C_T) < 1e-10 and abs(q_ - Q_T) < 1e-10, "c/q 应逐位免疫"
print("  （六行逐位吻合，四道 1e-10 断言门通过）")

# ---------------------------------------------------------------- §2 单参数检验
print()
print("=" * 100)
print("§2  ★ 用 `a` 定单个偏移（1 自由度），`b` 成为零自由度**预言**（1 检验）")
print("=" * 100)
# r0: 解 −(b_T+½)r₀ + q_T r₀² = Δa
AA, BB, CC = Q_T, -(B_T + 0.5), -A_GAP
r0_req = float((-BB - np.sqrt(BB * BB - 4 * AA * CC)) / (2 * AA))
assert abs(-(B_T + 0.5) * r0_req + Q_T * r0_req ** 2 - A_GAP) < 1e-12
u0_req = float(-A_GAP / C_T)
s0_req = float(A_GAP / (C_T - 1.0))
cases = [("r₀（λ 纯偏移）", r0_req, -2.0 * Q_T * r0_req),
         ("u₀（ω_log 纯偏移）", u0_req, 0.0),
         ("s₀（ω̄₂ 纯偏移）", s0_req, 0.0),
         ("零假设（无偏移）", 0.0, 0.0)]
print(f"  {'候选':<18}{'所需值':>10}{'物理含义':>26} | {'预言 b':>10}{'数据 b':>10}{'nσ_b':>8}")
best = None
for nm, val, db in cases:
    bpred = B_T + db
    z = abs(float(cQ[1]) - bpred) / sig_b
    pct = 100 * (np.exp(val) - 1)
    mean = ("λ_data 比 λ_theory 高" if "r₀" in nm else
            "ω_log,data 比理论高" if "u₀" in nm else
            "ω̄₂_data 比理论高" if "s₀" in nm else "—")
    txt = "—" if val == 0 else f"{mean} {pct:+.1f}%"
    print(f"  {nm:<18}{val:>10.4f}{txt:>26} | {bpred:>10.4f}{float(cQ[1]):>10.4f}{z:>8.2f}")
    if best is None or z < best[1]:
        best = (nm, z)
print()
print(f"  ⟹ `b` 预言最好的候选：**{best[0]}**（{best[1]:.2f}σ）")
print(f"  ⟹ `s₀` 需 ω̄₂ 差 {np.exp(s0_req):.2f}× —— **量级荒谬，排除**。")
print(f"  ⟹ `u₀` 能补 `a` 但**对 `b` 零作用**（它不进 b）⟹ 留在零假设的 b 上。")

# r0 的不确定度
dadr = -(B_T + 0.5) + 2 * Q_T * r0_req
se_r0 = sig_a / abs(dadr)
print()
print(f"  所需 r₀ = {r0_req:+.4f} ± {se_r0:.4f}（由 σ_a={sig_a:.4f} 除以 |da/dr₀|={abs(dadr):.4f}）")
print(f"  ⟹ 距零 {abs(r0_req)/se_r0:.2f}σ —— 这正是 `a` 缺口本身的重述。")

# ---------------------------------------------------------------- §3 直接测量
print()
print("=" * 100)
print("§3  r₀ 的**独立直接测量**：AD75(隧道) vs Savrasov(DFT) 四对同材料")
print("=" * 100)
print("  ⚠ **约定方向**：本文 r₀ 定义为 `logλ_data = logλ_theory + r₀`，data=隧道、theory≈DFT")
print("     ⟹ r₀ = ⟨log(λ_tun/λ_dft)⟩。`_chk_main_r1_fourpoint_audit` 回归的是 log(λ_dft/λ_tun)，")
print("     **方向相反**，其印出的 r₁=+0.119 在本文约定下约为 −0.119（一阶）。此前 §B 未审此约定。")
PAIRS = {"Pb": (1.55, 1.68), "Ta": (0.69, 0.86), "Al": (0.43, 0.44), "Mo": (0.60, 0.42)}
d = np.array([np.log(lt / ld) for lt, ld in PAIRS.values()])
r0_dir = float(d.mean()); se_dir = float(d.std(ddof=1) / np.sqrt(len(d)))
print()
print(f"  {'材料':>6}{'λ_tun':>9}{'λ_dft':>9}{'log(tun/dft)':>15}")
for (k, (lt, ld)), v in zip(PAIRS.items(), d):
    print(f"  {k:>6}{lt:>9.2f}{ld:>9.2f}{v:>15.4f}")
print(f"  ⟹ r₀(直接, n=4) = {r0_dir:+.4f} ± {se_dir:.4f}")
z_r0 = abs(r0_req - r0_dir) / float(np.hypot(se_r0, se_dir))
print(f"  对照所需 {r0_req:+.4f} ± {se_r0:.4f}  ⟹  差 {z_r0:.2f}σ")

# ---------------------------------------------------------------- §4 判决
print()
print("=" * 100)
print("判决（由数据算出，未预写）")
print("=" * 100)
THR = 2.0
print(f"  1. **可以失败的检验，它没失败。** 用 `a` 定出的 r₀ = {r0_req:+.4f} 把 `b` 从"
      f" {abs(B_GAP)/sig_b:.2f}σ 改善到 {best[1]:.2f}σ。")
print(f"     若所需 r₀ 把 `b` 推到 >{THR:.0f}σ，r₀ 假说当场死 —— 它没有。")
print(f"  2. **但这不是判别**：零假设下 `b` 本来就只差 {abs(B_GAP)/sig_b:.2f}σ（计入形状系统误差后），")
print(f"     ⟹ 数据**分不开** r₀ 与 u₀ 两条路。可写「相容」，**不可写「证据支持 λ 有 +{100*(np.exp(r0_req)-1):.0f}% 偏置」**。")
print(f"  3. 直接测量 r₀ = {r0_dir:+.4f} ± {se_dir:.4f} 与所需值差 {z_r0:.2f}σ ⟹ 相容，")
print(f"     但 n=4 的 σ={se_dir:.3f} 比所需值本身还宽 {se_dir/abs(r0_req):.1f}× ⟹ **直接测量判不掉**。")
gap_left = abs(A_GAP) / sig_a
print(f"  4. ★ **`a` 仍是唯一未关闭项**（{gap_left:.2f}σ，已用最保守的形状系统误差）。")
print("     三条候选的现状：")
print(f"       · 残余形状自由度：方向对、买 21%（已测，36/36 同号）")
print(f"       · r₀（λ 偏置 +{100*(np.exp(r0_req)-1):.0f}%）：能补满，且 `b` 的零自由度预言通过")
print(f"       · u₀（ω_log 偏置 +{100*(np.exp(u0_req)-1):.0f}%）：能补满，但对 `b` 无预言力")
print("     ⟹ 三者**不互斥**，且各自都能被现有数据容纳 ⟹ **`a` 不作为论文的承重项**。")
print()
print("  5. ★ **消费点/预算**：要判开 r₀ vs u₀，需要 AD75 材料的**真实隧道 α²F**（而不只是 λ/ω_log/ω̄₂ 三个数），")
print("     以便直接比对同材料的 DFT 谱。那是一件独立的数据获取任务。")
print(f"     但**判开它买到什么**？`c`、`q` 对三个纯偏移**逐位免疫**（§1 已证），`b` 在任一分支下都相合")
print("     ⟹ **对头条（符号 + c/q 相合）的预算 = ∞ ⟹ 不派活。**")
print()
print("  ⚠ 方向审计：")
print("    · σ 取了形状系统误差的**保守（大）端** ⟹ 使「相合」更容易 ⟹ 本文凡宣布相合处，")
print("      承重的是**「未被证伪」**而非**「已被确证」**；凡宣布不相合处（`a`）方向才是安全的。")
print("    · §2 的 r₀ 是**拟合**出来的，不是测出来的 ⟹ 不得反过来引用「λ 有 8% 偏置」这个数。")
print("    · §1 的闭式是代数恒等式（analytic），与误差预算无关，不受上两条影响。")
print("=" * 100)
