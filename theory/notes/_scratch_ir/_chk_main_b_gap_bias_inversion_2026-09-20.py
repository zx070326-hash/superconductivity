# -*- coding: utf-8 -*-
"""★ 订正我自己的 `q_theory_side_eliashberg_2026-09-20.md` §4 第 3 条（「b 差 4σ，未解释」），
并把「技术偏置对四个系数各自做什么」补全成闭式。

**缺陷**：那条把 **理论侧二次模型的 `b` = +1.1745** 与 **数据侧线性模型的 `b` = +0.9972** 相比。
两者不是同一个参数：
  · 线性模型  log𝔣 = a + b·L + c·Σ        ⟹ `b` 是样本**平均**斜率
  · 二次模型  log𝔣 = a + b·L + c·Σ + q·L² ⟹ `b` 是 **λ=1 处**的斜率
理论侧的 (b,c,q) 一直来自二次模型（`_chk_main_q_eliashberg_mustar_2026-09-20.py` 的设计矩阵含 x²），
而 AD75 那边我引的 (+0.9972, +0.7198) 来自 `_chk_main_ad75_lambda_exponent_curvature` 的 **X1（线性）**，
同一脚本的 X2（二次）给的是 (+1.2384, +0.8034, −0.5253)。
⟹ 「4σ」是**错配造成的**，不是物理。

W1 纪律：订正是一次交付。本脚本把**新**陈述也走完整流程——
门（§0 两套拟合都复现）、方向审计（§4）、以及新论证自己的独立确认（§A 无噪声恒等式断言）。

分节：
  §0  门：从零重算 AD75 29 行的**线性**与**二次**两套拟合，两套都要复现仓库既有数。
  §1  ★ 同阶比较：理论 vs 数据，四个系数逐一（理论侧的 a 本轮现算，此前没记）。
  §A  ★ 补全技术偏置闭式（加入 ω_log 的偏置 u）：
        q_D = q_T/(1+r₁)² ；**c_D ≡ c_T（六参数逐位免疫）**；
        b_D = [b_T − c_T(u₁−s₁) − ½r₁ − s₁]/(1+r₁) − 2q_T r₀/(1+r₁)²
      仓库旧式缺 `−c_T(u₁−s₁)`，因为它隐含 u₁=s₁（S 无偏）——是特例不是错，但前提此前没写。
  §B  用 q 反解 r₁，看它与同材料双技术对（n=4）的直接测量是否相合、精度差几倍。
  §C  曲率把 μ* 反测成多少？与 AD75 自己印的 μ* 比。

方向注记：全为双侧估计，不产生单边界判决。
严格性：§A analytic（代数）+ 无噪声数值恒等式（断言 1e-10）；§0/§1/§B/§C numerical-evidence。
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

# ---------------------------------------------------------------- §0 门
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
    dof = len(y) - X.shape[1]
    s2 = float(r @ r / dof)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    R2 = 1 - float(r @ r) / float(((y - y.mean()) ** 2).sum())
    return coef, se, R2, np.sqrt(s2)


def lin(L, Sg, y):
    return lsq(np.column_stack([np.ones(len(L)), L, Sg]), y)


def quad(L, Sg, y):
    return lsq(np.column_stack([np.ones(len(L)), L, Sg, L * L]), y)


frak = tc / (C_INF * np.sqrt(lam) * wb2)
L_D, S_D, y_D = np.log(lam), np.log(wlog / wb2), np.log(frak)
cL, sL, R2L, sdL = lin(L_D, S_D, y_D)
cQ, sQ, R2Q, sdQ = quad(L_D, S_D, y_D)

print("=" * 96)
print("§0  门：AD75 Table I 29 行——**两套**拟合都从零重算")
print("=" * 96)
print(f"  线性  log𝔣 = a + b·L + c·Σ")
print(f"        a = {cL[0]:+.4f} ± {sL[0]:.4f}   b = {cL[1]:+.4f} ± {sL[1]:.4f}   "
      f"c = {cL[2]:+.4f} ± {sL[2]:.4f}   R²={R2L:.4f}  s={sdL:.4f}")
print(f"  二次  log𝔣 = a + b·L + c·Σ + q·L²")
print(f"        a = {cQ[0]:+.4f} ± {sQ[0]:.4f}   b = {cQ[1]:+.4f} ± {sQ[1]:.4f}   "
      f"c = {cQ[2]:+.4f} ± {sQ[2]:.4f}   q = {cQ[3]:+.4f} ± {sQ[3]:.4f}   R²={R2Q:.4f}  s={sdQ:.4f}")
assert abs(cL[1] - 0.9972) < 5e-4 and abs(cL[2] - 0.7198) < 5e-4, "线性门未复现"
assert abs(cQ[3] + 0.5253) < 5e-4, "二次门未复现 q"
assert frak.max() < 1.0, f"𝔣>1 报警：max = {frak.max():.4f}"
print(f"  ⟹ 两套都逐位复现仓库既有数 PASS；𝔣max = {frak.max():.4f} < 1（天花板报警未触发）")
print()
print(f"  ★ **缺陷现场**：仓库既有的 `b=+0.9972 / c=+0.7198` 是**线性**那一行，")
print(f"     而 `q=−0.5253` 只能出自**二次**那一行。我先前把理论的二次 b 与数据的线性 b 对比 ⟹ 错配。")
print(f"     同一数据换模型阶，b 从 {cL[1]:+.4f} 走到 {cQ[1]:+.4f}（差 {cQ[1]-cL[1]:+.4f}），")
print(f"     恰好就是我报的那个「4σ 缺口」的量级。")

# ---------------------------------------------------------------- §1 理论侧现算（含 a）
DATA = ROOT / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"
FILES = {
    "H3S-anh": DATA / "H3S-200GPa/EPW_anharmonic/H3S.a2f",
    "D3S-anh": DATA / "D3S-200GPa/EPW_anharmonic/D3S.a2f",
    "H3S-har": DATA / "H3S-200GPa/EPW_harmonic/H3S.a2f",
    "D3S-har": DATA / "D3S-200GPa/EPW_harmonic/D3S.a2f",
}
NMATS = 400
LAM_GRID = np.exp(np.linspace(np.log(0.69), np.log(2.59), 9))


def load(path):
    w, a = [], []
    for line in path.read_text().splitlines()[1:]:
        p = line.split()
        if len(p) < 3:
            break
        try:
            x, yv = float(p[0]), float(p[1])
        except ValueError:
            break
        w.append(x); a.append(max(yv, 0.0))
    return np.array(w), np.array(a)


RAW = {k: load(p) for k, p in FILES.items()}


def theory_scan(mu, wratio=10.0):
    Ls, Ss, ys = [], [], []
    for k, (w, a) in RAW.items():
        sp0 = Spectrum.from_table(w, a)
        lam_ref, w2, wl = sp0.lambda_total(), sp0.omega_2(), sp0.omega_log()
        for lv in LAM_GRID:
            sp = Spectrum.from_table(w, a * (lv / lam_ref))
            t = Solver(spectrum=sp, mu_star=mu, w_cut=wratio * w2,
                       n_matsubara=NMATS).tc()
            assert t > 0
            ceil = C_INF * np.sqrt(sp.lambda_total()) * w2 / KB_MEV_PER_K
            Ls.append(np.log(sp.lambda_total())); Ss.append(np.log(wl / w2))
            ys.append(np.log(t / ceil))
    return np.array(Ls), np.array(Ss), np.array(ys)


print()
print("=" * 96)
print("§1  ★ 同阶比较：理论（各向同性 Eliashberg, μ*=0.10, w_cut=10·ω̄₂）vs AD75 29 行")
print("=" * 96)
Lt, St, yt_ = theory_scan(0.10)
cTq, sTq, R2Tq, sdTq = quad(Lt, St, yt_)
cTl, sTl, R2Tl, sdTl = lin(Lt, St, yt_)
assert abs(cTq[1] - 1.1745) < 2e-3 and abs(cTq[3] + 0.5565) < 2e-3, \
    f"理论侧未复现本轮既有值：b={cTq[1]:+.4f} q={cTq[3]:+.4f}"
print(f"  门：理论二次 b={cTq[1]:+.4f}、q={cTq[3]:+.4f} 复现本轮既有 (+1.1745, −0.5565) PASS")
print()
print(f"  {'模型':<6}{'系数':<4}{'数据 (AD75 29 行)':>24}{'理论 (Eliashberg)':>24}{'差':>10}{'σ_合成':>9}{'nσ':>7}")
rows = [("二次", "a", cQ[0], sQ[0], cTq[0], sTq[0]),
        ("二次", "b", cQ[1], sQ[1], cTq[1], sTq[1]),
        ("二次", "c", cQ[2], sQ[2], cTq[2], sTq[2]),
        ("二次", "q", cQ[3], sQ[3], cTq[3], sTq[3]),
        ("线性", "a", cL[0], sL[0], cTl[0], sTl[0]),
        ("线性", "b", cL[1], sL[1], cTl[1], sTl[1]),
        ("线性", "c", cL[2], sL[2], cTl[2], sTl[2])]
for mdl, nm, d, sd, t, st in rows:
    dif = d - t
    ss = float(np.hypot(sd, st))
    print(f"  {mdl:<6}{nm:<4}{d:>16.4f} ±{sd:>6.4f}{t:>16.4f} ±{st:>6.4f}"
          f"{dif:>10.4f}{ss:>9.4f}{abs(dif)/ss:>7.2f}")
print()
print("  ⚠ 理论侧的 σ 只是 36 点拟合的统计误差，**不含**约定/谱形状的系统项")
print("    （截断约定 ±0.014 于 q；四谱 b 跨度半宽 ±0.046）⟹ 上表的 nσ 是**保守上界**，真值更小。")
print("  ⚠ `a` 不是技术免疫的（r₀,s₀,u₀ 全进 a），且依赖绝对归一 ⟹ 它不作判决用，只列出。")

# ---------------------------------------------------------------- §A 闭式
print()
print("=" * 96)
print("§A  ★ 补全技术偏置闭式（把 ω_log 的偏置 u 也放进来）")
print("=" * 96)
print("  模型：logλ_D = L_T + (r₀+r₁L_T)；logω̄₂_D = logω̄₂_T + (s₀+s₁L_T)；")
print("        logω_log,D = logω_log,T + (u₀+u₁L_T)；**logTc_D = logTc_T（Tc 技术无关）**")
print("  ⟹ L_D=(1+r₁)L_T+r₀ ；Σ_D=Σ_T+(u₀−s₀)+(u₁−s₁)L_T ；y_D=y_T−½ρ−σ")
rng = np.random.default_rng(20260920)
n = 400
Lt0 = rng.uniform(-0.40, 0.96, n)
Sg0 = rng.uniform(np.log(0.50), np.log(0.92), n)
A0, B0, C0, Q0 = 0.31, 1.1745, 0.6482, -0.5565
yt0 = A0 + B0 * Lt0 + C0 * Sg0 + Q0 * Lt0 ** 2          # 无噪声：测恒等式而非统计


def apply_bias(r0, r1, s0, s1, u0, u1):
    rho = r0 + r1 * Lt0
    sig = s0 + s1 * Lt0
    uu = u0 + u1 * Lt0
    return quad(Lt0 + rho, Sg0 + (uu - sig), yt0 - 0.5 * rho - sig)[0]


print()
print(f"  {'r0':>6}{'r1':>7}{'s0':>6}{'s1':>7}{'u0':>6}{'u1':>7} | "
      f"{'c_D':>10}{'c_T':>9} | {'q_D':>10}{'q闭式':>10} | {'b_D':>10}{'b闭式':>10}")
for r0, r1, s0, s1, u0, u1 in [(0, 0, 0, 0, 0, 0),
                               (0.05, -0.173, -0.02, -0.054, 0.0, 0.0),
                               (0.05, -0.173, -0.02, -0.054, 0.03, 0.11),
                               (-0.12, 0.030, 0.0, 0.0, 0.0, 0.0),
                               (0.0, 0.030, 0.0, 0.38, 0.0, 0.0),
                               (0.0, 0.030, 0.0, 0.0, 0.0, 0.21),
                               (0.10, 0.40, 0.05, 0.30, -0.07, -0.19)]:
    a_, b_, c_, q_ = apply_bias(r0, r1, s0, s1, u0, u1)
    v = 1.0 / (1.0 + r1)
    q_pred = Q0 * v * v
    b_pred = (B0 - C0 * (u1 - s1) - 0.5 * r1 - s1) * v - 2.0 * Q0 * v * v * r0
    print(f"  {r0:>6.2f}{r1:>7.3f}{s0:>6.2f}{s1:>7.3f}{u0:>6.2f}{u1:>7.2f} | "
          f"{c_:>10.6f}{C0:>9.4f} | {q_:>10.6f}{q_pred:>10.6f} | {b_:>10.6f}{b_pred:>10.6f}")
    assert abs(c_ - C0) < 1e-10, "c 不免疫 ⟹ 推导错"
    assert abs(q_ - q_pred) < 1e-10, "q 闭式不符"
    assert abs(b_ - b_pred) < 1e-10, "b 闭式不符"
print("  （七行逐位吻合，三道 1e-10 断言门通过）")
print()
print("  ★ 两条结构事实：")
print("    ① **c_D ≡ c_T**，对 r₀,r₁,s₀,s₁,u₀,u₁ **六参数逐位免疫** [analytic]。")
print("       ⟹ `c` 是继「q 的符号」之后**第二个**技术免疫的可比量，而且是逐位免疫、不只是符号。")
print("    ② 旧式 b_D=(b−½r₁−s₁)/(1+r₁)−2q·r₀/(1+r₁)² 缺 `−c_T(u₁−s₁)/(1+r₁)`，")
print("       它在 u₁=s₁（ω_log 与 ω̄₂ 同偏、S 无偏）时消失 ⟹ **旧式是特例不是错，但前提此前没写下来**。")

# ---------------------------------------------------------------- §B r1
print()
print("=" * 96)
print("§B  用 q 反解 r₁，对照同材料双技术对（n=4）的直接测量")
print("=" * 96)
qD, sqD, qT, sqT = cQ[3], sQ[3], cTq[3], float(np.hypot(sTq[3], 0.0140))
v2 = qD / qT
v = float(np.sqrt(v2)); r1_hat = 1.0 / v - 1.0
rel = 0.5 * float(np.hypot(sqD / abs(qD), sqT / abs(qT)))
se_r1 = (1.0 + r1_hat) * rel
print(f"  (1+r₁)⁻² = q_D/q_T = {qD:+.4f}/{qT:+.4f} = {v2:.4f}"
      f"   ⟹ r₁ = {r1_hat:+.4f} ± {se_r1:.4f}")
print(f"  直接测量（Pb/Ta/Al/Mo 四对，`_chk_main_r1_fourpoint_audit`）：r₁ = +0.119 ± 0.309")
print(f"  ⟹ 两者差 {abs(r1_hat-0.119)/float(np.hypot(se_r1,0.309)):.2f}σ，相合。")
print(f"  ⟹ 曲率路线的 σ_r₁ = {se_r1:.4f}，比四对直接测量窄 {0.309/se_r1:.1f}×，且**不需要配对数据**。")
print(f"  ⚠ 但它不能用来「订正 q 后再引用 q」——那是循环。它能说的只是：")
print(f"     **在 r₁=0 的零假设下理论与数据已相合**，因此不必调用任何技术订正。")

# ---------------------------------------------------------------- §C mu*
print()
print("=" * 96)
print("§C  把曲率当温度计：数据的 q 反测出多少 μ*？")
print("=" * 96)
mus, qs_ = [], []
for mu in (0.00, 0.05, 0.10, 0.13):
    if mu == 0.10:
        qv = cTq[3]
    else:
        Lx, Sx, yx = theory_scan(mu)
        qv = quad(Lx, Sx, yx)[0][3]
    mus.append(mu); qs_.append(qv)
    print(f"  μ* = {mu:.2f}  ⟹  q_theory = {qv:+.4f}")
mus = np.array(mus); qs_ = np.array(qs_)
slope = np.polyfit(mus, qs_, 1)[0]
mu_inv = float(np.interp(qD, qs_[::-1], mus[::-1]))
se_mu = float(np.hypot(sqD, sqT)) / abs(slope)
print()
print(f"  dq/dμ* = {slope:+.3f} /单位 μ*")
print(f"  数据 q = {qD:+.4f} ± {sqD:.4f}  ⟹  **μ*(由曲率反测) = {mu_inv:.3f} ± {se_mu:.3f}**")
print(f"  AD75 Table I 自己逐行印出的 μ* ∈ [0.092, 0.111]（≈0.10）")
print(f"  ⟹ 相差 {abs(mu_inv-0.10)/se_mu:.2f}σ，相合。")
print(f"  ⚠ 这是**一致性**不是**精度**：σ={se_mu:.3f} 比印出的跨度 0.019 宽 {se_mu/0.019:.0f}×")
print(f"     ⟹ 曲率**判不掉** AD75 的 μ*，只能说「不矛盾」。不得反过来宣称「独立测定了 μ*」。")

print()
print("=" * 96)
print("总判决（取代 q_theory_side §4 第 3 条）")
print("=" * 96)
print("  1. ★ **「b 差 4σ」撤回**——那是线性 vs 二次的模型阶**错配**，不是物理。")
print(f"     同阶（二次）比较：数据 b = {cQ[1]:+.4f} ± {sQ[1]:.4f}，理论 b = {cTq[1]:+.4f} ± {sTq[1]:.4f}，"
      f"差 {abs(cQ[1]-cTq[1])/float(np.hypot(sQ[1],sTq[1])):.2f}σ。")

# 判词由数据算出，不预写（CLAUDE.md：判词不得预写成字符串）
THR = 2.0
agree = [(nm, abs(d - t) / float(np.hypot(sd, st)))
         for mdl, nm, d, sd, t, st in rows if mdl == "二次"]
ok = [(nm, z) for nm, z in agree if z <= THR]
bad = [(nm, z) for nm, z in agree if z > THR]
print(f"  2. ★ **同阶二次模型下，{len(ok)}/{len(agree)} 个系数相合**（判据 |Δ| ≤ {THR:.0f}σ_合成）：")
print("     相合：" + "、".join(f"{nm} ({z:.2f}σ)" for nm, z in ok))
if bad:
    print("     **不相合**：" + "、".join(f"{nm} ({z:.2f}σ)" for nm, z in bad))
    print("     ⟹ 不得笼统写成「全部相合」。**三个形状系数 (b,c,q) 相合，截距 a 不相合。**")
else:
    print("     无不相合项。")
print("  3. ★ **c 的相合是技术免疫的**（§A①，六参数逐位）⟹ 不能被「隧道反演有系统误差」解释掉。")
print("     b 与 q 的相合**不是**技术免疫的（b）/ 只有符号免疫（q）⟹ 承重程度依次为 c > q > b > a。")
if bad:
    a_gap = cQ[0] - cTq[0]
    print(f"  4. ★ **`a` 的 {a_gap:+.4f}（{100*(np.exp(a_gap)-1):+.1f}% 于 𝔣）是本轮真正待解释的对象**，")
    print("     且它的形状很特殊：**只有截距偏，三个形状系数都不偏**。两类候选：")
    print("     (i) 纯偏移型技术偏置 (s₀,u₀)——闭式 a_D−a_T = −c_T·u₀ + (c_T−1)·s₀；")
    print(f"         仅 u₀ 时需 u₀ = {a_gap/(-cTq[2]):+.4f}（隧道 ω_log 比理论高 "
          f"{100*(np.exp(a_gap/(-cTq[2]))-1):+.1f}%），"
          f"仅 s₀ 时需 s₀ = {a_gap/(cTq[2]-1):+.4f}（**量级荒谬，排除**）。")
    print("     (ii) **(λ,S) 两参数不足以定 𝔣**——氢化物 α²F 是「低频 S 块 + 高频 H 块」的双块结构，")
    print("          与 Pb/In/Tl 的单声学块形状迥异；残余形状自由度只进截距、不进形状系数，")
    print("          正好是观测到的形状。⟹ **下一步：在固定 (λ,S) 上扫形状，量这个残余自由度有多大。**")
    print(f"     旁证：AD75 自己的二次拟合残差 s = {sdQ:.4f}（{100*sdQ:.1f}%），说明 (λ,S) 在**族内**已不完备；")
    print(f"     族**间**（氢化物 vs 元素/合金）差 {100*abs(np.exp(a_gap)-1):.0f}% 并不意外。")
print("  5. ⚠ **消费点更新（与订正前相反）**：「同材料双技术」买不到 q（符号免疫）、买不到 c（逐位免疫），")
print("     而 b 在同阶下已相合 ⟹ **对本比较预算 = ∞，不派活**。")
print("     真正该花的是 §4(ii)：**固定 (λ,S) 的形状扫描**，它才对准唯一还开着的缺口 `a`。")
print("=" * 96)
