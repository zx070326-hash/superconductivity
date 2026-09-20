# -*- coding: utf-8 -*-
"""★ 接 `_chk_main_b_gap_bias_inversion_2026-09-20.py` §4(ii)：
**(λ, S) 这两个形状参数到底够不够定 𝔣？**

背景：同阶比较后，三个形状系数 (b,c,q) 全部相合（1.40σ / 0.68σ / 0.39σ），
**只有截距 a 差 4.2σ**（数据 −1.0003 vs 理论 −0.8707，即 𝔣 上 −12.2%）。
「只有截距偏、形状系数全不偏」这个形状指向一个具体机制：
**残余形状自由度**——在固定 (λ, S) 上仍可自由变形的那部分 α²F，
它把 𝔣 整体上下平移（进 a），却几乎不改变 𝔣 对 λ 的依赖（不进 b,q）。

本脚本直接量它：构造**同一 (λ, S)、不同形状**的谱族，解 Eliashberg，看 𝔣 散多少。

族（全部解析构造，逐点校验 λ/S 命中）：
  · 双 Einstein：λ₁=xλ 于 ω₁、λ₂=(1−x)λ 于 ω₂=ρω₁，ρ 由 S 反解（x 扫）。
    这正是氢化物「低频 S 块 + 高频 H 块」的最简模型，x 就是两块的权重比。
  · Lorentzian：宽度 Γ 由 S 反解（每个 S 只有一个形状 ⟹ 一个独立点）。
  · 矩形块：α²F 在 [ω_a, ω_b] 上为常数，ω_b/ω_a 由 S 反解（另一个独立点）。

判据（预先写死，不看结果再定）：
  若「固定 (λ,S) 上 𝔣 的形状散布」≥ 观测到的 12.2%，则 a 的缺口被该机制解释；
  若 ≤ 一半，则 a 的缺口另有出处，须另找。

方向注记：双侧估计，不产生单边界判决。
严格性：numerical-evidence（float64 Eliashberg + 有限形状族）。
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
MU, WRATIO, NMATS = 0.10, 10.0, 400
W2_TARGET = 50.0                       # ω̄₂ 固定（μ*>0 只依赖 w_cut/ω̄₂，此值不影响判决）
LAM_GRID = np.exp(np.linspace(np.log(0.69), np.log(2.59), 9))
A_GAP = -0.1296                        # 观测到的 a 缺口（数据 − 理论），见上游脚本 §1


def bisect(f, lo, hi, tol=1e-14, it=200):
    flo, fhi = f(lo), f(hi)
    assert flo * fhi < 0, f"区间不变号：f({lo})={flo:.4g}, f({hi})={fhi:.4g}"
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if hi - lo < tol * max(1.0, abs(mid)):
            break
    return 0.5 * (lo + hi)


# ------------------------------------------------------------ 形状族构造
def sp_twopeak(lam, S, x):
    """双 Einstein：权重比 x = λ₁/λ；ρ=ω₂/ω₁ 由 S 反解。"""
    def lnS(rho):
        return (1 - x) * np.log(rho) - 0.5 * np.log(x + (1 - x) * rho ** 2)
    rho = bisect(lambda r: lnS(r) - np.log(S), 1.0 + 1e-12, 1e6)
    w1 = W2_TARGET / np.sqrt(x + (1 - x) * rho ** 2)
    return Spectrum.two_peak(x * lam, w1, (1 - x) * lam, rho * w1)


def _table_sp(lam, w, dens):
    """按给定 α²F 形状（未归一）构造 Spectrum，并把 λ 调到目标值。"""
    sp = Spectrum.from_table(w, dens)
    return Spectrum.from_table(w, np.asarray(dens) * (lam / sp.lambda_total()))


def sp_block(lam, S):
    """矩形块：α²F 在 [w0, R·w0] 上常数；R 由 S 反解。"""
    def lnS(R):
        w = np.linspace(1.0, R, 20001)
        d = np.ones_like(w)
        wl = np.exp(np.trapezoid(d / w * np.log(w), w) / np.trapezoid(d / w, w))
        w2 = np.sqrt(np.trapezoid(d * w, w) / np.trapezoid(d / w, w))
        return np.log(wl / w2)
    R = bisect(lambda r: lnS(r) - np.log(S), 1.0 + 1e-6, 1e4)
    w = np.linspace(1.0, R, 20001)
    w = w * (W2_TARGET / np.sqrt(np.trapezoid(w, w) / np.trapezoid(1.0 / w, w)))
    return _table_sp(lam, w, np.ones_like(w))


def sp_lorentz(lam, S):
    """Lorentzian(ω₀, Γ)：Γ/ω₀ 由 S 反解。"""
    def build(g):
        w0 = 1.0
        w = np.linspace(max(1e-4, w0 - 12 * g), w0 + 12 * g, 40001)
        w = w[w > 1e-6]
        return w, 1.0 / ((w - w0) ** 2 + g ** 2)

    def lnS(g):
        w, d = build(g)
        wl = np.exp(np.trapezoid(d / w * np.log(w), w) / np.trapezoid(d / w, w))
        w2 = np.sqrt(np.trapezoid(d * w, w) / np.trapezoid(d / w, w))
        return np.log(wl / w2)
    g = bisect(lambda gg: lnS(gg) - np.log(S), 1e-4, 0.40)
    w, d = build(g)
    w = w * (W2_TARGET / np.sqrt(np.trapezoid(d * w, w) / np.trapezoid(d / w, w)))
    return _table_sp(lam, w, d)


def frak(sp):
    t = Solver(spectrum=sp, mu_star=MU, w_cut=WRATIO * sp.omega_2(),
               n_matsubara=NMATS).tc()
    assert t > 0
    return t / (C_INF * np.sqrt(sp.lambda_total()) * sp.omega_2() / KB_MEV_PER_K)


S_TARGET = 0.80          # AD75 支撑域 [0.5446, 0.8919] 的中部
SHAPES = ([(f"2峰 x={x:.2f}", lambda l, x=x: sp_twopeak(l, S_TARGET, x))
           for x in (0.15, 0.30, 0.45, 0.60, 0.75, 0.90)]
          + [("矩形块", lambda l: sp_block(l, S_TARGET)),
             ("Lorentz", lambda l: sp_lorentz(l, S_TARGET))])

print("=" * 96)
print(f"门：每个形状都必须命中目标 (λ, S)  —— S* = {S_TARGET}，λ* = 1.30")
print("=" * 96)
print(f"  {'形状':<14}{'λ':>10}{'S':>10}{'ω̄₂[meV]':>12}{'ω_log[meV]':>12}{'λ 误差':>10}{'S 误差':>10}")
for nm, mk in SHAPES:
    sp = mk(1.30)
    lv, sv = sp.lambda_total(), sp.omega_log() / sp.omega_2()
    print(f"  {nm:<14}{lv:>10.6f}{sv:>10.6f}{sp.omega_2():>12.4f}{sp.omega_log():>12.4f}"
          f"{abs(lv/1.30-1):>10.2e}{abs(sv/S_TARGET-1):>10.2e}")
    assert abs(lv / 1.30 - 1) < 1e-9, f"{nm} λ 未命中"
    assert abs(sv / S_TARGET - 1) < 2e-6, f"{nm} S 未命中：{sv:.8f}"
print("  ⟹ 八个形状全部在 (λ, S) 上逐位一致 ⟹ 它们的 𝔣 之差**纯粹**是残余形状自由度。")

print()
print("=" * 96)
print(f"结果：固定 (λ, S) 上 𝔣 的形状散布（μ*={MU}, w_cut={WRATIO:.0f}·ω̄₂, N={NMATS}）")
print("=" * 96)
hdr = "  " + f"{'形状':<14}" + "".join(f"{l:>9.2f}" for l in LAM_GRID) + f"{'a':>9}{'b':>8}{'q':>8}"
print(hdr)
recs = {}
for nm, mk in SHAPES:
    fs = np.array([frak(mk(l)) for l in LAM_GRID])
    X = np.column_stack([np.ones(len(fs)), np.log(LAM_GRID), np.log(LAM_GRID) ** 2])
    co, *_ = np.linalg.lstsq(X, np.log(fs), rcond=None)
    recs[nm] = (fs, co)
    print("  " + f"{nm:<14}" + "".join(f"{v:>9.4f}" for v in fs)
          + f"{co[0]:>9.4f}{co[1]:>8.4f}{co[2]:>8.4f}")

F = np.array([recs[nm][0] for nm, _ in SHAPES])
CO = np.array([recs[nm][1] for nm, _ in SHAPES])
print()
print("  " + f"{'跨形状 极差/中位':<14}"
      + "".join(f"{(F[:, j].max()-F[:, j].min())/np.median(F[:, j]):>8.1%} " for j in range(F.shape[1])))
for j, nm in enumerate(("a", "b", "q")):
    print(f"  {nm}: 跨形状范围 [{CO[:, j].min():+.4f}, {CO[:, j].max():+.4f}]，"
          f"极差 {CO[:, j].max()-CO[:, j].min():.4f}")

fmax = F.max()
print()
print(f"  𝔣 全局最大 = {fmax:.4f} {'<<< ⚠ 破天花板！' if fmax > 1 else '< 1（天花板报警未触发）'}")
assert fmax < 1.0, f"𝔣 > 1，必须当场停下：{fmax:.4f}"

a_spread = float(CO[:, 0].max() - CO[:, 0].min())
b_spread = float(CO[:, 1].max() - CO[:, 1].min())
q_spread = float(CO[:, 2].max() - CO[:, 2].min())

# ---------------------------------------------------------------- 标定
# 人工族不是真实材料的形状分布：它**含病态形状**（低频权重极小的双峰、极宽 Lorentz），
# 会高估散布；而我那四条氢化物谱彼此形状相近，会**低估**。
# 唯一的第三方标尺是 AD75 自己二次拟合的残差 s = 0.0470（log𝔣 上 4.7%）——
# 它就是「(λ,S) 在真实族内有多不完备」的直接测量。
# 标定链：在本族内测出 σ_coef / σ_log𝔣 的比，再乘 AD75 的 s。
S_AD75 = 0.0470
JREF = int(np.argmin(np.abs(LAM_GRID - 1.30)))
sd_logf = float(np.std(np.log(F[:, JREF]), ddof=1))
sd = CO.std(axis=0, ddof=1)
REG = [i for i, (nm, _) in enumerate(SHAPES) if nm not in ("2峰 x=0.15", "Lorentz")]
sd_reg = CO[REG].std(axis=0, ddof=1)
sd_logf_reg = float(np.std(np.log(F[REG, JREF]), ddof=1))

print()
print("=" * 96)
print("标定：把本族的散布折算成「AD75 真实族」的形状系统误差")
print("=" * 96)
print(f"  本族（全 8 形状）  σ(log𝔣@λ={LAM_GRID[JREF]:.2f}) = {sd_logf:.4f}"
      f"   σ(a)={sd[0]:.4f}  σ(b)={sd[1]:.4f}  σ(q)={sd[2]:.4f}")
print(f"  去掉两个病态形状（2峰 x=0.15 与 Lorentz，低频权重极小 ⟹ ω_log 被一个微弱低频峰拖住）：")
print(f"                     σ(log𝔣) = {sd_logf_reg:.4f}"
      f"   σ(a)={sd_reg[0]:.4f}  σ(b)={sd_reg[1]:.4f}  σ(q)={sd_reg[2]:.4f}")
print()
print(f"  第三方标尺：AD75 29 行二次拟合的残差 s = {S_AD75:.4f}")
print(f"  ⟹ 真实族的形状系统误差（按 σ_coef/σ_log𝔣 的比线性折算）：")
print(f"  {'':<8}{'全族比':>10}{'⟹ σ_sys':>11}{'':<4}{'去病态比':>10}{'⟹ σ_sys':>11}")
sys_all, sys_reg = {}, {}
for j, nm in enumerate(("a", "b", "q")):
    r1_, r2_ = sd[j] / sd_logf, sd_reg[j] / sd_logf_reg
    sys_all[nm], sys_reg[nm] = r1_ * S_AD75, r2_ * S_AD75
    print(f"  {nm:<8}{r1_:>10.3f}{sys_all[nm]:>11.4f}{'':<4}{r2_:>10.3f}{sys_reg[nm]:>11.4f}")

print()
print("=" * 96)
print("判决")
print("=" * 96)
print(f"  ★ 第 1 条（我原先的**预言被证伪**）：我预期「残余形状只进截距、不进形状系数」。")
print(f"     实测相反——极差 a {a_spread:.4f} < q {q_spread:.4f} < b {b_spread:.4f}。")
print(f"     ⟹ 「只有 a 偏」**不能**归因于残余形状自由度的特殊结构；该机制会把三个一起推。")
print()
print(f"  ★ 第 2 条（真正承重的产出，方向对我不利）：**我此前给理论侧 q 的误差棒太小。**")
print(f"     笔记里写的是 q = −0.5565 ± 0.0256(拟合) ± 0.0140(约定)，合成 ±"
      f"{float(np.hypot(0.0256, 0.0140)):.4f}。")
print(f"     但**形状系统误差**（用 AD75 自己的残差标定）是 ±{sys_reg['q']:.4f}（去病态）"
      f"~ ±{sys_all['q']:.4f}（全族），")
print(f"     即比原先合成的误差大 {sys_reg['q']/float(np.hypot(0.0256,0.0140)):.1f}–"
      f"{sys_all['q']/float(np.hypot(0.0256,0.0140)):.1f}×。")
print(f"     根因：我那四条谱**全是氢化物**、彼此形状相近（q 跨度仅 0.0127）⟹ 四谱跨度**不是**形状系统误差的估计，")
print(f"           它只测了「H3S/D3S × 谐/非谐」这一个很窄的形状邻域。")
q_tot = float(np.hypot(np.hypot(0.0256, 0.0140), sys_reg["q"]))
print(f"     ⟹ **订正后的理论预言：q = −0.5565 ± {q_tot:.4f}**（拟合 ⊕ 约定 ⊕ 形状）。")
print(f"       与数据 −0.5253 ± 0.0749 的差距从 0.39σ 变成 "
      f"{abs(-0.5253+0.5565)/float(np.hypot(0.0749, q_tot)):.2f}σ ——**结论（相合）不变，但精度声明必须下调**。")
print()
print(f"  ★ 第 3 条：`a` 的缺口 {A_GAP:+.4f} 对比形状系统误差 ±{sys_reg['a']:.4f}~±{sys_all['a']:.4f}")
za = abs(A_GAP) / sys_all["a"]
print(f"     ⟹ {za:.1f}σ_shape（用全族比）/ {abs(A_GAP)/sys_reg['a']:.1f}σ_shape（去病态）。")
if za < 2.0:
    print("     ⟹ **在形状系统误差之内，`a` 的缺口不再是一个异常**——它和三个形状系数一样，")
    print("        落在「氢化物谱 ≠ 元素/合金谱」的预期偏差里。原先的「4.2σ」只是因为分母漏了这一项。")
else:
    print("     ⟹ 即便计入形状系统误差，`a` 仍偏 ⟹ 另有机制，不得就此结案。")
print()
print("  ⚠ 承重边界（本节自己的弱点，写在前面）：")
print("    · 标定链假设 σ_coef/σ_log𝔣 这个比在人工族与真实族之间可迁移——**未验证**，")
print("      它是把 8 个解析形状的内部结构外推到 29 个真实材料。这是本节最弱的一环。")
print("    · AD75 的残差 s=0.0470 里既有形状不完备，也有实验误差与转录误差 ⟹ 它是**上界**，")
print("      方向上使 σ_sys 偏大 ⟹ 「相合」的判决因此是**保守**的（这个方向是安全的）。")
print("    · 反过来，用它来说「a 的缺口不异常」时方向就**不**安全（σ 偏大使异常变小）⟹")
print("      第 3 条只作**降级**用（把 4.2σ 降为「未判定」），不作**结案**用。")
print("=" * 96)
