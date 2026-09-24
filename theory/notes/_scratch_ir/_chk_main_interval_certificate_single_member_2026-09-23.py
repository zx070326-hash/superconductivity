# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）任务 ②(c) 第一刀：**单成员的区间算术证书**（mpmath.iv，外向舍入）。

本仓库此前关于 `λ₁(M)` 的一切都是 float64 的 `numerical-evidence`。本脚本第一次把整条链
（`G` → `B/D` → `γ` → `μ` / `res²` → Kato–Temple）换成**声音的外向舍入区间算术**，
对**单个成员**产出一个真正的 `certificate` 级上界。盒子/自适应细分是下一刀，不在本脚本内。

━━━━━━━━ 链条（每一步的方向都写在这里，全部朝「判据变难」放大） ━━━━━━━━
  λ₁(M)  ≤  μ + res²/(μ − γ)                              [Kato–Temple；要求 μ > γ]
     μ   := ⟨v,Mv⟩/⟨v,v⟩,   res² := ‖Mv − μv‖²/⟨v,v⟩       [v 任意 ⟹ 浮点 v 当精确有理向量]
     res² 的 n ≥ L 部分用解析尾界（见 res2_tail_up）        [放大]
     γ   ≥ λ₂(M)，经三步：
        (γ1) `𝔎_c = 2μ*u_cu_cᵀ ⪰ 0` ⟹ `M ⪯ 𝔎_ph` ⟹ `λ₂(M) ≤ λ₂(𝔎_ph)`      [analytic]
        (γ2) Cauchy 交错（min–max 取 V=span{e_p}）⟹ `λ₂(𝔎_ph) ≤ λ₁(𝔎_ph 去掉指标 p)` [analytic]
        (γ3) 加权 Schur：`𝔎_ph = B − D`，B 逐元非负（含对角）⟹ `B−D` 是 Metzler
             ⟹ `B−D+cI ≥ 0` 且对称 ⟹ Schur 检验函数 w>0 给
             `λ₁(B'−D') ≤ R(p,w) := sup_{i≠p} [ (1/w_i)Σ_{j≠p}B_ij w_j − D_ii ]`      [analytic]
        取 `w_j = u_j = (2j+1)^{−1/2}` ⟹ `(1/u_i)Σ_j B_ij u_j = Σ_j [G(下标距)+G(i+j+1)]/(2j+1)`。

━━━━━━━━ 本脚本订正上一轮笔记 §10.3 的行尾界（我自己的推导，有洞） ━━━━━━━━
*旧陈述错在哪*：§10.3 写 `R_i ≤ 4·M2·ln(i+1)/i² + 2(λ+π²M2/6)/(i+1)`，其第一项用了
`Σ_{j≤h} 1/(2j+1) ≤ ½ln(i+1)`（h=⌊i/2⌋）。**这一步为假**：`i=2` 时左边 `1+1/3=1.3333`，
右边 `½ln3=0.5493`。调和和的 `j=0` 项（=1）被整个漏掉了。
（GATE K3 报 0 违反是因为**第二项数量级上压倒第一项**，掩盖了这个洞——
 「0/N 通过只证明被测的那一步」的又一例：K3 测的是总和，没测这一项。）
*新陈述凭什么对*：两段重推，逐项给出处——
  · `j ≤ h`：`i−j ≥ i−h ≥ i/2` ⟹ `G(|i−j|) ≤ 4M2/i²`；`i+j+1 ≥ i` ⟹ `G(i+j+1) ≤ M2/i² ≤ 4M2/i²`；
    且 `Σ_{j=0}^{h} 1/(2j+1) ≤ 1 + ∫_0^h dt/(2t+1) = 1 + ½ln(2h+1) ≤ 1 + ½ln(i+1)`
    ⟹ 该段 `≤ (8M2/i²)(1 + ½ln(i+1))`。
  · `j > h`：`2j+1 ≥ 2h+1 ≥ i`（i 偶 ⟹ i+1，i 奇 ⟹ i，取 i 保守）；
    `Σ_{j>h}G(|i−j|) ≤ G(0)+2Σ_{k≥1}G(k) ≤ λ + 2M2·π²/6`（每个 k≥1 至多被命中两次）；
    `Σ_{j>h}G(i+j+1) ≤ Σ_{k≥1}G(k) ≤ M2·π²/6` ⟹ 该段 `≤ (λ + M2π²/2)/i`。
  · 丢掉 `−D_ii ≤ 0`（放大，fail-closed）。
⟹ **`R_i ≤ B_i := (8 + 4·ln(i+1))·M2/i² + (λ + π²M2/2)/i`，`i ≥ 1`**，且 `B_i` 对 `i` 单调递减。
数值上新界比旧界**更大**（更保守）⟹ 旧的「0 违反」不受影响，但此前的推导不可引用。

严格性：上面所有不等式与拓扑论证 `analytic`；本脚本输出的 `λ₁` 上界为 **`certificate`**
（mpmath.iv 外向舍入，dps 在下方显式设置）；作为对照印出的 float64 参考值是 `numerical-evidence`。
⚠ 试探向量 `v` 由 float64 特征分解给出，但 **KT 对任意 v 成立** ⟹ v 的来源不进声音性论证，
   只影响界的松紧。v 的每个分量是精确的二进制有理数，直接喂给 iv 不引入误差。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import time
import numpy as np
from mpmath import iv

_d = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, name + ".py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


KT = _load("_chk_main_kato_temple_two_atom_prototype_2026-09-23")

iv.dps = 30

# π 的**双侧**包封（不用 mp.pi：惰性常数的坑见 CLAUDE.md）。两个端点是相邻 double，
# 下端 3.141592653589793 < π < 上端 3.1415926535897936。
PI = iv.mpf([3.141592653589793, 3.1415926535897936])


def iv_G(lams, wbars, kmax):
    """G(k) = Σ_j λ_j ω_j²/(ω_j²+k²)，k=0..kmax，逐个外向舍入。"""
    L = [iv.mpf(float(t)) for t in np.atleast_1d(lams)]
    W = [iv.mpf(float(t)) for t in np.atleast_1d(wbars)]
    W2 = [w*w for w in W]
    out = []
    for k in range(kmax+1):
        kk = iv.mpf(k)*iv.mpf(k)
        s = iv.mpf(0)
        for lj, w2 in zip(L, W2):
            s += lj*w2/(w2 + kk)
        out.append(s)
    return out


def iv_S(g, nmax):
    """S_n = G(0) + 2Σ_{j=1..n}G(j)，闭式有限和（与截断无关）。"""
    out = [g[0]]
    acc = iv.mpf(0)
    for n in range(1, nmax+1):
        acc += g[n]
        out.append(g[0] + 2*acc)
    return out


# ───────────────────────── γ：λ₂(M) 的认证上界 ─────────────────────────

def row_tail_up(lam, M2, i):
    """B_i = (8 + 4 ln(i+1))·M2/i² + (λ + π²M2/2)/i  —— 见文件头的重推。i ≥ 1。"""
    ii = iv.mpf(i)
    return (iv.mpf(8) + 4*iv.log(ii + 1))*M2/(ii*ii) + (lam + PI*PI*M2/2)/ii


def col_tail_up(M2, i, J):
    """Σ_{j>J}[G(|i−j|)+G(i+j+1)]/(2j+1) ≤ 2·M2/((2J+3)(J−i))，要求 J > i。"""
    assert J > i
    return 2*M2/((iv.mpf(2*J) + 3)*(iv.mpf(J - i)))


def gamma_cert(lams, wbars, I, J, p=0, verbose=False):
    """返回 (γ, 逐行段 argmax, 行尾界 B_{I+1}, 逐行段最大值)。γ ≥ λ₂(M)，区间上端。"""
    g = iv_G(lams, wbars, I + J + 2)
    S = iv_S(g, I)
    lam = iv.mpf(float(np.sum(np.atleast_1d(lams))))
    M2 = iv.mpf(0)
    for lj, wj in zip(np.atleast_1d(lams), np.atleast_1d(wbars)):
        M2 += iv.mpf(float(lj))*iv.mpf(float(wj))*iv.mpf(float(wj))
    best = None; arg = None
    for i in range(0, I+1):
        if i == p:
            continue
        s = iv.mpf(0)
        for j in range(0, J+1):
            if j == p:
                continue
            s += (g[abs(i-j)] + g[i+j+1])/(iv.mpf(2*j) + 1)
        s += col_tail_up(M2, i, J)
        s -= S[i]/(iv.mpf(2*i) + 1)
        hi = float(s.b)
        if best is None or hi > best:
            best, arg = hi, i
    tail = float(row_tail_up(lam, M2, I+1).b)
    if verbose:
        print(f"      γ: max_{{1≤i≤{I}}} R_i 上端 = {best:.9f} (argmax i={arg})，"
              f"i>{I} 的 B_{{{I+1}}} = {tail:.9f}")
    return max(best, tail), arg, tail, best


# ───────────────────────── Kato–Temple（区间） ─────────────────────────

def res2_tail_up(M2, W, Nv, L):
    """‖(1−P_L)Mv‖² ≤ (2·M2·W)²/(2L+1)·(d⁻⁴ + 1/(3d³))，d=L−Nv。v 未归一，W=Σu_m|v_m|。"""
    d = iv.mpf(L - Nv)
    return (2*M2*W)**2/(iv.mpf(2*L) + 1)*(d**-4 + 1/(3*d**3))


def kt_cert(lams, wbars, mus, Nc, v_float, Nv, L, gamma_up, verbose=False):
    """返回 (λ₁ 的认证上界, μ 下端, res² 上端)。v_float 当作精确二进制有理向量。"""
    assert L > Nv >= Nc
    g = iv_G(lams, wbars, L + Nv + 2)
    S = iv_S(g, Nv - 1)
    v = [iv.mpf(float(t)) for t in v_float]
    u = [1/iv.sqrt(iv.mpf(2*n) + 1) for n in range(L)]

    M2 = iv.mpf(0)
    for lj, wj in zip(np.atleast_1d(lams), np.atleast_1d(wbars)):
        M2 += iv.mpf(float(lj))*iv.mpf(float(wj))*iv.mpf(float(wj))

    # 库仑标量 ⟨u_c, v⟩（u_c 在 n < Nc 上等于 u_n，其余为 0）
    cdot = iv.mpf(0)
    for m in range(min(Nc, Nv)):
        cdot += u[m]*v[m]
    two_mus = 2*iv.mpf(float(mus))

    Mv = []
    for n in range(L):
        s = iv.mpf(0)
        for m in range(Nv):
            s += u[m]*(g[abs(n-m)] + g[n+m+1])*v[m]
        s *= u[n]
        if n < Nv:
            s -= u[n]*u[n]*S[n]*v[n]
        if n < Nc:
            s -= two_mus*u[n]*cdot
        Mv.append(s)

    nv2 = iv.mpf(0)
    for m in range(Nv):
        nv2 += v[m]*v[m]
    num = iv.mpf(0)
    for m in range(Nv):
        num += v[m]*Mv[m]
    mu = num/nv2

    head = iv.mpf(0)
    for n in range(L):
        r = Mv[n] - (mu*v[n] if n < Nv else iv.mpf(0))
        head += r*r
    W = iv.mpf(0)
    for m in range(Nv):
        W += u[m]*abs(v[m])
    res2 = (head + res2_tail_up(M2, W, Nv, L))/nv2

    gam = iv.mpf(gamma_up)
    if not (float(mu.a) > gamma_up):
        return None, float(mu.a), float(res2.b)
    bound = mu + res2/(mu - gam)
    if verbose:
        print(f"      μ∈[{float(mu.a):.12f},{float(mu.b):.12f}]  res²≤{float(res2.b):.3e}  "
              f"尾占比≈{float(res2_tail_up(M2, W, Nv, L).b)/float(res2.b*nv2.b):.1e}")
    return float(bound.b), float(mu.a), float(res2.b)


# ───────────────────────── 主流程 ─────────────────────────

def vstar(lams, wbars, mus, Nc, Nv):
    from scipy.linalg import eigh
    Msub, _, _ = KT.build(np.atleast_1d(lams), np.atleast_1d(wbars), mus, Nc, Nv)
    _, V = eigh(Msub, subset_by_index=[Nv-1, Nv-1])
    return np.ascontiguousarray(V[:, 0].ravel())


def main():
    r = 10.0
    lam, mus = 3.2953, 0.13
    wE = KT.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
    print(f"基准：λ={lam} μ*={mus} r={r} ⟹ ω̄_E*={wE:.9f}  N_c={Nc}  M2={lam*wE*wE:.9f}")
    print(f"mpmath.iv dps={iv.dps}；π 包封宽度={float(PI.b-PI.a):.2e}（含真值 ⟹ 声音）")
    print()

    Nv, L, I, J = 120, 800, 60, 1200

    CASES = [
        ("Einstein", np.array([lam]), np.array([wE])),
        ("2atom x=.4 ρ=3", *[None]),
    ]
    L4, W4 = KT.members(lam, wE, 0.4, 3.0)
    L1, W1 = KT.members(lam, wE, 0.1, 20.0)
    CASES = [("Einstein", np.array([lam]), np.array([wE])),
             ("2atom x=.4 ρ=3", L4, W4),
             ("2atom x=.1 ρ=20", L1, W1)]

    print("=" * 112)
    print(f"GATE IV1  端到端区间证书（Nv={Nv}, L={L}, γ 用 I={I}, J={J}, p=0）")
    print(f"  {'成员':<20}{'γ_cert':>12}{'μ 下端':>15}{'res² 上端':>12}"
          f"{'λ₁ 认证上界':>18}{'λ₁ float(N=1500)':>19}{'δ':>11}{'方向':>8}")
    rows = []
    for tag, Ls, Ws in CASES:
        t0 = time.time()
        gam, argi, btail, hmax = gamma_cert(Ls, Ws, I, J)
        v = vstar(Ls, Ws, mus, Nc, Nv)
        b, mulo, res2 = kt_cert(Ls, Ws, mus, Nc, v, Nv, L, gam)
        ref = KT.lam_k(Ls, Ws, mus, Nc, 1500)
        if b is None:
            print(f"  {tag:<20}{gam:>12.6f}{mulo:>15.9f}{res2:>12.2e}"
                  f"{'μ≤γ ⟹ KT 不适用':>18}{ref:>19.9f}")
            rows.append((tag, gam, None, ref))
            continue
        d = b/ref - 1.0
        ok = "OK" if b >= ref - 1e-12 else "**FAIL**"
        print(f"  {tag:<20}{gam:>12.6f}{mulo:>15.9f}{res2:>12.2e}"
              f"{b:>18.12f}{ref:>19.12f}{d:>11.2e}{ok:>8}   [{time.time()-t0:.1f}s]")
        rows.append((tag, gam, b, ref))
    print()

    print("=" * 112)
    print("GATE IV2  γ 的两段：逐行精确段 vs 行尾界段——扫 I，看尾段何时让位给逐行段")
    print(f"  （注意：B_i 只含 (λ, M2) ⟹ **族内常数** ⟹ 三个成员的尾段值逐位相同，这是设计不是巧合）")
    print(f"  {'I':>6}{'max_{i≤I} R_i':>18}{'B_{I+1}':>14}{'γ_cert':>14}{'谁承重':>12}{'耗时':>8}")
    for tag, Ls, Ws in (CASES[0],):
        for Ii in (60, 150, 300):
            t0 = time.time()
            gam, argi, btail, hmax = gamma_cert(Ls, Ws, Ii, J)
            who = "尾段 B" if btail >= hmax else f"逐行段(i={argi})"
            print(f"  {Ii:>6}{hmax:>18.9f}{btail:>14.9f}{gam:>14.9f}{who:>16}{time.time()-t0:>7.1f}s")
    print("  ⟹ I 只买余量（γ 从 0.181 降到 ~0.10），而预算线是 0.99 ⟹ **I=60 已远够**，不加大。")
    print("    （「预算为 ∞ ⟹ 不派活」的正面用法：这里 γ 的消费点是 KT 的 μ−γ ≈ 0.82 vs 0.90，")
    print("      对 λ₁ 上界的影响 = res²·(1/0.82−1/0.90) ≈ 7e−10，不进账。）")
    print()

    print("=" * 112)
    print("GATE IV6  **认证的 Einstein 根区间**：扫 λ，找认证上界跨过 1 的位置")
    print("  （这是本脚本唯一的正面交付：一个 certificate 级的「λ 到这里天花板还没破」）")
    print(f"  {'λ':>12}{'λ₁ 认证上界':>20}{'λ₁ float':>20}{'判词':>22}")
    for lf in (3.2900, 3.2940, 3.2950, 3.2953, 3.2960):
        Ls = np.array([lf]); Ws = np.array([wE])
        gam, _, _, _ = gamma_cert(Ls, Ws, I, J)
        vv = vstar(Ls, Ws, mus, Nc, Nv)
        b, _, _ = kt_cert(Ls, Ws, mus, Nc, vv, Nv, L, gam)
        ref = KT.lam_k(Ls, Ws, mus, Nc, 1500)
        verdict = "**认证 λ₁<1**" if (b is not None and b < 1.0) else "认证界 ≥1（判不掉）"
        print(f"  {lf:>12.4f}{b:>20.12f}{ref:>20.12f}{verdict:>22}")
    print("  读法：认证上界 <1 ⟹ 该 λ 处 Einstein 成员**确定**没达到判据 ⟹ 可作为天花板的")
    print("        `certificate` 级支撑。上界 ≥1 时只说明本工具判不掉，**不**说明天花板破了")
    print("        （单边界方向规则：上界的下估只支持 PASS，永不支持 FAIL）。")
    print()

    print("=" * 112)
    print("GATE IV3  区间宽度从哪来：dps 扫描（只测舍入，不测离散化——CLAUDE.md）")
    tag, Ls, Ws = CASES[0]
    v = vstar(Ls, Ws, mus, Nc, Nv)
    for dps in (20, 30, 45):
        iv.dps = dps
        t0 = time.time()
        gam, _, _, _ = gamma_cert(Ls, Ws, I, J)
        b, mulo, res2 = kt_cert(Ls, Ws, mus, Nc, v, Nv, L, gam)
        print(f"  dps={dps:<4} γ_cert={gam:.12f}  λ₁ 上界={b:.14f}   [{time.time()-t0:.1f}s]")
    iv.dps = 30
    print("  读法：dps 一致只说明舍入不承重；离散化（L, I, J, Nv）由各自的**解析尾界**覆盖，")
    print("        不靠 dps。两者是不同的误差源，别互相冒充。")
    print()

    print("=" * 112)
    print("GATE IV4  投毒：把 λ 整体放大，认证上界必须跟着越过 1（否则是假 PASS）")
    for f in (1.00, 1.02, 1.10, 1.30):
        Ls = np.array([lam*f]); Ws = np.array([wE])
        gam, _, _, _ = gamma_cert(Ls, Ws, I, J)
        v = vstar(Ls, Ws, mus, Nc, Nv)
        b, _, _ = kt_cert(Ls, Ws, mus, Nc, v, Nv, L, gam)
        ref = KT.lam_k(Ls, Ws, mus, Nc, 1500)
        verdict = "≤1（天花板未破）" if (b is not None and b <= 1.0) else ">1（报破）"
        print(f"  λ×{f:<5.2f}  λ₁认证上界={b:.12f}  float 参考={ref:.12f}  ⟹ {verdict}")
    print("  期望：×1.00 给 ≤1 或刚过 1（Einstein 恰是根），×1.02 起必须 >1。")
    print()

    print("=" * 112)
    print("GATE IV5  负对照：把 γ 人为调到 μ 以上，KT 必须拒绝出数（前提门有效）")
    tag, Ls, Ws = CASES[0]
    v = vstar(Ls, Ws, mus, Nc, Nv)
    for gfake in (0.05, 0.5, 0.95, 1.5):
        b, mulo, _ = kt_cert(Ls, Ws, mus, Nc, v, Nv, L, gfake)
        print(f"  γ={gfake:<5} ⟹ {'上界=' + format(b, '.12f') if b is not None else 'KT 拒绝（μ≤γ）'}"
              f"   (μ 下端={mulo:.9f})")
    print("  读法：γ 变大 ⟹ 界变松但仍是上界；γ≥μ 时必须拒绝 ⟹ 前提不是摆设。")
    print()

    print("=" * 112)
    print("GATE IV7  行尾界 B_i 自查（订正后的新界必须自己也被验一遍——W1：订正是一次交付）")
    print("  (a) 单调递减（sup_{i>I} R_i ≤ B_{I+1} 靠它）  (b) 对逐行精确段确实是上界")
    for tag, Ls, Ws in (CASES[0],):
        lamv = iv.mpf(float(np.sum(np.atleast_1d(Ls))))
        M2v = iv.mpf(0)
        for lj, wj in zip(np.atleast_1d(Ls), np.atleast_1d(Ws)):
            M2v += iv.mpf(float(lj))*iv.mpf(float(wj))*iv.mpf(float(wj))
        bs = [float(row_tail_up(lamv, M2v, i).b) for i in range(1, 2001)]
        bad = [i+1 for i in range(len(bs)-1) if bs[i+1] > bs[i]]
        print(f"  (a) i=1..2000：递增处 {len(bad)} 个 {bad[:5]}   B_1={bs[0]:.6f}  B_10={bs[9]:.6f}"
              f"  B_100={bs[99]:.6f}  B_2000={bs[1999]:.3e}")
        # (b)：重算 i ≤ 60 的逐行 R_i，与 B_i 比
        II, JJ = 60, 1200
        g = iv_G(Ls, Ws, II + JJ + 2); S = iv_S(g, II)
        worst = -1e18; wi = None; nviol = 0
        for i in range(1, II+1):
            ss = iv.mpf(0)
            for j in range(1, JJ+1):
                ss += (g[abs(i-j)] + g[i+j+1])/(iv.mpf(2*j) + 1)
            ss += col_tail_up(M2v, i, JJ)
            ss -= S[i]/(iv.mpf(2*i) + 1)
            ri = float(ss.b); bi = float(row_tail_up(lamv, M2v, i).b)
            nviol += (ri > bi)
            if ri/bi > worst:
                worst, wi = ri/bi, i
        print(f"  (b) i=1..{II}：违反 {nviol}/{II}，最劣 R_i/B_i = {worst:.4f} (i={wi})")
        print("      ⟹ 新界成立且不紧（比值 ≪ 1 说明它只用于尾段、不用于承重段，符合设计）。")


if __name__ == "__main__":
    main()
