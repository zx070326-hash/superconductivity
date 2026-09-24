# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）任务 ②(c)：**盒上区间证书 + 自适应分支定界**（mpmath.iv，外向舍入）。

上一刀（`_chk_main_interval_certificate_single_member_2026-09-23.py`）对**单个点参数成员**
出了 certificate 级的 λ₁ 上界。本脚本把参数本身换成**区间**，对 (q,σ) 正方形做自适应细分，
目标是一句**对整族一致**的判词：

    sup_{ν ∈ 𝒞̄(λ, M2)} λ₁(M[ν]) < 1      ⟹  该 (λ, M2) 下天花板未破（certificate）

坐标与恒等式见 `_chk_main_qsigma_coords_and_deficit_identity_2026-09-23.py`（GATE P0–P3）：
  q := λ₂/λ、σ := λ₂ω₂²/M2，(q,σ) ∈ [0,1]²，**整块正方形就是紧化后的全部极点**
  （对角 σ=q = Einstein；q=0 = 曲线(A)；σ=1 = 曲线(B)；(0,1) = ν≡0）。
  G(k) = λω̄₂²[F(1−σ,1−q) + F(σ,q)]，F(a,b)=ab/(ω̄₂²a+k²b)，F 对 a、b 各自不减；
  G_E(k) − G(k) = λω̄₂⁴k²(σ−q)²/[(ω̄₂²+k²)·u·(W−u)]，u=ω̄₂²σ+k²q，W−u=ω̄₂²(1−σ)+k²(1−q)。

━━━━━━━━ 盒上的 G 包围：两种写法**取交** ━━━━━━━━
  · F-写法：每个变量在每项里只出现一次 ⟹ 该项的区间包围**精确**；在对角线附近是一阶宽。
  · Δ-写法：G ∈ [G_E − Δ_hi, G_E − Δ_lo]；跨对角线的盒子 Δ_lo = 0 ⟹ 上端**恰为 G_E**，
    且 Δ 对「离对角线距离」是二次的 ⟹ 宽度二阶。角点 (0,0)/(1,1) 处 u 或 W−u 可为 0
    ⟹ Δ_hi = ∞，那里由 F-写法接管（一阶宽，与该处真实的一阶变化相符）。
  两者都是声音上/下界 ⟹ **取交仍声音**，且各自补上对方的短板。

━━━━━━━━ 每个盒子的判决链（全部朝「判据变难」放大） ━━━━━━━━
  工具 1（无前提，粗）：λ₁(M) ≤ λ₁(𝔎_ph) ≤ R₀ := max_i[(1/u_i)(Bu)_i − D_ii]（加权 Schur，w=u）
  工具 2（有前提，精）：Kato–Temple  λ₁ ≤ μ_hi + res²_hi/(μ_lo − γ)，要求 μ_lo > γ
      γ ≥ λ₂(M)：库仑 PSD ⟹ M ⪯ 𝔎_ph ⟹ λ₂(M) ≤ λ₂(𝔎_ph)；Cauchy 交错去掉指标 p=0；
                 再对去指标算子用同一个加权 Schur。
      行尾 i>I 用 `B_i = (8+4ln(i+1))M2/i² + (λ+π²M2/2)/i`（§11.1 的重推，**族内常数**，
      因 𝒞̄ 的两条约束是 ≤ ⟹ 任何成员的 Σλ_j ≤ λ、Σλ_jω_j² ≤ M2，而 B_i 对两者都单调不减）。
      列尾 j>J 用 2M2/((2J+3)(J−i))；res² 的 n≥L 尾用 (2M2W)²/(2L+1)(d⁻⁴+1/(3d³))。
  KT 的方向：对盒内任一 ν，μ(ν) ≤ μ_hi 且 res²(ν)/(μ(ν)−γ) ≤ res²_hi/(μ_lo−γ) ⟹ 和 ≤ 上式。

⚠ 试探向量 v 取盒中心成员的 float64 特征向量，**不进声音性论证**（KT 对任意 v 成立），
   只影响松紧；其分量是精确二进制有理数，喂进 iv 不引入误差。

严格性：判词为 `certificate`（mpmath.iv 外向舍入 + 上列解析尾界）；作为对照印出的
float64 λ₁ 为 `numerical-evidence`。盒子数是**实测值**，不是估计（§10.7 的 UNCOMPUTED 由此清掉）。
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
QS = _load("_chk_main_qsigma_coords_and_deficit_identity_2026-09-23")

iv.dps = 25
PI = iv.mpf([3.141592653589793, 3.1415926535897936])      # 双侧包封，含真值


# ───────────────────────── 盒上的 G 包围 ─────────────────────────

def _F(a, b, w2, k2):
    """F(a,b) = ab/(ω̄₂²a + k²b)，a,b ≥ 0 的 float；a=0 或 b=0 时连续延拓为 0。返回 iv。"""
    if a <= 0.0 or b <= 0.0:
        return iv.mpf(0)
    A, B = iv.mpf(a), iv.mpf(b)
    return A*B/(w2*A + k2*B)


def G_box(lam, wE, box, kmax):
    """返回 [None] + [G(k) 的 iv 包围]_{k=1..kmax}，对盒内**每一个** (q,σ) 都有效。"""
    qlo, qhi, slo, shi = box
    lamv, w2 = iv.mpf(lam), iv.mpf(wE)*iv.mpf(wE)
    pref = lamv*w2
    dlo_, dhi_ = slo - qhi, shi - qlo                      # d = σ−q 的范围
    if dlo_ <= 0.0 <= dhi_:
        d2lo, d2hi = 0.0, max(dlo_*dlo_, dhi_*dhi_)
    else:
        d2lo = min(abs(dlo_), abs(dhi_))**2; d2hi = max(dlo_*dlo_, dhi_*dhi_)
    out = [None]
    for k in range(1, kmax+1):
        k2 = iv.mpf(k)*iv.mpf(k)
        # ── F-写法（逐项精确）
        hiF = pref*(_F(1.0-slo, 1.0-qlo, w2, k2) + _F(shi, qhi, w2, k2))
        loF = pref*(_F(1.0-shi, 1.0-qhi, w2, k2) + _F(slo, qlo, w2, k2))
        hi, lo = hiF.b, loF.a
        # ── Δ-写法（近对角二阶）
        GE = pref/(w2 + k2)
        ulo = w2*iv.mpf(slo) + k2*iv.mpf(qlo)
        uhi = w2*iv.mpf(shi) + k2*iv.mpf(qhi)
        vlo = w2*(1 - iv.mpf(shi)) + k2*(1 - iv.mpf(qhi))
        vhi = w2*(1 - iv.mpf(slo)) + k2*(1 - iv.mpf(qlo))
        num = pref*w2*k2/(w2 + k2)
        if float(uhi.a) > 0 and float(vhi.a) > 0:          # Δ 的**下**界要用 u、W−u 的上端
            dlo_iv = num*iv.mpf(d2lo)/(uhi*vhi)
            hi = min(hi, (GE - dlo_iv).b)
        if float(ulo.a) > 0 and float(vlo.a) > 0:          # Δ 的**上**界要用下端
            dhi_iv = num*iv.mpf(d2hi)/(ulo*vlo)
            lo = max(lo, (GE - dhi_iv).a)
        if lo < 0:
            lo = iv.mpf(0).a                                # G ≥ 0（正测度）
        out.append(iv.mpf([lo, hi]))
    return out


# ───────────────────────── 尾界（解析，族内常数） ─────────────────────────

def row_tail_up(lam, M2, i):
    """B_i = (8+4ln(i+1))M2/i² + (λ+π²M2/2)/i（§11.1 重推）。对 λ、M2 单调不减 ⟹ 用族常数合法。"""
    ii = iv.mpf(i)
    return (iv.mpf(8) + 4*iv.log(ii+1))*M2/(ii*ii) + (lam + PI*PI*M2/2)/ii


def col_tail_up(M2, i, J):
    assert J > i
    return 2*M2/((iv.mpf(2*J) + 3)*iv.mpf(J - i))


def res2_tail_up(M2, W, Nv, L):
    d = iv.mpf(L - Nv)
    return (2*M2*W)**2/(iv.mpf(2*L) + 1)*(d**-4 + 1/(3*d**3))


# ───────────────────────── 加权 Schur：R₀（不剥）与 γ（剥 p=0） ─────────────────────────

def schur_bounds(g, lam_iv, M2_iv, I, J, p=0):
    """一次双循环同时给 (R₀, γ)。R₀ ≥ λ₁(𝔎_ph) ≥ λ₁(M)；γ ≥ λ₁(去掉 p 的 𝔎_ph) ≥ λ₂(M)。
    对角用**抵消后**的写法（不含 G(0)）：[G(2i+1) − 2Σ_{j=1..i}G(j)]/(2i+1)。"""
    csum = [iv.mpf(0)]                                     # Σ_{j=1..i} G(j)
    for i in range(1, I+1):
        csum.append(csum[-1] + g[i])
    best0 = best1 = None
    for i in range(0, I+1):
        s_all = iv.mpf(0); s_peel = iv.mpf(0)
        for j in range(0, J+1):
            if j == i:
                continue
            t = (g[abs(i-j)] + g[i+j+1])/(iv.mpf(2*j) + 1)
            s_all += t
            if j != p:
                s_peel += t
        diag = (g[2*i+1] - 2*csum[i])/(iv.mpf(2*i) + 1)    # = [(Bu)_i/u_i 的 j=i 项] − D_ii
        ct = col_tail_up(M2_iv, i, J)
        r0 = float((s_all + diag + ct).b)
        best0 = r0 if best0 is None else max(best0, r0)
        if i != p:
            r1 = float((s_peel + diag + ct).b)
            best1 = r1 if best1 is None else max(best1, r1)
    tail = float(row_tail_up(lam_iv, M2_iv, I+1).b)
    return max(best0, tail), max(best1, tail)


# ───────────────────────── Kato–Temple（盒上） ─────────────────────────

def kt_box(g, lam_iv, M2_iv, mus, Nc, v_float, Nv, L, gamma_up):
    """返回 (bound_hi 或 None, μ 下端, μ 上端, res² 上端)。g 为盒上的 G 包围。"""
    v = [iv.mpf(float(t)) for t in v_float]
    u = [1/iv.sqrt(iv.mpf(2*n) + 1) for n in range(L)]
    csum = [iv.mpf(0)]
    for i in range(1, Nv):
        csum.append(csum[-1] + g[i])

    cdot = iv.mpf(0)
    for m in range(min(Nc, Nv)):
        cdot += u[m]*v[m]
    two_mus = 2*iv.mpf(float(mus))

    Mv = []
    for n in range(L):
        s = iv.mpf(0)
        for m in range(Nv):
            if m == n:
                continue
            s += u[m]*(g[abs(n-m)] + g[n+m+1])*v[m]
        s *= u[n]
        if n < Nv:                                          # 对角：G(0) 已抵消
            s += u[n]*u[n]*(g[2*n+1] - 2*csum[n])*v[n]
        if n < Nc:
            s -= two_mus*u[n]*cdot
        Mv.append(s)

    nv2 = iv.mpf(0); num = iv.mpf(0); W = iv.mpf(0)
    for m in range(Nv):
        nv2 += v[m]*v[m]; num += v[m]*Mv[m]; W += u[m]*abs(v[m])
    mu = num/nv2
    head = iv.mpf(0)
    for n in range(L):
        rr = Mv[n] - (mu*v[n] if n < Nv else iv.mpf(0))
        head += rr*rr
    res2 = (head + res2_tail_up(M2_iv, W, Nv, L))/nv2
    mlo, mhi, r2 = float(mu.a), float(mu.b), float(res2.b)
    if not (mlo > gamma_up):
        return None, mlo, mhi, r2
    return mhi + r2/(mlo - gamma_up), mlo, mhi, r2


# ───────────────────────── 分支定界 ─────────────────────────

def vstar(lam, wE, q, s, mus, Nc, Nv):
    from scipy.linalg import eigh
    L_, W_ = QS.mem_qs(lam, wE, q, s)
    Msub, _, _ = KT.build(L_, W_, mus, Nc, Nv)
    _, V = eigh(Msub, subset_by_index=[Nv-1, Nv-1])
    return np.ascontiguousarray(V[:, 0].ravel())


def certify(lam, wE, mus, Nc, target, I=60, J=400, Nv=40, L=160,
            maxdepth=40, maxbox=20000, verbose=True, log_every=100):
    """对 (q,σ) 正方形的**上三角** σ≥q 做自适应细分（下三角由对称覆盖）。

    每个盒子：一次加权 Schur 双循环同时给 (R₀, γ) ——
      · R₀ < target ⟹ 收（regime 1，无前提，远盒用）；
      · 否则若 μ_lo > γ ⟹ 走 Kato–Temple，b < target 则收（regime 2，近对角用）；
      · 都不收 ⟹ 沿较长边中点二分。
    γ **逐盒**算，不复用根盒值：根盒 γ=0.608 会把 KT 的适用下限抬到 0.61，而 R₀<1 只覆盖
    λ₁≲0.28 ⟹ 中间出现两工具皆不收的带，细分无法消除。逐盒 γ≈0.18（行尾主导）时两区重叠。

    对称性（analytic）：交换两个原子 ⟺ (q,σ)↦(1−q,1−σ) 给出**同一个测度**（GATE B2s）
    ⟹ 严格下三角 σ<q 的成员与某个上三角成员逐字相同；闭盒且中点二分 ⟹ 跳过 shi≤qlo 合法。
    """
    lam_iv = iv.mpf(lam)
    M2_iv = iv.mpf(lam)*iv.mpf(wE)*iv.mpf(wE)
    kmax = max(I + J + 2, L + Nv + 2)
    t0 = time.time()
    stack = [(0.0, 1.0, 0.0, 1.0, 0)]
    nleaf = n1 = n2 = nsplit = nskip = nproc = 0
    worst = -1e18; worst_box = None; fails = []; leaves = []
    gmax = -1e18; maxdep_seen = 0
    while stack:
        qlo, qhi, slo, shi, dep = stack.pop()
        if shi <= qlo:                                       # 严格下三角（含只触到对角一点）
            nskip += 1; continue
        # 限额必须数**已处理的盒子**：DFS 的栈深 ≤ maxdepth，用 nleaf+len(stack) 计数在
        # 投毒（永远收不了）时**永远触发不了**，会把深度 maxdepth 的整棵树遍历完。
        nproc += 1
        if nproc > maxbox:
            fails.append(("BUDGET", qlo, qhi, slo, shi, dep)); break
        maxdep_seen = max(maxdep_seen, dep)
        g = G_box(lam, wE, (qlo, qhi, slo, shi), kmax)
        r0, gam = schur_bounds(g, lam_iv, M2_iv, I, J)
        gmax = max(gmax, gam)
        acc = None; b = None
        if r0 < target:
            acc, b = 1, r0
        else:
            qc, sc = 0.5*(qlo+qhi), 0.5*(slo+shi)
            v = vstar(lam, wE, qc, sc, mus, Nc, Nv)
            b, mlo, mhi, r2 = kt_box(g, lam_iv, M2_iv, mus, Nc, v, Nv, L, gam)
            if b is not None and b < target:
                acc = 2
        if acc is not None:
            nleaf += 1; n1 += (acc == 1); n2 += (acc == 2)
            leaves.append((qlo, qhi, slo, shi, acc, b))
            if b > worst:
                worst, worst_box = b, (qlo, qhi, slo, shi, acc)
            if verbose and nleaf % log_every == 0:
                print(f"    …{nleaf} 叶 (深度≤{maxdep_seen}, 栈{len(stack)}, {time.time()-t0:.0f}s)",
                      flush=True)
            continue
        if dep >= maxdepth:
            fails.append(("DEPTH", qlo, qhi, slo, shi, dep)); continue
        nsplit += 1
        if (qhi - qlo) >= (shi - slo):
            qm = 0.5*(qlo + qhi)
            stack.append((qlo, qm, slo, shi, dep+1)); stack.append((qm, qhi, slo, shi, dep+1))
        else:
            sm = 0.5*(slo + shi)
            stack.append((qlo, qhi, slo, sm, dep+1)); stack.append((qlo, qhi, sm, shi, dep+1))
    return dict(nleaf=nleaf, n_R0=n1, n_KT=n2, nsplit=nsplit, worst=worst, gamma_max=gmax,
                worst_box=worst_box, fails=fails, maxdepth=maxdep_seen, nskip=nskip,
                nproc=nproc, leaves=leaves, secs=time.time()-t0)


def main():
    r, mus = 10.0, 0.13
    lam_star = 3.2953
    wE_star = KT.einstein_root(lam_star, mus, r, 600)
    M2 = lam_star*wE_star*wE_star                            # 固定 η/M ⟹ M2 固定
    print(f"固定 M2 = λ⟨ω²⟩ = {M2:.9f}（对应 η/M=9.172 在 T=300K 的松原单位）；μ*={mus}, r={r}")
    print(f"Einstein 根：λ*={lam_star}（ω̄₂={wE_star:.9f}），λ<λ* 时 Einstein 成员 λ₁<1")
    print(f"mpmath.iv dps={iv.dps}")
    print()

    def setup(lam):
        wE = float(np.sqrt(M2/lam)); Nc = int(np.floor(r*wE))
        return wE, Nc

    # ───────── GATE B0：盒上 G 包围确实**外包**盒内成员（含角点与随机内点） ─────────
    print("=" * 112)
    print("GATE B0  盒上 G 包围的外包性：盒内随机 (q,σ) 的真 G(k) 必须落在 [G_lo,G_hi] 内")
    lam = 3.29; wE, Nc = setup(lam)
    rng = np.random.default_rng(7)
    boxes = [(0.0, 1.0, 0.0, 1.0), (0.3, 0.35, 0.3, 0.35), (0.0, 0.02, 0.0, 0.02),
             (0.98, 1.0, 0.98, 1.0), (0.1, 0.2, 0.7, 0.9), (0.0, 0.05, 0.9, 1.0)]
    nbad = 0; nchk = 0; relw = []
    for bx in boxes:
        g = G_box(lam, wE, bx, 60)
        wid = max(float(g[k].b - g[k].a)/max(float(g[k].b), 1e-300) for k in range(1, 61))
        relw.append(wid)
        pts = [(bx[0], bx[2]), (bx[0], bx[3]), (bx[1], bx[2]), (bx[1], bx[3])]
        pts += [(rng.uniform(bx[0], bx[1]), rng.uniform(bx[2], bx[3])) for _ in range(40)]
        for (q, s) in pts:
            L_, W_ = QS.mem_qs(lam, wE, q, s)
            gg = QS.G_of(L_, W_, np.arange(1, 61))
            for k in range(1, 61):
                nchk += 1
                if not (float(g[k].a) - 1e-14 <= gg[k-1] <= float(g[k].b) + 1e-14):
                    nbad += 1
        print(f"  盒 q∈[{bx[0]:.2f},{bx[1]:.2f}] σ∈[{bx[2]:.2f},{bx[3]:.2f}]："
              f"最大相对宽度 {wid:.3e}")
    print(f"  ⟹ 外包违反 {nbad}/{nchk}")
    print()

    # ───────── GATE B1：单盒方向门（盒上界 ≥ 盒内 float64 λ₁） ─────────
    print("=" * 112)
    print("GATE B1  单盒端到端：盒上认证上界 vs 盒内成员的 float64 λ₁（方向必须是 ≥）")
    lam_iv = iv.mpf(lam); M2_iv = iv.mpf(lam)*iv.mpf(wE)*iv.mpf(wE)
    I, J, Nv, L = 25, 250, 50, 250
    kmax = max(I + J + 2, L + Nv + 2)
    print(f"  参数 I={I} J={J} Nv={Nv} L={L}；λ={lam} ω̄₂={wE:.9f} N_c={Nc}")
    print(f"  {'盒':<34}{'R₀':>10}{'γ':>10}{'KT 上界':>16}{'盒内 max λ₁(float)':>20}{'方向':>7}{'秒':>7}")
    for bx in [(0.3, 0.35, 0.3, 0.35), (0.3, 0.4, 0.35, 0.45), (0.0, 0.1, 0.0, 0.1),
               (0.1, 0.2, 0.7, 0.9), (0.45, 0.55, 0.45, 0.55)]:
        t0 = time.time()
        g = G_box(lam, wE, bx, kmax)
        r0, gam = schur_bounds(g, lam_iv, M2_iv, I, J)
        qc, sc = 0.5*(bx[0]+bx[1]), 0.5*(bx[2]+bx[3])
        v = vstar(lam, wE, qc, sc, mus, Nc, Nv)
        b, mlo, mhi, r2 = kt_box(g, lam_iv, M2_iv, mus, Nc, v, Nv, L, gam)
        mx = -1e18
        for _ in range(60):
            q, s = rng.uniform(bx[0], bx[1]), rng.uniform(bx[2], bx[3])
            L_, W_ = QS.mem_qs(lam, wE, q, s)
            mx = max(mx, KT.lam_k(L_, W_, mus, Nc, 600))
        for (q, s) in [(bx[0], bx[2]), (bx[1], bx[3])]:
            L_, W_ = QS.mem_qs(lam, wE, q, s)
            mx = max(mx, KT.lam_k(L_, W_, mus, Nc, 600))
        ok = "OK" if (b is None or b >= mx - 1e-12) and r0 >= mx - 1e-12 else "**FAIL**"
        bs = "μ≤γ" if b is None else f"{b:.9f}"
        print(f"  q[{bx[0]:.2f},{bx[1]:.2f}] σ[{bx[2]:.2f},{bx[3]:.2f}]          "
              f"{r0:>10.4f}{gam:>10.4f}{bs:>16}{mx:>20.9f}{ok:>7}{time.time()-t0:>7.1f}")
    print()
    # ───────── GATE B2s：对称性 (q,σ)↦(1−q,1−σ)（上三角裁剪的**前提**） ─────────
    print("=" * 112)
    print("GATE B2s  点反射 (q,σ)↦(1−q,1−σ) 必须给出**同一个测度** ⟹ 只认证上三角合法")
    nbad = 0; worst_d = 0.0; worst_l = 0.0
    for _ in range(400):
        q, s = rng.uniform(0, 1), rng.uniform(0, 1)
        L1, W1 = QS.mem_qs(lam, wE, q, s)
        L2, W2 = QS.mem_qs(lam, wE, 1.0-q, 1.0-s)
        a = sorted(zip(np.atleast_1d(L1), np.atleast_1d(W1)))
        b = sorted(zip(np.atleast_1d(L2), np.atleast_1d(W2)))
        d = max(abs(x[0]-y[0]) + abs(x[1]-y[1]) for x, y in zip(a, b)) if len(a) == len(b) else 9e9
        dl = abs(KT.lam_k(L1, W1, mus, Nc, 300) - KT.lam_k(L2, W2, mus, Nc, 300))
        worst_d = max(worst_d, d); worst_l = max(worst_l, dl)
        if d > 1e-12 or dl > 1e-12:
            nbad += 1
    print(f"  400 组随机 (q,σ)：原子表最劣偏差 {worst_d:.3e}，λ₁ 最劣偏差 {worst_l:.3e}"
          f"  ⟹ 违反 {nbad}/400")
    print()

    # ───────── GATE B2：整族自适应分支定界（**盒子数是这里测出来的**） ─────────
    print("=" * 112)
    print("GATE B2  整族分支定界：sup_{ν∈𝒞̄} λ₁(M[ν]) < 1 ？（certificate；盒子数为实测）")
    I2, J2, Nv2, L2_ = 60, 400, 40, 160
    print(f"  工具参数 I={I2} J={J2} Nv={Nv2} L={L2_}；判据 target=1.0（严格 <）")
    rows = []
    for lam in [3.20, 3.26, 3.29]:
        wE, Nc = setup(lam)
        lE = KT.lam_k([lam], [wE], mus, Nc, 600)
        print(f"  ── λ={lam}  ω̄₂={wE:.9f}  N_c={Nc}  λ₁(Einstein,float)={lE:.9f}"
              f"  余量 {1-lE:.3e}", flush=True)
        R = certify(lam, wE, mus, Nc, 1.0, I=I2, J=J2, Nv=Nv2, L=L2_, maxbox=40000)
        rows.append((lam, lE, R))
        print(f"     盒子 {R['nleaf']}（R₀ 收 {R['n_R0']}，KT 收 {R['n_KT']}）"
              f"，分裂 {R['nsplit']}，对称跳过 {R['nskip']}，最大深度 {R['maxdepth']}"
              f"，最劣上界 {R['worst']:.9f}，γ_max {R['gamma_max']:.4f}"
              f"，{R['secs']:.0f}s，失败盒 {len(R['fails'])}", flush=True)
        for f in R["fails"][:5]:
            print(f"       FAIL {f}")
    print()
    print(f"  {'λ':>7}{'λ₁(E,float)':>15}{'盒子数':>9}{'R₀/KT':>12}{'最劣认证上界':>18}"
          f"{'深度':>6}{'秒':>8}{'判词':>10}")
    for lam, lE, R in rows:
        ok = "PASS" if (not R["fails"] and R["worst"] < 1.0) else "FAIL"
        split = f"{R['n_R0']}/{R['n_KT']}"
        print(f"  {lam:>7.2f}{lE:>15.9f}{R['nleaf']:>9}{split:>12}"
              f"{R['worst']:>18.9f}{R['maxdepth']:>6}{R['secs']:>8.0f}{ok:>10}")
    print()


if __name__ == "__main__":
    main()
