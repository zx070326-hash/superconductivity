# -*- coding: utf-8 -*-
"""
主会话从零：细化 (a) 的**加强版**——把极值性从「抽样未证伪」换成**逐 k 的单调性导数**。

━━ 结构事实（本轮从 AD1975 Appendix A 末段读出，analytic 部分）━━
1. 固定 (λ, ω̄₂) 的族里，`G[0] = 2∫α²F/ω = λ` 是**常数**；只有 k≥1 的 G[k] 变。
2. 两矩都对 α²F **线性**，G 也对 α²F 线性 ⟹ 该族（及其 G 像）是**凸集**，Einstein 是极点。
3. 库仑项 `−2μ* u_c u_cᵀ` 与 α²F **完全无关** ⟹ **μ\* 只通过顶本征向量 v 进入极值性问题**。
4. AD1975 (A10)+(A13) 末段给的是 **E1**：`x/(x+a²)` 二阶导恒负 ⟹ Einstein 谱**逐 k 最大化 G[k]**。
   （这一步是 AD 的，analytic，μ\* 无关——本脚本数值复核。）
⟹ **E1 + 「λ₁ 对 {G[k]}_{k≥1} 单调不减」 ⟹ Einstein 极值性（对任意 μ\*）。**
   仓库既有判断「新颖性只在 E1⇒本征值」由此得到显式形式：缺的**正是**这条单调性。

━━ 导数（解析，本脚本对有限差分验证）━━
  M_ij = (G[|i−j|]+G[i+j+1])u_i u_j − δ_ij(L1_i + G[0]u_i²) − 2μ*(u_c)_i(u_c)_j
  L1_n = (Σ_{k=1}^{n} 2G[k])/(2n+1) = 2u_n² Σ_{k=1}^{n} G[k]
  令 w := v ⊙ u（v 为归一顶本征向量）。对 k ≥ 1：
      D[k] := ∂λ₁/∂G[k] = 2Σ_{n≥k} w_n w_{n−k}  +  Σ_{i+j=k−1} w_i w_j  −  2Σ_{n≥k} w_n²
  （Hellmann–Feynman，λ₁ 单重时成立；本脚本印谱隙确认单重。）

━━ 逐 k 单调性是**死路**，且**由库仑驱动**（本轮门 1 / 1b 实测）━━
  实测 D[10] = −4.97e−2 < 0（λ=2, μ*=0.13, r=10）⟹ **λ₁ 对 {G[k]} 不单调**。
  门 1b：同一 (λ,谱,r) 下 μ*=0 时 D[k]<0 的 k **一个都没有**（min D = +0.0），
  μ*=0.13 时有 18 个（k=6–311）⟹ **失效是库仑引起的，不是 Z 通道的固有性质**。
  （2026-09-23 订正：此前我写「与 μ* 无关」，被门 1b 直接推翻。）
  机制假说（门 1c 检验，不承重）：μ*>0 使能隙函数出现**节点** ⟹ v 带号 ⟹ w=v⊙u 带号
  ⟹ 交叉项 Σ w_n w_{n−k} 可为负。μ*=0 时 M+cI 逐元非负 ⟹ Perron ⟹ v≥0 ⟹ 交叉项非负。
  ⟹ 「E1 + 逐 k 单调 ⟹ 极值性」这条路在 μ*>0 处**不成立**——而 μ*>0 正是待证的那一块。
  ⟹ 正确的量是沿线段的**方向导数**，见下。

━━ 方向导数（本脚本的真正判词）━━
  线段 G_t := (1−t)G_spec + t·G_E 全程留在族内（凸性，第 2 条）⟹ 合法路径。
      Ψ(t) := dλ₁/dt = Σ_{k≥1} D[k]·(G_E[k] − G_spec[k])
  E1 保证括号**逐项 ≥ 0**；D[k] 变号 ⟹ 和的符号才是判词。
  **Ψ(t) ≥ 0 沿整条线段 ⟹ λ₁(Einstein) ≥ λ₁(spec)**，即该谱不优于 Einstein。

━━ 方向审计（承重）━━
· **单侧**：Ψ≥0 只给「不优于」（PASS 方向）。出现 Ψ<0 **不**等于反例——沿线段某点下降
  仍可能两端点顺序不变；须另行直接比 λ₁。⟹ Ψ<0 一律报「未判定」，不报 FAIL。
· 正对照承重：已知失效带（λ ≲ 2μ*）**必须**看到 Ψ<0；看不到 ⟹ 本诊断是盲的。
· 门 1 的有限差分误差是**绝对**的（≈2ε·λ₁/(2h)=2e−10，与 |D| 无关）⟹ 容差必须挂绝对轴：
  `|差| ≤ 1e−5·|D| + 3e−9`。用相对容差卡小 |D| 的点 = 要求差分分辨到噪声地板以下（**不是**解析式错）。

严格性：E1 与凸性 analytic（此处数值复核）；D[k] 的符号扫描 numerical-evidence。
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.linalg import eigh

KB = 8.617333262e-2
rng = np.random.default_rng(20260923)
WB = 150.0                      # ω̄₂，标度无关


def build(T, G, N, mus, wc):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    L1 = np.zeros(N)
    L1[1:] = 2.0 * np.cumsum(G[1:N]) / (2.0 * n[1:] + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    M = (G[np.abs(i - j)] + G[i + j + 1]) * u[i] * u[j] - np.diag(L1)
    M[n, n] -= G[0] * u * u
    if mus > 0.0:
        Nc = int(np.floor(wc / (2.0 * np.pi * T)))
        if Nc > 0:
            uc = np.where(n < Nc, u, 0.0)
            M -= 2.0 * mus * np.outer(uc, uc)
    return M, u


def top2(T, G, N, mus, wc):
    M, u = build(T, G, N, mus, wc)
    ev, V = eigh(M, subset_by_index=[N - 2, N - 1])
    return ev[1], ev[0], V[:, 1], u          # λ₁, λ₂, v₁, u


def Gk(kmax, T, aw):
    x2 = (2.0 * np.pi * T * np.arange(kmax, dtype=float)) ** 2
    G = np.zeros(kmax)
    for a, w in aw:
        G += 2.0 * a * w / (w * w + x2)
    return G


def Tc(aw, N, mus, wc, lo=0.02, hi=1200.0):
    f = lambda T: top2(T, Gk(2 * N + 2, T, aw), N, mus, wc)[0] - 1.0
    lo, hi = lo * KB, hi * KB
    if not (f(lo) > 0 > f(hi)):
        return np.nan
    for _ in range(46):
        m = np.sqrt(lo * hi)
        if f(m) > 0: lo = m
        else: hi = m
    return np.sqrt(lo * hi) / KB


def two_atom(lam, wbar, w1, wB):
    p, q = lam / 2.0, lam * wbar ** 2 / 2.0
    d = wB ** 2 - w1 ** 2
    if d <= 1e-12: return None
    a1 = w1 * (p * wB ** 2 - q) / d
    a2 = wB * (q - p * w1 ** 2) / d
    if a1 < -1e-14 or a2 < -1e-14: return None
    return [(max(a1, 0.0), w1), (max(a2, 0.0), wB)]


def Dk(v, u, N, kmax):
    """D[k] = ∂λ₁/∂G[k]，k = 1..kmax-1。"""
    w = v * u
    D = np.zeros(kmax)
    tail = np.cumsum(w[::-1] ** 2)[::-1]                # tail[n] = Σ_{m≥n} w_m²
    conv = np.convolve(w, w)                            # conv[s] = Σ_{i+j=s} w_i w_j
    for k in range(1, kmax):
        a = 2.0 * np.dot(w[k:], w[:N - k]) if k < N else 0.0
        b = conv[k - 1] if k - 1 < len(conv) else 0.0
        c = 2.0 * (tail[k] if k < N else 0.0)
        D[k] = a + b - c
    return D


# ── 门 1：D[k] 对有限差分（承重，解析式若抄错这里必炸）────────────────
print("=" * 104)
print("门 1：D[k] 解析式 vs 有限差分（Hellmann–Feynman）")
print("=" * 104)
N, MUS, R = 300, 0.13, 10.0
lam = 2.0
aw = two_atom(lam, WB, 0.6 * WB, 2.0 * WB)
T = Tc(aw, N, MUS, R * WB) * KB
G = Gk(2 * N + 2, T, aw)
l1, l2, v, u = top2(T, G, N, MUS, R * WB)
D = Dk(v, u, N, 2 * N + 1)
print(f"  λ₁ = {l1:.10f}   λ₂ = {l2:.10f}   谱隙 = {l1-l2:.3e}  （需远离 0 才能用 H–F）")
# 容差挂在**绝对**轴上：中心差分的舍入噪声 ≈ 2ε·λ₁/(2h) = 2e−16/1e−6 = 2e−10（与 |D| 无关），
# 截断项 O(h²·λ₁''') ≈ 1e−12。取 pad = 3e−9（≈15× 噪声地板）。
# 方向审计：这是**放松**验收门，必须确认它仍能抓错——解析式若抄错，误差是 O(|D|)，
# 在最小的受测点 k=150 上为 5.5e−6，是 pad 的 1800×；⟹ 放松后仍远不可能放行一个错公式。
FD_PAD, FD_REL = 3e-9, 1e-5
worst_abs, worst_ratio, ntested = 0.0, 0.0, 0
for k in [1, 2, 3, 5, 10, 50, 150, 299, 400]:
    h = 1e-6
    Gp = G.copy(); Gp[k] += h
    Gm = G.copy(); Gm[k] -= h
    fd = (top2(T, Gp, N, MUS, R * WB)[0] - top2(T, Gm, N, MUS, R * WB)[0]) / (2 * h)
    ae = abs(fd - D[k])
    tol = FD_REL * abs(D[k]) + FD_PAD
    worst_abs = max(worst_abs, ae); worst_ratio = max(worst_ratio, ae / tol); ntested += 1
    print(f"    k={k:4d}   解析 {D[k]:+.8e}   差分 {fd:+.8e}   绝对差 {ae:.2e}  /容差 {ae/tol:6.3f}"
          f"{'   ← D<0!' if D[k] < 0 else ''}")
assert ntested == 9 and worst_ratio < 1.0, "D[k] 解析式与差分不符，停"
print(f"  ⟹ {ntested}/9 个点全部通过，最劣 |差| = {worst_abs:.2e}，最劣 差/容差 = {worst_ratio:.3f} < 1\n")

# ── 门 1b：μ*=0 时 D[k] 是否也变号？（决定单调性死路是否与库仑无关）────
print("=" * 104)
print("门 1b：单调性失效与 μ* 有关吗？（若 μ*=0 也变号 ⟹ 与库仑无关，是 Z 通道的固有性质）")
print("=" * 104)
for mus_t in (0.0, 0.05, 0.13, 0.25):
    aw_t = two_atom(2.0, WB, 0.6 * WB, 2.0 * WB)
    Tt = Tc(aw_t, N, mus_t, R * WB) * KB
    Gt = Gk(2 * N + 2, Tt, aw_t)
    _, _, vt, ut = top2(Tt, Gt, N, mus_t, R * WB)
    Dt = Dk(vt, ut, N, 2 * N + 1)
    neg = np.where(Dt[1:] < 0)[0] + 1
    vs = vt * np.sign(vt[np.argmax(np.abs(vt))])       # 定向：主分量取正
    nnode = int(np.sum(vs[:-1] * vs[1:] < 0))          # 门 1c：v 的变号次数（能隙节点）
    print(f"  μ*={mus_t:.2f}  Tc={Tt/KB:7.2f} K   D[k]<0 的 k："
          f"{'无' if len(neg)==0 else f'{neg.min()}–{neg.max()}（共 {len(neg)} 个）':<22}"
          f"   min D = {Dt[1:].min():+.3e}   v 变号次数 = {nnode}   min v = {vs.min():+.3e}")
print("  门 1c：若「v 变号次数 = 0 ⟺ D 无负」逐行成立 ⟹ 支持节点机制假说（不承重，仅定向）。")
print()

# ── 门 2：E1（Einstein 逐 k 最大化 G[k]）数值复核 ─────────────────────
print("=" * 104)
print("门 2：E1 —— Einstein 谱是否逐 k 最大化 G[k]？（AD1975 (A10)(A13)，analytic；此处复核）")
print("=" * 104)
bad = 0; ntest = 0; worstE = 0.0
for lam in (1.5, 3.0, 5.0):
    GE = Gk(2 * N + 2, T, [(lam * WB / 2.0, WB)])
    for w1f in (0.05, 0.3, 0.7, 0.95):
        for w2f in (1.05, 1.5, 3.0, 20.0):
            a = two_atom(lam, WB, w1f * WB, w2f * WB)
            if a is None: continue
            Gs = Gk(2 * N + 2, T, a)
            d = GE[1:] - Gs[1:]
            ntest += 1
            if d.min() < -1e-13:
                bad += 1
            worstE = min(worstE, d.min())
print(f"  {ntest} 个 2 原子谱 × 全部 k≥1：E1 违反 {bad} 个，最劣 min(G_E−G_spec) = {worstE:.3e}")
print(f"  （另核 G[0]：Einstein 与任一同族谱皆 = λ，差 "
      f"{abs(Gk(2*N+2,T,[(3.0*WB/2,WB)])[0] - Gk(2*N+2,T,two_atom(3.0,WB,0.3*WB,3.0*WB))[0]):.2e}）\n")

# ── 门 3：正对照 + 判据区，沿线段扫方向导数 Ψ(t) 的符号 ────────────────
def seg_scan(lam, mus, r, N=N, nt=9):
    """沿 G_t=(1−t)G_spec+t·G_E 扫 Ψ(t)=Σ_k D[k](G_E[k]−G_spec[k])。
       返回 (Tc_E, Nc, minΨ, argmin, 最小谱隙, 点数, E1违反数, 端点直比的最大增益)。"""
    wc = r * WB
    E = [(lam * WB / 2.0, WB)]
    tE = Tc(E, N, mus, wc)
    if not np.isfinite(tE):
        return None
    T = tE * KB
    Nc = int(np.floor(wc / (2.0 * np.pi * T)))
    GE = Gk(2 * N + 2, T, E)
    lE = top2(T, GE, N, mus, wc)[0]
    mn, arg, gapmin, e1bad, gmax = np.inf, None, np.inf, 0, -np.inf
    psi1min, psi1arg, ghatmax, ghatarg = np.inf, None, -np.inf, None
    # 归一：两个端点统计量都除以方向长度 ‖ΔG‖₁，否则极值恒被近简并谱 (0.999,1.001) 取到
    # （那里 ΔG≈0 ⟹ Ψ(1)≈0、Δλ₁≈0，两列一致是恒等式而非独立证据）。
    specs = []
    for w1f in (0.02, 0.2, 0.5, 0.8, 0.95, 0.999):
        for w2f in (1.001, 1.05, 1.3, 2.0, 5.0, 50.0):
            a = two_atom(lam, WB, w1f * WB, w2f * WB)
            if a is not None:
                specs.append(((w1f, w2f), Gk(2 * N + 2, T, a)))
    for tagp, Gs in specs:
        dG = GE[1:] - Gs[1:]                      # E1：应逐项 ≥ 0
        if dG.min() < -1e-13:
            e1bad += 1
        nrm = float(np.abs(dG).sum())
        if nrm < 1e-12:
            continue                              # 完全简并方向：无信息，排除
        dl = top2(T, Gs, N, mus, wc)[0] - lE      # >0 即直接反例
        gmax = max(gmax, dl)
        if dl / nrm > ghatmax:
            ghatmax, ghatarg = dl / nrm, tagp
        for t in np.linspace(0.0, 1.0, nt):
            Gt = (1 - t) * Gs + t * GE
            l1, l2, v, u = top2(T, Gt, N, mus, wc)
            gapmin = min(gapmin, l1 - l2)
            psi = float(np.dot(Dk(v, u, N, 2 * N + 2)[1:], dG))
            if psi < mn:
                mn, arg = psi, (tagp, round(float(t), 2))
            if t == 1.0 and psi / nrm < psi1min:
                psi1min, psi1arg = psi / nrm, tagp
    return tE, Nc, mn, arg, gapmin, len(specs) * nt, e1bad, gmax, psi1min, psi1arg, ghatmax, ghatarg


print("=" * 104)
print("门 3：沿线段扫 min_t Ψ(t)，Ψ = Σ_k D[k]·(G_E[k]−G_spec[k])")
print("      Ψ≥0 全程 ⟹ λ₁(Einstein) ≥ λ₁(spec)（PASS 方向）；Ψ<0 ⟹ 未判定，非 FAIL")
print("=" * 104)
print(f"  {'区段':>8}{'μ*':>6}{'r':>4}{'λ':>7}{'Tc_E[K]':>9}{'N_c':>6}"
      f"{'min_t Ψ':>11}{'min Ψ̂(1)':>12}{'argmin Ψ̂':>16}{'max Δ̂λ₁':>12}{'argmax Δ̂':>16}{'最小谱隙':>10}{'E1违反':>7}")
pos_neg, pos_rows, crit_neg, crit_rows, crit_min, crit_ghat = 0, 0, 0, 0, np.inf, -np.inf
loc_rows, loc_bad, loc_min = 0, 0, np.inf
for tag, mus, r, lams in [("正对照", 0.30, 5.0, [0.54, 0.60, 0.66, 0.84]),
                          ("判据区", 0.13, 10.0, [1.5, 2.0, 3.0, 4.0, 5.0]),
                          ("判据区", 0.13, 6.0, [2.0, 5.0])]:
    for lam in lams:
        out = seg_scan(lam, mus, r)
        if out is None:
            print(f"  {tag:>8}{mus:6.2f}{r:4.0f}{lam:7.3f}   — Einstein 本身不超导"); continue
        tE, Nc, mn, arg, gap, npt, e1bad, gmax, p1, p1a, ghat, gha = out
        flag = " ⚠N_c>N" if Nc > N else ""
        print(f"  {tag:>8}{mus:6.2f}{r:4.0f}{lam:7.3f}{tE:9.2f}{Nc:6d}"
              f"{mn:11.3e}{p1:12.3e}{str(p1a):>16}{ghat:12.3e}{str(gha):>16}{gap:10.2e}{e1bad:7d}"
              f"{'  ←Ψ<0' if mn < 0 else ''}{'  ←Ψ̂(1)<0!' if p1 < 0 else ''}"
              f"{'  ←直接反例!' if gmax > 0 else ''}{flag}", flush=True)
        if Nc > N:
            continue
        loc_rows += 1; loc_bad += (p1 < 0); loc_min = min(loc_min, p1)
        if tag == "正对照":
            pos_rows += 1; pos_neg += (mn < 0)
        else:
            crit_rows += 1; crit_neg += (mn < 0)
            crit_min = min(crit_min, mn); crit_ghat = max(crit_ghat, ghat)

print()
print("-" * 104)
if pos_neg == 0:
    print(f"⚠⚠ VOID：{pos_rows} 个正对照行全部 Ψ≥0。已知失效带（λ≲2μ*）**应当**出现 Ψ<0")
    print("   ⟹ 本诊断对失效是盲的 ⟹ 判据区的 Ψ≥0 不作数。")
else:
    print(f"✓ 正对照亮了：{pos_neg}/{pos_rows} 行出现 Ψ<0 ⟹ 本诊断**不盲**。")
    print(f"⟹ 判据区：{crit_rows} 行，Ψ<0 的 {crit_neg} 行，全局 min_t Ψ = {crit_min:.3e}；"
          f"端点直比（归一）最大 Δ̂λ₁ = {crit_ghat:.3e}")
    if crit_neg == 0:
        print("   判词：**扫到的每条线段上 λ₁ 单调不减 ⟹ 这些端点谱皆不优于 Einstein**。")
        print("   ⚠ 覆盖面：线段族有限，Ψ≥0 是充分条件；**不写「已证极值性」**。")
    else:
        print("   判词：Ψ 这条**充分条件**在判据区失效 ⟹ 极值性在此**未判定**（非 FAIL）。")
        print("   ⟹ 须看「端点直比」列：其 ≤0 才说明端点上仍无反例。")
print("-" * 104)
print()
print("=" * 104)
print("门 4（局部极大）：Ψ(t=1) —— 从 Einstein 出发朝任一同族谱走，λ₁ 的方向导数")
print("      Ψ(1) ≥ 0 ∀方向 ⟹ **Einstein 是该族内的局部极大**（Carbotte RMP62 §II.H 自陈未证的那句）")
print("=" * 104)
print(f"  {loc_rows} 个合法行：Ψ(1)<0 的 {loc_bad} 行，全局 min Ψ(1) = {loc_min:.3e}")
if loc_bad == 0:
    print("  判词：**局部极大在扫到的全部方向上成立**（含正对照的失效带！）—— numerical-evidence。")
    print("  ⟹ 与门 3 合读：Einstein 是局部极大，但**不是**沿线段单调可达 ⟹ 失效（若有）只能是全局型。")
    print("  ⚠ 正对照在门 3 亮、在门 4 不亮 ⟹ 门 4 对「已知失效带」是**盲**的，其 PASS 不可单独承重。")
else:
    print("  判词：连局部极大都在部分方向失效 ⟹ 该处极值性**已被证伪**（须与端点直比列对照确认）。")
print("=" * 104)
