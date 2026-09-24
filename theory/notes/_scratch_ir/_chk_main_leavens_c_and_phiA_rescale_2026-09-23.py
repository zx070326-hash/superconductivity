# -*- coding: utf-8 -*-
"""
主会话从零独立重算（2026-09-23）。两件事，共用一套**我自己推导并搭起来**的算子
（不 import、不复制 `_chk_mustar_ceiling_roomT_2026-09-23.py` 或 2cf6342 上的脚本）：

(1) **Leavens 面积天花板常数 `c(μ*, r)`** —— 上一轮判词里唯一自陈没重算的一环。
    笔记 `ceiling/phiA_vs_leavens_roomT_crossing_2026-09-23.md` 写明：
    `c(0.13,10)=0.16574` 是从 memory 记的交叉常数 `λ*=(2C∞/c)²=4.862` **代数反解**出来的，
    本轮把它**真解一次**（对 Einstein 族做 `max_λ kTc/A`）。
    消费点：`λ_area`（300K 所需耦合的面积-天花板阈值）∝ `1/c²` ⟹ 判词对 c 的 11% 误差敏感
    （c 若 ≥0.1838 则面积阈值降到 (Φ-A) 阈值之下，胜负关系反转）。

(2) **(Φ-A) 的标度一致性**：定理原文（`coulomb/mustar_resolved_ceiling_PhiE_2026-09-21.md`
    第 442–443、490 行）把 `K_ph` 定为**无库仑、λ=1 归一**的 pencil，配 `μ̃ = μ*/λ`；
    该约定下 Tc 判据是 `λ₁ = 1/λ`。而 `_chk_mustar_ceiling_roomT_2026-09-23.py:84-94`
    把 `μ₁,μ₂` 乘了 λ（即换到 λ-标度约定）却仍用 `μ̃` 当秩一强度、且 `Δ` 未乘 λ、
    判据用 `gstar=1.0` ⟹ 两个约定混用。本脚本两版都算，并用**精确 λ₁** 做方向门。

━━━━━━━━━━ 算子的推导（analytic 部分，本脚本唯一的解析内容）━━━━━━━━━━
线性化各向同性 Eliashberg（Tc 处，Δ 通道）：
    Z_n Δ_n = πT Σ_{m∈ℤ} [λ(n−m) − μ*θ] Δ_m/|ω_m|,
    Z_n     = 1 + (πT/ω_n) Σ_{m∈ℤ} λ(n−m) sgn(ω_m),      ω_n=(2n+1)πT.
πT/|ω_m| = 1/(2m+1)。用 Δ 偶性折叠 m<0（|ω_{−m−1}|=ω_m、λ(n+m+1)）：
    RHS = Σ_{m≥0} [G(|n−m|)+G(n+m+1) − 2μ*1_{n,m<Nc}] Δ_m/(2m+1),
    Z_n = 1 + u_n² S_n,   S_n := Σ_{m≥0}[G(|n−m|) − G(n+m+1)],  u_n := (2n+1)^{−1/2}.
代入 Δ_n = y_n/u_n 并两边乘 u_n：
    Σ_m u_n u_m [G(|n−m|)+G(n+m+1)] y_m − u_n²S_n y_n − 2μ*(u_c u_cᵀ y)_n = y_n.
⟹ **对称**矩阵 M，判据 `λ₁(M) = 1`。Einstein 谱：G(k) = λ ϖ̄²/(ϖ̄²+k²)，ϖ̄ := ω_E/(2πT)。
M 对 λ 齐次（G ∝ λ，S ∝ λ），故 M = λK̂ − 2μ* u_cu_cᵀ，`λ₁(M)=1 ⟺ λ₁(K̂ − 2μ̃u_cu_cᵀ)=1/λ`。
（注：S_n 本脚本用**矩阵行和**精确求，未用 `G(0)+2Σ_{k≤n}G(k)` 的 cumsum 近似
 ——这是与仓库脚本相互独立的第二套实现，两者只应差 N 截断尾。）

面积天花板的换算（Einstein）：α²F = (λω_E/2)δ(ω−ω_E) ⟹ A := ∫α²F dω = λω_E/2；
在根处 kTc = ω_E/(2πϖ̄*) ⟹ **kTc/A = 1/(π λ ϖ̄*(λ))**，故 `c(μ*,r) = max_λ 1/(πλϖ̄*)`。

严格性：以上算子形式 analytic；本脚本一切数值 **numerical-evidence**（float64, eigh, 二分）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh

C_INF = 0.1827262477890      # kTc = C_inf sqrt(lam <w^2>) 的 lam->oo 常数（仓库既有）
HB = 64.65415134             # meV / sqrt(eV/A^2/amu):  eta/M = lam * wbar2^2 / HB^2
KB = 8.617333262e-2          # meV/K


# ────────────────────────────── 算子 ──────────────────────────────
def pieces(wbar, N):
    """返回 lam=1 归一的 K_hat（无库仑），以及 u。"""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    k = np.arange(2 * N + 2)
    G = wbar * wbar / (wbar * wbar + k * k)          # lam = 1
    Kp = G[np.abs(i - j)]
    Km = G[i + j + 1]
    S = (Kp - Km).sum(axis=1)                        # Z 的重整和（N 截断内精确）
    Khat = (Kp + Km) * np.outer(u, u) - np.diag(u * u * S)
    return Khat, u


def lam1_exact(wbar, lam, mus, r, N):
    """精确 lambda1(M)，M = lam*Khat - 2 mu* u_c u_c^T。判据 lam1 = 1。"""
    Khat, u = pieces(wbar, N)
    M = lam * Khat
    if mus > 0.0:
        Nc = int(np.floor(r * wbar))
        if Nc > 0:
            uc = np.where(np.arange(N) < Nc, u, 0.0)
            M -= 2.0 * mus * np.outer(uc, uc)
    return float(eigh(M, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])


def phiA(wbar, lam, mus, r, N, variant):
    """(Phi-A) 对 lambda1(M) 的上界，返回与判据 1.0 同标度的量（即 lam-标度）。

    variant='consistent': 严格按定理的齐次性，秩一强度 = mu*，gap = lam*gamma_hat，
                          ceiling = lam*mu1_hat - Delta，Delta = 2mu* c1^2/(1+S+)。
    variant='repo'      : 复刻 _chk_mustar_ceiling_roomT_2026-09-23.py:84-94 的混用约定。
    返回 None 表示前提 (PA) 不成立。"""
    Nc = int(np.floor(r * wbar))
    Khat, u = pieces(wbar, N)
    ev, V = eigh(Khat)
    m1h, m2h = ev[-1], ev[-2]
    uc = np.where(np.arange(N) < Nc, u, 0.0)
    c1s = float(V[:, -1] @ uc) ** 2
    un = float(uc @ uc)
    if variant == "consistent":
        nu, gam = mus, lam * (m1h - m2h)
    elif variant == "repo":
        nu, gam = mus / lam, lam * (m1h - m2h)
    else:
        raise ValueError(variant)
    den = gam - 2.0 * nu * un
    if den <= 0.0:
        return None
    sp = 2.0 * nu * (un - c1s) / den
    return lam * m1h - 2.0 * nu * c1s / (1.0 + sp)


# ─────────────────────── 求根：lam1(wbar) = 1 ───────────────────────
def roots_in_wbar(fn, gstar, lo, hi, ngrid=26, nbis=34):
    """粗扫定位**第一个**变号区间再二分。返回 (变号次数, 最小根)。
    最小根 = sup 规则（最高 Tc），与 (Phi-A) 笔记 §(P1) 一致。
    注：粗扫只为找括号；`multi_root_scan` 另行用细网格专门查多根/非单调。"""
    ws = np.linspace(lo, hi, ngrid)
    vals = np.empty(len(ws))
    for t, w in enumerate(ws):
        v = fn(w)
        vals[t] = np.nan if v is None else v - gstar
    sgn = np.sign(vals)
    cross = [t for t in range(len(ws) - 1)
             if np.isfinite(vals[t]) and np.isfinite(vals[t + 1]) and sgn[t] != sgn[t + 1]]
    if not cross:
        return 0, None
    a, b, fa = ws[cross[0]], ws[cross[0] + 1], vals[cross[0]]
    for _ in range(nbis):
        m = 0.5 * (a + b)
        fm = fn(m) - gstar
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return len(cross), 0.5 * (a + b)


def multi_root_scan(lam, mus, r, N, lo=0.05, hi=12.0, ngrid=240):
    """细网格专查 lambda1(wbar)-1 的变号次数与单调性（库仑饱和赝根的已知风险）。"""
    ws = np.linspace(lo, hi, ngrid)
    v = np.array([lam1_exact(w, lam, mus, r, N) - 1.0 for w in ws])
    ncross = int(np.sum(np.sign(v[:-1]) != np.sign(v[1:])))
    ndec = int(np.sum(np.diff(v) < 0))
    return ncross, ndec, len(ws) - 1


def c_of_lambda(lam, mus, r, N, lo=0.05, hi=12.0):
    """kTc/A = 1/(pi lam wbar*)。"""
    ncr, w = roots_in_wbar(lambda w: lam1_exact(w, lam, mus, r, N), 1.0, lo, hi)
    if w is None:
        return None, ncr
    return 1.0 / (np.pi * lam * w), ncr


def maximize_c(mus, r, N, lgrid):
    best = (-1.0, None, None)
    tab = []
    for lam in lgrid:
        c, ncr = c_of_lambda(lam, mus, r, N)
        tab.append((lam, c, ncr))
        if c is not None and c > best[0]:
            best = (c, lam, ncr)
    return best, tab


def refine_max(mus, r, N, l0, half, tol=2e-5):
    """黄金分割求 max_lam c(lam)（c 在峰附近单峰，已由粗扫确认）。"""
    gr = 0.5 * (np.sqrt(5.0) - 1.0)
    a, b = l0 - half, l0 + half
    x1, x2 = b - gr * (b - a), a + gr * (b - a)
    f = lambda l: (c_of_lambda(l, mus, r, N)[0] or -1.0)
    f1, f2 = f(x1), f(x2)
    while b - a > tol:
        if f1 < f2:
            a, x1, f1 = x1, x2, f2
            x2 = a + gr * (b - a)
            f2 = f(x2)
        else:
            b, x2, f2 = x2, x1, f1
            x1 = b - gr * (b - a)
            f1 = f(x1)
    l = 0.5 * (a + b)
    return l, f(l)


def phiA_unit(wbar, lam, mus, r, N):
    """定理**原文约定**（`mustar_resolved_ceiling_PhiE_2026-09-21.md` L442-443, L412）：
    K_ph 取 lam=1 归一的 K_hat，秩一强度取 mu~ = mu*/lam，前提 (PA) gamma_hat > 2 mu~||u||^2。
    返回 (ceiling_hat, gstar_hat=1/lam)。判据 lambda1(K_hat - 2mu~ uu^T) = 1/lam。
    这是**独立于** phiA(...,'consistent') 的第二条实现路径；两者必须满足
    lam * ceiling_hat == consistent（齐次性）——这正是本脚本的 GATE H。"""
    Nc = int(np.floor(r * wbar))
    Khat, u = pieces(wbar, N)
    ev, V = eigh(Khat)
    m1h, m2h = ev[-1], ev[-2]
    uc = np.where(np.arange(N) < Nc, u, 0.0)
    c1s = float(V[:, -1] @ uc) ** 2
    un = float(uc @ uc)
    mut, gamh = mus / lam, m1h - m2h
    den = gamh - 2.0 * mut * un
    if den <= 0.0:
        return None, 1.0 / lam
    sp = 2.0 * mut * (un - c1s) / den
    return m1h - 2.0 * mut * c1s / (1.0 + sp), 1.0 / lam


def gate_homogeneity(mus, r, N):
    """GATE H：定理原文约定 (lam=1 归一, 判据 1/lam) 与我用的 lam-标度约定必须逐点等价。
    若不等价，说明我对 (Phi-A) 的重标度读错了 —— 这条门是对**我自己这次订正**的独立复核。"""
    print("GATE H  齐次性：lam * [定理原文约定的 ceiling] == [我的 lam-标度 ceiling] ?")
    print(f"{'wbar':>6}{'lam':>6}{'lam*ceil_hat':>15}{'consistent':>13}{'|差|':>12}"
          f"{'gstar_hat':>12}{'判据等价':>10}")
    worst = 0.0
    for wbar in (0.5, 0.7, 1.0, 1.5, 2.2, 3.0):
        for lam in (2.0, 3.0, 4.0):
            ch, gs = phiA_unit(wbar, lam, mus, r, N)
            cc = phiA(wbar, lam, mus, r, N, "consistent")
            if ch is None or cc is None:
                print(f"{wbar:6.2f}{lam:6.1f}{'None':>15}{'None' if cc is None else '%.6f'%cc:>13}"
                      f"{'-':>12}{gs:12.6f}{'(PA)失效':>10}")
                continue
            d = abs(lam * ch - cc)
            worst = max(worst, d)
            # 判据等价性：ceil_hat >= 1/lam  <=>  lam*ceil_hat >= 1
            eq = ((ch >= gs) == (cc >= 1.0))
            print(f"{wbar:6.2f}{lam:6.1f}{lam*ch:15.9f}{cc:13.9f}{d:12.2e}"
                  f"{gs:12.6f}{'OK' if eq else '**NO**':>10}")
    print(f"        最劣绝对差 = {worst:.3e}  "
          f"{'(机器精度内等价 -> 我的重标度读法成立)' if worst < 1e-9 else '** 不等价, 我的读法有问题 **'}")


# ────────────────────────────── 主流程 ──────────────────────────────
def main():
    N = 400
    print("=" * 92)
    print("GATE A  mu*=0 盲测锚点：Leavens 1975 印刷值 c(0)=0.2309, omega_E/A=1.750")
    print("        （我的求解器此前独立给 0.230868 / 1.7488；本脚本是第三套独立实现）")
    (cb, lb, _), tab = maximize_c(0.0, 0.0, N, np.arange(0.7, 1.85, 0.05))
    lstar, cstar = refine_max(0.0, 0.0, N, lb, 0.06)
    print(f"        c(0) = {cstar:.6f}   (印刷 0.2309, 相对差 {abs(cstar-0.2309)/0.2309*100:.3f}%)")
    print(f"        argmax lam = {lstar:.5f} -> omega_E/A = 2/lam = {2.0/lstar:.4f} "
          f"(印刷 1.750, 相对差 {abs(2.0/lstar-1.750)/1.750*100:.3f}%)")
    for NN in (200, 800):
        _, cN = refine_max(0.0, 0.0, NN, lstar, 0.05)
        print(f"        N-收敛: N={NN:4d} -> c(0)={cN:.6f}  (N=400 差 {abs(cN-cstar):.2e})")
    print()

    print("=" * 92)
    print("STEP 1  上一轮唯一没重算的一环：c(mu*=0.13, r=10)")
    mus, r = 0.13, 10.0
    (cb2, lb2, ncr2), tab2 = maximize_c(mus, r, N, np.arange(0.7, 2.65, 0.05))
    l2, c2 = refine_max(mus, r, N, lb2, 0.06)
    print(f"        粗扫 argmax lam={lb2:.2f}  c={cb2:.6f}   (根的变号次数={ncr2})")
    print(f"        精化       lam={l2:.5f}  c(0.13,10) = {c2:.6f}")
    c_backed = 2.0 * C_INF / np.sqrt(4.862)
    print(f"        上一轮从 memory 反解的值            = {c_backed:.6f}"
          f"   (相对差 {abs(c2-c_backed)/c_backed*100:.2f}%)")
    for NN in (200, 800):
        _, cN = refine_max(mus, r, NN, l2, 0.05)
        print(f"        N-收敛: N={NN:4d} -> c={cN:.6f}  (N=400 差 {abs(cN-c2):.2e})")
    print("        多根/单调性细扫（库仑饱和赝根风险，memory: coulomb-saturation-creates-a-pseudo-root）：")
    for ll in (l2, 1.5, 2.5):
        nc, nd, nt = multi_root_scan(ll, mus, r, N)
        print(f"          lam={ll:.4f}: lambda1(wbar)-1 变号 {nc} 次, 下降步 {nd}/{nt}"
              f"  {'(单调增, 单根)' if (nc == 1 and nd == 0) else '** 非单调/多根, 需人工看 **'}")
    print()

    # ── 300K 阈值 ──
    eta_over_M = 9.172
    pref = HB * np.sqrt(eta_over_M)          # meV;  omega_E = pref/sqrt(lam)
    print("=" * 92)
    print(f"STEP 2  300K 阈值（eta/M={eta_over_M}, mu*={mus}, r={r}）")
    print("        面积天花板: kTc <= c*A = c*lam*omega_E/2 = c*pref*sqrt(lam)/2")

    def lam_area(c):
        return (2.0 * 300.0 * KB / (c * pref)) ** 2

    print(f"        用我重算的 c={c2:.6f}     -> lam_area = {lam_area(c2):.4f}")
    print(f"        用上一轮反解的 c={c_backed:.6f} -> lam_area = {lam_area(c_backed):.4f}"
          f"   (上一轮报 2.538)")
    print()

    print("=" * 92)
    print("STEP 3  (Phi-A) 方向门：ceiling >= 精确 lambda1 ？ 两个标度约定都测")
    print(f"{'wbar':>6}{'lam':>6}{'exact':>12}{'consistent':>13}{'repo':>12}"
          f"{'c_ok':>7}{'r_ok':>7}{'紧度 c/r':>10}")
    bad_c = bad_r = 0
    for wbar in (0.5, 0.7, 1.0, 1.5, 2.2, 3.0):
        for lam in (1.5, 2.0, 3.0, 4.0):
            ex = lam1_exact(wbar, lam, mus, r, N)
            gc = phiA(wbar, lam, mus, r, N, "consistent")
            gr = phiA(wbar, lam, mus, r, N, "repo")
            okc = (gc is not None) and (gc >= ex - 1e-9)
            okr = (gr is not None) and (gr >= ex - 1e-9)
            bad_c += (gc is not None) and (not okc)
            bad_r += (gr is not None) and (not okr)
            tight = ((gc - ex) / (gr - ex)) if (gc is not None and gr is not None
                                                and gr - ex > 1e-12) else float("nan")
            print(f"{wbar:6.2f}{lam:6.1f}{ex:12.6f}"
                  f"{('None' if gc is None else '%.6f' % gc):>13}"
                  f"{('None' if gr is None else '%.6f' % gr):>12}"
                  f"{'OK' if okc else 'VIOL':>7}{'OK' if okr else 'VIOL':>7}{tight:10.4f}")
    print(f"        方向违反：consistent {bad_c}/24, repo {bad_r}/24")
    print()

    print("=" * 92)
    print("STEP 4  (Phi-A) 的 300K lambda 阈值（两个约定）")

    def Tc_phiA(lam, variant):
        wE = pref / np.sqrt(lam)
        _, w = roots_in_wbar(lambda w: phiA(w, lam, mus, r, N, variant), 1.0, 0.05, 12.0)
        if w is None:
            return None
        return wE / (2.0 * np.pi * w) / KB

    out = {}
    for variant in ("repo", "consistent"):
        print(f"  -- {variant} --")
        print(f"{'lam':>7}{'Tc_ceil[K]':>13}")
        for lam in np.arange(1.5, 4.01, 0.25):
            t = Tc_phiA(lam, variant)
            print(f"{lam:7.2f}{('%.2f' % t) if t else 'None':>13}")
        lo, hi = 1.2, 9.0
        f = lambda l: (Tc_phiA(l, variant) or -1e9) - 300.0
        if f(lo) < 0 < f(hi):
            for _ in range(45):
                m = 0.5 * (lo + hi)
                if f(m) < 0:
                    lo = m
                else:
                    hi = m
            out[variant] = 0.5 * (lo + hi)
            print(f"  ==> lam_ceil({variant}) = {out[variant]:.4f}")
        else:
            out[variant] = None
            print(f"  ==> 扫描区间内无变号 (f(lo)={f(lo):.2f}, f(hi)={f(hi):.2f})")
    print()

    print("=" * 92)
    print("判词")
    la = lam_area(c2)
    print(f"  面积天花板 (Leavens, 我重算的 c={c2:.5f})  lam_area   = {la:.4f}")
    for v in ("repo", "consistent"):
        if out[v] is not None:
            rel = "面积更紧（(Phi-A) 输）" if la > out[v] else "(Phi-A) 更紧（赢）"
            print(f"  (Phi-A) [{v:10s}]                       lam_ceil  = {out[v]:.4f}   -> {rel}")


if __name__ == "__main__":
    main()
