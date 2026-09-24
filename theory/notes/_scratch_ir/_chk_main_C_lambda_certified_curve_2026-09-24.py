# -*- coding: utf-8 -*-
"""
v2（2026-09-24）：吸收团队长两条追加消息后的重写。

变更历史（诚实记录，不删旧版本的教训）：
  v1（本目录同名文件的更早版本，已被本文件覆盖）用「从最小偏移升序试探」，
  在 λ=1（弱耦合端，ω̄₂^E=2.2，N_c=22）的第一个 certify() 调用上就吃满了 22 分钟总预算——
  offset 越小越接近阈值、盒子数指数增长，从最难的一侧起步是设计错误。
  **教训**：证书搜索必须从「大偏移=容易过」起步，早停在第一次失败，不能反过来。

━━━━━━━━ 标度关系（复核见 GATE S1，末尾附代数推导，analytic） ━━━━━━━━
kTc = C(λ,μ*)·√M2^phys，  C(λ,μ*) := 1/(2π·ω̄₂_dimless(λ,μ*)·√λ)。
ω̄₂_dimless(λ,μ*) 是无量纲算子 λ₁=1 的根（build/Gk/kt_box 只吃 (λ_a,ω̄_a,μ*,N_c)，
无任何显式 T——GATE S1 逐源码签名核对）。

━━━━━━━━ 消费点（团队长两条订正后定稿） ━━━━━━━━
`theory/derivations/P4_conjecture2.md`（INDEX.md 第 307 行，方向 A，已并入主线）已经证了
    Tc < C_∞·√(λ⟨ω²⟩),  C_∞ = (1/2π)√𝔤(2) = 0.1827262477   （analytic，对全部 P、全部 μ*≥0，单边）
这条对 λ、μ* **都不依赖**（自由的、非锐的覆盖：库仑 PSD ⟹ 单调性免费给出 μ*≥0）。
本脚本产出的 C(λ,μ*) 是**分辨 λ 和 μ* 的**，其价值是**紧度比 C(λ,μ*)/C_∞**（越小越紧）
和**μ* 买到的锐化量 C(λ,0)−C(λ,μ*)**。

**一致性门（有牙，不是判词的浪漫版）**：
  · 任何认证**上界**侧 C ≥ C_∞ ⟹ 与仓库内已并入的 analytic 定理**直接冲突** ⟹ 立刻停、报告。
  · 任何认证**下界**侧（Rayleigh μ_lo>1 侧）C > C_∞ ⟹ 同样冲突 ⟹ 立刻停、报告。
脚本对每个点都做这两条检查，触发时把该点单独打印到 "CONSISTENCY ALARM" 段。

━━━━━━━━ 方向 ━━━━━━━━
GATE D（沿用 v1 结果，未变）：ω̄ 增大 ⟹ λ₁ 增大，全部 λ 一致 ⟹
  certify() 在 ω̄<ω̄^E 给 sup<1（PASS 方向）⟹ 真阈值 ω̄*≥该 ω̄ ⟹ C 的**认证上界**；
  Rayleigh μ_lo>1（FAIL 方向）在 ω̄>ω̄^E 处成立 ⟹ 真阈值 ω̄*≤该 ω̄ ⟹ C 的**认证下界**。

跑法：PYTHONIOENCODING=utf-8 OMP_NUM_THREADS=1 python -u <本脚本> > <日志> 2>&1
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import time
import math
import numpy as np
from mpmath import iv

_d = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, name + ".py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


BB = _load("_chk_main_box_branch_bound_certificate_2026-09-23")
KT = BB.KT

R = 10.0
C_INF = 0.1827262477                       # KAY Conjecture 2 常数 = 仓库 P4_conjecture2.md 的 C_∞（analytic）
LAMS_FULL = [1.0, 1.5, 2.0, 2.5, 3.0, 3.2953, 4.0, 5.0, 6.0, 8.0, 10.0]
LAMS_MUS0 = [2.0, 3.2953, 4.0, 6.0]         # μ*=0 spot-check（团队长：3–4 个 λ 够）
TIME_BUDGET_SECS = 42 * 60
MAXBOX = 600                                 # 单次 certify() 的盒预算上限（控速，非声音性相关）
OFFS_DESC = [0.30, 0.15, 0.06, 0.025, 0.01]  # 降序：先易后难，第一次失败即停


def C_of(lam, w):
    return 1.0 / (2.0 * math.pi * w * math.sqrt(lam))


def upper_cert(lam, mus, wE, t_deadline):
    """certify() 沿 OFFS_DESC 降序（易→难）试探，第一次失败立即停（早停控速）。"""
    best = None; tried = []
    for f in OFFS_DESC:
        if time.time() > t_deadline:
            break
        w = wE * (1.0 - f)
        Nc = int(math.floor(R * w))
        t0 = time.time()
        Rres = BB.certify(lam, w, mus, Nc, 1.0, I=60, J=400, Nv=40, L=160,
                           maxbox=MAXBOX, verbose=False)
        secs = time.time() - t0
        ok = (len(Rres["fails"]) == 0) and (Rres["worst"] < 1.0)
        tried.append(dict(f=f, w=w, Nc=Nc, ok=ok, worst=Rres["worst"],
                           nleaf=Rres["nleaf"], secs=secs))
        print(f"    [upper μ*={mus}] f={f:.3f} w̄={w:.6f} Nc={Nc} nleaf={Rres['nleaf']} "
              f"worst={Rres['worst']:.9f} ({secs:.0f}s) {'PASS' if ok else 'FAIL-stop'}", flush=True)
        if ok:
            best = dict(status="OK", w=w, Nc=Nc, C=C_of(lam, w), tried=tried)
        else:
            break
    if best is None:
        return dict(status="EXHAUSTED" if tried else "BUDGET", tried=tried)
    return best


def lower_cert(lam, mus, wE, t_deadline, Nv=80, pad=200):
    """Rayleigh μ_lo（快）同样降序早停。"""
    Lr = Nv + pad
    lam_iv = iv.mpf(lam)
    best = None; tried = []
    for f in OFFS_DESC:
        if time.time() > t_deadline:
            break
        w = wE * (1.0 + f)
        Nc = int(math.floor(R * w))
        t0 = time.time()
        M2_iv = iv.mpf(lam) * iv.mpf(w) * iv.mpf(w)
        v = BB.vstar(lam, w, 0.5, 0.5, mus, Nc, Nv)
        g = BB.G_box(lam, w, (0.5, 0.5, 0.5, 0.5), Lr + Nv + 2)
        _, mlo, mhi, _ = BB.kt_box(g, lam_iv, M2_iv, mus, Nc, v, Nv, Lr, -1e9)
        secs = time.time() - t0
        ok = mlo > 1.0
        tried.append(dict(f=f, w=w, Nc=Nc, ok=ok, mlo=mlo, secs=secs))
        print(f"    [lower μ*={mus}] f={f:.3f} w̄={w:.6f} Nc={Nc} mu_lo={mlo:.9f} "
              f"({secs:.1f}s) {'PASS' if ok else 'FAIL-stop'}", flush=True)
        if ok:
            best = dict(status="OK", w=w, Nc=Nc, C=C_of(lam, w), tried=tried)
        else:
            break
    if best is None:
        return dict(status="EXHAUSTED" if tried else "BUDGET", tried=tried)
    return best


def gate_S1():
    print("=" * 110)
    print("GATE S1  标度关系前提：build/Gk/kt_box 的入参只有 (λ_a, ω̄_a[无量纲], μ*, N_c)")
    import inspect
    for fn in (KT.Gk, KT.build, BB.kt_box):
        print(f"  {fn.__name__}{inspect.signature(fn)}")
    print("  ⟹ 无任何显式 T / 物理常数；标度关系 kTc=C(λ,μ*)√M2^phys 为代数恒等式（analytic）。")
    print()


def gate_mus0_Nc_independence():
    print("=" * 110)
    print("GATE N0  μ*=0 时 N_c 是否真的不进判据（不从代码结构推，实测）")
    lam, w = 4.0, 0.4775933
    vals = [KT.lam_k([lam], [w], 0.0, Nc, 900) for Nc in (0, 1, 3, 5, 10, 50)]
    same = all(abs(v - vals[0]) < 1e-13 for v in vals)
    print(f"  λ={lam} ω̄={w}: N_c∈{{0,1,3,5,10,50}} 给 λ₁∈{['%.15f'%v for v in vals]}")
    print(f"  ⟹ 逐位相同：{same}（μ*=0 时 2μ*=0 系数直接使库仑项消失，与 N_c 无关，确认）")
    print()


def gate_anchors():
    print("=" * 110)
    print("GATE A  复核团队长给的两个锚点（float64, Einstein 成员, N=900，独立算）")
    for lam, mus, tag in [(4.0, 0.0, "μ*=0"), (3.2953, 0.13, "μ*=0.13")]:
        wE = KT.einstein_root(lam, mus, R, 900)
        Nc = int(math.floor(R * wE))
        C = C_of(lam, wE)
        print(f"  λ={lam} {tag}: ω̄₂={wE:.7f} N_c={Nc} C={C:.7f} C/C∞={C/C_INF:.5f}")
    print()


def run_side(lam, mus, wE, t_deadline, side):
    if side == "upper":
        return upper_cert(lam, mus, wE, t_deadline)
    return lower_cert(lam, mus, wE, t_deadline)


def main():
    t_start = time.time()
    t_deadline = t_start + TIME_BUDGET_SECS
    gate_S1()
    gate_mus0_Nc_independence()
    gate_anchors()

    combos = [(lam, 0.13) for lam in LAMS_FULL] + [(lam, 0.0) for lam in LAMS_MUS0]
    results = {}
    alarms = []

    print("=" * 110)
    print("主循环：逐 (λ,μ*) 求 C 的认证上/下界（早停降序，先覆盖全表再各自尝试收紧）")
    for lam, mus in combos:
        wE = KT.einstein_root(lam, mus, R, 600)
        NcE = int(math.floor(R * wE))
        lE = KT.lam_k([lam], [wE], mus, NcE, 600)
        print(f"── λ={lam} μ*={mus}  ω̄₂^E={wE:.9f}  N_c^E={NcE}  λ₁(E,float)={lE:.9f}", flush=True)
        if time.time() > t_deadline:
            print("   [BUDGET] 总预算耗尽，本点及之后未算")
            results[(lam, mus)] = (wE, NcE, lE, dict(status="BUDGET_SKIP"), dict(status="BUDGET_SKIP"))
            continue
        up = upper_cert(lam, mus, wE, t_deadline)
        lo = lower_cert(lam, mus, wE, t_deadline)
        results[(lam, mus)] = (wE, NcE, lE, up, lo)
        if up.get("status") == "OK" and up["C"] >= C_INF:
            alarms.append(("upper", lam, mus, up))
            print(f"   !! CONSISTENCY ALARM: 上界 C={up['C']:.9f} >= C_INF={C_INF}")
        if lo.get("status") == "OK" and lo["C"] > C_INF:
            alarms.append(("lower", lam, mus, lo))
            print(f"   !! CONSISTENCY ALARM: 下界 C={lo['C']:.9f} > C_INF={C_INF}")
        print()

    print("=" * 110)
    print("汇总表")
    hdr = (f"  {'λ':>8}{'μ*':>6}{'ω̄₂^E(f)':>13}{'N_c^E':>7}"
           f"{'C_upper':>14}{'C_up/C∞':>10}{'w_up':>10}"
           f"{'C_lower':>14}{'C_lo/C∞':>10}{'w_lo':>10}")
    print(hdr)
    for (lam, mus), (wE, NcE, lE, up, lo) in results.items():
        def fmt(d):
            if d.get("status") == "OK":
                return f"{d['C']:>14.9f}{d['C']/C_INF:>10.5f}{d['w']:>10.6f}"
            return f"{d.get('status',''):>14}{'':>10}{'':>10}"
        print(f"  {lam:>8.4f}{mus:>6.2f}{wE:>13.9f}{NcE:>7}" + fmt(up) + fmt(lo))
    print()
    print(f"  总耗时 {time.time()-t_start:.0f}s（预算 {TIME_BUDGET_SECS}s）  C_INF={C_INF}")

    print()
    print("一致性告警清单：")
    if not alarms:
        print("  （无——全部认证点都在 C_INF 之下，与已并入的 P4_conjecture2.md 定理相容）")
    else:
        for side, lam, mus, d in alarms:
            print(f"  !! side={side} λ={lam} μ*={mus} C={d['C']:.9f} w={d['w']:.9f} Nc={d['Nc']}")

    print()
    print("未算/未达清单：")
    any_missing = False
    for (lam, mus), (wE, NcE, lE, up, lo) in results.items():
        for side, d in (("upper", up), ("lower", lo)):
            if d.get("status") != "OK":
                any_missing = True
                print(f"  λ={lam} μ*={mus} 侧={side} 状态={d.get('status')} "
                      f"已试偏移={[t['f'] for t in d.get('tried',[])]}")
    if not any_missing:
        print("  （无）")


if __name__ == "__main__":
    main()
