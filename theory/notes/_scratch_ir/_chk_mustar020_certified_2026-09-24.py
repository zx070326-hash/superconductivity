# -*- coding: utf-8 -*-
"""
派活：把 μ*=0.20 那一格从 float64 抬到双侧认证括号（team-lead 消息，2026-09-24）。

沿用 `_chk_main_C_lambda_certified_curve_2026-09-24.py` 的同一套管线
（BB.certify / BB.kt_box / BB.vstar / BB.G_box / KT.einstein_root / KT.lam_k），
只把 μ* 固定在 0.20、λ 取 [2.0, 3.2953, 4.0, 6.0]（与已有 μ*=0 spot-check 同一组 λ）。

━━━━━━━━ 消费点与预算（团队长的原话，抄录不改写） ━━━━━━━━
头条表「我们对 Leavens 的优势随 μ* 增大」逐 μ* 比值：1.0152/1.1853/1.2982/1.8156  【2026-09-24 订正：该列现为 1.7967，原 1.8156 的 λ_ceil 取了最大根】
（μ*=0/0.10/0.13/0.20）。μ*=0.20 这一列现在是 numerical-evidence。
预算：整列删掉 ⟹ 仍余三列，「随 μ* 单调升」定性结论不翻转，只是覆盖上限从 0.20 缩到 0.13。
⟹ 本活买的是覆盖范围（物理上限 μ*≈0.2），不买通过；优先级低于 Leavens 侧认证。
不许把交付写成「关闭了某个缺口」。

━━━━━━━━ 方向（与基线脚本相同，未新推导） ━━━━━━━━
GATE D：ω̄ 增大 ⟹ λ₁ 增大 ⟹
  certify() 在 w<ω̄^E 给 sup<1（PASS）⟹ 真阈值 ω̄*≥该 w ⟹ C 的**认证上界**；
  kt_box(...,gamma=-1e9) 在 w>ω̄^E 给 mu_lo>1（PASS）⟹ 真阈值 ω̄*≤该 w ⟹ C 的**认证下界**。

跑法：PYTHONIOENCODING=utf-8 OMP_NUM_THREADS=1 python -u <本脚本> > <日志> 2>&1
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import time
import math
from mpmath import iv

_d = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, name + ".py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


BB = _load("_chk_main_box_branch_bound_certificate_2026-09-23")
KT = BB.KT

R = 10.0
C_INF = 0.1827262477
MUS = 0.20
LAMS = [2.0, 3.2953, 4.0, 6.0]
TIME_BUDGET_SECS = 42 * 60
MAXBOX = 600
OFFS_DESC = [0.30, 0.15, 0.06, 0.025, 0.01]


def C_of(lam, w):
    return 1.0 / (2.0 * math.pi * w * math.sqrt(lam))


def upper_cert(lam, mus, wE, t_deadline):
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
            best = dict(status="OK", w=w, Nc=Nc, C=C_of(lam, w), tried=tried, last_pass_f=f)
        else:
            break
    if best is None:
        return dict(status="EXHAUSTED" if tried else "BUDGET", tried=tried)
    return best


def lower_cert(lam, mus, wE, t_deadline, Nv=80, pad=200):
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
            best = dict(status="OK", w=w, Nc=Nc, C=C_of(lam, w), tried=tried, last_pass_f=f)
        else:
            break
    if best is None:
        return dict(status="EXHAUSTED" if tried else "BUDGET", tried=tried)
    return best


def gate_Nc_no_jump():
    """派活坑 1：N_c = floor(R*w) 台阶——收紧前后（Einstein 点 vs 各 f 试探点）N_c 不许跳档
    未经报告地跨越导致的『换模型』。这里只做事后核对：在主循环里逐格打印 N_c，
    本 gate 单独在 Einstein 锚点上复核 floor 的邻域没有意外跳变（数值证据，非声音性）。"""
    print("=" * 110)
    print("GATE J  N_c=floor(R*w) 邻域跳档检查（Einstein 锚点 ±1e-6 扰动）")
    for lam in LAMS:
        wE = KT.einstein_root(lam, MUS, R, 600)
        NcE = int(math.floor(R * wE))
        Nc_lo = int(math.floor(R * wE * (1 - 1e-6)))
        Nc_hi = int(math.floor(R * wE * (1 + 1e-6)))
        print(f"  λ={lam}: ω̄₂^E={wE:.9f} N_c(E)={NcE} N_c(E·(1-1e-6))={Nc_lo} N_c(E·(1+1e-6))={Nc_hi}")
    print("  （本 gate 只核对 Einstein 点附近；各 f 格自身的 N_c 见主循环逐行打印）")
    print()


def gate_freedom_can_fail():
    """派活坑 2：任何自检门要先问『它可能失败吗』。μ*=0.20≠0 ⟹ 2μ*≠0 ⟹ N_c 代数上确实
    进入 λ_k（与 μ*=0 的退化情形不同），此门在本设定下**可能失败**，计入通过率合法。"""
    print("=" * 110)
    print("GATE F  自检门可判性：μ*=0.20 下 N_c 是否代数进入 λ_k（对照 μ*=0 的退化）")
    lam, w = 4.0, 0.4775933
    vals_mu020 = [KT.lam_k([lam], [w], 0.20, Nc, 900) for Nc in (0, 1, 3, 5, 10, 50)]
    same020 = all(abs(v - vals_mu020[0]) < 1e-13 for v in vals_mu020)
    vals_mu0 = [KT.lam_k([lam], [w], 0.0, Nc, 900) for Nc in (0, 1, 3, 5, 10, 50)]
    same0 = all(abs(v - vals_mu0[0]) < 1e-13 for v in vals_mu0)
    print(f"  μ*=0.20: N_c∈{{0,1,3,5,10,50}} → λ₁∈{['%.9f'%v for v in vals_mu020]}  逐位相同={same020}")
    print(f"  μ*=0.00(对照): 逐位相同={same0}")
    print(f"  ⟹ μ*=0.20 处该门『扫 N_c 看结果不变』{'不可能失败——不计入通过率' if same020 else '可能失败，本次通过——计入通过率合法'}")
    print()


def main():
    t_start = time.time()
    t_deadline = t_start + TIME_BUDGET_SECS
    gate_Nc_no_jump()
    gate_freedom_can_fail()

    results = {}
    alarms = []

    print("=" * 110)
    print(f"主循环：μ*={MUS} 固定，逐 λ∈{LAMS} 求 C 的认证上/下界")
    for lam in LAMS:
        wE = KT.einstein_root(lam, MUS, R, 600)
        NcE = int(math.floor(R * wE))
        lE = KT.lam_k([lam], [wE], MUS, NcE, 600)
        print(f"── λ={lam} μ*={MUS}  ω̄₂^E={wE:.9f}  N_c^E={NcE}  λ₁(E,float)={lE:.9f}", flush=True)
        if time.time() > t_deadline:
            print("   [BUDGET] 总预算耗尽，本点及之后未算")
            results[lam] = (wE, NcE, lE, dict(status="BUDGET_SKIP"), dict(status="BUDGET_SKIP"))
            continue
        up = upper_cert(lam, MUS, wE, t_deadline)
        lo = lower_cert(lam, MUS, wE, t_deadline)
        results[lam] = (wE, NcE, lE, up, lo)
        if up.get("status") == "OK" and up["C"] >= C_INF:
            alarms.append(("upper", lam, up))
            print(f"   !! CONSISTENCY ALARM: 上界 C={up['C']:.9f} >= C_INF={C_INF}")
        if lo.get("status") == "OK" and lo["C"] > C_INF:
            alarms.append(("lower", lam, lo))
            print(f"   !! CONSISTENCY ALARM: 下界 C={lo['C']:.9f} > C_INF={C_INF}")
        print()

    print("=" * 110)
    print("汇总表（含收紧深度：最后一次 PASS 的 f，未到 0.010 的显式标出）")
    hdr = (f"  {'λ':>8}{'μ*':>6}{'ω̄₂^E(f)':>13}{'N_c^E':>7}"
           f"{'C_upper':>14}{'C_up/C∞':>10}{'w_up':>10}{'f_up(last-PASS)':>16}"
           f"{'C_lower':>14}{'C_lo/C∞':>10}{'w_lo':>10}{'f_lo(last-PASS)':>16}")
    print(hdr)
    for lam, (wE, NcE, lE, up, lo) in results.items():
        def fmt(d):
            if d.get("status") == "OK":
                f_last = d.get("last_pass_f")
                flag = f"{f_last:.3f}" + ("" if f_last <= 0.01 + 1e-12 else "  <-- 未达0.010!")
                return f"{d['C']:>14.9f}{d['C']/C_INF:>10.5f}{d['w']:>10.6f}{flag:>16}"
            return f"{d.get('status',''):>14}{'':>10}{'':>10}{'':>16}"
        print(f"  {lam:>8.4f}{MUS:>6.2f}{wE:>13.9f}{NcE:>7}" + fmt(up) + fmt(lo))
    print()
    print(f"  总耗时 {time.time()-t_start:.0f}s（预算 {TIME_BUDGET_SECS}s）  C_INF={C_INF}")

    print()
    print("一致性告警清单：")
    if not alarms:
        print("  （无——全部认证点都在 C_INF 之下，与已并入的 P4_conjecture2.md 定理相容）")
    else:
        for side, lam, d in alarms:
            print(f"  !! side={side} λ={lam} C={d['C']:.9f} w={d['w']:.9f} Nc={d['Nc']}")

    print()
    print("未算/未达清单：")
    any_missing = False
    for lam, (wE, NcE, lE, up, lo) in results.items():
        for side, d in (("upper", up), ("lower", lo)):
            if d.get("status") != "OK":
                any_missing = True
                print(f"  λ={lam} 侧={side} 状态={d.get('status')} "
                      f"已试偏移={[t['f'] for t in d.get('tried',[])]}")
            elif d.get("last_pass_f", 0) > 0.01 + 1e-12:
                any_missing = True
                print(f"  λ={lam} 侧={side} 状态=OK 但收紧深度只到 f={d['last_pass_f']:.3f}（未达0.010）")
    if not any_missing:
        print("  （无）")

    print()
    print("严格性标级：本脚本产出的 C_upper/C_lower 为 certificate 级（mpmath.iv 外向舍入 +"
          " Kato-Temple/Schur 解析尾界，见基线脚本头注）；C(λ,μ*=0.20,E,float) 为 numerical-evidence。"
          " 本活的消费点：把头条表最后一列（μ*=0.20 比 1.7967，订正前印作 1.8156）从 numerical-evidence 抬到 certificate，"
          " 不改变『随 μ* 单调升』定性结论（该结论删去此列仍成立，见脚本头注的预算记录）。")


if __name__ == "__main__":
    main()
