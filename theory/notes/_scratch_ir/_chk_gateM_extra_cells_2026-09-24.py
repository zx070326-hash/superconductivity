# -*- coding: utf-8 -*-
"""
gateM-cells（2026-09-24）：补 4 格证书，关掉 GATE M（`C(λ)` 括号链单调）的 3 对未定。

派活方（team-lead）消息见交付记录；消费点、预算、可达性均由
`_chk_main_gateM_budget_and_feasibility_2026-09-24.py` 算好，本脚本**不重算预算**，
只执行派活清单里的 4 格，并复用基线脚本
`_chk_main_C_lambda_certified_curve_2026-09-24.py`（v2）里的 `upper_cert`/`lower_cert`
同款判定逻辑（`certify()` / `kt_box()`），**不改基线脚本**（CLAUDE.md §2）。

━━━━━━━━ 4 格清单（λ, μ*=0.13, 侧, 目标 f） ━━━━━━━━
  1. λ=3.2953  lower   f=0.0060   关 λ 3.0 → 3.2953
  2. λ=6.0     lower   f=0.0060   关 λ 5.0 → 6.0
  3. λ=8.0     upper   f=0.0100（若过，再试 f=0.0060）  关 λ 8.0 → 10.0（硬约束 f_up ≤ 0.011872）
  4. λ=10.0    lower   f=0.0013（若 λ=8 upper 到 0.0060，可放宽到 0.0042）

━━━━━━━━ 方向（与基线脚本一致，逐字抄） ━━━━━━━━
  上侧：certify() 在 w̄ < ω̄* 给 sup<1（PASS）⟹ 真阈值 ω̄*≥w̄ ⟹ C 的**认证上界**。
  下侧：kt_box(...,-1e9) 给的 Rayleigh 商在 w̄ > ω̄* 处 mu_lo>1 ⟹ 真阈值 ω̄*≤w̄ ⟹ C 的**认证下界**。

━━━━━━━━ FAIL-stop 纪律（CLAUDE.md 单边界方向规则） ━━━━━━━━
  `ok = (len(fails)==0) and (worst<1.0)`；若 worst<1 但 ok=False，只能是 fails 非空，
  ⟹ 该 worst 是**已判定盒子上的 max**，是真 sup 的**下估**，只支持 FAIL，不得读成「差一点过」。
  某格 FAIL 就如实报 FAIL 并说明卡在哪（叶数打满/超时），不降级参数硬凑 PASS。

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
MUS = 0.13
C_INF = 0.1827262477
MAXBOX = 600            # 与基线一致


def C_of(lam, w):
    return 1.0 / (2.0 * math.pi * w * math.sqrt(lam))


def upper_cell(lam, mus, wE, fs):
    """certify() 在给定 f 列表（升序尝试，逐个报告；不早停——我们要精确知道每个目标格的结果）。"""
    out = []
    for f in fs:
        w = wE * (1.0 - f)
        Nc = int(math.floor(R * w))
        t0 = time.time()
        Rres = BB.certify(lam, w, mus, Nc, 1.0, I=60, J=400, Nv=40, L=160,
                           maxbox=MAXBOX, verbose=False)
        secs = time.time() - t0
        ok = (len(Rres["fails"]) == 0) and (Rres["worst"] < 1.0)
        row = dict(f=f, w=w, Nc=Nc, ok=ok, worst=Rres["worst"], nleaf=Rres["nleaf"],
                   nsplit=Rres["nsplit"], fails=len(Rres["fails"]), secs=secs)
        out.append(row)
        print(f"    [upper λ={lam} μ*={mus}] f={f:.4f} w̄={w:.6f} Nc={Nc} "
              f"nleaf={Rres['nleaf']} nsplit={Rres['nsplit']} fails={len(Rres['fails'])} "
              f"worst={Rres['worst']:.9f} ({secs:.0f}s) {'PASS' if ok else 'FAIL'}", flush=True)
    return out


def lower_cell(lam, mus, wE, fs, Nv=80, pad=200):
    Lr = Nv + pad
    lam_iv = iv.mpf(lam)
    out = []
    for f in fs:
        w = wE * (1.0 + f)
        Nc = int(math.floor(R * w))
        t0 = time.time()
        M2_iv = iv.mpf(lam) * iv.mpf(w) * iv.mpf(w)
        v = BB.vstar(lam, w, 0.5, 0.5, mus, Nc, Nv)
        g = BB.G_box(lam, w, (0.5, 0.5, 0.5, 0.5), Lr + Nv + 2)
        _, mlo, mhi, _ = BB.kt_box(g, lam_iv, M2_iv, mus, Nc, v, Nv, Lr, -1e9)
        secs = time.time() - t0
        ok = mlo > 1.0
        row = dict(f=f, w=w, Nc=Nc, ok=ok, mlo=mlo, secs=secs)
        out.append(row)
        print(f"    [lower λ={lam} μ*={mus}] f={f:.4f} w̄={w:.6f} Nc={Nc} mu_lo={mlo:.9f} "
              f"({secs:.1f}s) {'PASS' if ok else 'FAIL'}", flush=True)
    return out


def best_ok(rows):
    """rows 里 ok=True 的、C 最紧（|C-w 距离最小，即|f|最小）的一条；无则 None。"""
    oks = [r for r in rows if r["ok"]]
    if not oks:
        return None
    return min(oks, key=lambda r: r["f"])


def main():
    print("=" * 112)
    print("GATE M 补格：4 格清单，逐格独立算，不早停（每条都要精确报告）")
    print("=" * 112)

    results = {}

    # ── 1. λ=3.2953 lower, f=0.0060 ──
    lam = 3.2953
    wE = KT.einstein_root(lam, MUS, R, 600)
    NcE = int(math.floor(R * wE))
    lE = KT.lam_k([lam], [wE], MUS, NcE, 600)
    print(f"\n── 格 1: λ={lam} μ*={MUS}  ω̄₂^E={wE:.9f}  N_c^E={NcE}  λ₁(E,float)={lE:.9f}")
    rows1 = lower_cell(lam, MUS, wE, [0.0060])
    results[("3.2953", "lower")] = (wE, rows1)

    # ── 2. λ=6.0 lower, f=0.0060 ──
    lam = 6.0
    wE = KT.einstein_root(lam, MUS, R, 600)
    NcE = int(math.floor(R * wE))
    lE = KT.lam_k([lam], [wE], MUS, NcE, 600)
    print(f"\n── 格 2: λ={lam} μ*={MUS}  ω̄₂^E={wE:.9f}  N_c^E={NcE}  λ₁(E,float)={lE:.9f}")
    rows2 = lower_cell(lam, MUS, wE, [0.0060])
    results[("6.0", "lower")] = (wE, rows2)

    # ── 3. λ=8.0 upper, f=0.0100 then 0.0060 (只在 0.0100 PASS 后才试) ──
    lam = 8.0
    wE = KT.einstein_root(lam, MUS, R, 600)
    NcE = int(math.floor(R * wE))
    lE = KT.lam_k([lam], [wE], MUS, NcE, 600)
    print(f"\n── 格 3: λ={lam} μ*={MUS}  ω̄₂^E={wE:.9f}  N_c^E={NcE}  λ₁(E,float)={lE:.9f}")
    rows3 = upper_cell(lam, MUS, wE, [0.0100])
    if rows3[-1]["ok"]:
        rows3 += upper_cell(lam, MUS, wE, [0.0060])
    else:
        print("    f=0.0100 未过（如实报 FAIL），跳过 f=0.0060 的尝试")
    results[("8.0", "upper")] = (wE, rows3)

    # ── 4. λ=10.0 lower ──
    lam = 10.0
    wE = KT.einstein_root(lam, MUS, R, 600)
    NcE = int(math.floor(R * wE))
    lE = KT.lam_k([lam], [wE], MUS, NcE, 600)
    print(f"\n── 格 4: λ={lam} μ*={MUS}  ω̄₂^E={wE:.9f}  N_c^E={NcE}  λ₁(E,float)={lE:.9f}")
    # 若格 3 的 f=0.0060 通过，λ=8 上侧收紧 ⟹ 下侧目标可放宽到 0.0042；否则用硬目标 0.0013
    lam8_best = best_ok(rows3)
    target4 = 0.0013
    if lam8_best is not None and lam8_best["f"] <= 0.0060 + 1e-12:
        target4 = 0.0042
        print(f"    （λ=8 上侧最紧 f={lam8_best['f']:.4f} ⟹ 下侧 λ=10 目标放宽到 f={target4}）")
    rows4 = lower_cell(lam, MUS, wE, [target4])
    results[("10.0", "lower")] = (wE, rows4)

    print()
    print("=" * 112)
    print("格结果汇总（best_ok = 该格实测通过的最紧 f）")
    print("=" * 112)
    for key, (wE, rows) in results.items():
        b = best_ok(rows)
        if b is None:
            print(f"  {key}: 全部尝试 FAIL — {[(r['f'], r['ok']) for r in rows]}")
        else:
            side = key[1]
            if side == "upper":
                w = wE * (1.0 - b["f"]); C = C_of(float(key[0]), w)
                print(f"  {key}: PASS f={b['f']:.4f} w̄={w:.6f} Nc={b['Nc']} "
                      f"worst={b['worst']:.9f} C_upper={C:.9f} C/C∞={C/C_INF:.5f}")
            else:
                w = wE * (1.0 + b["f"]); C = C_of(float(key[0]), w)
                print(f"  {key}: PASS f={b['f']:.4f} w̄={w:.6f} Nc={b['Nc']} "
                      f"mu_lo={b['mlo']:.9f} C_lower={C:.9f} C/C∞={C/C_INF:.5f}")

    # ── GATE M 整表重算：用新端点替换旧端点 ──
    print()
    print("=" * 112)
    print("GATE M 整表重算（μ*=0.13, r=10）：新端点替换基线表对应行")
    print("=" * 112)

    # 基线表 C_lower/C_upper（逐字抄自 C_lambda_certified_curve_2026-09-24.md §1，μ*=0.13 11 行）
    BASE_MUS013 = {
        1.0:    dict(C_lower=0.071626887, C_upper=0.085109595),
        1.5:    dict(C_lower=0.099972499, C_upper=0.101992146),
        2.0:    dict(C_lower=0.114376567, C_upper=0.116687205),
        2.5:    dict(C_lower=0.122918849, C_upper=0.125402058),
        3.0:    dict(C_lower=0.128397265, C_upper=0.130991149),
        3.2953: dict(C_lower=0.130720384, C_upper=0.133361200),
        4.0:    dict(C_lower=0.135034160, C_upper=0.137762123),
        5.0:    dict(C_lower=0.139078844, C_upper=0.141888517),
        6.0:    dict(C_lower=0.141812277, C_upper=0.144677171),
        8.0:    dict(C_lower=0.145685231, C_upper=0.150914957),
        10.0:   dict(C_lower=0.147582982, C_upper=0.150564456),
    }
    table = {k: dict(v) for k, v in BASE_MUS013.items()}

    def apply_new(key_lam, side, wE, rows):
        b = best_ok(rows)
        if b is None:
            return False
        lam = float(key_lam)
        if side == "upper":
            w = wE * (1.0 - b["f"])
            table[lam]["C_upper"] = C_of(lam, w)
        else:
            w = wE * (1.0 + b["f"])
            table[lam]["C_lower"] = C_of(lam, w)
        return True

    updated = []
    for key, (wE, rows) in results.items():
        lam_s, side = key
        if apply_new(lam_s, side, wE, rows):
            updated.append(key)

    lams_sorted = sorted(table)
    updated_lam_side = {(float(k[0]), k[1]) for k in updated}
    print(f"\n  本次实测更新的格：{updated}")
    print(f"\n  {'λ':>8}{'C_lower':>16}{' ':>2}{'C_upper':>16}{' ':>2}   (* = 本次更新)")
    for lam in lams_sorted:
        lo_mark = "*" if (lam, "lower") in updated_lam_side else " "
        up_mark = "*" if (lam, "upper") in updated_lam_side else " "
        print(f"  {lam:>8.4f}{table[lam]['C_lower']:>16.9f}{lo_mark:>2}"
              f"{table[lam]['C_upper']:>16.9f}{up_mark:>2}")

    print()
    print("  GATE M 逐对判定：C_lower(λ_{i+1}) > C_upper(λ_i) ？")
    pairs_ok = 0
    pairs_undet = []
    for i in range(len(lams_sorted) - 1):
        li, lj = lams_sorted[i], lams_sorted[i + 1]
        clo_j = table[lj]["C_lower"]
        cup_i = table[li]["C_upper"]
        est = clo_j > cup_i
        tag = "确立" if est else "未定"
        if est:
            pairs_ok += 1
        else:
            pairs_undet.append((li, lj))
        print(f"    λ {li:g} → {lj:g}:  C_lower({lj:g})={clo_j:.9f}  vs  C_upper({li:g})={cup_i:.9f}"
              f"   差={clo_j - cup_i:+.3e}   {tag}")

    print()
    print(f"  μ*=0.13: {pairs_ok}/{len(lams_sorted)-1} 对确立严格升；仍未定：{pairs_undet}")
    print(f"  （μ*=0 的 3/3 不受本次改动影响，未重算——本次 4 格全部 μ*=0.13）")


if __name__ == "__main__":
    main()
