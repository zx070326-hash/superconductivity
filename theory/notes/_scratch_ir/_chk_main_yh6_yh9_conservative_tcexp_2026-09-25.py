#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""主会话从零复核 hydride-m2-hunt 的 YH6/YH9 天花板检验（2026-09-25）。

被复核的交付：theory/notes/ceiling/hydride_ceiling_r_h3s_yh6_yh9_2026-09-25.md
本脚本不 import 交付方的脚本，输入只有 Ding 2025 Table I 的五行原始数字。

三件事：
 (A) 重算 ceiling / r，并把 Tc_exp 从 fail-open 的 220/240 换成保守的 224/243；
 (B) 把 r_calc（定理自检门）与 r_exp（外部效度）分开计数；
 (C) 量出 M2 跨方法跨度，按「两端是否共享出处」分组——检验副发现「≤2.6%」。

严格性：全部 numerical-evidence（λ/ν̄₂/Tc 均为外部已发表数）。
C_INF 本身是 analytic。ν̄₂ = <ω²>^{1/2} 是 analytic（见笔记 §2 订正框），本脚本不重证。
"""
import math

C_INF = 0.1827262477            # = sqrt(g(2))/(2*pi), analytic

# --- Ding, Chen, Shi, arXiv:2507.03383v2 (2025) TABLE I ------------------------
# (system, method, lambda(0), nu2_K, Tc_paper_K, Tc_paper_renorm_K or None, spectrum_source)
TABLE_I = [
    ("YH6", "DFPT",  2.534, 1399, 270, 248, "Ref[11] Troyan2021"),
    ("YH6", "SSCHA", 1.779, 1662, 245, 218, "Ref[11] Troyan2021"),
    ("YH6", "SPIA",  1.712, 1738, 212, None, "this work (PIMD+MLFF)"),
    ("YH9", "DFPT",  3.038, 1308, 252, 236, "Ref[20] Du2023 @260GPa"),
    ("YH9", "SPIA",  1.749, 1905, 240, None, "this work (PIMD+MLFF) @255GPa"),
]

# 实验 Tc：两套，方向相反
TCEXP_FAILOPEN = {"YH6": 220.0, "YH9": 240.0}   # Ding 正文约整转述；同时是 PIMD 模拟温度
TCEXP_CONSERV  = {"YH6": 224.0, "YH9": 243.0}   # Troyan2021 (=Ref[11]) / Kong2021 (=Ref[12] 标题)

NB3GE_R = 0.6260   # tier0_eta_ceiling_32_elements_2026-09-25.md 的最紧点（Tc_exp 类）


def ceiling_K(lam, nu2_K):
    """kTc <= C_inf * sqrt(M2), M2 = lam * nu2^2  =>  Tc_ceiling[K] = C_inf*sqrt(lam)*nu2[K]."""
    return C_INF * math.sqrt(lam) * nu2_K


def main():
    print("=" * 78)
    print("(A) ceiling 与 r —— 两套 Tc_exp 并列，方向已标")
    print("=" * 78)
    print(f"{'row':13s} {'lam':>6s} {'nu2(K)':>7s} {'ceil(K)':>8s} "
          f"{'r@220/240':>10s} {'r@224/243':>10s}")
    rows = []
    for sysname, meth, lam, nu2, tc, tcr, src in TABLE_I:
        c = ceiling_K(lam, nu2)
        r_fo = TCEXP_FAILOPEN[sysname] / c
        r_cons = TCEXP_CONSERV[sysname] / c
        rows.append((sysname, meth, lam, nu2, c, tc, tcr, r_fo, r_cons, src))
        print(f"{sysname+' '+meth:13s} {lam:6.3f} {nu2:7.0f} {c:8.2f} "
              f"{r_fo:10.4f} {r_cons:10.4f}")

    max_fo = max(r[7] for r in rows)
    max_cons = max(r[8] for r in rows)
    print(f"\n  max r_exp (fail-open 220/240) = {max_fo:.4f}")
    print(f"  max r_exp (保守    224/243) = {max_cons:.4f}   <-- 承重值")
    print("  方向：r 对 Tc_exp 单调增 => 取小的 Tc_exp 使「r<1」变易 => fail-open。")
    assert max_cons > max_fo, "保守值必须更紧，否则方向搞反了"
    print(f"  vs Nb3Ge {NB3GE_R:.4f} : "
          f"{'仍低（Nb3Ge 仍是最紧点）' if max_cons < NB3GE_R else '★超过，需重排'}")

    print()
    print("=" * 78)
    print("(B) 两类检验必须分开计数")
    print("=" * 78)
    # r_calc <= 1 是定理（同一 alpha^2F、mu*>=0 解 Eliashberg）=> 自由度 0，是门不是证据
    rc = []
    for sysname, meth, lam, nu2, c, tc, tcr, _, _, src in rows:
        rc.append(tc / c)
        if tcr is not None:
            rc.append(tcr / c)
    n_viol_calc = sum(1 for x in rc if x >= 1.0)
    n_viol_exp = sum(1 for r in rows if r[8] >= 1.0)
    print(f"  定理自检门 r_calc : {len(rc)} 行, 违反 {n_viol_calc}, 最劣 {max(rc):.4f}")
    print( "                      (r_calc<=1 由天花板定理保证 => 0 自由度，通过是门不是证据)")
    print(f"  外部效度   r_exp  : {len(rows)} 行 / 2 个体系, 违反 {n_viol_exp}, 最紧 {max_cons:.4f}")
    print( "                      (YH6 三行共享同一 Tc_exp => 独立自由度是体系数不是行数)")
    assert n_viol_calc == 0 and n_viol_exp == 0

    print()
    print("=" * 78)
    print("(C) 副发现「ceiling 对声子方法免疫，跨度<=2.6%」—— 按出处分组重测")
    print("=" * 78)
    M2 = {f"{s} {m}": lam * nu2 ** 2 for s, m, lam, nu2, *_ in TABLE_I}
    C = {f"{s} {m}": ceiling_K(lam, nu2) for s, m, lam, nu2, *_ in TABLE_I}

    def cmp(a, b, label):
        dm = abs(M2[b] / M2[a] - 1) * 100
        dc = abs(C[b] / C[a] - 1) * 100
        print(f"  {label:46s} dM2={dm:6.2f}%  dCeiling={dc:6.2f}%")
        return dc

    d1 = cmp("YH6 DFPT", "YH6 SSCHA", "YH6 DFPT vs SSCHA  [同一篇 Ref[11] 两条谱]")
    d2 = cmp("YH6 DFPT", "YH6 SPIA",  "YH6 DFPT vs SPIA   [Ref[11] vs 本文自算]")
    d3 = cmp("YH9 DFPT", "YH9 SPIA",  "YH9 DFPT vs SPIA   [Ref[20] vs 本文自算, 260 vs 255 GPa]")
    y6 = [C["YH6 DFPT"], C["YH6 SSCHA"], C["YH6 SPIA"]]
    print(f"\n  YH6 三方法 spread = {100*(max(y6)/min(y6)-1):.2f}%  <- 交付方引用的 2.6%")
    print(f"  但同一对方法在 YH9 上 = {d3:.1f}%，是 YH6 同对({d2:.1f}%)的 {d3/d2:.1f} 倍")
    print( "  => 跨度随「两端是否共享出处」单调变化：%.2f%% < %.2f%% < %.2f%%"
           % (d1, d2, d3))
    print( "  => 「对声子方法免疫」在 YH9 上未被复现；且 YH9 那一对带 260/255 GPa 压强混杂，")
    print( "     10.5%% 里方法与压强不可分 => 只能说「未复现」，不能说「已证伪」。")
    assert d3 > 2.6, "若此断言失败说明我读错了表"

    print()
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
