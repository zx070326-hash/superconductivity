"""
theory/solvers/examples/dos_scan.py

第二阶段功能 6（能量依赖 DOS）的第一版扫描图：对数型 van Hove 奇点（DOS.log_vhs）
相对费米面的位置 delta 在 [-3, 3]*omega_log 范围内扫描，输出 Tc 与同位素系数 alpha
（全缩放，mu_star_fixed=False，即 mu* 跟着 w_cut 走 Tolmachev 关系一起变）随 delta
的变化表。

这是 P1 §3.4 提到的 "s_DOS" 图的第一版数值实现（数值验证器角色，不是最终论文图）。

模型参数（团队指定）：
    Einstein 谱 lambda=1.84, omega0=93 meV, mu_star=0.16（w_cut=5*omega_log=465 meV）
    DOS.log_vhs(a=0.3, W=1000 meV)，delta/omega_log 取 [-3, 3] 之间 13 个等间距点

cut 的选择（团队未指定，本脚本的选择与理由）：
    cut=5 meV。理由：(1) 需要远小于 W=1000 meV，避免正则化尺度污染"背景能标"；
    (2) 需要与 pi*T（T~100-160K 时约 30-45 meV）同量级或更小，否则正则化本身会
    抹平最靠近费米面的几个 Matsubara 频率感受到的 vHs 形状；实测 cut=2 meV 和
    cut=5 meV 给出的 Tc 定性形状一致（数值差 <1%），说明结果对 cut 在这个量级内
    不敏感，取 5 meV 只是为了和 tests/test_eliashberg_iso.py 里用的同一个值保持
    一致，方便交叉核对。

运行：
    python theory/solvers/examples/dos_scan.py

注意：本脚本只打印结果，不做 assert；数值判定见 tests/test_eliashberg_iso.py
（"13. SPEC 必备功能 6"一节）。delta=0 处的结果初看反直觉（Tc 低于常数 DOS 基线，
而不是 SPEC 测试(c)字面预期的"更高"）——这不是 bug，是 N̄(0)≡N(ε)/N(0)|_{ε=0}≡1
这个归一化定义本身在 delta=0 处的必然推论，详见 README「能量依赖 DOS」一节与
tests 里 test_dos_vhs_pinned_exactly_at_ef_suppresses_both_tc_and_isotope_effect
的 docstring。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from eliashberg_iso import Spectrum, Solver, DOS  # noqa: E402

LAM = 1.84
W0_MEV = 93.0
MU_STAR = 0.16
W_CUT_RATIO = 5.0
DOS_A = 0.3
DOS_W_MEV = 1000.0
DOS_CUT_MEV = 5.0
DELTA_OVER_WLOG = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]


def main() -> None:
    sp = Spectrum.einstein(lam=LAM, w0=W0_MEV)
    wlog = sp.omega_log()
    w_cut = W_CUT_RATIO * wlog

    solver_const = Solver(sp, mu_star=MU_STAR, w_cut=w_cut)
    tc_const = solver_const.tc()
    alpha_const, _ = solver_const.isotope_coefficient(
        {"full": sp}, which="all", mu_star_fixed=False
    )

    print("=" * 70)
    print("DOS scan: log-type van Hove singularity vs distance from E_F")
    print(f"Einstein spectrum lambda={LAM}, omega_log={wlog:.4f} meV, "
          f"mu_star={MU_STAR}, w_cut={w_cut:.4f} meV")
    print(f"DOS.log_vhs(a={DOS_A}, W={DOS_W_MEV} meV, cut={DOS_CUT_MEV} meV)")
    print(f"Baseline (constant DOS): Tc={tc_const:.4f} K, alpha={alpha_const:.4f}")
    print("=" * 70)
    print(f"{'delta/wlog':>10} | {'delta(meV)':>10} | {'Tc(K)':>9} | "
          f"{'Tc/Tc0':>7} | {'alpha':>7}")
    print("-" * 70)

    rows = []
    for ratio in DELTA_OVER_WLOG:
        delta = ratio * wlog
        dos = DOS.log_vhs(delta=delta, a=DOS_A, W=DOS_W_MEV, cut=DOS_CUT_MEV)
        solver = Solver(sp, mu_star=MU_STAR, w_cut=w_cut, dos=dos, dos_n_eps=200)
        tc = solver.tc()
        alpha, _ = solver.isotope_coefficient({"full": sp}, which="all", mu_star_fixed=False)
        rows.append((ratio, delta, tc, alpha))
        print(f"{ratio:10.2f} | {delta:10.2f} | {tc:9.4f} | {tc / tc_const:7.4f} | {alpha:7.4f}")

    print()
    print(f"{'sanity check':<40}: alpha at delta=0 vs baseline: "
          f"{[r[3] for r in rows if r[0] == 0.0][0]:.4f} vs {alpha_const:.4f} "
          f"({'suppressed' if [r[3] for r in rows if r[0] == 0.0][0] < alpha_const else 'NOT suppressed'})")


if __name__ == "__main__":
    main()
