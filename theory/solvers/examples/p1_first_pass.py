"""
theory/solvers/examples/p1_first_pass.py

P1 (theory/derivations/P1_isotope_master_formula.md) 首轮数值检验：
对 lambda=2, omega_log=100 meV 的模型谱（Einstein 单模 + 双峰两种实现），
扫 mu_star = 0.10/0.13/0.20/0.25/0.30，输出数值 alpha（mu_star_fixed=False）
与 P1 式(3.2)（McMillan 闭式）的对照表；另外输出数值 d ln Tc/d ln lambda
与 McMillan 解析值的对照。

运行：
    python theory/solvers/examples/p1_first_pass.py

注意：本脚本只打印结果，不做 assert——数值判定在 tests/test_eliashberg_iso.py 里。
本脚本的输出被原样贴进 README.md 的"首轮结果"一节，运行本脚本的人和读 README 的人
看到的应该是同一张表。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from eliashberg_iso import Spectrum, Solver  # noqa: E402

W_CUT_RATIO = 5.0  # w_cut = 5 * omega_log，与 tests/test_eliashberg_iso.py 保持一致的约定
MU_STAR_LIST = [0.10, 0.13, 0.20, 0.25, 0.30]


def eq32_mcmillan_alpha(lam: float, mu_star: float) -> float:
    """P1 sec 3.1 eq (3.2): McMillan-based closed-form isotope coefficient.

    Deliberately re-derived here (not imported from eliashberg_iso.py) so this
    script compares the numeric solver against an independent transcription of
    the target formula, not against itself.
    """
    d = lam - mu_star * (1 + 0.62 * lam)
    return 0.5 * (1 - 1.04 * (1 + lam) * (1 + 0.62 * lam) * mu_star ** 2 / d ** 2)


def mcmillan_dlnTc_dlnlam(lam: float, mu_star: float, h: float = 1e-6) -> float:
    """McMillan Phi = exp[-1.04(1+lam)/D] 的 d ln(Phi)/d ln(lam)，中心差分对解析式求导。"""
    def ln_phi(l):
        d = l - mu_star * (1 + 0.62 * l)
        return -1.04 * (1 + l) / d
    return (ln_phi(lam + h) - ln_phi(lam - h)) / (2 * h) * lam


def alpha_table_for_spectrum(sp: Spectrum, label: str) -> None:
    print(f"\n--- {label}  (lambda={sp.lambda_total():.4f}, omega_log={sp.omega_log():.4f} meV) ---")
    print(f"{'mu*':>6} | {'alpha_num':>10} | {'alpha_eq32':>10} | {'diff':>8}")
    print("-" * 44)
    lam = sp.lambda_total()
    for mu_star in MU_STAR_LIST:
        solver = Solver(sp, mu_star=mu_star, w_cut=W_CUT_RATIO * sp.omega_log())
        alpha_num, _ = solver.isotope_coefficient({"full": sp}, which="all", mu_star_fixed=False)
        alpha_ref = eq32_mcmillan_alpha(lam, mu_star)
        diff = alpha_num - alpha_ref
        flag = "" if abs(diff) < 0.01 else "  (out of <0.01 tolerance)"
        print(f"{mu_star:6.2f} | {alpha_num:10.4f} | {alpha_ref:10.4f} | {diff:+8.4f}{flag}")


def lambda_sensitivity_table() -> None:
    print("\n--- d ln Tc / d ln lambda: numeric vs McMillan analytic (mu*=0.13) ---")
    print(f"{'lambda':>6} | {'s_num':>8} | {'s_mcmillan':>10} | {'diff':>8}")
    print("-" * 42)
    mu_star = 0.13
    for lam in (1.0, 2.0, 3.0):
        sp = Spectrum.einstein(lam=lam, w0=100.0)
        solver = Solver(sp, mu_star=mu_star, w_cut=W_CUT_RATIO * sp.omega_log())
        s_num = solver.log_sensitivity("lambda")
        s_ref = mcmillan_dlnTc_dlnlam(lam, mu_star)
        print(f"{lam:6.1f} | {s_num:8.4f} | {s_ref:10.4f} | {s_num - s_ref:+8.4f}")


def main() -> None:
    sp_einstein = Spectrum.einstein(lam=2.0, w0=100.0)
    sp_two_peak = Spectrum.two_peak(1.0, 50.0, 1.0, 200.0)

    print("=" * 60)
    print("P1 first-pass check: numeric alpha vs eq (3.2), lambda=2, omega_log=100 meV")
    print("Convention: w_cut = 5 * omega_log, mu_star_fixed=False (Tolmachev RG)")
    print("=" * 60)
    alpha_table_for_spectrum(sp_einstein, "Einstein single mode (lam=2, w0=100 meV)")
    alpha_table_for_spectrum(sp_two_peak, "Two-peak (lam1=1,w1=50 + lam2=1,w2=200 meV)")

    print()
    print("=" * 60)
    print("Secondary check: d ln Tc/d ln lambda, numeric solver vs McMillan analytic")
    print("=" * 60)
    lambda_sensitivity_table()


if __name__ == "__main__":
    main()
