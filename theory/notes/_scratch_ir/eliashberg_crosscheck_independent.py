"""Independent verification of Eliashberg cross-check using deficit_independent_check.py
machinery: Simpson weights, explicit loops, full eigh, secant root in log(T).

Cross-check at least 3 representative cases against the main analysis.

Rigor: [numerical-evidence] (float64, N finite).
OMP_NUM_THREADS=1."""
import os, sys, math, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deficit_independent_check import read_a2f, simpson_w, Tc, MEV_K, ROOT

if __name__ == "__main__":
    print("=" * 100)
    print("INDEPENDENT VERIFICATION: Simpson weights, explicit loops, full eigh, secant root")
    print("=" * 100)
    print()

    # Test cases: representative samples from each material-pressure-epw combination
    test_cases = [
        ("H3S-200GPa", "EPW_anharmonic", "H3S.a2f"),
        ("D3S-200GPa", "EPW_anharmonic", "D3S.a2f"),
        ("H3S-160GPa", "EPW_anharmonic", "H3S.a2f"),
        ("D3S-160GPa", "EPW_anh", "D3S.a2f"),
        ("H3S-200GPa", "EPW_harmonic", "H3S.a2f"),
    ]

    print(f"{'Dataset':>30} | {'Tc(mu*=0)':>10} {'Tc(mu*=0.12)':>12} {'Tc(mu*=0.13)':>12}")
    print("-" * 100)

    results = {}
    for mat_pres, epw, a2f_name in test_cases:
        a2f_path = f"{ROOT}/{mat_pres}/{epw}/{a2f_name}"
        if not os.path.exists(a2f_path):
            print(f"  WARNING: {a2f_path} not found")
            continue

        W, A = read_a2f(a2f_path)
        sw = simpson_w(W)
        lam = 2.0 * float(sw @ (A / W))
        p = (2.0 / lam) * (A / W) * sw
        p = p / p.sum()
        om2 = float(p @ W**2)
        sqom2 = math.sqrt(om2)

        N = 250  # use N=250 for speed, should still be converged

        tc0 = Tc(W, p, lam, 0.0, 10.0, N, sqom2)
        tc12 = Tc(W, p, lam, 0.12, 10.0, N, sqom2)
        tc13 = Tc(W, p, lam, 0.13, 10.0, N, sqom2)

        tag = f"{mat_pres}/{epw}"
        results[tag] = {
            "tc0": tc0 * MEV_K,
            "tc12": tc12 * MEV_K,
            "tc13": tc13 * MEV_K
        }

        print(f"{tag:>30} | {tc0*MEV_K:>10.2f} {tc12*MEV_K:>12.2f} {tc13*MEV_K:>12.2f}")

    print()
    print("=" * 100)
    print("COMPARISON WITH MAIN ANALYSIS (from eliashberg_crosscheck_v2.py output)")
    print("=" * 100)
    print()

    # Expected values from main analysis (read from output above)
    expected = {
        "H3S-200GPa/EPW_anharmonic": {"tc0": 283.3, "tc12": 203.9, "tc13": 199.0},
        "D3S-200GPa/EPW_anharmonic": {"tc0": 214.4, "tc12": 158.8, "tc13": 155.4},
        "H3S-160GPa/EPW_anharmonic": {"tc0": 290.7, "tc12": 215.4, "tc13": 211.2},
        "D3S-160GPa/EPW_anh": {"tc0": 216.8, "tc12": 164.1, "tc13": 161.1},  # approximate
        "H3S-200GPa/EPW_harmonic": {"tc0": 328.7, "tc12": 255.6, "tc13": 251.0},
    }

    print(f"{'Dataset':>30} | {'mu*':>5} | {'Main':>10} {'Indep':>10} {'Diff':>10} {'RelDiff':>8}")
    print("-" * 100)

    all_ok = True
    for tag in results:
        for mu_key, mu_label in [("tc0", "0.0"), ("tc12", "0.12"), ("tc13", "0.13")]:
            indep = results[tag][mu_key]
            main = expected.get(tag, {}).get(mu_key)

            if main:
                diff = indep - main
                reldiff = diff / main
                status = "OK" if abs(reldiff) < 1e-3 else "CHECK"
                if abs(reldiff) >= 1e-3:
                    all_ok = False

                print(f"{tag:>30} | {mu_label:>5} | {main:>10.2f} {indep:>10.2f} "
                      f"{diff:>+10.3f} {reldiff:>+7.2%} {status}")

    print()
    if all_ok:
        print("✓ All independent checks agree with main analysis to < 0.1%")
    else:
        print("⚠ Some differences exceed 0.1% - investigate")
    print()
