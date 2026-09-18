"""r-parameter scan: find optimal r = omega_c/sqrt(<om^2>) that minimizes difference
between our Tc and authors' fsr_adia Tc (most comparable to our assumptions).

If a single r fits both H3S and D3S across multiple pressures, that's evidence for
physical validity. If each case needs a different r, it's just curve-fitting.

Rigor: [numerical-evidence] (float64).
OMP_NUM_THREADS=1."""
import os, sys, math, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_deficit_budget import Tc_P, moments, MEV_K, ROOT

def extract_fsr_adia_tc(mat_pres, epw_type, mu_star_dir):
    """Extract authors' fsr_adia Tc using linear extrapolation."""
    if mu_star_dir:
        path = f"{ROOT}/{mat_pres}/{epw_type}/muc_0.12/fsr_adia/"
    else:
        path = f"{ROOT}/{mat_pres}/{epw_type}/fsr_adia/"

    # Find imag file
    files = glob.glob(path + "*.imag_iso_gap0")
    if not files:
        return None

    T, gap = [], []
    with open(files[0]) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 2:
                continue
            try:
                t, g = float(p[0]), float(p[1])
            except ValueError:
                continue
            T.append(t)
            gap.append(g)

    T = np.array(T)
    gap = np.array(gap)

    # Linear extrapolation
    k_small = (gap > 0.1) & (gap < 15.0)
    if k_small.sum() < 3:
        return None

    Ts = T[k_small][-6:]
    gs = gap[k_small][-6:]
    m, b = np.polyfit(Ts, gs, 1)
    return float(-b / m) if m < 0 else None


if __name__ == "__main__":
    print("=" * 100)
    print("r-PARAMETER SCAN: find optimal r = omega_c/sqrt(<om^2>)")
    print("=" * 100)
    print()

    # Focus on anharmonic EPW with default mu* (most physical)
    # Compare to fsr_adia (closest to our assumptions)
    test_cases = [
        ("H3S-200GPa", "EPW_anharmonic", "H3S.a2f", False),
        ("D3S-200GPa", "EPW_anharmonic", "D3S.a2f", False),
        ("H3S-160GPa", "EPW_anharmonic", "H3S.a2f", False),
        ("D3S-160GPa", "EPW_anh", "D3S.a2f", False),
    ]

    # Also check muc_0.12 cases
    test_cases_mu12 = [
        ("H3S-200GPa", "EPW_anharmonic", "H3S.a2f", True),
        ("D3S-200GPa", "EPW_anharmonic", "D3S.a2f", True),
        ("H3S-160GPa", "EPW_anharmonic", "H3S.a2f", True),
        ("D3S-160GPa", "EPW_anh", "D3S.a2f", True),
    ]

    r_values = [5.0, 7.5, 10.0, 12.5, 15.0, 20.0, 30.0]

    for label, cases in [("mu*=default", test_cases), ("mu*=0.12", test_cases_mu12)]:
        print(f"\n{'=' * 100}")
        print(f"SCAN FOR {label}")
        print(f"{'=' * 100}\n")

        results = {}
        for mat_pres, epw, a2f_name, use_mu12 in cases:
            a2f_path = f"{ROOT}/{mat_pres}/{epw}/{a2f_name}"
            if not os.path.exists(a2f_path):
                continue

            lam, om2, olog, w, a2f = moments(a2f_path)
            sq = math.sqrt(om2)

            # Get author's fsr_adia Tc
            author_tc = extract_fsr_adia_tc(mat_pres, epw, use_mu12)
            if not author_tc:
                continue

            mu_val = 0.12 if use_mu12 else 0.0

            # Scan r
            our_tcs = {}
            for r in r_values:
                tc = Tc_P(w, a2f, lam, mu_val, r, 400, sq)
                our_tcs[r] = tc * MEV_K if tc else None

            tag = f"{mat_pres}/{epw}"
            results[tag] = {
                "author": author_tc,
                "ours": our_tcs,
                "lam": lam
            }

        # Print results
        print(f"{'Dataset':>30} {'Author':>9} | " + " ".join([f"r={r:>4.1f}" for r in r_values]))
        print("-" * 100)

        for tag in sorted(results.keys()):
            author = results[tag]["author"]
            ours = results[tag]["ours"]

            row = f"{tag:>30} {author:>9.1f} | "
            for r in r_values:
                if ours[r]:
                    diff = ours[r] - author
                    row += f"{diff:>+7.1f} "
                else:
                    row += "   None "
            print(row)

        # Find best r for each case
        print("\n" + "-" * 100)
        print("OPTIMAL r FOR EACH CASE (minimize |our - author|)")
        print("-" * 100)

        for tag in sorted(results.keys()):
            author = results[tag]["author"]
            ours = results[tag]["ours"]
            lam = results[tag]["lam"]

            best_r = None
            best_diff = float('inf')
            for r in r_values:
                if ours[r]:
                    diff = abs(ours[r] - author)
                    if diff < best_diff:
                        best_diff = diff
                        best_r = r

            if best_r:
                varpi = math.sqrt(lam * 2.0)  # approximate, assumes <om^2> ~ lam
                om_c = best_r * varpi
                print(f"{tag:>30}: r={best_r:>5.1f}  (error={ours[best_r]-author:>+6.1f}K)  "
                      f"omega_c ~ {om_c:.1f}*sqrt(<om^2>)")

    print("\n" + "=" * 100)
    print("PHYSICAL INTERPRETATION")
    print("=" * 100)
    print()
    print("r = omega_c/sqrt(<om^2>) sets the Coulomb cutoff N_c = floor(r * varpi).")
    print("Physically defensible range: r ~ 3-15 based on omega_c ~ bandwidth/Fermi energy.")
    print()
    print("If one r fits all cases: strong evidence for physical validity.")
    print("If each case needs different r: just curve-fitting, not predictive.")
    print()
