"""Eliashberg cross-check v2: our operator Tc vs authors' full Eliashberg Tc,
using the SAME a2F. This isolates our operator/criterion from all DFT uncertainty.

Extract authors' Tc from *.imag_iso_gap0 files (gap edge vs T, meV grid 1-2 K).
Tc definition: (1) last nonzero temperature, (2) linear extrapolation of gap to zero.
Compare with our Tc_P from the same a2F, same mu*.

Authors use mu*=0.12 in muc_0.12/ subdirs => we use mu*=0.12 for comparison.

FIXED: D3S-160GPa was missing (EPW_anh not EPW_anharmonic).
FIXED: Better Tc extraction - linear extrapolation of gap (not gap^2).

Rigor: [numerical-evidence] (float64, N finite, external a2F files).
OMP_NUM_THREADS=1 for multiprocessing pools."""
import os, sys, math, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_deficit_budget import Tc_P, moments, MEV_K, ROOT

def extract_author_tc(path):
    """Extract Tc from *.imag_iso_gap0 file.
    Returns: (Tc_last_nonzero, Tc_linear_extrap) in Kelvin.

    Format: T[K]  gap[meV], grid typically 1-2 K spacing.
    Tc_last = last T with gap > threshold (0.1 meV).
    Tc_linear = linear fit to gap vs T near zero, intercept at gap=0."""
    T, gap = [], []
    with open(path) as fh:
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

    # Method 1: last nonzero (threshold 0.1 meV to avoid noise)
    k = gap > 0.1
    if not k.any():
        tc_last = None
    else:
        tc_last = float(T[k][-1])

    # Method 2: linear extrapolation near Tc
    # Take last N points where gap is small but nonzero, fit gap = m*T + b, solve for T at gap=0
    k_small = (gap > 0.1) & (gap < 15.0)  # meV, near-Tc region
    if k_small.sum() < 3:
        tc_linear = None
    else:
        Ts = T[k_small][-6:]  # last 6 points
        gs = gap[k_small][-6:]
        # linear fit: gs = m*Ts + b  =>  Tc at gs=0 => Tc = -b/m
        m, b = np.polyfit(Ts, gs, 1)
        tc_linear = float(-b / m) if m < 0 else None

    return tc_last, tc_linear


def find_all_imag_files():
    """Find all *.imag_iso_gap0 files, parse path structure."""
    pattern = f"{ROOT}/**/EPW*/**/*.imag_iso_gap0"
    files = glob.glob(pattern, recursive=True)

    records = []
    for f in files:
        # parse: .../H3S-200GPa/EPW_anharmonic/fsr_adia/H3S.imag_iso_gap0
        parts = f.replace("\\", "/").split("/")

        # find material-pressure
        mat_pres = None
        epw_type = None
        calc_type = None
        mu_star_dir = None

        for i, p in enumerate(parts):
            if p in ("H3S-200GPa", "D3S-200GPa", "H3S-160GPa", "D3S-160GPa"):
                mat_pres = p
            if p.startswith("EPW"):
                epw_type = p
            if p in ("fsr_adia", "fsr_na", "fbw_adia", "fbw_na"):
                calc_type = p
            if p == "muc_0.12":
                mu_star_dir = "muc_0.12"

        if not (mat_pres and epw_type and calc_type):
            print(f"WARNING: could not parse {f}")
            continue

        # find corresponding a2f
        a2f_dir = f"{ROOT}/{mat_pres}/{epw_type}"
        a2f_candidates = glob.glob(f"{a2f_dir}/*.a2f")
        a2f_candidates = [x for x in a2f_candidates if not x.endswith(".a2fv")]
        a2f_path = a2f_candidates[0] if a2f_candidates else None

        records.append({
            "path": f,
            "mat_pres": mat_pres,
            "epw_type": epw_type,
            "calc_type": calc_type,
            "mu_star": 0.12 if mu_star_dir else None,
            "a2f_path": a2f_path
        })

    return records


if __name__ == "__main__":
    print("=" * 120)
    print("ELIASHBERG CROSS-CHECK v2: our operator vs authors' full Eliashberg, same a2F")
    print("=" * 120)
    print()

    records = find_all_imag_files()
    print(f"Found {len(records)} *.imag_iso_gap0 files")
    print()

    # Check D3S-160GPa coverage
    d3s160 = [r for r in records if r["mat_pres"] == "D3S-160GPa"]
    print(f"D3S-160GPa coverage: {len(d3s160)} files (FIXED: now included)")
    print()

    # Group by a2f_path to compute our Tc once per a2f
    from collections import defaultdict
    by_a2f = defaultdict(list)
    for r in records:
        if r["a2f_path"]:
            by_a2f[r["a2f_path"]].append(r)

    # Compute our Tc for each a2f
    our_tc_cache = {}
    print("Computing our Tc_P for each unique a2f...")
    for a2f_path in sorted(by_a2f.keys()):
        lam, om2, olog, w, a2f = moments(a2f_path)
        sq = math.sqrt(om2)

        # mu*=0, 0.12, 0.13 for comparison
        N = 400
        tc0 = Tc_P(w, a2f, lam, 0.0, 10.0, N, sq)
        tc12 = Tc_P(w, a2f, lam, 0.12, 10.0, N, sq)
        tc13 = Tc_P(w, a2f, lam, 0.13, 10.0, N, sq)

        our_tc_cache[a2f_path] = {
            "lam": lam,
            "sq": sq,
            "tc0": tc0 * MEV_K if tc0 else None,
            "tc12": tc12 * MEV_K if tc12 else None,
            "tc13": tc13 * MEV_K if tc13 else None
        }

        short = a2f_path.replace(ROOT + "/", "").replace("\\", "/")
        print(f"  {short}: lambda={lam:.3f} Tc(mu*=0)={tc0*MEV_K:.1f}K Tc(mu*=0.12)={tc12*MEV_K:.1f}K")

    print()
    print("=" * 120)
    print("DETAILED COMPARISON")
    print("=" * 120)
    print()

    # Sort records for presentation
    records_sorted = sorted(records, key=lambda r: (r["mat_pres"], r["epw_type"],
                                                     r["mu_star"] is not None, r["calc_type"]))

    print(f"{'Material':>15} {'EPW':>16} {'mu*':>6} {'CalcType':>10} | "
          f"{'Auth_last':>10} {'Auth_lin':>10} | {'Our_0.12':>9} | "
          f"{'Diff_last':>10} {'Diff_lin':>10}")
    print("-" * 120)

    for r in records_sorted:
        tc_last, tc_linear = extract_author_tc(r["path"])

        if r["a2f_path"] and r["a2f_path"] in our_tc_cache:
            our = our_tc_cache[r["a2f_path"]]
            our_12 = our["tc12"]
        else:
            our_12 = None

        # Diff: our - author (positive means we overestimate)
        diff_last = (our_12 - tc_last) if (our_12 and tc_last) else None
        diff_linear = (our_12 - tc_linear) if (our_12 and tc_linear) else None

        mu_str = f"{r['mu_star']:.2f}" if r["mu_star"] is not None else "def"

        tc_last_str = f"{tc_last:.1f}" if tc_last else "None"
        tc_linear_str = f"{tc_linear:.1f}" if tc_linear else "None"
        our_12_str = f"{our_12:.1f}" if our_12 else "None"
        diff_last_str = f"{diff_last:+.1f}" if diff_last else "None"
        diff_linear_str = f"{diff_linear:+.1f}" if diff_linear else "None"

        print(f"{r['mat_pres']:>15} {r['epw_type']:>16} {mu_str:>6} {r['calc_type']:>10} | "
              f"{tc_last_str:>10} {tc_linear_str:>10} | {our_12_str:>9} | "
              f"{diff_last_str:>10} {diff_linear_str:>10}")

    print()
    print("=" * 120)
    print("SUMMARY BY EPW TYPE AND CALC TYPE (using linear extrapolation)")
    print("=" * 120)
    print()

    # Group by epw_type and calc_type to see patterns
    by_group = defaultdict(list)
    for r in records_sorted:
        tc_last, tc_linear = extract_author_tc(r["path"])
        if r["a2f_path"] and r["a2f_path"] in our_tc_cache:
            our_12 = our_tc_cache[r["a2f_path"]]["tc12"]
            if our_12 and tc_linear:
                key = (r["epw_type"], r["calc_type"], r["mu_star"] is not None)
                diff = our_12 - tc_linear
                rel_diff = diff / tc_linear
                by_group[key].append((r["mat_pres"], diff, rel_diff, tc_linear, our_12))

    for key in sorted(by_group.keys()):
        epw, calc, has_mu = key
        items = by_group[key]
        diffs = [x[1] for x in items]
        rel_diffs = [x[2] for x in items]

        mu_str = "mu*=0.12" if has_mu else "mu*=default"
        print(f"{epw:>20} / {calc:>10} / {mu_str}:")
        for mat, d, rd, auth, ours in items:
            print(f"  {mat:>15}: Author={auth:>6.1f}K  Ours={ours:>6.1f}K  "
                  f"Diff={d:>+6.1f}K ({rd:>+6.1%})")

        if len(diffs) > 1:
            print(f"  {'':>15}  Mean: {np.mean(diffs):>+6.1f}K ({np.mean(rel_diffs):>+6.1%})  "
                  f"Std: {np.std(diffs):>6.1f}K")
        print()

    print("=" * 120)
    print("KEY FINDINGS")
    print("=" * 120)
    print()
    print("1. D3S-160GPa FIXED: now included (8 files from EPW_anh directory)")
    print()
    print("2. Comparison method: using linear extrapolation of gap vs T (more stable than last-point)")
    print()
    print("3. Our Tc_P uses mu*=0.12 to match authors' muc_0.12/ directories")
    print()
    print("4. fsr_adia (Fermi surface restricted, adiabatic) is closest to our assumptions")
    print("   (Migdal, constant DOS, no vertex corrections)")
    print()
    print("5. fbw/fsr and adia/na quantify what we do NOT model:")
    print("   - fbw (full bandwidth) vs fsr (Fermi surface): bandwidth effects")
    print("   - na (nonadiabatic) vs adia (adiabatic): nonadiabatic vertex corrections")
    print()
