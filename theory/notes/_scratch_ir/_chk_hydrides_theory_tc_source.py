"""[numerical-evidence] Is Tc_predicted_K in hydrides_theory.csv the SOURCE's own
reported Tc, or is it Allen-Dynes recomputed by the harvester?

This matters: if it is recomputed AD, the column must not be quoted as "the
paper's Eliashberg prediction" -- most hydride papers solve full Migdal-Eliashberg
(or SCDFT), which differs from AD by 5-20%.  Rows where AD reproduces the value
to <1% are recomputations; rows that deviate carry independent information.

Uses each row's own mustar_used.  Run:
    OMP_NUM_THREADS=1 python _chk_hydrides_theory_tc_source.py
"""
import os, csv
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

CSV = os.path.join("..", "..", "data", "hydrides_theory.csv")


def tc_ad(lam, wlog_meV, w2_meV, mus):
    """AD 1975 Eq.34; omegas in meV -> Tc in K via 1 meV = 11.604518 K."""
    wlog, w2 = wlog_meV * 11.604518, w2_meV * 11.604518
    r = w2 / wlog
    core = lam - mus * (1.0 + 0.62 * lam)
    if core <= 0:
        return np.nan
    L1 = 2.46 * (1.0 + 3.8 * mus)
    L2 = 1.82 * (1.0 + 6.3 * mus) * r
    f1 = (1.0 + (lam / L1) ** 1.5) ** (1.0 / 3.0)
    f2 = 1.0 + (r - 1.0) * lam ** 2 / (lam ** 2 + L2 ** 2)
    return f1 * f2 * wlog / 1.2 * np.exp(-1.04 * (1.0 + lam) / core)


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


if __name__ == "__main__":
    rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
    print(f"{'material':>8} {'P':>5} {'treat':>11} {'lam':>6} {'w2/wlog':>8}"
          f" {'mu*':>5} {'Tc_csv':>7} {'Tc_AD':>7} {'dev':>8} {'tier':>10}")
    for x in rows:
        lam, wl, w2 = num(x["lambda"]), num(x["omega_log_meV"]), num(x["omega_2_meV"])
        mus, tc = num(x["mustar_used"]), num(x["Tc_predicted_K"])
        ad = tc_ad(lam, wl, w2, mus) if not any(
            np.isnan(v) for v in (lam, wl, w2, mus)) else np.nan
        dev = ad / tc - 1.0 if tc else np.nan
        r = w2 / wl if wl else np.nan
        d = f"{100*dev:+7.1f}%" if np.isfinite(dev) else "     n/a"
        print(f"{x['material']:>8} {x['pressure_GPa']:>5} {x['treatment']:>11}"
              f" {lam:6.3f} {r:8.3f} {mus:5.2f} {tc:7.1f}"
              f" {ad if np.isfinite(ad) else float('nan'):7.1f} {d} "
              f"{x['evidence_tier']:>10}")
