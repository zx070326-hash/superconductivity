"""[numerical-evidence] From-zero recomputation of the derived columns in
theory/data/ambient_hydride_candidates.csv.

Nobody has independently checked these.  Two derived columns are load-bearing
for the deficit-budget / filling-factor argument:

    ceiling_K        =?  0.1827 * sqrt(lambda) * omega_2      (Allen-Dynes
                         strong-coupling asymptote Tc -> 0.1827 sqrt(<w^2> lambda))
    filling_factor_f =?  Tc_AD / ceiling_K
    omega_log_meV    =?  omega_log_K / 11.604518
    omega_2_meV      =?  omega_2_K   / 11.604518
    omega2_over_omegalog =? omega_2_K / omega_log_K

The formula for ceiling_K is NOT stated in the CSV; it is inferred here and then
tested on every row.  If the inference is wrong the residuals will not be small.

Run: OMP_NUM_THREADS=1 python _chk_ambient_csv_recompute.py
"""
import os, csv
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

CSV = os.path.join("..", "..", "data", "ambient_hydride_candidates.csv")
K_PER_MEV = 11.604518
AD_ASYMPTOTE = 0.1827


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


if __name__ == "__main__":
    rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
    print(f"rows = {len(rows)}")
    cols = ["ceiling_K", "filling_factor_f", "omega_log_meV", "omega_2_meV",
            "omega2_over_omegalog"]
    worst = {c: (0.0, "") for c in cols}
    nbad = {c: 0 for c in cols}
    fvals, names = [], []

    for r in rows:
        lam, w2, wl = num(r["lambda"]), num(r["omega_2_K"]), num(r["omega_log_K"])
        tc = num(r["Tc_AD_K"])
        pred = {
            "ceiling_K": AD_ASYMPTOTE * np.sqrt(lam) * w2,
            "omega_log_meV": wl / K_PER_MEV,
            "omega_2_meV": w2 / K_PER_MEV,
            "omega2_over_omegalog": w2 / wl,
        }
        pred["filling_factor_f"] = tc / pred["ceiling_K"]
        for c in cols:
            got, exp = num(r[c]), pred[c]
            if np.isnan(got) or np.isnan(exp):
                continue
            rel = abs(got - exp) / max(abs(exp), 1e-12)
            if rel > worst[c][0]:
                worst[c] = (rel, r["material"])
            if rel > 2e-3:
                nbad[c] += 1
                if nbad[c] <= 3:
                    print(f"   MISMATCH {c:22s} {r['material']:12s}"
                          f" csv={got:.6g} recomputed={exp:.6g} rel={rel:.2e}")
        if not np.isnan(pred["filling_factor_f"]):
            fvals.append(pred["filling_factor_f"]); names.append(r["material"])

    print(f"\n{'column':>22} {'worst rel.diff':>15} {'at':>14} {'#>2e-3':>7}")
    for c in cols:
        print(f"{c:>22} {worst[c][0]:15.2e} {worst[c][1]:>14} {nbad[c]:7d}")

    f = np.array(fvals)
    o = np.argsort(f)
    print(f"\n== filling factor f = Tc_AD / (0.1827 sqrt(lam) omega_2), recomputed ==")
    print(f"   n={len(f)}  mean={f.mean():.4f}  sd={f.std(ddof=1):.4f}"
          f"  min={f.min():.4f}  max={f.max():.4f}  sd/mean={f.std(ddof=1)/f.mean():.3f}")
    print(f"   lowest 4 : " + ", ".join(f"{names[i]} {f[i]:.3f}" for i in o[:4]))
    print(f"   highest 4: " + ", ".join(f"{names[i]} {f[i]:.3f}" for i in o[-4:]))
