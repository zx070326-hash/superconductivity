"""[numerical-evidence] Is Tc_AD_K in ambient_hydride_candidates.csv just the
Allen-Dynes formula applied to the listed (lambda, omega_log, omega_2)?

The CSV has no mustar column.  If Tc_AD reproduces from AD at some single mu*,
then filling_factor_f = Tc_AD / (0.1827 sqrt(lam) omega_2) carries ZERO
information beyond (lambda, mu*, omega_2/omega_log) -- it is a deterministic
function of them, not an empirical regularity.  That changes what the spread in
f is allowed to mean.

Allen-Dynes 1975 Eq. (34):
    Tc = (f1 f2 omega_log / 1.2) exp[ -1.04(1+lam) / (lam - mu*(1+0.62 lam)) ]
    f1 = [1 + (lam/Lam1)^{3/2}]^{1/3},        Lam1 = 2.46 (1 + 3.8 mu*)
    f2 = 1 + (r-1) lam^2/(lam^2 + Lam2^2),    Lam2 = 1.82 (1 + 6.3 mu*) r,
    r  = omega_2/omega_log
Also compare the plain McMillan-form (f1=f2=1) to see which factor matters.

Run: OMP_NUM_THREADS=1 python _chk_ambient_tc_is_AD.py
"""
import os, csv
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

CSV = os.path.join("..", "..", "data", "ambient_hydride_candidates.csv")


def tc_ad(lam, wlog, w2, mus, use_f=True):
    r = w2 / wlog
    core = lam - mus * (1.0 + 0.62 * lam)
    if core <= 0:
        return 0.0
    f1 = f2 = 1.0
    if use_f:
        L1 = 2.46 * (1.0 + 3.8 * mus)
        L2 = 1.82 * (1.0 + 6.3 * mus) * r
        f1 = (1.0 + (lam / L1) ** 1.5) ** (1.0 / 3.0)
        f2 = 1.0 + (r - 1.0) * lam ** 2 / (lam ** 2 + L2 ** 2)
    return f1 * f2 * wlog / 1.2 * np.exp(-1.04 * (1.0 + lam) / core)


if __name__ == "__main__":
    rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
    lam = np.array([float(r["lambda"]) for r in rows])
    wl = np.array([float(r["omega_log_K"]) for r in rows])
    w2 = np.array([float(r["omega_2_K"]) for r in rows])
    tc = np.array([float(r["Tc_AD_K"]) for r in rows])
    name = [r["material"] for r in rows]

    print(f"{'mu*':>6} {'form':>10} {'median rel.dev':>15} {'max |rel.dev|':>14}"
          f" {'at':>12} {'#within 2%':>11}")
    best = None
    for mus in (0.08, 0.10, 0.12, 0.13, 0.15):
        for use_f, tag in ((True, "AD f1f2"), (False, "McMillan")):
            pred = np.array([tc_ad(a, b, c, mus, use_f)
                             for a, b, c in zip(lam, wl, w2)])
            dev = pred / tc - 1.0
            med, mx = np.median(dev), np.abs(dev).max()
            i = np.abs(dev).argmax()
            nw = int((np.abs(dev) < 0.02).sum())
            print(f"{mus:6.2f} {tag:>10} {med:15.4f} {mx:14.4f} {name[i]:>12}"
                  f" {nw:6d}/{len(tc)}")
            if best is None or mx < best[0]:
                best = (mx, mus, use_f, tag, dev, pred)

    mx, mus, use_f, tag, dev, pred = best
    print(f"\nbest: mu*={mus}, {tag}, max|rel.dev| = {mx:.4f}")
    o = np.argsort(-np.abs(dev))
    print("   worst 5 rows: " + ", ".join(
        f"{name[i]} csv={tc[i]:.1f} AD={pred[i]:.1f} ({100*dev[i]:+.1f}%)"
        for i in o[:5]))
