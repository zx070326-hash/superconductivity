"""[analytic identity + numerical-evidence] filling_factor_f is a closed-form
function of (lambda, r) ONLY -- it carries no material-specific information.

Since Tc_AD_K is verified (_chk_ambient_tc_is_AD.py) to be Allen-Dynes at
mu*=0.10 on all 48 rows, and ceiling = 0.1827 sqrt(lam) omega_2, the two
omega_2 factors cancel and omega_log drops out entirely:

    f = Tc_AD / (0.1827 sqrt(lam) omega_2)
      = f1(lam) f2(lam,r) / (1.2 * 0.1827 * sqrt(lam) * r)
        * exp[-1.04(1+lam)/(lam - mu*(1+0.62 lam))]                    (*)

i.e. f = F(lam, r) with NO free material input.  Test (*) against the CSV
column, then ask how much of the observed spread in f is explained by lam alone.

Run: OMP_NUM_THREADS=1 python _chk_filling_factor_is_closed_form.py
"""
import os, csv
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

CSV = os.path.join("..", "..", "data", "ambient_hydride_candidates.csv")
MUS = 0.10


def f_closed(lam, r, mus=MUS):
    """Eq. (*): filling factor from (lambda, r) alone."""
    L1 = 2.46 * (1.0 + 3.8 * mus)
    L2 = 1.82 * (1.0 + 6.3 * mus) * r
    f1 = (1.0 + (lam / L1) ** 1.5) ** (1.0 / 3.0)
    f2 = 1.0 + (r - 1.0) * lam ** 2 / (lam ** 2 + L2 ** 2)
    core = lam - mus * (1.0 + 0.62 * lam)
    return (f1 * f2 / (1.2 * 0.1827 * np.sqrt(lam) * r)
            * np.exp(-1.04 * (1.0 + lam) / core))


if __name__ == "__main__":
    rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
    lam = np.array([float(r["lambda"]) for r in rows])
    r = np.array([float(x["omega2_over_omegalog"]) for x in rows])
    fcsv = np.array([float(x["filling_factor_f"]) for x in rows])
    name = [x["material"] for x in rows]

    fcl = f_closed(lam, r)
    dev = fcl / fcsv - 1.0
    i = np.abs(dev).argmax()
    print(f"closed form (*) vs CSV column: max |rel.dev| = {np.abs(dev).max():.4f}"
          f" at {name[i]}   (#within 2%: {int((np.abs(dev)<0.02).sum())}/{len(fcsv)})")

    print(f"\nrange of inputs:  lam in [{lam.min():.2f},{lam.max():.2f}]"
          f"   r in [{r.min():.3f},{r.max():.3f}]")
    print(f"observed f:       mean={fcsv.mean():.4f} sd={fcsv.std(ddof=1):.4f}"
          f"  (sd/mean={fcsv.std(ddof=1)/fcsv.mean():.3f})")

    # how much of the spread does lambda alone explain?
    rb = float(np.median(r))
    f_lam_only = f_closed(lam, rb)
    res = fcsv - f_lam_only
    print(f"\nf at fixed r = median = {rb:.3f}, varying lam only:")
    print(f"   sd of residual f - F(lam, r_med) = {res.std(ddof=1):.4f}"
          f"   => lam alone explains "
          f"{100*(1 - res.var(ddof=1)/fcsv.var(ddof=1)):.1f}% of the variance")

    # NB: f turns out to be INCREASING in lam over this range (checked below).
    print(f"\n{'lam':>6} {'F(lam, r_med)':>14}")
    for a in (0.5, 0.8, 1.0, 1.4, 2.0, 3.0, 5.0):
        print(f"{a:6.2f} {f_closed(a, rb):14.4f}")
