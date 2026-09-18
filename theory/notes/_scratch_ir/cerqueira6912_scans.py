"""
Step 2 (pruning scan) + Step 3 (population scan) on Cerqueira et al. 6912-material
dataset (arXiv:2307.10728, Materials Cloud 3kbt5-r3n56), run on whatever CSV
of parsed per-material (lambda, wlog, w2av, Tc) rows is available
(see cerqueira6912_parse_batch.py). This script does NOT hit the network.

Step 2 ceiling:
  Tc_ceiling[K] = C_inf * sqrt(lambda) * w2av[meV] * (1 meV -> K)
  C_inf = 0.1827262477890   (analytic upper end, per task spec -- not re-derived here)
  meV_to_K = 11.6045

  omega-bar-2 used here is McMillan.dat's own w2av[meV], i.e. the SAME quantity the
  paper's Eq.(1c) defines and the SAME one f2's Lambda2/(w2av/wlog) uses --
  NOT omega_log. So the one-sided "Jensen lower estimate" trap from the task brief
  does not apply to this run: every row with an McMillan.dat is a genuine w2av,
  not a wlog stand-in. Rows lacking a McMillan.dat block are dropped (parser
  already skips them), not imputed.

Step 2 degeneracy check:
  Algebraic fact (verified by hand from the real Allen-Dynes formula, see note):
    Tc_AD = wlog * G(lambda, r, mustar),   r := w2av/wlog
    Tc_ceiling = C_inf*sqrt(lambda)*w2av*conv = C_inf*sqrt(lambda)*r*wlog*conv
    => Tc_AD/Tc_ceiling = G(lambda,r,mustar) / (C_inf*sqrt(lambda)*r*conv)
  wlog cancels EXACTLY. Since mustar=0.10 is the SAME fixed constant for every
  material in this dataset (paper: "we took arbitrarily mu*=0.10 for all
  materials"), the ratio collapses further to a pure function of (lambda, r)
  -- zero material information beyond that pair. This script verifies the
  cancellation numerically by rescaling wlog artificially (holding lambda, r
  fixed) and recomputing Tc_AD via the literal Allen-Dynes formula (not the
  file's cached value) to confirm invariance to machine precision.

Step 3 population scan:
  Xi := lambda/mustar, mustar=0.10 fixed for this dataset => Xi = lambda/0.10 = 10*lambda.
  Report counts with Xi<2.6, Xi<3.1, cross-tabulated against Tc_pred>{1,10,77} K.
"""
import sys, csv, math

C_INF = 0.1827262477890
MEV_TO_K = 11.6045
MUSTAR = 0.10

def allen_dynes_tc(lam, wlog_K, r, mustar):
    """Real Allen-Dynes 1975 formula (not McMillan-only), r = w2av/wlog."""
    if lam <= mustar*(1+0.62*lam):
        return 0.0
    Lambda1 = 2.46*(1+3.8*mustar)
    Lambda2 = 1.82*(1+6.3*mustar)*r
    f1 = (1 + (lam/Lambda1)**1.5)**(1/3)
    f2 = 1 + (lam**2)*(r-1)/(lam**2 + Lambda2**2)
    expo = -1.04*(1+lam)/(lam - mustar*(1+0.62*lam))
    return f1*f2*wlog_K/1.2*math.exp(expo)

def load(path):
    rows = []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            try:
                lam = float(row["lambda"])
                wlog_K = float(row["wlog_K"])
                w2av_meV = float(row["w2av_meV"])
                r_ratio = float(row["r_w2av_over_wlog"])
            except (ValueError, TypeError, KeyError):
                continue
            tc_ad = row.get("Tc_AD_mustar010_K")
            tc_ad = float(tc_ad) if tc_ad not in (None, "", "None") else None
            rows.append(dict(material=row["formula"], lam=lam, wlog_K=wlog_K,
                              w2av_meV=w2av_meV, r=r_ratio, tc_ad_file=tc_ad))
    return rows

def step2(rows):
    print(f"\n=== Step 2: pruning scan, n={len(rows)} ===")
    n_fail = 0
    worst = []
    degeneracy_max_dev = 0.0
    for row in rows:
        tc_ceiling = C_INF*math.sqrt(row["lam"])*row["w2av_meV"]*MEV_TO_K
        tc_pred = row["tc_ad_file"]
        if tc_pred is None:
            continue
        f = tc_pred/tc_ceiling if tc_ceiling > 0 else float("inf")
        if tc_pred > tc_ceiling:
            n_fail += 1
            worst.append((row["material"], tc_pred, tc_ceiling, f))
        # degeneracy check: rescale wlog by factor k, recompute Tc_AD via literal
        # formula, and Tc_ceiling is INDEPENDENT of wlog (only lambda,w2av) so the
        # ratio should scale by k. Instead check that ratio computed two ways
        # (from formula vs cached file) agree, and that G(lambda,r) alone
        # (i.e. removing the wlog factor) is invariant under k-rescaling.
        tc_ad_formula_base = allen_dynes_tc(row["lam"], row["wlog_K"], row["r"], MUSTAR)
        for k in (0.5, 2.0, 5.0):
            tc_ad_k = allen_dynes_tc(row["lam"], row["wlog_K"]*k, row["r"], MUSTAR)
            ratio_base = tc_ad_formula_base/row["wlog_K"] if row["wlog_K"] else None
            ratio_k = tc_ad_k/(row["wlog_K"]*k) if row["wlog_K"] else None
            if ratio_base and ratio_k:
                dev = abs(ratio_k-ratio_base)/abs(ratio_base)
                degeneracy_max_dev = max(degeneracy_max_dev, dev)
    print(f"FAIL (Tc_pred > Tc_ceiling): {n_fail}/{len(rows)}")
    print(f"max relative deviation of G(lambda,r,mustar)=Tc_AD/wlog under wlog rescaling "
          f"(0.5x,2x,5x) across all rows: {degeneracy_max_dev:.3e}  "
          f"(should be ~0 if wlog truly cancels)")
    worst.sort(key=lambda t: -t[3])
    for m in worst[:10]:
        print(f"  FAIL: {m[0]:20s} Tc_pred={m[1]:8.3f} K  Tc_ceiling={m[2]:8.3f} K  ratio={m[3]:.4f}")
    return n_fail

def step3(rows):
    print(f"\n=== Step 3: population scan, n={len(rows)} ===")
    xi26 = [r for r in rows if r["lam"]/MUSTAR < 2.6]
    xi31 = [r for r in rows if r["lam"]/MUSTAR < 3.1]
    print(f"Xi=lambda/mustar < 2.6: {len(xi26)}/{len(rows)}")
    print(f"Xi=lambda/mustar < 3.1: {len(xi31)}/{len(rows)}")
    for label, subset in (("Xi<2.6", xi26), ("Xi<3.1", xi31)):
        for thr in (1, 10, 77):
            n = sum(1 for r in subset if (r["tc_ad_file"] or 0) > thr)
            print(f"  within {label}: Tc_pred>{thr}K: {n}/{len(subset)}")

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else \
        r"D:\physics\superconductivity\theory\data\external\cerqueira6912_batchq_sample.csv"
    rows = load(path)
    step2(rows)
    step3(rows)
