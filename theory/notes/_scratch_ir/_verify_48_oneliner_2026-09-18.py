import csv, math

C_INF = 0.1827262477890  # (1/2pi)*sqrt(g(2)), upper CERT-D endpoint, analytic-labelled constant

rows = list(csv.DictReader(open("theory/data/ambient_hydride_candidates.csv", encoding="utf-8")))
fs = []
viol = 0
ratio_om = []
for r in rows:
    tc = float(r["Tc_AD_K"]); lam = float(r["lambda"])
    w2K = float(r["omega_2_K"]); wlogK = float(r["omega_log_K"])
    ceil_K = C_INF * math.sqrt(lam) * w2K
    f = tc / ceil_K
    fs.append(f)
    if f > 1:
        viol += 1
    ratio_om.append(w2K / wlogK)

print("n =", len(rows))
print("f: min %.4f  max %.4f  mean %.4f" % (min(fs), max(fs), sum(fs) / len(fs)))
print("violations f>1:", viol)
print("omega2/omegalog: min %.3f  max %.3f  mean %.3f" % (min(ratio_om), max(ratio_om), sum(ratio_om) / len(ratio_om)))
