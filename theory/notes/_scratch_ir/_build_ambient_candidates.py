# -*- coding: utf-8 -*-
"""把 48 个常压氢化物候选落成正式数据文件（此前只存在于探针脚本里）。
数据源：Cerqueira, Fang, Errea, Sanna, Marques, arXiv:2403.13496 TABLE I + II。
逐行由主会话从库内 PDF 的 pdftotext 全文转录。
单位 K 由脚本换算到 meV，不手敲。
"""
import csv, io, math, importlib.util, contextlib, os
K_PER_MEV = 11.604518
C_INF = 0.1827262477                      # (1/2pi)*sqrt(g(2))  [analytic]

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("h48", os.path.join(HERE, "_chk_ceiling_48_hydrides.py"))
m = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(m)

COLS = ["material","family","Tc_AD_K","omega_log_K","omega_2_K","lambda",
        "omega_log_meV","omega_2_meV","omega2_over_omegalog",
        "ceiling_K","filling_factor_f","column_assignment","notes"]
rows = []
for block, fam, ca in ((m.CLEAN[:16], "SM2-TM-H6 / Table I", "confirmed"),
                       (m.CLEAN[16:31], "Double perovskites", "confirmed"),
                       (m.CLEAN[31:], "Table II (misc.)", "confirmed"),
                       (m.PEROV, "Perovskites ABH3", "UNCERTAIN")):
    for nm, tc, wlog, w2, lam in block:
        ceil = C_INF * math.sqrt(lam) * w2
        rows.append({
            "material": nm, "family": fam, "Tc_AD_K": tc,
            "omega_log_K": wlog, "omega_2_K": w2, "lambda": lam,
            "omega_log_meV": round(wlog / K_PER_MEV, 2),
            "omega_2_meV": round(w2 / K_PER_MEV, 2),
            "omega2_over_omegalog": round(w2 / wlog, 3),
            "ceiling_K": round(ceil, 1),
            "filling_factor_f": round(tc / ceil, 4),
            "column_assignment": ca,
            "notes": ("pdftotext 把 (omega_2, lambda) 列拆到了另一行，列指派由主会话按顺序+正文旁证复原，"
                      "**不确定**；该块另有最坏配对上界 f<=0.9629" if ca == "UNCERTAIN" else ""),
        })

p = "D:/physics/superconductivity/theory/data/ambient_hydride_candidates.csv"
with io.open(p, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS); w.writeheader(); w.writerows(rows)

fs = [r["filling_factor_f"] for r in rows]
print(f"wrote {len(rows)} rows -> {p}")
print(f"f: min {min(fs):.4f}  max {max(fs):.4f}  mean {sum(fs)/len(fs):.4f}   violations(f>1): {sum(1 for x in fs if x>1)}")
den_min = min(C_INF*math.sqrt(l)*w for _,_,_,w,l in m.PEROV)
print(f"perovskite 最坏配对上界 f_max <= {max(t/den_min for _,t,_,_,_ in m.PEROV):.4f}")
