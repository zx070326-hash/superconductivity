"""rev3: compare the 4 server-recomputed cells (committed cell_job, kna=16) with the x4_run2.log lines;
and recount the X6 coverage 64/60/46 from the logged wbar_c table + the interval rule."""
import json, numpy as np
rows = [json.loads(l) for l in open("rev3_srv_cells.log", encoding="utf-8") if l.startswith("{")]
log = {(20.000, 20.115): (0.071, 0.454, {"0.0": 1.203, "0.1": 0.937, "0.13": 0.857, "0.26": 0.510, "0.33": 0.323}),
       (34.906, 35.107): (0.097, 0.515, {"0.0": 1.149, "0.1": 0.882, "0.13": 0.802, "0.26": 0.455, "0.33": 0.269}),
       (59.883, 60.000): (0.044, 0.504, {"0.0": 1.161, "0.1": 0.895, "0.13": 0.815, "0.26": 0.468, "0.33": 0.281})}
for r in rows:
    key = (round(r["a"], 3), round(r["b"], 3))
    print(f"cell [{r['a']:.4f},{r['b']:.4f}] Na={r['Na']} t={r['t']:.0f}s: cA={r['cA']:.3e} cB={r['cB']:.3e} Lw={r['Lw']:.4f} Lb={r['Lb']:.4f} slack={r['slack']:.4f} al={r['al']:.5f}  al-slack(=B+D+c)={r['al']-r['slack']:.5f}")
    for mt, (m, ks, case) in r["floors"].items():
        g = m - r["al"]
        ref = log.get(key, (None, None, {}))[2].get(mt)
        print(f"     mt={mt:<4} mu2={m:.5f} ks={ks} case={case} g={g:+.4f}" + (f"  log {ref:+.3f} {'match' if abs(g-ref)<6e-4 else 'MISMATCH'}" if ref is not None else ""))
    if key in log:
        print(f"     slack log {log[key][0]} vs {r['slack']:.3f}; al log {log[key][1]} vs {r['al']:.3f}")
# X6 recount from the logged table
WMIN = {0.0: 0.2, 0.1: 0.5391, 0.13: 0.6427, 0.26: 1.1686, 0.33: 1.5666}
MUT = (0.0, 0.1, 0.13, 0.26, 0.33)
tab10 = {0.4: (6.80, 12.90, 25.50, 39.40), 0.5: (4.04, 6.40, 10.20, 13.50), 0.6: (2.86, 4.04, 5.80, 7.10), 0.7: (2.22, 2.94, 3.88, 4.56),
         0.8: (1.84, 2.34, 2.94, 3.36), 1.0: (1.40, 1.68, 2.00, 2.20), 1.5: (0.96, 1.08, 1.22, 1.30), 2.0: (0.76, 0.84, 0.94, 0.98)}
tab3 = {0.4: (6.80, 13.10, 28.90, 51.20), 0.5: (4.04, 6.40, 10.90, 15.60), 0.6: (2.86, 4.04, 6.00, 7.80), 0.7: (2.22, 2.94, 4.00, 4.88),
        0.8: (1.84, 2.32, 2.98, 3.50), 1.0: (1.40, 1.68, 2.02, 2.26), 1.5: (0.96, 1.06, 1.20, 1.30), 2.0: (0.76, 0.84, 0.92, 0.98)}
c60 = c20 = c5 = 0
for tab in (tab10, tab3):
    for lam, ws in tab.items():
        for mus, wc in zip((0.0, 0.05, 0.10, 0.13), ws):
            mt = mus / lam
            key = min(m for m in MUT if m >= mt - 1e-12)
            c60 += WMIN[key] <= wc <= 60
            if key <= 0.26: c20 += WMIN[key] <= wc <= 20; c5 += WMIN[key] <= wc <= 5
print(f"X6 recount from logged wbar_c: {c60}/64 with upper 60, {c20} with upper 20, {c5} with upper 5  (log: 64/60/46)")
