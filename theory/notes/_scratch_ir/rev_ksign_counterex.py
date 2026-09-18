"""REVIEWER: hunt for a cell [a,b] where the SHIPPED test passes (k-signs of the best case >= 0)
but some other window case has a negative k AND the returned bound is violated inside the cell.
This is the concrete failure mode of checking only the best case.  Scans mutilde and a."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import rev_common as rc
from p4d_coulomb_certify import mu2_pencil_cell
found = []
tested = 0; shipped_ok_otherfail = 0
for mt in (0.05, 0.1, 0.13, 0.2, 0.26, 0.3, 0.4, 0.5, 0.7, 1.0):
    for a in np.geomspace(0.3, 5.0, 300):
        b = a * 1.05
        mu2, k1, k2, k3, case = mu2_pencil_cell(a, b, mt)
        shipped_ok = bool(k1 >= 0 and k2 >= 0 and k3 >= 0)
        if not shipped_ok: continue
        tested += 1
        other_fail = any(min(rc.ksigns(a, mt, c)[:3]) < 0 for c in (0, 1, 2) if c != case)
        if other_fail: shipped_ok_otherfail += 1
        tv = min(rc.lam1_2x2(w, mt, c) for w in np.linspace(a, b, 21) for c in (0, 1, 2))
        if tv < float(mu2) - 1e-13:
            found.append((mt, a, b, case, float(mu2), tv, other_fail))
print(f"cells with shipped k-sign test passing: {tested}; of those with ANOTHER case failing k-signs: {shipped_ok_otherfail}")
print(f"bound violations among shipped-passing cells: {len(found)}")
for f in found[:10]:
    print(f"  mt={f[0]} cell [{f[1]:.4f},{f[2]:.4f}] best case={f[3]} bound={f[4]:.6f} true min={f[5]:.6f} other_case_kfail={f[6]}")
