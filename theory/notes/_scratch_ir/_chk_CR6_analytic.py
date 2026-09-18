"""[analytic + numerical-evidence] C_R <= 6, proved, and the 12-cell cost of it.

Residues of Phi_s (sympy, exact):
    r(s+i) = -i/(4(s+i)),  r(s-i) = i/(4(s-i)),  r(-s+i) = i/(4(s-i)),
    r(-s-i) = -i/(4(s+i)),  r(i) = r(-i) = 1/(2(s^2+1))
=>  sum_j |r_j| = 4 * 1/(4 sqrt(s^2+1)) + 2 * 1/(2(s^2+1))
                = 1/sqrt(1+s^2) + 1/(1+s^2)
=>  s V2[Phi_s] <= 4 s sum|r_j| = 4s/sqrt(1+s^2) + 4s/(1+s^2) <= 4 + 2 = 6
    (first term < 4; second <= 2 by s^2+1 >= 2s).   ALL of it elementary.

Part 1 verifies the residue sum symbolically; part 2 runs the 12-cell criterion
at C_R=6 against the C_R=1.6 / 20 / inf columns.
Run: OMP_NUM_THREADS=1 python _chk_CR6_analytic.py
"""
import os, io, contextlib, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, sympy as sp

s_ = sp.symbols("s", positive=True)
claim = 1/sp.sqrt(1+s_**2) + 1/(1+s_**2)
terms = [sp.Abs(-sp.I/(4*s_+4*sp.I)), sp.Abs(sp.I/(4*(s_-sp.I))),
         sp.Abs(sp.I/(4*(s_-sp.I))), sp.Abs(-sp.I/(4*s_+4*sp.I)),
         sp.Abs(1/(2*(s_**2+1))), sp.Abs(1/(2*(s_**2+1)))]
tot = sp.simplify(sum(terms))
print("sum_j |r_j| =", tot)
print("  minus claim ->", sp.simplify(tot - claim))
f = 4*s_*claim
print("  4 s sum|r| =", sp.simplify(f))
print("  sup over s>0 (sympy maximum):", sp.N(sp.maximum(f, s_, sp.Interval.open(0, sp.oo)), 15))
print("  elementary bound: 4s/sqrt(1+s^2) < 4 ; 4s/(1+s^2) <= 2  =>  <= 6")
print("  4s/(1+s^2) <= 2 residual (s-1)^2 >= 0:",
      sp.simplify(2 - 4*s_/(1+s_**2) - (s_-1)**2/(1+s_**2)))

HERE = os.path.dirname(os.path.abspath(__file__))
def _load(name, alias):
    spc = importlib.util.spec_from_file_location(alias, os.path.join(HERE, name))
    m = importlib.util.module_from_spec(spc)
    with contextlib.redirect_stdout(io.StringIO()):
        spc.loader.exec_module(m)
    return m
gc = _load("_chk_grid_coarsen.py", "gc")
pf = _load("_chk_padR_full.py", "pf")
TS = gc.THETA_STAR

print("\n== 12-cell worst at C_R = 6 (the PROVED constant) ==")
CRS = [1.6, 6.0, 20.0, np.inf]
print(f"{'cell':>15} |" + "".join(f"{f'C_R={c:g}':>10}" for c in CRS))
worst = {c: -np.inf for c in CRS}
for lo, hi, al in pf.CELLS:
    row = {c: -np.inf for c in CRS}
    for w in np.geomspace(lo, hi, 11):
        w = float(w); bs = gc.base(w); h = bs["h"]
        Pb = float(gc.peel_all_betas(w).min())
        for c in CRS:
            Lf = np.inf if np.isinf(c) else pf.lam_full(w, c)
            row[c] = max(row[c], gc.ratio(bs, al, min(Pb, Lf, TS, 0.6*h)))
    for c in CRS:
        worst[c] = max(worst[c], row[c])
    print(f"[{lo:5.2f},{hi:5.2f}] |" + "".join(f"{row[c]:10.4f}" for c in CRS))
print(f"{'12-cell worst':>15} |" + "".join(f"{worst[c]:10.4f}" for c in CRS))
