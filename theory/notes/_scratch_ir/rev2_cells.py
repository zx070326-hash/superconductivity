"""rev2: rerun X0+X1 of the committed extension script via import and dump every
cell (a,b,Na,L,alpha_cell) plus per-mutilde (mu2,g,ksign,case) to a pickle for
the truth-side comparisons.  Claim side only -- no truth here."""
import sys, time, pickle, io
from contextlib import redirect_stdout
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
import p4_interval_extension as X
t0 = time.time()
buf = io.StringIO()
with redirect_stdout(buf):
    X.run_X0()
    lo_grid = X.ext_grid(X.LO_EXT, X.LO_BASE)
    X.run_X1()
out = buf.getvalue()
print(out)
# recompute the lower cells (run_X0 does not store them) -- cheap (Na=70)
lo_cells = X.alpha_side(lo_grid)
hi_cells = list(X._CELLS_HI)
def rows(cells):
    d = {}
    for mt in X.MUTILDES:
        d[mt] = X.cell_rows(cells, mt)
    return d
data = dict(lo=[(a, b, Na, float(L), float(al)) for a, b, Na, L, al in lo_cells],
            hi=[(a, b, Na, float(L), float(al)) for a, b, Na, L, al in hi_cells],
            lo_rows=rows(lo_cells), hi_rows=rows(hi_cells),
            G_LO=dict(X.G_LO), G_HI=dict(X.G_HI), PASS=X.PASS[0])
with open(r"D:\physics\superconductivity\theory\notes\_scratch_ir\rev2_cells.pkl", "wb") as f:
    pickle.dump(data, f)
print(f"dumped {len(lo_cells)} lo + {len(hi_cells)} hi cells; PASS={X.PASS[0]}; {time.time()-t0:.0f}s")
