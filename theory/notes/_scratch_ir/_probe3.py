import numpy as np, sys
sys.path.insert(0, 'theory/notes/_scratch_ir')
from probe_non_lorentz_numerics import KERNELS, h_curve

def hpp_nonuniform(xgrid, h):
    x0, x1, x2 = xgrid[:-2], xgrid[1:-1], xgrid[2:]
    h0, h1, h2 = h[:-2], h[1:-1], h[2:]
    Hpp = 2 * ((h2 - h1) / (x2 - x1) - (h1 - h0) / (x1 - x0)) / (x2 - x0)
    return Hpp, x1

xgrid = np.geomspace(1e-6, 20000.0, 1200)

names = [n for n in KERNELS if n.startswith("yukawa") or n.startswith("expcut")]
for name in names:
    maker = KERNELS[name]
    w = np.sqrt(xgrid)
    h = h_curve(maker, w)
    Hpp, xm = hpp_nonuniform(xgrid, h)
    pos = Hpp > 1e-8
    if pos.any():
        lo_cross = xm[pos].max()
        msg = ("%s: VIOLATION, H'' > 0 for x below ~%.4g; max=%.3e at x=%.3e"
               % (name, lo_cross, Hpp.max(), xm[np.argmax(Hpp)]))
    else:
        msg = "%s: no violation on x in [1e-6,2e4]; max H''=%.3e" % (name, Hpp.max())
    print(msg)
