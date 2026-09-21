import numpy as np, sys
sys.path.insert(0, 'theory/notes/_scratch_ir')
from probe_non_lorentz_numerics import KERNELS, h_curve, k_pow


def hpp_nonuniform(xgrid, h):
    x0, x1, x2 = xgrid[:-2], xgrid[1:-1], xgrid[2:]
    h0, h1, h2 = h[:-2], h[1:-1], h[2:]
    Hpp = 2 * ((h2 - h1) / (x2 - x1) - (h1 - h0) / (x1 - x0)) / (x2 - x0)
    return Hpp, x1


for gamma in (1.0, 1.5, 2.0, 2.2, 2.5, 3.0):
    maker = lambda w, g=gamma: (lambda k: k_pow(k, w, g))
    xgrid = np.geomspace(1e-6, 200.0, 800)
    w = np.sqrt(xgrid)
    h = h_curve(maker, w)
    Hpp, xm = hpp_nonuniform(xgrid, h)
    pos = Hpp > 1e-8
    label = "second derivative"
    if pos.any():
        lo_cross = xm[pos].max()
        msg = ("gamma=%s: VIOLATION found, %s > 0 for x below ~%.4g (w below ~%.4g); "
               "max value %.3e at x=%.3e") % (gamma, label, lo_cross, np.sqrt(lo_cross),
                                               Hpp.max(), xm[np.argmax(Hpp)])
    else:
        msg = "gamma=%s: no violation found on x in [1e-6,200]; max %s = %.3e" % (
            gamma, label, Hpp.max())
    print(msg)
