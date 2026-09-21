# rev7: where do the inner-sum pad profiles T_f, U_f attain their sup, and are they decreasing beyond the grid end (240)?
import sys, numpy as np
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
import p4_cii_m2 as m2
for nm in ("G", "phi", "kap"):
    for wG in (False, True):
        ss = [0.3, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 60.0, 120.0, 240.0, 500.0, 1000.0]
        v = [m2.abs_sig2(nm, s, wG) for s in ss]
        i = int(np.argmax(v))
        print(f"{nm}-{'U' if wG else 'T'}: sup on the sample at s={ss[i]} val={v[i]:.4f} (x1.5={1.5*v[i]:.4f});  s=240: {v[-3]:.4f}  500: {v[-2]:.4f}  1000: {v[-1]:.4f}")
