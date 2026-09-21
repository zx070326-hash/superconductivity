# Is H(x) := h(sqrt(x)) concave in x = wbar^2 ?  (would give k(P,T) <= h(wbar_rms): Einstein is extremal at fixed <w^2>)
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "derivations"))
from p4_interchange_uniformity import H_mat
def h(w, N): return np.linalg.eigvalsh(H_mat(N, w))[-1]
xs = np.concatenate([np.geomspace(1e-3, 1, 25), np.geomspace(1.2, 3600, 40)])
H = []
for x in xs:
    w = np.sqrt(x); N = 400 if w < 2 else int(min(2500, max(800, 40*w)))
    H.append(h(w, N))
H = np.array(H)
# second differences on the nonuniform grid
viol = 0
print("   x        wbar      H(x)     H'(x)     H''(x)")
for i in range(1, len(xs)-1):
    x0,x1,x2 = xs[i-1],xs[i],xs[i+1]; y0,y1,y2 = H[i-1],H[i],H[i+1]
    d1 = (y2-y0)/(x2-x0)
    d2 = 2*((y2-y1)/(x2-x1) - (y1-y0)/(x1-x0))/(x2-x0)
    if d2 > 0: viol += 1
    if i % 4 == 0 or d2 > 0: print(f"{x1:9.4g} {np.sqrt(x1):8.3f} {y1:9.5f} {d1:10.4e} {d2:11.3e} {'CONVEX!' if d2>0 else ''}")
print("convexity violations:", viol, "of", len(xs)-2)
# Jensen consequence test: two-atom P, compare k(P,T) vs h(wbar_rms)
print("\n two-atom P: k(P) vs h(rms) (must be <= if theorem true)")
rng = np.random.default_rng(0); worst = -1
for _ in range(40):
    w1, w2 = rng.uniform(0.2, 8, 2); p = rng.uniform(0.05, 0.95)
    N = 900
    K = p*H_mat(N, w1) + (1-p)*H_mat(N, w2); k = np.linalg.eigvalsh(K)[-1]
    rms = np.sqrt(p*w1*w1 + (1-p)*w2*w2); hr = h(rms, N)
    worst = max(worst, k/hr)
print(" max k(P)/h(rms) over 40 random 2-atom P =", worst)
