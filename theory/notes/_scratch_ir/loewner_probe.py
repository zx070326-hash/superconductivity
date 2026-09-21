# Loewner tests: (a) is x -> H(x) (x=wbar^2) operator-CONCAVE?  D = H(mid) - (H(x)+H(y))/2 >= 0 ?
#                (b) is H'(wbar) >= 0 (operator monotone in wbar)?  lam_min(H')
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "derivations"))
from p4_interchange_uniformity import H_mat
N = 600
print("(a) operator concavity in x=wbar^2: lam_min(D), D = H(sqrt((x+y)/2)) - [H(sqrt x)+H(sqrt y)]/2")
worst = np.inf
for (w1, w2) in [(0.1,0.3),(0.3,1.0),(0.5,2.0),(1.0,3.0),(2.0,8.0),(0.2,10.0),(5.0,30.0),(0.05,0.5)]:
    x, y = w1*w1, w2*w2; m = np.sqrt((x+y)/2)
    D = H_mat(N, m) - 0.5*(H_mat(N, w1) + H_mat(N, w2))
    e = np.linalg.eigvalsh(D); worst = min(worst, e[0])
    print(f"  w1={w1:5.2f} w2={w2:5.2f}: lam_min(D)={e[0]:+.3e}  lam_max(D)={e[-1]:+.3e}  diag_min={D.diagonal().min():+.3e}")
print("  worst lam_min(D) =", worst, "-> operator concave?" , worst > -1e-10)
print("(b) operator monotone in wbar: lam_min(H'(wbar)), lam_max")
for w in (0.1, 0.5, 1.0, 2.0, 5.0, 20.0):
    e = np.linalg.eigvalsh(H_mat(N, w, 1)); print(f"  w={w:5.1f}: lam_min={e[0]:+.4e}  lam_max={e[-1]:+.4e}")
