# Conjecture 2 reduction probe: is h(w) := lam1(H(w)) <= g2 * w^2 for all w>0, with g2 = rho(G(2)) in [1.31714,1.31820]?
# (convexity of lam_max + linearity of K(P,T)=int H(w)P(dw) would then give Conj. 2.)  Also: crude analytic tail
# ||H(w)|| <= 2*S1(w), S1 = sum_k w^2/(k^2+w^2) <= pi*w/2  -> H bound pi*w < g2 w^2 iff w > pi/g2.
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "derivations"))
from p4_interchange_uniformity import H_mat
g2_lo, g2_hi = 1.31714, 1.31820
print(f"pi/g2_lo = {np.pi/g2_lo:.4f}  (beyond this the trivial row-sum bound closes the inequality)")
print("   w     h(w)=lam1   g2_lo*w^2   ratio h/(g2 w^2)   2*S1(w)   pi*w")
for w in (0.05,0.1,0.2,0.3,0.5,0.7,1.0,1.4142,2.0,2.4,3.0,5.0,10.0):
    N = 1500 if w >= 1 else 400
    h = np.linalg.eigvalsh(H_mat(N, w))[-1]
    k = np.arange(1, 200000); S1 = np.sum(w*w/(k*k+w*w))
    print(f"{w:6.3f}  {h:9.5f}  {g2_lo*w*w:10.5f}  {h/(g2_lo*w*w):14.5f}   {2*S1:8.4f}  {np.pi*w:7.4f}")
