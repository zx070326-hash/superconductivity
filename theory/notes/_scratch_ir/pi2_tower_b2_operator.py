# -*- coding: utf-8 -*-
"""Independent b2 from real finite-wbar top eigenvalues h(wbar)=lambda_1(H(wbar)).
   R2(w) := (h - L - c0 - b1/L) * L^2  ->  b2 + b3/L + b4/L^2 + ...
   Compare to the scaling-integral prediction b2 = 2 kappa_2 - b1*(gamma+ln4-1/2) = -0.32707...
"""
import sys, numpy as np
sys.path.insert(0, "theory/derivations")
from p4_extremal_asymptotics import top_eig, k_brk
import math

gamma = 0.5772156649015328606
c0 = gamma + math.log(4) - 1.5
b1 = 49*math.pi**2/144 - 3
b2_pred = -0.3270714675337773   # from scaling integrals

def hval(w, mult):
    N = int(mult*w)
    return top_eig(N, k_brk(N, w, 0))[0]

rows=[]
for w in [100.,200.,400.,800.,1600.,3200.,6400.]:
    h60 = hval(w, 60)
    h120 = hval(w, 120)
    # Richardson in N to kill (w/N)^2-type truncation (60 vs 120): (4*h120-h60)/3
    hR = (4*h120 - h60)/3
    L = math.log(w)
    R2 = (hR - L - c0 - b1/L)*L*L
    R2_raw = (h120 - L - c0 - b1/L)*L*L
    rows.append((w, L, h120, hR, R2, R2_raw))
    print("w=%7.0f  L=%.4f  h=%.10f  R2(Rich)=%.6f  R2(N120)=%.6f" % (w,L,h120,R2,R2_raw))

# fit R2 = b2 + b3/L + b4/L^2  using last 4 points
import numpy as np
Ls = np.array([r[1] for r in rows[-5:]])
R2s= np.array([r[4] for r in rows[-5:]])
A = np.vstack([np.ones_like(Ls), 1/Ls, 1/Ls**2]).T
coef,*_ = np.linalg.lstsq(A, R2s, rcond=None)
print()
print("fit R2 = b2 + b3/L + b4/L^2  (last 5 pts):  b2=%.5f  b3=%.4f  b4=%.3f" % (coef[0],coef[1],coef[2]))
print("scaling-integral prediction b2 =", b2_pred)
print("difference =", coef[0]-b2_pred)
