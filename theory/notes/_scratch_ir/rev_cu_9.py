"""rev_cu_9: broad independent check of the three inequality steps (dense eigh),
over the WHOLE range kappa<1, not just the note's four cells."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "D:/physics/superconductivity/theory/derivations")
import numpy as np
from rev_cu_lib import RCell, uc

NN = 200
n1 = n2 = n3 = n4 = 0
w1 = w2 = w3 = 0.0
for w in (0.05,0.1,0.3,0.5,0.74,1.0,1.5,2.0,3.0,5.0,10.0,20.0,50.0):
    for Nc in (1,2,3,4,6,10,20,40,80):
        C = RCell(w, Nc, NN)
        u = uc(NN, Nc)
        for mut in np.linspace(0.005, 3.0, 60):
            kap = C.kappa(mut)
            if kap >= 1.0:
                continue
            e, V = np.linalg.eigh(C.A - 2*mut*np.outer(u, u))
            lam = float(e[-1]); v = V[:, -1]
            # Lemma 1 (iii)
            n1 += 1
            if C.mu1 - lam > 2*mut*C.c1**2 + 1e-12: w1 += 1
            # Lemma 1 (i): lam in (mu2, mu1) and simple
            n4 += 1
            if not (C.mu2 < lam < C.mu1) or (lam - float(e[-2])) <= 0: w2 += 1
            # Lemma 2
            cth = abs(float(C.V[:, -1] @ v))
            t = np.sqrt(max(0.0, 1-cth*cth))/max(cth, 1e-300)
            bd = 2*mut*C.c1*C.up/(C.gam - 2*mut*C.c1**2)
            n2 += 1
            if t > bd*(1+1e-9) + 1e-12: w3 += 1
            # Lemma 3
            q = float(v @ (C.Of @ v))
            lb = (1/(1+t*t))*(C.q1 - 2*C.b*t - C.m*t*t)
            n3 += 1
            if q < lb - 1e-11: print("   LEMMA3 VIOLATION", w, Nc, mut, q, lb)
print(f"Lemma 1(iii) mu_1-lam <= 2 mut c_1^2 : {n1} tests, {int(w1)} violations")
print(f"Lemma 1(i)   lam in (mu2,mu1), simple: {n4} tests, {int(w2)} violations")
print(f"Lemma 2      tan theta bound         : {n2} tests, {int(w3)} violations")
print(f"Lemma 3      quadratic lower bound   : {n3} tests, 0 violations unless printed above")
