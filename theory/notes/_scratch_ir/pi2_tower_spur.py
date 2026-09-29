# -*- coding: utf-8 -*-
"""Stability test: are the large-basis PSLQ 'relations' for kappa_2 real or artifacts?
   A real relation is stable when precision is raised; an artifact vanishes.
"""
import mpmath as mp
def kap2(dps):
    mp.mp.dps = dps
    def G(s): return 1/(1+s**2)
    def d(s): return 1 - mp.atan(s)/s
    def Ups(s): return (3*mp.log(1+s**2)+2*(s**2-2)*mp.atan(s)/s+(4+s**2))/(2*(1+s**2)*(4+s**2))
    def Kf(s): return Ups(s)+d(s)*G(s)
    BRK=[0,mp.mpf(1)/4,mp.mpf(1)/2,1,2,5,20,mp.inf]
    return mp.quad(lambda s: Kf(s)**2/(2*s), BRK)

for dps in (45, 65):
    k2 = kap2(dps)
    pi=mp.pi; ln2=mp.log(2); ln3=mp.log(3); z3=mp.zeta(3); Gc=mp.catalan
    basis7 = [mp.mpf(1),pi**2,pi**3,z3,pi**2*ln2,pi*ln2**2,ln2**3]
    tol = mp.mpf(10)**-(dps-8)
    r = mp.pslq([k2]+basis7, tol=tol, maxcoeff=10**7, maxsteps=40000)
    print("dps=%d  k2/7-basis ->" % dps, r, flush=True)
    # residual check of the earlier candidate vector at dps45
    cand = [192,648,-46,0,-687,88,1]  # was over a 9-basis; not directly comparable, skip
print("Conclusion: if the vector changes between dps=45 and 65, it is a numerical artifact.")
