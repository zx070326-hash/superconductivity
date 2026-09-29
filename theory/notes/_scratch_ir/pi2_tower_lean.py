# -*- coding: utf-8 -*-
"""Lean: recognize kappa_2 (weight 3) and compute kappa_3 (weight 4) -> b3, beta3.
   Flush output; modest maxsteps.
"""
import mpmath as mp
mp.mp.dps = 45

def G(s):   return 1/(1+s**2)
def d(s):   return 1 - mp.atan(s)/s
def Ups(s): return (3*mp.log(1+s**2) + 2*(s**2-2)*mp.atan(s)/s + (4+s**2)) \
                   / (2*(1+s**2)*(4+s**2))
def Kf(s):  return Ups(s) + d(s)*G(s)
def R(s,sg):
    return 2*s**2*sg**2*(5+s**2+sg**2)/((1+(s-sg)**2)*(1+(s+sg)**2)*(1+s**2)*(1+sg**2))

BRK = [0, mp.mpf(1)/4, mp.mpf(1)/2, 1, 2, 5, 20, mp.inf]
k2 = mp.quad(lambda s: Kf(s)**2/(2*s), BRK)
print("kappa_2 =", mp.nstr(k2, 40), flush=True)

pi=mp.pi; ln2=mp.log(2); ln3=mp.log(3); z3=mp.zeta(3); Gc=mp.catalan
zoo = {"1":mp.mpf(1),"pi^2":pi**2,"pi^3":pi**3,"zeta3":z3,
       "pi^2 ln2":pi**2*ln2,"pi^2 ln3":pi**2*ln3,"pi ln2^2":pi*ln2**2,
       "pi ln2 ln3":pi*ln2*ln3,"ln2^3":ln2**3,"ln2^2 ln3":ln2**2*ln3,
       "pi G":pi*Gc,"pi^3(a)":pi**3}
def tryp(subset, maxc=10**7):
    b=[k2]+[zoo[n] for n in subset]
    r=mp.pslq(b, tol=mp.mpf(10)**-33, maxcoeff=maxc, maxsteps=20000)
    print("  k2 /", subset, "->", r, flush=True)
tryp(["1","pi^2","zeta3"])
tryp(["1","pi^2","pi^3","zeta3"])
tryp(["1","pi^2","pi^3","zeta3","pi^2 ln2","pi ln2^2","ln2^3"])
tryp(["1","pi^2","pi^3","zeta3","pi^2 ln2","pi^2 ln3","pi ln2 ln3","ln2^2 ln3","pi G"])

print("\n--- kappa_3 (weight 4) ---", flush=True)
mp.mp.dps = 30
BRK2 = [0, mp.mpf(1)/3, 1, 3, 12, mp.inf]
def K2f(s):
    return mp.quad(lambda sg: R(s,sg)*Kf(sg)/(2*sg), BRK2) + d(s)*Kf(s)
k3 = mp.quad(lambda s: Kf(s)*K2f(s)/(2*s), BRK2)
print("kappa_3 =", mp.nstr(k3, 22), flush=True)

b1=49*pi**2/144-3; C=mp.euler+mp.log(4)-mp.mpf(1)/2; k1=b1/2
b2=2*k2-b1*C
b3=2*k3-4*k2*C+2*k1*(C**2-b1)
print("beta2=2k2 =", mp.nstr(2*k2,20), " b2 =", mp.nstr(b2,20), flush=True)
print("beta3=2k3 =", mp.nstr(2*k3,20), " b3 =", mp.nstr(b3,20), flush=True)
# is beta3 (weight 4) in Q[pi^2] (i.e. rational*pi^4+rational*pi^2+rational)?
pi4=pi**4
r=mp.pslq([2*k3, mp.mpf(1), pi**2, pi4], tol=mp.mpf(10)**-20, maxcoeff=10**6, maxsteps=20000)
print("  beta3 / {1,pi^2,pi^4} ->", r, "  (weight-4: in Q[pi^2]?)", flush=True)
r2=mp.pslq([2*k3, mp.mpf(1), pi**2, pi4, z3*ln2, pi**2*ln2**2], tol=mp.mpf(10)**-20, maxcoeff=10**6, maxsteps=20000)
print("  beta3 / {1,pi^2,pi^4,zeta3 ln2,pi^2 ln2^2} ->", r2, flush=True)
