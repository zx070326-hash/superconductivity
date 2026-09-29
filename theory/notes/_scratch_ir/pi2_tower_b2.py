# -*- coding: utf-8 -*-
"""b2 via the moment reduction  delta = 2 sum_{j>=1} kappa_j / y^j,  y = h+1 = 2B + delta.
   kappa_j = <tau, K^j tau>,  K = calR + (I - D),  (K tau)(s) = Ups_G(s) + d(s) G(s).
   Scaling limit:  <tau, f g> -> (1/2) int_0^oo f(s) g(s) ds/s   (since w_n = D/(2 s_n)).
"""
import mpmath as mp
mp.mp.dps = 50

def G(s):   return 1/(1+s**2)
def d(s):   return 1 - mp.atan(s)/s
def Ups(s): return (3*mp.log(1+s**2) + 2*(s**2-2)*mp.atan(s)/s + (4+s**2)) \
                   / (2*(1+s**2)*(4+s**2))
def Kf(s):  return Ups(s) + d(s)*G(s)          # (K tau)(s)

BRK = [0, mp.mpf(1)/4, mp.mpf(1)/2, 1, 2, 5, 20, mp.inf]

# kappa_1 = (1/2) int G * Kf / s   (must equal b1/2)
k1 = mp.quad(lambda s: G(s)*Kf(s)/(2*s), BRK)
# kappa_2 = (1/2) int Kf^2 / s
k2 = mp.quad(lambda s: Kf(s)**2/(2*s), BRK)

b1 = 49*mp.pi**2/144 - 3
gamma = mp.euler
C = gamma + mp.log(4) - mp.mpf(1)/2          # 2B = ln wbar + C ;  c0 = C - 1
b2 = 2*k2 - b1*C

print("kappa_1                =", mp.nstr(k1, 40))
print("b1/2                   =", mp.nstr(b1/2, 40))
print("kappa_1 - b1/2         =", mp.nstr(k1 - b1/2, 6), "  (sanity: K tau scaling fn correct)")
print()
print("kappa_2                =", mp.nstr(k2, 40))
print("2*kappa_2              =", mp.nstr(2*k2, 40))
print("C = gamma+ln4-1/2      =", mp.nstr(C, 40))
print("b1*C                   =", mp.nstr(b1*C, 40))
print("b2 = 2k2 - b1*C        =", mp.nstr(b2, 40))
print("   (b1 doc 3-pt fit had b2 ~ -0.320)")
print()

def show_pslq(name, x, basis, names, tol=mp.mpf(10)**-38, maxc=10**6):
    rel = mp.pslq([x] + basis, tol=tol, maxcoeff=maxc, maxsteps=10**5)
    print(name, "PSLQ over {", ",".join(names), "} ->", rel)
    if rel and rel[0] != 0:
        c0 = rel[0]
        expr = " + ".join("%s/%s*%s" % (-rel[i+1], c0, names[i]) for i in range(len(basis)) if rel[i+1]!=0)
        print("      =>", name, "=", expr)
    return rel

pi2 = mp.pi**2; pi4 = mp.pi**4; z3 = mp.zeta(3); ln2 = mp.log(2)
cat = mp.catalan; pi = mp.pi; pi3 = mp.pi**3

print("=== is kappa_2 in Q[pi^2] ? ===")
show_pslq("kappa_2", k2, [mp.mpf(1), pi2],                 ["1","pi^2"])
show_pslq("kappa_2", k2, [mp.mpf(1), pi2, pi4],            ["1","pi^2","pi^4"])
print("--- richer (weight<=3) basis ---")
show_pslq("kappa_2", k2, [mp.mpf(1), pi2, z3],             ["1","pi^2","zeta3"])
show_pslq("kappa_2", k2, [mp.mpf(1), pi2, z3, pi2*ln2, ln2**3, pi**3],
          ["1","pi^2","zeta3","pi^2 ln2","ln2^3","pi^3"])
print()
print("=== is b2 in Q[pi^2] ? (structural claim predicts NO: b2 = 2k2 - b1*(gamma+ln4-1/2)) ===")
show_pslq("b2", b2, [mp.mpf(1), pi2],                      ["1","pi^2"])
show_pslq("b2", b2, [mp.mpf(1), pi2, pi4],                 ["1","pi^2","pi^4"])
print("--- allow gamma (should catch the -b1*gamma term) ---")
show_pslq("b2", b2, [mp.mpf(1), pi2, gamma, gamma*pi2, ln2, ln2*pi2],
          ["1","pi^2","gamma","gamma pi^2","ln2","ln2 pi^2"])
print()
print("=== identify() free-form ===")
print("kappa_2 identify:", mp.identify(k2, ['pi**2','zeta(3)']))
print("b2      identify:", mp.identify(b2, ['pi**2']))
print("b2+b1*C identify:", mp.identify(b2 + b1*C, ['pi**2','zeta(3)']), " (=2 kappa_2)")
