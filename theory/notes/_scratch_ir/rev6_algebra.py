# rev6: algebraic checks of (M1) identity, diagonal bookkeeping, Theta(s) closed form, V_eta, Lambda_f
import sympy as sp
s, t, u, a, x = sp.symbols('s sigma u a x', positive=True)
G = lambda y: 1/(1+y**2)
R = 2*s**2*t**2*(5+s**2+t**2)/((1+(s-t)**2)*(1+(s+t)**2)*(1+s**2)*(1+t**2))
d = sp.simplify(sp.together(G(s-t)+G(s+t)-2*G(s)*G(t) - R))
print("[1] G(s-t)+G(s+t)-2G(s)G(t)-R(s,t) simplifies to:", d)
# diagonal: 2 tau_n^2 + w_n R(s_n,s_n) - Dg_nn  must equal  w_n G(2 s_n) - l_G(n)
w = sp.Symbol('w', positive=True); l = sp.Symbol('l')
Dg = w + l
lhs = 2*w*G(s)**2 + w*R.subs(t, s) - Dg
rhs = w*G(2*s) - l
print("[1b] diagonal bookkeeping lhs-rhs:", sp.simplify(lhs-rhs))
# R >= 0 obvious (all factors positive). R(s,0)=R(0,t)=0:
print("[1c] R(s,0), R(0,t):", sp.simplify(R.subs(t,0)), sp.simplify(R.subs(s,0)))
# Theta(s) = (1/2) int_0^inf R(s,t)/t dt  closed form
Theta_claim = (sp.log(1+s**2) + 2*s*sp.atan(s))/(2*(1+s**2))
I = sp.integrate(R/t, (t, 0, sp.oo))
print("[3] sympy integral:", sp.simplify(I/2))
print("[3] difference with claim:", sp.simplify(I/2 - Theta_claim))
import mpmath as mp
mp.mp.dps = 30
for sv in [0.3, 1.5, 4.0, 20.0]:
    f = lambda tt: 2*sv**2*tt**2*(5+sv**2+tt**2)/((1+(sv-tt)**2)*(1+(sv+tt)**2)*(1+sv**2)*(1+tt**2))/tt
    num = mp.quad(f, [0, sv/2, sv, 2*sv, mp.inf])/2
    cl = (mp.log(1+sv**2)+2*sv*mp.atan(sv))/(2*(1+sv**2))
    print(f"    s={sv}: numeric {num}  closed {cl}  diff {num-cl}")
# Theta* high precision
Th = lambda ss: (mp.log(1+ss**2)+2*ss*mp.atan(ss))/(2*(1+ss**2))
sst = mp.findroot(lambda ss: mp.diff(Th, ss), 1.5)
print("[3] s* =", sst, " Theta* =", Th(sst), " stationarity atan(s)(1-s^2)-s(ln(1+s^2)-2) =", mp.atan(sst)*(1-sst**2)-sst*(mp.log(1+sst**2)-2))
# Theta decreasing beyond 2.53? derivative sign
dTh = sp.diff(Theta_claim, s)
print("[3] Theta'(2.53) =", sp.N(dTh.subs(s, 2.53)), " Theta'(2.4)=", sp.N(dTh.subs(s,2.4)), " Theta'(100)=", sp.N(dTh.subs(s,100)))
# where does Theta' change sign? numerator of Theta'
num = sp.simplify(sp.numer(sp.together(dTh)))
print("[3] numerator of Theta':", sp.factor(num))
# V_eta = int_0^inf |eta''|, eta = u/(1+u^2)^2
eta = u/(1+u**2)**2
e2 = sp.factor(sp.diff(eta, u, 2))
print("[2] eta'' =", e2)
e1 = sp.diff(eta, u)
V = (e1.subs(u,0) - e1.subs(u,1)) + (sp.limit(e1, u, sp.oo) - e1.subs(u,1))
print("[2] V_eta =", sp.simplify(V))
Xi = a*x/(a**2+x)**2
print("[2] int_0^inf Xi da =", sp.integrate(Xi, (a, 0, sp.oo)))
# Lambda_f
phi = s**2/(1+s**2)**3; kap = s**2/(1+s**2)**2
Lphi = sp.Rational(1,8)*(sp.atan(s) + s*(s**2-1)/(1+s**2)**2)
Lkap = sp.Rational(1,2)*(sp.atan(s) - s/(1+s**2))
print("[4] Lambda_phi' - phi:", sp.simplify(sp.diff(Lphi, s)-phi), " Lambda_kap' - kap:", sp.simplify(sp.diff(Lkap,s)-kap))
print("[4] I[G phi]=", sp.integrate(G(s)*phi/s,(s,0,sp.oo))/2, " I[G kap]=", sp.integrate(G(s)*kap/s,(s,0,sp.oo))/2,
      " I[phi^2]=", sp.integrate(phi**2/s,(s,0,sp.oo))/2, " I[kap^2]=", sp.integrate(kap**2/s,(s,0,sp.oo))/2)
# (M2) identity check symbolically for f = kappa: R_f definition consistency is by definition; check R_f(s,0)=0
Rf = lambda f: f(s-t)+f(s+t)-2*f(s)*G(t)-2*f(t)*G(s)
fk = lambda y: y**2/(1+y**2)**2
print("[1d] R_kap(s,0) =", sp.simplify(Rf(fk).subs(t,0)), " R_kap(0,t)=", sp.simplify(Rf(fk).subs(s,0)))
# does R_kap have negative values?
import numpy as np
Rk = sp.lambdify((s,t), Rf(fk), 'numpy')
S,T = np.meshgrid(np.linspace(0.01,6,600), np.linspace(0.01,6,600))
print("[1d] min R_kap on [0,6]^2 =", Rk(S,T).min(), " min R_phi =", sp.lambdify((s,t), Rf(lambda y: y**2/(1+y**2)**3),'numpy')(S,T).min())
