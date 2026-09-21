# rev7: lemmas L1-L4, psi-jets, Theta* algebra, V_f'' closed forms, tail majorant derivation -- sympy / mpmath
import sympy as sp, mpmath as mp, numpy as np, time
t0 = time.time(); mp.mp.dps = 30
s, sg, t, a, b, D, y = sp.symbols('s sigma t a b D y', positive=True)
# ---- (L1) Peano: int_a^b h - D h(m) = int K h''  (check on generic polynomial degree 5 and on exp)
hfun = sp.Function('h')
m = (a+b)/2
for htest in [t**5 - 3*t**3 + t, sp.exp(-t), sp.cos(3*t)]:
    lhs = sp.integrate(htest, (t, a, b)) - (b-a)*htest.subs(t, m)
    h2 = sp.diff(htest, t, 2)
    rhs = sp.integrate((t-a)**2/2*h2, (t, a, m)) + sp.integrate((b-t)**2/2*h2, (t, m, b))
    print("[L1] Peano identity residual:", sp.simplify(lhs - rhs))
print("[L1] int_0^D K =", sp.simplify(sp.integrate((t)**2/2, (t, 0, D/2)) + sp.integrate((D-t)**2/2, (t, D/2, D))), " max K = D^2/8")
# ---- (L3): the cells [(k-1/2)D,(k+1/2)D], k=1..n, tile [D/2, s_n] -- trivial; check l_f identity numerically done in script. Check Lam_f closed forms:
G = 1/(1+y**2); KAP = y**2/(1+y**2)**2; PHI = y**2/(1+y**2)**3
for nm, f, Lam in [('G', G, sp.atan(y)), ('kap', KAP, (sp.atan(y) - y/(1+y**2))/2), ('phi', PHI, (sp.atan(y) + y*(y**2-1)/(1+y**2)**2)/8)]:
    print(f"[L3] Lam_{nm}' - f =", sp.simplify(sp.diff(Lam, y) - f), "  Lam(0)=", Lam.subs(y, 0))
    f2 = sp.simplify(sp.diff(f, y, 2))
    roots = [r for r in sp.solve(sp.numer(sp.together(f2)), y) if r.is_real and r > 0]
    roots = sorted([sp.nsimplify(r) for r in roots])
    f1 = sp.diff(f, y)
    if nm == 'G':
        V = -2*f1.subs(y, roots[0])
    else:
        V = 2*f1.subs(y, roots[0]) - 2*f1.subs(y, roots[1])
    Vnum = mp.quad(lambda yy: abs(sp.lambdify(y, f2, 'mpmath')(yy)), [0] + [float(r) for r in roots] + [mp.inf])
    print(f"      f''={f2}; positive roots {roots}; closed V={sp.nsimplify(V)} = {float(V):.13f}; mp.quad int|f''| = {mp.nstr(Vnum, 14)}; f''(0)={f2.subs(y,0)}")
    # sup |f''| = 2 ?  check the extrema of f'' 
    f3 = sp.numer(sp.together(sp.diff(f2, y)))
    crit = [float(r) for r in sp.solve(f3, y) if r.is_real and r > 0]
    vals = [abs(float(f2.subs(y, c))) for c in crit] + [abs(float(f2.subs(y, 0)))]
    print(f"      sup|f''| over critical pts {crit} -> {max(vals):.6f}")
    # Lam_f(s) <= s^3/3 for f in {kap, phi}
    if nm != 'G':
        print(f"      f <= y^2 ? min over y of y^2 - f at y=0.1,1,3:", [float((y**2 - f).subs(y, v)) for v in (0.1, 1, 3)], " (f<=y^2 <=> Lam<=y^3/3)")
# ---- (L4)
for r in (2, 3, 4):
    print(f"[L4] zeta_h({r}) = {mp.nstr(mp.nsum(lambda n: (n+mp.mpf(1)/2)**(-r), [0, mp.inf]), 20)}  vs (2^r-1)zeta(r) = {mp.nstr((2**r-1)*mp.zeta(r), 20)}",
          "  pi^2/2" if r == 2 else ("  7zeta(3)" if r == 3 else "  pi^4/6"), mp.nstr([mp.pi**2/2, 7*mp.zeta(3), mp.pi**4/6][r-2], 20))
# ---- psi-jets: A_j = S[G^j] exact, independent: mpmath nsum with high precision at w=8.4761 and 45 vs the recursion (own implementation)
def A_jet_own(w, jmax=6):
    w = mp.mpf(w); z = mp.mpc(mp.mpf(1)/2, w)
    jet = [mp.re((1j)**k * (mp.digamma(z) if k == 0 else mp.polygamma(k, z)))/2 for k in range(jmax+1)]
    jet[0] += (mp.euler + 2*mp.log(2))/2
    A = {1: jet[0]}; cur = jet
    for j in range(1, jmax):
        L = [(k*cur[k] + w*cur[k+1])/2 for k in range(len(cur)-1)]
        cur = [cur[k] - L[k]/j for k in range(len(L))]
        A[j+1] = cur[0]
    return A
mp.mp.dps = 25
for w in (8.4761, 45.0):
    A = A_jet_own(w)
    for j in (1, 2, 3, 6):
        direct = mp.nsum(lambda n: (1/(2*n+1)) * (w**2/(w**2 + (n+mp.mpf(1)/2)**2))**j, [0, mp.inf])
        print(f"[psi] w={w} A_{j}: jet {mp.nstr(A[j], 16)}  nsum {mp.nstr(direct, 16)}  diff {mp.nstr(A[j]-direct, 3)}")
# check the recursion symbolically: x d/dx G^j = j (G^j - G^{j+1}) with G = x/(x+a^2)
x, aa = sp.symbols('x a', positive=True); Gx = x/(x+aa**2)
for j in (1, 2, 5):
    print(f"[psi] x d/dx G^{j} - {j}(G^{j}-G^{j+1}) =", sp.simplify(x*sp.diff(Gx**j, x) - j*(Gx**j - Gx**(j+1))))
# ---- Theta*: Theta' = N/(1+s^2)^2 ; N' as stated
Th = (sp.log(1+s**2) + 2*s*sp.atan(s))/(2*(1+s**2))
N = 2*s + (1-s**2)*sp.atan(s) - s*sp.log(1+s**2)
print("[Theta] Theta' - N/(1+s^2)^2 =", sp.simplify(sp.diff(Th, s) - N/(1+s**2)**2))
Np = 2 + (1-3*s**2)/(1+s**2) - 2*s*sp.atan(s) - sp.log(1+s**2)
print("[Theta] N' - stated =", sp.simplify(sp.diff(N, s) - Np), "; N'(0) =", Np.subs(s, 0), "; N''=", sp.simplify(sp.diff(Np, s)), "(<0 for s>0)")
print("[Theta] N(1) =", float(N.subs(s, 1)), " N(2) =", float(N.subs(s, 2)))
mp.mp.dps = 40
Nf = lambda v: 2*v + (1-v*v)*mp.atan(v) - v*mp.log(1+v*v)
sstar = mp.findroot(Nf, 1.5); Thf = lambda v: (mp.log(1+v*v) + 2*v*mp.atan(v))/(2*(1+v*v))
print("[Theta] s* =", mp.nstr(sstar, 20), " Theta* =", mp.nstr(Thf(sstar), 20))
# ---- Tail majorant: Phi_f(s) <= M_f(s) for s in [4, 300] (own quadrature of Theta_f), and mu_f bounds
mp.mp.dps = 15
import scipy.integrate as si
Gn = lambda v: 1/(1+v*v); KAn = lambda v: v*v/(1+v*v)**2; PHn = lambda v: v*v/(1+v*v)**3
K2 = lambda v: (6*v**4 - 16*v*v + 2)/(1+v*v)**4; P2 = lambda v: (20*v**4 - 26*v*v + 2)/(1+v*v)**5
LK = lambda v: (np.arctan(v) - v/(1+v*v))/2; LP = lambda v: (np.arctan(v) + v*(v*v-1)/(1+v*v)**2)/8
def Theta_f(f, sv):
    g = lambda tt: (f(abs(sv-tt)) + f(sv+tt) - 2*f(sv)*Gn(tt) - 2*f(tt)*Gn(sv))/tt
    pts = sorted(set([1e-12, 0.5, 1.0, sv/2, sv-1, sv, sv+1, 2*sv, 4*sv]))
    tot = sum(si.quad(g, pts[i], pts[i+1], limit=300)[0] for i in range(len(pts)-1)) + si.quad(g, pts[-1], np.inf, limit=300)[0]
    return 0.5*tot
F0 = {'phi': np.pi/16, 'kap': np.pi/4}; MU = {'phi': lambda r: 26.7/r**6, 'kap': lambda r: 10.2/r**4}
FUP = {'phi': lambda v: 1/v**4, 'kap': lambda v: 1/v**2}
yy = np.geomspace(2.0, 1e4, 200001)
for nm, f, f2 in (('kap', KAn, K2), ('phi', PHn, P2)):
    # mu_f(r) = sup_{y>=r}|f''| <= MU(r)?  check on r in [3, 1000]
    absf2 = np.abs(f2(yy)); suff = np.maximum.accumulate(absf2[::-1])[::-1]
    rr = np.geomspace(3.0, 1000.0, 400); worst = max(suff[np.searchsorted(yy, r)]/MU[nm](r) for r in rr)
    print(f"[tail] {nm}: max_r sup_(y>=r)|f''| / MU(r) over r in [3,1000] = {worst:.4f}")
    # Phi_f <= M_f at a few (w, s)
    for w in (8.4761, 45.0):
        A = float(A_jet_own(w)[1]); Aj = A_jet_own(w); Sf = float(Aj[1]-Aj[2]) if nm == 'kap' else float(Aj[2]-Aj[3])
        Lf = LK if nm == 'kap' else LP
        for sv in (4.0, 10.0, 30.0, 100.0, 300.0):
            Phi = 2*A*f(sv) + Theta_f(f, sv) + 2*Gn(sv)*Sf + Lf(sv)/sv
            M = 2*A*FUP[nm](sv) + MU[nm](sv-1)/4 + 0.5*FUP[nm](sv)*(np.log(2)+np.log(sv/2)) + 0.5*FUP[nm](sv/2)*np.log(sv/2) + 3*F0[nm]/sv + 2*Sf/sv**2 + F0[nm]/sv
            print(f"       w={w} s={sv:5g}: Phi_f={Phi:.6e}  M_f={M:.6e}  ok={Phi<=M}")
print(f"[done {time.time()-t0:.0f}s]")
