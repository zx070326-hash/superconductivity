# rev6: recompute the absolute constants with mpmath (nested quad), identify
import mpmath as mp, time
mp.mp.dps = 20
G = lambda y: 1/(1+y**2)
PHI = lambda y: y**2/(1+y**2)**3
KAP = lambda y: y**2/(1+y**2)**2
LPHI = lambda s: (mp.atan(s) + s*(s*s-1)/(1+s*s)**2)/8
LKAP = lambda s: (mp.atan(s) - s/(1+s*s))/2
def Rf(f, s, t): return f(abs(s-t)) + f(s+t) - 2*f(s)*G(t) - 2*f(t)*G(s)
def RG(s, t): return 2*s*s*t*t*(5+s*s+t*t)/((1+(s-t)**2)*(1+(s+t)**2)*(1+s*s)*(1+t*t))
def I(u, pts=None):
    pts = pts or [0, 0.5, 1, 2, 5, 20, mp.inf]
    return mp.quad(lambda s: u(s)/s, pts)/2
def Ups(Rr, s, wG=True):
    g = (lambda t: G(t)*Rr(s,t)/t) if wG else (lambda t: Rr(s,t)/t)
    pts = sorted(set([0, s/2, s, 2*s, 5, 20])) + [mp.inf]
    return mp.quad(g, pts)/2
t0 = time.time()
out = {}
d = lambda s: 1 - mp.atan(s)/s
out['P1'] = I(lambda s: G(s)*PHI(s)); out['K1'] = I(lambda s: G(s)*KAP(s))
out['mPP'] = I(lambda s: PHI(s)**2); out['mKK'] = I(lambda s: KAP(s)**2)
out['cDel'] = I(lambda s: d(s)*G(s)**2); out['nD2'] = I(lambda s: d(s)**2*G(s)**2)
for nm, f, L in [('phi', PHI, LPHI), ('kap', KAP, LKAP)]:
    l = lambda s, L=L: L(s)/s
    out[nm+'_tDt'] = I(lambda s: l(s)*G(s)**2)
    out[nm+'_fDt'] = I(lambda s: f(s)*l(s)*G(s))
    out[nm+'_DD'] = I(lambda s: l(s)**2*G(s)**2)
print("1D constants (%.0fs):" % (time.time()-t0))
for k,v in out.items(): print(f"  {k:8s} {mp.nstr(v,16)}   identify: {mp.identify(v, ['pi','pi**2','log(2)','zeta(3)'])}")
# double integrals: GU for phi and kap (feed 0.1257572 and 0.1277654), qR, nR2, and the fU/lU/UU used in e_kappa/e_phi
mp.mp.dps = 15
t0 = time.time()
for nm, Rr in [('phi', lambda s,t: Rf(PHI,s,t)), ('kap', lambda s,t: Rf(KAP,s,t)), ('G', RG)]:
    U = lambda s, Rr=Rr: Ups(Rr, s, True)
    gu = I(lambda s: G(s)*U(s))
    print(f"  {nm}_GU  = {mp.nstr(gu,12)}   ({time.time()-t0:.0f}s)")
    if nm == 'G':
        nr2 = I(lambda s: U(s)**2); print(f"  nR2    = {mp.nstr(nr2,12)}   ({time.time()-t0:.0f}s)")
    else:
        f = PHI if nm=='phi' else KAP; L = LPHI if nm=='phi' else LKAP
        fu = I(lambda s: f(s)*U(s)); lu = I(lambda s: (L(s)/s)*G(s)*U(s)); uu = I(lambda s: U(s)**2)
        print(f"  {nm}_fU  = {mp.nstr(fu,12)}  {nm}_lU = {mp.nstr(lu,12)}  {nm}_UU = {mp.nstr(uu,12)}   ({time.time()-t0:.0f}s)")
print("0.1257572 claim: -(phi_GU + phi_tDt) ; 0.1277654 claim: -(kap_GU + kap_tDt)")
