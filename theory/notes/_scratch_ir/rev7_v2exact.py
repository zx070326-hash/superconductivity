# rev7: mpmath recomputation (no splines, no trapezoid) of the closed-form-integrand constants I[u], V2[u] used by M2, vs CV
import mpmath as mp, numpy as np, time, json
t0 = time.time(); mp.mp.dps = 20
CV = json.load(open(r"D:\physics\superconductivity\theory\notes\_scratch_ir\rev7_consts.json"))["CV"]
G = lambda s: 1/(1+s*s); KA = lambda s: s*s/(1+s*s)**2; PH = lambda s: s*s/(1+s*s)**3
LP = lambda s: (mp.atan(s)+s*(s*s-1)/(1+s*s)**2)/8; LK = lambda s: (mp.atan(s)-s/(1+s*s))/2
d = lambda s: 1 - mp.atan(s)/s
U = {"nD2": lambda s: d(s)**2*G(s)**2, "phi_tDt": lambda s: LP(s)/s*G(s)**2, "kap_tDt": lambda s: LK(s)/s*G(s)**2,
     "phi_fDt": lambda s: PH(s)*LP(s)/s*G(s), "kap_fDt": lambda s: KA(s)*LK(s)/s*G(s),
     "phi_DD": lambda s: (LP(s)/s)**2*G(s)**2, "kap_DD": lambda s: (LK(s)/s)**2*G(s)**2}
def I_and_V2(u):
    h = lambda s: u(s)/(2*s)
    I = mp.quad(h, [0, 0.5, 1, 2, 5, 20, mp.inf])
    h2 = lambda s: mp.diff(h, s, 2)
    xs = np.geomspace(1e-3, 300, 3000); vals = [float(h2(mp.mpf(x))) for x in xs]
    pts = [mp.mpf("1e-3")] + [mp.mpf(float(xs[i])) for i in range(len(xs)-1) if vals[i]*vals[i+1] < 0] + [mp.mpf(300)]
    V = mp.quad(lambda s: abs(h2(s)), pts) + mp.quad(lambda s: abs(h2(s)), [300, mp.inf]) + abs(mp.quad(h2, [0, mp.mpf("1e-3")]))
    return I, V, len(pts)-2
for k, u in U.items():
    I, V, nz = I_and_V2(u)
    print(f"  {k:8s}: I = {mp.nstr(I, 12)} (M2 {CV[k]:.10f}, diff {float(I)-CV[k]:+.2e})   V2 = {mp.nstr(V, 10)} (M2 {CV['V_'+k]:.10f}, rel diff {(float(V)-CV['V_'+k])/CV['V_'+k]:+.2e}, {nz} sign changes)   [{time.time()-t0:.0f}s]", flush=True)
h = lambda s: d(s)*G(s)**2/s/(2*s); I = mp.quad(h, [0, 1, 5, mp.inf]); h1 = lambda s: mp.diff(h, s, 1)
xs = np.geomspace(1e-3, 300, 2000); vals = [float(h1(mp.mpf(x))) for x in xs]; pts = [mp.mpf("1e-3")]+[mp.mpf(float(xs[i])) for i in range(len(xs)-1) if vals[i]*vals[i+1] < 0]+[mp.mpf(300)]
V1 = mp.quad(lambda s: abs(h1(s)), pts) + mp.quad(lambda s: abs(h1(s)), [300, mp.inf]) + abs(mp.quad(h1, [0, mp.mpf("1e-3")]))
print(f"  J1: I = {mp.nstr(I, 12)} (M2 {CV['J1']:.10f})  V1 = {mp.nstr(V1, 10)} (M2 {CV['V_J1']:.10f}); h(0+) = {mp.nstr(mp.limit(h, 0), 8)}  (V1 = 1/6 exactly iff h monotone: {nz} sign changes of h1 found: {len(pts)-2})")
print(f"[done {time.time()-t0:.0f}s]")
