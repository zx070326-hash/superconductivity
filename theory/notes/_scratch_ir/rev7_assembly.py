# rev7: independent re-assembly of the M2 criterion from the script CONSTANTS only (CV, TBAR, THETA*, the Schur grid arrays);
# own psi-jet A_j (mpmath), own formulas for every bound, own bisection, dense grid to 1e4, sensitivity (+-1% per constant, safety 1.5->1.0)
import sys, time, json, numpy as np, mpmath as mp
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
import p4_cii_m2 as m2
t0 = time.time(); m2.build_profiles(verbose=False); m2.build_constants(verbose=False); m2.build_sups(verbose=False)
print(f"built {time.time()-t0:.0f}s")
CV = dict(m2.CV); TB = dict(m2.TBAR); TH = m2.THETA_STAR; VF = dict(m2.V_FPP); SCH = m2.SCH; SS = m2._SS
Z3H = float(7*mp.zeta(3)); F0 = {'phi': np.pi/16, 'kap': np.pi/4}
mp.mp.dps = 25
def Aj(w, jmax=6):
    w = mp.mpf(w); z = mp.mpc(mp.mpf(1)/2, w)
    jet = [mp.re((1j)**k*(mp.digamma(z) if k == 0 else mp.polygamma(k, z)))/2 for k in range(jmax+1)]
    jet[0] += (mp.euler+2*mp.log(2))/2
    A = {1: jet[0]}; cur = jet
    for j in range(1, jmax):
        L = [(k*cur[k]+w*cur[k+1])/2 for k in range(len(cur)-1)]; cur = [cur[k]-L[k]/j for k in range(len(L))]; A[j+1] = cur[0]
    return {k: float(v) for k, v in A.items()}
def schur_top(nm, A, Sf, curv=1.5):
    d = SCH[nm]; F = 2*A*d['f'] + d['psi'] + 2*Sf*d['G']; F2 = np.abs(2*A*d['f2'] + d['psi2'] + 2*Sf*d['G2'])
    ds = np.diff(SS); cell = np.maximum(F[:-1], F[1:]) + ds**2*curv*np.maximum(F2[:-1], F2[1:])/8
    head = 2*A*d['f'][0] + max(0.0, d['psi'][0]) + 2*Sf
    s = 30.0; fu = (1/s**4, 1/(s/2)**4) if nm == 'phi' else (1/s**2, 1/(s/2)**2); mu = 26.7/(s-1)**6 if nm == 'phi' else 10.2/(s-1)**4
    tail = 2*A*fu[0] + mu/4 + 0.5*fu[0]*(np.log(2)+np.log(s/2)) + 0.5*fu[1]*np.log(s/2) + 3*F0[nm]/s + 2*Sf/s**2 + F0[nm]/s
    return max(float(cell.max()), head, tail)
def crit(w, C=CV, tb=TB, th=TH, safety_T=1.0, curv=1.5, ns_scale=1.0, detail=False):
    D = 1/w; q = D*D/8; A = Aj(w); A1, B = A[1], A[2]
    S_Gphi = A[3]-A[4]; S_Gkap = A[2]-A[3]; S_kk = A[2]-2*A[3]+A[4]; S_pp = A[4]-2*A[5]+A[6]; S_phi = A[2]-A[3]; S_kap = A[1]-A[2]
    rho = 2*B-1; g = rho - th - q*tb['G']*safety_T
    beta = (D/2 - np.arctan(D/2)) + q*VF['G']
    ND = C['nD2'] + q*C['V_nD2'] + 2*beta*(C['J1'] + D/2*C['V_J1']) + beta**2*Z3H/(2*D*D)
    NR = np.sqrt(C['nR2'] + q*C['V_nR2']) + q*np.sqrt(C['W_G'] + D/2*C['V_W_G'])
    Cr = (NR + np.sqrt(ND))/np.sqrt(B); sn = Cr/g
    def e1(nm):
        V = VF[nm]; return q*(np.log(1+V*w) + V/(V+D))
    def e2(nm):
        V = VF[nm]; return Z3H*D**4/576 + q*q*(4*np.log(1+V*w)+2)
    def padX(nm):
        return q*C['V_'+nm+'_GU'] + q*(C['GU_'+nm] + D/2*C['V_GU_'+nm])
    lo, hi = {}, {}
    for nm, Sgf in (('phi', S_Gphi), ('kap', S_Gkap)):
        pX = padX(nm); Yhi = C[nm+'_tDt'] + q*C['V_'+nm+'_tDt'] + e1(nm); Ylo = max(0.0, C[nm+'_tDt'] - q*C['V_'+nm+'_tDt'] - e1(nm))
        lo[nm] = 4*Sgf + (C[nm+'_GU'] - pX - Yhi)/B; hi[nm] = 4*Sgf + (C[nm+'_GU'] + pX - Ylo)/B
    Ft = lo['phi']; qk = max(abs(lo['kap']), abs(hi['kap']))
    e = {}
    for nm, Sff, Sgf in (('phi', S_pp, S_Gphi), ('kap', S_kk, S_Gkap)):
        V = lambda k: q*C['V_'+nm+'_'+k]
        nZ2 = 4*B*B*Sff - 4*B*(C[nm+'_fDt']-V('fDt')) + (C[nm+'_DD']+V('DD')) + 4*B*(C[nm+'_fU']+V('fU')) - 2*(C[nm+'_lU']-V('lU')) + (C[nm+'_UU']+V('UU'))
        tZ = 2*B*Sgf - (C[nm+'_tDt']+V('tDt')) + (C[nm+'_GU']-V('GU'))
        dZ = np.sqrt(e2(nm)) + q*np.sqrt(C['W_'+nm] + D/2*C['V_W_'+nm])
        e[nm] = (np.sqrt(max(0.0, nZ2 - max(0.0, tZ)**2/B)) + dZ)/np.sqrt(B)
    NS = {}
    for nm, Sf in (('phi', S_phi), ('kap', S_kap)):
        NS[nm] = ns_scale*(schur_top(nm, A1, Sf, curv) + q*tb[nm]*safety_T + D*D/3)
    Flo = (1-sn*sn)*Ft - 2*sn*e['phi'] - sn*sn*NS['phi']; Ehi = e['kap'] + sn*(qk + NS['kap'])
    r = Ehi**2/(g*Flo) if Flo > 0 else np.inf
    if detail:
        return r, dict(A=A1, B=B, g=g, Cr=Cr, s=sn, Ft=Ft, qk=qk, e_phi=e['phi'], e_kap=e['kap'], NP=NS['phi'], NK=NS['kap'], Flo=Flo, Ehi=Ehi)
    return r
def bis(fn, lo=2.0, hi=60.0, it=60):
    for _ in range(it):
        mid = (lo+hi)/2
        if fn(mid) > 1: lo = mid
        else: hi = mid
    return hi
print("\n[5] own assembly vs script bounds2:")
for w in (8.4761, 8.5, 10.0, 45.0, 1000.0):
    r, d = crit(w, detail=True); b = m2.bounds2(w)
    print(f"  w={w:8g}: own ratio {r:.6f}  script {b["ratio"]:.6f} | g {d["g"]:.5f}/{b["g"]:.5f} s {d["s"]:.5f}/{b["s"]:.5f} Ft {d["Ft"]:.5f}/{b["F_t"]:.5f} qk {d["qk"]:.5f}/{b["q_k"]:.5f} e_k {d["e_kap"]:.5f}/{b["e_kap"]:.5f} NP {d["NP"]:.5f}/{b["NP"]:.5f} NK {d["NK"]:.5f}/{b["NK"]:.5f} Flo {d["Flo"]:.5f}/{b["F_lo"]:.5f} Ehi {d["Ehi"]:.5f}/{b["E_hi"]:.5f}")
wA = bis(crit); print(f"\n[5] own wbar_A = {wA:.5f}")
grid = np.geomspace(wA*1.0001, 1e4, 3000); rr = np.array([crit(x) for x in grid])
print(f"[5] dense grid (3000 pts, [{grid[0]:.4f},1e4]): max ratio {rr.max():.5f}, monotone decreasing: {np.all(np.diff(rr) < 0)}, ratio(1e4)={rr[-1]:.5f}")
if not np.all(np.diff(rr) < 0):
    i = np.where(np.diff(rr) >= 0)[0]; print("    non-monotone at w =", grid[i][:10])
for x in (1e6, 1e10, 1e20):
    print(f"    ratio({x:g}) = {crit(x):.5f}")
print(f"[5] ratio just below/above: {crit(wA*0.999):.5f} / {crit(wA*1.001):.5f}")
print("\n[4] safety factors -> 1.0:")
for lab, kw in (("T-bar 1.5->1.0", dict(safety_T=1/1.5)), ("curv 1.5->1.0", dict(curv=1.0)), ("both", dict(safety_T=1/1.5, curv=1.0))):
    f = lambda w, kw=kw: crit(w, **kw); print(f"   {lab}: wbar_A = {bis(f):.5f}; ratio@45 = {f(45.0):.6f}")
C0 = {k: (0.0 if (k.startswith("V_")) else v) for k, v in CV.items()}
f0 = lambda w: crit(w, C=C0, tb={k: 0.0 for k in TB}); print(f"   (V2 and TBAR pads set to 0; eps, D^2/3, beta kept) wbar_A = {bis(f0):.5f}")
print("\n[4] sensitivity: ratio change for +-1 percent of each constant at w=wbar_A and 45  (ratio itself: %.5f / %.5f)" % (crit(wA), crit(45.0)))
rows = []
for k in CV:
    for sgn in (+1, -1):
        C = dict(CV); C[k] = CV[k]*(1+0.01*sgn)
        rows.append((k, sgn, crit(wA, C=C)-crit(wA), crit(45.0, C=C)-crit(45.0)))
for k in TB:
    for sgn in (+1, -1):
        tb = dict(TB); tb[k] = TB[k]*(1+0.01*sgn); rows.append(("TBAR_"+k, sgn, crit(wA, tb=tb)-crit(wA), crit(45.0, tb=tb)-crit(45.0)))
for sgn in (+1, -1):
    rows.append(("THETA*", sgn, crit(wA, th=TH*(1+0.01*sgn))-crit(wA), crit(45.0, th=TH*(1+0.01*sgn))-crit(45.0)))
    rows.append(("NS(both Schur sups)", sgn, crit(wA, ns_scale=1+0.01*sgn)-crit(wA), crit(45.0, ns_scale=1+0.01*sgn)-crit(45.0)))
rows.sort(key=lambda r: -abs(r[2]))
SG = {1: "+", -1: "-"}
for k, sgn, d1, d2 in rows[:16]:
    print(f"   {k:22s} {SG[sgn]}1%: d(ratio)@wA = {d1:+.5f}   @45 = {d2:+.5f}")
print(f"   (all other |d(ratio)@wA| < {abs(rows[16][2]):.1e})")
print("\n[4] how large an error in the single worst constant flips ratio>1 at w=8.5 / at w=45:")
for k, sgn, d1, d2 in rows[:4]:
    for w in (8.5, 45.0):
        def rat(eps):
            if k.startswith("TBAR_"):
                tb = dict(TB); tb[k[5:]] = TB[k[5:]]*(1+eps*sgn); return crit(w, tb=tb)
            if k == "THETA*": return crit(w, th=TH*(1+eps*sgn))
            if k.startswith("NS"): return crit(w, ns_scale=1+eps*sgn)
            C = dict(CV); C[k] = CV[k]*(1+eps*sgn); return crit(w, C=C)
        lo, hi = 0.0, 5.0
        if rat(hi) < 1:
            print(f"   {k} {SG[sgn]}: not even a 500% error flips at w={w}"); continue
        for _ in range(40):
            mid = (lo+hi)/2
            if rat(mid) > 1: hi = mid
            else: lo = mid
        print(f"   {k:22s} {SG[sgn]}: relative error {hi*100:.1f}% flips ratio>1 at w={w}")
print(f"[done {time.time()-t0:.0f}s]")
