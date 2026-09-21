# rev6: independent re-assembly of the criterion and wbar_A from the script's printed constants (CONST), own A,B, own Theta_f quadrature
import numpy as np, mpmath as mp
from scipy.integrate import quad
mp.mp.dps = 25
C = dict(P1=1/12, K1=1/8, mPP=1/80, mKK=1/24,
    cDel=0.05027541, nD2=0.00865336, qR=0.12892810, nR2=0.04925665,
    phi_tDt=0.03017391, phi_fDt=0.00465127, phi_DD=0.00174685, phi_fU=-0.01263179, phi_lU=-0.00455746, phi_UU=0.01422639, phi_GU=-0.09558297,
    kap_tDt=0.03736229, kap_fDt=0.01149276, kap_DD=0.00323943, kap_fU=-0.01913140, kap_lU=-0.00593825, kap_UU=0.01546955, kap_GU=-0.09040317)
TH_STAR = 0.6349287917881
G = lambda y: 1/(1+y*y); PHI = lambda y: y*y/(1+y*y)**3; KAP = lambda y: y*y/(1+y*y)**2
LPHI = lambda s: (np.arctan(s)+s*(s*s-1)/(1+s*s)**2)/8; LKAP = lambda s: (np.arctan(s)-s/(1+s*s))/2
def AB(w):
    z = mp.mpc(0.5, w); A = (mp.re(mp.digamma(z))+mp.euler+2*mp.log(2))/2
    return float(A), float(A + w/4*mp.im(mp.polygamma(1,z)))
def Theta_f(f, s):   # (1/2) int R_f(s,t)/t dt, signed
    g = lambda t: (f(abs(s-t))+f(s+t)-2*f(s)*G(t)-2*f(t)*G(s))/t
    pts = sorted(set([1e-12, s/2, s, 2*s, 5.0, 40.0]))
    tot = sum(quad(g, pts[i], pts[i+1], limit=200)[0] for i in range(len(pts)-1)) + quad(g, pts[-1], np.inf, limit=200)[0]
    return 0.5*tot
sgrid = np.concatenate([np.linspace(0.02, 3, 150), np.linspace(3.1, 12, 60), np.geomspace(13, 200, 30)])
TP = np.array([Theta_f(PHI, s) for s in sgrid]); TK = np.array([Theta_f(KAP, s) for s in sgrid])
def schur(A, f, T, L, If): return float(np.max(2*A*f(sgrid) + T + 2*G(sgrid)*If + L(sgrid)/sgrid))
PAD = lambda w: 0.25/w**2 + 2e-6
def crit(w, pad=True):
    A, B = AB(w); pd = PAD(w) if pad else 0.0
    rho = 2*B-1; Cr = (np.sqrt(C['nR2'])+np.sqrt(C['nD2']))/np.sqrt(B) + pd
    g = rho - TH_STAR - pd; s = Cr/g
    Ft = 1/3 + (C['phi_GU']-C['phi_tDt'])/B - pd
    qk = 0.5 + (C['kap_GU']-C['kap_tDt'])/B
    e = {}
    for nm, mff, If in [('phi', C['mPP'], C['P1']), ('kap', C['mKK'], C['K1'])]:
        nZ2 = 4*B*B*mff - 4*B*C[nm+'_fDt'] + C[nm+'_DD'] + 4*B*C[nm+'_fU'] - 2*C[nm+'_lU'] + C[nm+'_UU']
        tZ = 2*B*If - C[nm+'_tDt'] + C[nm+'_GU']
        e[nm] = np.sqrt((nZ2 - tZ*tZ/B)/B) + pd
    NP = schur(A, PHI, TP, LPHI, 1/8); NK = schur(A, KAP, TK, LKAP, 1/4)
    Flo = (1-s*s)*Ft - 2*s*e['phi'] - s*s*NP; Ehi = e['kap'] + s*(abs(qk)+NK)
    return Ehi**2/(g*Flo), dict(g=g, s=s, Flo=Flo, Ehi=Ehi, NP=NP, NK=NK, e_phi=e['phi'], e_kap=e['kap'], qk=qk)
for w in (7.776, 8.0, 45.0, 1000.0):
    r, d = crit(w); print(f"w={w}: ratio={r:.5f}  " + "  ".join(f"{k}={v:.5f}" for k,v in d.items()))
lo, hi = 2.0, 60.0
for _ in range(50):
    mid = (lo+hi)/2
    if crit(mid)[0] > 1: lo = mid
    else: hi = mid
print(f"independent wbar_A = {hi:.4f}")
lo, hi = 2.0, 60.0
for _ in range(50):
    mid = (lo+hi)/2
    if crit(mid, pad=False)[0] > 1: lo = mid
    else: hi = mid
print(f"wbar_A with pad removed = {hi:.4f}")
# sensitivity: which term dominates E_hi at w=7.78 and 45
for w in (7.776, 45.0):
    r, d = crit(w); print(f"w={w}: e_kap={d['e_kap']:.4f}  s*|q|={d['s']*abs(d['qk']):.4f}  s*NK={d['s']*d['NK']:.4f}   F: (1-s^2)Ft vs 2s e_phi={2*d['s']*d['e_phi']:.4f} s^2 NP={d['s']**2*d['NP']:.4f}")
