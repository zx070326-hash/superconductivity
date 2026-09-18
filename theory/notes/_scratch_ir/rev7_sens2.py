# rev7: (a) eps1 inflated (doc-form per-n bound sums to 1.36 eps1) -> wbar_A ; (b) eps2 vs sum of per-n bounds squared WITHOUT G^2 ;
# (c) own-quadrature Schur profile Phi_f vs script grid top ; (d) J1 head-term bug size
import sys, time, numpy as np, scipy.integrate as si
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
exec(open(r"D:\physics\superconductivity\theory\notes\_scratch_ir\rev7_asm_a.py").read())
# (a) wrap crit with an eps1 multiplier by monkeypatching VF? eps1 uses VF only through formula; emulate by scaling: define crit2
def crit_eps(w, mult1=1.0, mult2=1.0):
    D = 1/w; q = D*D/8; A = Aj(w); A1, B = A[1], A[2]
    S_Gphi = A[3]-A[4]; S_Gkap = A[2]-A[3]; S_kk = A[2]-2*A[3]+A[4]; S_pp = A[4]-2*A[5]+A[6]; S_phi = A[2]-A[3]; S_kap = A[1]-A[2]
    C = CV; tb = TB; th = TH
    rho = 2*B-1; g = rho - th - q*tb["G"]
    beta = (D/2 - np.arctan(D/2)) + q*VF["G"]
    ND = C["nD2"] + q*C["V_nD2"] + 2*beta*(C["J1"] + D/2*C["V_J1"]) + beta**2*Z3H/(2*D*D)
    NR = np.sqrt(C["nR2"] + q*C["V_nR2"]) + q*np.sqrt(C["W_G"] + D/2*C["V_W_G"])
    Cr = (NR + np.sqrt(ND))/np.sqrt(B); sn = Cr/g
    e1 = lambda nm: mult1*q*(np.log(1+VF[nm]*w) + VF[nm]/(VF[nm]+D))
    e2 = lambda nm: mult2*(Z3H*D**4/576 + q*q*(4*np.log(1+VF[nm]*w)+2))
    padX = lambda nm: q*C["V_"+nm+"_GU"] + q*(C["GU_"+nm] + D/2*C["V_GU_"+nm])
    lo, hi = {}, {}
    for nm, Sgf in (("phi", S_Gphi), ("kap", S_Gkap)):
        pX = padX(nm); Yhi = C[nm+"_tDt"] + q*C["V_"+nm+"_tDt"] + e1(nm); Ylo = max(0.0, C[nm+"_tDt"] - q*C["V_"+nm+"_tDt"] - e1(nm))
        lo[nm] = 4*Sgf + (C[nm+"_GU"] - pX - Yhi)/B; hi[nm] = 4*Sgf + (C[nm+"_GU"] + pX - Ylo)/B
    Ft = lo["phi"]; qk = max(abs(lo["kap"]), abs(hi["kap"]))
    e = {}
    for nm, Sff, Sgf in (("phi", S_pp, S_Gphi), ("kap", S_kk, S_Gkap)):
        V = lambda k: q*C["V_"+nm+"_"+k]
        nZ2 = 4*B*B*Sff - 4*B*(C[nm+"_fDt"]-V("fDt")) + (C[nm+"_DD"]+V("DD")) + 4*B*(C[nm+"_fU"]+V("fU")) - 2*(C[nm+"_lU"]-V("lU")) + (C[nm+"_UU"]+V("UU"))
        tZ = 2*B*Sgf - (C[nm+"_tDt"]+V("tDt")) + (C[nm+"_GU"]-V("GU"))
        dZ = np.sqrt(e2(nm)) + q*np.sqrt(C["W_"+nm] + D/2*C["V_W_"+nm])
        e[nm] = (np.sqrt(max(0.0, nZ2 - max(0.0, tZ)**2/B)) + dZ)/np.sqrt(B)
    NS = {nm: schur_top(nm, A1, Sf) + q*tb[nm] + D*D/3 for nm, Sf in (("phi", S_phi), ("kap", S_kap))}
    Flo = (1-sn*sn)*Ft - 2*sn*e["phi"] - sn*sn*NS["phi"]; Ehi = e["kap"] + sn*(qk + NS["kap"])
    return Ehi**2/(g*Flo)
print("[a] eps1/eps2 multipliers -> wbar_A, ratio@45:")
for m1_, m2_ in ((1.0, 1.0), (1.36, 1.0), (2.0, 1.0), (1.0, 2.0), (1.0, 5.0), (1.36, 2.0)):
    f = lambda w: crit_eps(w, m1_, m2_); print(f"    eps1 x{m1_}, eps2 x{m2_}: wbar_A = {bis(f):.4f}  ratio@45 = {f(45.0):.6f}  ratio@8.5 = {f(8.5):.6f}")
# (b) eps2 vs sum over n of the per-n bound squared, with and without G^2
print("\n[b] sum_n w_n [G^2] eps_f(s_n)^2 / eps2 (code per-n form):")
for w in (6.0, 8.5, 45.0, 1000.0, 1e5):
    D = 1/w; q = D*D/8; N = int(min(400*w, 2e7)); n = np.arange(N, dtype=float); sn = (n+0.5)*D; wn = 1/(2*n+1)
    for nm in ("phi", "kap"):
        V = VF[nm]; eps = (m2.LAM_F[nm](D/2) + q*np.minimum(2*n*D, V))/sn; e2 = Z3H*D**4/576 + q*q*(4*np.log(1+V*w)+2)
        print(f"    w={w:8g} {nm}: with G^2 {np.sum(wn/(1+sn*sn)**2*eps**2)/e2:.4f}   without {np.sum(wn*eps**2)/e2:.4f}")
# (c) own quadrature of Phi_f on [0.05,30] vs script top
Gn = lambda v: 1/(1+v*v); KAn = lambda v: v*v/(1+v*v)**2; PHn = lambda v: v*v/(1+v*v)**3
LK = lambda v: (np.arctan(v) - v/(1+v*v))/2; LP = lambda v: (np.arctan(v) + v*(v*v-1)/(1+v*v)**2)/8
def Theta_f(f, sv):
    g = lambda tt: (f(abs(sv-tt)) + f(sv+tt) - 2*f(sv)*Gn(tt) - 2*f(tt)*Gn(sv))/tt
    pts = sorted(set([1e-12, 0.5, 1.0, sv/2, sv, 2*sv, 4*sv, 40.0]))
    return 0.5*(sum(si.quad(g, pts[i], pts[i+1], limit=300)[0] for i in range(len(pts)-1)) + si.quad(g, pts[-1], np.inf, limit=300)[0])
print("\n[c] own-quadrature sup of Phi_f on 600 pts in [0.05,30] vs script grid top (chord+curv):")
sg = np.concatenate([np.linspace(0.05, 3.0, 400), np.linspace(3.02, 30, 200)])
TP = np.array([Theta_f(PHn, s) for s in sg]); TK = np.array([Theta_f(KAn, s) for s in sg])
for w in (8.4761, 45.0, 1000.0):
    A = Aj(w); ex = m2.exact_sums(w)
    for nm, f, T, L, Sf in (("phi", PHn, TP, LP, ex["S_phi"]), ("kap", KAn, TK, LK, ex["S_kap"])):
        Phi = 2*A[1]*f(sg) + T + 2*Gn(sg)*Sf + L(sg)/sg; i = int(np.argmax(Phi))
        top = schur_top(nm, A[1], Sf)
        print(f"    w={w:8g} {nm}: own max {Phi[i]:.6f} at s={sg[i]:.3f};  script top {top:.6f};  own/top = {Phi[i]/top:.6f}")
# (d) J1 head
s0 = m2._DS[0]; h0 = (1-np.arctan(s0)/s0)*Gn(s0)**2/(2*s0*s0)
print(f"\n[d] _int head term: s0={s0:.1e}, integrand of J1 at s0 = {h0:.6f} (limit 1/6); head added = {0.5*s0*h0:.2e}, true head ~ s0*h(0) = {s0/6:.2e}; J1 short by ~{s0/6 - 0.5*s0*h0:.2e} (M2 J1={CV['J1']:.8f}, mpmath 0.10704534)")
print(f"[done {time.time()-t0:.0f}s]")
