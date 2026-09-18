# rev7: OWN float64 truth (own KAY-II builder, rev6 lineage) vs EVERY M2 bound at w in {8.5,20,45,200} (N=40w),
# plus the true lattice double sums vs the padded intervals at w in {8.5,10,20,45,100}
import sys, time, numpy as np
from scipy.sparse.linalg import eigsh
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
import p4_cii_m2 as m2
t0 = time.time(); m2.build_profiles(verbose=False); m2.build_constants(verbose=False); m2.build_sups(verbose=False)
print(f"built {time.time()-t0:.0f}s")
G = lambda y: 1/(1+y*y); PH = lambda y: y*y/(1+y*y)**3; KA = lambda y: y*y/(1+y*y)**2
LP = lambda s: (np.arctan(s)+s*(s*s-1)/(1+s*s)**2)/8; LK = lambda s: (np.arctan(s)-s/(1+s*s))/2
def RG(s, t): return 2*s*s*t*t*(5+s*s+t*t)/((1+(s-t)**2)*(1+(s+t)**2)*(1+s*s)*(1+t*t))
def RF(f, s, t): return f(np.abs(s-t)) + f(s+t) - 2*f(s)*G(t) - 2*f(t)*G(s)
def opK(N, f):
    n = np.arange(N); sq = np.sqrt(2.0*n+1.0); Dd = np.abs(n[:,None]-n[None,:]); S = n[:,None]+n[None,:]+1
    M = (np.where(Dd>0, f[Dd], 0.0) + f[S])/(sq[:,None]*sq[None,:]); cs = np.concatenate(([0.0], np.cumsum(2.0*f[1:N]))); M[np.diag_indices(N)] -= cs/(2.0*n+1.0); return M
def tab(N, fun, w):
    k = np.arange(0, 2*N+2, dtype=float); f = np.zeros_like(k); f[1:] = fun(k[1:]/w); return f
def topk(M, k=3):
    if M.shape[0] <= 2500:
        e = np.linalg.eigh(M); return e[0][::-1][:k], e[1][:, ::-1][:, :k]
    e, V = eigsh(M, k=k, which="LA", tol=1e-12); o = np.argsort(e)[::-1]; return e[o], V[:, o]
def opnorm(M):
    if M.shape[0] <= 2500: return np.max(np.abs(np.linalg.eigvalsh(M)))
    return abs(eigsh(M, k=1, which="LM", tol=1e-10)[0][0])
ok = lambda c: "OK " if c else "VIOL"
print("\n[6] own truth (N=40w) vs M2 bounds2.  every line: bound / truth / side")
viol = 0
for w in (8.5, 20.0, 45.0, 200.0):
    N = int(40*w); n = np.arange(N); wn = 1/(2*n+1.0); sn = (n+0.5)/w; p = np.sqrt(wn)
    H = opK(N, tab(N, G, w)); OP = opK(N, tab(N, PH, w)); OK = opK(N, tab(N, KA, w))
    e, V = topk(H, 3); h, l2 = e[0], e[1]; v = V[:,0]*np.sign(V[0,0])
    tau = p*G(sn); B = tau@tau; th = tau/np.sqrt(B); rho = th@H@th; r = H@th - rho*th
    Ft = th@OP@th; yq = OP@th; e_phi = np.linalg.norm(yq - Ft*th); qk = th@OK@th; yk = OK@th; e_kap = np.linalg.norm(yk - qk*th)
    Fv = v@OP@v; y = OK@v; y -= (v@y)*v; Ev = np.linalg.norm(y); sin_ = np.sqrt(max(0, 1-(th@v)**2))
    nOP, nOK = opnorm(OP), opnorm(OK); gap = h-l2
    b = m2.bounds2(w); TQ = m2.THETA_STAR + m2.TBAR["G"]/(8*w*w)
    tests = [("rho_lo<=h", b["rho"], h, b["rho"] <= h), ("rho_lo<=rho", b["rho"], rho, b["rho"] <= rho), ("g_lo<=TRUE GAP h-l2", b["g"], gap, b["g"] <= gap),
             ("Theta*+qT>=lam2", TQ, l2, TQ >= l2), ("C_r>=||r||", b["Cr"], np.linalg.norm(r), b["Cr"] >= np.linalg.norm(r)),
             ("s>=sin", b["s"], sin_, b["s"] >= sin_), ("F_t lo<=F(t)", b["F_t"], Ft, b["F_t"] <= Ft), ("|q_k| bound>=", b["q_k"], abs(qk), b["q_k"] >= abs(qk)),
             ("e_phi>=", b["e_phi"], e_phi, b["e_phi"] >= e_phi), ("e_kap>=", b["e_kap"], e_kap, b["e_kap"] >= e_kap), ("NP>=||O[phi]||", b["NP"], nOP, b["NP"] >= nOP),
             ("NK>=||K||", b["NK"], nOK, b["NK"] >= nOK), ("F_lo<=F(v)", b["F_lo"], Fv, b["F_lo"] <= Fv), ("E_hi>=E(v)", b["E_hi"], Ev, b["E_hi"] >= Ev)]
    rb = b["ratio"]
    print(f" w={w:5g} N={N}: h={h:.5f} l2={l2:.5f} gap={gap:.5f} true ratio={Ev**2/(gap*Fv):.4f} M2 ratio={rb:.4f}")
    for nm, bb, tt, c in tests:
        viol += (not c); print(f"    {nm:22s} {bb:10.6f} {tt:10.6f} {ok(c)}")
    print(f"    [{time.time()-t0:.0f}s]")
print(f"[6] violations: {viol}")
print("\n[2] true lattice quantities vs M2 padded intervals: value, [lo,hi] or bound, and (true - center)/pad")
for w in (8.5, 10.0, 20.0, 45.0, 100.0):
    N = min(int(400*w), 8000); n = np.arange(N); wn = 1/(2*n+1.0); sn = (n+0.5)/w; p = np.sqrt(wn); D = 1/w; q = 1/(8*w*w)
    tau = p*G(sn); B = tau@tau; b = m2.bounds2(w, detail=True); ex = m2.exact_sums(w)
    UpsG = np.zeros(N); UpsP = np.zeros(N); UpsK = np.zeros(N); rowG = np.zeros(N); rowP = np.zeros(N); rowK = np.zeros(N)
    wG = wn*G(sn)
    for i0 in range(0, N, 400):
        i1 = min(N, i0+400); S = sn[i0:i1, None]; T = sn[None, :]
        rg = RG(S, T); rp = RF(PH, S, T); rk = RF(KA, S, T)
        UpsG[i0:i1] = rg@wG; UpsP[i0:i1] = rp@wG; UpsK[i0:i1] = rk@wG; rowG[i0:i1] = rg@wn
        rowP[i0:i1] = (PH(np.abs(S-T)) + PH(S+T))@wn; rowK[i0:i1] = (KA(np.abs(S-T)) + KA(S+T))@wn
    lG = 2*np.concatenate(([0.0], np.cumsum(G(np.arange(1, N)/w))))/(2*n+1); Dg = wn + lG
    lP = 2*np.concatenate(([0.0], np.cumsum(PH(np.arange(1, N)/w))))/(2*n+1); lK = 2*np.concatenate(([0.0], np.cumsum(KA(np.arange(1, N)/w))))/(2*n+1)
    out = []
    NR = np.sum(wn*UpsG**2); ND = np.sum(wn*G(sn)**2*(1-Dg)**2)
    out.append(("N_R=||R tau||^2", NR, m2.CV["nR2"], b["NR_hi"]**2, None)); out.append(("N_D=||(I-D)tau||^2", ND, m2.CV["nD2"], b["ND_hi"], None))
    for nm, Ups, lf, f, Sgf in (("phi", UpsP, lP, PH, ex["S_Gphi"]), ("kap", UpsK, lK, KA, ex["S_Gkap"])):
        X = np.sum(wn*G(sn)*Ups); Y = np.sum(wn*G(sn)**2*lf)
        pX = m2.padX(nm, w); Yp = q*m2.CV["V_"+nm+"_tDt"] + m2.eps1(nm, w); cX = m2.CV[nm+"_GU"]; cY = m2.CV[nm+"_tDt"]
        out.append((f"X_{nm}=<tau,R_f tau>", X, cX, (cX-pX, cX+pX), pX)); out.append((f"Y_{nm}=<tau,D_f tau>", Y, cY, (cY-Yp, cY+Yp), Yp))
        val = 4*Sgf + (X-Y)/B
        out.append((f"4S+(X-Y)/B [{nm}]", val, None, (b["out"][nm+"_lo"], b["out"][nm+"_hi"]), None))
        Z = p*(2*B*f(sn) - lf*G(sn) + Ups); nZ2 = Z@Z; tZ = tau@Z; ePerp = np.sqrt(nZ2 - tZ*tZ/B)/np.sqrt(B)
        out.append((f"||Z||^2 [{nm}]", nZ2, None, None, None)); out.append((f"<tau,Z> [{nm}] >= tZ_lo", tZ, None, b["e"][nm+"_tZ"], None))
        out.append((f"e_{nm}=||(I-P)Z||/sqrtB", ePerp, None, b["e"][nm], None))
        rows = rowP if nm == "phi" else rowK; NS = b["NS"][nm]; Sf = ex["S_phi"] if nm == "phi" else ex["S_kap"]
        out.append((f"Schur sup [{nm}] (true rows+l_f)", np.max(rows + lf), None, NS, None))
        out.append((f"  ..lattice Theta_f part sup_n", np.max(rows - 2*f(sn)*ex["A"] - 2*G(sn)*Sf), None, None, None))
    out.append(("sup_n R_G row sum", np.max(rowG), m2.THETA_STAR, m2.THETA_STAR + q*m2.TBAR["G"], q*m2.TBAR["G"]))
    print(f"  w={w:5g} (N={N}, s_max={sn[-1]:.0f}):")
    for nm, val, cen, bd, pad in out:
        if isinstance(bd, tuple):
            inside = bd[0] <= val <= bd[1]; extra = f"  (true-center)/pad = {(val-cen)/pad:+.3f}" if pad else ""
            print(f"     {nm:34s} {val:.7f} in [{bd[0]:.7f},{bd[1]:.7f}] {ok(inside)}{extra}")
        elif bd is None:
            print(f"     {nm:34s} {val:.7f}")
        else:
            up = ">= tZ_lo" not in nm; c = (val <= bd) if up else (val >= bd); rel = "<=" if up else ">="
            extra = f"  (true-center)/pad = {(val-cen)/pad:+.3f}" if pad else ""
            print(f"     {nm:34s} {val:.7f} {rel} {bd:.7f} {ok(c)}{extra}")
    print(f"    [{time.time()-t0:.0f}s]")
print(f"[done {time.time()-t0:.0f}s]")
