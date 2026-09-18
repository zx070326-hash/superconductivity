# rev7: build M2 once, dump every constant to JSON, and test the eps1/eps2 formulas against the SUM OF THE (L3) BOUNDS
import sys, os, json, time, numpy as np
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
import p4_cii_m2 as m2
t0 = time.time()
m2.build_profiles(verbose=False); m2.build_constants(verbose=False); m2.build_sups(verbose=False)
print(f"built {time.time()-t0:.0f}s")
out = dict(CV=m2.CV, TBAR=m2.TBAR, THETA_STAR=m2.THETA_STAR, S_STAR=m2.S_STAR, V_FPP=m2.V_FPP)
# Schur sups (top, head, tail) at a set of wbar
sch = {}
for w in (8.4761, 8.5, 10.0, 20.0, 45.0, 100.0, 1000.0):
    ex = m2.exact_sums(w); sch[str(w)] = {}
    for nm, Sf in (('phi', ex['S_phi']), ('kap', ex['S_kap'])):
        sch[str(w)][nm] = m2.schur_sup(nm, ex['A'], Sf, w)
    sch[str(w)]['ex'] = ex
out['schur'] = sch
json.dump(out, open(r"D:\physics\superconductivity\theory\notes\_scratch_ir\rev7_consts.json", "w"), indent=1)
print("CV:"); [print(f"  {k:12s} {v: .10f}") for k, v in m2.CV.items()]
print("TBAR:", m2.TBAR); print("THETA*", m2.THETA_STAR, m2.S_STAR)
# ---- eps1/eps2: the formula vs the sum of the per-n BOUND eps_f(s_n) (what a proof must dominate), and vs the true error
print("\n[eps] w | f | sum w G^2 eps_f(s_n) / eps1 | true / eps1 | sum w G^2 eps_f^2 / eps2 | true / eps2")
for w in (6.0, 8.0, 8.4761, 10.0, 20.0, 45.0, 100.0, 1000.0):
    D = 1/w; q = D*D/8; N = int(400*w); n, wn, sn = m2.lattice(N, w)
    for nm in ('phi', 'kap'):
        f = m2.FAM[nm][0]; V = m2.V_FPP[nm]
        lf = 2.0*np.concatenate(([0.0], np.cumsum(f(np.arange(1, N)/w))))/(2*n+1)
        e_true = np.abs(lf - m2.LAM_F[nm](sn)/sn)
        eps = (m2.LAM_F[nm](D/2) + q*np.minimum(2.0*n*D, V))/sn
        # also the docstring's literal form with s*sup|f''| (not (s-D/2))
        eps_doc = (m2.LAM_F[nm](D/2) + q*np.minimum(2.0*sn, V))/sn
        G2 = wn*m2.G(sn)**2
        S1b, S1t = float(np.sum(G2*eps)), float(np.sum(G2*e_true))
        S2b, S2t = float(np.sum(G2*eps**2)), float(np.sum(G2*e_true**2))
        S1d = float(np.sum(G2*eps_doc))
        print(f"  {w:8g} {nm:3s} | {S1b/m2.eps1(nm,w):.4f} (doc-form {S1d/m2.eps1(nm,w):.4f}) | {S1t/m2.eps1(nm,w):.4f} | {S2b/m2.eps2(nm,w):.4f} | {S2t/m2.eps2(nm,w):.4f}")
# ---- D^2/3 term: sup_n eps_f(s_n) vs D^2/3
print("\n[D^2/3] w | f | sup_n eps_f(s_n) / (D^2/3) | sup_n true / (D^2/3)")
for w in (8.0, 45.0, 1000.0):
    D = 1/w; q = D*D/8; N = int(400*w); n, wn, sn = m2.lattice(N, w)
    for nm in ('phi', 'kap'):
        f = m2.FAM[nm][0]; V = m2.V_FPP[nm]
        lf = 2.0*np.concatenate(([0.0], np.cumsum(f(np.arange(1, N)/w))))/(2*n+1)
        e_true = np.abs(lf - m2.LAM_F[nm](sn)/sn)
        eps = (m2.LAM_F[nm](D/2) + q*np.minimum(2.0*n*D, V))/sn
        print(f"  {w:6g} {nm} | {eps.max()/(D*D/3):.4f} | {e_true.max()/(D*D/3):.4f}")
print(f"[done {time.time()-t0:.0f}s]")
