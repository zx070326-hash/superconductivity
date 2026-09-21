import sys, numpy as np
exec(open(r"D:\physics\superconductivity\theory\notes\_scratch_ir\rev7_asm_a.py").read())
for w in (8.4761, 45.0):
    A = Aj(w); ex = m2.exact_sums(w)
    for nm, Sf in (("phi", ex["S_phi"]), ("kap", ex["S_kap"])):
        d = SCH[nm]; F = 2*A[1]*d["f"] + d["psi"] + 2*Sf*d["G"]; F2 = np.abs(2*A[1]*d["f2"] + d["psi2"] + 2*Sf*d["G2"])
        ds = np.diff(SS); corr = ds**2*1.5*np.maximum(F2[:-1], F2[1:])/8
        i = int(np.argmax(F)); print(f"w={w} {nm}: grid max F {F[i]:.6f} at s={SS[i]:.4f}; max|F2| on grid {F2.max():.3f}; max cell corr {corr.max():.2e}; corr at argmax {corr[min(i,len(corr)-1)]:.2e}; global-sup corr (ds^2 max|F2|/8) {(ds.max()**2*F2.max()/8):.2e}")
for w in (8.6, 8.8, 9.0, 9.5):
    print(f"ratio({w}) = {crit(w):.4f}")
