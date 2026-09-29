import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys, numpy as np
DER = "D:/physics/superconductivity/.claude/worktrees/agent-ae6c83eb3d884027f/theory/derivations"
sys.path.insert(0, DER)
from p4_concavity_q_route import O_ker, kernels, top_eigen
from p4_concavity_gap2_gates import h_lo_trial, lam2_hi, THETA_STAR


def karr(x, KX, which):
    k = np.arange(KX, dtype=float); k2 = k * k; d = k2 + x
    f = x / d; f[0] = 0.0
    if which == "f2": a = f * f
    elif which == "f3": a = f * f * f
    elif which == "f-f2": a = f - f * f
    elif which == "f-f3": a = f - f * f * f
    elif which == "g":
        a = 2.0 * x * k2 / (d ** 3); a[0] = 0.0
    a[0] = 0.0
    return a


def kL3(varr, karr_, N, KX):
    n = np.arange(N); u = 1.0 / np.sqrt(2 * n + 1.0); w = u * varr
    S = n[:, None] + n[None, :] + 1
    return float((karr_[np.minimum(S, KX - 1)] * np.outer(w, w)).sum())


def wschur(M, psi):
    return float(((np.abs(M) @ psi) / psi).max())


def perron(Mnn, iters=300):
    N = Mnn.shape[0]; psi = np.ones(N)
    for _ in range(iters):
        psi = np.abs(Mnn @ psi); psi /= np.linalg.norm(psi); psi += 1e-300
    return wschur(Mnn, psi)


def lam_max_dom(M, N, psi_kind="u"):
    di = np.diag(M).copy(); Mod = M.copy(); np.fill_diagonal(Mod, 0.0)
    Mdom = Mod + np.diag(np.maximum(di, 0.0))
    n = np.arange(N); u = 1.0 / np.sqrt(2 * n + 1.0)
    return wschur(Mdom, u), perron(Mdom)


def analyze(varpi, N=400):
    x = varpi ** 2; KX = 8 * N
    ker = kernels(x, KX=KX)
    A = O_ker(N, ker["f"]); B = O_ker(N, ker["f2"]); M = A - B
    h, v, l2 = top_eigen(A)
    E = float(v @ (O_ker(N, ker["g"]) @ v))
    q = float(v @ (B @ v)); m3 = float(v @ (O_ker(N, ker["f3"]) @ v))
    Var = float((B @ v) @ (B @ v) - q * q)
    Mmax = float(np.linalg.eigvalsh(M)[-1])
    hq = h - q; hm3 = h - m3
    hlo = h_lo_trial(x, KX); l2h = lam2_hi(x, N, KX); glo = hlo - l2h
    # cone captures
    cap = {}
    for nm, wh in [("f2", "f2"), ("f-f2", "f-f2"), ("f-f3", "f-f3"), ("g", "g")]:
        cap[nm] = kL3(v, karr(x, KX, wh), N, KX)
    Mmax_u, Mmax_p = lam_max_dom(M, N)
    return dict(varpi=varpi, x=x, h=h, q=q, m3=m3, E=E, Var=Var, Mmax=Mmax,
                hq=hq, hm3=hm3, glo=glo, cap=cap, Mmax_u=Mmax_u, Mmax_p=Mmax_p)


VARPI = [0.3, 0.5, 0.7, 1.0, 1.3, 1.5, 2.0, 3.0, 6.0, 10.0, 20.0]
print("cone-coercivity capture (kL3/true) for each kernel; and lam_max(M) upper-bound tightness")
print(f"{'w':>5} {'cap f-f2':>9} {'cap f-f3':>9} {'cap f2':>8} {'Mmax_u/M':>9} {'Mmax_p/M':>9}")
for w in VARPI:
    d = analyze(w); c = d["cap"]
    print(f"{w:5.2f} {c['f-f2']/d['hq']:9.4f} {c['f-f3']/d['hm3']:9.4f} {c['f2']/d['q']:8.4f} {d['Mmax_u']/d['Mmax']:9.4f} {d['Mmax_p']/d['Mmax']:9.4f}")

print()
print("FULLY-ANALYTIC closure attempt.  Mmax^hi=Perron(Mdom); cone lower bounds; kappa=0.66E (G2 analytic).")
print("op-norm-minus_an = Mmax_p^2 - (hq_cone)^2 ; qfree_an = glo*(hm3_cone) - glo^2/4")
print("RHS = glo*(q-m3) with (q-m3)>=x*kappa/2, kappa=0.66E.  ratio<1 closes.")
print(f"{'w':>5} {'opm_an/RHS':>10} {'qfree_an/RHS':>12} {'budget(glo,trueE)':>17}")
for w in VARPI:
    d = analyze(w); c = d["cap"]
    x = d["x"]
    hq_lo = c["f-f2"]           # cone lower bound on h-q (uses v; profile-analytic in principle)
    hm3_lo = c["f-f3"]          # cone lower bound on h-m3
    qm3_lo = x * (0.66 * d["E"]) / 2.0
    RHS = d["glo"] * qm3_lo
    opm_an = d["Mmax_p"] ** 2 - hq_lo ** 2
    qfree_an = d["glo"] * hm3_lo - d["glo"] ** 2 / 4.0
    budget = d["glo"] * (d["q"] - d["m3"]) / d["Var"]
    print(f"{w:5.2f} {opm_an/RHS:10.3f} {qfree_an/RHS if qfree_an>0 else -1:12.3f} {budget:17.2f}")
