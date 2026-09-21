import os, sys, time
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "derivations"))
from p4_interchange_uniformity import H_mat
from p4d_coulomb_uniqueness import uc_vec, thr_from, top_pencil, O_f

N = 600

def H_multi(NN, s, comps, order=0):
    out = np.zeros((NN, NN))
    for p, c in comps:
        out += p * H_mat(NN, c * s, order)
    return out

def O_g_multi(NN, s, comps):
    """O[g](s) = s dA_P/ds = sum_i p_i * O_f(c_i*s, NN).  Verified by finite
    difference: max|O_g_multi - s dA/ds| = 2.7e-11 (see probe_sign_check3.py)."""
    out = np.zeros((NN, NN))
    for p, c in comps:
        out += p * O_f(c * s, NN)
    return out

def six_constants_uniq(NN, s, comps, Nc):
    AP = H_multi(NN, s, comps, 0)
    Og = O_g_multi(NN, s, comps)
    e, V = np.linalg.eigh(AP)
    mu1, mu2 = float(e[-1]), float(e[-2])
    gam = mu1 - mu2
    uc = uc_vec(Nc, NN)
    c = V.T @ uc
    c1 = abs(float(c[-1]))
    up = float(np.linalg.norm(c[:-1]))
    G = V.T @ Og @ V
    q1 = float(G[-1,-1])
    b = float(np.linalg.norm(G[:-1,-1]))
    m = max(0.0, -float(np.linalg.eigvalsh(Og)[0]))
    thr = thr_from(gam, c1, up, q1, b, m)
    return dict(gam=gam, c1=c1, up=up, q1=q1, b=b, m=m, thr=thr, AP=AP, Og=Og, mu1=mu1, mu2=mu2)

comps = [(0.8, 0.3), (0.2, 4.0)]
ss = np.arange(0.05, 3.0 + 1e-9, 0.02)
Nc_list = [5, 8, 12]
mut_list = [0.13, 0.26, 0.33]

t0 = time.time()
lam = {(Nc, mt): np.empty(ss.size) for Nc in Nc_list for mt in mut_list}
uc_cache = {Nc: uc_vec(Nc, N) for Nc in Nc_list}
for i, s in enumerate(ss):
    AP = H_multi(N, s, comps, 0)
    for Nc in Nc_list:
        for mt in mut_list:
            lam1, _, _ = top_pencil(AP, uc_cache[Nc], mt)
            lam[(Nc, mt)][i] = lam1
print(f"scan done in {time.time()-t0:.1f}s")

n_violations = 0
n_correctly_excluded = 0
n_wrongly_included = 0
report = []
for Nc in Nc_list:
    for mt in mut_list:
        arr = lam[(Nc, mt)]
        d = np.diff(arr)
        neg_idx = np.where(d < -1e-10)[0]
        if len(neg_idx) == 0:
            continue
        n_violations += 1
        j = neg_idx[np.argmin(d[neg_idx])]
        s_at = ss[j]
        C = six_constants_uniq(N, s_at, comps, Nc)
        excluded = mt >= C['thr']
        report.append((Nc, mt, s_at, C['thr'], d[j], excluded))
        if excluded:
            n_correctly_excluded += 1
        else:
            n_wrongly_included += 1

print(f"\n{len(report)} non-monotone (Nc,mut) series found; checking threshold exclusion:\n")
print(f"{'Nc':>4} {'mut':>6} {'s_worst':>8} {'thr_U':>10} {'d(worst)':>12} {'excluded?':>10}")
for Nc, mt, s_at, thr, dworst, excl in report:
    print(f"{Nc:4d} {mt:6.3f} {s_at:8.3f} {thr:10.5f} {dworst:12.3e} {str(excl):>10}")

print(f"\nnon-monotone series: {n_violations}; correctly excluded by mutilde_*(P,T,Nc): {n_correctly_excluded}; WRONGLY INCLUDED: {n_wrongly_included}")
