"""直接对撞：它的 KT.lam_k 与我的 lam1，同 (lams,ws,mus,Nc,N)。
定位 1e-6 差是约定/实现差，还是我的重建错。"""
import os, importlib.util; os.environ["OMP_NUM_THREADS"]="1"
import numpy as np
from scipy.linalg import eigh
_d = os.path.join("theory","notes","_scratch_ir")
def _load(n):
    s = importlib.util.spec_from_file_location(n, os.path.join(_d, n+".py"))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
KT = _load("_chk_main_kato_temple_two_atom_prototype_2026-09-23")

def mine(lam, w, mus, Nc, nd):
    n = np.arange(nd); u = 1.0/np.sqrt(2.0*n+1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    k = np.arange(1, 2*nd+2); g = np.zeros(2*nd+2)
    g[1:] = lam*w*w/(w*w + k*k)
    off = g[np.abs(i-j)] + g[i+j+1]; np.fill_diagonal(off, 0.0)
    cum = np.concatenate([[0.0], np.cumsum(g[1:nd+1])])
    M = np.outer(u,u)*(off + np.diag(g[2*n+1] - 2.0*cum[:nd]))
    if mus>0 and Nc>0:
        uc = np.where(n<Nc,u,0.0); M = M - 2.0*mus*np.outer(uc,uc)
    return float(eigh(M, subset_by_index=[nd-1,nd-1], eigvals_only=True)[0])

print(f"{'lam':>8}{'w_lo':>11}{'Nc':>4}{'N':>6}{'KT.lam_k':>15}{'我的':>15}{'差':>12}")
for lam, w, ml in [(2.0,1.100033,1.021190594),(3.2953,0.731708,1.020081705),
                   (4.0,0.640491,1.021883712),(6.0,0.495741,1.025180308)]:
    for Nc in sorted({int(np.floor(10*w)), int(np.floor(10*w))-1}):
        for nd in (600, 1200):
            a = KT.lam_k([lam],[w],0.20,Nc,nd); b = mine(lam,w,0.20,Nc,nd)
            flag = ""
            if abs(a-ml) < 5e-9: flag = "  <== 与 mu_lo 相符!"
            print(f"{lam:8}{w:11.6f}{Nc:4d}{nd:6d}{a:15.9f}{b:15.9f}{a-b:+12.2e}{flag}")
    print(f"{'':>8}{'':>11}{'':>4}{'':>6}{'mu_lo=':>15}{ml:15.9f}")
