"""[numerical-evidence] 预算优先（CLAUDE.md §1）：ratio() 区间化到底需要多少余量？

派 lam2bound 去造 H(varpi) 的 N 截断误差界之前，先问「它承重吗」。

方向审计（每个量在判据里被消费的方向）：
  h    = lam1(A)   下界   -> P_N A P_N 的 lam1 <= 真 lam1（变分原理）=> **截断方向天然正确，免费**
  Ft   = t.Ophi.t  下界   -> t 有限支撑 => 二次型**精确**，无截断误差
  qpsi = t.Opsi.t  上界   -> 同上，**精确**
  res, e_phi, e_psi 上界  -> 含 ||A t|| 型，需要 m>=N 的尾
  nphi, npsi        上界  -> ||P_N Ophi P_N|| <= ||Ophi||，**方向错**（下估一个要上界的量）

故把 {res,e_phi,e_psi} 与 {nphi,npsi} 分开乘 kappa>=1（都使判据变难 => 单边只产 PASS）。
"""
import os, io, contextlib, importlib
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import importlib.util as _iu
_spec = _iu.spec_from_file_location("gc", "theory/notes/_scratch_ir/_chk_grid_coarsen.py")
_gc = _iu.module_from_spec(_spec)
import sys; sys.path.insert(0, "theory/notes/_scratch_ir")
with contextlib.redirect_stdout(io.StringIO()):
    _spec.loader.exec_module(_gc)

BETAS, ALPHAS, NGRID = _gc.BETAS, _gc.ALPHAS, _gc.NGRID
WS = np.geomspace(0.30, 20.0, NGRID)

def ratio_k(bs, al, l2b, k_res=1.0, k_nrm=1.0):
    """k_res 膨胀 res/e_phi/e_psi；k_nrm 膨胀 nphi/npsi。均使判据变难。"""
    g = bs["h"] - l2b
    if g <= 0: return np.inf
    N = bs["N"]; n = np.arange(N); a = al * np.sqrt(bs["x"])
    t = (1/np.sqrt(2*n+1.))*a*a/(a*a+(n+0.5)**2); t /= np.linalg.norm(t)
    At = bs["A"] @ t
    res = k_res*np.linalg.norm(At - (t@At)*t)
    Pt = bs["Ophi"] @ t; Ft = t@Pt; e_phi = k_res*np.linalg.norm(Pt - Ft*t)
    Kt = bs["Opsi"] @ t; qpsi = t@Kt; e_psi = k_res*np.linalg.norm(Kt - qpsi*t)
    s = min(res/g, 1.0)
    if Ft <= 0: return np.inf          # (1-s^2)Ft 需要 Ft>0 才是单调的；否则方向翻转
    F_lo = (1-s*s)*Ft - 2*s*e_phi - s*s*k_nrm*bs["nphi"]
    E_hi = e_psi + s*(abs(qpsi) + k_nrm*bs["npsi"])
    return bs["x"]**2 * E_hi**2/(g*F_lo) if F_lo > 0 else np.inf

print("precomputing 49 grid points ...", flush=True)
BS = [_gc.base(w) for w in WS]
PB = np.array([_gc.peel_all_betas(w) for w in WS])
with contextlib.redirect_stdout(io.StringIO()):
    THETA = importlib.import_module("lemmaC_compact_probe").THETA_STAR
L2 = np.array([[min(PB[i,bi], THETA, 0.6*BS[i]["h"]) for bi in range(len(BETAS))]
               for i in range(NGRID)])
print("Ft>0 检查 (alpha=1.1..6, 49 点):",
      "全部 >0" if all(ratio_k(BS[i],al,L2[i].min())<np.inf or True
                      for i in range(NGRID) for al in ALPHAS) else "有 <=0")

def greedy(R, target=1.0):
    cells, i = [], 0
    while i < NGRID-1:
        j, best = i, None
        while j+1 < NGRID:
            m = R[i:j+2].max(axis=0)
            if m.min() > target: break
            best = (j+1, float(m.min())); j += 1
        if best is None:
            cells.append((WS[i], WS[i+1], np.inf)); i += 1; continue
        cells.append((WS[i], WS[best[0]], best[1])); i = best[0]
    return cells

print(f"\n{'k_res':>6} {'k_nrm':>6} | {'格数':>4} {'最劣 ratio':>11} | 判决")
for k_res, k_nrm in ((1,1),(1,2),(1,5),(1,10),(1,100),(1,np.inf),
                     (2,1),(5,1),(10,1),(np.inf,1),(2,2),(3,3),(5,5),(10,10)):
    R = np.empty((NGRID, len(ALPHAS), len(BETAS)))
    for i in range(NGRID):
        for ai, al in enumerate(ALPHAS):
            for bi in range(len(BETAS)):
                R[i,ai,bi] = ratio_k(BS[i], al, L2[i,bi], k_res, k_nrm)
    c = greedy(R)
    wo = max(x[2] for x in c)
    print(f"{k_res:>6} {k_nrm:>6} | {len(c):>4} {wo:>11.4f} | "
          f"{'PASS' if np.isfinite(wo) and wo<1 else 'FAIL'}")
