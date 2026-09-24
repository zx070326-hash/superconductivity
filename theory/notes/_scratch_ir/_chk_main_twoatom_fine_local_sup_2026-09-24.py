"""细扫两原子族（固定 lam, M2=lam*w^2）看 sup 是否超过 Einstein 约 1e-6。
粗扫 N=200 定位，最优点用 N=600 复核。判据：sup - Einstein 是否 >= mu_lo - Einstein。"""
import os, importlib.util; os.environ["OMP_NUM_THREADS"]="1"
import numpy as np
_d = os.path.join("theory","notes","_scratch_ir")
def _load(n):
    s = importlib.util.spec_from_file_location(n, os.path.join(_d, n+".py"))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
KT = _load("_chk_main_kato_temple_two_atom_prototype_2026-09-23")

def f(lam, t, a, w, mus, Nc, nd):
    if t >= 1.0: return KT.lam_k([lam],[w],mus,Nc,nd)
    w2sq = (w*w - t*a*a)/(1.0-t)
    if w2sq <= 1e-14: return -9.0
    return KT.lam_k([lam*t, lam*(1-t)],[a, np.sqrt(w2sq)], mus, Nc, nd)

CASES = [(2.0,1.100033,0.20,1.021190594),(3.2953,0.731708,0.20,1.020081705),
         (4.0,0.640491,0.20,1.021883712),(3.2953,0.668045,0.13,1.011682278)]
for lam, w, mus, ml in CASES:
    Nc = int(np.floor(10*w))
    ein = KT.lam_k([lam],[w],mus,Nc,600)
    best = (-9.0, None, None)
    for t in np.concatenate([1-np.logspace(-4,-0.05,26), np.linspace(0.05,0.9,18)]):
        amax = w/np.sqrt(t)
        for a in np.concatenate([w*np.linspace(0.3,0.999,16),
                                 np.linspace(w*1.0001, amax*0.99999, 16)]):
            if a >= amax: continue
            v = f(lam, t, a, w, mus, Nc, 200)
            if v > best[0]: best = (v, t, a)
    _, t, a = best
    # 局部细化（N=200），再用 N=600 复核
    for scale in (1e-1, 1e-2, 1e-3):
        for tt in t + t*scale*np.linspace(-1,1,9):
            for aa in a + a*scale*np.linspace(-1,1,9):
                if not (0 < tt < 1) or aa <= 0: continue
                v = f(lam, tt, aa, w, mus, Nc, 200)
                if v > best[0]: best = (v, tt, aa)
        _, t, a = best
    sup600 = f(lam, t, a, w, mus, Nc, 600)
    gain, need = sup600 - ein, ml - ein
    print(f"lam={lam:<7} mu*={mus} w={w:.6f} Nc={Nc}")
    print(f"   Einstein(N=600)={ein:.9f}   最优两原子(N=600)={sup600:.9f}  (t={t:.5f}, w1/w={a/w:.4f})")
    print(f"   族内增益 {gain:+.2e}   mu_lo 需要的增益 {need:+.2e}   "
          f"{'==> 族内可达, mu_lo 合法' if gain >= need else ('==> mu_lo 无需增益(已<=Einstein)' if need<=0 else '** 族内未达, 仍缺 %.2e **'%(need-gain))}")
