"""Same external dps=120 enclosure audit, applied to v6 (see _chk_v4_enclosure.py).
Plus the F3 sweep: does v6's u_iv now enclose 1/sqrt(2n+1) for every n?
"""
import os, sys, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_deficit_cert_v6 import build_M_iv
from h3s_fill_factor import moments, ROOT

DPS_HI, DPS, N, T = 120, 50, 40, 1.0
g = [x for x in glob.glob(f"{ROOT}/H3S-200GPa/EPW_anharmonic/*.a2f") if not x.endswith(".a2fv")]
lam_f, om2_f, olog_f, w, a2f = moments(g[0])

mp.mp.dps = DPS_HI
W = [mp.mpf(float(x)) for x in w]      # EXACT binary, not repr()
Y = [mp.mpf(float(x)) for x in a2f]
n_pts = len(W)
tw = [None]*n_pts
tw[0] = (W[1]-W[0])/2; tw[-1] = (W[-1]-W[-2])/2
for i in range(1, n_pts-1): tw[i] = (W[i+1]-W[i-1])/2
raw = [(Y[i]/W[i])*tw[i] for i in range(n_pts)]
tot = mp.fsum(raw)
p = [raw[i]/tot for i in range(n_pts)]
xs = [(W[i]/(2*mp.pi*mp.mpf(float(T))))**2 for i in range(n_pts)]
def F_hi(k):
    if k == 0: return mp.mpf(1)
    k2 = mp.mpf(k)**2
    return mp.fsum([p[i]*xs[i]/(k2+xs[i]) for i in range(n_pts)])
Fh = [F_hi(k) for k in range(2*N+2)]
u_hi = [mp.mpf(1)/mp.sqrt(mp.mpf(2*n+1)) for n in range(N)]
L1_hi = [mp.mpf(0)]*N
acc = mp.mpf(0)
for nn in range(1, N):
    acc += 2*Fh[nn]; L1_hi[nn] = acc/mp.mpf(2*nn+1)
def A_hi(i, j):
    if i == j: return -L1_hi[i] + Fh[2*i+1]*u_hi[i]**2
    return (Fh[abs(i-j)] + Fh[i+j+1])*u_hi[i]*u_hi[j]

M_iv, u_iv, fker_iv, L1_iv, _, _ = build_M_iv(w, a2f, lam_f, 0.0, 10.0, T, N, DPS,
                                              Nc_fixed=None, poison=False)
def relw(iv):
    mp.mp.dps = DPS_HI
    a, b = mp.mpf(iv.a), mp.mpf(iv.b)
    m = (a+b)/2
    return float((b-a)/abs(m)) if m != 0 else float(b-a)

rows = ([('F', k, fker_iv[k], Fh[k]) for k in (1, 20, 80)]
        + [('L1', k, L1_iv[k], L1_hi[k]) for k in (1, 20, 39)]
        + [('A', ij, None, A_hi(*ij)) for ij in ((0,0), (0,1), (20,20), (20,21))])
print(f"{'quantity':<12}{'relwidth':>11}{'encl TRUE':>12}")
print("-"*35)
nfail = 0
for nm, key, iv, truth in rows:
    if iv is None:
        i, j = key
        iv = (-L1_iv[i] + fker_iv[2*i+1]*u_iv[i]**2 if i == j else
              (fker_iv[abs(i-j)] + fker_iv[i+j+1])*u_iv[i]*u_iv[j])
    mp.mp.dps = DPS_HI
    ok = mp.mpf(iv.a) <= truth <= mp.mpf(iv.b)
    nfail += (not ok)
    print(f"{nm+str(key):<12}{relw(iv):11.2e}{'YES' if ok else '*** NO ***':>12}")
print("-"*35)
print(f"v6 non-enclosing: {nfail}/10")

# F3 sweep on v6's ACTUAL u_iv (build_M_iv output), not a re-typed idiom --
# retyping the idiom would test my transcription, not v6's source.
M2_iv, u800_iv, _, _, _, _ = build_M_iv(w, a2f, lam_f, 0.0, 10.0, T, 800, DPS,
                                        Nc_fixed=None, poison=False)
bad = 0
first = []
for n in range(800):
    mp.mp.dps = 200
    t = mp.mpf(1)/mp.sqrt(mp.mpf(2*n+1))
    ok = mp.mpf(u800_iv[n].a) <= t <= mp.mpf(u800_iv[n].b)
    if not ok:
        bad += 1
        if len(first) < 6: first.append(n)
print(f"F3 sweep on v6's OWN u_iv, n=0..799: {bad}/800 non-enclosing "
      f"(v4's mp.sqrt gave 232/800)" + (f"  first: {first}" if first else ""))
