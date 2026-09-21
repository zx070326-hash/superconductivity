"""MY OWN recomputation of two agent claims (CLAUDE.md: never trust 'done').
(1) eliashberg-xcheck: authors' Tc vs ours, claimed -3.0% systematic, std 0.1 K.
(2) two-scale-two-band: claimed dCOU = 89.2% at w=0.8, freq_ratio=316.
numerical-evidence (float64)."""
import os, math, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_deficit_budget import Tc_P
from h3s_fill_factor import moments, MEV_K, ROOT

def authors_Tc(path):
    T, D = [], []
    for line in open(path):
        p = line.split()
        if len(p) >= 2:
            try: T.append(float(p[0])); D.append(float(p[1]))
            except ValueError: pass
    T, D = np.array(T), np.array(D)
    k = D > 0
    T, D = T[k], D[k]
    last = T[-1]
    # sqrt extrapolation: Delta^2 linear in T near Tc. Use last m points.
    out = {}
    for m in (2, 3, 5):
        t, d2 = T[-m:], D[-m:]**2
        a, b = np.polyfit(t, d2, 1)          # d2 = a T + b
        out[m] = -b/a
    return last, out

print("=== (1) authors' Tc extraction, MY extractor, muc_0.12 (mu*=0.12) ===")
print(f"{'case':>34} {'lastT':>7} {'sqrt2pt':>8} {'sqrt3pt':>8} {'sqrt5pt':>8} | {'ours':>7} {'diff%':>7}")
CASES = [("H3S-200GPa","EPW_anharmonic"), ("D3S-200GPa","EPW_anharmonic"),
         ("H3S-160GPa","EPW_anharmonic"), ("D3S-160GPa","EPW_anh")]
diffs = []
for tag, sub in CASES:
    gp = f"{ROOT}/{tag}/{sub}/muc_0.12/fsr_adia/*.imag_iso_gap0"
    g = glob.glob(gp)
    if not g:
        print(f"{tag+'/'+sub:>34}  MISSING {gp}"); continue
    last, ext = authors_Tc(g[0])
    ga = [x for x in glob.glob(f"{ROOT}/{tag}/{sub}/*.a2f") if not x.endswith(".a2fv")]
    lam, om2, ol, w, a2f = moments(ga[0])
    ours = Tc_P(w, a2f, lam, 0.12, 10.0, 400, math.sqrt(om2))*MEV_K
    d = (ours/ext[3]-1)*100
    diffs.append(d)
    print(f"{tag+'/'+sub:>34} {last:>7.1f} {ext[2]:>8.1f} {ext[3]:>8.1f} {ext[5]:>8.1f} | "
          f"{ours:>7.1f} {d:>+6.1f}%", flush=True)
print(f"  -> mean {np.mean(diffs):+.2f}%  std {np.std(diffs):.2f}%  (agent claimed -3.0%, std 0.1 K)")

print("\n=== (2) two-scale deficit: dCOU at w=0.8, freq_ratio=316 ===")
def two_scale(w, fr, npts=2):
    """P = w*delta(om1) + (1-w)*delta(om2), om2/om1 = fr. Return (w_arr,a2f_arr,lam,sqom2).
    a2F must satisfy P(dom)=(2/lam)a2F/om dom  =>  a2F_i = lam*p_i*om_i/2 for delta weights."""
    lam = 2.0
    om1, om2 = 1.0, float(fr)
    p = np.array([w, 1.0-w]); oms = np.array([om1, om2])
    a2f = lam*p*oms/2.0
    sqom2 = math.sqrt(float(p @ oms**2))
    return oms, a2f, lam, sqom2

for (w, fr) in [(0.8,316.0),(0.99,1e5),(0.5,316.0),(0.8,10.0)]:
    oms, a2f, lam, sq = two_scale(w, fr)
    t0 = Tc_P(oms, a2f, lam, 0.0,  10.0, 400, sq)
    t1 = Tc_P(oms, a2f, lam, 0.13, 10.0, 400, sq)
    vp = sq/(2*math.pi*t1); Nc = int(10*vp)
    print(f"  w={w:<5} fr={fr:<8.0f} lam={lam} sq<om2>={sq:9.2f} | Tc(0)={t0:.5f} Tc(.13)={t1:.5f} "
          f"dCOU={(1-t1/t0)*100:5.1f}%  varpi={vp:8.2f} Nc={Nc}", flush=True)
# Einstein reference at same lambda
oms, a2f, lam, sq = np.array([1.0]), np.array([1.0]), 2.0, 1.0
t0 = Tc_P(oms, a2f, lam, 0.0, 10.0, 400, sq); t1 = Tc_P(oms, a2f, lam, 0.13, 10.0, 400, sq)
print(f"  Einstein  lam=2.0                              | Tc(0)={t0:.5f} Tc(.13)={t1:.5f} "
      f"dCOU={(1-t1/t0)*100:5.1f}%  varpi={1/(2*math.pi*t1):8.2f} Nc={int(10/(2*math.pi*t1))}")
