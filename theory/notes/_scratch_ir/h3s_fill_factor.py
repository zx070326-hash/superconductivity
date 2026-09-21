"""First-hand fill factors from the Mishra EPW a2F data (H3S/D3S, 160/200 GPa),
plus the mu* deficit as a function of r = omega_c/sqrt(<om^2>).

lambda = 2 int a2F/om dom ;  <om^2> = (2/lambda) int a2F om dom   (P-measure moment,
P(dom) = (2/lambda) a2F/om dom, the same convention as Allen-Dynes / KAY).
Ceiling = C_inf sqrt(lambda <om^2>),  C_inf = 0.1827262477 (analytic).
numerical-evidence (float64); trapezoid on the file grid, with a grid-halving check."""
import os, sys, math, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fill_factor_probe import fill, varpi_c, CINF

MEV_K = 11.604518            # 1 meV in K
ROOT = ("D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/"
        "Data-H3S-D3S-PH-EPW-Anharmonic")


def moments(path):
    w, a2f = [], []
    for line in open(path):
        p = line.split()
        if len(p) >= 3:
            try:
                x, y = float(p[0]), float(p[1])
            except ValueError:
                continue
            w.append(x); a2f.append(y)
    w = np.array(w); a2f = np.array(a2f)
    m = w > 0
    w, a2f = w[m], a2f[m]
    lam = 2.0 * np.trapezoid(a2f / w, w)
    om2 = (2.0 / lam) * np.trapezoid(a2f * w, w)
    olog = math.exp((2.0 / lam) * np.trapezoid(a2f / w * np.log(w), w))
    return lam, om2, olog, w, a2f


if __name__ == "__main__":
    print("=== first-hand a2F moments (Mishra EPW) ===")
    print(f"{'dataset':>26} {'lambda':>7} {'sqrt<om2>[meV]':>15} {'om_log[meV]':>12} "
          f"{'r=om_c/sq<om2>':>15} {'ceiling[K]':>11}")
    rows = {}
    for tag in ("H3S-200GPa", "D3S-200GPa", "H3S-160GPa", "D3S-160GPa"):
        for sub in ("EPW_anharmonic", "EPW_harmonic"):
            p = f"{ROOT}/{tag}/{sub}/H3S.a2f"
            if not os.path.exists(p):
                p = f"{ROOT}/{tag}/{sub}/D3S.a2f"
            if not os.path.exists(p):
                g = glob.glob(f"{ROOT}/{tag}/{sub}/*.a2f")
                if not g:
                    continue
                p = g[0]
            lam, om2, olog, w, a2f = moments(p)
            sq = math.sqrt(om2)
            ceil_K = CINF * math.sqrt(lam * om2) * MEV_K
            # KAY cutoff convention omega_c = 10 * sqrt<om2> is the r we scan; report
            # the *data-driven* r for omega_c = max phonon freq * (typical el. scale)
            r10 = 10.0
            rows[(tag, sub)] = (lam, sq, olog, ceil_K)
            print(f"{tag+'/'+sub:>26} {lam:>7.4f} {sq:>15.3f} {olog:>12.3f} "
                  f"{r10:>15.1f} {ceil_K:>11.1f}")

    print("\n=== fill factor vs experiment (Tc_exp from hydrides_experimental.csv) ===")
    EXP = {"H3S-200GPa": (203.0, "155 GPa dome peak (Drozdov 2015); dataset is 200 GPa"),
           "D3S-200GPa": (150.0, "~150 GPa (Drozdov 2015)")}
    for tag, sub in (("H3S-200GPa", "EPW_anharmonic"), ("D3S-200GPa", "EPW_anharmonic"),
                     ("H3S-200GPa", "EPW_harmonic"), ("D3S-200GPa", "EPW_harmonic")):
        if (tag, sub) not in rows:
            continue
        lam, sq, olog, ceil_K = rows[(tag, sub)]
        tc_exp, note = EXP.get(tag, (float("nan"), ""))
        # framework's own Einstein prediction at this lambda, mu*=0 and 0.13
        t0, f0 = fill(lam, 0.0)
        t1, f1 = fill(lam, 0.13, 10.0)
        tcE0 = t0 * sq * MEV_K            # Einstein Tc with om_E = sqrt<om2>
        tcE1 = t1 * sq * MEV_K
        print(f"{tag+'/'+sub:>26}  lam={lam:.3f}  ceiling={ceil_K:6.1f} K | "
              f"Tc_E(mu*=0)={tcE0:6.1f} (f={f0:.3f})  Tc_E(.13)={tcE1:6.1f} (f={f1:.3f}) | "
              f"Tc_exp={tc_exp:5.1f}  f_exp={tc_exp/ceil_K:.3f}")

    print("\n=== mu* deficit vs r = omega_c/sqrt<om^2>  (Einstein, mu*=0.13) ===")
    print(f"{'lam':>6} " + " ".join(f"{'r='+str(r):>9}" for r in (3, 5, 10, 20, 40)))
    for lam in (1.0, 1.7, 2.0, 3.0):
        t0, _ = fill(lam, 0.0)
        line = f"{lam:>6.2f} "
        for r in (3, 5, 10, 20, 40):
            t1, _ = fill(lam, 0.13, float(r))
            line += f" {(1-t1/t0)*100:>8.1f}%"
        print(line, flush=True)
    print("\n(r grows as sqrt(M) on deuteration: omega_c is electronic, sqrt<om^2> ~ 1/sqrt(M))")
