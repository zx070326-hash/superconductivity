"""
Round 4: the shape of H_mu and the mu~ threshold for the Tolmachev hump.

 (I)  H_mu(varpi) at fixed N_c, plotted in the scaled variable s = varpi/N_c,
      for several mutilde  -> where does the N-shape appear?
 (J)  mutilde threshold for existence of the interior local maximum.

float64 => numerical-evidence.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, importlib.util, pathlib, sys

_here = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "chk1", _here / "_chk_hmu_concavity_2026-09-18.py")
chk1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(chk1)
OpFFT, f_mix, top_eigs = chk1.OpFFT, chk1.f_mix, chk1.top_eigs


def h_mu(varpi, mt, Nc, Nfac=20):
    N = max(400, 2 * Nc + 10, int(Nfac * varpi))
    kk = np.arange(2 * N, dtype=float)
    op = OpFFT(f_mix(kk, [1.0], [varpi ** 2]), N, mt=mt, Nc=Nc)
    return top_eigs(op, k=1)[0][0]


def has_hump(mt, Nc, npts=34, smax=4.0):
    """interior local max of varpi -> h_mu on (0.01*Nc, smax*Nc)?  returns
    (bool, s_at_max, H_max, s_at_min, H_min)."""
    s = np.exp(np.linspace(np.log(0.01), np.log(smax), npts))
    vs = s * Nc
    H = np.array([h_mu(v, mt, Nc) for v in vs])
    for i in range(1, npts - 1):
        if H[i] > H[i - 1] and H[i] > H[i + 1]:
            j = i
            while j + 1 < npts and H[j + 1] <= H[j]:
                j += 1
            return True, s[i], H[i], s[j - 1] if j > i else s[i], H[j - 1] if j > i else H[i]
    return False, np.nan, np.nan, np.nan, np.nan


def probe_shape():
    print("=" * 78)
    print("(I) shape of H_mu in the scaled variable s = varpi / N_c   (N_c = 200)")
    print("=" * 78)
    Nc = 200
    svals = np.array([0.02, 0.05, 0.1, 0.15, 0.25, 0.4, 0.5, 0.7, 1.0, 1.5,
                      2.0, 3.0, 5.0])
    mts = [0.0, 0.2, 0.4, 0.6, 1.0]
    print(f"  {'s=varpi/N_c':>12} " + " ".join(f"{('mt='+str(m)):>11}" for m in mts))
    tab = {}
    for s in svals:
        row = [h_mu(s * Nc, m, Nc) for m in mts]
        tab[s] = row
        print(f"  {s:12.3f} " + " ".join(f"{r:11.6f}" for r in row))
    print("\n  column-wise: mu*=0 is monotone increasing in s (h ~ ln varpi);")
    print("  large mutilde develops a local max then a dip (the Tolmachev crossover).")
    print("  Coulomb cost h(mu*=0) - h_mu :")
    print(f"  {'s':>12} " + " ".join(f"{('mt='+str(m)):>11}" for m in mts[1:]))
    for s in svals:
        row = tab[s]
        print(f"  {s:12.3f} " + " ".join(f"{row[0]-r:11.6f}" for r in row[1:]))
    print()


def probe_threshold():
    print("=" * 78)
    print("(J) mutilde threshold for the interior local maximum (bisection)")
    print("=" * 78)
    print(f"  {'N_c':>6} {'mt_hump':>10} {'s@max':>8} {'H_max':>9} {'s@min':>8} "
          f"{'H_min':>9}  (at mt = 1.5*mt_hump)")
    for Nc in (25, 50, 100, 200, 400):
        lo, hi = 0.02, 4.0
        if not has_hump(hi, Nc)[0]:
            print(f"  {Nc:6d}   no hump even at mutilde={hi}"); continue
        if has_hump(lo, Nc)[0]:
            print(f"  {Nc:6d}   hump already at mutilde={lo}"); continue
        for _ in range(9):
            mid = np.sqrt(lo * hi)
            if has_hump(mid, Nc)[0]:
                hi = mid
            else:
                lo = mid
        ok, sM, HM, sm, Hm = has_hump(1.5 * hi, Nc)
        print(f"  {Nc:6d} {hi:10.4f} {sM:8.4f} {HM:9.4f} {sm:8.4f} {Hm:9.4f}")
    print("\n  mutilde = mu*/lambda ;  a hump needs LARGE mutilde = WEAK coupling.")
    print()


if __name__ == "__main__":
    probe_shape()
    probe_threshold()
