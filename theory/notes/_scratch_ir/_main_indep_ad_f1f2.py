"""Does Allen-Dynes' f1*f2 prefactor close the residual drift in D_obs/D_AD?

Main session, from zero.  D_obs on the FULL 28-point grid is recomputed here with
the independent Einstein implementation (no agent numbers imported).

Allen-Dynes 1975 (PRB 12, 905), eq. (34):
    Tc = (f1 f2 / 1.20) w_log exp[ -1.04(1+lam) / (lam - mus(1+0.62 lam)) ]
    f1 = [1 + (lam/L1)^{3/2}]^{1/3},              L1 = 2.46 (1 + 3.8 mus)
    f2 = 1 + (w2bar/w_log - 1) lam^2/(lam^2+L2^2), L2 = 1.82 (1+6.3 mus)(w2bar/w_log)

KEY: on an Einstein spectrum w_log = w2bar = w_E, so w2bar/w_log = 1 and
     ** f2 == 1 identically **.  Only f1 survives -- and it does NOT cancel in
     D := 1 - Tc(mus)/Tc(0) because L1 carries mus.  So
     D_AD^{f1} = 1 - [f1(lam,mus)/f1(lam,0)] * exp(dE),   dE as in the bare case.
The 1/1.20 and w_log prefactors cancel in the ratio; that is why the bare
comparison was meaningful at all.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh

# ---------------------------------------------------------------- D_obs (ours)
def lam1(lam, mus, T, N=400, r=10.0):
    vp = 1.0 / (2.0 * np.pi * T)
    k = np.arange(2 * N + 1, dtype=float)
    F = vp**2 / (k**2 + vp**2)
    n = np.arange(N); u = 1.0 / np.sqrt(2.0 * n + 1.0)
    L1 = np.zeros(N); acc = 0.0
    for nn in range(1, N):
        acc += 2.0 * F[nn]                 # sum_{k=1}^{nn}, INCLUDES k=nn
        L1[nn] = acc / (2.0 * nn + 1.0)
    i, j = np.meshgrid(n, n, indexing='ij')
    M = (F[np.abs(i - j)] + F[i + j + 1]) * u[i] * u[j] - np.diag(L1)
    M[n, n] -= F[0] * u**2                 # L2 is strictly off-diagonal
    Nc = int(np.floor(r * vp))
    if mus > 0 and Nc > 0:
        uc = np.where(n < Nc, u, 0.0)
        M -= 2.0 * (mus / lam) * np.outer(uc, uc)
    return eigh(M, eigvals_only=True)[-1]

def Tc(lam, mus, N=400):
    f = lambda T: lam * lam1(lam, mus, T, N) - 1.0
    lo, hi = 1e-6, 5.0
    assert f(lo) > 0 > f(hi), (f(lo), f(hi))
    for _ in range(120):
        mid = np.sqrt(lo * hi)
        if f(mid) > 0: lo = mid
        else: hi = mid
    return np.sqrt(lo * hi)

# ------------------------------------------------------------------ AD 1975
def ad_exponent(lam, mus, enh=True):
    eff = mus * (1.0 + 0.62 * lam) if enh else mus
    return -1.04 * (1.0 + lam) / (lam - eff)

def f1(lam, mus):
    L1 = 2.46 * (1.0 + 3.8 * mus)
    return (1.0 + (lam / L1) ** 1.5) ** (1.0 / 3.0)

def f2(lam, mus, ratio=1.0):            # ratio = w2bar/w_log; Einstein => 1
    L2 = 1.82 * (1.0 + 6.3 * mus) * ratio
    return 1.0 + (ratio - 1.0) * lam**2 / (lam**2 + L2**2)

def D_ad(lam, mus, use_f1, enh=True, use_f2=True):
    pref = 1.0
    if use_f1: pref *= f1(lam, mus) / f1(lam, 0.0)
    if use_f2: pref *= f2(lam, mus) / f2(lam, 0.0)      # == 1 on Einstein
    return 1.0 - pref * np.exp(ad_exponent(lam, mus, enh)
                               - ad_exponent(lam, 0.0, enh))

def D_bcs(lam, mus):
    return 1.0 - np.exp(-1.0 / (lam - mus) + 1.0 / lam)

LAMS = [0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
MUS  = [0.05, 0.10, 0.13, 0.20]

print("sanity: f2 on Einstein (ratio=1) must be exactly 1 for every (lam,mus)")
bad = max(abs(f2(l, m) - 1.0) for l in LAMS for m in MUS)
print(f"  max|f2 - 1| = {bad:.3e}   ->  f2 drops out; only f1 can matter\n")

t0 = {l: Tc(l, 0.0) for l in LAMS}
rows = []
print(f"{'lam':>5}{'mus':>6}{'Tc(0)':>10}{'Tc(mus)':>10}{'D_obs%':>9}"
      f"{'/bare':>8}{'/bare+f1':>10}{'f1ratio':>9}"
      f"{'/noenh':>9}{'/noenh+f1':>11}{'/BCS':>9}")
for l in LAMS:
    for m in MUS:
        t1 = Tc(l, m); D = 1.0 - t1 / t0[l]
        rb = D / D_ad(l, m, use_f1=False)
        rf = D / D_ad(l, m, use_f1=True)
        rn = D / D_ad(l, m, use_f1=False, enh=False)
        rnf = D / D_ad(l, m, use_f1=True, enh=False)
        rows.append((l, m, D, rb, rf, rn, rnf, D / D_bcs(l, m)))
        print(f"{l:5.1f}{m:6.2f}{t0[l]:10.6f}{t1:10.6f}{100*D:9.3f}"
              f"{rb:8.3f}{rf:10.3f}{f1(l,m)/f1(l,0.0):9.4f}"
              f"{rn:9.3f}{rnf:11.3f}{D/D_bcs(l,m):9.2f}")

a = np.array([[r[3], r[4], r[5], r[6], r[7]] for r in rows])
print(f"\n{'baseline':<26}{'mean':>8}{'sd':>8}{'sd/mean':>9}"
      f"{'min':>8}{'max':>8}{'max/min':>9}")
for k, nm in enumerate(["AD bare exponent", "AD + f1 (f2=1)",
                        "AD no (1+0.62lam)", "AD no (1+0.62lam), +f1", "BCS"]):
    c = a[:, k]
    print(f"{nm:<26}{c.mean():8.3f}{c.std(ddof=1):8.3f}"
          f"{c.std(ddof=1)/c.mean():9.3f}{c.min():8.3f}{c.max():8.3f}"
          f"{c.max()/c.min():9.2f}")

print("\ndrift along lam at mus=0.13 (the thing f1 was supposed to fix):")
sel = [r for r in rows if abs(r[1] - 0.13) < 1e-12]
for lbl, idx in [("lam      ", 0), ("bare     ", 3), ("bare+f1  ", 4),
                 ("noenh    ", 5), ("noenh+f1 ", 6)]:
    fmt = "{:9.1f}" if idx == 0 else "{:9.3f}"
    print("  " + lbl + ":" + "".join(fmt.format(r[idx]) for r in sel))
