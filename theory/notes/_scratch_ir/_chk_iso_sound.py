"""Soundness checks for the PdH numbers.

(A) The Coulomb step cutoff Nc = floor(r*varpi(T)) makes g(T) := lam*lam1 - 1
    DISCONTINUOUS in T.  Bisection presumes a single crossing.  Here Nc ~ 200
    (vs 13 for H3S), so scan g finely and count sign changes / locate the
    LARGEST root (memory: Tc := sup over down-jumps).
(B) r-scan with mu* RUNNING per Anderson-Morel: 1/mu*(r) = 1/mu*(r0) - ln(r/r0),
    anchored mu*(r=10) = 0.085.  Fixed-mu* r-scans are not self-consistent.
(C) exact three-channel decomposition of alpha.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from _chk_iso_core import C_INF, frak_f, Tc_einstein, assemble_vec
from numpy.linalg import eigvalsh

M = {"H": 1.00782503, "D": 2.01410178, "T": 3.01604928}
LAM = {"H": 0.40, "D": 0.46, "T": 0.48}
WLOG = {"H": 405.0, "D": 304.0, "T": 257.0}
CM1_TO_MEV = 0.1239841984
MEV_TO_K = 11.604518121
MUS = 0.085


def g_of_T(lam, mus, T, N=800, r=10.0):
    vp = 1.0 / (2.0 * np.pi * T)
    k = np.arange(2 * N + 1, dtype=float)
    F = vp * vp / (k * k + vp * vp)
    Mx = assemble_vec(F, N)
    Nc = int(np.floor(r * vp))
    if mus > 0 and Nc > 0:
        u = 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)
        uc = np.where(np.arange(N) < min(Nc, N), u, 0.0)
        Mx = Mx - 2.0 * (mus / lam) * np.outer(uc, uc)
    return lam * eigvalsh(Mx)[-1] - 1.0, Nc


print("=" * 80)
print("(A) root multiplicity of g(T) around Tc  [lam=0.40 / 0.46, mu*=0.085, r=10]")
print("=" * 80)
for tag, lam in (("PdH", 0.40), ("PdD", 0.46)):
    T0 = Tc_einstein(lam, MUS, N=800)
    Ts = T0 * np.exp(np.linspace(-0.35, 0.35, 481))
    gs, ncs = [], []
    for T in Ts:
        a, b = g_of_T(lam, MUS, T)
        gs.append(a); ncs.append(b)
    gs = np.array(gs)
    sc = np.where(np.sign(gs[:-1]) != np.sign(gs[1:]))[0]
    print(f"{tag}: bisection Tc/wE = {T0:.9f}   Nc at Tc = {ncs[len(ncs)//2]}")
    print(f"     Nc range over the scan window: {min(ncs)} .. {max(ncs)}"
          f"   ({len(set(ncs))} distinct values, i.e. that many jump points)")
    print(f"     sign changes of g on the window: {len(sc)}  at T/Tc = "
          f"{', '.join(f'{Ts[i]/T0:.5f}' for i in sc)}")
    print(f"     g monotone decreasing? {'YES' if np.all(np.diff(gs)<0) else 'NO'}"
          f"   (n_increasing_steps = {int((np.diff(gs)>0).sum())} / {len(gs)-1})")
    print(f"     largest sign change  T/Tc = {Ts[sc[-1]]/T0:.6f}"
          f"   -> relative to bisection root: {100*(Ts[sc[-1]]/T0-1):+.4f}%")

print()
print("=" * 80)
print("(B) r-scan with mu* RUNNING (Anderson-Morel), anchor mu*(r=10)=0.085")
print("=" * 80)


def mus_run(r, r0=10.0, mus0=MUS):
    return 1.0 / (1.0 / mus0 - np.log(r / r0))


lnM = np.log(M["D"] / M["H"])
a_ceil_proxy = 0.5 * np.log((LAM["H"] * WLOG["H"] ** 2) /
                            (LAM["D"] * WLOG["D"] ** 2)) / lnM
print(f"{'r':>6}{'mu*(r)':>9}{'N':>6}{'f_H':>9}{'f_D':>9}"
      f"{'a_deficit':>11}{'a_budget':>10}{'Tc_H K':>9}")
for r, N in [(5.0, 800), (10.0, 800), (20.0, 1600), (40.0, 2400), (80.0, 3200)]:
    mu = mus_run(r)
    fH = frak_f(LAM["H"], mu, N=N, r=r)
    fD = frak_f(LAM["D"], mu, N=N, r=r)
    ad = -np.log(fD / fH) / lnM
    ceilH = C_INF * np.sqrt(LAM["H"]) * WLOG["H"] * CM1_TO_MEV * MEV_TO_K
    print(f"{r:6.1f}{mu:9.4f}{N:6d}{fH:9.5f}{fD:9.5f}{ad:+11.4f}"
          f"{a_ceil_proxy+ad:+10.4f}{fH*ceilH:9.3f}")

print()
print("=" * 80)
print("(C) exact channel decomposition of alpha(PdH->PdD), Einstein proxy, r=10")
print("=" * 80)
fH0, fD0 = frak_f(LAM["H"], 0.0, N=800), frak_f(LAM["D"], 0.0, N=800)
fH, fD = frak_f(LAM["H"], MUS, N=800), frak_f(LAM["D"], MUS, N=800)
a_lam = 0.5 * np.log(LAM["H"] / LAM["D"]) / lnM
a_w2 = np.log(WLOG["H"] / WLOG["D"]) / lnM
a_fin = -np.log(fD0 / fH0) / lnM
a_cou = -np.log((fD / fD0) / (fH / fH0)) / lnM
print(f"  a_ceiling  = a_lambda + a_<w2>        = {a_lam:+.4f} {a_w2:+.4f}"
      f" = {a_lam+a_w2:+.4f}")
print(f"  a_deficit  = a_FIN    + a_COU         = {a_fin:+.4f} {a_cou:+.4f}"
      f" = {a_fin+a_cou:+.4f}")
print(f"  a_budget                              = {a_lam+a_w2+a_fin+a_cou:+.4f}")
print(f"  (check vs direct)                     = "
      f"{-np.log((fD*np.sqrt(LAM['D'])*WLOG['D'])/(fH*np.sqrt(LAM['H'])*WLOG['H']))/lnM:+.4f}")
print()
print(f"  ln(M_D/M_H) = {lnM:.6f}  -> amplification: 1% error in the Tc RATIO"
      f" = {0.01/lnM:.4f} in alpha")
print(f"  to land alpha within +/-0.03 you need the Tc ratio to"
      f" {0.03*lnM*100:.2f}%")
