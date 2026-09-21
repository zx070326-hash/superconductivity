"""Two checks prompted by the team-lead's r-scan.

(1) ANCHOR INVARIANCE.  mu* is DEFINED AT a cutoff.  Scanning r at fixed mu*
    varies the cutoff while holding the pseudopotential -> an inconsistent pair.
    Anderson-Morel transport:  1/mu*(r) = 1/mu*(r0) - ln(r/r0).
    CLAIM: once the anchor pair (r0, mu*_0) is fixed, alpha is nearly r-independent,
    and the 17% "r span" is really the span over WHERE mu*=0.085 is anchored.
    Test: for each anchor r0, evaluate alpha at several r and see if it is flat,
    and whether its value equals the fixed-mu* value AT r = r0.

(2) ROUNDING BAND of Errea's 2-significant-figure Tc on the D->T step, computed
    separately for the harmonic row (47/34/30) and the SSCHA row (5.0/6.5/6.9).
    My note quotes +/-0.04 for the harmonic row; +/-0.04 is the SSCHA number.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.linalg import eigvalsh
from _chk_iso_core import C_INF, assemble_vec, f_einstein

M = {"H": 1.00782503, "D": 2.01410178, "T": 3.01604928}
LAM = {"H": 0.40, "D": 0.46, "T": 0.48}
WLOG = {"H": 405.0, "D": 304.0, "T": 257.0}
CM1_TO_MEV = 0.1239841984
MEV_TO_K = 11.604518121


def Tc(lam, mus, N, r, tol=1e-11):
    def g(T):
        vp = 1.0 / (2 * np.pi * T)
        Mx = assemble_vec(f_einstein(vp, N), N)
        Nc = int(np.floor(r * vp))
        assert Nc < N, f"Nc={Nc} >= N={N}  (raise N)"
        if mus > 0 and Nc > 0:
            u = 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)
            uc = np.where(np.arange(N) < Nc, u, 0.0)
            Mx = Mx - 2.0 * (mus / lam) * np.outer(uc, uc)
        return lam * eigvalsh(Mx)[-1] - 1.0
    lo, hi = 4e-3, 3.5e-2
    assert g(lo) > 0 > g(hi), (g(lo), g(hi))
    for _ in range(38):
        m = np.sqrt(lo * hi)
        if g(m) > 0:
            lo = m
        else:
            hi = m
    return np.sqrt(lo * hi)


def alpha_pair(tc, A, B):
    return -np.log(tc[B] / tc[A]) / np.log(M[B] / M[A])


def alpha_at(r, mus, N):
    tc = {k: Tc(LAM[k], mus, N, r) * WLOG[k] * CM1_TO_MEV * MEV_TO_K for k in "HDT"}
    return alpha_pair(tc, "H", "D"), alpha_pair(tc, "D", "T"), tc


def mus_run(r, r0, mus0=0.085):
    return 1.0 / (1.0 / mus0 - np.log(r / r0))


print("=" * 86)
print("(1) IS THE 17% SPAN 'r SENSITIVITY' OR 'ANCHOR SENSITIVITY'?")
print("=" * 86)
print("\n(1a) fixed mu*=0.085, scan r  -- the inconsistent scan (reproduces peer)")
print(f"{'r':>5}{'N':>7}{'Tc_H':>8}{'Tc_D':>8}{'Tc_T':>8}{'a(H->D)':>10}{'a(D->T)':>10}")
fixed = {}
for r, N in ((5., 700), (10., 700), (20., 1000), (40., 1800)):
    aHD, aDT, tc = alpha_at(r, 0.085, N)
    fixed[r] = aHD
    print(f"{r:5.0f}{N:7d}{tc['H']:8.2f}{tc['D']:8.2f}{tc['T']:8.2f}{aHD:+10.4f}{aDT:+10.4f}")

print("\n(1b) Anderson-Morel transport: fix an ANCHOR (r0, mu*=0.085), then scan r.")
print("     If the model is cutoff-consistent, each ROW should be FLAT,")
print("     and should equal the fixed-mu* value at r = r0 (bold column of 1a).")
print(f"{'anchor r0':>10}{'':>3}" + "".join(f"{'r='+str(int(r)):>12}" for r in (5, 10, 20, 40))
      + f"{'span':>9}{'fixed-mu* @r0':>15}")
for r0 in (5., 10., 20.):
    row = []
    for r, N in ((5., 700), (10., 700), (20., 1000), (40., 1800)):
        mu = mus_run(r, r0)
        aHD, _, _ = alpha_at(r, mu, N)
        row.append(aHD)
    print(f"{r0:10.0f}{'':>3}" + "".join(f"{v:+12.4f}" for v in row)
          + f"{max(row)-min(row):9.4f}{fixed[r0]:+15.4f}")

print("\n" + "=" * 86)
print("(2) ROUNDING BAND from Errea's 2-significant-figure Tc, D->T step")
print("=" * 86)
lnDT = np.log(M["T"] / M["D"])
for tag, (tD, tT, half) in {
        "harmonic row  Tc_D=34, Tc_T=30": (34.0, 30.0, 0.5),
        "SSCHA row     Tc_D=6.5, Tc_T=6.9": (6.5, 6.9, 0.05)}.items():
    a_c = -np.log(tT / tD) / lnDT
    a_hi = -np.log((tT - half) / (tD + half)) / lnDT     # most positive alpha
    a_lo = -np.log((tT + half) / (tD - half)) / lnDT     # most negative alpha
    print(f"  {tag}")
    print(f"     alpha_central = {a_c:+.4f}   band = [{a_lo:+.4f}, {a_hi:+.4f}]"
          f"   half-width = {max(a_hi-a_c, a_c-a_lo):.4f}")
print()
print("  harmonic control: my pure-w_log alpha(D->T) = +0.4515 (lam mass-independent)")
a_c = -np.log(30.0 / 34.0) / lnDT
a_hi = -np.log(29.5 / 34.5) / lnDT
print(f"     gap to Errea harm-ME central  = {0.4515 - a_c:.4f}")
print(f"     gap to the TOP of the rounding band = {0.4515 - a_hi:.4f}"
      f"   <- residual that rounding CANNOT absorb")
