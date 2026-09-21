"""Diagnostic for the team-lead discrepancy on alpha_COU.

Hypothesis: the peer put mu~ = mu* = 0.085 straight into
    M = O[F] - 2 mu~ uc uc^T
whereas the Tc criterion is  lam * lam1( O[F] - 2 (mu*/lam) uc uc^T ) = 1,
i.e. mu~ = mu*/lam  (= 0.2125 / 0.1848 / 0.1771 for PdH/PdD/PdT).

Run BOTH conventions through the identical machinery and compare against the
peer's reported Tc_H/Tc_D/Tc_T = 8.86 / 9.59 / 8.98 K and alpha = -0.115.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.linalg import eigvalsh
from _chk_iso_core import C_INF, assemble_vec, f_einstein, Tc_einstein

M = {"H": 1.00782503, "D": 2.01410178, "T": 3.01604928}
LAM = {"H": 0.40, "D": 0.46, "T": 0.48}
WLOG = {"H": 405.0, "D": 304.0, "T": 257.0}
CM1_TO_MEV = 0.1239841984
MEV_TO_K = 11.604518121
MUS = 0.085


def Tc_conv(lam, mus, N=900, r=10.0, divide_by_lam=True, tol=1e-13):
    """Tc/wE.  divide_by_lam=True -> mu~ = mu*/lam (the criterion).
                         =False -> mu~ = mu*      (the suspected peer reading)."""
    mut = (mus / lam) if divide_by_lam else mus

    def g(T):
        vp = 1.0 / (2.0 * np.pi * T)
        Mx = assemble_vec(f_einstein(vp, N), N)
        Nc = int(np.floor(r * vp))
        assert Nc < N, f"Nc={Nc} >= N={N}"
        if mut > 0 and Nc > 0:
            u = 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)
            uc = np.where(np.arange(N) < Nc, u, 0.0)
            Mx = Mx - 2.0 * mut * np.outer(uc, uc)
        return lam * eigvalsh(Mx)[-1] - 1.0

    lo, hi = 2.2e-3, 0.2      # lo chosen so that Nc = floor(r/(2 pi lo)) < N
    assert g(lo) > 0.0 > g(hi), (g(lo), g(hi))
    while hi / lo > 1.0 + tol:
        mid = np.sqrt(lo * hi)
        if g(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


def alpha(TcA, TcB, A, B):
    return -np.log(TcB / TcA) / np.log(M[B] / M[A])


print("=" * 84)
print("mu~ convention diagnostic   (N=900, r=10, mu*=0.085 fixed, w_log proxy)")
print("=" * 84)
res = {}
for tag, div in (("mu~ = mu*/lam  (criterion)", True),
                 ("mu~ = mu*      (suspected)", False)):
    print(f"\n--- {tag} ---")
    print(f"{'iso':>5}{'lam':>7}{'mu~':>9}{'Tc/wE':>11}{'Tc K':>9}{'Nc':>6}")
    tc = {}
    for k in "HDT":
        t = Tc_conv(LAM[k], MUS, divide_by_lam=div)
        wE = WLOG[k] * CM1_TO_MEV
        tcK = t * wE * MEV_TO_K
        tc[k] = tcK
        mut = (MUS / LAM[k]) if div else MUS
        print(f"{'Pd'+k:>5}{LAM[k]:7.2f}{mut:9.4f}{t:11.6f}{tcK:9.3f}"
              f"{int(np.floor(10.0/(2*np.pi*t))):6d}")
    res[div] = tc
    print(f"   alpha(H->D) = {alpha(tc['H'],tc['D'],'H','D'):+.4f}"
          f"   alpha(D->T) = {alpha(tc['D'],tc['T'],'D','T'):+.4f}")

print("\n" + "=" * 84)
print("peer reported: Tc = 8.86 / 9.59 / 8.98 K,  alpha(H->D) = -0.115")
print("mine, mu~=mu*  :  Tc = "
      + " / ".join(f"{res[False][k]:.2f}" for k in "HDT")
      + f" K,  alpha(H->D) = {alpha(res[False]['H'],res[False]['D'],'H','D'):+.4f}")
print("mine, mu~=mu*/lam: Tc = "
      + " / ".join(f"{res[True][k]:.2f}" for k in "HDT")
      + f" K,  alpha(H->D) = {alpha(res[True]['H'],res[True]['D'],'H','D'):+.4f}")

# alpha_COU under both conventions, same definition as _chk_iso_sound.py (C)
print("\n" + "=" * 84)
print("alpha_COU := alpha_budget(mu*) - alpha_budget(mu*=0)   [definition used in my note]")
lnMD = np.log(M["D"] / M["H"])
t0 = {k: Tc_conv(LAM[k], 0.0) for k in "HDT"}
ceil = {k: C_INF * np.sqrt(LAM[k]) * WLOG[k] * CM1_TO_MEV * MEV_TO_K for k in "HDT"}
a0 = alpha(t0["H"] * ceil["H"] / (C_INF * np.sqrt(LAM["H"]) * WLOG["H"] * CM1_TO_MEV * MEV_TO_K)
           * ceil["H"] / ceil["H"], 0, "H", "D") if False else None
tc0K = {k: t0[k] * WLOG[k] * CM1_TO_MEV * MEV_TO_K for k in "HDT"}
ab0 = alpha(tc0K["H"], tc0K["D"], "H", "D")
print(f"  alpha_budget(mu*=0)          = {ab0:+.4f}")
for div, tag in ((True, "mu~=mu*/lam"), (False, "mu~=mu*    ")):
    ab = alpha(res[div]["H"], res[div]["D"], "H", "D")
    print(f"  {tag}: alpha_budget = {ab:+.4f}   alpha_COU = {ab-ab0:+.4f}")
