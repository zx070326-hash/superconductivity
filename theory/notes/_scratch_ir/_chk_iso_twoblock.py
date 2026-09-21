"""PdH: general-P (two-block) spectrum + r-scan at adequate N.

Two-block model of alpha^2 F, built ONLY from numbers Errea 2013 prints in text:
  - "H-character optical modes have the largest contribution to lambda, between
     75% and 79% of the total depending on the isotope"           (p.4)
  - "The contribution of the low-energy Pd-character acoustic modes is similar
     for the three hydrides."                                      (p.4)
  - "At zone center the PdH optical modes are degenerate at 488 cm-1"  (p.3)
  - TABLE I: lam = 0.40/0.46/0.48, w_log = 405/304/257 cm-1 (SSCHA)

P(dw) = (2/lam)(a2F/w)dw is the probability measure; the lambda-fraction of a
block IS its P-weight.  Model: P = w_ac delta(w_ac) + w_op delta(w_op), with
  * lam_ac := lam * w_ac  held CONSTANT across isotopes ("similar for the three")
  * w_ac(H) fixed by requiring w_op to span the quoted 75-79%
  * w_ac itself (the Pd block frequency) a scan parameter, mass-independent
  * w_op(iso) then DETERMINED by the reported w_log:
        ln w_log = w_ac ln W_ac + w_op ln W_op
Gives <w^2> = w_ac W_ac^2 + w_op W_op^2 -- i.e. the missing second moment.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.linalg import eigvalsh
from _chk_iso_core import C_INF, assemble_vec, frak_f

CM1_TO_MEV = 0.1239841984
MEV_TO_K = 11.604518121
M = {"H": 1.00782503, "D": 2.01410178, "T": 3.01604928}
LAM = {"H": 0.40, "D": 0.46, "T": 0.48}
WLOG = {"H": 405.0, "D": 304.0, "T": 257.0}      # cm^-1, SSCHA
TC_ME = {"H": 5.0, "D": 6.5, "T": 6.9}
MUS = 0.085


def alpha(TcA, TcB, A, B):
    return -np.log(TcB / TcA) / np.log(M[B] / M[A])


# --------------------------------------------------- general-P Tc (my own)
def lam1_P(lam, mus, T, w, om, N=400, r=10.0):
    """w = P-weights (sum 1), om = frequencies in the SAME unit; T in that unit."""
    vp = np.asarray(om, float) / (2.0 * np.pi * T)
    k = np.arange(2 * N + 1, dtype=float)
    F = np.zeros_like(k)
    for wi, v in zip(w, vp):
        F += wi * v * v / (k * k + v * v)
    Mx = assemble_vec(F, N)
    if mus > 0.0:
        w2 = float(np.dot(w, np.asarray(om, float) ** 2))
        vprms = np.sqrt(w2) / (2.0 * np.pi * T)
        Nc = int(np.floor(r * vprms))
        if Nc > 0:
            u = 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)
            uc = np.where(np.arange(N) < min(Nc, N), u, 0.0)
            Mx = Mx - 2.0 * (mus / lam) * np.outer(uc, uc)
    return eigvalsh(Mx)[-1]


def Tc_P(lam, mus, w, om, N=400, r=10.0, tol=1e-12):
    g = lambda T: lam * lam1_P(lam, mus, T, w, om, N, r) - 1.0
    lo, hi = 1e-5 * max(om), 10.0 * max(om)
    assert g(lo) > 0.0 > g(hi), (g(lo), g(hi))
    while hi / lo > 1.0 + tol:
        mid = np.sqrt(lo * hi)
        if g(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


# --------------------------------------------------- two-block construction
def two_block(Wac, wop_H=0.75):
    """Wac in cm^-1.  Returns per-isotope (w_ac, w_op, W_op, <w^2>, rho)."""
    lam_ac = LAM["H"] * (1.0 - wop_H)          # constant absolute acoustic lambda
    out = {}
    for k in "HDT":
        wac = lam_ac / LAM[k]
        wop = 1.0 - wac
        # ln w_log = wac ln Wac + wop ln Wop
        lnWop = (np.log(WLOG[k]) - wac * np.log(Wac)) / wop
        Wop = np.exp(lnWop)
        w2 = wac * Wac ** 2 + wop * Wop ** 2
        out[k] = dict(wac=wac, wop=wop, Wop=Wop, w2=w2,
                      rho=w2 / WLOG[k] ** 2)
    return out


if __name__ == "__main__":
    print("=" * 86)
    print("PART A -- r-scan at adequate N (Nc can exceed N=400 at large r!)")
    print("=" * 86)
    print(f"{'r':>6}{'N':>6}{'Nc_H':>7}{'f_H':>9}{'f_D':>9}{'a_deficit':>11}{'a_budget':>10}")
    lnM = np.log(M["D"] / M["H"])
    a_ceil_proxy = 0.5 * np.log((LAM["H"] * WLOG["H"] ** 2) /
                                (LAM["D"] * WLOG["D"] ** 2)) / lnM
    for r, N in [(3.0, 800), (5.0, 800), (10.0, 800), (20.0, 1600), (40.0, 2400)]:
        from _chk_iso_core import Tc_einstein
        fH = frak_f(LAM["H"], MUS, N=N, r=r)
        fD = frak_f(LAM["D"], MUS, N=N, r=r)
        tH = Tc_einstein(LAM["H"], MUS, N=N, r=r)
        ncH = int(np.floor(r / (2 * np.pi * tH)))
        ad = alpha(fH, fD, "H", "D")
        print(f"{r:6.1f}{N:6d}{ncH:7d}{fH:9.5f}{fD:9.5f}{ad:+11.4f}"
              f"{a_ceil_proxy + ad:+10.4f}")
    print(f"   a_ceiling (proxy <w^2>=w_log^2) = {a_ceil_proxy:+.4f}")

    print()
    print("=" * 86)
    print("PART B -- two-block spectrum: supplies <w^2>, removes Einstein proxy")
    print("=" * 86)
    print(f"{'Wac':>6}{'iso':>5}{'w_ac':>7}{'w_op':>7}{'W_op':>8}"
          f"{'sqrt<w2>':>10}{'rho':>7}{'ceil K':>8}{'Tc_P K':>8}{'Tc_E K':>8}")
    results = {}
    for Wac in (100.0, 130.0, 160.0, 190.0):
        tb = two_block(Wac)
        res = {}
        for k in "HDT":
            b = tb[k]
            w = np.array([b["wac"], b["wop"]])
            om = np.array([Wac, b["Wop"]]) * CM1_TO_MEV
            rt2 = np.sqrt(b["w2"]) * CM1_TO_MEV
            ceil = C_INF * np.sqrt(LAM[k]) * rt2 * MEV_TO_K
            tcP = Tc_P(LAM[k], MUS, w, om, N=800) * MEV_TO_K
            tcE = frak_f(LAM[k], MUS, N=800) * ceil
            res[k] = dict(ceil=ceil, TcP=tcP, TcE=tcE, w2=b["w2"])
            print(f"{Wac:6.0f}{'Pd'+k:>5}{b['wac']:7.3f}{b['wop']:7.3f}"
                  f"{b['Wop']:8.1f}{np.sqrt(b['w2']):10.1f}{b['rho']:7.3f}"
                  f"{ceil:8.2f}{tcP:8.3f}{tcE:8.3f}")
        results[Wac] = res
        print(f"{'':6}{'-> alpha:':>12} "
              f"ceiling {alpha(res['H']['ceil'],res['D']['ceil'],'H','D'):+.4f}"
              f"   budget(2-block) {alpha(res['H']['TcP'],res['D']['TcP'],'H','D'):+.4f}"
              f"   budget(Einstein) {alpha(res['H']['TcE'],res['D']['TcE'],'H','D'):+.4f}"
              f"   [D->T 2-block {alpha(res['D']['TcP'],res['T']['TcP'],'D','T'):+.4f}]")
        print(f"{'':6}{'':12} s = rho_H/rho_D = "
              f"{two_block(Wac)['H']['rho']/two_block(Wac)['D']['rho']:.4f}"
              f"   d_BRD(H) = {100*(1-res['H']['TcP']/res['H']['TcE']):+.1f}%"
              f"   d_BRD(D) = {100*(1-res['D']['TcP']/res['D']['TcE']):+.1f}%")
        print()

    print("=" * 86)
    print("PART C -- two-block, mu* sweep (Wac=130 cm-1)")
    print("=" * 86)
    tb = two_block(130.0)
    print(f"{'mu*':>7}{'Tc_H K':>9}{'Tc_D K':>9}{'Tc_T K':>9}"
          f"{'a(H->D)':>10}{'a(D->T)':>10}")
    for mus in (0.0, 0.06, 0.085, 0.10, 0.13):
        tt = {}
        for k in "HDT":
            b = tb[k]
            w = np.array([b["wac"], b["wop"]])
            om = np.array([130.0, b["Wop"]]) * CM1_TO_MEV
            tt[k] = Tc_P(LAM[k], mus, w, om, N=800) * MEV_TO_K
        print(f"{mus:7.3f}{tt['H']:9.3f}{tt['D']:9.3f}{tt['T']:9.3f}"
              f"{alpha(tt['H'],tt['D'],'H','D'):+10.4f}"
              f"{alpha(tt['D'],tt['T'],'D','T'):+10.4f}")

    print()
    print("N-convergence spot check, two-block Wac=130, mu*=0.085:")
    b = tb["H"]
    for N in (400, 800, 1600):
        t = Tc_P(LAM["H"], MUS, np.array([b["wac"], b["wop"]]),
                 np.array([130.0, b["Wop"]]) * CM1_TO_MEV, N=N) * MEV_TO_K
        print(f"   N={N:5d}  Tc_P(PdH) = {t:.6f} K")
