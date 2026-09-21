"""PdH/PdD/PdT isotope exponent from the ceiling + deficit budget.

Inputs (Errea, Calandra, Mauri, PRL 111, 177002 (2013) = arXiv:1305.7123, TABLE I,
read off the in-repo pdftotext dump _pdftxt/2013_Ion_Errea_1305.7123.txt):
    harmonic : lam = 1.55 (all three), w_log = 205 / 150 / 125 cm^-1, Tc = 47/34/30 K
    SSCHA    : lam = 0.40/0.46/0.48,   w_log = 405 / 304 / 257 cm^-1, Tc = 5.0/6.5/6.9 K
    mu* = 0.085 (fixed).  Exp Tc: 9/11 (ref 7)  or  8/10 K (ref 8).

<w^2> is NOT reported anywhere for PdH.  w_log^2 is used as a PROXY; the sign of
the proxy error in the ISOTOPE RATIO is undetermined -> see sensitivity block.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from _chk_iso_core import C_INF, frak_f, Tc_einstein

CM1_TO_MEV = 0.1239841984          # h c / e * 100  (CODATA 2018)
MEV_TO_K = 11.604518121            # e / (1000 k_B)

M = {"H": 1.00782503, "D": 2.01410178, "T": 3.01604928}   # u, atomic masses

SSCHA = {"H": dict(lam=0.40, wlog_cm=405.0, Tc_me=5.0),
         "D": dict(lam=0.46, wlog_cm=304.0, Tc_me=6.5),
         "T": dict(lam=0.48, wlog_cm=257.0, Tc_me=6.9)}
HARM = {"H": dict(lam=1.55, wlog_cm=205.0, Tc_me=47.0),
        "D": dict(lam=1.55, wlog_cm=150.0, Tc_me=34.0),
        "T": dict(lam=1.55, wlog_cm=125.0, Tc_me=30.0)}
MUS = 0.085


def alpha(TcA, TcB, A, B):
    """alpha = -(ln Tc_B - ln Tc_A)/(ln M_B - ln M_A), Errea's TABLE I definition."""
    return -np.log(TcB / TcA) / np.log(M[B] / M[A])


def row(d, mus, r=10.0, N=400, rho=1.0):
    """rho = <w^2>/w_log^2 (proxy correction factor, 1.0 = pure proxy)."""
    w2bar = d["wlog_cm"] * CM1_TO_MEV * np.sqrt(rho)     # sqrt(<w^2>) in meV
    ceil_K = C_INF * np.sqrt(d["lam"]) * w2bar * MEV_TO_K
    f = frak_f(d["lam"], mus, N=N, r=r)
    return ceil_K, f, f * ceil_K


def table(data, mus, r=10.0, rho=None, label=""):
    rho = rho or {k: 1.0 for k in data}
    out = {}
    print(f"\n--- {label}  (mu*={mus}, r={r}) ---")
    print(f"{'iso':>4}{'lam':>7}{'w_log meV':>11}{'ceiling K':>11}"
          f"{'f':>9}{'Tc_bud K':>10}{'Tc_ME K':>9}")
    for k in data:
        c, f, tb = row(data[k], mus, r, rho=rho[k])
        out[k] = dict(ceil=c, f=f, Tc=tb)
        print(f"{'Pd'+k:>4}{data[k]['lam']:7.2f}"
              f"{data[k]['wlog_cm']*CM1_TO_MEV:11.3f}{c:11.3f}{f:9.5f}{tb:10.3f}"
              f"{data[k]['Tc_me']:9.2f}")
    return out


def alphas(out, data, A="H", B="D"):
    a_ceil = alpha(out[A]["ceil"], out[B]["ceil"], A, B)
    a_def = alpha(out[A]["f"], out[B]["f"], A, B)
    a_bud = alpha(out[A]["Tc"], out[B]["Tc"], A, B)
    a_me = alpha(data[A]["Tc_me"], data[B]["Tc_me"], A, B)
    return a_ceil, a_def, a_bud, a_me


if __name__ == "__main__":
    print("=" * 78)
    print("PdH ISOTOPE TEST -- ceiling + parameter-free deficit budget")
    print(f"C_inf = {C_INF},  proxy <w^2> = w_log^2")
    print("=" * 78)

    # ---- experimental alphas, recomputed from the raw Tc pairs -------------
    print("\n[EXP] alpha from raw Tc pairs (Errea TABLE I 'Exp.' columns):")
    for ref, (th, td) in {"ref[7] Stritzker-Buckel 1972": (9.0, 11.0),
                          "ref[8] Schirber-Northrup 1974": (8.0, 10.0)}.items():
        print(f"   {ref:32s} Tc {th}/{td} K  ->  alpha(H->D) = {alpha(th,td,'H','D'):+.4f}")
    print(f"   Errea SSCHA-ME   5.0/6.5 K -> alpha(H->D) = {alpha(5.0,6.5,'H','D'):+.4f}"
          f"   (paper prints -0.38)")
    print(f"   Errea SSCHA-ME   5.0/6.9 K -> alpha(H->T) = {alpha(5.0,6.9,'H','T'):+.4f}"
          f"   (paper prints -0.29)")
    print(f"   Errea harm-ME   47.0/34.0 K -> alpha(H->D) = {alpha(47.,34.,'H','D'):+.4f}"
          f"   (paper prints +0.47)")

    # ---- main table --------------------------------------------------------
    sa = table(SSCHA, MUS, label="SSCHA (anharmonic) alpha2F params")
    ha = table(HARM, MUS, label="harmonic params (control)")

    print("\n[ALPHA] SSCHA row:")
    print(f"{'pair':>8}{'a_ceiling':>11}{'a_deficit':>11}{'a_budget':>10}"
          f"{'a_SSCHA-ME':>12}{'a_exp':>9}")
    for A, B, aexp in [("H", "D", -0.29), ("D", "T", None), ("H", "T", None)]:
        ac, ad, ab, am = alphas(sa, SSCHA, A, B)
        e = f"{aexp:+.3f}" if aexp is not None else "  n/a"
        print(f"{'Pd'+A+'->'+B:>8}{ac:+11.4f}{ad:+11.4f}{ab:+10.4f}{am:+12.4f}{e:>9}")
    print("\n[ALPHA] harmonic control row (lam mass-independent):")
    for A, B in [("H", "D"), ("D", "T")]:
        ac, ad, ab, am = alphas(ha, HARM, A, B)
        print(f"{'Pd'+A+'->'+B:>8}{ac:+11.4f}{ad:+11.4f}{ab:+10.4f}{am:+12.4f}")

    # ---- residual amplification -------------------------------------------
    print("\n[AMPLIFICATION] Tc residuals -> alpha error")
    for A, B in [("H", "D")]:
        for tag, tex in [("exp ref[7] 9/11 K", (9.0, 11.0)),
                         ("exp ref[8] 8/10 K", (8.0, 10.0)),
                         ("SSCHA-ME 5.0/6.5 K", (5.0, 6.5))]:
            rA = sa[A]["Tc"] / tex[0] - 1.0
            rB = sa[B]["Tc"] / tex[1] - 1.0
            da = np.log((1 + rA) / (1 + rB)) / np.log(M[B] / M[A])
            print(f"   vs {tag:20s} residual H {100*rA:+7.1f}%  D {100*rB:+7.1f}%"
                  f"   -> d(alpha) = {da:+.4f}")

    # ---- sensitivity: mu* --------------------------------------------------
    print("\n[SENS-1] mu* (r=10, proxy):")
    print(f"{'mu*':>7}{'f_H':>9}{'f_D':>9}{'a_deficit':>11}{'a_budget':>10}"
          f"{'Tc_H K':>9}{'Tc_D K':>9}")
    for mus in (0.0, 0.04, 0.06, 0.085, 0.10, 0.13, 0.16):
        o = {k: dict(zip(("ceil", "f", "Tc"), row(SSCHA[k], mus))) for k in "HD"}
        ac, ad, ab, _ = alphas(o, SSCHA, "H", "D")
        print(f"{mus:7.3f}{o['H']['f']:9.5f}{o['D']['f']:9.5f}"
              f"{ad:+11.4f}{ab:+10.4f}{o['H']['Tc']:9.3f}{o['D']['Tc']:9.3f}")

    # ---- sensitivity: r ----------------------------------------------------
    print("\n[SENS-2] r = w_c/sqrt(<w^2>)  (mu*=0.085, proxy):")
    print(f"{'r':>6}{'Nc_H':>6}{'Nc_D':>6}{'f_H':>9}{'f_D':>9}"
          f"{'a_deficit':>11}{'a_budget':>10}{'Tc_H K':>9}")
    for r in (3.0, 5.0, 10.0, 20.0, 40.0):
        o = {k: dict(zip(("ceil", "f", "Tc"), row(SSCHA[k], MUS, r=r))) for k in "HD"}
        ac, ad, ab, _ = alphas(o, SSCHA, "H", "D")
        nc = {k: int(np.floor(r / (2 * np.pi * Tc_einstein(SSCHA[k]["lam"], MUS, r=r))))
              for k in "HD"}
        print(f"{r:6.1f}{nc['H']:6d}{nc['D']:6d}{o['H']['f']:9.5f}{o['D']['f']:9.5f}"
              f"{ad:+11.4f}{ab:+10.4f}{o['H']['Tc']:9.3f}")

    # ---- sensitivity: the proxy -------------------------------------------
    print("\n[SENS-3] proxy: s := (<w^2>/w_log^2)_H / (<w^2>/w_log^2)_D")
    print("   a_ceiling shifts by (1/2)ln(s)/ln(M_D/M_H); a_deficit unchanged")
    ac0 = alphas(sa, SSCHA, "H", "D")[0]
    ab0 = alphas(sa, SSCHA, "H", "D")[2]
    lnM = np.log(M["D"] / M["H"])
    print(f"{'s':>8}{'a_ceiling':>11}{'a_budget':>10}")
    for s in (0.70, 0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.30):
        print(f"{s:8.2f}{ac0+0.5*np.log(s)/lnM:+11.4f}{ab0+0.5*np.log(s)/lnM:+10.4f}")
    need = np.exp(2 * lnM * (-0.29 - ab0))
    print(f"   s needed to move a_budget to the experimental -0.29: s = {need:.4f}"
          f"   ({100*(need-1):+.1f}%)")

    # ---- sensitivity: lambda ----------------------------------------------
    print("\n[SENS-4] lambda (the driver).  Vary d(lam) = lam_D - lam_H at fixed lam_H:")
    print(f"{'lam_D':>8}{'f_D':>9}{'a_deficit':>11}{'a_budget':>10}")
    for lamD in (0.40, 0.42, 0.44, 0.46, 0.48, 0.50):
        dD = dict(SSCHA["D"]); dD["lam"] = lamD
        o = {"H": dict(zip(("ceil", "f", "Tc"), row(SSCHA["H"], MUS))),
             "D": dict(zip(("ceil", "f", "Tc"), row(dD, MUS)))}
        ac, ad, ab, _ = alphas(o, {"H": SSCHA["H"], "D": dD}, "H", "D")
        print(f"{lamD:8.2f}{o['D']['f']:9.5f}{ad:+11.4f}{ab:+10.4f}")
