"""Tighten the CE threshold by replacing Lemma 2's worst-gap step with the exact
resolvent moments of u_perp at the top eigenvalue.

Repo Lemma 2 (P4D_coulomb_extremal_threshold.md step 3) already uses the exact
rank-one resolvent v_P = (A_P - lam)^{-1} u_c, but then bounds

    tan Theta = (mu_1 - lam)/c_1 * sqrt( sum_{j>=2} c_j^2/(lam-mu_j)^2 )

by |lam - mu_j| >= gam_P - 2 mu~ c_1^2 for EVERY j>=2 and (mu_1-lam) <= 2 mu~ c_1^2.
Both steps throw away the spectral distribution of u_perp.  Exactly:

    mu_1 - lam = 2 mu~ c_1^2 / (1 + 2 mu~ S_perp),      S_perp = sum c_j^2/(lam-mu_j)
    tan Theta  = 2 mu~ c_1 sqrt(S_2) / (1 + 2 mu~ S_perp),  S_2 = sum c_j^2/(lam-mu_j)^2

Since lam <= mu_1, (mu_1-mu_j)/(lam-mu_j) = 1 + (mu_1-lam)/(lam-mu_j) <= 1 + kap~
with kap~ := 2 mu~ c_1^2/(gam_P - 2 mu~ c_1^2), so with the mu~-FREE constants

    S2_0     = || (mu_1 - A_P)^{-1} u_perp ||^2
    Sperp_0  = < u_perp, (mu_1 - A_P)^{-1} u_perp >

one gets the sound bound

    tan Theta <= 2 mu~ c_1 (1+kap~) sqrt(S2_0) / (1 + 2 mu~ Sperp_0)          (NEW)

versus the repo bound

    tan Theta <= 2 mu~ c_1 ||u_perp|| / (gam_P - 2 mu~ c_1^2)                 (OLD)

Both are upper bounds; NEW <= OLD iff gam_P sqrt(S2_0) <= ||u_perp|| (always, since
1/(mu_1-mu_j) <= 1/gam_P) times the denominator gain.  This script (a) checks both
against the EXACT tan Theta computed from the perturbed top vector, (b) reports the
new threshold at the physical anchors, in particular the three near-line rows.

Strictness: numerical-evidence (float64, N=400).  The NEW inequality itself is
analytic (derivation above); only its evaluation is numerical.
"""
import os, math, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.optimize import brentq

sys.argv = ["x"]
import p4d_coulomb_extremal_threshold as ce


def constants(C):
    """Six repo constants + the two new resolvent moments, from A_P and u_c."""
    ev, V = np.linalg.eigh(C.K_P)
    mu1, psi1 = ev[-1], V[:, -1]
    if psi1.sum() < 0:
        psi1 = -psi1
    uc = ce.uc_vec(C.Nc, C.N) if hasattr(ce, "uc_vec") else C.c_P
    c1 = float(uc @ psi1)
    uperp = uc - c1 * psi1
    up = float(np.linalg.norm(uperp))
    gam = mu1 - ev[-2]
    cj = V[:, :-1].T @ uc                      # j = 1..N-1 ascending, excludes top
    dj = mu1 - ev[:-1]                         # > 0
    S2_0 = float(np.sum(cj ** 2 / dj ** 2))
    Sp_0 = float(np.sum(cj ** 2 / dj))
    return dict(mu1=mu1, psi1=psi1, c1=c1, up=up, gam=gam, S2_0=S2_0, Sp_0=Sp_0,
                cj=cj, dj=dj,
                q1=C.q1_d, b=C.b_d, m=C.m_d, ev=ev, V=V, uc=uc)


def tmax(K):
    b, m, q1 = K["b"], K["m"], K["q1"]
    return (-b + math.sqrt(b * b + m * q1)) / m


def tan_old(K, mut):
    den = K["gam"] - 2 * mut * K["c1"] ** 2
    return np.inf if den <= 0 else 2 * mut * K["c1"] * K["up"] / den


def tan_new(K, mut):
    den = K["gam"] - 2 * mut * K["c1"] ** 2
    if den <= 0:
        return np.inf
    kap = 2 * mut * K["c1"] ** 2 / den
    return (2 * mut * K["c1"] * (1 + kap) * math.sqrt(K["S2_0"])
            / (1 + 2 * mut * K["Sp_0"]))


def tan_new2(K, mut):
    """Same chain but with the EXACT shifted second moment instead of the
    crude (1+kap~) factor:
        del := 2 mu~ c1^2/(1 + 2 mu~ Sperp_0)  >=  mu_1 - lam      (Sperp >= Sperp_0)
        lam >= mu_1 - del =: lam_lo
        S_2(lam) <= S_2(lam_lo) = sum_j c_j^2 / ((mu_1-mu_j) - del)^2
        tan Theta = (mu_1-lam) sqrt(S_2(lam))/c1 <= del sqrt(S_2(lam_lo))/c1
    Sound for del < gam_P; returns inf otherwise."""
    dl = 2 * mut * K["c1"] ** 2 / (1 + 2 * mut * K["Sp_0"])
    d = K["dj"] - dl                       # (mu_1 - mu_j) - del, need > 0
    if d.min() <= 0:
        return np.inf
    S2 = float(np.sum(K["cj"] ** 2 / d ** 2))
    return dl * math.sqrt(S2) / K["c1"]


def lam_secular(K, mut):
    """Exact lam = lam_1(M_P) as the Lemma 1(i) secular root in (mu_2, mu_1):
           2 mu~ sum_j c_j^2/(mu_j - lam) = 1.
    Independent of np.linalg.eigh on M_P (that is tan_exact's route)."""
    ev, cc = K["ev"], K["V"].T @ K["uc"]
    m1, m2 = ev[-1], ev[-2]
    w = (m1 - m2) * 1e-13 + 1e-300

    def phi(l):
        return 2 * mut * float(np.sum(cc ** 2 / (ev - l))) - 1.0
    return brentq(phi, m2 + w, m1 - w, xtol=1e-16, rtol=1e-15)


def tan_sec(K, mut):
    """EXACT tan Theta from Lemma 1 alone (no Davis-Kahan, no surrogate for lam):
           tan Theta = (mu_1 - lam) sqrt(S_2(lam)) / c_1,
    with lam the secular root.  Analytic identity; only its evaluation is numeric."""
    lam = lam_secular(K, mut)
    d = K["dj"] - (K["mu1"] - lam)          # = lam - mu_j > 0 for j >= 2
    if d.min() <= 0:
        return np.inf
    S2 = float(np.sum(K["cj"] ** 2 / d ** 2))
    return (K["mu1"] - lam) * math.sqrt(S2) / K["c1"]


def tan_exact(K, mut):
    """exact angle between the top vector of M_P = A_P - 2 mu~ u_c u_c^T and psi1"""
    M = np.diag(K["ev"]) - 2 * mut * np.outer(K["V"].T @ K["uc"], K["V"].T @ K["uc"])
    w, U = np.linalg.eigh(M)
    v = K["V"] @ U[:, -1]
    co = abs(float(v @ K["psi1"]))
    return math.sqrt(max(0.0, 1 - co * co)) / co


def thr(K, f, cap_kappa=True):
    t = tmax(K)
    if cap_kappa:
        hi = K["gam"] / (2 * K["c1"] ** 2) * 0.999999
    else:
        # tan_new2 only needs del < gam (so lam_lo > mu_2); del(mu~) increases to
        # c1^2/Sperp_0.  If that sup is <= gam the bound never diverges.
        sup = K["c1"] ** 2 / K["Sp_0"]
        if sup <= K["gam"]:
            hi = 10.0
        else:
            hi = brentq(lambda m: 2*m*K["c1"]**2/(1+2*m*K["Sp_0"]) - K["gam"]*0.999999,
                        1e-12, 1e6)
    if f(K, hi) < t:
        return hi
    return brentq(lambda m: f(K, m) - t, 1e-12, hi, xtol=1e-14, rtol=1e-13)


ANCH = [("LaH10 129", 3.62, 0.10, None), ("LaH10 163", 2.67, 0.10, 0.783),
        ("LaH10 214", 2.06, 0.10, None), ("LaH10 264", 1.73, 0.10, None),
        ("H3S 200", 1.84, 0.16, None), ("Nb", 1.23, 0.13, None),
        ("MgB2", 0.62, 0.13, None)]

if __name__ == "__main__":
    print(f"{'anchor':>10} {'r':>2} {'Nc':>3} {'mu~':>7} | "
          f"{'thr_OLD':>8} {'thr_N2':>8} {'thr_SEC':>8} | "
          f"{'m_OLD':>6} {'m_N2':>6} {'m_SEC':>6} {'v':>9} | "
          f"{'kapcap':>7} {'tmax':>7} {'old':>8} {'new2':>8} {'sec':>7} {'eigh':>7} "
          f"{'|sec-eigh|':>10} {'ord':>5}")
    worst = 0.0
    for name, lam, mus, vr in ANCH:
        vr = vr or 0.8709 / math.sqrt(lam)
        mut = mus / lam
        oc = 2 * math.pi * vr
        for r in (3, 10):
            Nc = max(1, int(r * vr))
            C = ce.CellExtremal(np.array([0.5, 0.5]), np.array([oc / 1.05, oc * 1.05]),
                                1.0, Nc, NN=400)
            K = constants(C)
            t_o = thr(K, tan_old)
            t_2 = thr(K, tan_new2)
            t_s = thr(K, tan_sec, cap_kappa=False)
            te, ts = tan_exact(K, mut), tan_sec(K, mut)
            to, t2, tm = tan_old(K, mut), tan_new2(K, mut), tmax(K)
            dd = abs(ts - te) / max(te, 1e-300)
            worst = max(worst, dd)
            ok = "OK" if ts <= t2 + 1e-12 <= to + 1e-12 else "**ORD**"
            v = ("Y" if mut < t_o else "n") + "->" + ("Y" if mut < t_2 else "n") \
                + "/" + ("Y" if mut < t_s else "n")
            print(f"{name:>10} {r:>2} {Nc:>3} {mut:>7.4f} | "
                  f"{t_o:>8.4f} {t_2:>8.4f} {t_s:>8.4f} | "
                  f"{t_o/mut:>6.3f} {t_2/mut:>6.3f} {t_s/mut:>6.3f} {v:>9} | "
                  f"{K['gam']/(2*K['c1']**2):>7.4f} {tm:>7.4f} "
                  f"{to:>8.4f} {t2:>8.4f} {ts:>7.4f} {te:>7.4f} {dd:>10.2e} {ok:>5}")
    print(f"\nworst |tan_sec - tan_eigh|/tan  = {worst:.3e}   "
          f"(two independent routes: secular root vs eigh of M_P)")

    # ---- does the secular route survive kappa >= 1? -----------------------
    # thr_SEC exceeds gam/(2c1^2) on several rows, i.e. it claims certification
    # in a regime where the repo's Lemma 1 hypothesis kappa < 1 fails.  Lemma
    # 1(iii) and Lemma 2 need kappa < 1; the secular route uses NEITHER -- only
    # 1(i)/(ii), whose proof (phi strictly increasing, phi -> -inf at mu_2+,
    # phi -> +inf at mu_1-) needs just c_1 != 0, c_2 != 0.  Test that claim
    # directly at mu~ ABOVE the cap.
    print("\nkappa >= 1 probe: lam_secular vs lam_1(eigh M_P), mu~ = 2 x gam/(2c1^2)")
    print(f"{'anchor':>10} {'r':>2} {'mu~':>8} {'kappa':>6} "
          f"{'|lam_sec-lam_eigh|':>19} {'|tan_sec-tan_eigh|/tan':>23}")
    for name, lam, mus, vr in ANCH:
        vr = vr or 0.8709 / math.sqrt(lam)
        oc = 2 * math.pi * vr
        for r in (3, 10):
            Nc = max(1, int(r * vr))
            C = ce.CellExtremal(np.array([0.5, 0.5]), np.array([oc / 1.05, oc * 1.05]),
                                1.0, Nc, NN=400)
            K = constants(C)
            cap = K["gam"] / (2 * K["c1"] ** 2)
            mt = 2.0 * cap                       # kappa = 2, well outside Lemma 1's hyp.
            ls = lam_secular(K, mt)
            M = np.diag(K["ev"]) - 2*mt*np.outer(K["V"].T @ K["uc"], K["V"].T @ K["uc"])
            le = float(np.linalg.eigvalsh(M)[-1])
            ts, te = tan_sec(K, mt), tan_exact(K, mt)
            print(f"{name:>10} {r:>2} {mt:>8.4f} {2*mt*K['c1']**2/K['gam']:>6.2f} "
                  f"{abs(ls-le):>19.3e} {abs(ts-te)/te:>23.3e}")
