"""LaH10 / LaD10 / CaH6: r-upper-bound from (lambda, omega_log), 2026-09-25.

Only lambda and omega_log are printed in the sources (no omega_bar_2).
Jensen / power-mean: omega_bar_2 = <w^2>^(1/2) >= omega_log under the SAME probability
measure 2 a2F/(w lambda) dw  ==> ceiling_log := C_inf*sqrt(lambda)*omega_log <= true ceiling
==> r_up := Tc_exp / ceiling_log >= true r.  One-sided: r_up <= 0.6260 implies true r <= 0.6260;
r_up > 0.6260 decides nothing.
Grade: numerical-evidence (external DFT inputs); the inequality itself is analytic.
"""
import math

C_INF = 0.1827262477            # sqrt(g(2))/(2 pi), analytic (repo)
KB_MEV = 0.08617333262          # meV/K
CM1_TO_K = 1.438776877          # K per cm^-1
TIGHTEST = 0.6260               # Nb3Ge


def ceil_K(lam, w_K):
    return C_INF * math.sqrt(lam) * w_K


# self-check: reproduce repo H3S-200 anharm ceiling 415.2 K (lam=1.721, nu2=149.26 meV)
chk = ceil_K(1.721, 149.26 / KB_MEV)
print(f"self-check H3S-200: ceiling={chk:.2f} K (repo 415.2) rel={chk/415.2-1:+.2e}")

# (label, P GPa, lambda, omega_log in K, Tc_calc, source)
rows = []
for P, lam, wl, tc in [(129, 3.62, 76.4, 230), (163, 2.67, 96.4, 225),
                       (214, 2.06, 115.5, 210), (264, 1.73, 126.6, 201)]:
    rows.append(("LaH10", P, lam, wl / KB_MEV, tc, "Errea2020 EDT-I (SCDFT Tc; ME 255.3/242.8/237.9/216.9)"))
for P, lam, wl, tc in [(159, 3.14, 63.5, 171), (210, 2.21, 81.7, 158), (260, 1.80, 92.2, 151)]:
    rows.append(("LaD10", P, lam, wl / KB_MEV, tc, "Errea2020 EDT-I (SCDFT Tc)"))
for P, lam, wl, tc in [(150, 2.75, 684, 238), (200, 1.99, 851, 220), (250, 1.71, 933, 210)]:
    rows.append(("CaH6", P, lam, wl * CM1_TO_K, tc, "Jeon2111.10797 (iso ME, mu*=0.13)"))

TCEXP = {
    "LaH10": [("Drozdov 250 (170 GPa)", 250.0), ("Somayazulu 260 (188-196 GPa)", 260.0),
              ("Somayazulu upper 280 (onset, 190-200 GPa)", 280.0)],
    "LaD10": [("Drozdov 168", 168.0)],
    "CaH6": [("Ma 215 (172 GPa)", 215.0), ("Li 210 (160-190 GPa)", 210.0)],
}
for lab, P, lam, wK, tc, src in rows:
    c = ceil_K(lam, wK)
    print(f"{lab:6s} P={P:3d} lam={lam:.2f} wlog={wK:7.1f} K  ceiling_log={c:6.1f} K  Tc_calc/ceil={tc/c:.4f}")
    for name, t in TCEXP[lab]:
        r = t / c
        flag = "<=0.626" if r <= TIGHTEST else ">0.626 (undecided)"
        print(f"      r_up[{name}] = {r:.4f}  {flag}")
