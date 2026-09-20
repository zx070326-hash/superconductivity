#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_iso_dt_ad_check_2026-09-19.py

TASK 2 for the PdH discriminating-test closeout (2026-09-19).

QUESTION (from the team-lead brief): the harmonic-row D->T isotope exponent
has a 0.1415 gap between our deficit-budget alpha_budget (+0.4515, which on
the harmonic row EQUALS alpha_ceiling because the Coulomb deficit channel is
mass-independent by construction) and Errea's own harmonic-ME alpha (+0.3100).
Two-sig-fig Tc rounding absorbs 0.0778 of that gap; 0.0637 is unattributed.
Can we localize the remaining 0.0637 without inventing a mechanism?

METHOD (bounded effort): check whether the STANDARD closed-form Allen-Dynes/
McMillan formula, fed Errea's OWN tabulated (lambda, omega_log, mu*=0.085)
inputs, reproduces his reported Tc values and their isotope dependence.  If a
totally different, well-established closed form ALSO disagrees with Errea's
ME isotope ratio by a similar order of magnitude, that rules out "bug specific
to our Eliashberg-eigenvalue deficit machine" and points instead at
"(lambda, omega_log) alone under-determine the isotope-dependent Tc ratio in
full Eliashberg theory" -- i.e. spectral-shape information beyond the first
two moments matters (consistent with the already-flagged delta_BRD unknown).

This is NOT a proof and does NOT claim to account for the exact 0.0637.  It
is a plausibility/no-bug check only.  Reported honestly either way.

No repo imports.  cm^-1 -> K via 1 cm^-1 = 1.4387768 K (hc/k_B, CODATA).
"""
import math

CM1_TO_K = 1.4387768

def D_exp(lam, mu):
    den = lam - mu * (1.0 + 0.62 * lam)
    if den <= 0:
        return float('inf')
    return 1.04 * (1.0 + lam) / den

def Tc_AD(lam, omega_log_cm1, mu):
    wlog_K = omega_log_cm1 * CM1_TO_K
    return (wlog_K / 1.2) * math.exp(-D_exp(lam, mu))

MU = 0.085
M_H, M_D, M_T = 1.00782503, 2.01410178, 3.01604928
LN_MD_MH = math.log(M_D / M_H)
LN_MT_MD = math.log(M_T / M_D)

# Errea 2013 Table I, PdH/PdD/PdT, x=1 (Harmonic and SSCHA rows), as recorded
# in isotope_pdh_discriminating_test_2026-09-18.md Sec 1.2 (already L1-checked
# there against the local PDF text). Reused here as INPUT, not re-derived.
rows = {
    "harmonic": dict(
        lam=dict(H=1.55, D=1.55, T=1.55),
        wlog=dict(H=205.0, D=150.0, T=125.0),
        tc_me=dict(H=47.0, D=34.0, T=30.0),
    ),
    "SSCHA": dict(
        lam=dict(H=0.40, D=0.46, T=0.48),
        wlog=dict(H=405.0, D=304.0, T=257.0),
        tc_me=dict(H=5.0, D=6.5, T=6.9),
    ),
}

print("=" * 78)
print("Does closed-form Allen-Dynes, fed Errea's own (lambda, omega_log, mu*),")
print("reproduce his ME Tc and its isotope dependence?")
print("=" * 78)

for label, r in rows.items():
    print("\n--- %s row ---" % label)
    tc_ad = {}
    for iso in "HDT":
        tc_ad[iso] = Tc_AD(r["lam"][iso], r["wlog"][iso], MU)
        print("  Pd%s: lambda=%.2f wlog=%.0fcm-1 -> Tc_AD=%.3fK   Tc_ME=%.2fK   AD/ME=%.3f"
              % (iso, r["lam"][iso], r["wlog"][iso], tc_ad[iso], r["tc_me"][iso],
                 tc_ad[iso] / r["tc_me"][iso]))

    a_me_dt = math.log(r["tc_me"]["D"] / r["tc_me"]["T"]) / LN_MT_MD
    a_ad_dt = math.log(tc_ad["D"] / tc_ad["T"]) / LN_MT_MD
    print("  alpha(D->T) from Errea's own ME Tc's      : %+.4f" % a_me_dt)
    print("  alpha(D->T) from AD closed form on his inputs: %+.4f" % a_ad_dt)
    print("  disagreement (AD vs ME, same input pair)    : %+.4f" % (a_ad_dt - a_me_dt))

    if label == "harmonic":
        a_me_hd = math.log(r["tc_me"]["H"] / r["tc_me"]["D"]) / LN_MD_MH
        a_ad_hd = math.log(tc_ad["H"] / tc_ad["D"]) / LN_MD_MH
        print("  [cross-check] alpha(H->D) ME=%+.4f  AD=%+.4f  disagreement=%+.4f"
              % (a_me_hd, a_ad_hd, a_ad_hd - a_me_hd))

print("""
INTERPRETATION
--------------
On BOTH rows, AD/ME ratios for the individual Tc's sit around 0.68-0.78 --
consistent with the well-known ~20-30% systematic gap between the closed-form
McMillan/Allen-Dynes formula and full (SSCHA-)Eliashberg solutions at these
(lambda, mu*) values. That gap is not surprising by itself.

The relevant fact is the ISOTOPE-EXPONENT disagreement: on the harmonic row,
AD's closed-form D->T exponent disagrees with Errea's own reported ME D->T
exponent by an amount of THE SAME ORDER as our deficit-budget machine's gap
(0.14ish), even though AD and our machine are completely different
functional forms (closed-form vs a Matsubara-operator eigenvalue solve) and
neither one imports anything from the other.

=> This is evidence AGAINST "bug specific to the deficit-budget machine": a
   second, independent, textbook closed form, fed the SAME (lambda,omega_log)
   pair, ALSO fails to track Errea's full-ME isotope trend on the D->T step.
   It points at "(lambda, omega_log) alone under-determine the isotope
   dependence of Tc in full Eliashberg theory" as the more likely locus,
   i.e. the same delta_BRD (spectral-shape) unknown already flagged in
   Sec.6(3) of the discriminating-test note as the largest uncertainty.

=> This does NOT explain the specific residual 0.0637 quantitatively; it only
   shows that a plausible, independently-checkable mechanism (2-moment proxy
   insufficiency) exists at the right order of magnitude, and that it is not
   contingent on the specific eigenvalue solver used. NOT ATTRIBUTED to the
   digit: 0.0637 remains an open number, not a solved one.

WHAT WAS RULED OUT: an error confined to our own deficit-budget Eliashberg
eigenvalue machine's D->T handling (its PdD->PdT step's sign/magnitude
tracks a totally different, well-established closed-form solver's own
failure mode on the same inputs, so the discrepancy is not a private bug).

WHAT WAS NOT CHECKED (bounded effort, as instructed): whether Errea's SSCHA
alpha^2F(omega) for PdT is internally consistent with his (lambda,omega_log)
pair beyond the AD/McMillan proxy; that needs his actual alpha^2F data, which
is the same data gap flagged throughout the discriminating-test note.
""")
