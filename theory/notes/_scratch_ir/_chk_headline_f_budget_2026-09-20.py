"""
Task C (budget) + Task B (decomposition) recompute for the deficit-factor headline note.
Independent from _chk_fpred_zero_param_2026-09-18.py: this script only consumes the
ALREADY-PUBLISHED numbers in deficit_budget_h3s_2026-09-17.md Table 1 and
h3s_tc_exp_pressure_audit_2026-09-18.md's corrected anchors (184.04 K / 145.5 K), and
recomputes the telescoping decomposition + degrees-of-freedom accounting from scratch.
No new Tc_P solve here (that would just re-import deficit_independent_check.py; this
script's job is bookkeeping/consistency, not a fourth code path for Tc_P).
"""
import math

# rows: name, lam, ceiling_K, TcE0, TcP0, TcP13, Tc_exp_legacy, Tc_exp_corrected(or None)
rows = [
    ("H3S-200anharm", 1.7210, 415.21, 327.44, 283.34, 198.974, 203.0, 184.04),
    ("D3S-200anharm", 1.7838, 298.58, 237.79, 214.42, 155.358, 150.0, 145.5),
    ("H3S-160anharm", 1.8760, 400.42, 323.09, 290.73, 211.156, 203.0, None),  # no digitized 160GPa anchor
]

print(f"{'row':16s} {'f_pred':>8s} {'f_obs(legacy)':>14s} {'f_obs(corr)':>12s} {'resid_legacy':>13s} {'resid_corr':>11s}")
for name, lam, ceil, tce0, tcp0, tcp13, texp_l, texp_c in rows:
    fpred = tcp13/ceil
    fobs_l = texp_l/ceil
    resid_l = (fpred-fobs_l)/fobs_l*100
    if texp_c is not None:
        fobs_c = texp_c/ceil
        resid_c = (fpred-fobs_c)/fobs_c*100
        print(f"{name:16s} {fpred:8.4f} {fobs_l:14.4f} {fobs_c:12.4f} {resid_l:12.2f}% {resid_c:10.2f}%")
    else:
        print(f"{name:16s} {fpred:8.4f} {fobs_l:14.4f} {'N/A':>12s} {resid_l:12.2f}% {'N/A':>10s}")

print()
print("=== Task C: budget degradation ===")
# Degrees of freedom check: H3S-200 and H3S-160 share the same compound (H3S) and same
# Mishra EPW workflow; D3S-200 shares everything with H3S-200 except isotope mass.
# => true independent "materials families" represented = 1 (H3S/D3S system only).
print("Independent materials families in current non-circular sample: 1 (H3S/D3S, single EPW/SSCHA workflow).")
print("H3S-200 vs H3S-160: same compound, same workflow -> not independent (pressure variant only).")
print("H3S vs D3S: same electronic structure, same workflow, differ only in nuclear mass -> systematic")
print("  errors (DFT functional, pseudopotential, alpha2F broadening model, mu* convention) are shared.")
print()
print("Degrading the claim: even if we allow f in [0.3,0.7] (task's stress-test window), the headline")
print("'f~0.5 stable across materials' requires >=2 INDEPENDENT families (different compounds, ideally")
print("different DFT codes/groups) to mean anything -- with 1 family, any sample size N (adding more")
print("pressures or isotopes of the SAME compound) cannot buy independence. Budget for sample count is")
print("therefore *unbounded* in the wrong direction: the quantity that is actually scarce is FAMILY count,")
print("not row count. This matches CLAUDE.md 'count degrees of freedom before N independent tests'.")
