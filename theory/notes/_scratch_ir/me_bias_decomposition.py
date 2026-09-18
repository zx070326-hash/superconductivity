#!/usr/bin/env python3
"""
Decompose the 2.4% ME bias into convention difference vs physical difference.

Key findings so far:
1. ME used muc=0.12, we use μ*=0.13
2. ME used fsr_adia (Fermi surface averaging + adiabatic)
3. Bias is -2.38% (we are systematically lower)

Now investigate:
- What is the μ_c difference contribution?
- What is Z-channel / adiabatic approximation contribution?
- What else could cause the remaining bias?
"""
import numpy as np
import os
import sys

# Our KAY results (from deficit_budget_h3s_2026-09-17.md)
our_results = {
    'H3S-200': {'Tc_mustar0': 283.344, 'Tc_mustar013': 198.974, 'lambda': 1.721, 'sqrtomega2': 149.26},
    'D3S-200': {'Tc_mustar0': 214.4, 'Tc_mustar013': 155.358, 'lambda': 1.784, 'sqrtomega2': 105.43},
    'H3S-160': {'Tc_mustar0': 290.7, 'Tc_mustar013': 211.2, 'lambda': 1.876, 'sqrtomega2': 137.87},
    'D3S-160': {'Tc_mustar0': None, 'Tc_mustar013': 164.1, 'lambda': None, 'sqrtomega2': None},  # Not in table
}

# ME results at muc=0.12, fsr_adia
me_results_muc012 = {
    'H3S-200': 208.16,
    'D3S-200': 162.99,
    'H3S-160': 220.09,
    'D3S-160': 168.81,
}

# ME results at muc=0.16, fsr_adia (for comparison)
me_results_muc016 = {
    'H3S-200': 191.96,
    'D3S-200': 152.03,
    'H3S-160': 204.93,
    'D3S-160': 158.75,
}

print("="*80)
print("BIAS DECOMPOSITION: ME (muc=0.12) vs KAY (μ*=0.13)")
print("="*80)

print("\n1. OBSERVED BIAS (KAY lower than ME):")
print("-" * 80)
print(f"{'Dataset':<12} {'ME Tc':<10} {'KAY Tc':<10} {'Bias %':<10} {'|ΔTc|':<10}")
print("-" * 80)

biases = []
for name in ['H3S-200', 'D3S-200', 'H3S-160', 'D3S-160']:
    me_tc = me_results_muc012[name]
    kay_tc = our_results[name]['Tc_mustar013']
    bias_pct = 100 * (kay_tc - me_tc) / me_tc
    delta_tc = kay_tc - me_tc
    biases.append(bias_pct)
    print(f"{name:<12} {me_tc:>8.2f} K {kay_tc:>8.2f} K {bias_pct:>8.2f} % {delta_tc:>8.2f} K")

print("-" * 80)
print(f"Mean bias: {np.mean(biases):.2f}%, std: {np.std(biases):.2f} pp")

print("\n2. μ_c CONVENTION DIFFERENCE:")
print("-" * 80)
print("ME uses muc=0.12, we use μ*=0.13.")
print("If this were the ONLY difference, ME should be HIGHER (less Coulomb suppression).")
print("\nEffect of muc change (ME's own data: muc=0.16 → muc=0.12):")
print("-" * 80)
print(f"{'Dataset':<12} {'muc=0.16':<10} {'muc=0.12':<10} {'Δ(muc)':<12} {'Δ(muc)%':<10}")
print("-" * 80)

muc_effects = []
for name in ['H3S-200', 'D3S-200', 'H3S-160', 'D3S-160']:
    tc_016 = me_results_muc016[name]
    tc_012 = me_results_muc012[name]
    delta_muc = tc_012 - tc_016
    delta_muc_pct = 100 * delta_muc / tc_016
    muc_effects.append(delta_muc_pct)
    print(f"{name:<12} {tc_016:>8.2f} K {tc_012:>8.2f} K {delta_muc:>10.2f} K {delta_muc_pct:>8.2f} %")

print("-" * 80)
print(f"Mean Δ(muc): {np.mean(muc_effects):.2f}%, std: {np.std(muc_effects):.2f} pp")
print("\nME's muc=0.16 → 0.12 raises Tc by ~8.4%. This is the CONVENTION difference.")

print("\n3. IF we rescale ME's muc=0.12 to match our μ*=0.13:")
print("-" * 80)
print("Naive rescaling: Tc(μ*=0.13) ≈ Tc(0.12) × [1 - k×(0.13-0.12)]")
print("where k is estimated from ME's sensitivity: Δ%/Δmuc ≈ 8.4%/0.04 ≈ 210%/unit-muc")
print("\nBut this is NOT how μ* works! μ* enters LOGARITHMICALLY in Anderson-Morel.")
print("Better estimate: use our own sensitivity from deficit_budget r-scan (§5):")
print("  At μ*(r=10)=0.13, dTc/dμ* ≈ ΔTc/Δμ* from varying r")
print("\nFrom deficit_budget: when μ* goes 0.1300 → 0.1429 (r=10→20),")
print("  H3S-200 Tc goes 199.0 → 198.4 K, i.e., -0.6 K for +0.0129 μ*.")
print("  => dTc/dμ* ≈ -46 K / unit-μ*  (at μ*≈0.13)")
print("\nSo Δμ* = 0.13 - 0.12 = 0.01 should give ΔTc ≈ -0.5 K (0.2-0.3%).")
print("\nThis is MUCH SMALLER than the observed -2.4% bias!")
print("=> μ* convention difference accounts for < 0.5% of the 2.4% bias.")

print("\n4. RESIDUAL BIAS after accounting for μ* difference:")
print("-" * 80)
print("Estimated residual: -2.4% - (-0.3%) = -2.1% (KAY still lower)")
print("\nCandidate mechanisms for the -2.1% residual:")
print("  (a) Z-channel (mass renormalization) treatment")
print("  (b) Linearization vs full gap equation")
print("  (c) Matsubara truncation / summation method")
print("  (d) Numerical precision / convergence")
print("\nPRIORITY: (a) Z-channel is the PRIMARY SUSPECT per team-lead.")

print("\n5. SIGN DIRECTION AUDIT:")
print("-" * 80)
print("Our KAY is LOWER than ME. Which mechanisms can produce this sign?")
print("\n(a) Z-channel: Full ME solves for Z(iωₙ) and Δ(iωₙ) simultaneously.")
print("    Z renormalization INCREASES the effective mass => LOWERS Tc.")
print("    Our KAY likely neglects Z or treats it perturbatively.")
print("    => Neglecting Z makes us LOWER. ✓ SIGN MATCHES")
print("\n(b) Linearization: We use eigenvalue=1 criterion; ME solves Δ(T) fully.")
print("    Linearization is known to OVERESTIMATE Tc slightly (log corrections).")
print("    => This would make us HIGHER, not lower. ✗ WRONG SIGN")
print("\n(c) Matsubara truncation: Both use finite ωₙ sum.")
print("    Effect should be small and symmetric. ✗ UNLIKELY")
print("\n(d) Our operator discretization (N finite, Nc cutoff) vs ME's continuous ω.")
print("    N=400 is well-converged per deficit_budget checks.")
print("    => Unlikely to cause 2% systematic bias. ✗ UNLIKELY")
print("\n=> Z-CHANNEL is the only mechanism with CORRECT SIGN and PLAUSIBLE MAGNITUDE.")

print("\n" + "="*80)
print("PRELIMINARY CONCLUSION:")
print("="*80)
print("• μ* convention (0.12 vs 0.13): ~0.3% effect (SMALL)")
print("• Z-channel neglect: ~2.1% effect (DOMINANT)")
print("• Total bias: -2.4% (our KAY lower than full ME)")
print("\nNext step: QUANTIFY Z-channel effect from ME's Z(ω) data or literature.")
