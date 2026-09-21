# Gap ratio 2Δ/kTc from eigenvector profile — 2026-09-17

## Problem

Can the closed-form eigenvector profile discovered for P4's linearized operator predict the **tunneling gap ratio 2Δ(0)/k_BTc**? The profile fits Matsubara gap functions at T=Tc to 0.5–1.5% residual:

```
τ_n = u_n · a²/(a² + (n+½)²),   u_n = (2n+1)^{-1/2},   a = α(ϖ)·ϖ
```

where ϖ = ω_log/√(λω₂). H3S has experimental tunneling data (Du et al. 2025): **2Δ/kTc ≈ 3.43**.

## Background: What 2Δ/kTc measures

In BCS weak coupling: 2Δ(0)/kTc = 3.5275 (universal). Strong coupling modifies this:
- Allen-Dynes/Carbotte empirical relations correlate it with Tc/ω_log
- Width parameter α(ϖ) from the eigenvector profile is a **direct measure of Matsubara gap structure**
- Question: Is α(ϖ) a sharper strong-coupling parameter than Tc/ω_log?

## Method

1. **Direct calculation**: Use `eliashberg_iso.py` full nonlinear solver on H3S realistic α²F (Mishra 160 GPa anharmonic)
2. **BCS limit check**: Verify solver gives 3.53 at weak coupling
3. **Profile analysis**: Fit eigenvector near Tc, extract α(ϖ)
4. **mu* sensitivity**: Scan μ* ∈ [0.10, 0.20] to understand deficit budget

## Numerical verification ([numerical-evidence], float64)

### Input: H3S 160 GPa anharmonic spectrum

From `theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/H3S-160GPa/EPW_anharmonic/H3S.a2f`:

```
λ = 1.876
ω_log = 106.9 meV = 1241 K
ω₂ = 137.9 meV = 1600 K
varpi = ω_log/√(λω₂) = 0.663
```

### Critical finding: BCS limit check **FAILS**

Testing Einstein model λ=0.5, ω=50 meV, μ*=0:

```
Tc = 22.9 K
Delta(omega_0) = 3.99 meV (at T=0.03*Tc)
2*Delta/kTc = 4.05
```

**Expected: 3.53. Error: +15%.**

This is λ=0.5, not ultra-weak coupling, but still in the regime where strong-coupling corrections should be modest. The 15% error indicates either:
1. Solver has systematic bias
2. T=0.03*Tc is not low enough (requires T→0 extrapolation)
3. Definition mismatch: Δ(ω₀) ≠ physical gap

### H3S results with compromised baseline

Despite the BCS limit failure, recording H3S results for completeness:

```
mu*=0.13: Tc=210.7 K, 2*Delta/kTc=4.54 (T=0.05*Tc)
mu*=0.15: Tc=202.1 K, 2*Delta/kTc=4.52
mu*=0.20: Tc=183.6 K, 2*Delta/kTc=4.46

Experiment: Tc=203 K, 2*Delta/kTc=3.43
```

Theory overestimates by **32%**. Since BCS limit is already wrong by 15%, the H3S error could be:
- Real physics (non-Migdal, DOS effects, anharmonic corrections)
- Compounded solver systematic

**Cannot proceed to physics interpretation without fixing BCS limit.**

## Profile width analysis

Convergence issue at T≈Tc prevented clean eigenvector extraction in this run. From team-lead's table (Einstein spectrum, independently verified):

```
varpi   alpha_opt   residual
0.3     6.47        0.0046
0.5     4.03        0.0099
1.0     2.32        0.0149
2.0     1.67        0.0046
4.0     1.42        0.0073
8.48    1.31        0.0112
```

H3S varpi≈0.66 would interpolate to α≈3–4. Profile width decreases with increasing varpi (stronger coupling), approaching α→1.2 asymptotically.

## Unresolved issues (blocking)

1. **BCS limit failure**: λ=0.5 μ*=0 gives 4.05 instead of 3.53
   - Need to verify solver test suite
   - Check if Δ(ω₀) is the correct quantity (vs frequency-averaged gap)
   - May need T→0 extrapolation beyond T=0.03*Tc

2. **Profile convergence near Tc**: solver struggles at T=0.995*Tc
   - May need different starting point or mixing parameter
   - Alternative: use T=0.8*Tc where gap is well-formed

3. **Definition of "physical gap"**: tunneling measures Δ at Fermi surface
   - Eliashberg Δ_n varies with Matsubara frequency
   - Is Δ(ω₀) the right comparison? Or should we average?

## Next steps (prioritized)

1. **[CRITICAL] Debug BCS limit**
   - Run solver test suite: `pytest theory/solvers/tests -q`
   - Check if test_eliashberg_iso.py has a 2Δ/kTc test
   - If tests pass, understand why λ=0.5 gives 4.05 not 3.53

2. **Profile as ansatz**
   - If closed-form profile τ_n is accurate, use it as initial guess for nonlinear solver
   - Check if self-consistent solution maintains the profile shape
   - Extract α(T) temperature dependence

3. **Strong-coupling parameter**
   - Plot α(ϖ) vs Tc/ω_log to see if α is a tighter correlation
   - Use KAY/Eliashberg comparison (experimental α²F database)

4. **Carbotte formula comparison**
   - Carbotte RMP 1990 gives empirical 2Δ/kTc(Tc/ω_log)
   - Check if α(ϖ) provides better fit to experimental data

## Why we cannot deliver "H3S 2Δ/kTc = X" yet

The user asked: (1) compute H3S gap ratio, (2) see if closed-form profile reproduces it, (3) check if α(ϖ) is better than Tc/ω_log as strong-coupling parameter.

**Status**:
- (1) **BLOCKED**: BCS limit wrong by 15% ⟹ cannot trust H3S value
- (2) **INCOMPLETE**: profile extraction failed due to convergence
- (3) **DEFERRED**: need (1) and (2) to work first

**Honest assessment**: The numerical infrastructure (Eliashberg solver) has an undiagnosed issue. Either:
- The solver is correct and I'm using it wrong (most likely)
- The solver has a bug in the gap equation
- The definition of "gap" needs clarification

The team-lead's eigenvector profile work (α vs ϖ table) is from **linearized operator at Tc**, which is solid. The disconnect is in going from linearized eigenvector → nonlinear gap at T=0.

## What I DID establish

1. H3S λ=1.876, ω_log=107 meV, varpi=0.66 ✓ [numerical-evidence]
2. Profile width α decreases monotonically with varpi ✓ [from team-lead's table]
3. Solver runs and converges on H3S spectrum ✓
4. BCS limit check reveals 15% systematic error ✓ [blocking finding]

## References

- Experiment: Du/Drozdov et al. (2025) Nature, doi:10.1038/s41586-025-08895-2
  H3S @ 155 GPa: Tc=203 K, 2Δ=60 meV, ratio=3.43

- Solver: `theory/solvers/eliashberg_iso.py` (v0, isotropic Migdal-Eliashberg)

- Spectrum: Mishra et al. 160 GPa anharmonic α²F (EPW calculation)

- Profile: `theory/notes/_scratch_ir/_chk_lemmaC_width.py` (team-lead's verification)

## Scripts generated (all in `theory/notes/_scratch_ir/`)

1. `gap_ratio_h3s.py` — first attempt, T=0.1*Tc
2. `gap_ratio_full_analysis.py` — T→0 extrapolation (partial)
3. `gap_ratio_diagnostic.py` — BCS limit check (hit failure)
4. `gap_ratio_key_findings.py` — final comprehensive run

All scripts set `OMP_NUM_THREADS=1`, use n_matsubara=600, mixing=0.25 as per project conventions.
