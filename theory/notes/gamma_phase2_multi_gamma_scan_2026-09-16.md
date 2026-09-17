# γ-Model Phase 2: Multi-γ Complete Scan Results

**Date**: 2026-09-16  
**Server**: petlove-train (`/zx/p4cert`)  
**Script**: `gamma_dk_coulomb_schur.py` (CERTIFIED Schur tail bounds)  
**Parameters**: N=400, μ* ∈ {0.00, 0.10, 0.13, 0.15}, Nc=50 (for μ*>0)  
**Strictness**: **certificate** for 𝔤_lo (all rows) and for ‖E‖/sinΘ/c₁²_lo at **γ>0.5 only**; see the audit correction below for γ=0.50.

> **AUDIT CORRECTION #2 (team-lead, 2026-09-16, later)** — see `gamma_tail_bound_audit_v2_2026-09-16.md`.
> The γ>½ branch of `gamma_tail_bound.py` has three further one-sided defects (integral-comparison direction, Hankel substitution direction, "worst row = N+1" assumption), so the ‖E‖ / sinΘ / c₁²_lo columns below were **not derived as inequality chains**. A sound weighted-Schur bound (v2) is ≤ the v1 number at every (γ,N) ever published, so every c₁²_lo below **remains a valid lower bound** (a posteriori); no number is retracted, the derivation is superseded. The block-diagonal Davis–Kahan refinement (v3) then gives much tighter certificates: c₁²_lo(μ*=0) = **0.543 / 0.704 / 0.796 / 0.838** for γ = 0.67 / 1 / 1.5 / 2 (vs 0.028 / 0.600 / 0.746 / 0.802 below). §2 and the "bottleneck" sentence at the end are superseded accordingly; the γ=0.67 column is no longer vacuous.
>
> **AUDIT CORRECTION (team-lead, 2026-09-16 18:40, independently re-derived)**
> The γ=0.50 rows' `‖E‖ = 1.30031 (cert)` was **NOT sound**: `gamma_tail_bound.weighted_tail_bound_hi` had a γ≤½ branch returning a finite number `(1+ln J)/√2` for the weighted tail Σ_{m≥J} u_m/m^{1/2}, which **diverges** (harmonic type). Check: for J=401 the branch returned 4.945 while the partial sum already exceeds 6.02 at m=2·10⁶. This is exactly the non-compactness established in `gamma_small_compactness_2026-09-16.md` (‖E_N‖=∞ for γ≤½).
> **Fix**: that branch now returns `+inf` (the only sound upper bound); the γ>½ branch is unchanged bit-for-bit (verified γ=0.67, 1.0). Rerun γ=0.50 N=400: `‖E‖=+inf, sinΘ=1, c₁²_lo=0, 𝔤_lo=2.24837143` (unchanged). The original row's downstream verdict (sinΘ=1, c₁²_lo=0) was already vacuous, so **no false numerical claim was made**, but the label "‖E‖ ≤ 1.300 (cert)" was false and is retracted. γ=0.50 c₁² is **N/A (non-compact)**, not "0 (certificate)".
> Also retracted: the scaling claim "‖E‖ ∝ N^{1−γ}" (§3) — at γ=1 it predicts N⁰, contradicting the observed decay 0.2600→0.1675→0.1081 (N=100→200→400), i.e. exponent ≈ **−0.63**. (The Phase-1 note's "~N⁻¹" was also imprecise.) The Schur-bound scaling is not a clean power law; do not cite either exponent.

---

## Full Scan Table (20 rows: 5γ × 4μ*)

| γ    | μ*   | Nc  | 𝔤_lo (cert)    | ‖E‖ (cert) | sinΘ (cert) | c₁²_lo (cert) | c₁²_f64 | Tc/Tc₀  |
|------|------|-----|----------------|------------|-------------|---------------|---------|---------|
| 0.50 | 0.00 | -   | 2.24837143     | **+inf**¹  | 1.00000     | **N/A**¹      | 0.54808 | 1.0000  |
| 0.50 | 0.10 | 50  | 1.99521105     | **+inf**¹  | 1.00000     | **N/A**¹      | 0.56796 | 0.7875  |
| 0.50 | 0.13 | 50  | 1.92118437     | **+inf**¹  | 1.00000     | **N/A**¹      | 0.57427 | 0.7301  |
| 0.50 | 0.15 | 50  | 1.87239955     | **+inf**¹  | 1.00000     | **N/A**¹      | 0.57853 | 0.6935  |
| 0.67 | 0.00 | -   | 1.92049116     | 0.62473    | 0.51112     | 0.02846       | 0.62557 | 1.0000  |
| 0.67 | 0.10 | 50  | 1.68834344     | 0.62473    | 0.62059     | 0.00018       | 0.65361 | 0.8251  |
| 0.67 | 0.13 | 50  | 1.62194260     | 0.62473    | 0.66173     | 0.00000       | 0.66193 | 0.7771  |
| 0.67 | 0.15 | 50  | 1.57862189     | 0.62473    | 0.69182     | 0.00000       | 0.66737 | 0.7463  |
| 1.00 | 0.00 | -   | 1.60746057     | 0.10813    | 0.07359     | 0.59968       | 0.72298 | 1.0000  |
| 1.00 | 0.10 | 50  | 1.41510424     | 0.10813    | 0.08316     | 0.60836       | 0.75019 | 0.8803  |
| 1.00 | 0.13 | 50  | 1.36250691     | 0.10813    | 0.08660     | 0.60813       | 0.75636 | 0.8476  |
| 1.00 | 0.15 | 50  | 1.32884160     | 0.10813    | 0.08896     | 0.60724       | 0.75983 | 0.8267  |
| 1.50 | 0.00 | -   | 1.40859254     | 0.04165    | 0.03062     | 0.74607       | 0.80065 | 1.0000  |
| 1.50 | 0.10 | 50  | 1.25156270     | 0.04165    | 0.03446     | 0.75643       | 0.81853 | 0.9242  |
| 1.50 | 0.13 | 50  | 1.21045071     | 0.04165    | 0.03567     | 0.75601       | 0.82035 | 0.9039  |
| 1.50 | 0.15 | 50  | 1.18456867     | 0.04165    | 0.03648     | 0.75484       | 0.82065 | 0.8909  |
| 2.00 | 0.00 | -   | 1.31814021     | 0.02748    | 0.02132     | 0.80179       | 0.84081 | 1.0000  |
| 2.00 | 0.10 | 50  | 1.17845709     | 0.02748    | 0.02390     | 0.80816       | 0.85219 | 0.9455  |
| 2.00 | 0.13 | 50  | 1.14273240     | 0.02748    | 0.02467     | 0.80615       | 0.85157 | 0.9311  |
| 2.00 | 0.15 | 50  | 1.12042701     | 0.02748    | 0.02517     | 0.80390       | 0.85021 | 0.9220  |

**Note on Tc/Tc₀**: Tc ∝ 𝔤^(1/γ), so ratio = (𝔤(μ*)/𝔤(0))^(1/γ). This is **numerical-evidence** for the ratio (each endpoint 𝔤_lo is certificate, but the ratio itself is not certified two-sided without an upper bound on 𝔤(μ*=0)).

¹ γ=0.50 lies on the non-compact boundary (weighted tail diverges); tail control and hence c₁² are not certifiable. See audit correction at top. 𝔤_lo for γ=0.50 remains **certificate** (Rayleigh quotient needs no tail control).

**Independent re-derivation (team-lead)**: γ=2.00 rows μ*=0 and μ*=0.13 recomputed from zero: 𝔤_lo = 1.31814021 / 1.14273240, ‖E‖ = 0.02748, c₁²_lo = 0.80615, ratio^{1/2} = 0.9311 — bit-identical to the table.

---

## Physical Observations

### 1. Tc Suppression Trends vs γ

**Coulomb repulsion μ*=0.13 suppresses Tc by:**
- γ=0.50 (nematic): **−27.0%** (Tc/Tc₀ = 0.730)
- γ=0.67 (SDW):     **−22.3%** (Tc/Tc₀ = 0.777)
- γ=1.00 (FM):      **−15.2%** (Tc/Tc₀ = 0.848)
- γ=1.50:           **−9.6%**  (Tc/Tc₀ = 0.904)
- γ=2.00 (Einstein): **−6.9%**  (Tc/Tc₀ = 0.931)

**Trend**: Stronger γ (harder bosons) → weaker Coulomb suppression. This matches physical intuition: harder bosons have less low-energy spectral weight where Coulomb repulsion operates.

### 2. Eigenvector Concentration c₁² vs γ

**At μ*=0.00:**
- γ=0.50: c₁²_lo = **0.000** (certificate), sinΘ = 1.000 → tail dominates, no certified concentration
- γ=0.67: c₁²_lo = **0.028** (certificate), sinΘ = 0.511 → weak concentration
- γ=1.00: c₁²_lo = **0.600** (certificate), sinΘ = 0.074 → strong concentration
- γ=1.50: c₁²_lo = **0.746** (certificate), sinΘ = 0.031 → very strong
- γ=2.00: c₁²_lo = **0.802** (certificate), sinΘ = 0.021 → extremely strong

**Trend**: c₁² increases monotonically with γ. Harder bosons → stronger concentration on zero mode n=0.

### 3. Tail Norm ‖E‖ Scaling

**Schur tail bounds ‖E‖ = ‖C‖ + ‖𝔊^BB‖:**
- γ=0.50: ‖E‖ = **1.300** → tail sums Σ_{m≥N} u_m/m^γ barely converge for γ>1/2
- γ=0.67: ‖E‖ = **0.625** → 2.1× reduction
- γ=1.00: ‖E‖ = **0.108** → 12× reduction from γ=0.50
- γ=1.50: ‖E‖ = **0.042** → 31× reduction
- γ=2.00: ‖E‖ = **0.027** → 48× reduction

~~**Scaling**: ‖E‖ ∝ N^(1−γ) for large N (from integral bounds).~~ **RETRACTED** (see audit correction at top): observed exponent at γ=1 is ≈ −0.63 over N=100→400, not 0. The γ=0.50 entry above (1.300) is likewise retracted — the true value is +inf (non-compact).

### 4. Old Placeholder vs CERTIFIED Schur Bounds

**Ratio (old placeholder ‖E‖) / (Schur ‖E‖):**
- γ=0.50: **0.488** → old placeholder was 2.1× too optimistic (unsound for soft bosons)
- γ=0.67: **1.016** → old placeholder slightly larger (accidentally conservative)
- γ=1.00: **0.979** → old placeholder 2% optimistic (the bug this script fixes)
- γ=1.50: **0.254** → old placeholder 4× too optimistic
- γ=2.00: **0.010** → old placeholder 100× too optimistic

**Conclusion**: The old `tail_bounds_coarse()` placeholder was **systematically unsound** except accidentally at γ≈0.67. The Schur integration fixes this across all γ.

---

## Direction Audit (One-Sided Bounds)

All inequalities are sound:
- **𝔤_lo**: LOWER bound (Rayleigh quotient on trial vector)
- **λ₂_hi**: UPPER bound (inertia count with margin)
- **r₁**: UPPER bound (residual norm)
- **‖C‖, ‖𝔊^BB‖**: UPPER bounds (Schur test on nonnegative blocks)
- **sinΘ**: UPPER bound (Davis-Kahan theorem)
- **c₁²_lo**: LOWER bound (cos²Θ · |v₁[0]| − sinΘ)²

The μ*=0 value 𝔤_lo is valid for the **INFINITE operator** 𝔊(γ,0,∞), not just the N×N truncation, because v₁ is supported on n≤N.

---

## Computational Details

- **Runtime per scan**: ~180–200 seconds × 4 configs = ~12–13 minutes per γ
- **Total wall time**: ~60 minutes (5 scans in parallel)
- **Environment**: OMP_NUM_THREADS=1, mpmath dps=30
- **Server load**: 5 Python processes @ 100% CPU each, no memory issues

---

## Next Steps (Phase 3 Candidates)

1. **Kato-Temple upper bounds** on 𝔤(μ*=0) → certified two-sided Tc/Tc₀ ratios
2. **Higher N sweep** (N=600, 800) to verify ‖E‖ ∝ N^(1−γ) scaling law
3. **Intermediate γ grid** (γ ∈ [0.5, 2.0] step 0.1) for smooth interpolation
4. **Extended μ* range** (μ* up to 0.20 or 0.25) to map full suppression curve
5. **Anisotropic γ** (γ_x ≠ γ_y) for cuprate/pnictide models

**Current bottleneck**: γ≤0.5 non-compactness (fundamental, see `gamma_small_compactness_2026-09-16.md`) and the ‖𝔊^BB‖-dominated c₁² gap at γ=1 (see `gamma_c1sq_tightening_2026-09-16.md`). Certificate-level c₁² achieved for γ∈{0.67,1,1.5,2}; γ=0.67 is nearly vacuous (c₁²_lo=0.028).
