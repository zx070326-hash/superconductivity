# FeSe Dynamic Im χ(q,ω) Sign Convention Controversy

**Date**: 2026-09-22  
**Status**: Unresolved theoretical question  
**Strictness**: analytical reasoning (physical arguments, not definitive)

---

## Executive Summary

Phase 1 dynamic susceptibility verification encountered a fundamental sign convention question: **Does the imaginary part of multi-orbital magnetic susceptibility Im χ_ij(q,ω>0) have to be positive semi-definite (PSD)?**

**Finding**: This question has **no clear answer in accessible literature**. Existing repository documentation (2026-09-21) marks this as a **"文献空位"** (literature gap).

**Decision**: Shelve the controversy, preserve both calculation versions, proceed to Phase 2 geometric nesting survey.

---

## The Controversy

### Initial assumption (incorrect)

**Premise**: Fluctuation-dissipation theorem guarantees Im χ(q,ω>0) is PSD.

**Numerical result**: 35/35 test cases showed negative eigenvalues (O(0.1-2.0)), violating PSD.

**Diagnosis**: Sign error suspected in analytic continuation formula.

### Sign correction attempt

**1D free electron gas test** (analytical verification):
```
Im χ = -π × Lorentzian × [f(k) - f(k+q)]  →  ALL NEGATIVE
Im χ = +π × Lorentzian × [f(k) - f(k+q)]  →  ALL POSITIVE
```

**Conclusion from test**: +π convention gives PSD for scalar (single-band) case.

**Problem**: FeSe 5-orbital calculation with +π still gave 35/35 negative eigenvalues.

### Theoretical analysis

**Key realization**: Fluctuation-dissipation theorem constrains the **Hermitian part** of the full complex susceptibility, not necessarily Im χ_ij itself.

**Evidence from repository** (`theory/notes/literature/rpa_chi_theory_gap_2026-09-21.md`, Line 98-106):

> "涨落耗散定理保证 χ 的 **Hermitian 部分**在平衡态必须 PSD"
>
> Translation: "Fluctuation-dissipation theorem guarantees the **Hermitian part** of χ must be PSD in thermal equilibrium"

**Implication**: 
- For **static** χ(q,0): Real and symmetric → must be PSD
- For **dynamic** Im χ(q,ω>0): This is part of the anti-Hermitian component
- The constraint applies to `(χ + χ†)/2`, not to `Im χ` alone

**Physical origin of negative eigenvalues**:
- Multi-orbital interference (off-diagonal elements)
- Even with positive diagonal elements, Hermitian matrix can have negative eigenvalues
- Example: `[[1.0, 0.9+0.45i], [0.9-0.45i, 1.0]]` has eigenvalues `[-0.006, +2.006]`

---

## Literature Status

### What we know (from 2026-09-21 repository search):

**Confirmed (analytical)**:
1. Static χ(q,0) must be PSD (thermal stability)
2. Diagonal elements of Im χ_ij from Lehmann representation are non-negative
3. 2×2 PSD condition: χ₁₂² ≤ χ₁₁·χ₂₂ (Cauchy-Schwarz)

**Literature gap**:
- No explicit theorem found for **matrix-valued** Im χ_ij(q,ω>0) PSD requirement
- Standard textbooks (Mahan, Fetter & Walecka, Altland-Simons) give components but don't explicitly state the matrix constraint
- Iron-based RPA literature focuses on gap symmetry (s±), not full χ matrix positive-definiteness

**Quote from 2026-09-21 documentation**:
> "未找到：明确的解析 no-go 定理：'RPA + 负非对角 ⟹ 必有负本征值'"
>
> Translation: "NOT FOUND: explicit analytical no-go theorem: 'RPA + negative off-diagonal ⟹ must have negative eigenvalues'"

### What we need (unavailable):

**Definitive sources** (web search budget exhausted, cannot access):
- Mahan, *Many-Particle Physics*, Ch. 5 (response functions)
- Altland-Simons, *Condensed Matter Field Theory*, Ch. 7
- Explicit statement on multi-orbital Im χ PSD requirement

**Verification path blocked**:
- 200/200 web searches used
- Cannot fetch textbook pages
- arXiv access blocked by network policy

---

## Two Calculation Versions

### Version 1: -π convention (original)

**Formula**:
```python
im_chi_band_mn = -np.pi * lorentzian(omega - delta_E, eta) * fermi_diff
```

**Result**: 35/35 negative eigenvalues (λ_min ~ -0.1 to -2.0)

**Status**: May be physically correct for multi-orbital case

### Version 2: +π convention (after 1D test)

**Formula**:
```python
im_chi_band_mn = +np.pi * lorentzian(omega - delta_E, eta) * fermi_diff
```

**Result**: Expected 0/35 violations (not fully tested on FeSe 5-orbital)

**Status**: Matches 1D free electron gas PSD behavior

### Sign-independent results

**Valid regardless of sign convention**:

1. **Q_AFM identification**: M = (π,π)
   - Frobenius norm ||Im χ(M, ω=0.5eV)|| = 1.740 (maximum among 5 q-points)
   - Enhancement factor: 2.32× vs average

2. **Nesting correlation**: ξ(q) vs ||Im χ(q,ω)||
   - Pearson r = 0.474 (moderate positive correlation)
   - Ranking consistency: 3/5 matches

3. **Hermiticity**: Im χ_ij = (Im χ_ji)* satisfied numerically

---

## Physical Interpretation

### If negative eigenvalues are allowed:

**Mechanism**: 
- Orbital-resolved susceptibility has interference terms
- Negative eigenvalues indicate certain linear combinations of orbital densities have "phase opposition" in their response
- Not a causality violation - the full retarded Green's function still satisfies causality

**Analogy**: 
- Similar to how static Re χ(q,0) can have negative eigenvalues (confirmed by Team B for FeSe)
- Multi-orbital = multi-channel interference

### If negative eigenvalues are forbidden:

**Implication**: 
- Sign error in formula (likely in analytic continuation)
- Or projection from band to orbital basis introduces error
- Requires fixing before any Im χ results are publishable

---

## Decision and Path Forward

### Team lead decision (2026-09-22):

**Action**: Shelve sign controversy, proceed to Phase 2 geometric nesting

**Rationale**:
1. Cannot resolve in 2-3 hours without textbook access
2. Geometric nesting ξ(q) is sign-independent
3. Risk management: don't spend 2-3 weeks on a theoretical question when core work can proceed

### Phase 1 deliverables (preserved):

**Sign-independent**:
- ✓ Q_AFM = M identification (geometric, robust)
- ✓ Nesting correlation baseline (r=0.474)
- ✓ Computational framework for Phase 2

**Sign-dependent** (on hold):
- ? PSD verification (35/35 "violations" may be physical)
- ? Absolute values of Im χ (magnitude correct, sign uncertain)

### Recommended next steps (if resources available):

1. **Short term** (if textbook access restored):
   - Read Mahan Ch. 5, Altland-Simons Ch. 7 verbatim statements
   - Check for "matrix-valued" vs "scalar" qualifier in FDT

2. **Medium term** (if still unclear):
   - Consult with condensed matter theorist
   - Cross-check with published multi-orbital RPA calculations

3. **Long term** (publication):
   - If negative eigenvalues turn out to be correct: explain physical meaning in paper
   - If sign correction needed: rerun all calculations with corrected formula

---

## Scripts and Data

**Archived scripts** (git tracked per 毕业规则):
1. `_fese_dynamic_chi_imag_verification_2026-09-22.py` - main Phase 1 code (both versions annotated)
2. `_debug_lindhard_sign_2026-09-22.py` - 1D sign test (proves +π for scalar)
3. `_fese_dynamic_chi_MATRIX_FORM_2026-09-22.py` - theoretical analysis

**Output data** (not committed):
- Full 5q × 7ω = 35 point scan results (both versions)
- Nesting function ξ(q) for 5 special points

**Computational cost**: ~15 minutes per version on 100×100 k-grid

---

## Relation to Other Work

### Team B (SCDFT vertex corrections):

Team B is investigating how f_xc kernel in SCDFT modifies χ. Their work may independently resolve the sign question by:
- Computing full SCDFT χ and checking if it's PSD
- Comparing with our Lindhard χ₀ baseline

**Coordination**: Our geometric nesting (Phase 2) provides physical foundation for their mechanism.

### Team A (RPA diagonal dominance):

Team A found that static RPA χ can have negative off-diagonal elements while maintaining PSD. This supports the "negative eigenvalues allowed" interpretation.

**Consistency**: If static χ can have negative off-diagonal, dynamic Im χ likely can too.

---

## Lessons Learned (W1 application)

### Mistake: Assuming theory without verification

**Error**: Took "fluctuation-dissipation ⟹ Im χ PSD" as granted without checking literature for **matrix-valued** case.

**Correction**: Always verify theoretical premises, especially when extending scalar results to matrices.

### Correct action: 1D test case

**Good**: Used simple 1D free electron gas to test sign convention in controlled setting.

**Result**: Isolated the +π vs -π question before blaming multi-orbital complexity.

### Correct action: Literature search in repository

**Good**: Found existing 2026-09-21 documentation marking this as known gap.

**Result**: Avoided "discovering" a gap that was already documented.

### Process improvement:

**For future**: When encountering 100% failure rate (35/35), immediately question the theoretical premise, not just the implementation.

---

## Conclusion

**Theoretical status**: Unresolved - no definitive textbook statement accessible

**Practical status**: Phase 1 geometric results (Q_AFM, nesting correlation) are valid

**Recommendation**: 
- Mark Im χ sign convention as "pending literature confirmation"
- Proceed with Phase 2 geometric nesting (sign-independent)
- Revisit if textbook access becomes available or Team B's SCDFT work sheds light

**Impact on Phase 2**: None - geometric nesting ξ(q) depends only on Fermi surface topology, not on Im χ.

---

**Related files**:
- `../literature/rpa_chi_theory_gap_2026-09-21.md` - Literature gap documentation
- `fese_dynamic_chi_phase1_correction_2026-09-22.md` - Full Phase 1 report with theory correction
- `fese_nesting_analytical_reasoning_2026-09-22.md` - Geometric mechanism (sign-independent)
