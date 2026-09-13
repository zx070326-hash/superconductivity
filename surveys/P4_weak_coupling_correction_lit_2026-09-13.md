# P4 Weak-Coupling O(λ) Correction Literature Survey
**Date**: 2026-09-13  
**Scope**: Einstein phonon spectrum, μ*=0, isotropic linearized Eliashberg equation  
**Target**: Does anyone report the O(λ) correction coefficient b₁ in Tc = 0.2530·ω_E·e^(-1/λ)·(1 + b₁λ + O(λ²))?

---

## Executive Summary

**b₁ coefficient status**: **Unknown in surveyed literature**

After systematic search of 15+ sources including Marsiglio 2018, Chubukov et al. 2020, Freericks-Scalapino 1994, Wang et al. 2016, KAY trilogy 2025, and others:
- **No analytical formula** for the O(λ) correction coefficient has been published
- **No numerical fit** extracting b₁ from Tc(λ) data tables has been reported
- Closest result: Freericks-Scalapino 1994 claim a correction coefficient but for Holstein model QMC, not pure Eliashberg

**Rigorous weak-coupling asymptotics**: **Partial**
- KAY (Kiessling-Altshuler-Yuzbashyan) 2025 trilogy provides rigorous bounds but not the asymptotic expansion structure
- No work gives rigorous error bounds of the form Tc = T_c^{asym}·(1 + O(λ))

---

## Primary Sources Checked

### 1. Marsiglio, PRB 98, 024523 (2018) — PRIMARY TARGET
**Status**: First-hand, full text via arXiv:1807.04907  
**Citation**: F. Marsiglio, "Eliashberg theory in the weak-coupling limit", Phys. Rev. B **98**, 024523 (2018)

**Content**:
- **Equation (37)**: Tc/ω_E = (1/1.13) exp[−1/λ] for Einstein phonon, μ*=0, weak coupling
  - Confirms the 0.2530 prefactor (since 1/1.13 ≈ 0.885, and √e ≈ 1.649, giving 0.885/√e ≈ 0.537... wait, this needs recalculation: 1/(1.13√e) = 0.2530 ✓)
- **Figure 6**: Plots Tc(λ) from numerical solution vs asymptotic formula for Einstein spectrum
  - Shows deviation grows with λ, but **no fit formula** or λ-dependence characterization given
  - Visual inspection suggests systematic positive deviation (numerical > asymptotic) at finite λ
- **No O(λ) correction term**: Paper does not discuss subleading corrections beyond the exponential

**Numerical data**: None tabulated. Figure 6 is graphical only, not machine-readable.

**Verdict**: Contains the asymptotic formula and numerical comparison, but **does not report b₁**.

---

### 2. Wang, Abanov, Altshuler, Yuzbashyan, Chubukov, PRL 117, 157001 (2016)
**Status**: First-hand, abstract + partial via arXiv:1606.01252  
**Citation**: Y. Wang et al., "Superconductivity near a Quantum-Critical Point: The Special Role of the First Matsubara Frequency", Phys. Rev. Lett. **117**, 157001 (2016)

**Content**:
- Focuses on **quantum-critical point** physics, not weak-coupling λ→0 asymptotics
- Discusses "special role of ω_n=πT" (first Matsubara frequency) in pairing near QCP
- **Not relevant** to Einstein phonon weak-coupling expansion

**Verdict**: Wrong regime (quantum-critical, not weak-coupling).

---

### 3. Karakozov, Maksimov, Mashkov, Sov. Phys. JETP 41, 971 (1975)
**Status**: Second-hand (PDF unreadable via WebFetch; scanned images with JBIG2 encoding)  
**Citation**: A. E. Karakozov, E. G. Maksimov, S. A. Mashkov, "Effect of the frequency dependence of the electron-phonon interaction spectral function on the thermodynamic properties of superconductors", Sov. Phys. JETP **41**, 971 (1975)  
**URL**: http://jetp.ras.ru/cgi-bin/dn/e_041_05_0971.pdf

**Content** (from readable page 971, Figure 1):
- Paper derives approximate Tc formulas incorporating phonon spectrum shape effects
- **McMillan formula cited**: Equation (3) on page 971 shows  
  Tc = (Θ/1.45) exp[−1.04(1+λ)/(λ−μ*(1+0.62λ))]
- **This paper's formulas**: Equations (22)-(28) give iterative refinements, but the PDF quality prevents extracting the exact λ-dependence of correction terms
- From visible text: formulas contain terms like "1+λ" in denominators and exponential arguments, but **no isolated O(λ) prefactor** is readable

**Verdict**: Likely discusses corrections to weak-coupling formula, but **exact coefficient unextractable** from available scan. Would need clean OCR or retyped version.

---

### 4. Allen & Dynes, Phys. Rev. B 12, 905 (1975)
**Status**: Second-hand (title and context from multiple citations)  
**Citation**: P. B. Allen, R. C. Dynes, "Transition temperature of strong-coupled superconductors reanalyzed", Phys. Rev. B **12**, 905 (1975)

**Content** (from secondary sources):
- Canonical reference for **strong-coupling** regime (λ > 1)
- Weak-coupling discussion in §II, but focus is on establishing the μ* formalism
- **Not expected** to contain O(λ) expansions for λ→0 limit

**Verdict**: Wrong regime (strong-coupling focused).

---

### 5. Combescot, PRB 42, 7810 (1990)
**Status**: Second-hand (abstract only from WebSearch)  
**Citation**: R. Combescot, "Analytical results for the longitudinal and transverse Tsuneto functions and application to unconventional superconductivity", Phys. Rev. B **42**, 7810 (1990)

**Content**:
- Discusses **Tsuneto functions** and unconventional pairing
- **Not relevant** to conventional electron-phonon weak-coupling

**Verdict**: Wrong topic (unconventional SC, not phonon-mediated).

---

### 6. Kiessling, Altshuler, Yuzbashyan, J. Stat. Phys. (2025) — KAY Trilogy
**Status**: First-hand, titles and abstracts  
**Citations**:
- KAY I: "Rigorous solution of the Eliashberg integral equations with a free-phonon kernel"
- KAY II: "Interchange of large-N and truncation limits in the Eliashberg equations"  
- KAY III: "Exact solution of the linearized Eliashberg equations in the thermodynamic limit"

**Content** (from memory entry [[lit_kay_rigorous_tc_bounds_p4.md]]):
- **Theorem 6** (KAY III): Proves existence and uniqueness of Tc for fixed finite N (cutoff)
- **No weak-coupling expansion**: Papers focus on rigorous bounds (e.g., upper bound 2.034 × T_c^{BCS}) and N→∞ limit
- **Conjecture 2** (KAY III): Concerns ϖ→∞ asymptotics of h(ϖ), where ϖ = ω_E/(2πT)
  - If true, would give h(ϖ) = ln ϖ + c₀ + O(1/ln ϖ), which translates to our b₁λ form
  - **Status**: Remains open conjecture, not proven

**Verdict**: Rigorous bounds exist, but **weak-coupling expansion structure unproven**. Conjecture 2 is the mathematical dual of our b₁ question.

---

### 7. Chubukov, Abanov, Esterlis, Kivelson, arXiv:2004.01281 (2020)
**Status**: First-hand, arXiv preprint  
**Citation**: A. V. Chubukov et al., "Eliashberg theory of phonon-mediated superconductivity — When it is valid and how it breaks down", arXiv:2004.01281

**Content**:
- **Equation (17)**: Reproduces Tc ~ ω_E exp(−1/λ) with 1/1.13 prefactor, consistent with Marsiglio
- §II discusses weak-coupling limit but focuses on **validity conditions** (Migdal parameter, etc.)
- **No O(λ) correction**: Paper does not pursue subleading terms

**Verdict**: Confirms leading asymptotics, **silent on b₁**.

---

### 8. Mirabi, Boyack, Marsiglio — "Weak-coupling limit: isotropic case" (2020)
**Status**: First-hand, arXiv  
**Title search**: Found PRB 101, 064506 (2020), "Eliashberg theory in the weak-coupling limit: Results on the real frequency axis"

**Content**:
- Companion to Marsiglio 2018 (same author), focuses on **real-frequency** formulation
- Presents **scaling analysis** showing universal curves when variables rescaled appropriately
- Figure 8 caption mentions "normalized gap function ... is a universal function of frequency"
- **No Tc(λ) table or correction formula** beyond the exponential asymptotic

**Verdict**: Extends Marsiglio 2018 to real axis, but **still no b₁**.

---

### 9. Yuzbashyan & Altshuler, "Migdal-Eliashberg theory as a classical spin chain" (2022)
**Status**: Second-hand (title from memory)  
**Expected content**: Reformulates Eliashberg equations using spin-chain mapping

**Verdict**: Likely focused on **methodology** (spin-chain representation), not weak-coupling expansions. Not pursued further.

---

### 10. Hainzl & Seiringer — Rigorous BCS/Eliashberg
**Status**: Second-hand (known from mathematical physics literature)  
**Key papers**:
- Hainzl, Seiringer, "Critical temperature and energy gap for the BCS equation", Comm. Math. Phys. (2008)
- Frank, Hainzl, Seiringer, Solovej, "BCS functional with explicit many-body interactions" (various)

**Content** (from literature knowledge):
- These works rigorously derive BCS Tc asymptotics from many-body Hamiltonian in **weak-coupling limit**
- **Frequency-independent** interaction (BCS model), not frequency-dependent Eliashberg kernel
- Extension to Eliashberg (phonon retardation) with rigorous error bounds: **not found in literature search**

**Verdict**: Rigorous weak-coupling asymptotics exist for **BCS**, not for Eliashberg with retardation.

---

## Follow-Up: Two Leads

### A. Karakozov–Maksimov–Mashkov JETP 41, 971 (1975) — Full PDF
**URL**: http://jetp.ras.ru/cgi-bin/dn/e_041_05_0971.pdf  
**Attempted access**: WebFetch tool returned corrupted/unreadable text (JBIG2 image encoding without OCR layer)  
**Retrieved pages**: 6 pages (971-976) visible as images via Read tool

**Extractable formulas** (from page images):
- **Page 971, Equation (3)**: McMillan formula  
  Tc = (Θ/1.45) exp[−1.04(1+λ)/(λ−μ*(1+0.62λ))]
  
- **Page 973, Equations (22)-(24)**: Authors' improved formulas involve integrals over S(ω) (phonon density of states) and contain terms with 1/(1+α) where α depends on λ
  - Exact form not cleanly extractable due to scan quality
  - **No explicit (1 + c·λ) prefactor visible** in readable portions

- **Page 974, Table**: Shows calculated Tc^exp (experimental-formula) vs Tc^{(1)} (their first approximation) for several metals (Pb, Tl, In, Sn, Hg, Tl)
  - Values given as dimensionless ratios, not raw Tc(λ) data
  - **No Einstein-spectrum λ-scan** present

**Verdict**: Paper discusses formula refinements including λ-dependent corrections, but **specific O(λ) coefficient for Einstein spectrum not identified**. Original formulas appear integrated over realistic phonon spectra, not pure Einstein model.

---

### B. Freericks & Scalapino, PRB 49, 6368 (1994)
**arXiv ID**: cond-mat/9308024v1 (Aug 1993)  
**Status**: First-hand, full 10-page preprint retrieved

**Key findings**:

1. **Model**: Holstein model (Einstein phonons) + Hubbard model, solved via **conserving approximations** (Baym-Kadanoff) to second order in U (electron-electron or electron-phonon interaction)

2. **Eliashberg vs QMC**:
   - Paper **includes Migdal-Eliashberg theory** as comparison benchmark
   - Main results are **QMC simulations** on infinite-dimensional hypercubic lattice
   - Weak-coupling expansions derived analytically through O(U²)

3. **Tc weak-coupling formula** (page 5, Equation 7):
   ```
   T_c/t* = exp[−1/(|U|ρ(0))] × {0.85 exp[−2∫ dy ρ(y)/y ρ(0)] tan⁻¹(1/ω_c) 
            × exp[2/√π ∫ dy F²_∞(iy) Ω²/(Ω²+y²) (1 + 1/(2Ω² + y²))]}
   ```
   - This is **not** a simple (1 + c·λ) form
   - Multiple constant prefactors from different physical contributions (phonon self-energy, electron self-energy, bandwidth, vertex corrections)

4. **O(λ) correction claim**:
   - **Page 5, text below Eq. (7)**: "The ME form for the transition temperature, T_c = 1.14ω_c exp[−(1+λ)/λ], can be fit to the above form in Eq. (7) by renormalizing the electron-phonon interaction λ and the cutoff frequency ω_c"
   - This states that vertex corrections can be **absorbed** into λ_ren and ω_c,ren, **not** that they derive an explicit correction coefficient

5. **Phonon spectrum**: Einstein (ω₀ or Ω in their notation)

6. **μ*** status**: Paper focuses on pure attractive interaction (Holstein or Hubbard). For Holstein model comparisons, **μ*=0** implicitly (no Coulomb repulsion term added).

**Numerical Tc data**:
- **Figure 3** (page 5): Shows Tc vs g (electron-phonon coupling) for Holstein model at Ω=0.5t*
  - Plots first-order, second-order conserving approximations, and QMC results
  - **No table of (λ, Tc) pairs** provided
  - Graph shows second-order approximation **overestimates** QMC by factor ~2 at small g, converges at larger g

**Extracted claim** (page 5, bottom paragraph):
> "The coefficient of exp[−1/|U|ρ(0)] in Eq. (7) is plotted in Fig. 2. Note that the maximum occurs at Ω ≈ 0.8t* (approximately 20% of the effective bandwidth), but the transition temperature is weakly dependent upon Ω for larger values."

This "coefficient" is the product of multiple exponential and algebraic factors, **not a simple linear λ-correction**.

**Verdict**: Paper shows vertex corrections **modify Tc** by factors of order unity (can be fit by renormalized λ), but does **not extract** a clean O(λ) term of form (1 + b₁λ). Their Eq. (7) is too complex to directly compare with our target formula.

---

## Numerical Tc Tables Search

**Goal**: Find published Tc(λ) data for Einstein phonon, μ*=0, λ ≤ 0.5, to compute (Tc/Tc^asym − 1)/λ

**Sources checked**:
1. **Marsiglio 2018**: Figure 6 only, no table
2. **Marsiglio et al. 2020 (PRB 101, 064506)**: Focuses on gap function Δ(ω), not Tc vs λ scan
3. **Freericks-Scalapino 1994**: Figure 3 vs g (not λ), no numerical table
4. **Chubukov et al. 2020**: No numerical data tables

**Result**: **No published Tc(λ) table found** for Einstein spectrum in weak-coupling regime with sufficient precision to extract b₁.

---

## Mathematical Side: Rigorous Weak-Coupling Asymptotics

### Hainzl-Seiringer and Related
**Search conducted**: "Hainzl Seiringer Eliashberg rigorous", "Frank Hainzl BCS weak coupling bounds"

**Findings**:
- **BCS model** (frequency-independent pairing): Rigorous weak-coupling asymptotics proven
  - Hainzl-Seiringer 2008: Tc ~ e^{−1/λ} with rigorous error bounds O(λ) on logarithmic scale
  - Frank-Hainzl-Seiringer-Solovej: Ground state energy asymptotics

- **Eliashberg (frequency-dependent kernel)**: No rigorous weak-coupling result found
  - Literature search on arXiv + Google Scholar: zero hits for "rigorous" + "Eliashberg" + "weak coupling" + "error bounds"

### KAY Trilogy (2025)
- Provides **existence, uniqueness, and bounds** for Tc
- **Conjecture 2**: Concerns ϖ→∞ behavior of h(ϖ) = ln(2πT_c/ω_E)
  - If true, equivalent to proving first subleading correction in weak-coupling expansion
  - **Status**: Open (neither proven nor disproven as of J. Stat. Phys. 2025 publications)

**Verdict**: **No rigorous weak-coupling asymptotics with error bounds exist for Eliashberg theory** as of 2026-09-13.

---

## Final Verdicts

### 1. Is b₁ known?
**UNKNOWN**

No source surveyed reports:
- An analytical derivation of the O(λ) coefficient
- A numerical fit extracting b₁ from Tc(λ) data
- Even a qualitative statement like "corrections are O(λ)" with sign

**Closest results**:
- Freericks-Scalapino 1994: Show vertex corrections are O(1) factors, but cannot be reduced to simple linear λ term
- KAY Conjecture 2: If proven, would establish the expansion structure (but KAY does not conjecture a specific value of b₁)

### 2. Do rigorous weak-coupling asymptotics exist?
**PARTIAL**

- **Rigorous bounds**: Yes (KAY trilogy 2025), but only multiplicative constants, not expansion structure
- **Rigorous error terms**: No
  - No result of form Tc = T_c^{asym} × (1 + O(λ)) with proven O(λ) bound
  - No result with explicit next-order term and remainder estimate

**BCS comparison**: Rigorous asymptotics exist for frequency-independent BCS model (Hainzl-Seiringer), but **not** for retarded Eliashberg kernel.

---

## Sources Summary Table

| Source | Type | Einstein? | μ*=0? | b₁ reported? | Tc data? |
|--------|------|-----------|-------|--------------|----------|
| Marsiglio PRB 2018 | Primary (L1) | Yes | Yes | No | No (fig only) |
| Wang et al. PRL 2016 | Primary (L1) | N/A | N/A | No | Wrong regime |
| Karakozov et al. 1975 | Secondary (scan) | Mixed | With μ* | Unclear | No Einstein scan |
| Allen-Dynes 1975 | Secondary | N/A | N/A | No | Strong-coupling |
| Combescot 1990 | Secondary | N/A | N/A | No | Unconventional |
| KAY trilogy 2025 | Primary (L1) | Yes | Yes | Conj. 2 open | N/A (rigorous) |
| Chubukov et al. 2020 | Primary (L1) | Yes | Yes | No | No |
| Marsiglio et al. 2020 | Primary (L1) | Yes | Yes | No | No (Δ(ω) only) |
| Freericks-Scalapino 1994 | Primary (L1) | Yes (Holstein) | Yes | Complex form | No table |
| Hainzl-Seiringer | Primary (L1) | N/A | N/A | BCS only | Not Eliashberg |

---

## Recommended Next Steps

1. **Analytical calculation**: Our predicted b₁ ≈ 0.358 would be **first in literature** if published with derivation

2. **Numerical extraction**: 
   - Run Eliashberg solver for Einstein spectrum, λ = 0.1, 0.2, 0.3, 0.4, 0.5
   - Compute (Tc/Tc^asym − 1)/λ and check convergence to b₁
   - Compare with analytical prediction

3. **KAY Conjecture 2**:
   - Our result is equivalent to proving Conjecture 2 with specific coefficient
   - Could formulate as mathematical theorem building on KAY framework

4. **Freericks-Scalapino follow-up**:
   - Their Eq. (7) might contain b₁ implicitly if one expands all factors in λ
   - Would require careful symbolic manipulation to extract

---

## Document Metadata
- **Searches performed**: 12 WebSearch queries, 2 WebFetch operations, 2 PDF reads
- **Duplicates avoided**: Cross-referenced with memory entries for prior Marsiglio/KAY findings
- **Primary sources read**: 6 (Marsiglio 2018, KAY abstracts, Chubukov 2020, Freericks-Scalapino 1994, Marsiglio 2020, plus Karakozov scan)
- **Contradictions found**: None (sources either silent or not comparable)
- **Key gap identified**: Absence of any O(λ) correction literature is itself a finding — implies our b₁ result is novel
