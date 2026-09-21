# RPA Diagonal Dominance Literature Review

**Date**: 2026-09-22  
**Task**: Literature review on susceptibility matrix structure in RPA  
**Research questions**:
1. Does Lindhard χ₀ have diagonal dominance at the Fermi surface?
2. What do field theory textbooks say about χ₀ structure?
3. Are there results about matrix structure preservation under Dyson resummation?
4. Any mention of "diagonal dominance" or "Gershgorin" in RPA context?

---

## Executive Summary

**Key finding**: No explicit diagonal dominance theorem found in the literature for RPA susceptibility matrices in orbital space.

**Evidence landscape**:
- Lindhard χ₀ structure is well-studied, but **intraband vs interband** distinction is better established than diagonal dominance in orbital basis
- RPA matrix inversion procedures are standard, but focus on **convergence with basis size** rather than diagonal dominance
- Local field effects and screening show that **off-diagonal terms are significant** and must be retained for accuracy
- **Gershgorin theorem not mentioned** in RPA susceptibility context (search returned general Gershgorin material only)

**Assessment**: Diagonal dominance preservation is **likely NOT a standard result**. The main obstacle is that off-diagonal (interorbital) susceptibility elements can be substantial, especially near Fermi surface hot spots with strong orbital mixing.

---

## 1. Lindhard Function Structure

### 1.1 Basic Properties

The Lindhard function χ₀(q,ω) describes density-density response of non-interacting fermions ([Lindhard theory - Wikipedia](https://en.wikipedia.org/wiki/Lindhard_theory)). It can be obtained in closed form and is fundamental to RPA.

**Key structure**:
- At q=0: connects to thermodynamic compressibility via sum rule ([arXiv:cond-mat/0209059](https://arxiv.org/pdf/cond-mat/0209059))
- Exhibits branch points related to Fermi surface topology ([arXiv:cond-mat/0604496](https://arxiv.org/html/cond-mat/0604496v1))
- Shows singular behavior at wavevectors matching Fermi surface nesting ([APL 120, 111901](https://pubs.aip.org/aip/apl/article/120/11/111901/2833165/Fermi-surface-nesting-and-the-Lindhard-response))

### 1.2 Intraband vs Interband Contributions

**Critical distinction** ([arXiv:2212.03192](https://arxiv.org/html/2212.03192v1), [arXiv:2504.03882](https://arxiv.org/html/2504.03882v1)):
- **Intraband** contributions dominate at small q (long wavelength limit)
- **Interband** contributions become important for:
  - Finite momentum transfers
  - Optical response
  - Systems with strong orbital mixing near Fermi surface

In flat band systems, intraband conductivity **vanishes entirely**, leaving only interband contributions. This suggests that "diagonal dominance" is basis-dependent and not universal.

### 1.3 Multiband Structure

**Important caveat** ([arXiv:1207.1106](https://arxiv.org/html/1207.1106v3)):

> "The irreducible susceptibility is **not simply proportional to the Lindhard function** because it's affected by matrix elements involving wavefunctions. Both Fermi surface nesting **and** orbital overlap must be considered together."

This means that even if single-band Lindhard has favorable structure, multiband generalization introduces **orbital matrix elements** that can destroy simple dominance patterns.

---

## 2. RPA Matrix Structure

### 2.1 Matrix Formulation

The multiorbital RPA uses matrix formulation ([arXiv:1606.06189](https://arxiv.org/html/1606.06189v2)):
- Treats spin, charge, and orbital fluctuations systematically
- Matrix structure couples different orbital channels
- Off-diagonal terms represent inter-orbital interactions

**Key equation**: χ = χ₀ + χ₀ V χ (Dyson series)

where χ₀ is bare susceptibility, V is interaction kernel, χ is dressed susceptibility.

### 2.2 Dielectric Matrix Inversion

Standard RPA procedure ([PRB 12, 1487 (1975)](http://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.12.1487/fulltext)):
1. Calculate full RPA dielectric matrix ε(G,G') with multiple reciprocal lattice vectors
2. **Invert the matrix**: ε⁻¹(G,G')
3. Extract macroscopic dielectric function from **first diagonal element** of inverse

**Key observation**: The procedure requires **full matrix inversion**, suggesting off-diagonal elements are non-negligible. If diagonal dominance held, simpler approximations would suffice.

### 2.3 Local Field Effects

Local field corrections describe how particle correlations modify the restoring force ([arXiv:nucl-th/0003044](https://arxiv.org/html/nucl-th/0003044v1), [arXiv:cond-mat/0104229](https://arxiv.org/html/cond-mat/0104229v5)):

- **Off-diagonal terms** in ε(G,G') represent local field effects
- These are **essential** for accurate dielectric response in real materials
- Diamond, silicon studies show local fields significantly affect results

**Implication**: Diagonal approximation (neglecting off-diagonal ε) gives **qualitatively wrong** results. This suggests diagonal dominance is weak or absent.

---

## 3. Matrix Structure Preservation Under Dyson Resummation

### 3.1 Matrix Dyson Equation

Searched for structure preservation theorems in context of matrix Dyson equation (MDE).

**Findings** ([arXiv:1604.08188](https://arxiv.org/html/1604.08188v5), [arXiv:2402.16376](https://arxiv.org/html/2402.16376v2)):
- Stability analysis of MDE is well-developed (Ajanki-Erdős-Krüger)
- **Spectral dominance approach** exists for large random matrices
- BUT: these results are for **random matrix theory**, not RPA susceptibility

**No explicit theorem found** stating: "If χ₀ is diagonally dominant, then χ = (1 - χ₀V)⁻¹χ₀ is also diagonally dominant."

### 3.2 Why Structure Might Not Preserve

The Dyson series χ = χ₀(1 + Vχ₀ + Vχ₀Vχ₀ + ...) involves **matrix products**. Even if χ₀ is diagonally dominant, Vχ₀V can amplify off-diagonal elements, especially if V itself has large off-diagonal terms (as Coulomb interaction does in orbital basis).

---

## 4. Gershgorin Circle Theorem in RPA Context

**Search results**: No papers found applying Gershgorin theorem to RPA susceptibility matrices.

General Gershgorin content found ([Wikipedia](https://en.wikipedia.org/wiki/Gershgorin_circle_theorem), [arXiv:1605.07239](https://arxiv.org/html/1605.07239v2)), but **no connection to RPA** in electronic structure literature.

**Interpretation**: This absence suggests either:
1. Diagonal dominance doesn't hold, so Gershgorin is not applicable
2. The RPA community hasn't explored this angle (gap in literature)

Given that local field effects are known to be important (Section 2.3), option (1) seems more likely.

---

## 5. Related Theoretical Frameworks

### 5.1 Kukkonen-Overhauser Method

The KO approach uses **local field factors** G±(q) to capture exchange-correlation beyond RPA ([arXiv:cond-mat/9607178](https://arxiv.org/html/cond-mat/9607178v1), [arXiv:cond-mat/0211442](https://arxiv.org/html/cond-mat/0211442v1)):

χ⁻¹ = χ₀⁻¹ - V[1 + G(q)]

The local field factor G(q) is wavevector-dependent and represents **non-local correlations**. Its q-dependence shows that simple diagonal approximations fail.

### 5.2 Quantum Geometry Contributions

Recent work shows magnetic susceptibility decomposes into band dispersion + quantum geometric contributions ([arXiv:2602.14511](https://arxiv.org/pdf/2602.14511.pdf)):

> "In materials like LaFeAsO, dominant magnetic fluctuations originate from the **geometric contribution** rather than traditional Fermi surface nesting."

This suggests susceptibility structure is controlled by quantum geometry, not just diagonal density of states.

---

## 6. Textbook Content (Attempted)

### 6.1 Altland-Simons "Condensed Matter Field Theory"

**Chapter 7: Response Functions** ([Cambridge Core](https://www.cambridge.org/core/books/condensed-matter-field-theory/response-functions/8BC9D6C9F169FDD218A9B7878370A85F))

Content confirmed to cover Lindhard function. Unfortunately, full text not accessible through WebFetch (enterprise restrictions).

**Known content** (from secondary sources):
- Surveys response function formalism
- Lindhard function as key example of density-density response
- Connects correlation functions to experimental observables

**Speculation**: Standard textbook treatment likely covers single-band Lindhard in detail but may not address multiband diagonal dominance question.

### 6.2 Mahan "Many-Particle Physics"

Referenced multiple times in literature ([arXiv:1208.3210](https://arxiv.org/html/1208.3210v1)) as standard reference for Lindhard susceptibility.

**Content**: Lindhard function as baseline for plasmons, screening, correlation effects in interacting systems.

**Not accessible**: Full chapter content not retrieved.

---

## 7. Computational Evidence

### 7.1 Convergence Requirements

Modern RPA implementations ([arXiv:2307.12029](https://arxiv.org/html/2307.12029v1)):
- Sub-quadratic scaling achieved through **locality and sparsity** patterns
- Can handle 1000+ atom systems

**Key point**: Computational efficiency comes from **real-space locality** (exponential decay with distance), NOT from diagonal dominance in momentum/orbital space.

### 7.2 Basis Set Convergence

RPA correlation energies converge slowly with plane-wave cutoff ([arXiv:0902.0889](https://arxiv.org/html/0902.0889v2)). This **slow convergence** suggests that many reciprocal lattice vectors (off-diagonal blocks in ε(G,G')) contribute significantly.

---

## 8. Assessment: Is Diagonal Dominance Provable?

**Answer: Unlikely in general orbital basis.**

### 8.1 Evidence Against

1. **Local field effects are essential**: Diamond/silicon studies show diagonal approximation fails
2. **Orbital matrix elements matter**: Irreducible susceptibility ≠ bare Lindhard in multiband case
3. **Interband transitions**: Can dominate over intraband at finite q
4. **No theorem in literature**: 40+ years of RPA calculations, but no diagonal dominance cited
5. **Full matrix inversion required**: Standard procedure doesn't exploit diagonal structure

### 8.2 Where It Might Hold

Possible special cases:
- **q→0 limit**: Intraband dominates, but this is band-diagonal not orbital-diagonal
- **Well-separated bands**: If Fermi surface involves only one band, other bands give small interband χ₀
- **Specific basis choices**: Natural orbitals or Wannier functions might enhance diagonal character

### 8.3 Main Obstacle

The **Coulomb interaction V** is itself **not diagonally dominant** in orbital basis (off-diagonal Coulomb matrix elements are substantial). Even if bare χ₀ had favorable structure, the Dyson series χ = χ₀ + χ₀Vχ₀ + ... would amplify off-diagonals through matrix products.

---

## 9. Alternative Approaches

If diagonal dominance proof is blocked, what might work instead?

### 9.1 Perron-Frobenius / M-Matrix Theory

Look for **positivity structure** rather than diagonal dominance:
- Does χ or χ⁻¹ have M-matrix structure (off-diagonals ≤ 0)?
- Analogy: our Theorem GF uses positivity, not dominance

### 9.2 Resolvent Bounds

Use **spectral bounds** instead of element-wise dominance:
- ‖χ‖ ≤ ‖χ₀‖ · something involving ‖V‖
- Related to Matrix Dyson Equation stability analysis

### 9.3 Physical Constraints

Use sum rules and Ward identities:
- Compressibility sum rule constrains q→0
- Gauge invariance constrains current-current response
- Maybe these imply useful inequalities?

---

## Bibliography

### Core RPA Theory
- [Altland & Simons, "Condensed Matter Field Theory", Ch.7](https://www.cambridge.org/core/books/condensed-matter-field-theory/response-functions/8BC9D6C9F169FDD218A9B7878370A85F) — Response functions (not fully accessible)
- [Lindhard Theory - Wikipedia](https://en.wikipedia.org/wiki/Lindhard_theory) — Basic formalism
- [arXiv:1211.6873](https://arxiv.org/html/1211.6873v2) — RPA applied to solids, molecules, graphene interfaces

### Multiorbital RPA
- [arXiv:1606.06189](https://arxiv.org/html/1606.06189v2) — Matrix RPA for multiorbital Hubbard model
- [arXiv:1212.0210](https://arxiv.org/html/1212.0210v1) — Two-orbital model for BiS₂ superconductors
- [arXiv:1301.3254](https://arxiv.org/html/1301.3254v1) — Two-particle self-consistent method for multiorbital systems

### Matrix Structure & Local Fields
- [PRB 12, 1487 (1975)](http://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.12.1487/fulltext) — RPA dielectric matrix for diamond with local fields
- [arXiv:nucl-th/0003044](https://arxiv.org/html/nucl-th/0003044v1) — Momentum conservation and local field corrections
- [arXiv:cond-mat/0104229](https://arxiv.org/html/cond-mat/0104229v5) — Dynamical local field and sum rules

### Lindhard Function Details
- [arXiv:1111.5337](https://arxiv.org/html/1111.5337v1) — Lindhard function in d dimensions
- [arXiv:1208.3210](https://arxiv.org/html/1208.3210v1) — Lindhard and RPA in cuprates
- [arXiv:cond-mat/9806170](https://arxiv.org/pdf/cond-mat/9806170v1) — Charge/spin fluctuations in d-wave SC

### Intraband/Interband Structure
- [arXiv:2212.03192](https://arxiv.org/html/2212.03192v1) — Conductivity in flat bands (intraband vanishes)
- [arXiv:2504.03882](https://arxiv.org/html/2504.03882v1) — Quantum geometry in flat bands
- [arXiv:2602.14511](https://arxiv.org/pdf/2602.14511.pdf) — Magnetic fluctuations from quantum geometry

### Matrix Dyson Equation
- [arXiv:1604.08188](https://arxiv.org/html/1604.08188v5) — Stability of matrix Dyson equation
- [arXiv:2402.16376](https://arxiv.org/html/2402.16376v2) — Spectral dominance approach to large random matrices

### Kukkonen-Overhauser & Local Fields
- [arXiv:cond-mat/9607178](https://arxiv.org/html/cond-mat/9607178v1) — Charge and spin response with local field factors
- [arXiv:cond-mat/0211442](https://arxiv.org/html/cond-mat/0211442v1) — Local fields in 2D electron gas
- [arXiv:2607.20403](https://arxiv.org/pdf/2607.20403) — Exact exchange local field factor

### Sum Rules & Long-Wavelength Limits
- [arXiv:cond-mat/0209059](https://arxiv.org/pdf/cond-mat/0209059) — Compressibility sum rule
- [arXiv:1301.2568](https://arxiv.org/html/1301.2568v2) — Constraints on response theories
- [Wikipedia: Thomas-Fermi screening](https://en.wikipedia.org/wiki/Thomas%E2%80%93Fermi_screening) — q→0 limit

### Computational RPA
- [arXiv:2307.12029](https://arxiv.org/html/2307.12029v1) — Sub-quadratic scaling RPA for periodic systems
- [arXiv:0902.0889](https://arxiv.org/html/0902.0889v2) — Efficient RPA correlation energy calculation
- [VASP Wiki: ACFDT/RPA](https://www.vasp.at/wiki/index.php/ACFDT/RPA_calculations) — Practical implementation

### Fermi Surface & Nesting
- [APL 120, 111901](https://pubs.aip.org/aip/apl/article/120/11/111901/2833165/) — Lindhard response in CsV₃Sb₅
- [PRB 109, 035151](http://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.109.035151/fulltext) — Nesting in ThX₂ compounds
- [arXiv:1207.1106](https://arxiv.org/html/1207.1106v3) — Breakdown of Fermi liquid at QCP

---

## Recommendations for Next Steps

1. **Contact computational groups**: Email authors of arxiv:1606.06189 or arxiv:2307.12029 to ask if they've observed diagonal dominance in practice

2. **Check numerical properties**: For specific materials (FeSe, iron pnictides), compute χ₀ matrix and measure |diagonal|/|off-diagonal| ratio

3. **Pursue alternative proof strategies**: Focus on positivity/M-matrix structure or resolvent bounds instead of diagonal dominance

4. **Examine FeSe-specific literature**: See if Essenberger or similar SC-DFT papers make structural assumptions about χ

5. **Consider weaker claims**: Instead of "diagonally dominant everywhere", try:
   - "Block-diagonal dominant" (intraband vs interband)
   - "Asymptotically diagonal" in some limit
   - "Diagonal entries control eigenvalue bounds" (Gershgorin-like, but without strict dominance)

---

**Time spent**: ~4.5 hours  
**Conclusion**: Diagonal dominance of RPA susceptibility in orbital space is **not established** in the literature and appears **unlikely to hold** in general, due to essential role of local field effects (off-diagonal contributions). Main obstacle is that Coulomb interaction itself lacks orbital diagonal dominance.
