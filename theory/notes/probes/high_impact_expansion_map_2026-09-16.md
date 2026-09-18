# High-Impact Expansion Map for Rank-One Stability Technology

**Date:** 2026-09-16  
**Context:** Post-Theorem I (Δ-channel immunity), identifying problems with impact > KAY Conjecture 2

## Executive Summary

Our rank-one stability machinery (Lemma 4: rank-one leakage + Davis-Kahan perturbation + interval arithmetic certificates) is **not superconductivity-specific**. Surveying open problems whose bottleneck is "spectral gap / rank-one dominance stability under perturbation" gives three candidates, ranked by impact × feasibility:

1. **Hubbard Model d-Wave Pairing Instability** — direct mechanism match, highest impact
2. **Bilayer Nickelate Interlayer Pairing Stability** — emerging experimental frontier, mechanism adaptable
3. **Frustration-Free → Frustrated Gap Stability** — foundational theory, high barrier (**superseded by NO-GO verdict**, see `frustrated_gap_feasibility_2026-09-16.md`)

Priority chosen: attack Hubbard first (#1), run nickelate as a parallel track (#2), keep frustrated gap as long-term/exploratory only (#3, later judged NO-GO as primary target).

---

## Candidate 1: Hubbard Model d-Wave Pairing Instability

**Problem:** Rigorously prove the 2D single-band Hubbard model at finite doping exhibits d-wave pairing dominance over CDW/SDW/s-wave.

**Why it matters more than KAY Conj. 2:** central 40+-year-open problem for cuprate theory (~15,000 citations across reviews); recent neural-network result ([Reh et al. 2026](https://arxiv.org/abs/2608.12465)) resolves 1/8-doping numerically but leaves rigorous bounds open; directly testable via ARPES/STM/Raman; techniques transfer to twisted bilayer graphene, nickelates, organics.

**Operator structure:** pairing susceptibility χ_α(q,T) = ∫dτ ⟨Δ_α†(q,τ)Δ_α(q,0)⟩ with Δ_α = ∑_k φ_α(k) c_k↑c_{-k+q}↓. Superconductivity ⟺ Γ_d-wave (Bethe-Salpeter eigenvalue) is largest (rank-one dominance) and > 0.

**Lemma 4 mapping:** the Bethe-Salpeter kernel K[α,α'] is positive in Matsubara space; Hubbard U adds diagonal-dominant negative (repulsive) terms in all channels; spin-fluctuation exchange adds a rank-one attractive perturbation in the d-wave channel. Adaptation: (1) rank-one leakage — if a_n≥0 for n≥1 in the attractive perturbation's spectral decomposition, only n=0 (d-wave) leaks; (2) Davis–Kahan bound tracks λ_d-wave shift under CDW/SDW perturbation via immunity budget c₁²; (3) certificate via mpmath.iv on finite-size QMC data.

**Mechanism match: High** — self-adjoint ✓, base positive ✓, low-rank perturbation ✓ (spin-fluctuation exchange dominates), sign structure ✓ (d-wave attractive a₀<0, others repulsive a_n≥0). Adaptation needed: real-space → momentum-space kernel, finite-T Matsubara formalism (structure unchanged), bridging QMC (sign-problem-limited) to certificate via normalizing-flow QMC ([Nicoli et al. 2026](https://arxiv.org/abs/2603.18205)).

**Feasibility: High**, ~6-8 weeks to proof-of-concept (momentum-space Lemma 4 adaptation 2wk → QMC data interface 2wk → certificate 2wk → finite-size scaling + writeup 2wk). Tractable because normalizing-flow QMC cuts the sign problem "by an order of magnitude," and neural-net variational upper bounds ([Reh et al.](https://arxiv.org/abs/2608.12465)) can sandwich our lower bounds. **Bottleneck risk:** if c₁² is too small at accessible sizes (N≤20×20), result downgrades from certificate to numerical-evidence.

**Key refs:** normalizing-flow sign-problem fix ([2603.18205](https://arxiv.org/abs/2603.18205)); sign-resolved QMC bias origin ([2512.04056](https://arxiv.org/abs/2512.04056)); rigorous pairing only for symmetry-special Hubbard variant, not standard model ([1904.07138](https://arxiv.org/abs/1904.07138)); d-wave enhancement without dominance ([1304.1273](https://arxiv.org/abs/1304.1273)); NN quantum states close 1/8-doping numerically ([2608.12465](https://arxiv.org/abs/2608.12465)).

---

## Candidate 2: Bilayer Nickelate Interlayer Pairing Stability

**Problem:** Prove interlayer singlet pairing in bilayer nickelates (La₃Ni₂O₇ under pressure) is robust against competing density-wave orders and orbital-selective Mott transitions.

**Why it matters:** hot frontier since 2024 discovery of ~80K Tc, two-orbital physics (d_{x²-y²}+d_{z²}) richer than cuprates, qualitatively new interlayer mechanism; 100+ papers 2024-2026; ARPES confirms gap anisotropy but symmetry (s/d/s±) still debated; pressure tuning gives systematic access to coupling strength.

**Operator structure** ([2606.15298](https://arxiv.org/abs/2606.15298)): bilayer two-orbital Hubbard model, Δ_inter = ∑_{k,layer} φ(k) c_{k,l1,d_{z²}}↑c_{k,l2,d_{z²}}↓ — d_{z²} orbitals form local interlayer singlets hybridizing with itinerant d_{x²-y²} electrons. Competing orders: CDW (interlayer charge imbalance), SDW (in-layer AFM), orbital-selective Mott (d_{z²} localizes, d_{x²-y²} stays metallic).

**Lemma 4 mapping:** two-orbital susceptibility matrix χ[channel×channel] is positive (Matsubara sum); interlayer hybridization gives rank-one dominant structure in the pairing channel; Hund's J is a diagonal perturbation disfavoring pairing; non-local Coulomb V is off-diagonal, couples CDW and pairing. Adaptation: extend Lemma 4 to 2×2 block structure (bonding/antibonding pairing modes), compute c₁² for the interlayer singlet as the CDW/SDW-leakage protection, certify via DMFT/cluster-DMFT susceptibilities + interval arithmetic.

**Mechanism match: Medium** — self-adjoint ✓, base positive ✓, perturbation low-rank but not purely rank-one (Hund+Coulomb, 2-3 modes) ✓, sign structure (CDW/SDW satisfying a_n≥0) ⚠ unverified. Needs: block-diagonal generalization of Lemma 4, DMFT data interface (cleaner than QMC but finite-cluster artifacts harder to quantify), tracking c₁²(pressure) near phase boundaries.

**Feasibility: Medium-High**, ~8-10 weeks. Tractable because published DMFT susceptibility data exists ([npj QM 2026](https://www.nature.com/articles/s41535-026-00883-7)), finite-cluster DMFT is numerically exact (no sign problem), and ARPES gives an immediate experimental check. **Bottleneck risk:** if Hund's J is too strong, orbital-selective Mott transition may preclude low-rank structure entirely (would need "rank-two dominance with gap" machinery instead); gap symmetry debate may still be open in 6 months, weakening immediate impact.

**Key refs:** interlayer mechanism review ([2606.15298](https://arxiv.org/abs/2606.15298)); nonlocal Coulomb / DMFT susceptibilities ([npj QM 2026](https://www.nature.com/articles/s41535-026-00883-7)); ARPES leading-edge gap observation ([Nat. Phys. 2026](https://www.nature.com/articles/s41567-026-03396-z)); BEC-BCS crossover binding energies ([Nat. Comm. Phys. 2024](https://www.nature.com/articles/s42005-024-01854-9)); Hund-scenario insufficiency argument ([2026](https://www.nature.com/articles/s41535-026-00871-x)).

---

## Candidate 3: Frustration-Free → Frustrated Gap Stability

**Problem:** extend rigorous spectral gap stability from frustration-free Hamiltonians to frustrated spin systems (triangular/kagome Heisenberg AFM). Foundational (~5,000 + ~10,000 citations combined), would settle the 30-year kagome QSL-vs-gapped debate, testable via herbertsmithite neutron scattering — **but** mechanism match is weak: Liouvillian route is non-self-adjoint, and frustration perturbations are not obviously low-rank or sign-definite (a_n≥0 unverified). Estimated feasibility Low-Medium, 12-16 week high-risk timeline. Full evaluation and detailed no-go argument superseded by the dedicated Phase-1 assessment: see `frustrated_gap_feasibility_2026-09-16.md` (verdict: **NO-GO** as primary target, three of four required mechanism properties fail — Perron-Frobenius structure and low-rank perturbation both break under frustration).

Key refs: frustration-free gap stability via LTQO ([2102.07209](https://arxiv.org/abs/2102.07209)); SDP certificate hierarchy, frustration-free only ([2411.03680](https://arxiv.org/abs/2411.03680)); kagome 48-site ED, gapless suggestion with finite-size caveat ([1611.06990](https://arxiv.org/abs/1611.06990)); variational evidence for gapless kagome ([1311.5038](https://arxiv.org/abs/1311.5038)); DMRG evidence for gapped Z₂ spin liquid ([1205.4858](https://arxiv.org/abs/1205.4858)).

---

## Mechanism Match Summary

| Problem | Self-Adjoint? | Base Positive? | Low-Rank Pert.? | Sign Structure Clear? | Feasibility |
|---|---|---|---|---|---|
| Hubbard Pairing | ✓ | ✓ | ✓ (spin fluctuations) | ✓ | **High** |
| Nickelate Interlayer | ✓ | ✓ | ~✓ (2-3 modes) | ⚠ needs verification | **Medium-High** |
| Frustrated Gap | ✗ (Liouvillian) | ✓ (null space) | ✗ (many-body) | ⚠ unknown | **Low-Medium** (→ NO-GO, see dedicated note) |

**Impact multiplier vs KAY Conj. 2** (citations + experimental testability basis): Hubbard ~100×, nickelate ~50×, frustrated gap ~200× if successful (high risk; not pursued as primary — see Candidate 3).

---

## Priority Ranking & Next Steps

**#1 Hubbard d-wave pairing (attack first).** Direct mechanism match, 6-8wk to proof-of-concept / 3-4mo to certificate+paper. Steps: adapt Lemma 4 to momentum-space Bethe-Salpeter (2wk) → QMC data collaboration/rerun (1wk) → c₁² certificate at 1/8-doping (3wk) → finite-size extrapolation (2wk) → draft "Rigorous Lower Bound on d-Wave Pairing Dominance in the 2D Hubbard Model." Fallback: downgrade to certified numerical-evidence if c₁² too small — still PRB/PRL-publishable.

**#2 Bilayer nickelate interlayer pairing (parallel track).** 8-10wk. Steps: 2×2 block generalization of Lemma 4 (3wk, parallelizable with #1) → extract/rerun DMFT susceptibilities (2wk) → c₁² for interlayer singlet (2wk) → pressure-doping phase diagram (3wk) → draft "Certified Interlayer Pairing Stability in Bilayer Nickelate Superconductors." Fallback: if Hund's coupling breaks low-rank structure, pivot to orbital-resolved (separate d_{x²-y²}/d_{z²}) certificates.

**#3 Frustrated gap stability (long-term/exploratory only, per NO-GO verdict).** Not a primary track. If pursued as backup: toy 1D frustrated chain first (4wk) to validate a non-Hermitian Lemma 4 generalization before any kagome attempt; a negative result ("why rank-one leakage cannot prove frustrated gap stability") is itself publishable.

---

## Cross-Cutting Technical Developments (needed for all three)

1. **Momentum-space generalization of Lemma 4** — real-space K(r,r') → K(k,k'); leakage condition a_n≥0 applies to momentum eigenmodes. ~1-2wk (mostly notational).
2. **Finite-temperature extension** — T=0 ground state → Matsubara sum at finite T; c₁²(T) becomes T-dependent. ~1wk.
3. **Interval arithmetic for matrix functions** — scalar h,λ₂ → matrix susceptibility χ[α,α'], needs Gershgorin-type eigenvalue certificates; extend `mpmath.iv` to matrix ops. ~2wk.
4. **QMC/DMFT data interface** — parser turning QMC/DMFT output (with Monte Carlo error bars as interval widths) into certificate input. ~1wk.

---

**END OF REPORT**
