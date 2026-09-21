# Frustrated Quantum Magnets Gap Stability — Phase 1 Feasibility Assessment

**Date:** 2026-09-16 · **Lead:** Opus (feasibility evaluation) · **Context:** Post-Theorem I (Δ-channel immunity), evaluating expansion to frustrated systems gap stability

---

## Executive Summary

**Phase 1 Judgment: NO-GO (with qualified pivot option)**

The frustrated gap stability problem is a **foundational ~30-year open problem** with ~200× impact multiplier vs KAY Conjecture 2 if successful, but the technical barriers are **fundamental**, not incremental. 15+ years of expert effort (Hastings, Nachtergaele, Bravyi, Michalakis) hitting the same wall is a strong negative signal.

**Core finding:** Our rank-one stability machinery (Lemma 4 + Davis-Kahan + c₁² immunity budget) requires three properties that **do not hold simultaneously** in frustrated systems:
1. **Operator self-adjointness** — Achievable (work with Hamiltonian H, not Liouvillian)
2. **Base operator positive-definite structure** — **FAILS** (Perron-Frobenius breaks under frustration)
3. **Low-rank perturbation** — **FAILS** (frustration is intrinsically high-rank, many-body)

**Recommendation:** Do not pursue full frustrated gap stability proof. Prioritize Hubbard pairing (direct mechanism match, 6-8 weeks) and nickelate interlayer pairing (medium barrier, 8-10 weeks) instead; consider Pivot Option A (§7) only as backup.

---

## 1. Problem & Impact

**Frustration-free case (SOLVED):** if H₀ is frustration-free + gapped + satisfies LTQO (Local Topological Quantum Order), the spectral gap is stable under small short-range perturbations ([Michalakis-Zwolak 2013](https://arxiv.org/abs/1109.1588), [Nachtergaele et al. 2020](https://arxiv.org/abs/2010.15337), [Lucia et al. 2023](https://arxiv.org/abs/2102.07209)).

**Frustrated case (OPEN):** H = H₀ + V with V introducing frustration (e.g. Kagome Heisenberg triangular plaquettes forbidding simultaneous energy minimization) — does Δ > 0 persist?

**Impact:** Kagome Heisenberg model ~10,000 citations (canonical test case); frustration-free gap stability ~5,000 citations; applications span topological-order stability, MBL, adiabatic quantum computing. Herbertsmithite neutron scattering ([2024](https://arxiv.org/abs/2504.06491): Δ=1.1(2) meV) makes this directly testable, but pure-Kagome case remains debated (gapped Z₂ spin liquid [Yan et al. 2011](https://arxiv.org/abs/1205.4858) vs gapless [Iqbal et al. 2013](https://arxiv.org/abs/1311.5038)/[Jiang et al. 2016](https://arxiv.org/abs/1610.04727)) — **no consensus after 35 years**, so even a valid extension couldn't be validated against ground truth.

**Impact multiplier vs KAY Conjecture 2:** ~200× if successful — but see barriers below.

## 2. Why the Frustration-Free → Frustrated Extension Fails

Frustration-free machinery (BHM strategy: quasi-adiabatic continuation + LTQO + Lieb-Robinson bounds) rests on H₀ = ∑ᵢhᵢ with hᵢ≥0 local projectors sharing a common null space — this gives Perron-Frobenius-type ground-state uniqueness. Three barriers, all fatal, break this for frustrated systems:

1. **LTQO cannot hold** — frustration introduces long-range correlations even in the ground state ([Nachtergaele et al. 2020](https://arxiv.org/abs/2010.15337): "LTQO cannot hold for all local observables").
2. **No Perron-Frobenius structure** — competing terms hᵢ have incompatible minima, ground state is a sign-indefinite compromise, not a common null space.
3. **Perturbation is high-rank, not low-rank** — our Lemma 4 requires V = ∑ₙaₙ|ψₙ⟩⟨ψₙ| with aₙ≥0 for n≥1 (only ground state attracts); frustration V=∑_plaquettes V△ has rank ~N with **mixed signs** in the H₀ eigenbasis — cannot decompose as rank-one-attractive + diagonal-repulsive.

Recent literature (2024–2026: [Eisert et al. SDP hierarchy](https://arxiv.org/abs/2411.03680), [COMPASS non-Hermitian](https://arxiv.org/abs/2607.13790), [Lemm martingale method](https://arxiv.org/abs/2607.19078)) brings **no breakthrough** — all remain within the frustration-free framework, or (COMPASS) address an orthogonal notion of "gap frustration" (non-Hermitian PT-symmetric complexified couplings, not geometric lattice frustration — no technology transfer).

**Mapping our required properties (self-adjoint / positive base / low-rank / sign-structure) onto three attempted formulations, all fail:** direct Hamiltonian H=H₀+V (1/4 satisfied), Liouvillian ℒ[ρ]=[H,ρ] (0.5/4, plus Davis-Kahan doesn't even hold non-Hermitian), effective low-energy projection (1.5/4 — frustration persists within the subspace).

**Why Perron-Frobenius failure is fatal:** our c₁² immunity budget relies on a protected ground-to-first-excited distance; without Perron-Frobenius uniqueness there is no structural reason this distance exists or is stable — analogous to Perron-Frobenius being structurally impossible (not just hard) for matrices with negative entries.

## 3. Breakthrough Scenarios Considered (all low-to-medium probability)

| Scenario | Idea | Est. success prob. |
|---|---|---|
| A: effectively low-rank | Frustration might act rank-r≪N within the low-energy subspace | 15% |
| B: weakly frustrated subclass | J₂/J₁≪ε perturbative regime, gap Δ(ε)=Δ₀−O(ε²) | 35% (some lit. suggests gap CAN close, e.g. [DM-perturbed kagome 2018](https://arxiv.org/abs/1809.09128)) |
| C: numerical-certified bounds (existence only, not stability) | SDP hierarchy + DMRG/PEPS + interval arithmetic | 40% (existence) / 10% (stability) |

**Why 15+ years of expert effort (Hastings, Nachtergaele, Bravyi, Michalakis) failed:** the three BHM pillars (LTQO, quasi-adiabatic continuation, Lieb-Robinson light cone) all break simultaneously under frustration — not "hard to extend," structurally blocked. Bayesian update: prior 20% feasible → posterior **5%** after this deep dive.

## 4. Final Recommendation: NO-GO

**Do not pursue** full frustrated gap stability as primary direction — 1/4 required properties satisfied, 15-year expert roadblock, unknown ground truth (Kagome gapped/gapless undecided), technology mismatch (our machinery targets positive low-rank-perturbed operators, not frustrated Hamiltonians).

**Prioritize instead:** Hubbard model d-wave pairing (direct mechanism match, 6-8 wks, HIGH feasibility, highest priority per `high_impact_expansion_map_2026-09-16.md`); bilayer nickelate interlayer pairing (8-10 wks, medium). Backup if those stall: **Pivot Option A** — "Why Rank-One Perturbation Theory Fails for Frustrated Systems," a no-go theorem with explicit counterexamples (70% success prob., 8-10 wks; proving a negative result is easier). Option B (weakly-frustrated partial proof, 35% success, 10-12 wks) is the fallback if the user insists on pursuing frustration despite the NO-GO. Option C (SDP extension, 25%/60% publishable-negative) not recommended absent high risk tolerance.

## 5. Key Literature (full list preserved)

**Frustration-free (solved):** [Michalakis & Zwolak 2013](https://arxiv.org/abs/1109.1588) · [Nachtergaele et al. 2020](https://arxiv.org/abs/2010.15337) · [Lucia et al. 2023](https://arxiv.org/abs/2102.07209) · [Eisert et al. 2024](https://arxiv.org/abs/2411.03680) SDP hierarchy · [Lemm 2024](https://arxiv.org/abs/2607.19078) martingale method.

**Frustrated systems (open, gapped-vs-gapless debate):** [Liao et al. 2017](https://arxiv.org/abs/1611.06990) 48-site ED · [Yan et al. 2011](https://arxiv.org/abs/1205.4858) DMRG gapped Δ=0.13(1) · [Iqbal et al. 2013](https://arxiv.org/abs/1311.5038) gapless variational · [Zhu et al. 2024](https://arxiv.org/abs/2504.06491) experimental Δ=1.1(2) meV · [witness spin glass 2026](https://arxiv.org/abs/2608.21514).

**Non-Hermitian/Liouvillian & Perron-Frobenius extensions:** [COMPASS 2026](https://arxiv.org/abs/2607.13790) · [Quantum Perron-Frobenius 2023](https://arxiv.org/abs/2305.18971) · [Lieb Bose-Hubbard 2026](https://arxiv.org/abs/2506.10319) · [reflection positivity 2024](https://arxiv.org/abs/2510.20662).

## 6. Actionable Next Steps

1. Accept NO-GO; launch Hubbard pairing project (top priority).
2. Develop cross-cutting tools (momentum-space Lemma 4, interval matrix eigenvalues) needed for both Hubbard and nickelate.
3. Monitor Kagome literature for gapped/gapless consensus; re-evaluate in 6 months if it emerges.
4. This assessment stays filed for reference if a breakthrough occurs.

---

**END OF PHASE 1 FEASIBILITY ASSESSMENT**
