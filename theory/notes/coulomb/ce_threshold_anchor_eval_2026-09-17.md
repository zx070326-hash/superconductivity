# CE threshold evaluation at physical anchor points (2026-09-17)

## Task

Evaluate Theorem CE threshold μ̃*^CE at 7 physical material anchor points to determine if physical μ̃ = μ*/λ falls below the extremality threshold.

## Method

- **Theorem CE formula**: μ̃*^CE := 𝔱^d γ_P / (2 c₁^P ‖u_⊥^P‖ + 2 𝔱^d (c₁^P)²)
  - Six constants: (γ_P, c₁^P, ‖u_⊥^P‖, q₁^d, b^d, m^d) evaluated on A_P=𝔎(P,T)
  - Implementation: `CellExtremal` class from `p4d_coulomb_extremal_threshold.py`

- **Material anchor points**: 7 materials, 2 r values each (r=3, r=10), total 14 points
  - LaH10: 129/163/214/264 GPa (λ=3.62/2.67/2.06/1.73, μ*=0.10)
  - H₃S: 200 GPa (λ=1.84, μ*=0.16)
  - Nb: (λ=1.23, μ*=0.13)
  - MgB₂: (λ=0.62, μ*=0.13)

- **ϖ_rms proxy**:
  - LaH10 163 GPa: measured ϖ_rms = 0.783 (from ω_rms=1202.6 K, Tc^E=244.4 K)
  - Other 6 points: analytic floor ϖ_rms ≥ 0.8709/√λ
  - Floor is conservative: μ*>0 suppresses Tc → raises ϖ → bound remains valid

- **P proxy**: Einstein 2-atom equal-weight proxy
  - weights = [0.5, 0.5]
  - omegas = [ω_c/1.05, ω_c×1.05] where ω_c = 2π ϖ_rms T
  - spread = 1.05 (near case, most conservative)

- **Computation**: N=400, float64, OMP_NUM_THREADS=1

## Results

### RETRACTION (2026-09-17): First table had N_c bug

The first evaluation (ce_threshold_seven_materials.py) reported r=3 and r=10 as having **identical** thresholds despite different N_c values, which is impossible (the threshold decreases monotonically in N_c). **The root cause was never located.** The explanation previously written here ("module-scope N / NN=N shadowing") is retracted: it is not what `CellExtremal.__init__` does, and it would affect the truncation N, not N_c. The script as it now stands on disk passes N_c correctly and, when rerun by team-lead, reproduces the corrected table below to 6 digits — so the defect was either in a since-overwritten version of the script or in the hand-transcription of its output into the first table. Which of the two is unknown.

**Corrected table below matches team-lead's independent recomputation** (same class, explicit parameter passing verified) and a fresh rerun of the script.

### Summary table (N=400, CORRECTED)

| Material        | lam  | mu*  |  r | varpi |  Nc | mutilde | mu*_CE   | margin | verdict | note     |
|-----------------|------|------|----|-------|-----|---------|----------|--------|---------|----------|
| LaH10 129 GPa   | 3.62 | 0.10 |  3 | 0.4577|   1 | 0.0276  | 0.0955   | 3.46   | PASS    | floor    |
| LaH10 129 GPa   | 3.62 | 0.10 | 10 | 0.4577|   4 | 0.0276  | 0.0572   | 2.07   | PASS    | floor    |
| LaH10 163 GPa   | 2.67 | 0.10 |  3 | 0.7830|   2 | 0.0375  | 0.1675   | 4.47   | PASS    | measured |
| LaH10 163 GPa   | 2.67 | 0.10 | 10 | 0.7830|   7 | 0.0375  | 0.1111   | 2.96   | PASS    | measured |
| LaH10 214 GPa   | 2.06 | 0.10 |  3 | 0.6068|   1 | 0.0485  | 0.1504   | 3.10   | PASS    | floor    |
| LaH10 214 GPa   | 2.06 | 0.10 | 10 | 0.6068|   6 | 0.0485  | 0.0800   | 1.65   | PASS    | floor    |
| LaH10 264 GPa   | 1.73 | 0.10 |  3 | 0.6620|   1 | 0.0578  | 0.1715   | 2.97   | PASS    | floor    |
| LaH10 264 GPa   | 1.73 | 0.10 | 10 | 0.6620|   6 | 0.0578  | 0.0912   | 1.58   | PASS    | floor    |
| H3S 200 GPa     | 1.84 | 0.16 |  3 | 0.6421|   1 | 0.0870  | 0.1638   | 1.88   | PASS    | floor    |
| H3S 200 GPa     | 1.84 | 0.16 | 10 | 0.6421|   6 | 0.0870  | 0.08711  | 1.001  | PASS    | floor    |
| Nb              | 1.23 | 0.13 |  3 | 0.7855|   2 | 0.1057  | 0.1682   | 1.59   | PASS    | floor    |
| Nb              | 1.23 | 0.13 | 10 | 0.7855|   7 | 0.1057  | 0.1115   | 1.05   | PASS    | floor    |
| MgB2            | 0.62 | 0.13 |  3 | 1.1060|   3 | 0.2097  | 0.2224   | 1.06   | PASS    | floor    |
| MgB2            | 0.62 | 0.13 | 10 | 1.1060|  11 | 0.2097  | 0.1544   | 0.736  | **FAIL**| floor    |

**Pass rate: 13/14 (92.9%)**

margin := μ̃*_CE / μ̃  (ratio of threshold to physical value)
- PASS: μ̃ < μ̃*_CE (physical value below threshold)
- FAIL: μ̃ ≥ μ̃*_CE (physical value at or above threshold)

### LaH10 verdict

**All four LaH10 pressure points PASS the CE threshold at both r=3 and r=10.**

- Margins: **1.577× to 4.473×** (r=3 gives 2.97–4.47×, r=10 gives 1.58–2.97×). The minimum is LaH10 264 GPa r=10 (1.577), followed by 214 GPa r=10 (1.648) — not 2.07 as previously stated here.
- LaH10 163 GPa uses measured ϖ_rms=0.783; others use floor proxy
- Physical μ̃ ranges from 0.0276 (129 GPa) to 0.0578 (264 GPa)
- Thresholds μ̃*_CE: r=3 gives 0.096–0.172, r=10 gives 0.057–0.111

LaH10 is certified by Theorem CE at both tested r values.

### Critical observation: H₃S 200 GPa r=10 near-threshold

**H₃S 200 GPa at r=10 has margin 1.001** (μ̃=0.08696, μ̃*_CE=0.08711), essentially at the threshold boundary. This is the tightest margin among all PASS points and suggests the theorem's certification boundary is physically meaningful, not artificially loose.

### MgB₂ r=10 failure investigation

MgB₂ at r=10 is the **ONLY failure** among 14 points. Investigation:

#### N convergence (MgB₂ r=10, ϖ_rms=1.106, Nc=11)

| N   | μ*_CE    | margin | delta      | verdict |
|-----|----------|--------|------------|---------|
| 200 | 0.154357 | 0.7362 | ---        | FAIL    |
| 400 | 0.154356 | 0.7362 | -0.000001  | FAIL    |
| 800 | 0.154356 | 0.7362 | -0.000000  | FAIL    |

**Threshold is converged at N=400.** No indication of rising threshold with N.

#### ϖ_rms sensitivity (N=400, Nc = ⌊r × factor × ϖ_floor⌋)

| factor | varpi  | Nc | mu*_CE   | margin | verdict |
|--------|--------|--- |----------|--------|---------|
| 1.00   | 1.1060 | 11 | 0.154356 | 0.7362 | FAIL    |
| 1.10   | 1.2166 | 12 | 0.168423 | 0.8032 | FAIL    |
| 1.20   | 1.3273 | 13 | 0.181623 | 0.8662 | FAIL    |
| 1.30   | 1.4379 | 14 | 0.193971 | 0.9251 | FAIL    |
| 1.50   | 1.6591 | 16 | 0.215780 | 1.0291 | **PASS**|

**The failure persists until ϖ_rms is 50% above the floor.** If true MgB₂ ϖ_rms is significantly higher than the floor proxy (ϖ_floor = 0.8709/√0.62 = 1.106), the material could pass.

#### Interpretation: one-sided bound direction

The floor proxy ϖ_rms ≥ 0.8709/√λ is a **lower bound** on true ϖ_rms:
- Using floor → conservative Nc
- Lower Nc → affects six constants → affects μ̃*_CE

Direction audit needed: Does lower ϖ_rms → lower or higher μ̃*_CE?

From sensitivity table: **higher ϖ → higher Nc → higher μ̃*_CE**

Therefore:
- Floor ϖ gives **lower bound** on μ̃*_CE
- FAIL from a lower-bound threshold is **NOT conservative**
- This FAIL cannot be dismissed as one-sided artifact

**Conclusion**: MgB₂ r=10 is a **genuine failure** of Theorem CE certification. The theorem does not certify this point. The physical μ̃=0.2097 exceeds the computed threshold μ̃*_CE=0.154 even at converged N and even accounting for floor-proxy conservatism.

### Comparison with team-lead's three-constant version

Team-lead reported: Q(P) ≤ Q(E) holds 42/42 at physical anchors, but CE threshold only certifies 27/42.

My six-constant evaluation: 13/14 PASS (92.9%), with MgB₂ r=10 as sole failure.

The 13/14 pass rate is consistent with team-lead's statement that CE threshold is incomplete (not all Q(P)≤Q(E) points are certified by the threshold formula).

**Cross-check**: Within task scope (7 materials × 2 r values), my implementation identifies the same gap between exact operator result (should pass) and threshold certification (doesn't always certify).

## What couldn't be delivered

1. **Three-constant vs six-constant 10% tolerance check**: Team-lead's three-constant version outputs were not provided in sufficient detail to perform point-by-point comparison. I evaluated the six-constant formula as implemented in the repository.

2. **Measured ϖ_rms for 6/7 materials**: Only LaH10 163 GPa has measured ω_rms. Floor proxy 0.8709/√λ used for others is rigorous (analytic) but may be loose, especially in low-λ regime (MgB₂).

3. **Root cause of MgB₂ failure**: Identified as genuine (not N-truncation, not floor artifact based on sensitivity direction), but did not diagnose which of the six constants or which structural assumption (H1)-(H5) is responsible.

## Rigor classification

- **Threshold formula (Theorem CE statement)**: analytic (Lemmas 0–4 in `P4D_coulomb_extremal_threshold.md` v0.2)
- **All numerical evaluations (μ̃*_CE values)**: **numerical-evidence**
  - float64 + finite N=400
  - Six constants not interval-certified
  - Gate G1-G4 use float64 at N=300 baseline
  - N convergence confirms 400 is adequate for 6-digit threshold stability
- **ϖ_rms floor proxy**: **analytic** (inequality 0.8709/√λ is proven)
- **Floor proxy direction**: **numerical-evidence** (higher ϖ → higher Nc → higher threshold verified numerically on test cases)
- **LaH10 PASS verdicts**: **numerical-evidence** (thresholds computed at numerical-evidence level)
- **MgB₂ r=10 FAIL verdict**: **numerical-evidence** (threshold value is numerical-evidence)
- **MgB₂ FAIL genuineness**: **numerical-evidence** (N-convergence and direction audit are numerical)

**Correction (2026-09-17)**: Initial rigor note incorrectly stated Theorem CE thresholds as certificate-level. The theorem **statement** (Lemmas 0-4) is analytic, but any concrete threshold **value** is numerical-evidence because the six constants are evaluated via float64 eigendecomposition at finite N without interval certification.

## Files

- Anchor evaluation: `theory/notes/_scratch_ir/ce_threshold_seven_materials.py`
- MgB₂ investigation: `theory/notes/_scratch_ir/mgb2_fail_investigation.py`
- Implementation: `theory/derivations/p4d_coulomb_extremal_threshold.py` (CellExtremal class)
- Theorem statement: `theory/derivations/P4D_coulomb_extremal_threshold.md` (v0.2)
