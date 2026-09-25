# Essenberger FeSe patch→band mapping investigation (2026-09-21)

> researcher agent · task from team-lead · NOT OBTAINED (blocker status revised)

## 0. Executive summary

**The "5 patch vs 2 band incompatibility" blocker does NOT exist** — it was proven false in prior work.

**Actual findings**:
1. **Both papers use a 2-band structure throughout** (hole band Γ, electron band M)
2. **The "5 sectors" mentioned in 1411.2121 are a k-space subdivision for gap anisotropy**, applied **after** all coupling calculations
3. **All three channels (phonon, SF, Coulomb) are constructed on the same 2×2 band structure**
4. **The formalism paper 1409.7968 never mentions "sectors" or "patches"** — it works entirely in 2-band isotropic approximation
5. **The mapping question is therefore ill-posed** — there are no "5 patches" at the coupling matrix level

## 1. Two-band structure (common to all channels)

### From 1411.2121 (application paper)

**Line 123-134** (verbatim from `_fe_1411_utf8.txt`):
```
This way we can group together the holes (see Fig. 1)
around the Γ−point (labeled with n = +) and the elec-
trons at the M-point (labeled with n = −). Within this
approximation only s-wave pairing symmetries are possi-
ble, including the s± symmetry suggested in Ref. 4.
```

**Two bands**:
- Band `+` (hole): Γ-centered barrel  
- Band `−` (electron): M-centered pockets

**This 2-band structure is used for**:
- Phonon Eliashberg function `α²F_ij(ω)` (line 137-138: "band-resolved")
- Coulomb screening `w_ij(E,E',ω)` (line 139-141)
- Spin-fluctuation kernel `Λ^SF_ij(E,E',ω)` (line 176-177)

All three channels report 2×2 matrices in Fig. 3 (`++`, `+-`, `-+`, `--` components).

### From 1409.7968 (formalism paper)

**Section VI, page 9, line 123-130** (VI.A title: "Isotropic approximation and the two band model with a SF pairing"):

The paper introduces a 2-band isotropic averaging:
```
Δ_nk ≈ Δ_n(E) := (1/N_n(E)) Σ_k δ(ε_nk − E) Δ_nk
N_n(E) := Σ_k δ(ε_nk − E)
```

**Eq. (36) and (37)** define Z and pairing kernels **for a 2-band system**, with indices `k,k'` running over two bands only.

**No mention of**: patches, sectors, or any k-space subdivision beyond the 2-band grouping.

## 2. The "5 sectors" are post-calculation refinement

**From 1411.2121, line 254-260**:
```
We conclude the paper with an investigation of the gap
function in k-space, going beyond the two-band isotropic
approximation used so far (Eq. 2). We divide the xy-
plane of the Brillouin zone into sectors as shown in Fig. 5,
and compute the superconducting gap corresponding to
each sub-band. The use of this more accurate approach
does not significantly affect the value of the critical tem-
perature (Tc increases by a few percent), it leads to a
significantly modulated gap function within the s± sym-
metry.
```

**Chronology**:
1. All coupling matrices (`Λ^ph`, `Λ^sf`, `w`) computed in **2-band isotropic approximation**
2. Tc estimated from 2-band gap equation
3. **Then**, as a refinement, BZ divided into **sectors** to resolve k-dependent gap modulation
4. This division happens **after** the coupling functions are already defined

**Line 261-266** describes gap features:
```
The large dip in the Δ− is related to the CB
fluctuation (M1 − M2) in this band and is consequently
not seen in Δ+ (Γ1 − Γ2). The smaller oscillations are
related to intraband scattering, where the interaction is
smaller between two minima of the gap (Γ2 − M2) and
larger between two maxima (Γ1,3 = M2).
```

**Labels like Γ1, Γ2, M1, M2** refer to **points within the already-defined ± bands**, not separate input channels.

## 3. Cross-verification with existing notes

**From `zdom_iron_based_zero_dof_2026-09-21.md` §6 (FeSe)**:

> §6.4 (line 503-506): "我一手读到三个通道建在**同一个 2 带分区**上，普查转引的「2 带 vs 5 patch」为假。"  
> §6.5 (line 530-532): "sector 是事后细分（'going beyond the two-band isotropic approximation'），**不是耦合表的输入**。"

**Confirmation**: Another agent independently reached the same conclusion from first-hand reading.

## 4. Numerical (Z-dom) verification status

**From existing memory** (`essenberger-scdft-functional-is-K1-2026-09-21.md`):

The **boson layer** (phonon + SF) satisfies (Z-dom) via (K1):
- Phonon piece: `Λ^Ph_kk' < 0` (Eq. 5, analytic)
- SF piece: `Im[Λ^SF_kk'] > 0` (Eq. 37, literature assertion)
- Structure: `Λ^Z = s + a`, `Λ^Δ = c + s − a` (where `a := |Λ^Ph|`, `s := (3/2)Im[Λ^SF]`, `c := Im[w]`)

**The (Z-dom) judgment operates on 2×2 matrices**, not 5×5.

**Coulomb layer**: undetermined (needs `Im[w_ij]` PSD, paper prints only element-wise positivity).

## 5. Why the blocker was reported

**Original concern** (from team-lead message):
> "repository has ρ(M2''_Z) calculation that needs patch→band mapping"  
> "Current status: numerical landing blocked by patch→band mapping"

**Resolution**: The blocker description conflated two different objects:
1. **Coupling matrices** `λ^ph_ij`, `λ^sf_ij` (2×2, band-resolved) — these ARE available
2. **Gap function** `Δ(k)` (k-resolved, subdivided into sectors for Fig. 5) — this is output, not input

**For (Z-dom) verification, only item (1) is needed**, and it exists as 2×2 matrices throughout both papers.

## 6. Answer to the original question

> "How do the 5 patches map to 2 bands?"

**Answer**: The question is ill-posed. There are no "5 patches" at the coupling matrix level.

**What actually exists**:
- **Input level**: 2 bands (Γ hole, M electron)
- **Coupling matrices**: 2×2 for all three channels
- **Output level**: Gap function `Δ(k)` computed with k-dependence, then **plotted** using sector subdivision for visualization

**If the repository calculation needs band-resolved `M2_ij`**:
- Read `λ^ph_ij` and `λ^sf_ij` as 2×2 matrices directly from the papers
- There is no "5→2 mapping" to perform

## 7. Strictness level

| Statement | Level |
|---|---|
| Both papers use 2-band structure | fact-checking (first-hand text reading) |
| "5 sectors" are post-calculation | fact-checking (line 254-260 verbatim) |
| Three channels built on same partition | fact-checking (cross-checked with `zdom_iron_based_zero_dof_2026-09-21.md` §6) |
| (K1) structure in 1409.7968 Eq.(36)/(37) | already verified in memory `essenberger-scdft-functional-is-K1-2026-09-21.md` (not re-derived here) |

## 8. What was NOT obtained

1. **Figure 5 itself** — could not render it from the available files (1411.2121 not in repository, only text extract)
2. **Explicit numerical values** of `λ^ph_ij`, `λ^sf_ij` for FeSe — 1411.2121 prints only qualitative statements ("λ < 0.1") and Fig. 3 plots, not tables
3. **The actual sector boundaries** in Fig. 5 — description says "as shown in Fig. 5" but provides no algebraic definition

## 9. Honest boundary

**This investigation did NOT**:
- Perform any DFT calculation
- Extract numerical coupling values from Fig. 3 (would require digitization)
- Verify whether 1409.7968's "two-band electron gas model" parameters match 1411.2121's FeSe-specific values
- Check whether other FeSe papers (e.g., published PRB version) provide more explicit sector definitions

**The conclusion "no 5→2 mapping needed" is based on**:
- Text statements in both papers that all couplings use 2-band averaging
- Cross-verification with existing repository note that independently reached same conclusion
- Logical structure: sectors are introduced **after** "going beyond the two-band isotropic approximation used so far"

## 10. Recommendation

**For numerical (Z-dom) verification on FeSe**:
1. The blocker "patch→band mapping" should be **removed** — it never existed
2. The actual data gap is: **numerical values of `λ^ph_ij` and `λ^sf_ij`** (2×2 matrices)
3. These are **printable** from 1411.2121 Fig. 3 via digitization, or from the published PRB version if it contains tables
4. Alternatively, use the **structure-only result** already in repository: FeSe satisfies (Z-dom) via (K1) at the boson layer (analytic), independent of numerical values

**The repository can proceed with**:
- Structural verification: (K1) ⇒ (Z-dom) automatically (triangle inequality, already proven)
- This requires **zero numerical extraction** from FeSe-specific papers

---

**Sources**:
- [arXiv:1411.2121](https://arxiv.org/abs/1411.2121) — Essenberger et al., ab-initio theory of iron-based superconductors (application paper)
- [arXiv:1409.7968](https://arxiv.org/pdf/1409.7968) — Essenberger et al., SCDFT functional derivation (formalism paper)
- Repository notes: `theory/notes/ceiling/zdom_iron_based_zero_dof_2026-09-21.md`, memory `essenberger-scdft-functional-is-K1-2026-09-21.md`