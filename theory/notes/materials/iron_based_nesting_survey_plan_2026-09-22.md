# Iron-Based Superconductor Nesting Survey Plan

**Date**: 2026-09-22  
**Team**: Team C Phase 2  
**Status**: In progress  
**Target**: Geometric nesting universality across 5 iron-based materials

---

## Objective

Establish the universality of Q_AFM nesting in iron-based superconductors and investigate potential correlation with Tc.

### Primary goals:
1. Compute Fermi surface topology for 5 representative materials
2. Calculate nesting function ξ(q) across Brillouin zone
3. Identify Q_AFM for each material
4. Test correlation: Tc vs ξ(Q_AFM)

### Deliverables:
- Fermi surface topology comparison (5-panel figure)
- ξ(q) heatmaps for each material
- Q_AFM summary table
- Tc correlation analysis (scatter plot + fit)

---

## Material Selection

| Material | Family | Tc (K) | Fermi surface | Status |
|----------|--------|--------|---------------|--------|
| LaFeAsO | 1111 | 26 | Γ hole + M electron | Parameter search |
| BaFe₂As₂ | 122 | ~38 (doped) | Same topology | Parameter search |
| LiFeAs | 111 | 18 | Same topology | Parameter search |
| FeSe | 11 | 8 | Same topology | ✓ Have (Team B) |
| FeTe₀.₅Se₀.₅ | 11 (doped) | 14 | Modified | Parameter search |

**Rationale**: These 5 materials span different families and Tc ranges while sharing the basic Fe-pnictide/chalcogenide structure.

---

## Technical Approach

### Step 1: Tight-binding model construction (Week 1)

**For each material, need**:
- 5-orbital tight-binding parameters (or simplified 2-3 orbital)
- Hopping integrals: t_ij between d-orbitals (d_xz, d_yz, d_xy, d_x²-y², d_3z²-r²)
- Crystal structure parameters (lattice constant, Fe-As/Se bond angles)

**Parameter sources** (priority order):
1. **Kuroki et al. PRL 101, 087004 (2008)** - LaFeAsO canonical
2. **Graser et al. arXiv:1003.0133 (2010)** - BaFe₂As₂ 3D tight-binding
3. **Borisenko et al. PRL 105, 067002 (2010)** - LiFeAs ARPES + fit
4. **Team B existing** - FeSe 5-orbital model
5. **Zhang et al. PRB 80, 014506 (2009)** - FeTe₀.₅Se₀.₅

**Fallback**: If 5-orbital parameters unavailable, use 2-orbital simplified model (d_xz, d_yz only) which captures essential Fermi surface topology.

### Step 2: Fermi surface calculation (Week 1)

**Method**: 
- Diagonalize H(k) on dense k-grid (200×200 for 2D BZ)
- Find Fermi level μ (chemical potential at given filling)
- Extract Fermi surface contours: {k : ε_α(k) = μ for any band α}

**Output**:
- Fermi surface plot for each material
- Identify Γ-centered hole pockets and M-centered electron pockets
- Measure pocket sizes (relevant for nesting quality)

### Step 3: Nesting function computation (Week 1-2)

**Definition**:
```
ξ(q) = Σ_{α,β} ∫ dk δ(ε_α(k)) δ(ε_β(k+q))
```

**Numerical implementation**:
```python
def compute_nesting(q, epsilon_k, mu, window=0.05):
    """
    Count k-points where both k and k+q are on Fermi surface
    
    Parameters:
    - q: wavevector (qx, qy)
    - epsilon_k: band energies on k-grid
    - mu: chemical potential
    - window: Fermi surface window (eV)
    
    Returns:
    - ξ(q): normalized nesting degree
    """
    on_fermi_k = (abs(epsilon_k - mu) < window)
    on_fermi_kq = (abs(epsilon_kq - mu) < window)
    return sum(on_fermi_k * on_fermi_kq) / N_k
```

**Output**:
- ξ(q) heatmap over full BZ (100×100 grid sufficient)
- Identify global maximum → Q_AFM candidate

### Step 4: Q_AFM identification (Week 2)

**Expected Q_AFM positions** (based on Fermi surface geometry):
- **Stripe order**: Q_AFM ≈ (π, 0) or (0, π) - connects Γ hole to M electron pockets
- **Checkerboard order**: Q_AFM ≈ (π, π) - less common in Fe-pnictides

**Verification**:
1. Check ξ(Q_AFM) >> ξ(other q) (enhancement factor > 2×)
2. Compare with experimental spin fluctuation peak positions (neutron scattering)
3. Consistency check: Q_AFM should connect dominant Fermi surface features

**Hypothesis**: All 5 materials have Q_AFM near (π,0) due to universal Γ-M topology.

### Step 5: Tc correlation analysis (Week 2)

**Candidate correlation observables**:
1. **Raw nesting**: ξ(Q_AFM)
2. **Normalized nesting**: ξ(Q_AFM) / ⟨ξ(q)⟩_BZ
3. **Pocket matching**: Size ratio of Γ hole / M electron pockets
4. **Fermi surface area**: Total Fermi surface length

**Statistical tests**:
- Pearson correlation coefficient (linear)
- Spearman rank correlation (monotonic)
- Fit models: Tc = a × ξ^b (power law) or Tc = a + b×ξ (linear)

**Success criteria**:
- **Strong correlation**: |r| > 0.8 → potential PRL
- **Moderate correlation**: 0.5 < |r| < 0.8 → interesting but needs more materials
- **Weak correlation**: |r| < 0.5 → nesting alone doesn't determine Tc (other factors dominate)

---

## Computational Details

**Tools**:
- Python + NumPy (tight-binding diagonalization)
- Matplotlib (Fermi surface + heatmap visualization)
- SciPy (correlation analysis)

**Reuse from Phase 1**:
- `build_5orbital_hamiltonian()` from Team B
- `compute_nesting_function()` from Phase 1 FeSe script
- k-grid integration framework

**Performance**:
- Single material: ~30 min (200×200 FS + 100×100 ξ scan)
- Total: ~2-3 hours for all 5 materials

**Output structure**:
```
theory/notes/_scratch_ir/
├── _iron_based_nesting_survey_2026-09-22.py  (main script)
├── _tb_params_LaFeAsO.py   (parameter files)
├── _tb_params_BaFe2As2.py
├── _tb_params_LiFeAs.py
├── _tb_params_FeTeSe.py
└── results/
    ├── fermi_surface_5panel.png
    ├── nesting_heatmaps_5panel.png
    └── tc_correlation.png
```

---

## Timeline

**Week 1** (Days 1-7):
- Day 1-2: Literature extraction (tight-binding parameters)
- Day 3-4: Script framework + LaFeAsO test case
- Day 5-7: Compute 5 Fermi surfaces + ξ(q) for 3 materials

**Week 2** (Days 8-14):
- Day 8-9: Complete remaining 2 materials
- Day 10-11: Q_AFM identification + validation
- Day 12-13: Tc correlation analysis
- Day 14: Documentation + figure preparation

**Week 3** (Days 15-21, if needed):
- Integration with Team B findings
- Joint notebook: nesting mechanism + SCDFT vertex corrections
- Decision on publication strategy

---

## Risk Management

### Risk 1: Parameters unavailable for some materials
**Mitigation**: 
- Use 2-orbital simplified model (d_xz, d_yz)
- Mark as "simplified topology" in results
- Still captures essential Γ-M nesting

### Risk 2: Tc correlation weak
**Mitigation**:
- Universality of Q_AFM ≈ (π,0) still valuable result
- Weak correlation suggests other factors (phonons, Coulomb) dominate
- Paper focus shifts from "prediction" to "mechanism validation"

### Risk 3: Computational bottleneck
**Mitigation**:
- Reduce k-grid density (200×200 → 100×100)
- Parallelize material loop (if possible)
- Use petlove-train server for heavy lifting

---

## Success Metrics

**Minimum viable result** (guaranteed):
- 5 Fermi surface topologies computed
- Q_AFM identified for each material
- Qualitative statement: "Q_AFM universality confirmed"

**Target result** (likely):
- Tc vs ξ(Q_AFM) correlation coefficient
- Quantitative nesting enhancement factors
- Publishable in PRB as "Fermi surface nesting universality"

**Stretch goal** (if strong correlation):
- Tc = f(ξ(Q_AFM)) predictive formula
- PRL submission as "Geometric origin of Tc variation"

---

## Next Actions

1. Extract LaFeAsO parameters from Kuroki PRL 2008 (next immediate step)
2. Code framework for tight-binding diagonalization
3. Reproduce FeSe Fermi surface as validation
4. Proceed to full 5-material scan

---

**Status**: Ready to start parameter extraction  
**Blocker**: None (web search budget exhausted but parameters likely available in existing notes or local PDFs)  
**Owner**: Team C (fese-nesting-universality agent)
