"""
Team-lead's high-priority follow-up (2026-09-22): does the masking-radius bug found
in the BaFe2As2 fine-grid rescan (see _chk_nesting_localmax_scan_2026-09-22.py) also
contaminate the Q_AFM claims for the OTHER 4 materials (LaFeAsO, LiFeAs, FeSe,
FeTeSe) reported by _iron_based_nesting_survey_2026-09-22.py
(-> materials/iron_based_Q_AFM_universality_2026-09-22.md,
    materials/iron_based_tight_binding_params_2026-09-22.md)?

STEP 1 (confirmed by reading the code, not the notes' prose):
  scan_nesting_bz() is called ONCE per material inside a single shared loop
  (run_survey(), line ~197), always with the SAME default q_exclude_radius=0.15
  (units of pi) -- never overridden per material.  So yes: same masking logic /
  same code path for all 5 materials.

  IMPORTANT DIFFERENCE from the BaFe2As2 fine-grid bug: that scan used q_grid=128
  (spacing ~0.0156pi) with an exclusion radius of 3 grid points (~0.047pi), so the
  disk covered several neighboring grid points and a genuine radius-in-grid-points
  sweep was meaningful.  THIS script uses q_grid=15 over [0,pi], spacing
  pi/14 ~= 0.2244pi -- LARGER than the 0.15pi exclusion radius itself.  So the
  mask here can only ever remove the single exact node q=(0,0); it is not a
  multi-point disk on this coarse grid.  The relevant question is therefore
  different: is the reported argmax (the first ring of surviving grid points)
  still sitting on the smooth shoulder/tail of the broad Gamma self-peak, i.e.
  does the argmax value keep dropping / does the location keep moving outward
  as we additionally exclude the first ring, second ring, etc. (in GRID-POINT
  units, not pi units, since pi-unit radius is coarser than the grid itself)?

Diagnostic: for each material, build the full xi_map on this script's own
(coarse) grid, then sweep exclusion radius in units of GRID POINTS
(0,1,2,3,...,7 covers up to half of q_grid=15), each time masking a disk of
that many grid-index steps around (0,0) (Euclidean distance in *index* units,
not radians) and recording argmax location + value.  Monotonic retreat
(location moving outward, value dropping) as radius grows = contaminated by
the same bug.  Stabilization (same location/value across several radii) =
diagnostic did not detect an issue (NOT the same as "confirmed clean" --
one-sided bound direction, per CLAUDE.md).
"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['PYTHONIOENCODING'] = 'utf-8'

import numpy as np
import importlib.util

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "_iron_based_nesting_survey_2026-09-22",
    os.path.join(SCRIPT_DIR, "_iron_based_nesting_survey_2026-09-22.py"))
survey = importlib.util.module_from_spec(spec)
# prevent the module's own __main__ block / run_survey() print-spam from firing on import
spec.loader.exec_module(survey)

MATERIALS = ['LaFeAsO', 'BaFe2As2', 'LiFeAs', 'FeSe', 'FeTeSe']

print("=" * 100)
print("Per-material masking-radius-in-GRID-POINTS diagnostic")
print("(reusing the survey script's own build_2orbital_hamiltonian / "
      "compute_nesting_function_2orbital / material params verbatim)")
print("=" * 100)

results_table = []

for material in MATERIALS:
    params, Tc = survey.get_material_hamiltonian_params(material)

    # rebuild the SAME xi_map the production script builds (q_grid=15, k_grid=100,
    # fermi_window=0.08 -- exact args used inside run_survey())
    q_grid = 15
    k_grid = 100
    fermi_window = 0.08
    qx_arr = np.linspace(0, np.pi, q_grid)
    qy_arr = np.linspace(0, np.pi, q_grid)
    xi_map = np.zeros((q_grid, q_grid))
    for i, qx in enumerate(qx_arr):
        for j, qy in enumerate(qy_arr):
            xi_map[i, j] = survey.compute_nesting_function_2orbital(
                (qx, qy), params, fermi_window, k_grid)

    # old reported value (radius=0.15pi in RADIAN units, effectively "exclude only (0,0)")
    q_exclude = 0.15 * np.pi
    QX, QY = np.meshgrid(qx_arr, qy_arr, indexing='ij')
    mask_old = (QX**2 + QY**2) < q_exclude**2
    xi_masked_old = xi_map.copy()
    xi_masked_old[mask_old] = -np.inf
    idx_old = np.unravel_index(np.argmax(xi_masked_old), xi_masked_old.shape)
    q_old = (qx_arr[idx_old[0]] / np.pi, qy_arr[idx_old[1]] / np.pi)
    xi_old = xi_map[idx_old]

    print(f"\n[{material}]  old reported (radius=0.15pi, radian units): "
          f"Q_AFM=({q_old[0]:.3f}pi,{q_old[1]:.3f}pi)  xi={xi_old:.4f}  "
          f"n_masked_nodes={mask_old.sum()}")

    # radius sweep in GRID-POINT (index) units: r=0 excludes only (0,0); r=1 excludes
    # the 3x3 index-neighborhood; etc. up to r=7 (half of q_grid=15)
    ii, jj = np.meshgrid(np.arange(q_grid), np.arange(q_grid), indexing='ij')
    seq = []
    for r in range(0, 8):
        mask_r = (ii**2 + jj**2) <= r**2  # distance from index (0,0), i.e. from q=(0,0)
        xm = xi_map.copy()
        xm[mask_r] = -np.inf
        if np.all(~np.isfinite(xm)):
            seq.append((r, None, None, mask_r.sum()))
            continue
        idx = np.unravel_index(np.argmax(xm), xm.shape)
        q_loc = (qx_arr[idx[0]] / np.pi, qy_arr[idx[1]] / np.pi)
        xi_val = xi_map[idx]
        seq.append((r, q_loc, xi_val, mask_r.sum()))

    print(f"  radius-sweep (grid-point units, r=0..7):")
    for r, q_loc, xi_val, nmask in seq:
        if q_loc is None:
            print(f"    r={r}: ALL MASKED (nmask={nmask})")
        else:
            print(f"    r={r}: argmax=({q_loc[0]:.3f}pi,{q_loc[1]:.3f}pi)  "
                  f"xi={xi_val:.4f}  nmask={nmask}")

    # classify: monotonic retreat = value strictly decreasing AND location moving
    # outward (index distance from origin non-decreasing) across the sweep, while
    # still changing (not already stable at r=0)
    finite_seq = [(r, q, v, n) for r, q, v, n in seq if q is not None]
    values = [v for _, _, v, _ in finite_seq]
    locs = [q for _, q, _, _ in finite_seq]
    # location changes at all across r=0..7?
    unique_locs = []
    for l in locs:
        if not unique_locs or unique_locs[-1] != l:
            unique_locs.append(l)
    n_transitions = len(unique_locs) - 1
    # monotonic value decrease check (strictly non-increasing, with at least one
    # strict decrease) up to the point of stabilization
    strictly_dropping = all(values[k] >= values[k+1] - 1e-9 for k in range(len(values)-1))
    ever_drops = any(values[k] > values[k+1] + 1e-9 for k in range(len(values)-1))

    if n_transitions >= 1 and ever_drops:
        verdict = "HIT (argmax location/value changes with radius -- retreat detected)"
    elif n_transitions == 0:
        verdict = "NOT DETECTED (argmax stable at r=0..7 -- cannot rule out a peak coincidentally near Gamma)"
    else:
        verdict = "AMBIGUOUS (location changes but value non-monotonic -- needs manual look)"

    print(f"  n_location_transitions(r=0..7) = {n_transitions}, "
          f"monotonic_nonincreasing = {strictly_dropping}")
    print(f"  VERDICT: {verdict}")

    # stabilized "true" value = value at the largest radius tested (r=7, ~half BZ)
    r_final, q_final, xi_final, _ = finite_seq[-1]

    results_table.append({
        'material': material,
        'q_old': q_old, 'xi_old': xi_old,
        'seq_locs': unique_locs,
        'q_final': q_final, 'xi_final': xi_final,
        'verdict': verdict,
    })

print("\n" + "=" * 100)
print("SUMMARY TABLE")
print("=" * 100)
print(f"{'Material':<12} {'old Q_AFM':<16} {'old xi':<10} {'stabilized(r=7)':<18} "
      f"{'xi(r=7)':<10} {'verdict'}")
for r in results_table:
    print(f"{r['material']:<12} "
          f"({r['q_old'][0]:.2f},{r['q_old'][1]:.2f})pi   "
          f"{r['xi_old']:<10.4f} "
          f"({r['q_final'][0]:.2f},{r['q_final'][1]:.2f})pi      "
          f"{r['xi_final']:<10.4f} "
          f"{r['verdict']}")
