#!/usr/bin/env python3
"""
Verify batch-q parsed output against cerqueira6912_batchq_sample.csv
Report any discrepancies in lambda, wlog, w2av (must match to full precision).
"""
import sys
import csv
from pathlib import Path

def load_csv(fpath):
    """Load CSV into list of dicts."""
    with open(fpath, 'r', encoding='utf-8') as f:
        return list(csv.DictReader(f))

def compare_float(val1_str, val2_str, field, mat_dir, tol=1e-15):
    """Compare two float strings, return True if match within tolerance."""
    if not val1_str or not val2_str:
        return val1_str == val2_str, f"{mat_dir}: {field} empty mismatch: '{val1_str}' vs '{val2_str}'"

    try:
        v1 = float(val1_str)
        v2 = float(val2_str)
        if abs(v1 - v2) > tol * max(abs(v1), abs(v2), 1):
            return False, f"{mat_dir}: {field} mismatch: {v1:.17g} vs {v2:.17g} (diff={v1-v2:.3e})"
        return True, None
    except ValueError as e:
        return False, f"{mat_dir}: {field} parse error: '{val1_str}' vs '{val2_str}' ({e})"

def main():
    ref_path = Path("D:/physics/superconductivity/theory/data/external/cerqueira6912_batchq_sample.csv")
    test_path = Path("batch-q-parsed.csv")

    if not ref_path.exists():
        print(f"ERROR: Reference file not found: {ref_path}", file=sys.stderr)
        sys.exit(1)

    if not test_path.exists():
        print(f"ERROR: Test file not found: {test_path}", file=sys.stderr)
        sys.exit(1)

    ref_data = load_csv(ref_path)
    test_data = load_csv(test_path)

    print(f"Reference: {len(ref_data)} rows (excluding header)")
    print(f"Test:      {len(test_data)} rows (excluding header)")

    if len(ref_data) != len(test_data):
        print(f"ERROR: Row count mismatch!", file=sys.stderr)
        sys.exit(1)

    # Build lookup by material_dir
    ref_lookup = {r['material_dir']: r for r in ref_data}
    test_lookup = {r['material_dir']: r for r in test_data}

    fields_to_check = ['lambda', 'wlog_K', 'wlog_meV', 'w2av_K', 'w2av_meV', 'r_w2av_over_wlog']
    test_field_map = {
        'r_w2av_over_wlog': 'r'
    }

    mismatches = []

    for mat_dir in sorted(ref_lookup.keys()):
        if mat_dir not in test_lookup:
            mismatches.append(f"{mat_dir}: MISSING in test output")
            continue

        ref_row = ref_lookup[mat_dir]
        test_row = test_lookup[mat_dir]

        for field in fields_to_check:
            test_field = test_field_map.get(field, field)
            ref_val = ref_row.get(field, '')
            test_val = test_row.get(test_field, '')

            match, msg = compare_float(ref_val, test_val, field, mat_dir)
            if not match:
                mismatches.append(msg)

    if mismatches:
        print(f"\n{len(mismatches)} MISMATCHES FOUND:")
        for m in mismatches[:20]:  # Show first 20
            print(f"  {m}")
        if len(mismatches) > 20:
            print(f"  ... and {len(mismatches)-20} more")
        sys.exit(1)
    else:
        print(f"\nVERIFICATION PASSED: All {len(ref_data)} materials match reference to full precision.")

if __name__ == '__main__':
    main()
