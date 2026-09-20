#!/usr/bin/env python3
"""
Cerqueira 6912 data parser: extract lambda, wlog, w2av, Tc from McMillan.dat and Eliashberg.dat
Outputs CSV: batch, material_dir, formula, mcid, lambda, wlog_K, wlog_meV, w2av_K, w2av_meV, r,
             Tc_ME_mu000, Tc_ME_mu010, Tc_ME_mu013, Tc_McM_mu010, Tc_AD_mu010, n_ME_numeric
"""
import sys
import os
import re
from pathlib import Path

def parse_mcmillan_dat(fpath):
    """Parse McMillan.dat for lambda, wlog, w2av, and Tc values at mu*=0.10 grid point."""
    try:
        with open(fpath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except Exception as e:
        return None

    # Extract header values
    lam = wlog_K = wlog_meV = w2av_K = w2av_meV = None
    for line in lines:
        if '# lambda' in line:
            lam = float(line.split('=')[1].strip())
        elif '# wlog[K]' in line:
            wlog_K = float(line.split('=')[1].strip())
        elif '# wlog[mev]' in line or '# wlog[meV]' in line:
            wlog_meV = float(line.split('=')[1].strip())
        elif '# w2av[K]' in line:
            w2av_K = float(line.split('=')[1].strip())
        elif '# w2av[meV]' in line:
            w2av_meV = float(line.split('=')[1].strip())

    if lam is None:
        return None

    # Extract Tc at mu*=0.10
    Tc_McM_mu010 = Tc_AD_mu010 = None
    for line in lines:
        if line.strip().startswith('0.100'):
            parts = line.split()
            if len(parts) >= 3:
                try:
                    val_mcm = parts[1].strip()
                    val_ad = parts[2].strip()
                    Tc_McM_mu010 = None if val_mcm == '-----' else float(val_mcm)
                    Tc_AD_mu010 = None if val_ad == '-----' else float(val_ad)
                except:
                    pass
            break

    r = w2av_K / wlog_K if (wlog_K and wlog_K > 0) else None

    return {
        'lambda': lam,
        'wlog_K': wlog_K,
        'wlog_meV': wlog_meV,
        'w2av_K': w2av_K,
        'w2av_meV': w2av_meV,
        'r': r,
        'Tc_McM_mu010': Tc_McM_mu010,
        'Tc_AD_mu010': Tc_AD_mu010,
    }

def parse_eliashberg_dat(fpath):
    """Parse Eliashberg.dat for mu* grid and Tc values at mu*=0.00, 0.10, 0.13."""
    try:
        with open(fpath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except Exception as e:
        return None

    mu_grid = []
    Tc_vals = []
    n_numeric = 0

    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                mu = float(parts[0])
                val_str = parts[1].strip()
                if val_str == '-----':
                    Tc_val = None
                else:
                    Tc_val = float(val_str)
                    n_numeric += 1
                mu_grid.append(mu)
                Tc_vals.append(Tc_val)
            except:
                pass

    # Extract specific mu* values
    Tc_ME_mu000 = Tc_ME_mu010 = Tc_ME_mu013 = None
    for mu, tc in zip(mu_grid, Tc_vals):
        if abs(mu - 0.0) < 1e-6:
            Tc_ME_mu000 = tc
        elif abs(mu - 0.10) < 1e-6:
            Tc_ME_mu010 = tc
        elif abs(mu - 0.13) < 1e-6:
            Tc_ME_mu013 = tc

    mu_step = mu_grid[1] - mu_grid[0] if len(mu_grid) > 1 else None
    mu_max = mu_grid[-1] if mu_grid else None

    return {
        'Tc_ME_mu000': Tc_ME_mu000,
        'Tc_ME_mu010': Tc_ME_mu010,
        'Tc_ME_mu013': Tc_ME_mu013,
        'n_ME_numeric': n_numeric,
        'mu_step': mu_step,
        'mu_max': mu_max,
    }

def extract_formula_mcid(material_dir):
    """Extract formula and mcid from directory name: <formula>_<mcid>."""
    parts = material_dir.rsplit('_', 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return material_dir, ''

def collect_batch(batch_root, batch_name):
    """Collect all materials in a batch directory."""
    batch_path = Path(batch_root)
    if not batch_path.exists():
        return []

    results = []
    for mat_dir in sorted(batch_path.iterdir()):
        if not mat_dir.is_dir():
            continue

        mcmillan_file = mat_dir / 'McMillan.dat'
        eliashberg_file = mat_dir / 'Eliashberg.dat'

        if not mcmillan_file.exists():
            continue

        mc_data = parse_mcmillan_dat(mcmillan_file)
        if mc_data is None:
            continue

        eli_data = parse_eliashberg_dat(eliashberg_file) if eliashberg_file.exists() else {}

        formula, mcid = extract_formula_mcid(mat_dir.name)

        row = {
            'batch': batch_name,
            'material_dir': mat_dir.name,
            'formula': formula,
            'mcid': mcid,
            **mc_data,
            **eli_data,
        }
        results.append(row)

    return results

def main():
    if len(sys.argv) < 2:
        print("Usage: python cerqueira_collect.py <root_dir> [<batch_name>]", file=sys.stderr)
        print("  root_dir: directory containing material subdirectories", file=sys.stderr)
        print("  batch_name: optional batch identifier (default: inferred from path)", file=sys.stderr)
        sys.exit(1)

    root_dir = sys.argv[1]
    batch_name = sys.argv[2] if len(sys.argv) > 2 else Path(root_dir).name

    results = collect_batch(root_dir, batch_name)

    if not results:
        print("No valid materials found.", file=sys.stderr)
        sys.exit(1)

    # Print CSV header
    fields = ['batch', 'material_dir', 'formula', 'mcid', 'lambda', 'wlog_K', 'wlog_meV',
              'w2av_K', 'w2av_meV', 'r', 'Tc_ME_mu000', 'Tc_ME_mu010', 'Tc_ME_mu013',
              'Tc_McM_mu010', 'Tc_AD_mu010', 'n_ME_numeric']
    print(','.join(fields))

    # Print data rows
    for row in results:
        vals = []
        for f in fields:
            v = row.get(f, None)
            if v is None:
                vals.append('')
            elif isinstance(v, float):
                vals.append(f'{v:.17g}')
            else:
                vals.append(str(v))
        print(','.join(vals))

    # Print diagnostic info to stderr
    print(f"\n[INFO] Processed {len(results)} materials from batch '{batch_name}'", file=sys.stderr)
    if results:
        eli_data = [r for r in results if r.get('n_ME_numeric', 0) > 0]
        print(f"[INFO] Eliashberg data available: {len(eli_data)}/{len(results)} materials", file=sys.stderr)
        if eli_data:
            mu_steps = [r.get('mu_step') for r in eli_data if r.get('mu_step') is not None]
            mu_maxs = [r.get('mu_max') for r in eli_data if r.get('mu_max') is not None]
            if mu_steps:
                print(f"[INFO] mu* step: {min(mu_steps):.3f} to {max(mu_steps):.3f}", file=sys.stderr)
            if mu_maxs:
                print(f"[INFO] mu* max: {min(mu_maxs):.3f} to {max(mu_maxs):.3f}", file=sys.stderr)

if __name__ == '__main__':
    main()
