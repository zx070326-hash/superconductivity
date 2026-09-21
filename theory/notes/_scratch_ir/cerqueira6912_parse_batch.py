"""
Parse McMillan.dat / Eliashberg.dat / lambda files from an extracted Cerqueira
et al. (arXiv:2307.10728, Materials Cloud 3kbt5-r3n56) batch directory into a CSV.

Source of truth for units/definitions (read directly off McMillan.dat headers,
cross-checked against arXiv:2307.10728 Eq.(1a-c) text fetched 2026-09-18):
  lambda        = 2*int[a2F(w)/w] dw                (dimensionless)
  wlog [K,meV]  = exp[ (2/lambda) int (a2F(w)/w) ln(w) dw ]
  w2av [K,meV]  = sqrt[ (2/lambda) int a2F(w)*w dw ]   <- this IS omega-bar-2
                  (paper calls it <omega^2>^(1/2); McMillan.dat literally
                  prints "w2av[K]" / "w2av[meV]" per material, so no manual
                  a2F.dos integration is needed -- contrary to the pre-download
                  assumption in omega2_data_sources_2026-09-17.md Q2(a), which
                  only inspected files_description.md, not actual file contents)
  mu* grid      = 0.000, 0.020, ..., step 0.02 (McMillan.dat / Eliashberg.dat)
  Tc columns    = TcMcMillan, TcAllenDynes (McMillan.dat), TcEliashberg (Eliashberg.dat)
                  '-----' means solver did not converge / no solution at that mu*
                  (typically because lambda - mu*(1+0.62*lambda) <= 0)

Each material directory name is "<Formula>_<materials-cloud-id>".

Output columns:
  material, mcid, lambda, wlog_K, wlog_meV, w2av_K, w2av_meV,
  Tc_AD_mustar010_K, Tc_McM_mustar010_K, Tc_Eliashberg_mustar010_K
"""
import sys, os, re, csv, glob

def parse_mcmillan(path):
    d = {}
    tc_rows = []
    with open(path, errors="replace") as f:
        for line in f:
            line = line.strip()
            m = re.match(r"#\s*lambda\s*=\s*([\-0-9.Ee+]+)", line)
            if m: d["lambda"] = float(m.group(1))
            m = re.match(r"#\s*wlog\[K\]\s*=\s*([\-0-9.Ee+]+)", line)
            if m: d["wlog_K"] = float(m.group(1))
            m = re.match(r"#\s*w2av\[K\]\s*=\s*([\-0-9.Ee+]+)", line)
            if m: d["w2av_K"] = float(m.group(1))
            m = re.match(r"#\s*wlog\[mev\]\s*=\s*([\-0-9.Ee+]+)", line, re.I)
            if m: d["wlog_meV"] = float(m.group(1))
            m = re.match(r"#\s*w2av\[mev\]\s*=\s*([\-0-9.Ee+]+)", line, re.I)
            if m: d["w2av_meV"] = float(m.group(1))
            m = re.match(r"^([\-0-9.]+)\s+(\-\-\-\-\-|[\-0-9.Ee+]+)\s+(\-\-\-\-\-|[\-0-9.Ee+]+)$", line)
            if m:
                mustar = float(m.group(1))
                tcm = None if m.group(2) == "-----" else float(m.group(2))
                tcad = None if m.group(3) == "-----" else float(m.group(3))
                tc_rows.append((mustar, tcm, tcad))
    return d, tc_rows

def parse_eliashberg(path):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, errors="replace") as f:
        for line in f:
            line = line.strip()
            m = re.match(r"^([\-0-9.]+)\s+(\-\-\-\-\-|[\-0-9.Ee+]+)$", line)
            if m:
                mustar = float(m.group(1))
                tce = None if m.group(2) == "-----" else float(m.group(2))
                rows.append((mustar, tce))
    return rows

def nearest(rows, target, idx_val=0):
    best = None
    for row in rows:
        if best is None or abs(row[0] - target) < abs(best[0] - target):
            best = row
    return best

def main(batch_dir, out_csv):
    mat_dirs = sorted(glob.glob(os.path.join(batch_dir, "*")))
    n_ok, n_skip = 0, 0
    with open(out_csv, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["material_dir","formula","mcid","lambda","wlog_K","wlog_meV",
                    "w2av_K","w2av_meV","r_w2av_over_wlog",
                    "Tc_McM_mustar010_K","Tc_AD_mustar010_K","Tc_Eliashberg_mustar010_K"])
        for md in mat_dirs:
            base = os.path.basename(md)
            mcm_path = os.path.join(md, "McMillan.dat")
            eli_path = os.path.join(md, "Eliashberg.dat")
            if not os.path.exists(mcm_path):
                n_skip += 1
                continue
            try:
                hdr, tc_rows = parse_mcmillan(mcm_path)
                eli_rows = parse_eliashberg(eli_path)
            except Exception as e:
                n_skip += 1
                continue
            if "lambda" not in hdr or "w2av_K" not in hdr or "wlog_K" not in hdr:
                n_skip += 1
                continue
            row010 = nearest(tc_rows, 0.10) if tc_rows else None
            eli010 = nearest(eli_rows, 0.10) if eli_rows else None
            if base.count("_") >= 1:
                formula, mcid = base.rsplit("_", 1)
            else:
                formula, mcid = base, ""
            r = hdr["w2av_K"]/hdr["wlog_K"] if hdr.get("wlog_K") else None
            w.writerow([base, formula, mcid,
                        hdr.get("lambda"), hdr.get("wlog_K"), hdr.get("wlog_meV"),
                        hdr.get("w2av_K"), hdr.get("w2av_meV"), r,
                        row010[1] if row010 else None,
                        row010[2] if row010 else None,
                        eli010[1] if eli010 else None])
            n_ok += 1
    print(f"parsed {n_ok} materials, skipped {n_skip} (incomplete/missing McMillan.dat)")

if __name__ == "__main__":
    batch_dir = sys.argv[1] if len(sys.argv) > 1 else \
        r"D:\physics\superconductivity\theory\data\external\cerqueira_raw\extracted\batch-q"
    out_csv = sys.argv[2] if len(sys.argv) > 2 else \
        r"D:\physics\superconductivity\theory\data\external\cerqueira6912_batchq_sample.csv"
    main(batch_dir, out_csv)
