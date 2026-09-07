"""
从 Mishra et al. Figshare 数据集的 *.imag_iso_gap0（能隙边 vs T）文件提取各处理下的 Tc，
构造 H3S/D3S 的 α_H 分解表：谐波/非谐 × FSR/FBW × 绝热/非绝热 × μ*=0.16/0.12 × 200/160 GPa。
Tc 取最后一个能隙非零的温度（EPW 在能隙消失处停止扫描，网格 1–2 K）。
运行：python theory/derivations/p1_mishra_tc_table.py
"""
import glob
import os
import sys
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parents[1] / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"
LN2 = np.log(2)


def tc_of(path: str) -> float:
    d = np.loadtxt(path, comments="#")
    d = d[None, :] if d.ndim == 1 else d
    T, g = d[:, 0], d[:, 1]
    nz = T[g > 0.02]
    return float(nz.max()) if len(nz) else float("nan")


def collect():
    out = {}
    for f in glob.glob(str(ROOT / "**" / "*.imag_iso_gap0"), recursive=True):
        rel = os.path.relpath(f, ROOT).replace("\\", "/")
        parts = rel.split("/")
        mat_p = parts[0]                       # H3S-200GPa
        kind = "anh" if "anh" in parts[1].lower() else "har"
        mu = "0.12" if "muc_0.12" in rel else "0.16"
        treat = parts[-2]                      # fsr_adia ...
        mat, P = mat_p.split("-")
        out[(mat, P, kind, mu, treat)] = tc_of(f)
    return out


tcs = collect()
treats = ["fsr_adia", "fsr_na", "fbw_adia", "fbw_na"]
print("Mishra et al. (Figshare 10.6084/m9.figshare.31136635) EPW 各向同性 ME，wscut=2.5 eV，μ* 固定（同位素不变）")
print("fsr=费米面限制(常 DOS)  fbw=全带宽(能量依赖 DOS)  adia=Migdal  na=一阶顶角修正\n")
for P in ["200GPa", "160GPa"]:
    for kind in ["har", "anh"]:
        for mu in ["0.16", "0.12"]:
            rows = []
            for t in treats:
                h = tcs.get(("H3S", P, kind, mu, t)); d = tcs.get(("D3S", P, kind, mu, t))
                if h is None or d is None:
                    continue
                rows.append((t, h, d, -np.log(d / h) / LN2))
            if not rows:
                continue
            print(f"=== {P}  {kind}  μ*={mu} ===")
            print(f"{'treatment':10s} {'Tc_H (K)':>9} {'Tc_D (K)':>9} {'α_H':>7}")
            for t, h, d, a in rows:
                print(f"{t:10s} {h:9.1f} {d:9.1f} {a:7.3f}")
            r = {t: a for t, _, _, a in rows}
            if all(k in r for k in treats):
                print(f"  δα_nad(FSR) = {r['fsr_na'] - r['fsr_adia']:+.3f}   δα_nad(FBW) = {r['fbw_na'] - r['fbw_adia']:+.3f}"
                      f"   δα_DOS(adia) = {r['fbw_adia'] - r['fsr_adia']:+.3f}   δα_DOS(na) = {r['fbw_na'] - r['fsr_na']:+.3f}"
                      f"   总(fbw_na − fsr_adia) = {r['fbw_na'] - r['fsr_adia']:+.3f}")
            print()
