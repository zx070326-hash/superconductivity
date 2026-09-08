r"""
P3 Task C (implementer numerics, 2026-09-09) — record envelope: hydride
(P, Tc, omega_max) points from theory/data/hydride_tc_record_envelope.csv
against the Semenok stability bounds (0.20/0.32 omega_max) and the Esterlis
bound (0.08-0.12 omega_D), all converted to the omega_max convention using
the document's own order-of-magnitude ratios (P3_tc_max_pressure.md Sec.1
line 36): omega_max/omega_D ~ 2.5  =>  Esterlis 0.08-0.12 omega_D
                                        ~ 0.032-0.048 omega_max.

omega_max values are NOT a CSV column (the CSV only carries omega_log/omega_D
in meV); they are hand-extracted from the CSV's own `note` field, which itself
quotes them from Semenok et al. 2407.12922 Table S1. Each value below is
tagged with the CSV row it came from -- see theory/data/hydride_tc_record_envelope.csv
for the full citation chain. Rows with no quoted omega_max are excluded from
this plot/table (honest gap, not filled in).

Run: python theory/derivations/p3_task_C_envelope.py   (from repo root)
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

K_PER_MEV = 1 / 0.08617333262  # 1 meV = 11.6045 K

# (material, P_GPa, Tc_K, omega_max_K, source note / caveat)
# All omega_max values quoted verbatim from the CSV `note` column (Semenok
# Table S1), matched to the row with the corresponding (lambda, omega_log).
DATA = [
    ("H3S",          157, 190,   2704, "Semenok Tbl S1 row; Tc=exp@157GPa, lambda/omega=theory@200GPa (CSV row 4 caveat)"),
    ("LaH10(Semenok)",214, 245,   3145, "Semenok Tbl S1 alt row for 214 GPa point (CSV row 9 note); Tc differs from Errea theory 237.9K"),
    ("(La,Y)H10",     183, 253,   3119, "CSV row 11"),
    ("YH6",           166, 224,   2450, "CSV row 12"),
    ("YH9",           201, 243,   3167, "CSV row 14"),
    ("CaH6",          172, 215,   2877, "CSV row 15"),
    ("(La,Ce)H9-10",  123, 190,   2894, "CSV row 21"),
    ("ThH10",         170, 161,   3122, "CSV row 19"),
    ("ThH9",          150, 146,   3362, "CSV row 20"),
]

# Esterlis bound is native in omega_D; document's own order-of-magnitude
# conversion omega_max/omega_D ~ 2.5 (P3_tc_max_pressure.md Sec.1 line 36).
OMEGA_MAX_OVER_OMEGA_D = 2.5
ESTERLIS_LO_D, ESTERLIS_HI_D = 0.08, 0.12
ESTERLIS_LO_MAX = ESTERLIS_LO_D / OMEGA_MAX_OVER_OMEGA_D
ESTERLIS_HI_MAX = ESTERLIS_HI_D / OMEGA_MAX_OVER_OMEGA_D
SEMENOK_STABLE = 0.20
SEMENOK_METASTABLE = 0.32

print("=" * 90)
print("P3 record envelope: Tc/omega_max ratio vs Semenok (0.20/0.32) and Esterlis (converted) bounds")
print("=" * 90)
print(f"Esterlis 0.08-0.12 omega_D -> {ESTERLIS_LO_MAX:.4f}-{ESTERLIS_HI_MAX:.4f} omega_max "
      f"(using doc's omega_max/omega_D~{OMEGA_MAX_OVER_OMEGA_D})")
print()
hdr = f"{'material':<16}{'P(GPa)':>7}{'Tc(K)':>7}{'wmax(K)':>9}{'Tc/wmax':>9}{'vs 0.20':>9}{'vs 0.32':>9}{'vs Esterlis':>13}"
print(hdr)
rows_out = []
for name, P, Tc, wmax, note in DATA:
    ratio = Tc / wmax
    margin_stable = ratio / SEMENOK_STABLE - 1       # >0 means VIOLATES stable bound
    margin_meta = ratio / SEMENOK_METASTABLE - 1     # >0 means VIOLATES metastable bound
    if ratio < ESTERLIS_LO_MAX:
        est_flag = "below range"
    elif ratio > ESTERLIS_HI_MAX:
        est_flag = "ABOVE range"
    else:
        est_flag = "in range"
    print(f"{name:<16}{P:>7}{Tc:>7}{wmax:>9}{ratio:>9.4f}{margin_stable*100:>+8.1f}%{margin_meta*100:>+8.1f}%{est_flag:>13}")
    rows_out.append((name, P, Tc, wmax, ratio))

print()
print("Reading: 'vs 0.20'/'vs 0.32' show the % by which Tc/omega_max exceeds (positive) or falls short of")
print("(negative) the Semenok stable/metastable line. 'vs Esterlis' shows whether Tc/omega_max falls inside")
print("the converted Esterlis band -- note ALL points here are far ABOVE the Esterlis band once converted to")
print("omega_max units, consistent with P3 Sec.1's point that Esterlis (defined on omega_D) looks much tighter")
print("than Semenok (defined on omega_max) only because of the frequency-convention gap, not because it is a")
print("tighter physical constraint in the same units.")

print()
print("Which line does each material sit closest to (of the two Semenok lines)?")
for name, P, Tc, wmax, ratio in rows_out:
    d_stable = abs(ratio - SEMENOK_STABLE)
    d_meta = abs(ratio - SEMENOK_METASTABLE)
    closest = "stable (0.20)" if d_stable < d_meta else "metastable (0.32)"
    verdict = "EXCEEDS metastable (0.32) -- would violate even the looser bound!" if ratio > SEMENOK_METASTABLE else (
        "between stable and metastable" if ratio > SEMENOK_STABLE else "below stable bound")
    print(f"  {name:<16} Tc/wmax={ratio:.3f}  closest to {closest:<20} [{verdict}]")

# --- plot ---
out_dir = os.path.join(os.path.dirname(__file__), "p3_out")
os.makedirs(out_dir, exist_ok=True)
fig, ax = plt.subplots(figsize=(8, 5.5))
names = [r[0] for r in rows_out]
Ps = [r[1] for r in rows_out]
ratios = [r[4] for r in rows_out]
ax.scatter(Ps, ratios, color="#c0392b", zorder=5, s=50)
for n, p, r in zip(names, Ps, ratios):
    ax.annotate(n, (p, r), textcoords="offset points", xytext=(5, 4), fontsize=8)
ax.axhline(SEMENOK_STABLE, color="#2980b9", linestyle="--", label=f"Semenok stable (0.20 omega_max)")
ax.axhline(SEMENOK_METASTABLE, color="#8e44ad", linestyle="--", label=f"Semenok metastable (0.32 omega_max)")
ax.axhspan(ESTERLIS_LO_MAX, ESTERLIS_HI_MAX, color="#27ae60", alpha=0.2,
           label=f"Esterlis 0.08-0.12 omega_D -> {ESTERLIS_LO_MAX:.3f}-{ESTERLIS_HI_MAX:.3f} omega_max (converted)")
ax.set_xlabel("Pressure (GPa)")
ax.set_ylabel("Tc / omega_max")
ax.set_title("Hydride record envelope: Tc/omega_max vs Semenok & Esterlis bounds\n"
              "(omega_max hand-extracted from Semenok Table S1 via CSV notes, see script header)")
ax.legend(fontsize=8, loc="upper right")
ax.set_ylim(0, max(ratios) * 1.3)
fig.tight_layout()
out_path = os.path.join(out_dir, "p3_record_envelope_omega_max.png")
fig.savefig(out_path, dpi=150)
print(f"\nSaved plot: {out_path}")
