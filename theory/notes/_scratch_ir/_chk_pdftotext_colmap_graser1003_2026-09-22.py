"""Reproducible "naive positional column mapping" for Graser et al. arXiv:1003.0133
Table I / Table II, as read directly out of the raw pdftotext -layout text stream
(no page-image information used here at all).

Purpose: this is the FIRST independent instrument in the two-instrument dual-check
required by CLAUDE.md Sec.2 ("抽 PDF 定承重字符时，两种仪器都要跑"). It reconstructs,
mechanically (character-position nearest-match, exactly what any naive whitespace-
aligned table reader would do), what value pdftotext's layout mode assigns to each
(row, column) cell. Cross-referencing this against the second instrument (200dpi
page-image reads, see theory/notes/literature/graser1003.0133_5orbital_TB_parameters_2026-09-22.md
Sec.4/5) is what pinned down the systematic column-misalignment bug documented in
that note's Sec.5.1 (5 mismatched cells in Table I, 23 in Table II, all caused by
sparse-column values shifting up to the previous non-blank row).

Run: PYTHONIOENCODING=utf-8 python _chk_pdftotext_colmap_graser1003_2026-09-22.py
"""
import os
import re

_HERE = os.path.dirname(__file__)
_TXT = os.path.join(_HERE, "graser1003.0133.txt")

NUM_RE = re.compile(r"[−-]?\d+\.\d+")


def header_cols(line):
    """Character offset of each column label 'i=<label>' or 'i = <label>' in a header line."""
    cols = []
    for m in re.finditer(r"i\s*=\s*([a-zA-Z]+)", line):
        cols.append((m.start(1), m.group(1)))
    return cols


def row_values(line, cols):
    """Greedy nearest-column assignment WITHOUT repeats: process tokens left-to-right
    (they appear in increasing character position, same order as the header columns),
    and for each token pick the nearest column that hasn't been used yet. A plain
    per-token argmin (no no-repeat constraint) creates spurious ties/collisions when
    two adjacent tokens are equidistant from the same header label -- that is a defect
    of a *sloppier* naive parser, not a faithful reproduction of the pdftotext
    misalignment itself, so it is corrected here."""
    tokens = [(m.start(), float(m.group(0).replace("−", "-"))) for m in NUM_RE.finditer(line)]
    remaining = list(cols)
    out = {}
    for pos, val in tokens:
        if not remaining:
            break
        col = min(remaining, key=lambda c: abs(c[0] - pos))
        remaining.remove(col)
        out[col[1]] = val
    return out


def main():
    with open(_TXT, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Table I: header on the "timm ..." line, data rows "m=1", "m=3", "m=4", "m=5"
    t1_header_idx = next(i for i, l in enumerate(lines) if l.strip().startswith("timm"))
    t1_cols = header_cols(lines[t1_header_idx])
    print("T1 header cols:", t1_cols)
    print()
    print("Table I naive positional mapping (pdftotext -layout text stream only):")
    for i in range(t1_header_idx + 1, t1_header_idx + 6):
        l = lines[i]
        m = re.match(r"\s*m\s*=\s*(\d+)", l)
        if not m:
            continue
        vals = row_values(l, t1_cols)
        print(f"  m={m.group(1)} {vals}")

    print()
    # Table II: header on the "timn ..." line, data rows "mn = 12" etc.
    t2_header_idx = next(i for i, l in enumerate(lines) if l.strip().startswith("timn"))
    t2_cols = header_cols(lines[t2_header_idx])
    print("T2 header cols:", t2_cols)
    print()
    print("Table II naive positional mapping (pdftotext -layout text stream only):")
    for i in range(t2_header_idx + 1, t2_header_idx + 9):
        l = lines[i]
        m = re.match(r"\s*mn\s*=\s*(\d+)", l)
        if not m:
            continue
        vals = row_values(l, t2_cols)
        print(f"  mn={m.group(1)} {vals}" if vals else f"  mn={m.group(1)}  (no numeric tokens -- blank line)")


if __name__ == "__main__":
    main()
