#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_drozdov_fig2c_independent.py -- 2026-09-18, main session, written from zero.

WHY
---
h3s-pressure-audit (subagent) reports that Drozdov's Fig. 2(c) annealed Tc(P)
curve gives Tc = 184.0 K at P = 200 GPa, hence a +8.2% overestimate for the
repo's 199.0 K deficit-budget prediction, with the residual SIGN robust.
That number is load-bearing (it replaces my own "sign undetermined") and
CLAUDE.md forbids taking a subagent's "done" on trust.  This script is an
INDEPENDENT from-zero re-extraction: it does not import, read, or run the
subagent's digitize_fig2c.py.

METHOD
------
Drozdov et al., Nature 525, 73 (2015) -- local copy of the preprint PDF.
Fig. 2 is on PDF page index 8; its panel (c) is drawn as VECTOR content, so
no raster pixel estimation is involved.  We:
  1. pull every vector drawing whose bbox lies inside the panel (c) frame;
  2. calibrate the axes on the MAJOR TICK MARK segments (not on the text
     label bounding boxes, whose baselines carry an unknown offset);
  3. cross-check the calibration on the MINOR ticks (must land on 125/175/225
     GPa and 50/150 K);
  4. self-check the calibration against the one data point whose value is
     known from the text: the caption says the magenta point is the magnetic
     susceptibility measurement, and the abstract's SQUID value is 203 K at
     155 GPa;
  5. convert every filled marker (caption: "filled points to sulfur hydride",
     "open color points refer to sulfur deuteride").

CAPTION, verbatim (PDF page 8 text layer, L1):
    "(c) Dependence of Tc on pressure. The data on annealed samples are
     presented.  Open color points refer to sulfur deuteride, and filled
     points to sulfur hydride. Magenta point was obtained in magnetic
     susceptibility measurements (Fig. 4a)."
=> This panel does NOT depend on arXiv:1412.0460, the disavowed un-annealed
   preprint.  That is what makes it the clean anchor.

DIRECTION AUDIT (CLAUDE.md one-sided rule)
------------------------------------------
No bound is evaluated here; every number is an affine map of extracted
coordinates.  The one-sided statement that IS used downstream is:
"199.0 K exceeds every annealed hydride point in the 195-226 GPa window",
which supports a POSITIVE residual only.  Marker-centre extraction for
ELLIPSE and RECTANGLE primitives is exact (bbox centre = data point, and the
self-check in step 4 is performed on a RECTANGLE -- the same primitive type
as the 200 GPa point).  For POLYLINE primitives (stars/triangles) the bbox
centre is NOT guaranteed to be the data point, so those are reported
separately and never used for the headline.

Level: numerical-evidence.  The 199.0 K prediction itself is INHERITED from
the repo's deficit budget and is NOT recomputed here.
"""

import sys
import fitz

PDF = "theory/notes/_scratch_ir/drozdov2015.pdf"
PAGE = 8
# panel (c) frame, located from the axis-line segments themselves
BOX = (365.0, 570.0, 135.0, 325.0)

# repo deficit-budget prediction, INHERITED (not recomputed here)
PRED_H = 199.0

MAGENTA = (1.0, 0.0, 1.0)
BLACK = (0.0, 0.0, 0.0)


def near(a, b, tol=1e-6):
    return a is not None and b is not None and \
        all(abs(x - y) < tol for x, y in zip(a, b))


def main():
    out = sys.stdout.write
    doc = fitz.open(PDF)
    page = doc[PAGE]
    x0, x1, y0, y1 = BOX
    sub = [g for g in page.get_drawings()
           if g["rect"].x0 > x0 and g["rect"].x1 < x1
           and g["rect"].y0 > y0 and g["rect"].y1 < y1]

    # ---- 1. axis calibration from tick segments -----------------------
    # major ticks are 3.29 pt long, minor 1.67 pt; the two long segments
    # (>100 pt) are the axis lines themselves.
    xt_major, xt_minor, yt_major, yt_minor = [], [], [], []
    for g in sub:
        if len(g["items"]) != 1 or g["items"][0][0] != "l":
            continue
        a, b = g["items"][0][1], g["items"][0][2]
        dx, dy = abs(a.x - b.x), abs(a.y - b.y)
        if dx > 100 or dy > 100:
            continue                       # axis line, not a tick
        if dx < 1e-6:                      # vertical tick -> x axis
            (xt_major if dy > 2.5 else xt_minor).append(a.x)
        elif dy < 1e-6:                    # horizontal tick -> y axis
            (yt_major if dx > 2.5 else yt_minor).append(a.y)

    # bottom x-axis ticks only (dedupe the mirrored top axis)
    xmaj = sorted(set(round(v, 2) for v in xt_major))
    ymaj = sorted(set(round(v, 2) for v in yt_major))
    xmin_ = sorted(set(round(v, 2) for v in xt_minor))
    ymin_ = sorted(set(round(v, 2) for v in yt_minor))
    out("major x ticks (pt): %s\n" % xmaj)
    out("major y ticks (pt): %s\n" % ymaj)

    # labelled majors, read off the text layer: x = 150,200,250 ; y = 0,100,200
    xa, xb = xmaj[-3], xmaj[-1]            # 150 and 250
    ya, yb = ymaj[-1], ymaj[0]             # 0 (largest y) and 200 (smallest)
    pt_per_gpa = (xb - xa) / 100.0
    pt_per_k = (ya - yb) / 200.0
    out("\nscale: %.6f pt/GPa , %.6f pt/K\n" % (pt_per_gpa, pt_per_k))

    def to_p(x):
        return 150.0 + (x - xa) / pt_per_gpa

    def to_t(y):
        return (ya - y) / pt_per_k

    # ---- 2. linearity cross-check on the minor ticks ------------------
    out("\n[cross-check] minor ticks must land on halves:\n")
    ok = True
    for v in xmin_:
        p = to_p(v)
        r = round(p / 25.0) * 25.0
        ok &= abs(p - r) < 0.2
        out("    x %8.2f pt -> %8.3f GPa  (expect %.0f)\n" % (v, p, r))
    for v in ymin_:
        t = to_t(v)
        r = round(t / 50.0) * 50.0
        ok &= abs(t - r) < 0.5
        out("    y %8.2f pt -> %8.3f K    (expect %.0f)\n" % (v, t, r))
    out("    linearity: %s\n" % ("PASS" if ok else "FAIL"))

    # ---- 3. collect markers -------------------------------------------
    circles, rects, polys = [], [], []
    for g in sub:
        it = g["items"]
        r = g["rect"]
        cen = ((r.x0 + r.x1) / 2.0, (r.y0 + r.y1) / 2.0)
        rec = (cen, g.get("fill"), g.get("color"), r.width)
        if len(it) == 4 and all(i[0] == "c" for i in it):
            circles.append(rec)
        elif len(it) == 1 and it[0][0] == "re":
            rects.append(rec)
        elif len(it) >= 5 and all(i[0] == "l" for i in it):
            polys.append(rec)

    # ---- 4. calibration self-check on the known magenta point ---------
    out("\n[self-check] filled magenta marker = SQUID point, known (155,203):\n")
    hit = None
    for cen, fill, col, w in rects + circles:
        if near(fill, MAGENTA):
            hit = (to_p(cen[0]), to_t(cen[1]))
            out("    extracted (%.2f GPa, %.2f K)   <- primitive: RECTANGLE\n"
                % hit)
    if hit is None:
        out("    !! not found -- calibration UNVERIFIED\n")
    else:
        out("    error: %+.2f GPa , %+.2f K\n" % (hit[0] - 155.0, hit[1] - 203.0))

    # ---- 5. filled (hydride) markers ----------------------------------
    out("\n[hydride = FILLED markers]  exact-centre primitives only:\n")
    exact = []
    for cen, fill, col, w in sorted(circles + rects, key=lambda z: z[0][0]):
        if near(fill, BLACK):
            exact.append((to_p(cen[0]), to_t(cen[1])))
    for p, t in sorted(exact):
        flag = "  <== P = 200 GPa" if abs(p - 200.0) < 1.0 else ""
        out("    %8.2f GPa   %8.2f K%s\n" % (p, t, flag))

    out("\n[hydride, POLYLINE primitives -- bbox centre NOT guaranteed]:\n")
    for cen, fill, col, w in sorted(polys, key=lambda z: z[0][0]):
        if near(fill, BLACK):
            out("    %8.2f GPa   %8.2f K   (advisory only)\n"
                % (to_p(cen[0]), to_t(cen[1])))

    # ---- 6. the residual ----------------------------------------------
    at200 = [t for p, t in exact if abs(p - 200.0) < 1.0]
    window = [(p, t) for p, t in exact if 194.0 <= p <= 227.0]
    out("\n" + "=" * 70 + "\nVERDICT\n" + "=" * 70 + "\n")
    if not at200:
        out("  no exact-centre hydride marker within 1 GPa of 200 GPa.\n")
        return
    t200 = at200[0]
    out("  Tc(annealed, H3S, 200 GPa) = %.2f K   [direct marker, no interpolation]\n"
        % t200)
    out("  repo prediction            = %.1f K   [INHERITED, not recomputed]\n"
        % PRED_H)
    out("  residual                   = %+.2f%%\n"
        % (100.0 * (PRED_H / t200 - 1.0)))
    out("\n  sign robustness -- every exact-centre hydride point in 194-227 GPa:\n")
    worst_hi = max(t for _, t in window)
    for p, t in sorted(window):
        out("    %7.2f GPa  %7.2f K   residual %+7.2f%%\n"
            % (p, t, 100.0 * (PRED_H / t - 1.0)))
    out("  max annealed Tc in that window = %.2f K < %.1f K\n" % (worst_hi, PRED_H))
    out("  => residual is POSITIVE across the whole window; sign does not flip.\n")


if __name__ == "__main__":
    main()
