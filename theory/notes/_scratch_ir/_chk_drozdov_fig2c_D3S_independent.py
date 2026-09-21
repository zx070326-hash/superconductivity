#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_drozdov_fig2c_D3S_independent.py -- 2026-09-18, main session, from zero.

WHY
---
oneliner-fix 的 f_obs 用了 D3S(200 GPa, annealed) = 145.5 K 作为承重锚点。
那个数是 `h3s_tc_exp_pressure_audit_2026-09-18.md` 第 65 行报的「145/146 两个相邻
标记」的中点，出自**另一个会话**的 digitize_fig2c.py。我上一轮只独立复核了 H3S 侧
三个点（189.99 / 184.04 / 170.04 K），**D3S 侧我没核过**。

而我上一轮恰好确立了一条筛选纪律：只有 RECTANGLE 与闭合 BEZIER（椭圆）图元的
bbox 中心才精确等于数据点；POLYLINE（星/三角）不是。D3S 是「空心彩色点」，
必须过同一道检验才能承重。

本脚本是从零的独立再提取：不 import / 不读 / 不运行 digitize_fig2c.py。
定标与 `_chk_drozdov_fig2c_independent.py` 同法（刻度线段定标 + 洋红已知点自检），
但本文件只关心**空心彩色**（deuteride）标记。

图注逐字（PDF 第 8 页文本层，L1）：
    "Open color points refer to sulfur deuteride, and filled points to
     sulfur hydride."
=> 判别规则：fill 为空（或与描边色不同的白）且 color 非黑 ⟹ deuteride。

Level: numerical-evidence（矢量坐标的仿射映射；无界、无方向判决）。
"""

import sys
import fitz

PDF = "theory/notes/_scratch_ir/drozdov2015.pdf"
PAGE = 8
BOX = (365.0, 570.0, 135.0, 325.0)
MAGENTA = (1.0, 0.0, 1.0)
BLACK = (0.0, 0.0, 0.0)


def near(a, b, tol=1e-6):
    return a is not None and b is not None and \
        all(abs(x - y) < tol for x, y in zip(a, b))


def main():
    out = sys.stdout.write
    page = fitz.open(PDF)[PAGE]
    x0, x1, y0, y1 = BOX
    sub = [g for g in page.get_drawings()
           if g["rect"].x0 > x0 and g["rect"].x1 < x1
           and g["rect"].y0 > y0 and g["rect"].y1 < y1]

    xt_major, yt_major = [], []
    for g in sub:
        if len(g["items"]) != 1 or g["items"][0][0] != "l":
            continue
        a, b = g["items"][0][1], g["items"][0][2]
        dx, dy = abs(a.x - b.x), abs(a.y - b.y)
        if dx > 100 or dy > 100:
            continue
        if dx < 1e-6 and dy > 2.5:
            xt_major.append(a.x)
        elif dy < 1e-6 and dx > 2.5:
            yt_major.append(a.y)
    xmaj = sorted(set(round(v, 2) for v in xt_major))
    ymaj = sorted(set(round(v, 2) for v in yt_major))
    xa, xb = xmaj[-3], xmaj[-1]
    ya, yb = ymaj[-1], ymaj[0]
    pt_per_gpa = (xb - xa) / 100.0
    pt_per_k = (ya - yb) / 200.0
    to_p = lambda x: 150.0 + (x - xa) / pt_per_gpa
    to_t = lambda y: (ya - y) / pt_per_k

    circles, rects, polys = [], [], []
    for g in sub:
        it, r = g["items"], g["rect"]
        cen = ((r.x0 + r.x1) / 2.0, (r.y0 + r.y1) / 2.0)
        rec = (cen, g.get("fill"), g.get("color"), r.width, r.height)
        if len(it) == 4 and all(i[0] == "c" for i in it):
            circles.append(("ELLIPSE", rec))
        elif len(it) == 1 and it[0][0] == "re":
            rects.append(("RECT", rec))
        elif len(it) >= 5 and all(i[0] == "l" for i in it):
            polys.append(("POLY", rec))

    out("[self-check] 洋红已知点 (155 GPa, 203 K):\n")
    for kind, (cen, fill, col, w, h) in rects + circles:
        if near(fill, MAGENTA):
            p, t = to_p(cen[0]), to_t(cen[1])
            out("    (%.2f GPa, %.2f K)  误差 %+.2f / %+.2f   图元=%s\n"
                % (p, t, p - 155.0, t - 203.0, kind))

    out("\n[全部标记：按 fill/color 分类，找空心彩色 = deuteride]\n")
    out("  %-8s %10s %10s  %-22s %-22s %6s\n"
        % ("图元", "P/GPa", "Tc/K", "fill", "stroke", "宽/pt"))
    d3s_exact, d3s_poly = [], []
    for kind, (cen, fill, col, w, h) in sorted(
            circles + rects + polys, key=lambda z: z[1][0][0]):
        p, t = to_p(cen[0]), to_t(cen[1])
        filled_black = near(fill, BLACK)
        is_open_color = (fill is None) and (col is not None) and not near(col, BLACK)
        if not is_open_color:
            continue
        out("  %-8s %10.2f %10.2f  %-22s %-22s %6.2f\n"
            % (kind, p, t, str(fill), str(tuple(round(c, 3) for c in col)), w))
        (d3s_exact if kind in ("RECT", "ELLIPSE") else d3s_poly).append((kind, p, t))

    out("\n" + "=" * 72 + "\n判决\n" + "=" * 72 + "\n")
    at200 = [(k, p, t) for k, p, t in d3s_exact if abs(p - 200.0) < 1.5]
    out("  精确中心图元（RECT/ELLIPSE）在 200±1.5 GPa 内的 deuteride 点：\n")
    for k, p, t in sorted(at200, key=lambda z: z[1]):
        out("      %s  %.2f GPa  %.2f K\n" % (k, p, t))
    if d3s_poly:
        out("  POLYLINE 图元（bbox 中心不保证=数据点，仅供参考，不承重）：\n")
        for k, p, t in sorted(d3s_poly, key=lambda z: z[1]):
            out("      %s  %.2f GPa  %.2f K\n" % (k, p, t))
    if at200:
        ts = sorted(t for _, _, t in at200)
        out("\n  ⟹ 我独立读到的 D3S@200GPa = %s K\n"
            % " / ".join("%.2f" % t for t in ts))
        out("  ⟹ 对照 oneliner-fix 用的 145.5 K（另一会话 145/146 的中点）:"
            " 差 %+.2f K\n" % (sum(ts) / len(ts) - 145.5))
    else:
        out("  ⟹ 200 GPa 附近**没有**精确中心的 deuteride 图元 ⟹ 145.5 K 不可承重\n")


if __name__ == "__main__":
    main()
