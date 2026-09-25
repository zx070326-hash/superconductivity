#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""主会话自查：我 2026-09-25 写下的 w2/w_log 标定是否挑了最有利的点？
消费点：ceiling/hydride_omegalog_corpus_extension_2026-09-25.md §3。
"""
import math, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
MEV_TO_K, C_INF, NB3GE_R, CM = 11.604518, 0.1827262477, 0.6260, 0.12398419843320026

# 仓库 hydrides_theory.csv 中**同时**有 w_log 与 w2 的全部行（同一次 alpha2F 积分）
CAL = [("H3S harm", 107.165, 129.32), ("H3S anh ", 109.995, 149.26),
       ("D3S harm",  79.562,  92.97), ("D3S anh ",  85.535, 105.43)]
print("== 全部可用标定点（不是只有 H3S）==")
rats = []
for tag, wl, w2 in CAL:
    r = w2 / wl; rats.append(r)
    print("  %s  w_log=%8.3f  w2=%8.3f  ratio=%.4f" % (tag, wl, w2, r))
lo, hi = min(rats), max(rats)
print("  标定范围 = [%.4f, %.4f]   （我笔记只引了 H3S 两点，并取了 max=%.4f）" % (lo, hi, hi))

ROWS = [("LaH10", 129, 3.62, 76.4, 1, 250.0), ("LaH10", 163, 2.67, 96.4, 1, 250.0),
        ("LaH10", 214, 2.06, 115.5, 1, 250.0), ("LaH10", 264, 1.73, 126.6, 1, 250.0),
        ("LaH10", 163, 2.67, 96.4, 1, 260.0), ("LaD10", 159, 3.14, 63.5, 1, 168.0),
        ("LaD10", 210, 2.21, 81.7, 1, 168.0), ("LaD10", 260, 1.80, 92.2, 1, 168.0),
        ("CaH6", 150, 2.75, 684.0, CM, 215.0), ("CaH6", 200, 1.99, 851.0, CM, 215.0),
        ("CaH6", 250, 1.71, 933.0, CM, 215.0)]
print()
print("== 每行：夺冠所需阈值，以及 4 个标定点各自给出的 r_true 判决 ==")
print("  %-6s %4s %8s %9s | %s" % ("体系","P","r_up","需比值","4 个标定点下 r_true（>0.626 = 夺冠）"))
n_flip = 0
for name,P,lam,wl,conv,tc in ROWS:
    wlK = wl*conv*MEV_TO_K
    r_up = tc/(C_INF*math.sqrt(lam)*wlK)
    need = r_up/NB3GE_R
    cells=[]; flips=0
    for rat in rats:
        rt = r_up/rat
        cells.append("%.3f%s" % (rt, "*" if rt>NB3GE_R else " "))
        if rt>NB3GE_R: flips+=1
    if flips: n_flip+=1
    print("  %-6s %4d %8.4f %9.4f | %s   (%d/4 夺冠)" % (name,P,r_up,need," ".join(cells),flips))
print()
print("== 判决 ==")
print("  夺冠所需阈值的最大值 = %.4f" % max(r_up/NB3GE_R for name,P,lam,wl,conv,tc in ROWS
      for r_up in [tc/(C_INF*math.sqrt(lam)*wl*conv*MEV_TO_K)]))
print("  标定范围 [%.4f, %.4f] **跨过** 1.2957 ⟹ 标定判不了。" % (lo,hi))
print("  有 %d 行在至少一个标定点下夺冠 ⟹ 我笔记 §3『若与 H3S 同量级则全部不夺冠』**方向取反了**：" % n_flip)
print("  比值**取大** ⟹ r_true 变小 ⟹ 偏向『不夺冠』。我取 max=1.3570 是单边取值，")
print("  却把它陈述成支持『不夺冠』的证据 —— 违反单边界方向规则。")
