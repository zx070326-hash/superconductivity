# -*- coding: utf-8 -*-
"""实测 theory/notes 下笔记数（排除 _scratch_ir/、INDEX.md、CLAUDE.md）。"""
import os, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
n = 0
per = {}
for r, d, f in os.walk('theory/notes'):
    rr = r.replace(os.sep, '/')
    if '_scratch_ir' in rr.split('/'):
        continue
    for x in f:
        if not x.endswith('.md'):
            continue
        if rr == 'theory/notes' and x in ('INDEX.md', 'CLAUDE.md'):
            continue
        n += 1
        per[rr] = per.get(rr, 0) + 1
print('total', n)
for k in sorted(per):
    print('  %-28s %d' % (k, per[k]))
