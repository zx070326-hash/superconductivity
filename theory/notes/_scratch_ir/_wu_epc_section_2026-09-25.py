# -*- coding: utf-8 -*-
"""逐字导出 Wu et al. 的 EPC 段落（按原文顺序），并定位参考文献 46/48/49 与 'second moment'。
只抽取、不推断。"""
import re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

RAW = open('_wu_PMC11549371_europepmc.xml', encoding='utf-8').read()

# 1) 保留段落边界：<p> -> 换行标记
s = RAW.replace('</p>', '\n@@P@@\n').replace('</title>', '\n@@T@@\n')
s = re.sub(r'<[^>]+>', ' ', s)
# 去掉 LaTeXML 的 documentclass 噪声，保留 $$...$$ 里的公式
s = re.sub(r'\\documentclass\[[^\]]*\]\{minimal\}', ' ', s)
s = re.sub(r'\\usepackage\{[^}]*\}', ' ', s)
s = re.sub(r'\\setlength\{[^}]*\}\{[^}]*\}', ' ', s)
s = s.replace('\\begin{document}', ' [TEX] ').replace('\\end{document}', ' [/TEX] ')
s = re.sub(r'[ \t]+', ' ', s)
lines = [re.sub(r'\s+', ' ', x).strip() for x in s.split('\n')]
paras = [x for x in lines if x and x not in ('@@P@@', '@@T@@')]

print('#' * 78)
print('# EPC 段落（含 "EPC strength" / "Allen" / "lambda" 的连续段，按原文顺序）')
print('#' * 78)
hit = [i for i, p in enumerate(paras)
       if ('Allen' in p and 'tau' in p.lower()) or 'EPC strength' in p
       or 'lambda \\langle' in p or '10}^{4}' in p]
rng = set()
for i in hit:
    for j in range(max(0, i - 1), min(len(paras), i + 2)):
        rng.add(j)
for j in sorted(rng):
    print('\n--- [para %d] ---' % j)
    print(paras[j][:2600])

print('\n' + '#' * 78)
print('# "second moment" 所在句')
print('#' * 78)
flat = re.sub(r'\s+', ' ', ' '.join(paras))
for m in re.finditer('second moment', flat):
    print(' ...', flat[max(0, m.start() - 500):m.start() + 400], '...')

print('\n' + '#' * 78)
print('# 参考文献 Allen / 编号 46,48,49 线索')
print('#' * 78)
for m in re.finditer(r'Allen[^|]{0,240}', flat):
    print(' *', m.group(0)[:240])
print()
for k in ['Supplementary Table 1', 'Supplementary Information']:
    for m in re.finditer(re.escape(k), flat):
        print(' #', flat[max(0, m.start() - 420):m.start() + 260])
        print()
