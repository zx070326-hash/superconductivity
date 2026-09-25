# -*- coding: utf-8 -*-
"""从 EuropePMC 全文 XML 逐字抽取 Wu et al. (LaH10 超快) 的承重字符串。
只做抽取与逐字打印，不做任何推断。"""
import re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

RAW = open('_wu_PMC11549371_europepmc.xml', encoding='utf-8').read()

# 去标签 -> 纯文本流（保留顺序）
txt = re.sub(r'<[^>]+>', ' ', RAW)
txt = re.sub(r'[ \t]+', ' ', txt)
txt = re.sub(r'\s*\n\s*', ' ', txt)
txt = re.sub(r'\s+', ' ', txt)

sents = re.split(r'(?<=[.]) ', txt)

def show(tag, keys, need_digit=True, limit=40, width=400):
    print('=' * 78)
    print('==', tag)
    print('=' * 78)
    n = 0
    for s in sents:
        low = s.lower()
        if any(k.lower() in low for k in keys):
            if need_digit and not re.search(r'\d', s):
                continue
            n += 1
            if n > limit:
                print('   ... (truncated)')
                break
            print(' *', s.strip()[:width])
            print()
    if n == 0:
        print('   [no match]')

# (a) 该样品自己的 Tc / 压强
show('A. 本样品的 Tc 与压强', ['superconducting transition', 'transition temperature',
                               ' Tc ', 'T c ', 'critical temperature', 'GPa'], limit=25)

# (b) lambda<Omega^2> 的两次出现 + "averaged value" 的另一条来源
show('B. lambda<Omega^2> / 2.3 / 2.4 / averaged', ['10 4', '104 (meV)', 'averaged value',
                                                    'we obtain', 'we estimate'], limit=25)

# (c) tau_fast 拟合的模型假设：<Omega>, <Omega^3>, Einstein, single mode
show('C. 拟合所用的矩 <Omega>/<Omega^3> 与模型假设',
     ['Omega 3', 'Omega 2', 'single mode', 'Einstein', 'moment', 'assume', 'Debye'], limit=30)

# (d) Allen / second moment / convention ref
show('D. 约定与 Allen 关系', ['Allen', 'convention', 'logarithmic', 'lifetime', 'two-temperature',
                              'electron-phonon coupling constant'], limit=30)

print('=' * 78)
print('== E. 关键子串的出现次数（逐字）')
print('=' * 78)
for k in ['2.3 ', '2.4 ', '96.4', '2.58', '2.67', '165', '163', '250', '260',
          'Allen', 'Omega log', 'second moment', 'Supplementary Table 1']:
    print('  %-22s x %d' % (repr(k), txt.count(k)))
