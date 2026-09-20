# -*- coding: utf-8 -*-
"""agent 的 deficit_decomposition_analysis §1.2 把 AD75 的 f2 写成
       f2 = 1 + [ (r-1) lam^2/(lam^2+L2^2) ]^(1/2)     <-- 多了一个开方
   AD75 Eq.(35) 实为
       f2 = 1 +   (r-1) lam^2/(lam^2+L2^2)
   本脚本判断：这一处是否就是它报 13.2% 误差（我算 3.85%）的原因。
"""
import numpy as np
exec(open("theory/notes/_scratch_ir/_chk_main_lambda_exponent_artifact_2026-09-20.py",
          encoding="utf-8").read().split('print("=" * 78)')[0])

def ad(lam, wlog, wb2, mus, sqrt_bug):
    L1 = 2.46*(1+3.8*mus); r = wb2/wlog; L2 = 1.82*(1+6.3*mus)*r
    f1 = (1+(lam/L1)**1.5)**(1/3)
    inner = (r-1)*lam**2/(lam**2+L2**2)
    f2 = 1 + (np.sqrt(np.abs(inner)) if sqrt_bug else inner)
    return f1*f2*(wlog/1.20)*np.exp(-1.04*(1+lam)/(lam-mus*(1+0.62*lam)))

for bug in (False, True):
    t = ad(lam, wlog, wb2, mus, bug)
    rel = np.abs(t/Tc-1)
    tag = "带开方 bug" if bug else "AD75 原式  "
    print(f"{tag}: |Tc_AD/Tc_exp-1| 均值 {rel.mean()*100:6.2f}%  "
          f"std {rel.std()*100:5.2f}%  最大 {rel.max()*100:6.2f}% ({name[int(rel.argmax())]})")
print()
print("agent 报的是: 均值 13.2% +- 4.1%, 最大 23.4% (b-Ga)")
print()
t_bug = ad(lam, wlog, wb2, mus, True); t_ok = ad(lam, wlog, wb2, mus, False)
print(f"{'材料':<18}{'Tc_exp':>8}{'AD原式':>9}{'err%':>8}{'AD+bug':>9}{'err%':>8}")
for i in np.argsort(-np.abs(t_bug/Tc-1))[:8]:
    print(f"{name[i]:<18}{Tc[i]:8.2f}{t_ok[i]:9.2f}{(t_ok[i]/Tc[i]-1)*100:8.2f}"
          f"{t_bug[i]:9.2f}{(t_bug[i]/Tc[i]-1)*100:8.2f}")
