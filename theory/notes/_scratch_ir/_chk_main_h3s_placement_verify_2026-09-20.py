# -*- coding: utf-8 -*-
"""主会话从零复核 h3s_on_ad75_trend_2026-09-20.md（2026-09-20）。
承重问题不是算术（算术我手算已对上），而是：判决对两处仓库已标记为有争议的输入是否稳健。
  (a) Tc_exp(H3S)=203K 在 155GPa，而 lambda/wbar2 在 200GPa。
      主会话 09-18 已矢量数字化 Errea Fig.2c：200GPa annealed = 184.04 K。
  (b) w_log(H3S): 仓库 alpha2F 给 109.995，Errea PRL2015 Table I 给 92.86（差 18%）。
  (c) D3S 的 w_log=77.778 由 1/sqrt(2) 精确伸缩推出，但 Errea 同时报了 H3S/D3S 两个
      w_log：73.3/92.86 = 0.7894，**不是** 0.7071。
严格性：numerical-evidence。
"""
import numpy as np

C_INF, MEV_K = 0.1827262477, 11.604518
B0, B1, B2 = -1.0691, 0.9972, 0.7198   # AD75 29-row regression
SCAT = 0.0794                           # residual scatter

def ceiling_K(lam, w2_mev): return C_INF*np.sqrt(lam)*w2_mev*MEV_K
def f_pred(lam, shape):     return np.exp(B0 + B1*np.log(lam) + B2*np.log(shape))

print("=== reproduce agent's two rows ===")
for name, lam, wlog, w2, tc in (("H3S-200anh",1.721,109.995,149.26,203.0),
                                ("D3S-200anh",1.784, 77.778,105.43,150.0)):
    ce = ceiling_K(lam,w2); fo = tc/ce; fp = f_pred(lam,wlog/w2)
    print(f"  {name}: ceiling={ce:7.2f}K  f_obs={fo:.4f}  f_pred={fp:.4f}  resid={(fo/fp-1)*100:+6.2f}%")

print("\n=== (c) Errea Table I reports BOTH w_log: is 1/sqrt(2) scaling supported? ===")
print(f"  Errea anh: w_log(D3S)/w_log(H3S) = 73.3/92.86 = {73.3/92.86:.4f}")
print(f"  Errea harm: 68.5/107.165(repo) n/a ; Errea harm pair 68.5/90.4 = {68.5/90.4:.4f}")
print(f"  exact uniform scaling 1/sqrt(2)  = {1/np.sqrt(2):.4f}")
print(f"  repo wbar2 ratio 105.43/149.26   = {105.43/149.26:.4f}   <- this one DOES scale")
print(f"  ==> w_log deviates {(73.3/92.86)/(1/np.sqrt(2))-1:+.1%} from uniform scaling; wbar2 only {(105.43/149.26)/(1/np.sqrt(2))-1:+.2%}")

print("\n=== (a)+(b) H3S verdict sensitivity: 2 contested inputs, 4 corners ===")
print(f"  band = +-{SCAT*100:.2f}%")
rows=[]
for tcl,tc in (("Tc=203 (155GPa, mismatched)",203.0),("Tc=184.04 (200GPa annealed, Fig2c)",184.04)):
    for wl,wlog in (("w_log=109.995 (repo a2F)",109.995),("w_log=92.86 (Errea TabI)",92.86)):
        ce=ceiling_K(1.721,149.26); fo=tc/ce; fp=f_pred(1.721,wlog/149.26); r=fo/fp-1
        ok = "PASS" if abs(r)<SCAT else "**FAIL**"
        rows.append(abs(r)<SCAT)
        print(f"  {tcl:36s} {wl:26s} f_obs={fo:.4f} f_pred={fp:.4f} resid={r*100:+7.2f}%  {ok}")
print(f"  ==> {sum(rows)}/4 corners pass  ==> verdict is {'ROBUST' if all(rows) else 'UNDETERMINED'}")

print("\n=== (c) D3S verdict sensitivity: derived vs Errea-ratio w_log ===")
for lbl,wlog in (("1/sqrt2 scaling (agent)",0.707107*109.995),("Errea ratio 0.7894",(73.3/92.86)*109.995)):
    ce=ceiling_K(1.784,105.43); fo=150.0/ce; fp=f_pred(1.784,wlog/105.43); r=fo/fp-1
    print(f"  {lbl:26s} w_log={wlog:6.2f} shape={wlog/105.43:.4f} resid={r*100:+7.2f}%  {'PASS' if abs(r)<SCAT else '**FAIL**'}")

print("\n=== D3S independence (my analytic invariance) ===")
print(f"  f_obs(H3S)={203.0/ceiling_K(1.721,149.26):.4f}   f_obs(D3S)={150.0/ceiling_K(1.784,105.43):.4f}"
      f"   differ {(150.0/ceiling_K(1.784,105.43))/(203.0/ceiling_K(1.721,149.26))-1:+.2%}")
print(f"  lambda 1.721 vs 1.784 differ {1.784/1.721-1:+.2%}  (exact scaling requires 0%)")
print("  ==> under exact harmonic uniform scaling f is INVARIANT ==> D3S is not an independent point;")
print("      its whole content is the SSCHA deviation from uniform scaling (few %).")
