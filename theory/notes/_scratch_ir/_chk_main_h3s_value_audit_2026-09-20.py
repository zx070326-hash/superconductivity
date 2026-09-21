# -*- coding: utf-8 -*-
"""H3S 落点线的「还值不值得推」审计（主会话，2026-09-20，从零）。
只做三件事，全部是算术核对，不产生新物理：
 (1) 用 AD75 29 行回归系数正算 H3S/D3S 的 f_pred，复核 closure note §3 的 0.4736/0.5260；
 (2) 把它与望远镜预算的 f_pred（0.4792/0.5203，Eliashberg 全解）比 —— 两条路差多少；
 (3) 复核谐波行的 z 是否由 λ 主导（77% 之说）。
严格性：numerical-evidence（float64 算术，输入为两份笔记印出的数）。
★ 方向：本文不含单边界判决，全是双侧比较。
"""
import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
A, B, C = -1.0691, 0.9972, 0.7198          # closure note §3 的回归系数
S = 0.0794                                  # 无偏残差标准差 (n-p=26)
rows = [  # name, lam, shape, f_obs, f_pred_telescope(None=无)
    ("H3S-200 非谐", 1.7210, 0.7369, 0.44325, 0.4792),
    ("D3S-200 非谐", 1.7838, 0.8113, 0.48731, 0.5203),
    ("H3S-200 谐波", 2.2929, 0.8287, 0.44325, 0.6045),
    ("D3S-200 谐波", 2.2938, 0.8557, 0.48731, 0.6111),
]
print(f"{'体系':>14}{'λ':>8}{'shape':>8}{'f_pred(回归)':>13}{'f_obs':>9}"
      f"{'偏离%':>9}{'f_pred(望远镜)':>15}{'两路差%':>9}")
for nm, lam, sh, fo, ft in rows:
    fp = np.exp(A + B*np.log(lam) + C*np.log(sh))
    dev = 100*(fo-fp)/fp
    d2 = 100*(ft-fp)/fp
    print(f"{nm:>14}{lam:>8.4f}{sh:>8.4f}{fp:>13.4f}{fo:>9.4f}{dev:>8.2f}%"
          f"{ft:>15.4f}{d2:>8.2f}%")
print()
# (3) 谐波 vs 非谐的 f_pred 之比，拆成 λ 与 shape 两项
for tag, (l1,s1), (l2,s2) in (("H3S",(2.2929,0.8287),(1.7210,0.7369)),
                              ("D3S",(2.2938,0.8557),(1.7838,0.8113))):
    dl, ds = B*np.log(l1/l2), C*np.log(s1/s2)
    print(f"  {tag}: log(f_har/f_anh) = {dl+ds:+.5f}  "
          f"= λ 项 {dl:+.5f} ({100*dl/(dl+ds):.1f}%) + shape 项 {ds:+.5f} "
          f"({100*ds/(dl+ds):.1f}%)   比值 {np.exp(dl+ds):.4f}")

# (4) 「任何单调 Tc(λ) 闭式都会拒绝谐波行」—— 这句话我自己要核，不许直接写进笔记
def tc_ad(lam, wlog, w2, mus):
    den = lam - mus*(1.0+0.62*lam)
    if den <= 0: return 0.0
    L1 = 2.46*(1.0+3.8*mus); r = w2/wlog
    L2 = 1.82*(1.0+6.3*mus)*r
    f1 = (1.0+(lam/L1)**1.5)**(1.0/3.0)
    f2 = 1.0 + (r-1.0)*lam**2/(lam**2+L2**2)
    return f1*f2*wlog/1.20*np.exp(-1.04*(1.0+lam)/den)
def tc_mcm(lam, wlog, mus):
    den = lam - mus*(1.0+0.62*lam)
    return 0.0 if den <= 0 else wlog/1.45*np.exp(-1.04*(1.0+lam)/den)
print("\n  AD75 Eq.(34) / McMillan 直接算 H3S（μ*=0.13），对照实验 184.04 K：")
for tag, lam, wlog, w2 in (("非谐",1.7210,1276.44,1732.12), ("谐波",2.2929,1243.60,1500.64)):
    a, m = tc_ad(lam,wlog,w2,0.13), tc_mcm(lam,wlog,0.13)
    print(f"    {tag}: AD = {a:7.2f} K ({100*(a-184.04)/184.04:+6.1f}%)   "
          f"McMillan = {m:7.2f} K ({100*(m-184.04)/184.04:+6.1f}%)")
