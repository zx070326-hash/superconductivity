# -*- coding: utf-8 -*-
"""诊断：c(lam):=max_w-root(kTc/A) 在 Einstein 族上对 lam 的锯齿结构。

背景：`_chk_main_leavens_c_and_phiA_rescale_2026-09-23.py` 的 `maximize_c`/`refine_max`
（粗网格 + 黄金分割局部精修）在寻 `c(mus=0.13, r=10)` 的 sup 时，不同种子会收敛到不同的
局部齿顶，导致两次独立复核给出 0.0704% 的差（0.165741 vs 0.165624）。

根因：`Nc=floor(r*wbar)` 是 wbar 的阶梯函数 ⟹ `lam1_exact(wbar,lam,mus,r,N)` 在
`wbar=k/r`（k 整数）处有跳变 ⟹ 求根 `lam1(wbar*)=1` 给出的 `wbar*(lam)`（随 lam 单调降）
连续扫过这些阈值时，`c(lam)=1/(pi*lam*wbar*(lam))` 跟着出现一串小的下跳 = 锯齿。

本脚本独立复核这个机制：
  GATE S1 — 在两个相邻 lam（跨过 wbar=1.20, k=12 分支边界）上直接打印 wbar*、Nc、c，
            逐点确认跳变位置与方向。
  GATE S2 — 对 lam 在 [0.30, 3.00] 步长 0.005 的网格上扫 c(lam)，报告最高的若干齿峰，
            确认 k=12 分支（lam≈1.60 附近）在本脚本独立复核中同样是目前扫到的最高齿。

用途：仅作为「台阶机制被逐点抓到」的数值证据（numerical-evidence 级）。
不作为 c(mus,r) 全局 sup 的认证——team-lead 已用广义本征值分支枚举给出解析结构
（lam*(w)=1/theta_max(w)，每支光滑、齿顶在分支左端 w=k/r），本脚本与该结论方向一致，
不重复其结果。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import numpy as np

_d = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, name + ".py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


LV = _load("_chk_main_leavens_c_and_phiA_rescale_2026-09-23")

MUS, R = 0.13, 10.0


def gate_s1(N=200):
    print("GATE S1  跨过 wbar=1.20 (k=12) 分支边界的逐点跳变", flush=True)
    for lam in [1.5990, 1.5995, 1.6000, 1.6004, 1.6010, 1.6050, 1.6100]:
        ncr, wstar = LV.roots_in_wbar(
            lambda w: LV.lam1_exact(w, lam, MUS, R, N), 1.0, 0.05, 12.0)
        if wstar is None:
            print(f"  lam={lam:.4f}  无根")
            continue
        Nc = int(np.floor(R * wstar))
        c = 1.0 / (np.pi * lam * wstar)
        print(f"  lam={lam:.4f}  wbar*={wstar:.8f}  Nc={Nc:3d}  c={c:.8f}", flush=True)
    print()


def gate_s2(N=200, lo=0.30, hi=3.00, step=0.005, top=8):
    print(f"GATE S2  lam in [{lo},{hi}] 步长 {step}，N={N}，报告最高 {top} 个齿峰", flush=True)
    lams = np.arange(lo, hi + 1e-9, step)
    cs = np.full(len(lams), np.nan)
    for i, lam in enumerate(lams):
        c, ncr = LV.c_of_lambda(lam, MUS, R, N)
        cs[i] = c if c is not None else np.nan
    order = np.argsort(-np.nan_to_num(cs, nan=-1.0))
    for idx in order[:top]:
        print(f"  lam={lams[idx]:.4f}  c={cs[idx]:.8f}", flush=True)
    j = int(np.nanargmax(cs))
    print(f"  本扫描全局 max: lam={lams[j]:.4f}  c={cs[j]:.8f}")
    print("  （已知：此值随 N 单调下降，见 team-lead 的 N=200..800 收敛表；"
          "N=200 处的读数不是 sup，只是同一齿在此截断下的欠收敛值。）")
    print()


if __name__ == "__main__":
    gate_s1()
    gate_s2()
