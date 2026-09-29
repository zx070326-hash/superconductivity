# -*- coding: utf-8 -*-
"""
主会话 2026-09-25：三约束天花板 Tc_B(λ, ω̄₁, ω̄₂; μ*) 在 Mishra H₃S/D₃S 六条 EPW α²F 上
（`theory/data/external/mishra_h3s_d3s/`，读法沿用 `_chk_crossover_allen_dynes_moments_2026-09-21.py`）。

目的：
  (i) 真实**全谱**的 Tc_true（直接解线性化 Eliashberg，同一 μ*、同一 ω_c）必须 ≤ Tc_B（一致性门：
      全谱 ν 属于 𝒞̄_ρ 且 λ,A,M2 恰等 ⟹ 若 Tc_true > Tc_B，单原子归约为假）——**这是 μ*>0 下
      单原子归约在非原子测度上的第一个检验**（【R′】只测了 2/3 原子）。
  (ii) 比较 H₃S 上 Tc_B vs 旧的 min(Leavens, λ-分辨 Einstein)（记忆：非谐 306.0 vs 327.4 K）。
  (iii) Tc_true / Tc_B = 该材料离「给定三矩后的最优谱」多远——这是**谱形工程**的余量。

μ* = 0.13（仓库 H₃S 约定之一），ω_c = 10·ω̄₂（r=10，与认证曲线同约定）；另报 μ*=0。
严格性：numerical-evidence（float64；全谱由梯形求积离散 ⟹ 本身即一个多原子测度，故 (i) 是精确的
同一对象比较，不受求积误差影响——求积只改变「是哪个测度」，不改变「它属于 𝒞̄_ρ」）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import glob
import importlib.util
import math
import numpy as np
from multiprocessing import Pool

_d = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("P", os.path.join(_d, "_chk_main_shape_factor_mustar_predictor_2026-09-25.py"))
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
ROOT = os.path.join(_d, "..", "..", "data", "external", "mishra_h3s_d3s", "Data-H3S-D3S-PH-EPW-Anharmonic")
KB = 0.0861733262  # meV/K
C_INF = 0.1827262477
LEAV = {0.0: 0.230868, 0.13: 0.165741}   # c(μ*)：本仓库从零解（记忆 leavens-1975），μ*=0.13 仅 r=10 约定


def load(path):
    w, a = [], []
    for line in open(path):
        s = line.split()
        if len(s) < 2: continue
        try: x, y = float(s[0]), float(s[1])
        except ValueError: continue
        if x > 0: w.append(x); a.append(max(y, 0.0))
    return np.array(w), np.array(a)


def job(p):
    w, a = load(p)
    dw = np.diff(w); wt = np.concatenate([[dw[0]/2], (dw[:-1]+dw[1:])/2, [dw[-1]/2]])   # 梯形权
    nu = a*wt                                   # 离散测度（meV）
    lam = float((2*nu/w).sum()); A = float(nu.sum()); M1 = float((nu*w).sum())
    w1 = 2*A/lam; w2 = math.sqrt(2*M1/lam); rho = w1/w2
    # 离散测度的 λ 份额与频率（ω̄₂≡1）
    af = (2*nu/w)/lam; wf = w/w2
    keep = af > 1e-12; af, wf = af[keep], wf[keep]
    out = [os.path.relpath(p, ROOT).replace(os.sep, "/")[:30], lam, w1, w2, rho]
    for mus in (0.0, 0.13):
        wc = 10.0
        tt = P.tsup(lam, af, wf, mus, wc)
        tb, wst = P.tB(lam, rho, mus, wc)
        tE = P.tsup(lam, [1.], [1.], mus, wc)
        K = w2/(2*math.pi)/KB
        leav = LEAV[mus]*A/KB
        out += [tt*K, tb*K, tE*K, leav]
    return out


def main():
    paths = sorted(glob.glob(os.path.join(ROOT, "*", "EPW_*", "*.a2f")))
    with Pool(12) as Pl: res = Pl.map(job, paths)
    for mus_i, mus in enumerate((0.0, 0.13)):
        print("="*118); print(f"μ*={mus}, ω_c=10ω̄₂   [K]"); print("="*118)
        print(f"{'spectrum':30s} {'λ':>6s} {'ρ₁':>6s} {'Tc_true':>8s} {'Tc_B':>8s} {'Tc_E(λ)':>8s} {'Leavens':>8s} "
              f"{'min(旧)':>8s} {'B/min':>7s} {'true/B':>7s} {'门':>4s}")
        for r in res:
            tag, lam, w1, w2, rho = r[:5]
            tt, tb, tE, lv = r[5+4*mus_i: 9+4*mus_i]
            mn = min(tE, lv) if mus == 0.0 else tE   # μ*>0 时 Leavens c(0.13) 是 r=10 约定的 Einstein 数，同类
            if mus > 0: mn = min(tE, lv)
            gate = "OK" if tt <= tb*(1+1e-9) else "**X**"
            print(f"{tag:30s} {lam:6.3f} {rho:6.3f} {tt:8.1f} {tb:8.1f} {tE:8.1f} {lv:8.1f} {mn:8.1f} "
                  f"{tb/mn:7.4f} {tt/tb:7.4f} {gate:>4s}")
    print("\n注：μ*=0.13 行的 Leavens 列用 c(0.13, r=10)=0.165741（ω_c 约定为 10ω_E，与此处 10ω̄₂ 在 Einstein 上一致，"
          "在展开谱上**不同**⟹ 该列只作参照）。")


if __name__ == "__main__":
    main()
