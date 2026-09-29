# -*- coding: utf-8 -*-
"""
主会话 2026-09-25：接 `_chk_main_shape_factor_mustar_predictor_2026-09-25.py` 的自由度警报。

警报：AD75 Table I 的 μ* 是 McMillan–Rowell 反演**拟合**出来的（对齐超导态数据 Δ₀，
而 Δ₀ 与 Tc 强相关）⟹ 印刷 μ* 携带 Tc 信息 ⟹ 「加 μ* 后 CV 5.69%」可能是 Tc 借 μ* 回流。
诚实的预测器必须用**不看超导态数据**的 μ*。三档对照（同一 29 行、同一输入信息类）：

  (a) 印刷 μ*（样本内上限，信息回流）
  (b) 常数 μ* = 0.10、ω_c = ω_ph（发现场景下唯一可用的做法）
  (c) 常数 μ* = 0.13
对每档同时算 AD75 公式（f₂ 无开方，仓库已修）——**同信息类**：AD 用 (λ, ω_log, ω̄₂, μ*)，
我们用 (λ, ω̄₁, ω̄₂, μ*)；两者都吃 3 个谱量 + μ*。

报：Tc_pred = c·Tc_B，c 由留一（LOO）拟合 ⟹ 样本外 |误差|；AD 同样允许一个 LOO 标度 c_AD
（否则对 AD 不公平——它有自己的系统偏置）。另报 AD 裸公式（c_AD=1）。
第四因子试探：r_B 对 ρ_log := ω_log/ω̄₂ 的线性回归（LOO），问 ω_log 作为**预测因子**
（F1 已证它进不了**上界**）能否再吃掉残差。

严格性：全部 numerical-evidence。拟合常数 c 不是定理；Tc_B 的上界性质不依赖 c。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import math
import numpy as np
from multiprocessing import Pool

_d = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("P", os.path.join(_d, "_chk_main_shape_factor_mustar_predictor_2026-09-25.py"))
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
T1 = P.T1


def tc_ad(lam, wlog, w2, mus):
    L1 = 2.46*(1+3.8*mus); L2 = 1.82*(1+6.3*mus)*(w2/wlog)
    f1 = (1+(lam/L1)**1.5)**(1/3)
    f2 = 1 + (w2/wlog - 1)*lam*lam/(lam*lam + L2*L2)
    return f1*f2*wlog/1.2*math.exp(-1.04*(1+lam)/(lam - mus*(1+0.62*lam)))


def job(args):
    row, mus = args
    m, wl, w1, w2, lam, wph, mus_pr, tc, eta = row
    mu = mus_pr if mus is None else mus
    tb, w = P.tB(lam, min(w1/w2, 1.0), mu, wph/w2)
    return m, tc/(tb*w2/(2*math.pi)), tc/tc_ad(lam, wl, w2, mu), wl/w2


def loo_scale(r):
    """Tc_pred_i = mean_{j≠i}(r_j) · base_i ⟹ 误差 = Tc_pred/Tc − 1 = mean_{j≠i}/r_i − 1
    （与仓库 AD 门同口径：|Tc_pred/Tc_exp − 1|）。"""
    r = np.asarray(r); n = len(r)
    return np.array([((r.sum()-r[i])/(n-1))/r[i] - 1 for i in range(n)])


def loo_linear(r, x):
    r = np.asarray(r); x = np.asarray(x); e = []
    for i in range(len(r)):
        mk = np.arange(len(r)) != i
        b = np.polyfit(x[mk], r[mk], 1)
        e.append(np.polyval(b, x[i])/r[i] - 1)
    return np.array(e)


def main():
    arms = [("(a) 印刷 μ*", None), ("(b) μ*=0.10", 0.10), ("(c) μ*=0.13", 0.13)]
    with Pool(18) as Pl:
        out = {name: Pl.map(job, [(row, mu) for row in T1]) for name, mu in arms}
    print(f"{'档':14s} {'CV r_B':>8s} {'CV r_AD':>8s} | {'B·LOO 均/最大':>16s} | {'AD·LOO 均/最大':>16s} | "
          f"{'AD 裸 均/最大':>15s} | {'B+ρ_log LOO':>14s}")
    for name, _ in arms:
        rb = np.array([x[1] for x in out[name]]); ra = np.array([x[2] for x in out[name]])
        rl = np.array([x[3] for x in out[name]])
        eb = np.abs(loo_scale(rb)); ea = np.abs(loo_scale(ra)); e0 = np.abs(1/ra - 1)
        el = np.abs(loo_linear(rb, rl))
        print(f"{name:14s} {rb.std(ddof=1)/rb.mean()*100:7.2f}% {ra.std(ddof=1)/ra.mean()*100:7.2f}% | "
              f"{eb.mean()*100:6.2f}% / {eb.max()*100:6.2f}% | {ea.mean()*100:6.2f}% / {ea.max()*100:6.2f}% | "
              f"{e0.mean()*100:5.2f}% / {e0.max()*100:6.2f}% | {el.mean()*100:5.2f}%/{el.max()*100:5.2f}%")
        print(f"{'':14s} mean r_B = {rb.mean():.4f}  max r_B = {rb.max():.4f}（{'**>1 越界**' if rb.max() > 1 else '≤1'}）"
              f"  corr(r_B, ρ_log) = {np.corrcoef(rb, rl)[0,1]:+.3f}")
    print("\n读法：AD 裸公式在 (a) 档应复现仓库已核的 3.85% / 12.62%（转录 + f₂ 实现门）。")


if __name__ == "__main__":
    main()
