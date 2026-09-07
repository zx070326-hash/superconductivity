# 各向同性 Eliashberg 求解器 · 规格

目的：为 `theory/derivations/P1_*` 与 P2 提供数值验证器。不是生产代码，是可信、可测、小而透明的研究工具。

## 语言与依赖
Python 3.12，仅 numpy / scipy（可选 matplotlib 画图）。不依赖 QE/EPW。单文件模块 `eliashberg_iso.py` + `tests/`。

## 物理方程（虚轴 Matsubara，各向同性，常 DOS 版本）

ω_n = πT(2n+1)。
  Z_n = 1 + (πT/ω_n) Σ_m λ(n−m) · ω_m / √(ω_m² + Δ_m²)
  Z_n Δ_n = πT Σ_m [ λ(n−m) − μ*(ω_c) θ(ω_c − |ω_m|) ] · Δ_m / √(ω_m² + Δ_m²)
  λ(n−m) = 2 ∫₀^∞ dω ω α²F(ω) / (ω² + (ω_n − ω_m)²)

Tc：线性化后核矩阵最大本征值 = 1，对 T 二分。也提供 T < Tc 的非线性解（给 Δ(T)、能隙比 2Δ₀/k_BTc 用——P6 需要）。

## 必备功能
1. **任意 α²F(ω)**：表格输入（ω 网格 + 值），内部三次样条或梯形积分；提供 Einstein、Lorentzian、双峰等模型谱生成器。输出 λ、ω_log、ω₂、Allen–Dynes Tc（含 f₁f₂）供对照。
2. **μ*(ω_c) 一致重整化**：内部约定 μ*(ω_c) 与 Matsubara 截止 ω_c 绑定；提供换算 1/μ*(ω_c₁) − 1/μ*(ω_c₂) = ln(ω_c₂/ω_c₁)，以及从 (μ, E_F) 生成 μ*(ω_c)。所有输出必须注明用的是哪个约定。
3. **质量缩放**：`scale_mass(a2F, s)` → α²F_s(ω) = α²F(ω·s)，s = √(M/M₀)（保持 λ 不变）。**部分缩放**：α²F = Σ_i α²F_i，只缩放指定分量（H 投影谱）。
4. **同位素系数**：`isotope_coefficient(a2F_parts, which, mu_star_spec, ...)` 用有限差分 −Δ ln Tc / Δ ln M，μ* 按 spec 决定是否随 ω_c 一起变（这是 P1 §3.1 的核心开关：`mu_star_fixed=True/False`）。
5. **对数敏感度**：`log_sensitivity(param)` 数值算 ∂ ln Tc / ∂ ln(参数)，参数包括 λ（按比例缩放 α²F 幅度）、ω_log（缩放 ω 轴）、μ*、ω_c、以及下面 6、7 的模型参数。
6. **能量依赖 DOS**（第二阶段）：Mitrović–Carbotte / Pickett 推广（Z、Δ、χ 三方程，对 ε 积分，N(ε) 表格输入）。先实现台阶型与对数 vHs 型模型 N(ε)。
7. **模型顶角核钩子**（第二阶段，P2 用）：允许 λ(n−m) → λ(n−m)·[1 + V(n, m; 参数)]，V 由调用方传入可调用对象。

## 测试（必须通过才算 v0）
- Einstein 谱 λ=1, ω_E=50 meV, μ*=0.1：Tc 与 Allen–Dynes 公式相差在文献已知的偏差范围内（<10%），且 Matsubara 截止与 T 收敛性有测试
- 强耦合 Einstein 谱 λ=2、3：Tc/ω_E 单调、与 Allen–Dynes 渐近 √λ 行为一致
- μ*=0、谐波、全缩放：数值 α = 0.500 ± 0.002（这是 P1 §1 的基本检验）
- μ*≠0、`mu_star_fixed=False`：数值 α 与 P1 式 (3.2) 在 λ=1–2、μ*≤0.15 下相差 <0.01；`mu_star_fixed=True` 时回到 0.5
- 部分缩放：两分量谱，只缩放一个，α_1 + α_2 = 0.5（Rainer–Culetto 和规则）
- 若本机有 Julia + IsoME（Kogler 2025），加一个交叉对照测试；没有就跳过并注明

## 接口草案
```python
from eliashberg_iso import Spectrum, Solver
sp = Spectrum.from_table(w, a2F)          # 或 Spectrum.einstein(lam=2, w0=0.15)
solver = Solver(sp, mu_star=0.13, w_cut=10*sp.w_max, n_matsubara=None)  # None → 自动按 w_cut/T
Tc = solver.tc()                           # K
Tc_AD = sp.allen_dynes(mu_star=0.13)       # 对照
alpha = solver.isotope_coefficient(parts=[sp_H, sp_host], scale="H", mu_star_fixed=False)
s_lam = solver.log_sensitivity("lambda")
```

## 交付
- `theory/solvers/eliashberg_iso.py`
- `theory/solvers/tests/test_eliashberg_iso.py`（pytest）
- `theory/solvers/README.md`：约定、收敛参数建议、已知局限
- 一个示例脚本 `examples/p1_first_pass.py`：对 λ=2、ω_log=100 meV 的模型谱，扫 μ*=0.1–0.3，输出 α 表，与 P1 §3.1 表对照
