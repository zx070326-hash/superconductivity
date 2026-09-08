"""
theory/derivations/p2_schrodi_fsr.py

Schrödi, Oppeneer, Aperis, "Full-bandwidth Eliashberg theory of
superconductivity beyond Migdal's approximation", Phys. Rev. B 102, 024503
(2020) = arXiv:1911.12872v2 (cond-mat.supr-con, "Dated: October 29, 2021"
on the v2 PDF header). PDF: library/pdf/D1/2020_Fabian_Schr_di_1911.12872.pdf.

本文件实现该文 Sec. II B "Fermi-surface restricted, isotropic approximation
to the non-adiabatic equations" 的 Eq. (19)-(24)（PDF 第 4-5 页），即：各向
同性、费米面限制（FSR）、含"第一顶角修正"（crossed diagram, Fig. 1(b)）的
Eliashberg 方程。与本仓库现有求解器 theory/solvers/eliashberg_iso.py 完全
独立——只 import 它的 Spectrum / Solver / KB_MEV_PER_K，不修改、不复用其
私有方法。

============================================================================
0. 逐字抄录（2026-09-08，对照 pdftotext -layout 与页面渲染图逐项复核两遍）
============================================================================

--- PDF 第 4 页，Sec. II B 末尾 ---

"We focus on the second term on the right-hand side of Eq. (4). Following
the same procedure as in deriving Eqs. (7-9) we obtain the following
expression for the non-adiabatic part of the self-energy:"

    Σ̂_m^(iso) = T² P Σ_{m',m''} λ_{m-m'} λ_{m'-m''}
                 × ρ̂₃ ĝ_{m'} ρ̂₃ ĝ_{m''} ρ̂₃ ĝ_{m''-m'+m} ρ̂₃                (19)

"where we use the label (iso) for 'isotropic'. The constant P is defined
as"

    P = (1/N₀³) Σ_{k,k',k''} δ(ξ_k) δ(ξ_k') δ(ξ_k'') δ(ξ_{k''-k'+k})       (20)

"and reflects momentum conservation of processes at the Fermi level."

--- PDF 第 5 页 ---

"The propagator ĝ_m is obtained by integrating the electron Green's
function over energy,"

    ĝ_{k,m} = ∫_{-∞}^{∞} dε Ĝ_{k,m}(ε)                                    (21)

"and the result is furthermore assumed to be isotropic, ĝ_{k,m} ≃ ĝ_m.
We get"

    ĝ_m = -π · (iω_m ρ̂₀ + Δ_m^(iso) ρ̂₁) / sqrt(ω_m² + (Δ_m^(iso))²)       (22)

"With the above considerations, the resulting Eliashberg equations are"

    Z_m^(iso) = 1 + (πT/ω_m) Σ_{m'} λ_{m-m'} γ_{m'}^(ω)
        + P (π³T²/ω_m) Σ_{m',m''} λ_{m-m'} λ_{m'-m''}
          × ( -γ_{m'}^(ω) γ_{m''}^(Δ) γ_{m''-m'+m}^(Δ)
              +γ_{m'}^(Δ) γ_{m''}^(ω) γ_{m''-m'+m}^(Δ)
              -γ_{m'}^(Δ) γ_{m''}^(Δ) γ_{m''-m'+m}^(ω)
              -γ_{m'}^(ω) γ_{m''}^(ω) γ_{m''-m'+m}^(ω) )                   (23)

    Δ_m^(iso) = (πT/Z_m^(iso)) Σ_{m'} λ_{m-m'} γ_{m'}^(Δ)
        + P (π³T²/Z_m^(iso)) Σ_{m',m''} λ_{m-m'} λ_{m'-m''}
          × ( -γ_{m'}^(Δ) γ_{m''}^(Δ) γ_{m''-m'+m}^(Δ)
              -γ_{m'}^(ω) γ_{m''}^(ω) γ_{m''-m'+m}^(Δ)
              +γ_{m'}^(ω) γ_{m''}^(Δ) γ_{m''-m'+m}^(ω)
              -γ_{m'}^(Δ) γ_{m''}^(ω) γ_{m''-m'+m}^(ω) )                   (24)

    [2026-09-08 更正：第 3、4 项的 Δ/ω 位置最初抄反了（曾误写成
    +γ_{m'}^(Δ)γ_{m''}^(ω)γ^(ω) 与 -γ_{m'}^(ω)γ_{m''}^(Δ)γ^(ω)）。经与 PDF 第 5
    页 6x 放大渲染图重新逐字核对（见下方复核记录），确认正确顺序是第 3 项
    (ω,Δ,ω) 系数 +1、第 4 项 (Δ,ω,ω) 系数 -1，现已按此改正上式与下面 §1b 的
    全部推导、kv_delta() 实现、CSV 结果。Eq (23) 未受影响（两次核对一致）。]

"Above we use γ_m^(ω) = ω_m/sqrt(ω_m²+(Δ_m^(iso))²) and
γ_m^(Δ) = Δ_m^(iso)/sqrt(ω_m²+(Δ_m^(iso))²)."

记号换算：λ_{m-m'} ≡ N₀V_{q,m-m'}（原文 Eq. 28-29，各向同性 FSR 耦合核，双动量
费米面平均），与本仓库 eliashberg_iso.py 的 λ(ν) 是同一个物理量：
λ_{m-m'} = spectrum.lambda_of_nu(2πT(m-m'))。ω_m = πT(2m+1) 是费米子 Matsubara
频率，m ∈ ℤ（原文未做正频折叠，Eq. 23-24 的 m,m',m'' 取遍全体整数——这与
eliashberg_iso.py 折叠前的"物理方程"一致，见该文件模块 docstring 第 14-18 行）。
Δ_m^(iso) 已经是 φ_m/Z_m（Eq. 24 两项都显式带 1/Z_m^(iso) 前因子），与
eliashberg_iso.Solver 里 "Z_nΔ_n = πTΣ[...]" 的 Δ_n 是同一个量——本文件的线性
化方程直接写成 "Z_m^(0)Δ_m = ..." 形式，与该求解器可以逐位对照（§2）。

我自己复核记录（2026-09-08 已按 coordinator 转达的 reviewer 高清 PDF 核对更
正）：Z 方程 (23) 括号内四项的 (γ_{m'},γ_{m''},γ_{m''-m'+m}) 类型序列是
(ω,Δ,Δ)→(-1)，(Δ,ω,Δ)→(+1)，(Δ,Δ,ω)→(-1)，(ω,ω,ω)→(-1)——这一项经两次独立
高清核对（我自己 + reviewer）均一致，未改动；Δ 方程 (24) 是 (Δ,Δ,Δ)→(-1)，
(ω,ω,Δ)→(-1)，(ω,Δ,ω)→(+1)，(Δ,ω,ω)→(-1)（此前误抄成 (Δ,ω,ω)→(+1)、
(ω,Δ,ω)→(-1)，即第 3、4 项的 (ω,Δ,ω)/(Δ,ω,ω) 类型与系数被对调了；已用 PDF
第 5 页 6x 放大渲染图重新逐字核对确认现在这个版本正确）。两个方程各自的四项
恰好对应"3 个 γ 里 Δ 的个数为偶数（0 或 2，共 4 种排列）"（Z 方程）与"奇数
（1 或 3，共 4 种排列）"（Δ 方程），与 C(3,0)+C(3,2)=4、C(3,1)+C(3,3)=4 吻合，
是一个内部一致性检查（8 种 (ω/Δ)³ 排列被两个方程各占一半，无遗漏无重复）——
注意这条一致性检查本身不能区分"(ω,Δ,ω)/(Δ,ω,ω) 谁配 +1 谁配 -1"，因为两种
分配都满足"偶数个 Δ→Z 方程，奇数个 Δ→Δ 方程"，所以它没能在第一次转录时
发现这个符号/类型错配——真正发现问题靠的是外部 reviewer 逐字比对高清图像。

============================================================================
1. 线性化推导（Δ→0，本文件的核心推导，供 reviewer 复核）
============================================================================

Δ→0 时 γ_m^(ω) → sign(ω_m) ≡ s_m（偶函数，O(Δ²) 修正，无 O(Δ) 项）；
γ_m^(Δ) → Δ_m/|ω_m|（奇函数，O(Δ) 领头阶，无 O(Δ⁰) 项）。

--- 1a. Z 方程：只保留 Δ=0 的项 ---

(23) 括号内四项里，含 γ^(Δ) 的三项在 Δ→0 时全部为零（每个 γ^(Δ) 至少贡献一个
O(Δ) 因子），只剩 (ω,ω,ω) 项存活：

    Z_m^(0) = 1 + (πT/ω_m) Σ_{m'} λ_{m-m'} s_{m'}
              - P (π³T²/ω_m) E(m)                                        (Z-lin)

    E(m) ≡ Σ_{m',m''} λ_{m-m'} λ_{m'-m''} s_{m'} s_{m''} s_{m''-m'+m}

（"-"号来自 (23) 里 (ω,ω,ω) 项本身的系数 -1。）第一项是标准 Migdal Z 方程；
第二项是顶角对 Z 通道的领头阶修正，与 Δ 无关（这是"Z 通道顶角修正"这个说法
的精确含义：它在 Tc 处就已经是非零的领头阶效应，不需要 Δ 展开）。

--- 1b. Δ 方程：保留 O(Δ¹) 的项，折成对 Δ_k 的线性核 ---

（本节 2026-09-08 已按更正后的 Eq (24) 重新推导；下面三项的类型/系数是修正
后的版本：(ω,ω,Δ)→-1 未变，(ω,Δ,ω)→+1、(Δ,ω,ω)→-1 与最初版本相比互换了。）

(24) 括号内四项里，(Δ,Δ,Δ) 是 O(Δ³)，线性化时丢弃；其余三项各含恰好一个
γ^(Δ)，是 O(Δ¹)：

    (ω,ω,Δ) 系数 -1：  -s_{m'} s_{m''} · Δ_{m''-m'+m}/|ω_{m''-m'+m}|
    (ω,Δ,ω) 系数 +1：  +s_{m'} · [Δ_{m''}/|ω_{m''}|] · s_{m''-m'+m}
    (Δ,ω,ω) 系数 -1：  -[Δ_{m'}/|ω_{m'}|] · s_{m''} · s_{m''-m'+m}

把这三项各自的"自由求和指标 + 被 Δ 携带的指标"做换元，统一折成对外部指标 m
与"Δ 所在指标" k 的双指标核。定义（与之前一致，未改动）

    C(a,b) ≡ Σ_j λ_{a-j} · s_j · s_{j+b-a}
    D(a,b) ≡ Σ_j λ_{a-j} · λ_{j-b} · s_j · s_{a+b-j}

逐项代数（三项分别记为 T2/T3/T4，对应上面 (ω,ω,Δ)/(ω,Δ,ω)/(Δ,ω,ω)）：

  T2 (ω,ω,Δ)，系数 -1：Δ 携带指标是 m''-m'+m。令 k=m''-m'+m，即
  m''=k-m+m'，代入得 m'-m''=m-k，故 λ_{m-m'}λ_{m'-m''}=λ_{m-m'}λ_{m-k}
  （λ_{m-k} 与求和变量 m' 无关，可提出）：
      T2(m) = Σ_k λ_{m-k} · [Σ_{m'} λ_{m-m'} s_{m'} s_{m'+(k-m)}] · Δ_k/|ω_k|
            = Σ_k λ_{m-k} · C(m,k) · Δ_k/|ω_k| · (-1)
  即对 K^V_Δ[m,k] 的贡献为 **-λ_{m-k}·C(m,k)**。

  T3 (ω,Δ,ω)，系数 +1：Δ 携带指标就是求和变量 m'' 本身，令 k=m''：
      T3(m) = Σ_k [Σ_{m'} λ_{m-m'}λ_{m'-k} s_{m'} s_{m+k-m'}] · Δ_k/|ω_k|
            = Σ_k D(m,k) · Δ_k/|ω_k|
  （Σ_{m'}λ_{m-m'}λ_{m'-k}s_{m'}s_{m+k-m'} 恰好是 D(a=m,b=k) 的定义式，
  j=m'。）即对 K^V_Δ[m,k] 的贡献为 **+D(m,k)**。

  T4 (Δ,ω,ω)，系数 -1：Δ 携带指标是 m' 本身，令 k=m'，
  λ_{m-m'}=λ_{m-k} 与求和变量 m'' 无关可提出：
      T4(m) = Σ_k λ_{m-k} · (-1) · [Σ_{m''} λ_{k-m''} s_{m''} s_{m''+(m-k)}] · Δ_k/|ω_k|
            = Σ_k λ_{m-k} · (-1) · C(k,m) · Δ_k/|ω_k|
  （Σ_{m''}λ_{k-m''}s_{m''}s_{m''+(m-k)} 恰好是 C(a=k,b=m) 的定义式。）
  即对 K^V_Δ[m,k] 的贡献为 **-λ_{m-k}·C(k,m)**。

三项相加，得到修正后的顶角线性核：

    K^V_Δ[m,k] = -λ_{m-k}·[C(m,k) + C(k,m)] + D(m,k)                      (KV-lin)

（对照最初版本的错误结果 λ_{m-k}·[C(k,m)-C(m,k)]-D(m,k)：由于 T3/T4 的类型
被对调，原推导把 D(m,k) 的符号和 C(m,k)/C(k,m) 的相对符号都搞反了——不是简单
的整体变号，两个公式的函数形式本身不同，数值上也确实给出不同结果，见下面
kv_delta() 与 tests 的更新。）

线性化后完整的 φ=Z·Δ 方程（"Z_mΔ_m=φ_m"，与 eliashberg_iso.py 的
"Z_nΔ_n=πTΣK2[n,m]Δ_m/ω_m" 同一形式）：

    Z_m^(0) Δ_m = πT Σ_k [λ_{m-k} - μ*(ω_c)θ(ω_c-|ω_k|)] Δ_k/|ω_k|
                  + P π³T² Σ_k K^V_Δ[m,k] · Δ_k/|ω_k|                    (Δ-lin)

（注意：标准电声项里 λ_{m-k} 前面**没有** sign(ω_k) 因子——γ_k^(Δ)→Δ_k/|ω_k|
本身就是"无符号"的，与 γ_k^(ω)→sign(ω_k)（Z 方程用，因为其分母是 ω_m 不是
|ω_m|）不同。首次实现时误将 s_k 乘进了这一项，用 eliashberg_iso.py 模块
docstring 第 16 行的物理方程 "Z_nΔ_n=πTΣ_m[λ(n-m)-μ*θ(ω_c-|ω_m|)]Δ_m/√(...)"
核对（该式电声项前也没有符号因子）后发现并改正。这与 Eq (24) 转录错误是两个
独立的 bug，分别修正。)

Tc：(Δ-lin) 除以 Z_m^(0) 是标准线性本征值问题 Δ_m=Σ_k M[m,k]Δ_k，
M[m,k] = {πT[λ_{m-k}s_k-μ*θ_k] + Pπ³T²K^V_Δ[m,k]} / (|ω_k| Z_m^(0)（第一个
下标 m 对应的值）)；对 T 用 brentq 找最大实本征值=1 处，与
eliashberg_iso.Solver._max_eigen_linearized 是同一套算法思想（矩阵在"未折
叠"的全体 Matsubara 指标 m∈[-N,N) 上直接构造，不做该文件的正负频折叠——因为
s_m=sign(ω_m) 本身破坏了折叠前对称性的显然性，直接在全轴上验证 P=0 时与折叠
版求解器逐位一致，比重新证明折叠在顶角项下仍成立更省时间、更不容易出错）。

对称性检查（供 tests 使用，也是这份推导本身的自洽性检查）：折叠假设
Δ_{-1-n}=Δ_n、Z_{-1-n}=Z_n 对应位置反射 i → M-1-i（因为 idx[M-1-i] =
(M-1-i)-N = -1-idx[i]）。手动代数验证（换元 m'=-1-a, m''=-1-b，用 λ 偶性
λ(x)=λ(-x) 与 s_{-1-x}=-s_x）：
    E(-1-m) = -E(m)                （E 本身是**反对称**的，不是对称的！）
    S(-1-m) = -S(m)，其中 S(m)=Σ_{m'}λ(m-m')s_{m'}（标准 Migdal 项的求和部分）
而 ω_{-1-m}=-ω_m 本身也反号，两个反号相消，两项各自除以 ω_m 后都变回对称：
    Z_m^(0) = 1+(πT/ω_m)S(m)-P(π³T²/ω_m)E(m) 在 m→-1-m 下不变（这才是折叠
    技巧 Z_{-1-n}=Z_n 成立所需要的不变量，E(m) 单独不是）。
K^V_Δ[m,k] 在同时反射 (m,k)→(M-1-m,M-1-k) 下直接不变（不需要额外的 1/ω 因子
抵消，因为它已经是 Δ_k/|ω_k| 前面的系数，|ω_k| 本身在反射下不变）——这条不变
性对修正后的 (KV-lin) 重新做了手动代数验证（换元 j=-1-j' 分别代入 C(a,b) 与
D(a,b) 的定义，用 λ 偶性与 s_{-1-x}=-x 的两次变号相消，可证 C(-1-a,-1-b)=
C(a,b) 与 D(-1-a,-1-b)=D(a,b) 对任意 a,b 都成立，从而 K^V_Δ 在联合反射下不变
这条结论在旧公式和新公式下都成立，不是巧合而是 C、D 各自的对称性——数值验证
见 tests/test_p2_schrodi_fsr.py（对 E(m) 检验反对称、对 Z_m^(0) 与 K^V_Δ
检验对称）。若验证失败说明线性化推导或实现有错。

============================================================================
2. 关于常数 P（Eq. 20）——诚实的局限说明
============================================================================

Eq. (20) 的 P 是原文特定紧束缚晶格模型（Eq. 25-27，含 t^(1)/t^(2) 近邻/次近邻
跳跃、维度参数 c_i、化学势 μ）在费米面上的四点相空间和，量纲 [P]=[能量]⁻¹
（由 Eq. 19 的 T²·P·(...) 必须匹配 Σ̂ 的能量量纲反推，(...)全是无量纲的 γ、λ
乘积）。原文正文没有给出这个量对任意"3D 常数 DOS"体系的闭式解析值——它依赖
具体能带（Eq. 25-27 的 t^(1)、t^(2)、c_i、μ），Appendix A 只说明了数值求和
的截断/延拓技巧，没有报告 P 的具体数值；正文唯一定量陈述是"1D 时 P 使方程
无法求解"（第 8 页，见下方 grep 记录），说明 P 在低维可以趋于 0/病态，但没给
3D 数值。

本文件因此**不**声称精确复现 Schrödi 3D 晶格的 P 值，而是：
    (i) 用量纲分析给出一个数量级占位值 P_dimensional_estimate(E_F) = 1/(4·E_F)
        （3D 自由电子费米面相空间的标准数量级，O(1) 前因子未定；仅用于把 P
        的扫描锚定在物理上合理的量级，不是 Eq. 20 的精确计算）；
    (ii) 把 P 作为自由扫描参数（{0, P_3D, 2·P_3D}），因为线性化后的 (Z-lin)/
        (KV-lin) 在 P 里严格线性（矩阵元素是 P 的一次函数），"Z 通道 vs Δ 通道
        谁主导、净效应升高还是压低 Tc"这个定性问题的答案不依赖 P 的具体数值
        （只要 P>0 且不太大以至于破坏线性化前提）——这正是本次数值实验要回答
        的问题，用相对扫描规避了对 P 绝对值的依赖。
E_F 取本仓库 P2_vertex_sign_criterion.md §5.2 对 H3S 200 GPa 反推的 E_F^eff≈
0.96 eV（DOS 缺损静态泛函独立给出，见该文档）作为量级锚点，P_3D≈1/(4×960meV)
≈2.6e-4 meV⁻¹。

============================================================================
3. 与 eliashberg_iso.py 的关系
============================================================================

本文件只 import Spectrum、Solver、KB_MEV_PER_K（公开 API），不 import、不调用
任何下划线开头的私有方法，不修改该文件。SchrodiFSRSolver 是一个独立实现，
在 P=0（或 vertex_in_Z=vertex_in_Delta=False）时应给出与 Solver.tc() 一致的
Tc（因为退化到同一组标准 Migdal 线性化方程，只是矩阵在全轴而非折叠轴上构造），
用作正确性的基准测试（tests §1）。
"""
from __future__ import annotations

import dataclasses
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Tuple

import numpy as np
from scipy import optimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "solvers"))
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# H3S 200 GPa 的 E_F^eff（meV），来自 P2_vertex_sign_criterion.md §5.2（DOS 缺损
# 静态泛函反推值 0.96 eV），仅用作 P 的数量级锚点，见模块 docstring §2。
E_F_EFF_H3S_MEV = 960.0


def p3d_dimensional_estimate(e_f_mev: float) -> float:
    """P 的量纲占位估计 P ~ 1/(4 E_F)（meV⁻¹），见模块 docstring §2。

    **不是** Eq. (20) 对任何具体晶格的精确计算——原文没有给出适用于任意 3D
    常数 DOS 体系的闭式 P，这里只保证量纲正确（[P]=能量⁻¹）与量级合理
    （自由电子费米面相空间，O(1) 前因子未定）。
    """
    if e_f_mev <= 0:
        raise ValueError("e_f_mev 必须 > 0")
    return 1.0 / (4.0 * e_f_mev)


# ==========================================================================
# 顶角核的数值构造（§1 推导的直接实现）
# ==========================================================================

def _toeplitz(lam_pos: np.ndarray, m: int) -> np.ndarray:
    """Λ[i,j] = λ(idx[i]-idx[j]) = lam_pos[|i-j|]（idx 是连续整数，值差=位置差）。"""
    idx = np.arange(m)
    return lam_pos[np.abs(idx[:, None] - idx[None, :])]


def _sign_array(n_pos: int) -> np.ndarray:
    """s[i] = sign(ω_{idx[i]})，idx[i]=i-n_pos，s[i]=+1 若 idx[i]>=0 否则 -1。"""
    idx = np.arange(-n_pos, n_pos)
    return np.where(idx >= 0, 1.0, -1.0)


def vertex_kernels(lam_pos: np.ndarray, n_pos: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """构造 Λ（标准 Migdal Toeplitz 核）、C(a,b)、D(a,b)、E(a)，见模块 docstring §1。

    位置索引 a,b,j ∈ [0, m)（m=2·n_pos）对应 Matsubara 指标 idx[i]=i-n_pos。
    关键简化：s(整数值) = sign(2·value+1) 对任意整数值都可以直接求值（不需要
    落在 [0,m) 位置窗口内），只有 λ 需要经由 lam_pos 表（λ 是连续函数的数值
    积分，s 是解析的符号函数）——这让 C、D 的换元求和不必担心指标越界截断。

        C(a,b) = Σ_j λ(idx[a]-idx[j])·s(idx[j])·s(idx[j]+idx[b]-idx[a])
        D(a,b) = Σ_j λ(idx[a]-idx[j])·λ(idx[j]-idx[b])·s(idx[j])·s(idx[a]+idx[b]-idx[j])
        E(a)   = Σ_{a'} λ(idx[a]-idx[a'])·s(idx[a'])·C(a',a)

    复杂度 O(m³)（对 a 的 Python 循环，内部向量化 O(m²)），m~100-200 时数值上
    很快（见 tests 里的计时注记）。
    """
    m = 2 * n_pos
    s = _sign_array(n_pos)
    lam = _toeplitz(lam_pos, m)  # Λ[j,b]，与 a 无关，只建一次
    j = np.arange(m)
    b = np.arange(m)

    C = np.empty((m, m))
    D = np.empty((m, m))
    for a in range(m):
        w_a = lam_pos[np.abs(a - j)] * s[j]                 # 长度 m：λ(idx[a]-idx[j])·s(idx[j])
        val_C = j[:, None] + b[None, :] - a - n_pos          # idx[j]+idx[b]-idx[a]，[j,b] 网格
        sign_C = np.where(val_C >= 0, 1.0, -1.0)
        C[a, :] = np.sum(w_a[:, None] * sign_C, axis=0)
        val_D = a + b[None, :] - j[:, None] - n_pos           # idx[a]+idx[b]-idx[j]，[j,b] 网格
        sign_D = np.where(val_D >= 0, 1.0, -1.0)
        D[a, :] = np.sum(w_a[:, None] * lam * sign_D, axis=0)

    E = np.sum(lam * s[None, :] * C.T, axis=1)  # E[a] = Σ_{a'} Λ[a,a']·s[a']·C[a',a]
    return lam, C, D, E


def kv_delta(lam: np.ndarray, C: np.ndarray, D: np.ndarray) -> np.ndarray:
    """K^V_Δ[m,k] = -λ_{m-k}·[C(m,k)+C(k,m)] + D(m,k)，见 (KV-lin)。

    2026-09-08 更正：此前（基于错误转录的 Eq 24 第 3、4 项）给出的公式是
    λ_{m-k}·[C(k,m)-C(m,k)]-D(m,k)，与现在这个公式在函数形式上不同（不只是
    整体变号），见模块 docstring §1b 的完整逐项重新推导。
    """
    return -lam * (C + C.T) + D


# ==========================================================================
# SchrodiFSRSolver：线性化 Tc 求解器
# ==========================================================================

@dataclass
class SchrodiFSRSolver:
    """Schrödi Eq. (23)-(24) 的线性化（Δ→0）Tc 求解器，全轴（未折叠）Matsubara
    指标。P=0（或两个 vertex_in_* 都 False）时退化为标准各向同性 FSR Migdal
    方程，应与 eliashberg_iso.Solver.tc() 逐位一致（<0.1 K，见 tests）。

    参数
    ----
    spectrum, mu_star, w_cut : 与 eliashberg_iso.Solver 同名参数同义。
    P : Eq. (20) 的相空间常数（meV⁻¹），见模块 docstring §2。
    vertex_in_Z, vertex_in_Delta : 顶角修正是否分别作用于 Z 通道 (Z-lin)、
        Δ 通道 (KV-lin)。两者独立开关，用于本次数值实验的核心对照
        （只开 Z / 只开 Δ / 两者都开）。
    oversample, margin : 与 Solver._n_for_T 同名参数同义、同公式（复刻公式，
        不 import 私有方法），决定给定 T 下的正频 Matsubara 数 n_pos。
    n_pos_cap : 顶角核 O(n_pos³) 计算量的上限（None 表示不封顶，走
        Solver._n_for_T 同款上限 20000／2；P!=0 时的生产扫描建议显式设置
        一个较小的值，如 60-100，见模块 docstring 与 __main__ 脚本）。
    """
    spectrum: Spectrum
    mu_star: float
    w_cut: float
    P: float = 0.0
    vertex_in_Z: bool = True
    vertex_in_Delta: bool = True
    oversample: float = 3.0
    margin: int = 30
    n_pos_cap: Optional[int] = None

    def _n_pos_for_T(self, t_mev: float) -> int:
        n = int(np.ceil(self.oversample * self.w_cut / (2 * np.pi * t_mev))) + self.margin
        n = max(n, 4)
        cap = self.n_pos_cap if self.n_pos_cap is not None else 10000
        return min(n, cap)

    def max_eigen(self, t_k: float) -> float:
        t_mev = KB_MEV_PER_K * t_k
        n_pos = self._n_pos_for_T(t_mev)
        m = 2 * n_pos
        idx = np.arange(-n_pos, n_pos)
        wn = np.pi * t_mev * (2 * idx + 1)
        abs_wn = np.abs(wn)
        s = np.where(idx >= 0, 1.0, -1.0)

        lam_pos = self.spectrum.lambda_of_nu(2 * np.pi * t_mev * np.arange(m))

        if self.P != 0.0 and (self.vertex_in_Z or self.vertex_in_Delta):
            lam, C, D, E = vertex_kernels(lam_pos, n_pos)
        else:
            lam = _toeplitz(lam_pos, m)
            E = None

        theta = np.clip((self.w_cut - abs_wn) / (2 * np.pi * t_mev) + 0.5, 0.0, 1.0)

        z0 = 1.0 + (np.pi * t_mev / wn) * (lam @ s)
        if self.P != 0.0 and self.vertex_in_Z:
            z0 = z0 - self.P * (np.pi ** 3 * t_mev ** 2 / wn) * E

        rhs = np.pi * t_mev * (lam - self.mu_star * theta[None, :])
        if self.P != 0.0 and self.vertex_in_Delta:
            rhs = rhs + self.P * (np.pi ** 3 * t_mev ** 2) * kv_delta(lam, C, D)

        mat = rhs / (abs_wn[None, :] * z0[:, None])
        eigvals = np.linalg.eigvals(mat)
        i = np.argmax(eigvals.real)
        return float(eigvals.real[i])

    def tc(self, t_bracket_k: Optional[Tuple[float, float]] = None, xtol: float = 1e-8) -> float:
        def f(ln_t_k: float) -> float:
            return self.max_eigen(np.exp(ln_t_k)) - 1.0

        if t_bracket_k is None:
            try:
                t_guess = self.spectrum.allen_dynes(self.mu_star)
                t_lo, t_hi = t_guess / 8.0, t_guess * 4.0
            except ValueError:
                wlog = self.spectrum.omega_log()
                t_lo, t_hi = 1e-2 * wlog / KB_MEV_PER_K, 0.5 * wlog / KB_MEV_PER_K
        else:
            t_lo, t_hi = t_bracket_k
        lo, hi = np.log(t_lo), np.log(t_hi)
        f_lo, f_hi = f(lo), f(hi)
        tries = 0
        while f_hi > 0 and tries < 30:
            hi += np.log(2.0)
            f_hi = f(hi)
            tries += 1
        tries = 0
        while f_lo < 0 and tries < 30:
            lo -= np.log(2.0)
            f_lo = f(lo)
            tries += 1
        if f_lo < 0 or f_hi > 0:
            raise RuntimeError(
                f"无法为 Tc 定位二分区间 (f_lo={f_lo:.4g}, f_hi={f_hi:.4g})"
            )
        ln_tc = optimize.brentq(f, lo, hi, xtol=xtol)
        return float(np.exp(ln_tc))


# ==========================================================================
# 数据加载（与 p1_mishra_a2f.py 的 load_epw_a2f 同一读法，独立复制一份，不 import
# 该文件——避免把本文件的正确性绑定到另一个 p1_*.py 脚本的顶层执行副作用上）
# ==========================================================================

def load_epw_a2f(path: Path) -> Spectrum:
    """EPW .a2f：第一行表头，之后 'w[meV] a2f int' 直到 'Integrated el-ph coupling' 行。"""
    w, a = [], []
    for line in path.read_text().splitlines()[1:]:
        parts = line.split()
        if len(parts) < 3:
            break
        try:
            x, y = float(parts[0]), float(parts[1])
        except ValueError:
            break
        w.append(x)
        a.append(max(y, 0.0))
    return Spectrum.from_table(np.array(w), np.array(a), name=path.stem)


DATA = ROOT / "data" / "external" / "mishra_h3s_d3s" / "Data-H3S-D3S-PH-EPW-Anharmonic"
MISHRA_FILES = {
    "H3S": DATA / "H3S-200GPa/EPW_anharmonic/H3S.a2f",
    "D3S": DATA / "D3S-200GPa/EPW_anharmonic/D3S.a2f",
}

W_CUT_MEV = 2500.0   # ω_c = 2.5 eV，与 EPW 一致（任务规格，用于 H3S/D3S 对照）
MU_STAR = 0.16
# Einstein 交叉检验谱的 w0=60 meV 远小于 H3S/D3S 的 ω_log(~110 meV)；固定 ω_c=2.5 eV
# 对它意味着 w_cut/T_mev 比值远大于 H3S 情形，导致自动 Matsubara 窗口（oversample*
# w_cut/2πT + margin）膨胀到数千，稠密本征值分解 O(n_pos^3) 不可行。ω_c 的物理作用
# 只是"远高于声子谱之上，μ*(ω_c) 重整化截止"，与 EPW 数据集绑定的具体数值（2.5 eV）
# 无关——Einstein 交叉检验单独用 ω_c=10·w0（标准收敛选择），只用于验证线性化/顶角核
# 实现本身（P=0 基准、小 P 线性度），不参与与 EPW δ_H/δ_D 目标值的对照。
W_CUT_EINSTEIN_RATIO = 10.0
CONFIGS = [
    ("Z_only", True, False),
    ("Delta_only", False, True),
    ("both", True, True),
]


def run_scan(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    spectra = {name: load_epw_a2f(p) for name, p in MISHRA_FILES.items()}
    einstein = Spectrum.einstein(lam=1.0, w0=60.0, name="einstein_test")
    einstein_D = einstein.scale_mass(np.sqrt(2.0))
    w_cut_einstein = W_CUT_EINSTEIN_RATIO * 60.0

    e_f_mev = E_F_EFF_H3S_MEV
    p3d = p3d_dimensional_estimate(e_f_mev)

    print(f"P_3D（量纲估计，E_F={e_f_mev:.0f} meV）= {p3d:.6g} meV^-1\n")

    systems = {
        "H3S": (spectra["H3S"], W_CUT_MEV),
        "D3S": (spectra["D3S"], W_CUT_MEV),
        "Einstein_H": (einstein, w_cut_einstein),
        "Einstein_D": (einstein_D, w_cut_einstein),
    }

    rows = []
    n_pos_cap = 200  # 见性能实测：H3S/D3S 需要 n_pos~100-120，200 留有余量且仍 <1s/次

    baselines = {}
    for name, (sp, wc) in systems.items():
        tc0 = SchrodiFSRSolver(sp, MU_STAR, wc, P=0.0).tc()
        baselines[name] = tc0
        print(f"{name}: Tc^0 (P=0) = {tc0:.3f} K   lambda={sp.lambda_total():.3f}  omega_log={sp.omega_log():.2f} meV  w_cut={wc:.0f} meV")
        rows.append(dict(system=name, P_label="0", P_value=0.0, config="none",
                          vertex_in_Z=False, vertex_in_Delta=False, Tc_K=tc0,
                          delta=0.0))

    print()
    for name, (sp, wc) in systems.items():
        tc0 = baselines[name]
        for p_label, p_val in [("P3D", p3d), ("2P3D", 2.0 * p3d)]:
            for cfg_name, vz, vd in CONFIGS:
                solver = SchrodiFSRSolver(
                    sp, MU_STAR, wc, P=p_val,
                    vertex_in_Z=vz, vertex_in_Delta=vd, n_pos_cap=n_pos_cap,
                )
                tc_v = solver.tc(t_bracket_k=(0.3 * tc0, 3.0 * tc0))
                delta = float(np.log(tc_v / tc0))
                rows.append(dict(system=name, P_label=p_label, P_value=p_val, config=cfg_name,
                                  vertex_in_Z=vz, vertex_in_Delta=vd, Tc_K=tc_v, delta=delta))
                print(f"{name:12s} P={p_label:5s} cfg={cfg_name:10s} Tc={tc_v:8.3f} K  delta={delta:+.4f}")

    import csv
    csv_path = out_dir / "schrodi_tc_scan.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n写出 {csv_path}")

    print("\n=== delta_D/delta_H（P3D/2P3D）对照 EPW FSR 目标 (H -0.21, D -0.15) ===")
    for p_label in ["P3D", "2P3D"]:
        for cfg_name, _, _ in CONFIGS:
            dh = next(r["delta"] for r in rows if r["system"] == "H3S" and r["P_label"] == p_label and r["config"] == cfg_name)
            dd = next(r["delta"] for r in rows if r["system"] == "D3S" and r["P_label"] == p_label and r["config"] == cfg_name)
            ratio = dd / dh if dh != 0 else float("nan")
            print(f"P={p_label:5s} cfg={cfg_name:10s} delta_H={dh:+.4f} delta_D={dd:+.4f} delta_D/delta_H={ratio:.3f}")

    # ------------------------------------------------------------------
    # 小 P（线性区）诊断：P3D 处顶角项系数 P·π³T² 与 λ 同量级（见模块 docstring
    # §2 对 P_3D 局限性的说明），不是受控小微扰，Z/Δ 通道谁主导的结论在 P3D 处
    # 可能已经偏离线性响应区。这里额外在 P→0 的严格线性区（P=1e-6, 1e-5 meV^-1，
    # 数值验证：Tc 对 P 的响应在此区间与 P 严格成正比，见 tests）单独提取通道
    # 权重与符号——这个结论不依赖 P_3D 绝对值是否准确，只要 P>0。
    # ------------------------------------------------------------------
    print("\n=== 线性区（P->0）通道分解，符号与相对权重不依赖 P_3D 绝对值 ===")
    lin_rows = []
    p_lin_values = [1e-6, 1e-5]
    for name, (sp, wc) in systems.items():
        tc0 = baselines[name]
        for p_val in p_lin_values:
            for cfg_name, vz, vd in CONFIGS:
                solver = SchrodiFSRSolver(
                    sp, MU_STAR, wc, P=p_val,
                    vertex_in_Z=vz, vertex_in_Delta=vd, n_pos_cap=n_pos_cap,
                )
                tc_v = solver.tc(t_bracket_k=(0.7 * tc0, 1.5 * tc0))
                delta = float(np.log(tc_v / tc0))
                lin_rows.append(dict(system=name, P_value=p_val, config=cfg_name,
                                      vertex_in_Z=vz, vertex_in_Delta=vd, Tc_K=tc_v, delta=delta,
                                      delta_over_P=delta / p_val))
                print(f"{name:12s} P={p_val:.0e} cfg={cfg_name:10s} delta={delta:+.6f}  delta/P={delta/p_val:+.4g}")

    csv_lin_path = out_dir / "schrodi_linear_regime.csv"
    with open(csv_lin_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(lin_rows[0].keys()))
        writer.writeheader()
        writer.writerows(lin_rows)
    print(f"写出 {csv_lin_path}")


if __name__ == "__main__":
    run_scan(ROOT / "derivations" / "p2_out")
