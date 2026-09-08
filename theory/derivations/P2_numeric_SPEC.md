# P2 数值任务规格：一阶电声顶角函数 P(q, ω_q; N(ε)) 的数值积分

来源公式：Grimaldi, Pietronero, Scattoni, Eur. Phys. J. B 10, 247 (1999), arXiv:cond-mat/9806012，**Eq. (6)**（本地 PDF `library/pdf/P2/1999_C_Grimaldi_cond-mat-9806012.pdf` 第 3 页，请用 Read 工具打开核对）。
目的：为 P2（`theory/derivations/P2_vertex_sign_criterion.md`）提供理论侧的顶角函数，检验三件事：(1) 在 H3S 参数区顶角平均 P̄ 是否为负；(2) P̄ 对声子频率的标度是否 ∝ ω₀（即同位素比 P̄_D/P̄_H ≈ 1/√2）；(3) H3S 的能量依赖 DOS 是否把 |P̄| 削弱约 2/3。这三个数已由第一性原理数据给出（P̄_H ≈ −0.095，比值 0.66–0.72，DOS 削弱 ~2/3），任务是看领头阶理论能否复现。

## 1. 公式

一阶顶角修正（自由电子传播子，Matsubara 求和已做）：

P(k+q, k; iω_k+iω_q, iω_k) = (ω₀/2) Σ_{k'} g²_{k−k'} · 1/(ε_{k'} − ε_{k'+q} + iω_q) ×
  [ (f(ε_{k'}) + n(−ω₀)) / (ε_{k'} + ω₀ − iω_k)
  − (f(ε_{k'}) + n(ω₀))  / (ε_{k'} − ω₀ − iω_k)
  − (f(ε_{k'+q}) + n(−ω₀)) / (ε_{k'+q} + ω₀ − i(ω_k+ω_q))
  + (f(ε_{k'+q}) + n(ω₀))  / (ε_{k'+q} − ω₀ − i(ω_k+ω_q)) ]

f 费米函数，n 玻色函数，n(−ω₀) = −1 − n(ω₀)。T→0：f = θ(−ε)，n(ω₀) = 0，n(−ω₀) = −1。取外线频率 ω_k = 0（或 πT，T 很小）。ω₀ 为 Einstein 声子频率。

**必须复现的解析极限**（原文 Eq. (11)、(19)、(20) 及正文；常数 DOS、结构无关耦合、半填充 E = 2E_F）：
- 静态极限（先 ω_q → 0 再 q → 0）：P_s = −ω₀/(ω₀ + E_F)
- 动态极限（先 q → 0 再 ω_q → 0）：P_d = +E_F/(ω₀ + E_F)
- 一电子极限 P_pol(q=0, ω_q=0) = λ·(E/2)/(ω₀ + E)（Eq. (11)；λ 的定义见 §2）
- 定性：固定 Q，P 随 ω_q 变化，**v_F q/ω_q < 1 时为正（动态区）、> 1 时为负（静态区）**——与 P_d > 0、P_s < 0 自洽；原文第 3 页写反了（第 1 页引言是对的），implementer 已用解析极限交叉验证（本行初版按原文第 3 页抄错，已勘误）

## 2. 各向同性模型（把 k' 求和变成能量-角度积分）

- 电子：各向同性带，能量 ε 从 −E_F（带底）到 E − E_F（带顶），DOS N(ε)；k' 求和 → ∫dε' N(ε') ∫₋₁¹ d(cosθ)/2，其中 θ 是 k' 与 q 的夹角
- ε_{k'+q}：线性化 ε_{k'+q} = ε' + v_F q cosθ（近费米面；q 用 Q ≡ q/(2k_F) 无量纲化，**v_F q = 4E_F·Q**——自由电子 v_F = 2E_F/k_F 与 q = 2k_F Q 的乘积；本规格初版误写 2E_F·Q，implementer 已按原文修正）。可选更精确的 ε_{k'+q} = ε' + ħ²(2k'q cosθ + q²)/2m（自由电子），两者在 Q ≲ 0.5 应一致，做一次对照
- 耦合 g²_{k−k'}：结构无关常数 g²，用 λ = 2N(0)g²/ω₀ 归一化（原文 g² = 2γ²/ω₀；请自己核对 λ 的定义使 Eq. (11) 成立），最终 P 以 λ 为单位输出（即输出 P/λ）
- 常数 DOS：N(ε) = N(0) 在 [−E_F, E_F]，半填充
- H3S DOS：`theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/H3S-200GPa/band/H3S.dos`（列：E(eV)、DOS(states/eV)、积分 DOS；E_F = 17.619 eV；带底约 −7.6 eV，即 E_F − 带底 ≈ 25 eV），N̄(ε) = N(ε)/N(E_F)，ε 相对 E_F；k_F、v_F 用自由电子等效值（n = 积分 DOS 到 E_F ≈ 9 电子/胞，胞体积从 `scf.in` 读；或简单取 E_F^free 使 N_free(E_F) = N(E_F)——写清你的选择）

## 3. 要算的量

A. **常数 DOS 验证**：P/λ 在 (Q, ω_q) 网格上，Q ∈ [0.01, 1.5]、ω_q ∈ [0, 5ω₀]，ω₀/E_F ∈ {0.01, 0.1, 0.5}；表格 + 图（matplotlib，存 `theory/derivations/p2_out/`）；核对 §1 的四个极限。
B. **Q 平均的顶角核** P_V(ω_q; Q_c) ≡ (1/Q_c)∫₀^{Q_c} dQ · P(Q, ω_q)（PGS 的"耦合在 Q < Q_c 内均匀"模型；若你认为该用 Q dQ 权重请说明并两种都给），Q_c ∈ {0.1, 0.3, 0.5, 1.0}，ω_q 网格同上。这是求解器 `vertex` 钩子的输入形式：λ(n−m) → λ[1 + λ P_V(ω_n − ω_m; Q_c)]。
C. **H3S 参数**：ω₀ = 110 meV（Mishra 非谐 ω_log），E_F 按 §2 的等效定义（写清数值），常数 DOS 与 H3S N̄(ε) 两种；输出 P_V(ω_q; Q_c) 表，以及**一个标量** P̄ ≡ P_V 在 ω_q ∈ [0, 3ω₀] 上按 Einstein 声子传播子 D(ω_q) = ω₀²/(ω_q² + ω₀²) 加权的平均（这近似对应进入 Eliashberg 核的有效值）。
D. **同位素**：ω₀ → 110/√2 meV（D3S），其余不变，重复 C；报告 P̄_D/P̄_H 对每个 Q_c。
E. **DOS 效应**：P̄(H3S DOS)/P̄(常数 DOS) 对每个 Q_c。

## 4. 交付

- `theory/derivations/p2_vertex_numeric.py`：函数 `vertex_P(Q, wq, w0, EF, dos: callable, T=0, wk=0)`、`vertex_PV(wq, Qc, ...)`、命令行入口跑 A–E 并把表/图写到 `theory/derivations/p2_out/`
- 单元测试 `theory/derivations/tests/test_p2_vertex.py`：四个解析极限（相对误差 <2%）
- 在文件顶部 docstring 写清所有约定（Σ_q 的符号、g² 与 λ 的关系、ε_{k'+q} 的近似、积分范围与节点数、T 与 ω_k 的取值）
- 报告：A 的极限核对表；C–E 的数字；你对"领头阶理论能否复现 P̄_H ≈ −0.095、比值 0.707、DOS 削弱 2/3"的判断；任何与原文公式不一致之处
- 提交（`-c user.name="p2-vertex" -c user.email="p2-vertex@local"`），**只 add 自己的文件**（显式路径），仓库多 agent 共享
- 不要改 `theory/solvers/` 与 `theory/derivations/P*.md`

## 5. 已知的物理张力（供你判断结果时参考）

Schrödi, Oppeneer, Aperis PRB 102, 024503 (2020) 的自洽全带宽解发现顶角函数**不是**小 q 峰化的，且非绝热参数越大越正、越各向同性——与 PGS 一阶自由传播子图像可能只在弱非绝热区一致。本任务只做一阶自由传播子（PGS 的层级），结果与第一性原理数据不符不算失败，是要记录的信息。
