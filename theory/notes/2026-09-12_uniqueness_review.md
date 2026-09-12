# P4 · Conjecture 2 唯一性半句（T↦𝔨(P,T) 严格递减：引理 A、B + 锥不变性主论证）· 对抗核验报告
日期：2026 年 9 月 12 日。审对象：分支 `p4-conj2` commit 5f05f83 —— `theory/derivations/P4_uniqueness_strategy.md`（§2.A 引理 A、B 证明，§3.A 主论证）、`theory/derivations/p4_uniqueness_lemmas.py`（37 检查）；上游 `P4_conjecture2.md` §1/§8、`P4_einstein_extremal.md` §1.C。立场：假设定理是错的，逐层找反例。
真值侧全部用**自写**代码（`_scratch_ir/rev8_lemmas.py`、`rev8_main.py`，按逐元定义 𝒪[f]_{nm}=−δ_{nm}2u_n²Σ_{k≤n}f(k)+(1−δ_{nm})f(|n−m|)u_nu_m+f(n+m+1)u_nu_m 从零构造，不 import 任何 derivations 模块）。`git diff 916296b 5f05f83 --stat`：只有策略文件（+51/−4）与新脚本（+331）两个文件变动。

裁决总表：

| 靶 | 裁决 | 一句话 |
|---|---|---|
| 1 引理 A（§2.A） | **CONFIRMED** | 三组分解逐指标重推无误（L≤L′ 下对角 n<L、m=n−k 支 n−k<L≤L′ 自动、L3 支 n∈[max(0,k−L′),min(L,k))）；自写暴力 k,L,L′≤300：min Q=0，Q=0 ⟺ L+L′≤k **逐格精确**，L+L′>k 时 min 3.33e-3；三组闭式 vs 矩阵部分和 max diff 5.3e-15；Q_{LL}(δ₁) 最小值恰为 1 |
| 2 引理 B（§2.A） | **CONFIRMED** | g(n)、R(n)、情形 I/II.a/II.b.1/II.b.2 逐行重推无误、情形穷尽；II.b.1 凸性恒等式与 II.b.2 的 F 恒等式用**自由符号 u_m**（sympy，K=0..5）验证为纯代数恒等式（不用 u 的任何性质）；c_m 凸/递减、d_m 凸/递减、(\*)、2Ku_Ku_{K+1}<1 均核过；自写暴力 K≤300、L≠L′≤300：min W=0，W=0 ⟺ K=0 ∨ L≥L′+K+1 **逐格精确**，L<L′ 且 K≥1 时 min 1.43e-10（K=1，(298,300)）；末行条件 (Me_{L′})_n≥0（n≥L′）逐格成立；δ₃ 对照 (2,1) = −0.4472 复现 |
| 3 §3.A 步骤 1（半群保锥） | **CONFIRMED** | 阶梯基坐标 a_L(y)=y_{L−1}−y_L（L<N）、a_N=y_{N−1} 下 C_N=非负象限；截断 P_N𝔎P_N 的系数矩阵 B 的非对角元恰是 W_{L,L′}（L<N）与 (Me_{L′})_{N−1}（L=N，由情形 I ≥0）；Abel 分解 F=Σ_K(F(K)−F(K+1))1_{k≤K} 逐元合法（每个矩阵元是 F 的有限线性组合，Σa_K=F(1)<∞，F(K)→0 由控制收敛）⟹ B 为非负权 Metzler 和 ⟹ e^{tB}≥0 逐元；N→∞ 强收敛 + C 闭成立。自写 N=600：随机 P 的 𝔎_T 阶梯基非对角最小值 +0.0、末行最小 +3.3e-8；e^{t𝔎}v 对 120 组 (v∈C, t∈{0.1,1,5,25}) 全部留在 C；对照 δ₃ 核 e^{𝒪}e_1 出锥（max diff +0.474） |
| 4 §3.A 步骤 2–3 | **DOWNGRADE（两处证明文字须改，结论不变）** | (a) **C−C 稠密的构造写错**：「x=(M·1_{[0,n₀)})−(M·1_{[0,n₀)}−x)」中第二项非负但**不必非增**（反例 x=(1,0,1,0)，M=10：Me₃−x=(9,10,9,0)）。正确构造：x=Σ_{L≤n₀}b_Le_L（b_L=x_{L−1}−x_L），拆 b=b⁺−b⁻，两半都在 C。(b) 「Σ_kf(k)Σa_La_{L′}(abs Q) 收敛因 ⟨v′,𝒪[f]v′⟩ 有限——用 Tonelli」不成立为论证（(n,m) 二重和含负号的 L1 项，有限≠绝对收敛）；正确做法：⟨v,𝒪[f]v⟩=lim_N⟨P_Nv,𝒪[f]P_Nv⟩（有界二次型连续），P_Nv 是有限 Abel 和，有限和内 ≥a_{L₀}²f(1)Q_{L₀L₀}(δ₁)（N>L₀），取极限。其余：s-lim e^{t(𝔎−λ₁)}=E({λ₁}) 由谱定理 + 控制收敛对任何有界自伴成立、不需孤立，只需 λ₁ 是本征值（紧 + λ₁≥𝔎₀₀=F(1)>0，后者文中未写但显然）；P₁C⊂E₁∩C、P₁C≠{0}、Rayleigh 比较、f≥0 与 f(1)>0、严格性均成立；不需 P 紧支、不需 λ₁ 单重、不需 Krein–Rutman 属实 |
| 5 数值对拍 | **CONFIRMED** | N=600 自写：40 组随机 1–11 原子 P × 5 个 T（160 对 T<T′）λ₁(𝔎_T)>λ₁(𝔎_{T′}) 零违反，min diff 3.1e-5；三组极端 P（单原子 T∈[1e-4,100]、两原子 ω=0.01/100 权 0.5/0.5 与 0.999/0.001）严格递减；顶本征向量 54/54 在 C（含 T=1e-4、T=100、极端两原子；tol 1e-10，最大违反 0）；⟨v,𝒪[f]v⟩ 对 300 组随机 (v∈C, f≥0 非单调) min +7.6e-5；对照非锥向量 e₅ 给 −0.18 |
| 6 复现 | **CONFIRMED** | `python p4_uniqueness_lemmas.py`：ALL 37 CHECKS PASSED，N=500，26.3 s；diff 只触及两个预期文件 |
| 7 记账 | **CONFIRMED（一处措辞）** | μ\*=0 限定、「Conjecture 2 完整」、Einstein h 于 (0,∞) 严格递增均无过度声称：KAY II Conj. 2 陈述（survey 行 36 逐字）本身无 μ\* 项；「完整」= 界（§1 已审）+ 存在（引理 X）+ 唯一（本文），恰是「There is a critical temperature … bounded above by」的全部内容。库仑束 Perron 向量 58/75 出锥的探针数字**本轮未复核**（只是限定范围的陈述，不承重）。措辞：§3 首段「λ₁′>0」与 §3.A「λ₁ 是本征值（紧性）」应补一句 λ₁≥𝔎₀₀=F_T(1)>0 |

---

## 1. 引理 A · CONFIRMED

逐元定义 𝒪[δ_k]_{nm}=u_nu_m[1_{|n−m|=k}+1_{n+m+1=k}]−2u_n²1_{n=m,n≥k}，L≤L′，矩形 n<L、m<L′：
- m=n+k：n<L 且 n+k<L′ ⟹ n<min(L,L′−k) ✓；m=n−k：n≥k、n<L，且 m=n−k<L≤L′ 自动 ✓；L3：m=k−1−n∈[0,L′) ⟺ k−L′<n≤k−1，并 n<L ⟹ n∈[max(0,k−L′),min(L,k)) ✓；对角 (n,n) 在矩形内 ⟺ n<L=min ✓。
- min(L,L′−k)≥L−k（L≤L′）⟹ m=n+k 支拆成 Σ_{n<L−k}u_nu_{n+k}=Σ_{j=k}^{L−1}u_{j−k}u_j 与余下 [max(0,L−k),min(L,L′−k)) = 组 2；前者与 m=n−k 支、L1 合成组 1 = 2Σ_{j=k}^{L−1}u_j(u_{j−k}−u_j) ✓。L≤k 时组 1 与「前者」同时为空，一致。
- 取等：组 1 空 ⟺ L≤k；组 2 空（在 L≤k 下）⟺ L′≤k；组 3 空（在 L≤k 下）⟺ L≤k−L′。合取 ⟺ L+L′≤k ✓。

自写暴力（`rev8_lemmas.log`）：k,L,L′≤300 全部 2.7e7 格 min Q=0、取等集 == {L+L′≤k} 逐格精确、其余 min 3.3334e-3；闭式 vs 矩阵抽样 max diff 5.3e-15；min_L Q_{LL}(δ₁)=1.000000。

## 2. 引理 B · CONFIRMED

行和 g(n)=(Me_{N′})_n 的推导：L2 窗口 m∈[max(0,n−K),min(N′,n+K+1)) 去掉 m=n（当 n<N′）贡献 u_n[S(·,·)−u_n1_{n<N′}]；L3 m<min(N′,K−n) 贡献 u_nS(0,min(N′,K−n))；对角 −2min(n,K)u_n²1_{n<N′}；合并 −u_n²(1+2min(n,K))1_{n<N′} ✓。
- 情形 I（n≥N′）：两行无对角，R(n)=S(max(0,n−K),N′)+S(0,min(N′,K−n)) 两片随 n 非增，u_n 递减、R≥0 ⟹ g 非增 ✓；R(n)=0 ⟺ 两片皆空 ⟺ n−K≥N′ ✓；g(n)=u_nR(n)≥0 即末行条件 ✓。
- II.a（n≤K−1，n+1<N′）：两行下限均为 0，R(n+1)−R(n)=u_{n+K+1}1_{n+K+1<N′}−u_{K−n−1}1_{K−n−1<N′}；K−n−1≤n+K+1 ⟹ 若前者示性为 1 则后者亦为 1，差 ≤0 ✓；W=u_nR(n)−u_{n+1}R(n+1)≥d_nR(n)≥d_nu_0 ✓（m=0 落在 L3 片内因 K−n≥1、N′≥2）。
- II.b（n≥K）：两行对角同为 2K+1，L3 片空。II.b.1：sympy 自由符号验证 g(n)−g(n+1)=Σ_{j=1}^K[c_{n−j}+c_{n+j}−2c_n] 为**代数恒等式**（K=0..5），只在 c 凸处用 u 凸、|u′|=(2x+1)^{−3/2} 凸 ✓。II.b.2：W(N′) 表达式与 F=Σ_{m=n−K}^n c_m+u_nu_{n+1}−(2K+1)c_n 同样自由符号验证 ✓；N′→N′+1 增量 d_nu_{N′} 在 n+2≤N′≤n+K 成立（N′=n+K+1 之后已属 II.b.1，不需要）✓；c 递减 ⟹ F≥u_nu_{n+1}−Kc_n=u_nu_{n+1}(1−2Ku_nu_{n+1})，2Ku_Ku_{K+1}=2K/√((2K+1)(2K+3))<1 ✓（K≤300 数值最大 0.9967）。
- 穷尽：L>L′ / L<L′；后者 n≤K−1 / n≥K；n≥K 下 N′≥n+K+2 / n+2≤N′≤n+K+1 ✓。

自写暴力：K≤300、L≠L′≤300 min W=0，取等集 == {K=0 ∨ L≥L′+K+1} 逐格精确，L<L′、K≥1 最小 1.4266e-10（K=1，(298,300)），末行 (Me_{L′})_n（n≥L′）min 0，δ₃ 对照 −0.4472。

## 3. 步骤 1 · CONFIRMED

「阶梯基 Metzler」确是所需条件：C_N 在坐标 a_L(y)=y_{L−1}−y_L（L<N）、a_N=y_{N−1} 下是非负象限，e^{tB}=lim(I+tB/n)^n 对 Metzler B 逐元非负，故保象限。截断矩阵 P_N𝔎P_N 的第 L′ 列系数只涉及 (𝔎e_{L′})_{L−1}、(𝔎e_{L′})_L（L≤N−1 时两者都在截断内，等于无穷矩阵的 W_{L,L′}）与 (𝔎e_{L′})_{N−1}（L=N），后者由情形 I ≥0——§2 推论 B′ 与 §3.A 的说法一致。Abel 和 Σ_K a_K𝒪[1_{k≤K}]：每个矩阵元是 F 的有限组合（L1 对角 Σ_{k≤n}、L2/L3 单项），Σ_{K≥k}a_K=F(k)−lim F=F(k)（⟦K⟧→0 有界，控制收敛），故逐元与 Σ_K 交换无条件；B(F)=Σa_KB(K) 非对角 ≥0 ✓。N→∞：‖e^{t𝔎^{(N)}}P_Nv−e^{t𝔎}v‖≤e^{t‖𝔎‖}‖P_Nv−v‖+‖(e^{t𝔎^{(N)}}−e^{t𝔎})v‖→0（幂级数逐项强收敛 + 范数尾一致）✓。数值见总表。

## 4. 步骤 2–3 · DOWNGRADE（证明文字两处）

1. **C−C 稠密（§3.A 步骤 2）构造错误，CONFIRMED**：文中取 x=(M·1_{[0,n₀)})−(M·1_{[0,n₀)}−x)，第二项 (M−x_j)_{j<n₀} 非增 ⟺ x 在 [0,n₀) 上非减，一般不成立；`rev8` 反例 x=(1,0,1,0)、M=10 给 (9,10,9,0)。**修复**（一行）：有限支撑 x=Σ_{L≤n₀}b_Le_L，b_L=x_{L−1}−x_L，x=Σb_L⁺e_L−Σb_L⁻e_L∈C−C。结论「C−C ⊃ 有限支撑序列 ⟹ 稠密 ⟹ P₁C≠{0}」仍成立。
2. **Tonelli 句（步骤 3）论证不成立但结论成立，CONFIRMED**：从 Σ_{n,m}v_n𝒪_{nm}v_m 重排到 Σ_{L,L′}a_La_{L′}Q_{LL′} 需绝对收敛；「⟨v′,𝒪[f]v′⟩ 有限」不蕴含它（L1 项为负）。直接修法：二次型连续 ⟹ ⟨v,𝒪[f]v⟩=lim_N⟨P_Nv,𝒪[f]P_Nv⟩，P_Nv=Σ_{L<N}a_Le_L+v_{N−1}e_N 有限和，⟨P_Nv,𝒪[f]P_Nv⟩=Σ_{L,L′≤N}a_La_{L′}Q_{LL′}(f)，Q_{LL′}(f)=Σ_{k≤L+L′−1}f(k)Q_{LL′}(δ_k)≥0（有限和），故 ≥a_{L₀}²f(1)Q_{L₀L₀}(δ₁)（N>L₀）；极限保 ≥。（或者：𝒪_abs[f]=L1+L2+L3 有界（`P4_conjecture2.md` §6 的 Rsum 界），⟨v,𝒪_abs v⟩<∞ 给绝对收敛——也行，但要引用。）
3. 谱投影：s-lim_{t→∞}e^{t(A−λ₁)}=E({λ₁}) 对任何有界自伴 A、λ₁=max spec 成立（‖(e^{t(A−λ₁)}−E{λ₁})v‖²=∫_{μ<λ₁}e^{2t(μ−λ₁)}d⟨v,E_μv⟩→0，控制收敛），**不需孤立**；需要的只是 E({λ₁})≠0，即 λ₁ 为本征值：紧 + λ₁>0。λ₁>0 文中两处（§3「λ₁′>0」、§3.A）未给理由，一句即够：λ₁≥⟨e₁,𝔎e₁⟩=𝔎₀₀=F_T(1)u_0²=F_T(1)>0（n=0 无 L1 项）。紧性：F(k)≤⟨ϖ²⟩/k² ⟹ Hilbert–Schmidt（或引 KAY II Prop 2）。
4. 其余核过：P₁v=lim e^{−tλ₁}e^{t𝔎}v∈C（正标量、闭）且 ∈ran P₁；v′≠0 ⟹ 某 a_L>0（否则 v′ 常值且 ∈ℓ² ⟹ 0）；𝒪[F_T]−𝒪[F_{T′}]=𝒪[F_T−F_{T′}]（逐元线性、两者有界）；f≥0 由 ⟦k⟧ 关于 ϖ 递增、f(1)>0 由 P(ω>0)>0；Rayleigh λ₁(𝔎_T)≥⟨v′,𝔎_Tv′⟩。全程未用单重性、Krein–Rutman、紧支 ✓。Einstein 特例即 F=⟦·⟧(ϖ)，得 h 于 (0,∞) 严格递增 ✓。

## 5. 数值 · CONFIRMED
见总表行 5；`rev8_main.log`。补充：极端 P（ω=0.01 与 100 权 0.999/0.001，T=0.5）的 λ₁≈4.0e-3、λ₂≈−8e-7，顶本征向量 min 7.1e-5>0、diff max −3.1e-7<0，锥内且严格。

## 6. 复现 · CONFIRMED
`python theory/derivations/p4_uniqueness_lemmas.py`：ALL 37 CHECKS PASSED (N=500, 26.3 s)。`git diff 916296b 5f05f83 --stat`：`P4_uniqueness_strategy.md` +51/−4、`p4_uniqueness_lemmas.py` +331，无其他文件。

## 7. 总评

**「对所有 P∈𝒫（P((0,∞))=1、⟨ω²⟩<∞），T↦𝔨(P,T) 在 μ\*=0 下严格递减」——PROVED。** 引理 A、B 的初等证明逐指标重推无误，关键恒等式（II.b.1 凸性和、II.b.2 的 F）是自由符号下的代数恒等式，取等集合与我方独立暴力（≤300）逐格精确一致；主论证的骨架（阶梯基 Metzler ⟹ 半群保锥 ⟹ 顶本征空间含锥向量 ⟹ Rayleigh 比较 + 推论 A′ 严格性）每一步成立，且确实不需要单重性、Krein–Rutman、紧支或 λ₁ 孤立。须修改的是两处**证明文字**而非结论：步骤 2 的 C−C 稠密构造（写错，一行 Abel 拆分修复）与步骤 3 的 Tonelli 句（改为有限截断取极限）；另补 λ₁>0 的一句理由。修改后，配合引理 X 与 §1 定理，KAY II Conjecture 2 在 μ\*=0 下完整成立（存在、唯一、T_c<T̃_c），Einstein 的 h 在全 (0,∞) 严格递增亦随之成立；Stage M 证书降为推论。无 CHALLENGED。

复现脚本（`theory/notes/_scratch_ir/`）：`rev8_lemmas.py` / `rev8_lemmas.log`（自写 𝒪 构造；引理 A/B 暴力 ≤300、闭式抽样、取等集、末行、δ₃ 对照、自由符号恒等式、凸性；1.4 s）、`rev8_main.py` / `rev8_main.log`（N=600：160 对 T<T′ 单调、极端 P、54 顶本征向量入锥、120 组半群保锥、δ₃ 出锥对照、阶梯基 Metzler、300 组二次型正性；26 s）。原脚本复跑 26 s。
