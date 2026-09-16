# Kondo 屏蔽的"对数发散 → 不敏感性"探针

*建立：2026-09-16，任务来自 team-lead 的 Kondo 普适性探索委托。*
*性质说明：本笔记是文献 + 解析论证，不含仓库自己的证书/数值管线；结论一律标
**[heuristic]**（标准教科书/文献论证，未在本仓库重新证明）或 **[literature]**（转引实验/数值结果）。
不与 `P4_*` 系列的 **[analytic]/[certificate]/[numerical-evidence]** 三级混淆。*

---

## 0. 结论先行

**Kondo 系统确实展现与 BCS 项目核心命题同构的"对数发散 → 不敏感性"现象，但边界条件更窄、且有一个已知的反例（van Hove 奇点）把它精确地击穿——这个反例本身很有信息量，因为它标出了不敏感性成立所需要的确切前提。**

具体地：

1. **成立部分**：一阶（poor-man's scaling）和精确 Bethe ansatz 层面，$T_K$ 的指数结构只依赖 $\rho(E_F)$（费米面处的态密度**值**），不依赖 $N(\varepsilon)$ 在带内的**形状**——只要 $\rho$ 在从 $D$ 标度到 $T_K$ 的整个 RG 窗口内近似为常数。这与 BCS 的 Cooper 对数 $\int N(\varepsilon)\tanh(\varepsilon/2T)/\varepsilon\,d\varepsilon\approx N(0)\ln(\omega_D/T)$ 只看 $N(E_F)$ 是同一个数学机制。
2. **前因子非普适**：连接微观参数（$U,\Gamma,\varepsilon_d$，或者 cutoff 方案）到 $T_K$ 的 $O(1)$ 前因子是**方案依赖**的——Haldane 近似公式与精确 Bethe ansatz/NRG 结果在前因子上有历史性的分歧，指数本身分毫不差。这与 BCS 项目里 "λ、⟨ω²⟩ 是相关方向，前因子/高阶矩是无关方向" 的 RG 语言完全平行。
3. **反例**：当费米面恰好落在（或极靠近）van Hove 奇点上，$N(\varepsilon)$ 在 RG 窗口内本身发散，不敏感性**按构造失效**——NRG 计算显示 $1/\ln T_K$ 从标准的 $\propto J$（线性）变成 $\propto\sqrt J$，即 DOS 的形状渗透进了指数本身，不再只是前因子修正。这精确标出了"不敏感性"定理需要的前提：**费米面附近 DOS 局部正则（不发散），且相关能标 $T_K\ll D$ 有真实的尺度分离**——与本仓库 Theorem E 需要 $\varpi\gg1$（真实声子谱窄、尺度分离）的前提是同一句话的两个物理实现。
4. **实验分散度的来源与预期一致**：同一宿主中不同杂质种类（Cr/Mn/Fe/Co/Ni in Cu）$T_K$ 跨越好几个数量级，但这**不是**因为不同 3d 元素让 $\rho(E_F)$（宿主的态密度）显著变化——noble metal host 的 $\rho(E_F)$ 变化通常不到 2 倍——而是因为**耦合常数 $J$**（由局域能级 $\varepsilon_d$、$U$、杂化 $\Gamma$ 决定）逐元素剧烈变化，指数放大了这个差异。这恰好印证"不敏感性"命题的精确内容：**对谱形不敏感，但对耦合常数本身极度敏感**——后者进指数，前者不进（在正则 DOS 前提下）。
5. **Kondo lattice 打破简单图像**：相干温度 $T_{\rm coh}$（Kondo lattice 的第二个能标）不是单杂质 $T_K$ 的简单推广，它对晶格细节、RKKY、维度、掺杂高度敏感——这是"多体自洽反馈"重新引入形状依赖的例子，与 BCS 从弱耦合 Theorem E（严格）过渡到强耦合 Eliashberg 修正（形状开始起微弱作用）是同一个物理直觉的镜像。

---

## 1. 理论构造：Anderson 杂质 → Kondo 温度

### 1.1 模型与 Schrieffer–Wolff

Anderson 杂质模型：
$$
H=\sum_{k\sigma}\varepsilon_k c^\dagger_{k\sigma}c_{k\sigma}+\varepsilon_d\sum_\sigma n_{d\sigma}+Un_{d\uparrow}n_{d\downarrow}+\sum_{k\sigma}V_k\big(c^\dagger_{k\sigma}d_\sigma+{\rm h.c.}\big).
$$
在局域矩区（$-\varepsilon_d,\ \varepsilon_d+U\gg\Gamma:=\pi\rho(E_F)V^2$）通过 Schrieffer–Wolff 变换积掉电荷涨落，得到 Kondo 模型的有效交换：
$$
J=2V^2\Big[\frac1{-\varepsilon_d}+\frac1{\varepsilon_d+U}\Big]\ \xrightarrow{\ \varepsilon_d=-U/2\ }\ J=\frac{8V^2}{U}.
$$
**[heuristic，标准 Schrieffer-Wolff]**。关键点：$J$ 由**局域、高能**参数（$U,\Gamma,\varepsilon_d$）决定，$\rho(E_F)$ 只通过 $\Gamma=\pi\rho V^2$ 这一个数进入——没有任何 $N(\varepsilon)$ 形状信息留在 $J$ 里。这一步已经暗示了后面的不敏感性：**微观模型把"谱形"和"耦合常数"从一开始就分离开**，正如 BCS 中 $\alpha^2F(\omega)$ 的信息只通过 $\lambda=2\int\alpha^2F(\omega)/\omega\,d\omega$ 这类矩进入 leading order。

### 1.2 Poor-man's scaling：为什么只看 $\rho(E_F)$

一阶标度方程 $dJ/d\ln D=-2\rho J^2$ 的推导过程：每降低一次 cutoff $D\to D-|dD|$，积掉的虚过程（电子/空穴被局域矩散射到能量 $D$ 附近再散射回来）贡献一个 $\rho(D)\,dD/D$ 的对数项。**只要 $\rho(\varepsilon)\approx\rho(E_F)$ 在整个降标度路径 $D\to T_K$ 上近似常数，$\rho(D)$ 就可以提到积分外**，给出
$$
T_K=D\exp\!\Big(-\frac1{\rho J}\Big)\times[\text{prefactor}].
$$
——这与 BCS 中 $\chi_{\rm pair}(0)=\int N(\varepsilon)\tanh(\varepsilon/2T)/\varepsilon\,d\varepsilon\approx N(0)\ln(\omega_D/T)$（Cooper 对数）把 $N(\varepsilon)$ 提到积分外的条件**是同一个假设**：相关能标远小于 DOS 变化的特征能标。两者都是"费米面对数积分对全局谱形的不敏感性"的具体实现，机制上是同一件事：
$$
\int_{T_K\ (\text{或 }T_c)}^{D\ (\text{或 }\omega_D)}\!\!\rho(\varepsilon)\,\frac{d\varepsilon}\varepsilon\ \approx\ \rho(E_F)\ln\frac D{T_K}\quad\text{当 }\rho\text{ 在窗口内近似常数}.
$$

### 1.3 精确 Bethe ansatz：前因子非普适，指数普适

对称 Anderson 模型的精确（Bethe ansatz）Kondo 标度：
$$
T_K \sim \sqrt{U\Gamma}\,\exp\!\Big(-\frac{\pi U}{8\Gamma}\Big)
$$
（Haldane 近似式，形式上），指数 $-\pi U/(8\Gamma)$ 与精确 Bethe ansatz/Wilson NRG 结果**完全一致**，但 $O(1)$ 前因子历史上有过修正——Haldane 1978 的原始前因子被后续 NRG（Krishna-murthy–Wilkins–Wilson 1980）和精确 Bethe ansatz（Tsvelick–Wiegmann, Andrei–Furuya–Lowenstein 1983）确认存在数值偏差 **[literature]**。教训完全平行于本仓库的教训：**指数（leading log 系数）是不敏感性/普适性成立的部分；$O(1)$ 前因子是 cutoff 方案（sharp vs Lorentzian vs Anderson 模型自身的能级结构）依赖的部分，从来没有人声称前因子普适**。这正是 BCS 项目里 "$\lambda,\langle\omega^2\rangle$ 是 relevant 方向，Jensen gap 是 bounded-but-nonzero 的 irrelevant 修正" 的逐字翻译。

来源：一般弱耦合渐近同时给出 $T_K/D\sim(\rho J)^n\exp[-1/(\rho J)]$ 型修正（$n$ 依赀阶数和 cutoff 方案变化于 $1/2$ 与 $2$ 之间），指数本身在所有方案下一致 **[literature, Bethe ansatz 综述]**。

---

## 2. Wilson NRG：比 BCS 更强的普适性——以及它在哪里断裂

### 2.1 单参数标度：$T/T_K$ 普适函数

Wilson 的数值重整化群给出一个比"$T_K$ 只依赖 $\rho J$"更强的陈述：**一旦 $T_K$ 固定，Kondo/Anderson 模型 Kondo 区的全部低能热力学量都坍缩成 $T/T_K$ 的普适函数** **[literature]**：
$$
T\chi(T)=\varphi(T/T_K),\qquad w:=4T_K\chi(0)\approx0.4128,\qquad R_W:=\frac{4\pi^2}{3}\frac{\chi_{\rm imp}(0)}{\gamma_{\rm imp}}=2\ (\text{精确 Bethe ansatz}).
$$
这比 BCS 项目的命题**强得多**：BCS 项目证的是"$T_c$ 这一个数对谱形不敏感"（Theorem E），而 Kondo 的单参数标度是"**所有**低能可观测量的**整条曲线**对模型的微观细节不敏感，只认 $T_K$"。原因是 Kondo 不动点是一个**精确可解、无长程序参量自洽反馈**的量子杂质不动点（等价于一个零维/边界共形场论问题），没有 BCS 那种"序参量 $\Delta$ 反过来重整化自身能隙方程"的自洽循环——这是 Kondo 普适性比 BCS 普适性"更干净"的结构原因，与本探针任务要求对比的"单粒子 vs 多体"区别直接相关（见 §3）。

### 2.2 断裂点一：van Hove 奇点——指数本身被 DOS 形状渗透

这是本探针**最重要的反例**，因为它精确标出不敏感性定理成立所需要的前提。$t$–$t'$ 方晶格 Kondo 模型（NRG 研究，arXiv:1109.1700）在 $t'\to0$ 时费米面精确落在对数发散的 van Hove 奇点上：

- 结果：**$1/\ln T_K$ 从标准的线性依赖 $\rho J$（远离奇点，$t'$ 较大，DOS 局部近似常数）交叉到 $\sqrt{\rho J}$ 依赖（$t'\to0$，费米面在奇点上）** **[literature]**。
- 物理解读：§1.2 的推导要求 $\rho(\varepsilon)$ 在整条 RG 路径上近似常数才能把它提出积分号；van Hove 奇点让 $\rho(\varepsilon)\sim\ln(1/|\varepsilon|)$ 在窗口内自身发散，积分 $\int\rho(\varepsilon)d\varepsilon/\varepsilon$ 不再简化为 $\rho(E_F)\ln(D/T_K)$，而是产生一个**双对数**结构，指数从 $-1/(\rho J)$ 变成 $\sim-1/\sqrt J$ 型——**这不是前因子修正，是指数的函数形式本身改变**。
- 同时报道：$t'\to0$ 处 Fermi 液体行为本身破缺（准粒子描述失效），暗示这里连"用 $T_K$ 做唯一标度"这件事的前提（存在唯一交叉能标+局域费米液体基态）都动摇了 **[literature]**。

**与 BCS 项目的镜像**：这精确对应本仓库对 Theorem E 的诚实边界——弱耦合级数、Lemma C 凹性都需要 $\varpi\gtrsim$ 某阈值（尺度分离），且团队自己发现"极端谱形（接近 delta 函数集中在带边、或强烈非 Einstein 的形状）"是压力测试 Theorem E 会不会推广到 $\mu^*>0$ 的地方。van Hove 奇点是 Kondo 版本的"病态谱形"：它不是"形状稍微偏离 Einstein"，而是"DOS 本身在相关能标发散"，恰好对应"什么时候不敏感性命题的前提被违反"这一类问题的一个独立确证。

### 2.3 断裂点二：输运性质不是普适的

热力学量（$\chi,\gamma,S$）在 Kondo 区普适于 $T/T_K$，但**输运**（电导、电阻率的温度依赖）一般**不是** $T/T_K$ 的普适函数——依赖于势散射相位、粒子-空穴不对称度等额外参数 **[literature]**。这也是一个"哪些量继承不敏感性、哪些不继承"的清晰划分，与 BCS 项目区分"$T_c$（可以证明不敏感）vs $2\Delta/T_c$ 之类比值（强耦合下依赖更多谱信息）"是同一类划分。

---

## 3. 机制对比：单粒子对数（Kondo）vs 多体 Matsubara 对数（BCS）

这是委托任务明确要求回答的问题。两个对数在**数学结构**上同源，在**物理层级**上不同：

| | Kondo $T_K$ | BCS $T_c$ |
|---|---|---|
| 对数的物理来源 | $T$-矩阵微扰级数里，局域矩把导带电子虚散射到 $E_F$ 附近的中间态（粒子或空穴），对中间态能量求和给出 $\int d\varepsilon/\varepsilon$ | Cooper 对易感受率 $\chi_{\rm pair}(0)=T\sum_{\omega_n}\int N(\varepsilon)\,G(\varepsilon,\omega_n)^2$，对 Matsubara 频率/能量求和给出 $\ln(\omega_D/T)$ |
| 问题的"体积" | 量子杂质：0 维（杂质自旋）+ 连续谱库；等价于边界 CFT / 单通道散射问题 | 真正的体系统：整个费米海参与配对，宏观序参量 $\Delta$ |
| 是否存在自洽反馈 | 否——RG 流是单方向的（$J$ 流向强耦合，无序参量反过来修改流），Kondo 不动点是**精确可解**的（Bethe ansatz） | 是——能隙方程 $\Delta=\lambda\int\cdots\tanh(E/2T)/E\,d\varepsilon$ 本身含 $\Delta$，$T_c$ 是这个自洽方程的**线性化**（$\Delta\to0$ 极限）临界点 |
| 求和/积分的变量 | 单个中间态能量 $\varepsilon$（单粒子传播子的极点） | 玻色型 Matsubara 频率求和 + 费米型能量积分（两个费米子传播子的乘积，配对涨落是玻色型的） |
| 不敏感性覆盖范围 | 覆盖**整条低能标度曲线**（Wilson 单参数标度），因为不动点本身对微扰是渐近自由/渐近强耦合、无竞争相 | 只覆盖 $T_c$ 这一个数（线性化间隙方程的本征值问题），$2\Delta/T_c$ 等衍生比值在强耦合下**会**重新依赖谱形（Eliashberg 阶数修正） |
| 已知反例 | van Hove 奇点（本笔记 §2.2）：DOS 发散 → 指数函数形式改变 | 强耦合、极端非 Einstein 谱形（宽/多峰）：Jensen gap 不再严格为零，只是"结构性小"（`paper_tc_insensitivity_core_statement_2026-09-16.md` §3） |

**关键判断**：委托任务问"单粒子 vs 多体的区别是否影响不敏感性"——回答是：**不影响不敏感性命题本身成立的条件（两者都归结为同一个"费米面附近对数积分只认局部 DOS 值"的引理），但影响不敏感性覆盖的"范围"**。Kondo 因为没有自洽反馈、且是可积不动点，把不敏感性从"一个能标"提升到"整条普适曲线"；BCS 因为间隙方程的自洽性，只能在线性化（$T_c$ 本身）层面拿到严格定理，非线性衍生量（$\Delta(0)$、比值、强耦合修正）需要额外工作（本仓库现有 Theorem E 到 Eliashberg 修正的过渡正是这条线）。这也解释了为什么 Kondo 的"反例"（van Hove）比 BCS 已知的"反例"（强耦合谱形修正）更"硬"——前者是把指数本身改变的**结构性**违反，后者是"仍然是小修正，只是不再精确为零"的**量级性**违反。

---

## 4. 实验/数值证据

### 4.1 同一宿主、不同杂质：分散度由 $J$ 主导，不由 $\rho$ 形状主导

3d 过渡金属杂质（Cr, Mn, Fe, Co, Ni）掺入 Cu、Au 等贵金属宿主，$T_K$ 跨越多个数量级（经典例子：Cu:Mn 的电阻极小值出现在 mK–K 量级，Cu:Fe 在 10–30 K 量级）**[literature]**。贵金属宿主的 $\rho(E_F)$ 在这些体系间变化通常不到 2 倍，远不足以解释指数级的 $T_K$ 差异；差异的主导来源是**局域参数** $\varepsilon_d,U,\Gamma$（进而 $J$）随 3d 壳层填充数系统变化——这正是"进指数的量（$J$）敏感，不进指数的量（$\rho$ 的形状）不敏感"这条分界线的直接实验确认 **[literature]**。

一个更尖锐的确认：**同一元素、不同局域环境**——Co 杂质本体（bulk Cu 内部）与 Co 在 Cu$_2$N/Cu(100) 表面配对时，$T_K$ 从体相的更高值降到 $\approx2.6\,$K **[literature]**。宿主 3d 电子态密度的"全局形状"没有发生数量级变化，改变的是**局域杂化 $\Gamma$**（配体环境改变了 $V_k$，从而改变 $\Gamma=\pi\rho V^2$ 里的 $V$，不是改变 $\rho$ 本身的形状）——再次确认真正敏感的自由度是耦合常数而非谱形。

同时也有一个警示性的反例（同一批文献报道）：用简化标度公式 $T_K\approx0.4D_0\exp(\cdots)$ 对 Co in Cu 做粗估会比 QMC/实验高出一个数量级以上 **[literature]**——这提醒"前因子/cutoff 定义非普适"这条边界在实际材料估算里是**实打实的误差来源**，不是纯理论洁癖；与本仓库反复强调"承重数值结论要独立重算，不能只信公式的渐近形式"是同一条纪律的物理版本。

### 4.2 Kondo lattice：相干温度打破简单图像

CeCoIn$_5$ / CeRhIn$_5$（同构 "115" 家族，仅过渡金属层不同：Co vs Rh）在环境压力下行为迥异——CeCoIn$_5$ 是超导体（$T_c\approx2.3\,$K）且更靠近量子临界点，CeRhIn$_5$ 常压反铁磁（$T_N\approx3.8\,$K），需要 $\gtrsim1.15\,$GPa 压力才出现超导 **[literature]**。差异的公认机制是 **c–f 杂化强度 $V$ 不同**（Rh/Co 的 d 电子结构不同导致 Ce-4f 与配体轨道的重叠不同），这在 Doniach 相图上表现为两个化合物处在不同的 $J_{\rm RKKY}/J_{\rm Kondo}$ 位置——压力通过压缩晶格增大 $V$，从而指数地推高 Kondo 关联能标、抑制 RKKY 主导的磁序。这仍然是**耦合常数（$V$，进而 $J$）驱动**的效应，与单杂质情形一致；但相干温度 $T_{\rm coh}$（heavy Fermi liquid 出现的温度）本身**不是**单杂质 $T_K$ 公式的简单推广——它由晶格周期性、RKKY 关联长度、掺杂/无序等**額外**自由度共同决定，这些自由度在单杂质图像里根本不存在。**这是"多体自洽反馈重新引入形状/细节依赖"的具体例子**：单杂质 Kondo 是可积、无竞争序参量的问题，Kondo lattice 一旦存在 RKKY 竞争和潜在磁序，就重新变成一个类似 BCS 的"多体自洽问题"，不敏感性的强版本（Wilson 单参数标度）不再自动成立。

---

## 5. 回答委托任务的具体问题

1. **$T_K$ 是否只依赖 $\rho,J$、对 $N(\varepsilon)$ 的能量依赖不敏感？**
   — 是，在"$\rho$ 在 RG 窗口内局部正则（不发散）"这一前提下，指数层面严格成立（poor-man's scaling 解析推导 + Bethe ansatz 精确解交叉确认）；前因子层面不严格（scheme 依赖，Haldane vs NRG 历史分歧即为例证）；**van Hove 奇点是已知的、指数层面本身被击穿的反例**，标出了这条不敏感性命题需要的确切前提（局部正则 DOS + 尺度分离 $T_K\ll D$）。
2. **Kondo lattice 相干温度对晶格细节的敏感性？**
   — 敏感。$T_{\rm coh}$ 依赖 RKKY、晶格周期性、维度、掺杂/无序，不是单杂质 $T_K$ 公式的直接推广；CeCoIn$_5$/CeRhIn$_5$ 的差异是耦合常数（杂化 $V$）驱动而非 DOS 形状驱动，但相干温度本身的**存在性和数值**依赖晶格层面的额外结构。
3. **数值/实验证据：不同杂质 $T_K$ 比值 vs 理论预期？**
   — 一致：贵金属宿主 $\rho(E_F)$ 变化很小，3d 杂质间 $T_K$ 的数量级差异主要由局域 $J$（$\varepsilon_d,U,\Gamma$）驱动，指数放大了这个差异；同元素不同局域环境（bulk vs 表面配对）差异同样由 $\Gamma$ 驱动，不由宿主 DOS 全局形状驱动。

**总评**：Kondo 系统是本仓库"对数发散 ⟹ 谱形不敏感性"这一大类现象的**第二个独立物理实现**（第一个是 BCS/Eliashberg 本身），机制上同源（费米面附近对数积分只认局域 DOS 值），覆盖范围因"是否存在自洽反馈"而不同（Kondo 更强，覆盖整条标度曲线；BCS 只严格覆盖 $T_c$ 本身），且都有一个结构清晰、有物理意义的"反例边界"（Kondo：van Hove 奇点；BCS：强耦合下的 $O(\text{Jensen gap})$ 修正）。这支持把"对数发散→不敏感性"当作一个**跨系统的一般性现象**来讨论，但任何推广陈述都必须显式带上"DOS/谱在相关能标局部正则 + 存在真实尺度分离"这条前提——这条前提在两个系统里都是不可省略的。

---

## 参考来源（均为 [literature]，未逐一验证原始论文全文，仅转引搜索摘要/公认教科书结果）

- Haldane, "Scaling Theory of the Asymmetric Anderson Model", PRL 1978（近似 $T_K$ 标度公式及其前因子历史订正）。
- Krishna-murthy, Wilkins, Wilson, PRB 21, 1003/1044 (1980)（NRG 对 Anderson 模型的数值解，普适标度曲线、Wilson ratio）。
- Tsvelick & Wiegmann, Adv. Phys. 32, 453 (1983); Andrei, Furuya, Lowenstein, Rev. Mod. Phys. 55, 331 (1983)（Bethe ansatz 精确解，$T_K$ 精确前因子与普适比值 $R=2$）。
- arXiv:1109.1700, "The Kondo effect in the presence of Van Hove singularities: A numerical renormalization group study"（$1/\ln T_K$ 从线性到 $\sqrt J$ 的交叉，Fermi 液体行为在 $t'\to0$ 破缺）。
- CeCoIn$_5$/CeRhIn$_5$ 相关：Nature Communications 2023（CeRhIn$_5$ 电荷退局域化交叉）、PNAS 2017（Kondo lattice 量子临界标度）、以及若干 DFT+DMFT Doniach 相图重建工作（c-f 杂化随压力增大、CeRhIn$_5$ 压力相图 $T_N$ 先升后降、$\gtrsim1.15\,$GPa 出现超导）。
- 3d 杂质在 Cu 中的 $T_K$：多篇 DMFT/NRG 对比综述（含 Co in Cu vs Cu$_2$N/Cu(100) 表面配对 $T_K\approx2.6\,$K 的对比数据）。
