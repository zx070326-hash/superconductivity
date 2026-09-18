# μ*>0 Einstein 极值性：独立复核推翻「ϖ≲0.33 诚实边界」，并给出一条无条件路线

**日期**：2026-09-16 · **分支** p4-conj2（未提交） · **作者**：team-lead 亲自从零重算
**触发**：`thmE-mustar` 的障碍图 `thmE_mustar_obstruction_map_2026-09-16.md`
**脚本**（全部本文件新建，未提交）：
`_scratch_ir/probe_mustar_concavity_Nc.py`、`probe_E4_delta_k_mustar.py`、
`probe_extremality_2atom_exhaustive.py`、`probe_convexity_bound.py`
算子在每个脚本里**从零重建**，不 import `spectral_efficiency.py`（后者 import 即跑全表）。

---

## 0. 结论摘要

| # | 陈述 | 级别 |
|---|------|------|
| 1 | `spectral_efficiency.py` 的 Coulomb 窗口缺 N_c 截断（Defect A）——确认 | 源码事实 |
| 2 | 修正方向**是**确定的：用物理 N_c 后凸性窗口**变宽变强**，不是变窄 | numerical-evidence |
| 3 | 5 点差分不是噪声（我怀疑的第三个缺陷**不成立**，撤回） | numerical-evidence |
| 4 | **「ϖ≲0.30–0.33 极值性可能失效」是伪命题**：定理 E 的**结论**在该窗口内成立 | numerical-evidence |
| 5 | 一般形式的核单调性（E4/引理 M）在 μ*>0 下**是假的**，不只是未证 | numerical-evidence |
| 6 | (CVX)：𝔨(P) ≤ ⟨h_μ⟩_P，**无条件、对一切 μ*、一切 ϖ** | analytic（一行） |
| 7 | 我上一条消息提的「端点弦下界」**已存在于仓库** `P4_tc_sandwich.md` 引理 A | 订正 |

---

## 1. Defect A：确认，且方向与 agent 所说相反

`spectral_efficiency.py:358`（与 `:402`）：

```python
u = 1.0 / np.sqrt(2.0 * np.arange(Nuse) + 1.0)
M = H - 2 * mu_tilde * np.outer(u, u)
```

全库 grep `Nc|N_c` 在该文件**零命中** ⇒ 纯粹漏写，没有被移到调用侧。Coulomb 窗口 = 数值截断
`Nuse = max(200, 5ϖ)`，而物理 N_c = ω_c/2πT ≈ 6–45。

‖u‖² = Σ_{n<N}1/(2n+1) = ½ln N + ½(γ_E+ln4) + O(N⁻²)：
N=200 → 3.6309、N_c=45 → 2.8851、N_c=6 → 1.8782。
即微扰强度被**放大 1.26–1.93 倍**（agent 报的「1.7–1.9×」只覆盖 N_c 小端）。

agent 说「修正方向不明显」。我算了，**方向完全确定**。H_μ''(x)，μ̃=0.05，dx=x/50：

| x | ϖ | N_c=6 | N_c=15 | N_c=45 | N_c=200（原脚本） |
|---|---|-------|--------|--------|------------------|
| 0.0100 | 0.100 | +2.878 | +2.513 | +1.717 | +0.785 |
| 0.0900 | 0.300 | +4.213 | +1.768 | +0.719 | +0.0042 |
| 0.1089 | 0.330 | +3.707 | +1.442 | +0.499 | −0.125 |
| 0.1600 | 0.400 | +1.699 | +0.610 | +0.024 | −0.379 |
| 0.2500 | 0.500 | −0.208 | −0.242 | −0.416 | −0.578 |
| 1.0000 | 1.000 | −0.302 | −0.291 | −0.284 | −0.280 |

N_c=200 一列复现了源笔记的 ϖ≈0.32 过零点 ✓（我的重建与被审对象是同一个算子）。
换成物理 N_c：过零点推到 **ϖ≈0.47（N_c=6）/ 0.45（N_c=15）/ 0.40（N_c=45）**，且凸性强度大 5–10 倍。
**bug 一直在掩盖 H″ 失效的严重程度。**

## 2. 我自己怀疑的第三个缺陷不成立（撤回）

担心 5 点模板在小 x 处除以 dx²≈10⁻⁶ 会淹没在 float64 相消噪声里。实测 dx = x/200、x/50、x/10
三档（跨 20 倍）：H″ 稳定到 4–5 位（x=0.1089 处 6.0925e-3 / 6.0925e-3 / 6.0945e-3）。
**源笔记的 H″ 数字是真的。** 记下来：这是一次我的误判。

## 3. Defect B：确认，并换成正确的检验

源笔记用 ϖ=1 与 √17 做「两原子 vs Einstein」，两点相距极远，与局域凹性无关。
正确做法是**合并检验**：x_± = x(1∓d)，权重各 ½，⟨x⟩=x 精确，Jensen 间隙 ≈ ½H″(x)x²d²。

gap/d² 在 d = 0.4→0.05 上收敛到 3–4 位，且**在所有 (x, μ̃, N_c) 处为负**，包括窗口内部：

| x | μ̃=0.05,N_c=6 | μ̃=0.05,N_c=45 | μ̃=0.13,N_c=6 | μ̃=0（对照） |
|---|---|---|---|---|
| 0.0900 | −4.353e-3 | −5.237e-3 | −2.480e-3 | −7.494e-3 |
| 0.1089 | −6.904e-3 | −7.603e-3 | −3.702e-3 | −1.043e-2 |

**H″>0 与「Einstein 仍然胜出」同时成立** ⇒ H″ 不是定理 E 的判据（见 §5）。

## 4. 定理 E 在 μ*>0 的**穷举**检验：结论成立

**穷举性论证**：M(P) = 𝒪[F_P] − 2μ̃|u_c⟩⟨u_c| 对 P **仿射**（𝒪 线性、Coulomb 项与 P 无关），
故 𝔨(P)=λ₁(M(P)) 是仿射泛函的上确界 ⇒ **对 P 凸**。凸函数在紧凸集
{P≥0 : 质量=1, 均值=x̄, supp⊂[0,X]} 上的**最大值取在端点**（Bauer 极大原理），而两条线性约束
割出的集合的端点是**至多两点支撑**的测度。⇒ 扫遍匹配 (质量, 均值) 的两原子测度即**穷举**。

> **归属**：凸性这一条**仓库已有**——`spectral_efficiency.py:368–370` 原文
> 「lam1(M) is convex in P (sup of affine). [analytic]」。仓库没做的是由它推出的两个推论：
> Bauer 归约（本节）与 (CVX) 界（§6）。不要把凸性本身算作新东西。

`probe_extremality_2atom_exhaustive.py`，N=400，w∈[0.02,0.98]×25 点 × x₁ 几何网格 40 点：

**best excess := max_P 𝔨(P) − 𝔨(δ_x̄)，全部为负**，7 组 (μ̃,N_c) × 7 个 x̄ = 49 格，零反例：

| μ̃ / N_c | x̄=0.01 | 0.09 | 0.1089 | 0.25 | 1.0 | 9.0 | 100 |
|---|---|---|---|---|---|---|---|
| 0.05 / 6 | −5.6e-10 | −8.1e-8 | −1.3e-7 | −6.9e-7 | −3.0e-6 | −4.8e-6 | −2.5e-7 |
| 0.13 / 6 | −5.2e-10 | −4.6e-8 | −6.9e-8 | −4.0e-7 | −2.9e-6 | −4.7e-6 | −2.6e-7 |
| 0.26 / 6 | −5.1e-10 | −3.8e-8 | −5.5e-8 | −2.5e-7 | −1.9e-6 | −4.5e-6 | −2.6e-7 |
| 0.13 / 45 | −1.1e-9 | −8.2e-8 | −1.2e-7 | −4.9e-7 | −2.6e-6 | −4.6e-6 | −2.4e-7 |

最大值点恒为 x₁≈x₂≈x̄（退化到 Einstein 点本身），excess 随原子合并 →0⁻。这正是
δ 为严格极大点的签名。**ϖ=0.30/0.33（x̄=0.09/0.1089）在窗口正中，照样没有反例。**

**判决**：`pairing_duality_theorem_2026-09-16.md` §5.2–5.3 的「诚实边界」应当**替换而非收紧**。
它描述的是一个**充分条件**的失效，不是定理的失效。论文里不能写成「μ*>0 时极值性在
ϖ≲0.33 可能不成立」。

*网格保留*：w/x₁ 网格是有限的；但 §3 的合并检验独立地在 d→0 上确认了退化点邻域的负号，
两者合起来才构成证据。仍为 `numerical-evidence`。

## 5. 真正的障碍：一般核单调性在 μ*>0 下是**假**的

定理 E 的链 E1+E4 中，E1（逐通道凹性）是 **μ\* 无关的恒等事实**：
f_x(k)=x/(k²+x)，∂²f/∂x² = −2k²/(k²+x)³ < 0（k≥1），f_x(0)=1 与 x 无关。Coulomb 项根本不碰 f。
所以 μ*>0 的全部困难压在 E4（引理 M）上。由变分原理，E4 只需对 M 的顶向量 v 有

  **(\*) ⟨v, 𝒪[δ_k] v⟩ ≥ 0，∀k≥1**（δ_k = 第 k 通道的指示；锥条件是它的**充分**条件）

`probe_E4_delta_k_mustar.py` 直接在**真实 Perron 向量**上测 (\*)：

| μ̃ / N_c | x̄=0.01 | 0.09 | 0.1089 | 0.25 | 1.0 | 9.0 | 100 | 变号数 |
|---|---|---|---|---|---|---|---|---|
| 0 / — | +2.2e-5 | +2.4e-5 | +2.4e-5 | +2.8e-5 | +4.4e-5 | +2.0e-4 | +1.7e-3 | 0 |
| 0.05 / 6 | **−0.443** | **−0.447** | **−0.419** | **−0.179** | +5.2e-5 | +2.3e-4 | +1.9e-3 | 1–2 |
| 0.13 / 45 | **−0.164** | **−0.159** | **−0.157** | **−0.142** | **−0.089** | **−0.027** | **−5.2e-3** | 1–2 |

μ*=0 时全正（锥成立）；μ*>0 时在 k=4–6 处负到 **−0.44，与 λ₁ 本身同量级**（λ₁≈0.02–0.5）。
由 Hellmann–Feynman，这意味着抬高单个通道 f(4) 会**降低** λ₁。

**所以这不是「锥条件未证」——(\*) 被决定性地证伪。** agent 关于「锥不变性路线不可救」
（μ_Metz ~ ϖ̄²/N_c³，低于物理 μ̃ 2–4 个数量级）的结论方向正确，但理由要改写：
问题不在锥的尺寸，而在一般 δ 的单调性本身就不真。

**而结论（§4）仍然成立** ⇒ 证明必须利用「物理 δ 张成的锥远小于 {δ≥0}」——
δ = F_P − f_x̄ 是 Stieltjes 型核的 Jensen 亏损，不是任意非负序列。这是 μ*>0 定理 E 的
**新的、准确的开放问题陈述**，取代原先的「恢复锥不变性」。

## 6. (CVX)：一条无条件、对一切 μ* 成立的上界

由 §4 的凸性，对 P = Σᵢwᵢδ_{xᵢ}（Σwᵢ=1）直接得

  **𝔨(P) ≤ Σᵢ wᵢ 𝔨(δ_{xᵢ}) = ⟨h_μ⟩_P**  ……（CVX）

其中 h_μ(x) := λ₁(𝒪[f_x] − 2μ̃|u_c⟩⟨u_c|) 用**同一个** μ̃ 与 N_c。

**严格性：analytic。** 证明一行，不需要锥、不需要 Perron 结构、不需要核单调性，
**对一切 μ*≥0、一切 ϖ、一切 N_c 成立**——正是现在缺的那块。

数值验证（`probe_convexity_bound.py`，36 组随机多原子谱，2–5 峰，x∈[10⁻³,400]）：
零违反，最小 slack +4.9e-5。并且在 30/36 组里 (CVX) 比定理 E 的界**更紧**；
剩下 6 组恰好落在 H_μ 凸的小 ϖ 区（此时 ⟨h⟩>h(x̄)，与 §1 的表一致——两者互为印证）。
实用界是 **min((CVX), 定理 E)**。

**物理读法**：任何谱的配对本征值不超过其各个 Einstein 模**单独作用**时的谱加权平均
⇒ **配对通道里声子模之间没有相长干涉**。这是一句可引用的物理陈述，且不依赖 μ*。

**与定理 E 的分工**：(CVX) 更紧但需要完整的 P；定理 E 压缩到两个矩 (λ,⟨ω²⟩) 但需要引理 C。
链条 𝔨(P) ≤ ⟨h_μ⟩_P ≤ h_μ(√x̄) 的第二步就是 H_μ 凹性——**所以源笔记算 H″ 并没有算错对象，
只是它是这条替代路线的判据，不是 E1+E4 那条路线的判据。** 我在 §3 的措辞据此修正。

## 7. 订正：端点弦下界不是新的

我上一条消息把「T_c(P) ≥ 端点两原子测度的 T_c」当作可以免费拿到的新下界。
**仓库已有**：`P4_tc_sandwich.md` §1 引理 A 下界-2，用 ψ(x) ≥ (x/x_max)ψ(x_max)（原点到 x_max 的弦），
给出 𝔨 ≥ (⟨ω²⟩/ω_max²)·h(ϖ_max)；下界-1 给 𝔨 ≥ h(ϖ_min)。
我的版本（[x_min,x_max] 端点弦）**严格更紧**，因为该弦在 [x_min,x_max] 上压过原点弦，
但它是对已有引理 A 的**收紧**，不是新定理，而且同样**依赖引理 M ⇒ 只在 μ\*=0 有效**。

μ*>0 下现成可用的下界只有平凡的变分式：取任一单位 v₀（例如 Einstein 顶向量）
𝔨(P) ≥ ⟨v₀,𝒪[F_P]v₀⟩ − 2μ̃⟨v₀,P_c v₀⟩。这条 `analytic`、μ*-稳健，但松。

## 8. 未决：γ-model 脚本的 μ* 归一化差一个 2

对称化把负频镜像折进来时，声子项得到 f(|n−m|)+f(n+m+1)（两项），而 Coulomb 项得到
**2μ̃ u_n u_m**（一个因子 2）。`spectral_efficiency.py` 用的是 2；
`_scratch_ir/gamma_dk_coulomb.py` 用的是 `H[n][m] -= mu*u_n*u_m`（因子 1）。
若两者意在同一个 μ*，则 γ-model 表里「μ\*=0.13 ⇒ Tc 压制 −15.2%」实际对应 μ̃=0.065，
**物理 μ̃=0.13 处的压制要更大**。已向 `thmE-mustar` 索取因子 2 的出处行号。
**在此结清之前，γ-model 的 μ* 数值标签不得进论文表。**

---

## 9. KAY II / III 源级核验（2026-09-17，我自己抽取 PDF 文本，非转引 agent）

`lit-mu-qcp` 的文献判决自评 L1.5 且承认 WebFetch 误解析。我把两篇 PDF 下到
`_scratch_ir/KAY2.pdf`、`KAY3.pdf`，用 `pdftotext -layout` 抽文本后逐条 grep。**源级事实**：

1. **KAY II Conjecture 2（eq.31）原文**：「There is a critical temperature Tc(λ,P) > 0 which for
   all P ∈ 𝒫 and λ > 0 is bounded above by Tc(λ,P) < T̄c(λ,P), with T̄c(λ,P) := (1/2π)√(𝔤(2) λ⟨ω²⟩)」。
   = 仓库的定理 1。**不含极值性，不含 μ\*。**
2. **KAY II Remark 3 原文**：「The asymptotic behavior … was first stated in [AD], see their eq.(22).
   Apparently they did not notice that this is presumably an upper bound on Tc for all λ > 0.」
   ⇒ AD eq.(22) 只给**渐近式**；「对一切 λ 成立的上界」才是 KAY 的新内容。
3. **`grep -ic "coulomb|pseudopotential"`：KAY II = 0，KAY III = 0。**
   整个 KAY 系列 μ*≡0。⇒ 本仓库 μ*>0 的一切结果（(CVX)、§4 穷举、Anderson–Morel）**不在其覆盖内**。
4. **KAY II 无 Einstein 极值性陈述**：`einstein` 仅 3 处命中——定义非色散模型、
   Outlook 说 part III 将**特化**到该模型、参考文献条目。无「δ 测度最大化 Tc」。
   KAY III Conjecture 1 是 Conjecture 2 的 Einstein 特化，同样无极值性。
5. **KAY II §9.1 对 AD 逐条审计**（eq.19/21/22/23），既给信用（Remark 3）也直接纠错
   （「We also found that it is not correct when [AD] claim…」）——**通篇未提 AD 有任何极值性主张**。

**仍未核**：~~AD 1975 原文~~ → **已核，见 §9b。**

## 9b. AD 1975 原文核到了（2026-09-17，L1，我自己读的 PDF）

`lit-mu-qcp` 给出了可用路线 `harvest.aps.org/v2/journals/articles/10.1103/<DOI>/fulltext`。
我自己下载 `_scratch_ir/AD1975.pdf`（18 页）与 `CARB.pdf`（RMP 62,1027），`pdftotext -layout` 抽文本后逐段读。
**agent 的三条承重断言全部在源头证实**，且原文比转述更锋利。

**(i) AD p.909 正文，逐字**（双栏 OCR，左栏重组）：

> 「…in the sense that T_c is maximized under the conditions of **fixed λ, μ\*, and ⟨ω²⟩**. … At smaller values of λ the shape does play a role.
> **It is clear from inspection of Eq. (15) that the eigenvalue ρ (and thus T_c) is maximized by making the parameters λ(n) as large as possible.**
> … we can then ask which shape of g(ω) (for fixed second moment ⟨ω²⟩) maximizes λ(n)/λ.
> **It is proved in Appendix A that the Einstein spectrum optimizes the integral (29) and thus is the preferred phonon spectrum.**」

**(ii) AD Appendix A，逐字**：

> 「given λ, μ\*, and ⟨ω^n⟩, for an arbitrary value n ≥ 2, the maximum T_c consistent with these parameters is obtained
> **when (and only when)** α²F is an Einstein spectrum.」
> 而实际证明只有一句：「The integrating factor **x/(x+a) has a negative second derivative** for all positive values of x.
> Thus (A4) is immediately proved using theorem (A10).」

**(A4) 就是我们的 E1。** f(a)=∫dx g(x)·x/(x+a)，取 a=k² 即 f_x(k)=x/(k²+x) 的逐通道凹性——同一行论证。
⇒ **AD 附录 A 证的恰好是 E1，一字不多。** E1 → 本征值（= 我们的 E4 / 引理 M）是
「clear from inspection」的**断言**，全文无证明。

**(iii) 判决（这是本项目的关键定位）**：
- **陈述**不新：AD 已在**完全相同的约束集**（且含 μ\*、且含任意矩 n≥2、且带 iff）写下。
- **E1 的证明**也不新：AD 附录 A 同一行凹性。
- **新的是 E1 ⇒ 本征值那一步**（引理 M）。μ\*=0 时那是我们的贡献。
- **μ\*>0 时 AD 的理由是错的，不只是未证**：§5 实测 ⟨v,𝒪[δ_k]v⟩ = **−0.44**（μ̃=0.05,N_c=6），
  抬高单个 λ(n) **降低** ρ。而 AD 恰恰把 μ\* 列进了他们的约束集。
  ⇒ 「AD 在他们自己声称的 μ\*≠0 区域里，用了一条可证伪的理由」。结论仍真（§4 穷举），路线必须换。

**(iv) Carbotte RMP 1990 §II.H p.1045 逐字**（我读到的）：
> 「We can prove that (2.48) represents a **local maximum** for any shape of α²F(ω) with a **fixed A**. …
> Thus we have maximized T_c/A [with the] optimum spectrum (2.49) for arbitrary A. **It follows** [that] the inequality
> T_c ≤ A c(μ\*) (2.50) holds. … This [was] first established by **Leavens (1975)** using somewhat different arguments.」
- 约束是**固定面积 A**，**不是**固定 λ+⟨ω²⟩ ⇒ Carbotte (2.50) **不是**我们的定理 E，是另一个约束流形。
- 同样的逻辑缺口：只证 local maximum，然后「It follows」跳到全局。p.1096 他自己承认：
  「Our experience with realistic spectra would lead us to believe that this is, in fact, an absolute maximum.」
⇒ **field-level 陈述成立**：Leavens(1975) → AD(1975) → Carbotte(1990)，δ 谱极值性一直是
「局部变分已证 + 全局靠经验/数值」。任何约束流形上的**全局**陈述都未被证明。

**(v) C_∞ 的数值先例**（我独立算过）：Carbotte eq.(2.29) p.1042「T_c = 0.2584 … = 0.183√λ ω_E」。
δ 谱 λ=2A/ω_E ⇒ 系数 0.2584/√2 = 0.1827164，其舍入区间 [0.1826810, 0.1827517]
**包含** KAY 的 C_∞ = 0.1827262477。⇒ C_∞ 在 1990 年已被数值知道到 4 位。引入 C_∞ 时应引 (2.29) 为数值先例。

**(vi) 未核、需软化**：Carbotte 把渐近式**并列**归给 Owen & Scalapino, Physica 55, 691 (1971)。
该文正文不可得（DOI 转 linkinghub 桩，pre-split Physica 无免费扫描）。
⇒ 论文里不得写「AD 1975 是最早出处」，须写「the earliest explicit statement we are aware of」。

## 10. 下一步（按价值排序）

1. **把 (CVX) 写成定理**并查新颖性（一行证明的东西必须查）。它单独就填上了「μ\*>0 无定理」的洞。
2. **Bauer 归约**：把定理 E 从无穷维测度问题归约到 **2 参数** (w,x₁) 不等式——μ\* 无关、完全严格。
   这是目前看得见的唯一一条能把 μ*>0 定理 E 做成 `analytic` 的路。
3. 重写 `pairing_duality_theorem_2026-09-16.md` §5.2–5.4 的「诚实边界」。
4. 结清 §8 的因子 2。
5. 端点弦收紧写进 `P4_tc_sandwich.md` 引理 A（μ*=0）。
