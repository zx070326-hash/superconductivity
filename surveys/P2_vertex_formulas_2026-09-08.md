# P2 支撑调研：Pietronero–Grimaldi 顶角函数与非绝热 Eliashberg 方程的可编码形式

调研 agent：p2-vertex-formulas（sonnet / researcher）　日期：2026-09-08。原始产出，未手改（两条回传合并）。

## 前提性发现
PRB 52,10516 (1995)、PRB 52,10530 (1995)、PRL 75,1158 (1995)、Cappelluti–Pietronero PRB 53,932 (1996) 经多次 arXiv API 查询确认**均不在 arXiv**。以下公式来自能在 arXiv 找到的后续论文（journal-ref 已核实），是转录/再推导而非原始出处逐字复刻。

## 1. 顶角函数 P_V、P_C
最佳替代：**cond-mat/9806012 = Grimaldi, Pietronero, Scattoni, Eur. Phys. J. B 10, 247 (1999)**，"Physical Origin of the Electron-Phonon Vertex Correction"。不用 P_V/P_C 命名，分解为 **P_pol**（晶格极化）+ **P_ex**（电子-空穴激发）：

```
P = P_pol + P_ex
P_pol(k+q,k;iω_k+iω_q,iω_k) = (ω₀/2) Σ_k' g²_{k−k'} × { [f(ε_k')+n(ω₀)]/[(ε_k'−ω₀−iω_k)(ε_k'+q−ω₀−i(ω_k+ω_q))]
                                                        + [1+n(ω₀)−f(ε_k')]/[(ε_k'+ω₀−iω_k)(ε_k'+q+ω₀−i(ω_k+ω_q))] }
P_ex 含 [f(ε_k'+q)−f(ε_k')]/(ε_k'+q−ε_k'−iω_q)
```
静态极限（先 ω_q→0 再 q→0）：P_ex → −λ；一般 **P_s = −ω₀/(ω₀+E_F)**（负）。动态极限（先 q→0）：P_ex → 0；**P_d = +E_F/(ω₀+E_F)**（正）。反绝热极限 ω₀→∞：P → P_ex。变量 Q = q/(2k_F)。
交叉验证：Schrödi 2020 Sec. IIIC 转述同一结论（静态负、动态正）。**警告**：公式经 WebFetch 对 arXiv HTML 自动提取，未逐字核对 PDF。

## 2. 顶角修正如何进入 Z、φ
**cond-mat/9710159 = Grimaldi, Cappelluti, Pietronero, Europhys. Lett. 42, 667 (1998)**：
```
Z(iω_n) = 1 + (πT/ω_n) Σ_m λ_Z(iω_n,iω_m,Q_c) · [ω_m/√(ω_m²+Δ_m²Z_m²)] · (2/π)arctan(E/(2|ω_m|Z(iω_m)))
λ_Z(iω_n,iω_m,Q_c) = λ[1 + λ P_V(iω_n,iω_m,Q_c)]
静态极限：Z₀ ≈ 1 + λ − λ²(π/4)(ω₀/(Q_c² E_F))
```
Q_c ≡ q_c/(2k_F)；E = 电子带宽。**φ/Δ 方程的显式系数没找到**。cond-mat/0007007（Botti, Cappelluti, Grimaldi, Pietronero, IJMPB 14, 2976 (2000)，疑为 PRB 66, 054532 (2002) 的会议前身）明说 λ_Z 与 λ_Δ "contain different nonadiabatic corrections and have therefore different functional forms"，但未展开。

## 3. McMillan 型 Tc(λ_eff, μ*_eff) 闭式
**没找到**。综述 cond-mat/0006325 未给。建议主会话自己用 λ_Z 展开代入 AD 做零阶。

## 4. 判据 Q_c*(m)
**没找到显式公式或数值**。只有定性："More correlated the system smaller q_c"；小 Q_c → 顶角趋正增强 Tc，动量结构各向同性 → 趋负。**这是真正的开放推导任务。**

## 5. Cappelluti–Pietronero PRB 53,932 (vHs + 非绝热)
完全找不到内容。

## 6. EPL 42,667 的同位素公式
```
极化子极限（Holstein，对照）：m* ≈ m·exp(g²/ω₀²)；α_{m*} = −½(g²/ω₀²)
非绝热顶角机制：α_{m*} ≈ −½(λ²/Z₀)(π/4)(ω₀/(Q_c² E_F))
算例 λ₀=1，ω₀/E_F=0.05：Q_c=0.4 → α_{m*} = −0.29；Q_c=0.2 → −0.6
```
|α_{m*}| ∝ ω₀/E_F（线性）、∝ 1/Q_c²。注意是 m* 的同位素系数，不是 Tc 的。另给 vHs 型 DOS 模型 N(ε) = −N₀ ln|2ε/E|。

## 7a. Schrödi, Oppeneer, Aperis PRB 102, 024503（本地 PDF 全文核对）
完整自能（Eq. 4）：Migdal 项 + 一阶顶角（交叉）项 ∝ |g|⁴。Z、χ、φ 三通道方程（Eq. 7–9）各有不同的 3×3 顶角矩阵 P^(Z)、P^(χ)、P^(φ)（Eq. 11–13）。
**FSR 各向同性近似（Eq. 19、23–24）**：
```
Σ̂_m^(iso) = T² P Σ_{m'm''} λ_{m−m'} λ_{m'−m''} ρ̂₃ĝ_{m'}ρ̂₃ĝ_{m''}ρ̂₃ĝ_{m''−m'+m}ρ̂₃
P = (1/N₀³) Σ δ(ξ_k)δ(ξ_k')δ(ξ_k'')δ(ξ_{k''−k'+k})   ——纯相空间常数
Z_m^(iso) = 1 + (πT/ω_m)Σ λ_{m−m'}γ^{(ω)}_{m'} + Pπ³T²/ω_m Σ λ_{m−m'}λ_{m'−m''}[四项 γ 组合]
```
**作者原话**：该近似"severely underestimates the effect of the first vertex correction and is therefore not valid"——只能当代码骨架/定性。
**符号条件（原文）**：取决于维度 + α=ω_ph/E_F + 耦合强度；3D α≪1 可忽略，3D α~1 可正可负（依 λ），2D 即使 α≪1 也不可忽略。**自洽解出的顶角函数不是小 q 峰化的（Fig. 7），α 越大越正、越各向同性**——与 PGS "小 Q_c 增强"直觉相反，"Q_c" 在强 α 区可能失去意义。

## 7b. Mishra 2026 (arXiv:2608.14949) 的 A_λ、R_V
正文显式：λ^V = 4∫dω∫dω' α²F^V(ω,ω')/(ωω')；η = ω_log/E_F^eff（E_F^eff 由 vHs 附近拟合）。A_λ、R_V 正文无定义式；用 Table I 反推：**A_λ ≡ (λ_har − λ_anh)/λ_har**（H3S 0.35 ✓）；**R_V ≡ λ^V_anh/λ_anh**（YH6 0.14 ✓、YH9 0.14 ✓、PdH 0.08 ✓）。反推，非抄录。

## 编码建议
不建议用单一乘子 λ(n−m)[1+V] 同时改 Z 和 φ——文献表明两通道结构不同。可行路径：(1) 只做 Z 通道用 EPL 42,667 的 λ_Z，φ 置零，先核量级；(2) Schrödi FSR 各向同性式作骨架/自检，不作定量；(3) 用 cond-mat/9806012 的 P_pol/P_ex 一般式配模型核 V(Q) 数值积分——唯一能回答 Q_c*(m) 的路径，工作量相当于一次小型全带宽实现。

## 文献 JSON
```json
[
  {"title":"Isotope effect on m* in high-Tc materials due to the breakdown of Migdal's theorem","first_author":"C. Grimaldi","authors_short":"Grimaldi, Cappelluti, Pietronero","year":1998,"venue":"Europhys. Lett. 42, 667","arxiv_id":"cond-mat/9710159","doi":"10.1209/epl/i1998-00303-0","tier":"core","one_line":"λ_Z=λ[1+λP_V(Qc)]、Z0静态近似、α_m*闭式","why":"P2 items 2,6"},
  {"title":"Physical Origin of the Electron-Phonon Vertex Correction","first_author":"C. Grimaldi","authors_short":"Grimaldi, Pietronero, Scattoni","year":1999,"venue":"Eur. Phys. J. B 10, 247","arxiv_id":"cond-mat/9806012","doi":"10.1007/s100510050852","tier":"core","one_line":"顶角函数分解为P_pol+P_ex；静态极限负、动态极限正","why":"P2 item 1"},
  {"title":"Nonadiabatic theory of the superconducting state","first_author":"M. Botti","authors_short":"Botti, Cappelluti, Grimaldi, Pietronero","year":2000,"venue":"Int. J. Mod. Phys. B 14, 2976","arxiv_id":"cond-mat/0007007","doi":"10.1142/S0217979200003198","tier":"important","one_line":"疑似PRB66,054532(2002)会议前身；λ_Z≠λ_Δ结构式","why":"P2 item 2"},
  {"title":"A survey of nonadiabatic superconductivity in cuprates and fullerides","first_author":"E. Cappelluti","authors_short":"Cappelluti, Grimaldi, Pietronero, Strässler","year":2000,"venue":"arXiv / proceedings","arxiv_id":"cond-mat/0006325","doi":null,"tier":"reference","one_line":"综述，只有定性陈述，无Qc*(m)判据","why":"确认 items 3,4 文献缺失"}
]
```
