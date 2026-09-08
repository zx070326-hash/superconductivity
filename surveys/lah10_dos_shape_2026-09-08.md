# LaH10 费米面附近 DOS 形状（文献读图）

调研 agent：lah10-dos-shape（sonnet / researcher）　日期：2026-09-08。原始产出摘录。WebSearch 配额耗尽，全用 WebFetch + 本地 PDF 读图。

## 来源
1. **Errea et al. Nature 2020 / arXiv:1907.11916，Ext. Data Fig. 10**（本地 PDF p.19）：Fm-3m 100/150/200/250/300 GPa 五条 total DOS 叠加，横轴 ±2 eV。紧贴 E_F 有窄高尖峰（~1.3–1.4 states/eV/f.u.），两侧 ±0.3–0.5 eV 浅谷（~0.75–0.8）。**分辨率不足以判定尖峰在 E_F 哪一侧**。图注：畸变结构在 150 GPa 使 N(E_F) 降 23%（C2）/ 11%（R-3m）。
2. **Liu, Wang, Yi, Kim, Kim, Cho, PRB 99, 140501(R) (2019) / arXiv:1811.08548，300 GPa**：Dirac 节线态在 L 点分裂出双峰 vHs，ΔE ≈ 90 meV；Fig. 2(b) 插图 + SM Fig. S5（横轴 −0.4～+0.2 eV，分辨率较好）：**峰在 E_F 下方 −0.10～−0.15 eV，E_F 在峰的高能侧下坡**——与 H3S 拓扑相同。读图：N̄(−100) ≈ 1.12、N̄(+100) ≈ 0.94、N̄(−200) ≈ 1.09、N̄(+200) ≈ 0.88（像素级，中低置信度）。**比 H3S 平坦得多**（H3S N̄(+100) = 0.76）。
3. Kruglov et al. arXiv:1810.01113：只有 N(E_F) 单点值（Fm-3m 170 GPa 11.2、200 GPa 10.6 states/f.u./Ry），无形状。
4. Quan, Ghosh, Pickett arXiv:1906.02695：文字说 vHs "at E_F"（与 1/2 的"偏移 90–150 meV"有张力），图不可读。

## 未找到
PNAS 2017、PRL 2017、PRB 98 100102(R) 无 arXiv；无任何数值 DOS 文件（Materials Cloud/Zenodo/SI）。

## 综合
三个来源同意：Fm-3m LaH10 的 E_F 附近有 vHs，E_F 大概率在其高能侧下坡（与 H3S 同向）；峰的陡峭度两个来源差 30–40%，且精细图都在 300 GPa。**若定性成立，I(ω₀) < 1 且比 H3S 更接近 1**——顶角修正与 H3S 同号（压制）但更弱。

## 主会话补算（粗模型，`p2_dos_deficit` 同一泛函）
用来源 2 的读图值构造 N̄(ε)（远处取 1）：LaH10 P_s/λ ≈ −0.03～−0.04（ω₀ = 96–110 meV）；同一粗模型下 H3S ≈ −0.07。即 LaH10 的静态顶角压制约为 H3S 的一半，**符号相同**。
