# 同位素系数：经原文核实的数据（curated）

来源：`surveys/P1_round2_pdf_extraction_2026-09-08.md`（researcher 读本地 PDF，页码 + 原句）。本表只收**原文直接给出或可由原文表格直接算出**的量；`hydrides_experimental.csv` 为一轮 WebSearch 摘要级原始表，其中 LaD10 (180 K / 0.46)、YD9 行已被本表推翻，请勿再用。

| system | isotope pair | P_H / P_D (GPa) | Tc_H / Tc_D (K) | alpha_exp | alpha_theory (method, μ* convention) | source | note |
|---|---|---|---|---|---|---|---|
| H3S | H/D | >170 区间 | 203 / ~150 | **0.3**（作者：P>170 GPa 区间比较） | 0.35（Errea 2015 SSCHA + iso-ME，μ*=0.16 fixed）；谐波 0.45 | Drozdov 2015 p.3；Errea 2015 Table I | 两点比 203/150 给 0.44，作者不用两点比 |
| LaH10 (fcc) | H/D | ~170 / ? | 250 / 168（正文）或 160（Ext.Data Fig.5） | **0.57**（作者，fcc 相） | 0.43（Errea 2020 SSCHA + aniso-ME，μ*=0.1 fixed，163/159 GPa）；AD 0.42；SCDFT 0.40 | Drozdov 2019 正文；Errea 2020 Ext.Data Table I | 文献中流传的 0.46–0.47 / 180 K **不在原文**；LaD10 压力与相待抽取 |
| LaH10 (P4/nmm) | H/D | ? / 142 | 215 / 140 (LaD10.6) | **0.62**（作者） | — | Drozdov 2019 | 竞争相，化学计量 LaD10.6 |
| YH6 | H/D | 166 / 172 | 224 / 170 | **0.40**（作者） | 谐波 iso-ME μ*(6Ry)=0.10–0.15：Tc 261–272；SSCHA：236–247；对上 Tc 与 α 需 μ*(6Ry)=0.19–0.22；SCDFT Tc=160，μ=0.187。**无 YD6 SSCHA** | Troyan 2021 p.5, p.7–8 | μ*(6 Ry) 约定问题见推导稿 §4.3 |
| YH6 (Kong) | H/D | 237 / 170 | 227 / 170 | （作者未算；压差 67 GPa） | — | Kong 2021 p.4 | 不可用于 α |
| YH9 | — | — | 243 / **无 D 数据** | **无** | — | Kong 2021 全文 | 唯一"YD9"来源为已撤稿 Snider |
| PdH | H/D | 常压 | 9 / 11（ref 7）；8 / 10（ref 8） | **−0.29 / −0.32** | 谐波 0.47；SSCHA **−0.38**（μ*=0.085 fixed） | Errea 2013 Table I | 最干净的反常同位素案例 |
| PdH | H/T | 常压 | — | — | 谐波 0.41；SSCHA −0.29 | Errea 2013 Table I | 无实验 α |

## SSCHA 两套谱的 (λ, ω_log, Tc)（用于两通道分解，推导稿 §3.2d）

| system | P (GPa) | λ_har | ω_log,har (meV) | Tc_har (K) | λ_anh | ω_log,anh (meV) | Tc_anh (K) | μ* | source |
|---|---|---|---|---|---|---|---|---|---|
| H3S | 200 | 2.64 | 90.4 | 250 | 1.84 | 92.86 | 194.0 | 0.16 fixed | Errea 2015 Table I |
| D3S | 200 | 2.64 | 68.5 | 183 | 1.87 | 73.3 | 152.0 | 0.16 fixed | 同上 |
| LaH10 | 163 | — | — | — | 2.67 | 96.4 | 242.8 (aniso-ME) | 0.1 | Errea 2020 Ext.Data Table I |
| LaD10 | 159 | — | — | — | 3.14 | 63.5 | 180.4 (aniso-ME) | 0.1 | 同上 |
| LaH10 | 129 / 214 / 264 | — | — | — | 3.62 / 2.06 / 1.73 | 76.4 / 115.5 / 126.6 | 255.3 / 237.9 / 216.9 | 0.1 | 同上 |
| LaD10 | 210 / 260 | — | — | — | 2.21 / 1.80 | 81.7 / 92.2 | 172.9 / 157.9 | 0.1 | 同上 |
| PdH | 0 | 1.55 | 25.42 | 47 | 0.40 | 50.21 | 5.0 | 0.085 fixed | Errea 2013 Table I |
| PdD | 0 | 1.55 | 18.60 | 34 | 0.46 | 37.69 | 6.5 | 0.085 fixed | 同上 |
| PdT | 0 | 1.55 | 15.50 | 30 | 0.48 | 31.86 | 6.9 | 0.085 fixed | 同上 |
| YH6 | 165 | 2.24 | 80.06 | 261–272 | 1.71 | 114.87 | 236–247 | μ*(6Ry)=0.10–0.15 | Troyan 2021 p.7–8 |

## 待抽取
- LaD10 168 K 转变对应的压力与相鉴定（Drozdov 2019 Fig. 3 / Ext. Data）
- Villa-Cortés 2022 结果 section：常 DOS 版 α(P)、el-ph/el-el 分量逐点值
- Troyan/Akashi 的 μ*(6 Ry) 定义：库仑核是否作用到 6 Ry 内全部 Matsubara 频率
