# D1 · 超越 Migdal-Eliashberg 的电声超导 Tc 计算

层级：L2 → L4　状态：调研中（2026-09-07）

## 范围

标准路径 DFPT → α²F(ω) → Allen-Dynes / 各向同性 ME 的四个前提分别被破坏时怎么办：

| 前提 | 破坏场景 | 修正方向 |
|---|---|---|
| Migdal 定理（ω_ph/E_F ≪ 1） | 氢化物 ω_ph ~ 0.2 eV | 顶角修正、非绝热理论（Pietronero–Grimaldi；Schrodi–Oppeneer–Aperis 2020） |
| N(ε) ≈ N(0) | H3S 的 van Hove 奇点落在声子能量尺度内 | full-bandwidth / energy-dependent DOS ME（Sano et al. 2016） |
| 费米面各向同性 | MgB2 双能隙；氢化物各向异性 | 各向异性 ME（EPW） |
| Allen-Dynes 拟合区（λ ≲ 1.5） | 氢化物 λ ~ 2–3 | 数值解 ME；SISSO/ML 修正公式（Xie et al. 2022） |

附：电声矩阵元 g 本身的精度（泛函依赖、GW 修正、SOC）。

## 核心问题

对 λ~2–3、ω_ph 高、DOS 剧变的体系，**哪一层近似贡献了最大的 Tc 误差**？哪些修正已有成熟工具，哪些还是纸面理论？

## 现状摘要

（待 surveys/D1 归档后填写）

## 可切入的具体问题

（待填）

## 关键文献

见 `library/index/D1.json`；PDF 在 `library/pdf/D1/`。

## 与其他方向的交叉

- D2：μ* 与截止频率绑定，数值解 ME 时必须一起处理
- D3：非谐修正后的 α²F 才是各向异性/顶角修正的正确输入
- D4：Allen-Dynes 的系统误差直接决定高通量排序的可靠性
