# BaFe₂As₂（Graser 真实 5-轨道参数）精细网格嵌套函数扫描

**日期**：2026-09-22
**任务**：Team C 优先级 1——用真实 Graser arXiv:1003.0133 5-轨道参数替换
`iron_based_tight_binding_params_2026-09-22.md`（Plan B）的手选简化 2-轨道模型，
重做嵌套度精细扫描并精确定位 Q_AFM。
**严格性级别**：`numerical-evidence`（float64，单套网格 + η 扫描做鲁棒性检验，
未做 certificate 级区间包围）。

**这篇笔记取代（不删除）** `iron_based_tight_binding_params_2026-09-22.md` /
`iron_based_nesting_survey_plan_2026-09-22.md` 中"BaFe₂As₂"那一行——那两篇明确
标注了"手选 W_scale，不是文献参数"；本篇用的是页图核验过的真实 Graser 参数
（见 [../literature/graser1003.0133_5orbital_TB_parameters_2026-09-22.md](../literature/graser1003.0133_5orbital_TB_parameters_2026-09-22.md)）。
**材料仍是 BaFe₂As₂(122族)，不是 FeSe(11族)**——拓扑结论按既有 caveat 可转移，
定量数值不可转移，见该文献笔记第 7 节。

---

## 结论摘要

1. **μ（6 e/Fe 填充）= 0.0109 eV**（N=128 网格），与独立实现（essenberger 的
   `_graser_find_mu_fermi_surface_2026-09-22.py`，N=120）重跑后给出 0.0127 eV
   **一致**（差异在网格分辨率量级内，非分歧——此前笔记误记的 −0.0656 eV 是
   参数订正前的陈旧缓存值，已排查清楚）。
2. **网格分辨率**：N=128×128，格距 2π/128 = 0.0156π，**优于** team-lead/Team C
   要求的 0.05π 目标（原 Plan B 用的是 q_grid=15、格距 ~0.07π 的粗网格）。
3. **X=(π,0)/Y=(0,π) 处确有增强**：ξ(X)/⟨ξ⟩ = 1.68×（η=0.03 eV），随 η 增大单调
   减弱，符合"展宽抹平嵌套信号"的物理预期。
4. **新发现，此前未报告**：全局非平凡极大值**不在** X/Y，而在对角方向
   q≈(0.51–0.62)π·(1,1)（依赖 η，见下表），其 ξ 值与 X/Y **相当或略高**
   （η=0.03 时 1.70× vs X 的 1.68×；η=0.01 时 2.95× vs X 的 2.98×，两者始终接近、
   谁更高随 η 摆动，不构成稳定排序）。**M=(π,π) 始终是三者中最弱的**，
   在 η≳0.1 eV 时增强降到 ≈1×（即与 BZ 平均无法区分）。
5. **一处方法论 bug 已发现并修复**：原始扫描脚本用半径 3 格
   （~0.047π）的圆盘遮蔽 q=0 自相关峰后取"全局 argmax"，该半径不足以盖住
   Γ 峰的平滑尾部——真实测出的"排除 Γ 后 argmax"其实仍是 Γ 峰的肩部（随遮蔽
   半径单调后退：68337→57609→48935→36698，直到半径 ≥0.125π 才稳定在真实
   的对角局域极大值 33224.86）。这是"平凡自相关峰"陷阱的一个新变种：
   `iron_based_tight_binding_params_2026-09-22.md` 记录过的旧 bug 是"没排除
   q=0"，这次是"排除半径不够大，仍混进了 Γ 峰的尾巴"——两者同属一类
   陷阱但不是同一个 bug，值得单独记录供复核。

## η 鲁棒性表（team-lead 若要复核，先看这张表）

| η (eV) | ξ(Γ) | ξ(X)/⟨ξ⟩ | ξ(M)/⟨ξ⟩ | 对角峰 q/π | 对角峰/⟨ξ⟩ |
|---|---|---|---|---|---|
| 0.01 | 202322.6 | 2.980 | 1.302 | 0.620 | 2.949 |
| 0.02 | 118575.9 | 2.030 | 1.193 | 0.510 | 2.044 |
| 0.03 | 85442.6 | 1.677 | 1.186 | 0.510 | 1.703 |
| 0.05 | 57266.2 | 1.371 | 1.149 | 0.510 | 1.392 |
| 0.08 | 41934.1 | 1.189 | 1.072 | 0.510 | 1.197 |
| 0.15 | 31162.2 | 1.055 | 0.988 | 0.500 | 1.054 |

**读法**：X 与对角峰在每个 η 下几乎并列第一（差距 ≤3%），M 恒为最弱、且随 η
增大最快跌到"无增强"。**Q_AFM 不能被本次扫描唯一定位在 X/Y=(π,0)**——
真实图像是"X/Y 与对角方向 q≈(0.5–0.6)π 附近同等强度地增强"，M 处反而最弱，
这与 Plan B 简化模型报告的"Q_AFM≈(0,π)"部分吻合（X/Y 确实增强）但补充了一个
Plan B 粗网格（q_grid=15，格距 0.07π）分辨率不足以看到的对角结构。

## 方法

- H(k)：`_scratch_ir/_graser1003_hamiltonian_2026-09-22.py`（essenberger 构建，
  已核验用的是订正后 Table I/II 数值，见文献笔记第 5.1 节 + 我发给
  essenberger 的订正确认）。
- μ：二分法解 6 e/Fe 填充（`find_mu`，`_graser1003_nesting_scan_2026-09-22.py`）。
- 嵌套函数：ξ(q) = Σ_k A(k)A(k+q)，A(k)=Σ_n Lorentzian(ε_n(k)−μ,η)（标准
  Lorentzian 展宽 joint spectral weight 代理），FFT 圆周自相关实现，
  O(N²logN)。
- η 鲁棒性：`_chk_nesting_eta_robustness_2026-09-22.py`（本次新增，已入库）。
- 自相关峰遮蔽半径鲁棒性：`_chk_nesting_localmax_scan_2026-09-22.py`
  （本次新增，已入库），确认对角峰位置对遮蔽半径 ≥0.125π 后稳定。

## 与费米面验收门的关系

X/Y 处的增强与"X/Y 确认为电子口袋"（见文献笔记第 8 节第 2 条，team-lead
2026-09-22 预算裁决）方向一致——电子口袋与 Γ 空穴袋之间的嵌套矢量正是
(π,0)/(0,π)。对角峰 q≈(0.5–0.6)π·(1,1) 连接的是哪两组费米面特征，本次未做
逐点费米面重建来回答，留作待办（不影响本文已交付的 ξ(q) 数值本身）。

## 待办 / 已知局限

1. 对角峰的物理归属（连接 Γ 的哪个口袋与哪个口袋）未做逐点核实。
2. 仅 float64 单精度网格 + 有限 η 扫描，未做 certificate 级区间验证——
   本文档标级 `numerical-evidence`，不得当作 `analytic`/`certificate` 引用。
3. 未做其余 4 种材料（LaFeAsO/LiFeAs/FeSe/FeTe0.5Se0.5）的真实参数版本——
   `iron_based_tight_binding_params_2026-09-22.md` 的简化模型结果对那 4 种
   材料仍是唯一交付物，标级不变（`numerical-evidence`，手选参数版本）。

## 相关笔记

- 参数来源+两仪器核验：[../literature/graser1003.0133_5orbital_TB_parameters_2026-09-22.md](../literature/graser1003.0133_5orbital_TB_parameters_2026-09-22.md)
- 被取代（保留存档）的简化模型：[iron_based_tight_binding_params_2026-09-22.md](iron_based_tight_binding_params_2026-09-22.md)、[iron_based_nesting_survey_plan_2026-09-22.md](iron_based_nesting_survey_plan_2026-09-22.md)、[iron_based_Q_AFM_universality_2026-09-22.md](iron_based_Q_AFM_universality_2026-09-22.md)
- 静态 χ₀ 规范依赖（不影响本文——本文只用能带能量，不构造 χ₀）：
  [../misc/essenberger_chi_sign_structure_2026-09-22.md](../misc/essenberger_chi_sign_structure_2026-09-22.md)
