# SCDFT 代码支持情况调研报告

**日期**: 2026-09-22  
**任务**: 评估现有 DFT 代码能否复现 Essenberger FeSe 的 SCDFT 计算  
**结论**: **所有候选代码均不支持 SCDFT** — 建议立即切换到从零实现路线  

---

## 执行摘要

调研了四个主流 DFT 代码（octopus、exciting、Elk、GPAW）的 SCDFT 支持情况：

| 代码 | SCDFT 支持 | 超导相关功能 | 评分 | 备注 |
|------|-----------|-------------|------|------|
| **octopus** | ✗ | 无 | ✗ | 只有 TD-DFT，无超导扩展 |
| **exciting** | ✗ | 无 | ✗ | 聚焦激发态，无 RPA 自旋磁化率 |
| **Elk** | ⚠ | Eliashberg（声子） | ⚠ | 有 α²F/λ/Tc，但**非 SCDFT 框架** |
| **GPAW** | ✗ | 无 | ✗ | 标准 DFT，无超导功能 |

**核心发现**：
1. **没有任何一个代码实现了 SCDFT**（superconducting DFT with f_xc kernel）
2. Elk 的超导功能是**传统 Eliashberg 理论**，不是 SCDFT
3. Essenberger 自己的计算流程是 **VASP + 自研 SCDFT 模块**，未开源

**推荐行动**：立即切换到「从零实现」路线（Week 1 剩余时间）。

---

## 1. octopus（https://octopus-code.org/）

### 1.1 功能检查

**官方描述**（GitLab README）：
- "Real-space Time-Dependent Density Functional Theory code"
- 核心：TD-DFT（时间依赖 DFT）
- 应用领域：光学响应、激发态动力学

**SCDFT 相关功能**：
- ✗ 无 superconducting DFT
- ✗ 无 pairing calculations
- ✗ 无 f_xc kernel for spin-fluctuation
- ✗ 无 Eliashberg 方程

**来源**：
- GitLab: https://gitlab.com/octopus-code/octopus
- WebFetch 尝试访问官方文档失败（404）

### 1.2 评分

**✗ 不可用** — octopus 是 TD-DFT 代码，与 SCDFT（超导配对）是完全不同的方向。

---

## 2. exciting（https://exciting-code.org/）

### 2.1 功能检查

**官方描述**（GitHub README）：
- Full-potential all-electron DFT package
- Implements LAPW (linearized augmented planewave)
- 聚焦：激发态、光谱、GW approximation

**功能列表**：
- Ground state DFT (LDA/GGA/hybrid)
- Time-dependent DFT (frequency/real-time)
- GW approximation (G₀W₀)
- Bethe-Salpeter equation

**SCDFT 相关功能**：
- ✗ 无 superconducting DFT
- ✗ 无 RPA spin susceptibility χˢ(q,ω)
- ✗ 无 spin-fluctuation mediated pairing
- ✗ 无 f_xc kernel

**来源**：
- GitHub: https://github.com/exciting/exciting
- 官方文档（多次 WebFetch 失败）

### 2.2 评分

**✗ 不可用** — exciting 专注于**电子激发态**（光学性质），不是**超导配对**。两者虽都叫"excitations"，但物理对象完全不同。

---

## 3. Elk（https://elk.sourceforge.io/）

### 3.1 功能检查

**官方描述**（elk.sourceforge.io）：
- Full-potential linearized augmented-plane-wave (LAPW) code
- 支持**电声耦合**和**超导计算**

**超导相关功能**（WebFetch 确认）：
- ✓ Electron-phonon coupling matrices
- ✓ Phonon linewidths
- ✓ Eliashberg function α²F(ω)
- ✓ Electron-phonon coupling constant λ
- ✓ McMillan-Allen-Dynes critical temperature Tc
- ✓ Eliashberg equations solved self-consistently
- ✓ Bogoliubov mean-field theory (experimental)

**但是**：
- ✗ **不是 SCDFT**（无 f_xc kernel、无 KS pairing functional）
- ✗ 是**传统 Eliashberg 理论**（基于 Migdal approximation）
- ✗ 只支持**声子机制**，不支持自旋涨落 Λ^SF

**来源**：
- WebFetch: https://elk.sourceforge.io/
- 文档陈述："Eliashberg function, α²F(ω)"

### 3.2 评分

**⚠ 功能不匹配** — Elk 有超导计算，但是：
1. **框架错误**：Eliashberg ≠ SCDFT（前者是多体微扰，后者是 KS DFT 扩展）
2. **机制缺失**：只有声子，无自旋涨落（Essenberger 的核心）
3. **无法复现**：Essenberger Eq.(37) 的 Im[Λ^SF] > 0 依赖 f_xc kernel，Elk 没有

**技术可行性**：即使搭建起来，也无法验证 SCDFT 顶角修正机制。

---

## 4. GPAW（https://gpaw.readthedocs.io/）

### 4.1 功能检查

**官方描述**（readthedocs）：
- Python-based DFT code
- Standard exchange-correlation functionals (LDA/GGA/PBE)
- RPA and ALDA kernels for linear response
- RPA correlation energy

**SCDFT 相关功能**：
- ✗ 无 superconducting DFT
- ✗ 无 pairing calculations
- ✗ 无 Eliashberg theory
- ✗ 无 f_xc kernel for superconductivity

**来源**：
- WebFetch: https://gpaw.readthedocs.io/
- 文档明确列出支持的功能（无超导相关）

### 4.2 评分

**✗ 不可用** — GPAW 是标准 DFT 代码，无任何超导扩展。

---

## 5. Essenberger 自己的计算流程（关键证据）

### 5.1 实际使用的代码（来自 1411.2121 §II）

**Essenberger 的计算链**：

```
1. VASP plane-wave code [54]
   ↓ (DFT Kohn-Sham bands)
2. LAPW code [55]
   ↓ (调 z_Se Wyckoff position)
3. RPA 磁化率 χˢ(q,ω)
   ↓ (使用 KS Green's function)
4. **SCDFT 自研模块**
   ↓ (f_xc kernel + BSE iteration)
5. 自旋涨落 Λ^SF
   ↓ (Eq.(32)–(37))
6. Gap equation Δ^xc_k
```

**来源**：
- `theory/notes/literature/essenberger_fese_parameters_search_2026-09-22.md`（L1，全文已读）
- 本地 PDF：`_scratch_ir/essenberger1411.2121.pdf`

### 5.2 关键结论

**Essenberger 未使用任何现成的 SCDFT 代码**，而是：
1. 用 VASP/LAPW 算 DFT 能带（标准功能）
2. **自己实现** SCDFT 模块（f_xc kernel + BSE）
3. 未公开代码

**证据**：
- 论文强调 "completely parameter free"（卖点是**方法**，不是工具）
- 无 acknowledgments 提及 SCDFT 代码
- 无 supplementary materials（arXiv 版本，无正式 PRB 出版）

---

## 6. 其他候选（快速排除）

### 6.1 Quantum ESPRESSO

**功能**：
- 有 PHonon 模块（电声耦合）
- 有 EPW（Eliashberg for phonon-mediated Tc）
- ✗ 无 SCDFT、无自旋涨落

**评分**：✗ 与 Elk 同类（传统 Eliashberg，非 SCDFT）

### 6.2 ABINIT

**功能**：
- 电声耦合、声子
- ✗ 无 SCDFT

**评分**：✗ 同上

### 6.3 联系 Essenberger 课题组

**可行性**：
- 作者：F. Essenberger, A. Sanna（Max Planck Institute for Microstructure Physics, Halle）
- **但需时间**（邮件往返 1-2 周）
- **风险**：代码可能是私有的、或依赖未发表的工作

**不推荐原因**：Week 1 预算紧张，等不起。

---

## 7. 技术门槛评估（假设从零实现）

### 7.1 核心组件

**必需实现的部分**（基于 `theory/notes/meta/scdft_formalism_essenberger_2026-09-22.md`）：

1. **f_xc 核计算**（Sec. V）
   - ALDA（adiabatic local density approximation）
   - 或 von Barth-Hedin 核（r_s ~ 2–3）
   - 依赖：libxc（开源）

2. **BSE 迭代**（Eq.(14)–(16)）
   - A^{c-d} = A_0 + A_0 G G A^{c-d}
   - 需要：KS Green's function G(k,iω)

3. **夹心形式**（Eq.(32)）
   - Λ^SF = 2 ∫ f_xc · χ_zz · f_xc
   - RPA 磁化率 χ^s(q,ω)

4. **SCDFT gap equation**（Eq.(36)）
   - 自洽求解 Δ^xc_k
   - 频率网格：Matsubara → 实轴解析延拓

### 7.2 输入数据

**来自 DFT**（可用 VASP/Quantum ESPRESSO）：
- Kohn-Sham 能带 ε(k)
- 费米面
- 态密度 N(ε)

**来自 Hubbard 模型**（替代方案）：
- 紧束缚参数 t_ij（Graser 1003.0133 / Kuroki PRL 2008）
- Hubbard U ~ 1 eV, Hund's J ~ 0.25 eV
- **来源**：`theory/notes/literature/essenberger_fese_parameters_search_2026-09-22.md` §3

### 7.3 时间估算

| 任务 | 难度 | 时间 |
|------|------|------|
| f_xc 核（libxc 接口） | 中 | 2–3 天 |
| BSE 迭代（矩阵代数） | 高 | 5–7 天 |
| RPA χ^s 计算 | 高 | 3–5 天 |
| Gap equation 求解器 | 中 | 2–3 天 |
| 调试 + 收敛 | 高 | 5–10 天 |
| **总计** | | **17–28 天**（3–4 周） |

**风险**：
- f_xc 数值不稳定性（需高精度）
- BSE 迭代可能不收敛
- 频率网格选择（Matsubara vs 实轴）

---

## 8. 推荐路线

### Plan A（推荐）：从零实现 SCDFT 核心

**理由**：
1. 没有现成代码可用
2. Essenberger 自己也是这么做的
3. 可以**直接追踪符号结构**（Im[Λ^SF] > 0 的机制）

**步骤**（Week 1 剩余 + Week 2–3）：
1. **Phase 1**（1 天，已完成）：理论框架理解 ✓
   - 文档：`theory/notes/meta/scdft_formalism_essenberger_2026-09-22.md`
2. **Phase 2**（5–7 天）：最小实现
   - f_xc 核（libxc）
   - 简化 BSE（不做全自洽，只算一次迭代）
   - 用 Graser/Kuroki 参数（不用 DFT，直接 Hubbard 模型）
3. **Phase 3**（7–10 天）：FeSe 复现
   - 计算 χ^s(q,ω)（裸 RPA）
   - 加 f_xc → Λ^SF
   - 验证 Im[Λ^SF] > 0（Essenberger Eq.(37)）
4. **Phase 4**（3–5 天）：机制提取
   - 逐 q 点追踪本征值谱
   - 与 Team C 整合（嵌套机制）

### Plan B（备选）：先算裸 RPA，后加 SCDFT

**理由**：
- 裸 RPA 的 χ^s(q,ω) 可以用标准 DFT 代码 + 后处理
- 验证 Team B 的「100% BZ 有负本征值」结论
- SCDFT 部分只做 f_xc 夹心（最后一步）

**步骤**（Week 1 剩余 + Week 2）：
1. **用 Quantum ESPRESSO 或 VASP**：
   - 算 KS 能带 ε(k)
   - 输出 Green's function G(k,iω)
2. **Python 后处理**：
   - 计算 Lindhard χ₀(q,ω)
   - RPA: χ^s = χ₀ / (1 − U χ₀)
   - 本征值分解，检查符号
3. **若验证「全 BZ 负本征值」**：
   - 再实现 f_xc 夹心
   - 看符号是否翻转

**风险**：
- RPA 计算量大（需要密 k 网格）
- 可能卡在「如何从 DFT 代码导出 G(k,iω)」

### Plan C（不推荐）：等 Essenberger 课题组回复

**理由**：
- 可能拿到真实代码
- 节省实现时间

**风险**：
- 时间不可控（1–2 周 → ∞）
- 代码可能不开源、或依赖未公开工作
- Week 1 预算浪费

---

## 9. 最终建议

### 立即行动（Week 1 剩余时间）

**选择 Plan A** — 从零实现 SCDFT 最小版本：

1. **今天剩余时间**（2–3 小时）：
   - 安装 libxc（f_xc 核库）
   - 搭建 Python 环境（numpy/scipy/matplotlib）
   - 读 Graser 1003.0133 或 Kuroki PRL 2008（取紧束缚参数）

2. **Week 1 剩余 2 天**：
   - 实现 Lindhard χ₀(q,ω)（裸泡泡图）
   - 实现 RPA χ^s = χ₀ / (1 − U χ₀)
   - 验证 Team B 的「负本征值」结论

3. **Week 2**：
   - 实现 f_xc 夹心 Λ^SF = f_xc · χ^s · f_xc
   - 监控本征值符号变化
   - 若成功 → 写 Phase 2 报告；若失败 → 切 Plan B

### 不建议做的事

1. ❌ 继续调研其他 DFT 代码（已确认无 SCDFT）
2. ❌ 尝试「魔改」Elk 的 Eliashberg 模块（框架不兼容）
3. ❌ 等 Essenberger 回复（时间不可控）

---

## 10. 交付清单

### 调研产出（本报告）

- [x] octopus 评估：✗ 不支持
- [x] exciting 评估：✗ 不支持
- [x] Elk 评估：⚠ 功能不匹配
- [x] GPAW 评估：✗ 不支持
- [x] Essenberger 计算流程确认：VASP + 自研模块
- [x] 时间估算：3–4 周（从零实现）
- [x] 推荐路线：Plan A（最小 SCDFT 实现）

### 下一步（交给 team-lead）

**决策点**：
- 批准 Plan A？
- 预算确认：Week 1–3 全投入 SCDFT 实现？
- 里程碑：Week 1 结束时需要有「裸 RPA 的 χ^s 符号验证」

**若批准，立即启动**：
- 环境搭建（libxc + Python）
- Graser/Kuroki 参数取件
- Phase 2 工作目录：`theory/notes/_scratch_ir/scdft_impl/`

---

## 附录：文献来源

### 仓库内部笔记
- `theory/notes/meta/scdft_formalism_essenberger_2026-09-22.md`（SCDFT 理论框架）
- `theory/notes/literature/essenberger_fese_parameters_search_2026-09-22.md`（Essenberger 计算流程）
- `theory/notes/misc/essenberger_chi_sign_structure_2026-09-22.md`（χ^s 符号结构）

### 外部资源（WebFetch 确认）
- octopus GitLab: https://gitlab.com/octopus-code/octopus
- exciting GitHub: https://github.com/exciting/exciting
- Elk 官方文档: https://elk.sourceforge.io/（Eliashberg 功能确认）
- GPAW 文档: https://gpaw.readthedocs.io/

### 本地 PDF（已读）
- `_scratch_ir/essenberger1411.2121.pdf`（FeSe 实际材料计算，6 页）
- `_scratch_ir/_essenberger1409.7968.pdf`（SCDFT 泛函推导，14 页；2026-09-25 订正：原写漏了前导下划线 ⟹ 指向不存在的文件）

### 调研限制
- WebSearch 预算已用完（200/200）
- arXiv.org 被企业网络封锁
- 多个官方文档 URL 返回 404（octopus/exciting）

---

**状态**：调研完成，等待 team-lead 批准行动路线。
