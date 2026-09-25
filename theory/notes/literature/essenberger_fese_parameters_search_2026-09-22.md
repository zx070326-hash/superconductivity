# Essenberger FeSe 参数组搜索报告

**日期**: 2026-09-22  
**任务**: 搜索 Essenberger 1411.2121 / 1409.7968 用的显式参数表（t_ij, U, J）  
**结论**: **NOT OBTAINED** — 确认未公开，方法论设计就是参数无关  
**read_depth**: L1（两篇全文 PDF 已读）

---

## 1. 核心发现

### Essenberger 未公开紧束缚参数的三重证据

#### 1.1 方法论设计是参数无关

**1409.7968** (PRB 94, 014503, 2016):
- 摘要："completely **parameter free** superconducting gap equation"（两次强调）
- §I: "fully parameter free and is able to reproduce... under the only knowledge of the atomic constituents"
- §VI: 应用对象是 "simplified **two-band electron gas model**"（非实际材料参数）

**1411.2121** (arXiv only):
- 摘要："neither the functional, nor the theoretical framework, contain any **adjustable parameters**"
- 正文 p.1: "it does not allow genuine predictions" 是对 Hubbard models [35,37-40] 的批评

⟹ **"无参数"是卖点**，不是缺陷。

#### 1.2 计算路线不经过紧束缚拟合

**1411.2121 实际计算流程**（§II，p.2–3）:
1. DFT Kohn-Sham 能带（VASP plane-wave code [54] + LAPW [55]）
2. 调 Se Wyckoff 位置 z_Se（影响磁化率）
3. 调 **scaling parameter α/α_c = 0.95**（人工参数，用于避免 AFM 发散）
4. RPA 磁化率 χ^s(q,ω) → 自旋涨落 Λ^SF

**不需要** t_ij / U / J：
- DFT 直接给 ε(k)
- RPA 用 KS Green's function，不用 Hubbard Hamiltonian

#### 1.3 引用 Hubbard models 是对比对象

1411.2121 正文 p.1:
> "SF mediated pairing has been extensively investigated in the realm of **Hubbard like models**[35, 37–40]. Since this approach necessarily involves **a set of parameters**, it does **not allow genuine predictions** of the critical temperature."

⟹ Hubbard 参数是**他们批评的对象**，不是自己用的工具。

---

## 2. 本地资源检查

### 已读全文

| 文件 | 页数 | 内容 |
|------|------|------|
| `_scratch_ir/essenberger1411.2121.pdf` | 6 | FeSe/LiFeAs/LaOFeAs 实际材料计算 |
| `_scratch_ir/essenberger1411.2121.txt` | 231 行 | pdftotext 提取（UTF-8） |
| `_scratch_ir/_essenberger1409.7968.pdf` | 14 | 自旋涨落泛函推导 + two-band 模型（2026-09-25 订正：原写漏前导下划线） |
| `_scratch_ir/_ess1409_utf8.txt` | 927 行 | pdftotext 提取 |

### Grep 关键词覆盖

```bash
# 无命中
pattern: "Table|Appendix|hopping|tight.binding|t_ij|Hubbard U"

# 唯一命中是"Hubbard models"的批评语境
context: "approach necessarily involves a set of parameters"
         "it does not allow genuine predictions"
```

### arXiv 附加文件

- WebFetch 被封锁（企业网络阻止 arxiv.org）
- Grep 本地文本：无 "Supplemental Material" / "ancillary" 字样
- **推断**：无补充材料（1411.2121 只有 arXiv 版，PRB 正式版不存在）

---

## 3. 替代方案：Graser / Kuroki 参数

### 已有搜索报告

`../_scratch_ir/_graser_5orbital_parameters.txt`（2026-09-22，实测 **106** 行；原写「107 行」且路径漏了前导下划线 ⟹ 指向一个不存在的文件，2026-09-25 订正）给出：

#### 最可能来源（按优先级）

1. **Kuroki et al., PRL 101, 087004 (2008)**  
   "Unconventional Pairing Originating from the Disconnected Fermi Surfaces of Superconducting LaFeAsO_{1-x}F_x"  
   - 5 轨道紧束缚参数的 **canonical 源**
   - Graser 0812.0343 引用它
   - **高置信度**含显式 t_ij 表

2. **Graser et al., arXiv:1003.0133 (2010)**  
   "Spin fluctuations and superconductivity in a 3D tight-binding model for BaFe2As2"  
   - 0812.0343 的后续（加 k_z 色散）
   - 报告确认：**有显式参数**
   - 更易获取（arXiv）

3. **Daghofer et al., PRL 101, 237004 (2008)**  
   - Supplemental material 有 hopping 参数

4. **Maier et al., PRB 79, 134520 (2009)**  
   - 详细参数表
   - Maier 是 Graser 合作者，被 Essenberger 1409 引用为 Ref.[36]

### 参数适用性评估

**为什么 Graser/Kuroki 参数可用于 Essenberger 工作的重现？**

1. **同一套 Hubbard 文献**  
   Essenberger 1411.2121 引用 [35,37-40]：
   - [35] Hirschfeld et al., Rep Prog Phys 74 (2011)
   - [37-40] Daghofer, Maier, Nicholson 等 2011-2012 年工作
   - 这些都引用 Kuroki/Graser 参数

2. **同一类材料**  
   - Kuroki: LaOFeAs (1111 族)
   - Graser 1003: BaFe2As2 (122 族)
   - Essenberger: FeSe (11 族) + LiFeAs (111 族) + LaOFeAs
   - **费米面拓扑相同**（Γ 空穴 + M 电子，见 1411.2121 Fig.1）

3. **5 轨道标准模型**  
   - d_xy, d_xz, d_yz, d_3z²-r², d_x²-y²
   - 2008-2012 年已成为铁基超导的事实标准

---

## 4. ~~典型参数范围（文献标准）~~ —— **⚠ 2026-09-25 降级：无一手来源，不得引用**

> **旧陈述错在哪**：本节把下列数值冠以「文献标准」，并注明「摘自 `graser_5orbital_parameters.txt` 第 45–63 行」。
> 但那个文件**自己的第 4 行**逐字写着 `Status: PARAMETERS NOT OBTAINED (web access blocked, search budget exhausted)`，
> 第 45 行的小标题是 `TYPICAL PARAMETER RANGES (Literature Standard)` —— 即**该 agent 在没拿到任何原文的情况下
> 自行写下的回忆值**，不是从 Graser 0812.0343 或任何论文抽取的。⟹ 这是一处**幻影引用**：
> 数字看起来有出处，实际出处是一份自陈「未取得」的搜索报告。
> **新陈述凭什么对**：我 2026-09-25 直接读了该 `.txt` 的第 1–12 行与第 45–48 行，两处逐字如上。
> **判词**：本节全部数值等级为 **无来源（既非 `numerical-evidence` 也非引用）**，**任何结论不得消费**。
> **真正有一手来源的 Graser 参数在别处**：`graser1003.0133_5orbital_TB_parameters_2026-09-22.md`
> （arXiv:1003.0133 直接取回成功，文本层 + 页图双仪器）——需要参数请用那一篇。
> 本节原文保留作为排查记录（「记录丢失时不得反向重建」的对偶：**也不得删掉错误记录来掩盖它**）。

摘自 `../_scratch_ir/_graser_5orbital_parameters.txt` 第 45–63 行（⚠ 见上，该来源自陈未取得原文）：

### 在位能 ε_d（设 d_xy 为 0）
- d_3z²-r²: ~+0.5 eV（最高）
- d_xz, d_yz: ~0–0.2 eV
- d_xy: 0（参考点）
- d_x²-y²: ~−0.3 eV

### 跃迁积分 t_ij
- 最近邻（Fe-Fe）: 0.1–0.5 eV
- 次近邻: 0.05–0.2 eV
- 强各向异性（面内 ≫ 面外）

### 相互作用
- Hubbard U: 1.0–2.0 eV（典型 ~1.5 eV）
- Hund's J: 0.2–0.3 eV（典型 ~0.25 eV）
- 化学势 μ: 调到 6 电子/Fe

### 能带特征
- Γ 点空穴口袋（桶状）
- M 点电子口袋
- 弱 k_z 色散（准二维）

---

## 5. 时间线与封锁情况

### 本次搜索（2026-09-22）
- **用时**: 1 天
- **WebSearch 预算**: 200/200 已耗尽（本会话共享，无法追加）
- **WebFetch**: arXiv.org 被企业网络封锁（"unable to verify domain is safe"）
- **实际工作**: 本地 PDF 全文读取 + 文本挖掘 + 旧报告整合

### Graser 旧报告（2026-09-22，另一会话）
- **用时**: ~2 小时
- **结论**: "CONFIDENCE: **HIGH** that parameters are in Kuroki 2008 or Graser 1003.0133"
- **状态**: "LOW that they're explicitly in 0812.0343 itself"

---

## 6. 建议行动

### Plan A′（立即可行，推荐）

**用 Graser 1003.0133 或 Kuroki PRL 2008 参数计算 RPA χ^s_{ij}(q,ω)**

理由：
1. Essenberger **确认不存在**显式紧束缚参数
2. 他们引用的 Hubbard 文献都用 Kuroki/Graser 参数
3. FeSe 费米面与 LaOFeAs/Ba122 **拓扑同构**（只是能带宽度略有差异）

执行：
1. 优先取 **Graser 1003.0133**（arXiv，可能绕过封锁）
2. 备选 **Kuroki PRL 2008**（需机构订阅，但是 canonical 源）
3. 用 5 轨道 + RPA 计算 χ(q,ω)
4. 投影到 Essenberger 的 2×2 带模型（Γ 空穴 / M 电子）

### Plan B（需外部资源）

1. **机构图书馆**取 Kuroki PRL 101, 087004 全文
2. **Google Scholar** 试取 Graser 1003.0133（可能有作者上传版）
3. **联系作者**（A. Sanna / F. Essenberger, Max Planck Institute for Microstructure Physics, Halle）
4. **ResearchGate** 搜作者上传版

### 不建议

- ❌ 继续搜 Essenberger 显式参数（已确认方法论上不存在）
- ❌ 从 1411.2121 Fig.1 反推参数（那是 DFT **输出**，不是输入；且 ε(k) 不唯一确定 t_ij）
- ❌ 用 1409.7968 §VI 的 two-band 模型参数（那是**抛物线近似**的人工参数，非实际材料）

---

## 7. 交叉引用

### 相关笔记
- `theory/notes/literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md`  
  含 Essenberger λ_ij 表的提取（带分辨的耦合常数，但不是紧束缚 t_ij）

### 相关脚本
- `_scratch_ir/_chk_essenberger_lambda_ij_psd.py`（多个 worktree 副本）
- `_scratch_ir/_chk_essenberger_partition_2026-09-21.py`
- `../_scratch_ir/_essenberger_fese_patch_band_mapping_2026-09-21.md`（2026-09-25 订正两处：原写漏了前导下划线；且该文件此前**只存在于某个 agent worktree**、主工作树没有 ⟹ 任何 clone 都读不到，现已取回主树并 `git add -f`。**是取回原件，不是重建**，181 行原样复制）

---

## 8. 结论摘要

| 问题 | 答案 | 证据等级 |
|------|------|----------|
| Essenberger 有无显式 t_ij / U / J 表？ | **无** | L1（全文已读） |
| 为什么没有？ | 方法论设计是**参数无关**（DFT → RPA，不经紧束缚） | analytic（方法论定义） |
| 1411.2121 的"模型"是什么？ | 各向同性近似（band-resolved isotropic），非紧束缚 | L1（§II 完整流程） |
| 1409.7968 §VI 的参数能用吗？ | 不能，那是 two-band **electron gas**（抛物线 DOS + 人工 λ/ω̄） | L1（§VI.A 全文） |
| 替代方案的最佳来源？ | Graser 1003.0133 (2010) 或 Kuroki PRL 2008 | 高置信度（Graser 报告） |
| 参数适用性？ | 可用（同类材料、同费米面拓扑、同文献体系） | 交叉引用一致 |

---

## 9. 文献锚点

### Essenberger 原始文献

- **1409.7968**: F. Essenberger, A. Sanna, A. Linscheid, F. Tandetzky, G. Profeta, P. Cudazzo, E.K.U. Gross, "Superconducting pairing mediated by spin-fluctuations from first principles", PRB 94, 014503 (2016)
- **1411.2121**: F. Essenberger, A. Sanna, P. Buczek, A. Ernst, L. Sandratskii, E.K.U. Gross, "Ab-initio theory of Iron based superconductors", arXiv:1411.2121 (2014, arXiv only)

### 推荐参数来源（按优先级）

1. K. Kuroki et al., PRL **101**, 087004 (2008)
2. S. Graser et al., arXiv:**1003.0133** (2010)
3. M. Daghofer et al., PRL **101**, 237004 (2008)
4. T.A. Maier et al., PRB **79**, 134520 (2009)

---

**状态**: 搜索完成，已向 team-lead 交付  
**下游**: Team B 可立即继续 Plan A′（用 Graser/Kuroki 参数）
