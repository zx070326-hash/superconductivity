# FeSe 费米面嵌套机制分析 - Phase 1 准备工作总结

**日期**：2026-09-22  
**任务**：Team C - 追踪 FeSe 反铁磁矢量 Q 的作用（费米面嵌套机制）  
**状态**：Phase 1 完成，等待 Team B 数值结果

---

## 一、背景理解

### 1.1 FeSe 费米面结构（from Essenberger 1411.2121 Fig. 1）

- **Band +（空穴）**：Γ 点 (0,0) 中心的圆柱形费米面
- **Band −（电子）**：M 点 (π,0) 中心的圆柱形费米面
- **2×2 带结构**：整篇论文使用此简化（L.254 "two-band isotropic approximation used so far"）

### 1.2 物理图像（from 研究计划 L.164-170）

```
铁基超导自旋涨落峰化在 Q_AFM ≈ (π,0)
不同轨道费米面被 Q_AFM "连接"（嵌套）
⟹ χ^s_{αβ}(Q_AFM) 都看到同一个峰
⟹ 可能都是正的（相干增强，无干涉抵消）
```

### 1.3 Essenberger 已有的证据（Fig. 4c）

- **已画出嵌套函数沿 Γ-M 线**（normalized nesting function）
- 文字："nesting for the SF pairing (see Fig. c) and leads to an increased interband pairing"
- 这正是我们要计算的对象 ξ_{αβ}(q)

---

## 二、理论框架

### 2.1 嵌套函数定义

```
ξ_{αβ}(q) = ∫ dk δ(ε_α(k)) δ(ε_β(k+q))
```

**物理意义**：
- 衡量波矢 q 将轨道 α 费米面上的态散射到轨道 β 费米面上的相空间体积
- 高嵌套度 → 强散射 → χ^s_{αβ}(q) 大
- 低嵌套度 → 弱散射 → χ^s_{αβ}(q) 小（可能为负？）

### 2.2 核心假说（待检验）

**正面假说**：χ^s_{αβ}(q) > 0（逐元）⟺ 所有轨道对在 Q_AFM ≈ (π,0) 附近都有高嵌套度

**反面假说**（如果 Team B 发现负非对角）：χ^s_{αβ}(q) < 0 的 q 点 ⟺ ξ_{αβ}(q) ≈ 0（嵌套度接近零）

### 2.3 与 RPA 理论的联系

从 Lindhard 函数：
```
χ₀_{αβ}(q,ω) = ∑_k [f(ε_α(k)) - f(ε_β(k+q))] / [ω - ε_β(k+q) + ε_α(k) + iη]
```

静态极限 ω=0、零温下：
```
Re χ₀_{αβ}(q,0) ~ ξ_{αβ}(q)  （嵌套度主导）
```

RPA 增强后：
```
χ^RPA_{αβ}(q) = [χ₀(q) @ inv(1 - U @ χ₀(q))]_{αβ}
```

如果 χ₀ 对角主导 + RPA 保持对角主导 → 所有 χ^RPA_{αβ} 同号

---

## 三、技术方案

### 3.1 已准备的工具

脚本：`theory/notes/_scratch_ir/_fese_nesting_framework.py`

**核心函数**：
1. `SimpleFeSeFermiSurface`：简化 2 带紧束缚模型（占位参数）
2. `compute_nesting_function(q, α, β)`：计算单点嵌套度
3. `scan_brillouin_zone()`：生成完整嵌套度热图
4. `check_special_points()`：检验 Q_AFM、Γ、M、X 点

**输出**：
- 3 个嵌套度热图：ξ_{++}(q), ξ_{−−}(q), ξ_{+−}(q)
- 与 Team B 的 χ^s 符号图逐点比较

### 3.2 需要 Team B 提供的数据

1. **χ^s_{ij}(q,ω) 符号结构**：
   - 在 (q_x, q_y) 平面上哪些区域 Im χ^s_{12}(q) > 0 / < 0
   - 关键问题：**非对角元是否存在负值区域**

2. **可选（如果可获得）**：
   - Essenberger 使用的紧束缚参数（跃迁积分 t_ij、化学势 μ）
   - 从 Kemper NJP 2010 或 Maier PRB 79 获取

### 3.3 分析流程（Phase 2）

```
收到 Team B 数值后：
  │
  ├─ 情形 A：Im χ^s_{ij} 全为正
  │    └─ 计算所有 q 点的 ξ_{αβ}(q)
  │    └─ 证明：Q_AFM 附近所有轨道对都高嵌套
  │    └─ 结论：费米面全局嵌套 → 逐元正
  │
  └─ 情形 B：存在负非对角
       └─ 识别 Im χ^s_{12}(q) < 0 的 q 区域
       └─ 计算这些 q 点的 ξ_{+−}(q)
       └─ 检验：负区域 ⟺ 低嵌套度
       └─ 结论：费米面几何限制 → 符号结构
```

---

## 四、与其他队的依赖关系

### 4.1 Team B（essenberger-chi-reconstruction）

**CRITICAL 依赖**：我的工作完全依赖 Team B 先给出符号结构

- Team B 任务：从 Essenberger 参数组计算 χ^s_{ij}(q,ω)
- Team B 产出：符号图 + 数值表
- 我的输入：Team B 的符号图
- 我的产出：物理机制解释（费米面嵌套 → 符号结构）

### 4.2 Team A（rpa-diag-dominance-proof）

**独立并行**：Team A 在证明一般 RPA 定理，我在解释 FeSe 特例

- 如果 Team A 成功：一般定理 + FeSe 验证 = 双重确认
- 如果 Team A 受阻：FeSe 经验仍有价值（材料特异性）
- 两队结果互不依赖

---

## 五、已读文献

1. **Essenberger PRB 94, 014503 (2016) = arXiv:1411.2121**
   - 读 pp.1-4（Fig. 1 费米面、Fig. 3 配对函数、Fig. 4 嵌套+能隙、Fig. 5 k 分辨能隙）
   - 关键：Fig. 4c 嵌套函数已存在

2. **Essenberger PRB 90, 214504 (2014) = arXiv:1409.7968**
   - 读 pp.1-3（SCDFT 形式论、粒子-空穴传播子）
   - 这是 1411.2121 的理论基础

3. **研究计划**：`theory/notes/meta/rpa_offdiag_positivity_research_plan_2026-09-22.md`
   - 路线 2II（L.162-191）明确了我的任务范围

4. **理论缺口**：`theory/notes/literature/rpa_chi_theory_gap_2026-09-21.md`
   - 确认 RPA 层次无 no-go 定理
   - 0/23000 数值证据（负非对角 + PSD 未见）

---

## 六、当前状态与下一步

### 6.1 Phase 1 完成标志

- ✓ 背景材料全部读取
- ✓ 物理图像理解清晰
- ✓ 技术框架已搭建
- ✓ 依赖关系明确
- ✓ 消息已发送给 Team B

### 6.2 等待事项

**等待 Team B 交付**（预计 1 周）：
1. χ^s_{ij}(q,ω) 符号结构普查
2. 符号图（热图或散点图）
3. 关键问题答案：**非对角元是否逐元为正？**

### 6.3 Phase 2 启动条件

收到 Team B 消息后立刻启动：
1. 运行 `_fese_nesting_framework.py` 生成嵌套度热图
2. 与 Team B 符号图定量关联（逐点或逐区域）
3. 撰写物理机制分析笔记
4. 交付时间：1-2 周

---

**Phase 1 准备工作完成。处于就绪状态，等待 Team B。**
