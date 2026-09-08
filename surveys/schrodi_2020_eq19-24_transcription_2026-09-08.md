# Schrödi, Oppeneer, Aperis PRB 102, 024503 (2020) Eq. (19)–(24) 逐字转录（researcher）

日期：2026-09-08。PDF `library/pdf/D1/2020_Fabian_Schr_di_1911.12872.pdf`，页码 = PDF 页码。

## 方程（p.4–5）
- (19) Σ̂_m^(iso) = T² P Σ_{m',m''} λ_{m−m'} λ_{m'−m''} ρ̂₃ ĝ_{m'} ρ̂₃ ĝ_{m''} ρ̂₃ ĝ_{m''−m'+m} ρ̂₃
- (20) P = (1/N₀³) Σ_{k,k',k''} δ(ξ_k)δ(ξ_k')δ(ξ_k'')δ(ξ_{k''−k'+k}) ≥ 0（"momentum conservation of processes at the Fermi level"）
- (21) ĝ_{k,m} = ∫dε Ĝ_{k,m}(ε)；(22) ĝ_m = −π (iω_m ρ̂₀ + Δ_m ρ̂₁)/√(ω_m² + Δ_m²)
- (23) **Z 方程**：Z_m = 1 + (πT/ω_m)Σ λ_{m−m'} γ^ω_{m'} + P(π³T²/ω_m) Σ_{m'm''} λ_{m−m'}λ_{m'−m''} ( −γ^ω' γ^Δ'' γ^Δ''' + γ^Δ' γ^ω'' γ^Δ''' − γ^Δ' γ^Δ'' γ^ω''' − γ^ω' γ^ω'' γ^ω''' )
- (24) **Δ 方程**：Δ_m = (πT/Z_m)Σ λ_{m−m'} γ^Δ_{m'} + P(π³T²/Z_m) Σ λλ ( −γ^Δ' γ^Δ'' γ^Δ''' − γ^ω' γ^ω'' γ^Δ''' **+ γ^ω' γ^Δ'' γ^ω''' − γ^Δ' γ^ω'' γ^ω'''** )　←第 3、4 项按 PDF 第 5 页 6× 放大核对更正（初版把 Δ/ω 位置对调；reviewer 与 p2-schrodi implementer 两次独立核对，2026-09-08）
- γ^ω_m = ω_m/√(ω_m²+Δ_m²)，γ^Δ_m = Δ_m/√(ω_m²+Δ_m²)；' = m'，'' = m''，''' = m''−m'+m
- λ_{m−m'} = λ Ω²/(Ω² + ω²_{m−m'})（拼接自 Eq. 5、28、29，推断）

两式各"3 负 1 正"，但不是 ω↔Δ 互换。

## 线性化（主会话）
Δ→0 时 γ^ω → sgn(ω)，γ^Δ → Δ/|ω|：
- Z 方程只剩第 4 项：−P(π³T²/ω_m)Σλλ s's''s'''（s = sgn ω）——m', m'' 靠近 m 时三个符号相同，乘积 = sgn(ω_m) → 顶角**降低** Z（与 Schrödi Table I 小 α 时 λ_Z^(2) < 0 一致）
- Δ 方程线性项（按更正后的 Eq. 24）：−s's''Δ'''/|ω'''| + s's'''Δ''/|ω''| − s''s'''Δ'/|ω'|——低频同符号区 ≈ −Δ/|ω|：顶角**削弱配对核**（reviewer：⟨s's''s'''⟩ 在 m=0 只有 0.44、m=1 为 0.73，m≥5 变号，"≈−Δ/|ω|" 只在最低两个 Matsubara 频率成立，是启发不是受控展开）
- 净效应：Schrödi p.12 "the FSR isotropic approximation ... seems to **systematically produce a negative correction**"——Δ 通道的削弱压过 Z 通道的减小；**`p2_schrodi_fsr.py` 数值确认**：Z-only δ/P > 0、Δ-only δ/P < 0 且约 2× Z、两通道净负，符号与 EPW FSR 一致

## 数值（一带模型）
- Table I（2D，λ=1.5，p.10）：λ_Z^(2) = −0.55 / −0.05 / +0.41（α = 0.05 / 0.1 / 1）——顶角对 Z 的贡献随 α 变号
- Fig. 6(c)（完整各向异性解）：Tc 75 K（α=0.05，压低）/ 120 K（0.1，增强）/ 135 K（1，增强），绝热 ≈100 K
- p.9："the FSR isotropic approximation severely underestimates the effect of the first vertex correction and is therefore not valid"
- p.12："the more isotropic the resulting vertex function is, the more positive its effect on the pairing becomes. Interestingly, these results are opposite to those of Ref. [28]."

## 与 Mishra 2507 的可比性（researcher 告诫）
Mishra 的 FSR 只是常数 DOS，保留完整各向异性 g；Schrödi 的 FSR-isotropic 还把耦合各向同性平均掉。Mishra 的 FSR 顶角压低（246→198，19.5%）比 FBW（238→226，5%）**更大**，与"FSR 低估顶角"的说法相反——两个 FSR 不是一回事。
