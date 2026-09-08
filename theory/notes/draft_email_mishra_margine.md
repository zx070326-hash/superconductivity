# 致 Mishra / Margine 组的数据请求（草稿，待用户决定是否发送）

收件人：Shashi B. Mishra (mshashi125@gmail.com), Elena R. Margine (rmargine@binghamton.edu)
主题：Request for isotope (deuteride) vertex-correction data and DOS from arXiv:2608.14949 / 2510.25902

---

Dear Dr. Mishra and Prof. Margine,

We have been working with the Figshare dataset accompanying your Annalen der Physik paper on H3S/D3S (10.6084/m9.figshare.31136635) and found it unusually valuable: because it contains both isotopes under the four treatments (FSR/FBW × adiabatic/vertex), it allows a complete decomposition of the hydrogen isotope coefficient α_H into anharmonic, vertex, DOS, and Coulomb contributions. From your data we obtain, at 200 GPa (anharmonic, μ* = 0.16): α_H = 0.337 (FSR-adia) → 0.252 (FSR-na) → 0.313 (FBW-na), i.e. the vertex correction lowers α_H by ≈0.085 and the energy-dependent DOS cancels about two thirds of that. Combined with the μ* cutoff renormalization, the four-factor total (0.29–0.38 between 160 and 200 GPa) matches the experimental α_H ≈ 0.3 of Drozdov et al. We also find that the vertex effect on ln Tc scales as δ_D/δ_H ≈ 0.66–0.72, consistent with the leading-order M^{-1/2} expectation.

We would like to test whether these findings generalize beyond H3S, and would be very grateful if you could share, for the compounds of arXiv:2608.14949:

1. For YH6 (and, if available, YH9 or LaBeH8): the corresponding **deuteride** calculations (α²F(ω) and the gap-vs-T scans for FSR/FBW × adia/na), analogous to D3S in the Figshare set. Even a single pressure would suffice.
2. For the hydrides in Table I: the **FSR** (Fermi-surface-restricted) adiabatic and vertex-corrected Tc values, in addition to the FBW values published, since we find the FSR quantity is the one that isolates the vertex effect from its interplay with the DOS.
3. The numerical **N(ε)** within ±1 eV of ε_F (e.g. the EPW `.dos` files) for each compound, to compare with a simple DOS-based estimate of the static vertex correction, P_s/λ = (ω₀/2)∫dε N̄(ε)/(|ε|+ω₀)² − 1, which reproduces your H3S vertex suppression and its pressure dependence rather well and which we would like to test across your material set (YH9, with ε_F in a DOS valley, is the most interesting case).

We will of course cite the relevant papers and datasets and would be happy to share our analysis scripts and results.

With best regards,
[用户署名 / 单位]

---

## 备注（给用户）
- 第 3 条里的公式是我们 P2 v0.3 的结果，尚在 reviewer 核验中；发信前等核验结果，若被推翻则删去该句
- 数字来自 `theory/derivations/P1_isotope_master_formula.md` §4.1 与 `P2_vertex_sign_criterion.md` §1
- 若不想透露判据，可只提第 1、2 条
