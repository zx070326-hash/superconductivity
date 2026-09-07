# P1 三轮支撑：原文补读（LaD10 相/压力、Villa-Cortés 结果、Errea 2020 约定、Troyan μ*(6Ry)）

调研 agent：p1-data-round2（sonnet / researcher，第三轮）　日期：2026-09-08
方法：Read 本地 PDF，页码 + 英文原句。原始产出，未手改（两条回传合并）。

## CSV

```csv
task,material,isotope,pressure_GPa,quantity,value,convention,source_file,page,quote_or_note
1,LaH10,H,150,structure,"Fm-3m，XRD 确认",NA,P1/2019_A_P_Drozdov_1812.01561.pdf,p.4,"The structure of the superconducting phase with Tc~250K at 150GPa was determined as face-centered cubic (fcc) lattice Fm3m"
1,LaD10,D,NA,structure,"Fm-3m，仅为假设，未经 XRD 独立确认",NA,P1/2019_A_P_Drozdov_1812.01561.pdf,"p.4 + Ext.Data Fig.5图注","正文: ""We can naturally suppose that...the deuteride with maximum Tc (168K) have the same crystal structure, that is the Fm3m phase""；Ext.Data Fig.5图注: ""we did not determine the structure for the phase with the highest Tc ~160K"""
1,"LaD10.6(P4/nmm)",D,142,structure,"P4/nmm，XRD 确认",NA,P1/2019_A_P_Drozdov_1812.01561.pdf,Fig.3c,"the La+D2 mixture with a Tc~140K at 142GPa"；XRD 精修给 LaD10.6
1,LaH10/LaD10,H&D,NA,alpha_pair_fcc,"Tc_H=250K vs Tc_D=168K → α=0.57（压力不匹配，仅按'同结构'假设配对，非同压对比）",NA,P1/2019_A_P_Drozdov_1812.01561.pdf,p.4,"the hydride with the maximum Tc (250K) and the deuteride with maximum Tc(168K)...Under this assumption, we obtained the isotope coefficient α of about 0.57...for the fcc phase"
1,LaH10/LaD10,H&D,NA,alpha_pair_P4nmm,"Tc_H=215K vs Tc_D=140K → α=0.62",NA,P1/2019_A_P_Drozdov_1812.01561.pdf,p.4,"the hydride with Tc=215K and the deuteride with Tc=140K should have a common P4/nmm structure...α of about...0.62"
1,LaD10,D,"139 / 145 / 152（插图读值，低置信度）",Tc_exp,"约163-165 / 158-160 / 150 K（从 Ext.Data Fig.4b 小插图估读，非精确数字）",NA,P1/2019_A_P_Drozdov_1812.01561.pdf,Ext.Data Fig.4b,"Superconducting transitions at different pressures (pressure determined from the D2 vibron scale)"；插图轴范围 140-165K
1,LaD10,D,NA,Tc_qualitative_internal_mismatch,"正文: ~168K, ~140K；Ext.Data Fig.5图注: ~160K, 140K, ~70K",NA,P1/2019_A_P_Drozdov_1812.01561.pdf,"正文 vs Ext.Data Fig.5图注","""three transitions in lanthanum deuterides with Tc~160K,140K and ~70K"""
1,LaD10,D,NA,Tc_70K_point_pressure,NOT_FOUND,NA,P1/2019_A_P_Drozdov_1812.01561.pdf,全文已通读两遍,"氘化物~70K 转变全文没有给出对应压力"
2,H3S,H&D,162,alpha_el-ph,"0.42（N(ε)）/ 0.45（N(0)）",N(ε) 与 N(0) 两套均给出,P1/2022_S_Villa-Cort_s_2103.12055.pdf,p.8,"α^el-ph becomes positive, reaching a total value of 0.42(0.45) at 162 GPa and 0.39(0.41) at 250 GPa, when N(ε)(N(0)) is used"
2,H3S,H&D,250,alpha_el-ph,"0.39（N(ε)）/ 0.41（N(0)）",同上,P1/2022_S_Villa-Cort_s_2103.12055.pdf,p.8,同上引文
2,H3S,H&D,162,alpha_el-el,"-0.09（N(ε)）/ -0.07（N(0)）",同上,P1/2022_S_Villa-Cort_s_2103.12055.pdf,p.9,"While for 162 GPa its value is quite marginal, -0.09(-0.07) for Nε(N0)"
2,H3S,H&D,250,alpha_el-el,"-0.13（N(ε)）/ -0.12（N(0)）",同上,P1/2022_S_Villa-Cort_s_2103.12055.pdf,p.9,"it improves as pressure increases, reaching -0.13(-0.12) with Nε(N0)"
2,H3S,note,NA,alpha_Total_N0_vs_Neps,"N(0) 版本比 N(ε) 版本在 250GPa 处最多高估 24%",NA,P1/2022_S_Villa-Cort_s_2103.12055.pdf,p.8,"α^Total with N(0) overestimates the previous discussed α values by an amount as large as 24% at 250 GPa"
2,H3S,note,NA,mu_star_convention,"μ* 不是固定值，是在每个压力点各自独立拟合的自由参数",逐点拟合非固定,P1/2022_S_Villa-Cort_s_2103.12055.pdf,p.3,"Once ωc is set, the only adjustable parameter left is μ*, which is fitted at each pressure of interest."
2,H3S,note,NA,alpha_pointwise_full_series,"NOT_TABULATED——逐点数值只以 Fig.12(a)(b)(c) 曲线形式给出",NA,P1/2022_S_Villa-Cort_s_2103.12055.pdf,"Fig.12(a)(b)(c)","需要数字化图像"
3,LaH10/LaD10,note,163/159,mu_star_Mc_AD列,"μ*=0.1 固定（经验值）",μ*=0.1 fixed,D3/2020_Ion_Errea_1907.11916.pdf,p.8 正文Methods,"we also include the values obtained with McMillan equation and Allen-Dynes modified McMillan equation (μ*=0.1)"
3,LaH10/LaD10,note,NA,mu_star_MEani列,"无经验 μ*；用 k 平均静态（ω=0）RPA 屏蔽库仑相互作用，完全从头计算",ab initio静态RPA,D3/2020_Ion_Errea_1907.11916.pdf,p.7,"the screened Coulomb interaction... is approximated by k-averaged static Coulomb interaction within the random phase approximation"
3,LaH10/LaD10,note,NA,mu_star_SCDFT列,"无经验 μ*；用能量依赖的 μ(ε,ε') 屏蔽库仑核，完全从头计算",ab initio能量依赖,D3/2020_Ion_Errea_1907.11916.pdf,p.7,"Eq.(3) allows us to calculate Tc completely ab initio, without introducing an empirical μ* parameter"
3,LaH10/LaD10,note,~160,alpha_maintext,"理论 α=0.43（160GPa 附近）；正文称与'experimentally reported α of 0.46'一致",NA,D3/2020_Ion_Errea_1907.11916.pdf,p.4,"the isotope coefficient α=−[lnTc(LaD10)−lnTc(LaH10)]/ln2 is close to 0.5 (0.43 around 160 GPa), the expected value in BCS theory, and it is in agreement with the experimentally reported α of 0.46"
3,LaH10,note,NA,LaD10_agreement_statement,"作者称 LaD10 Fm-3m 相的计算与'实验点'一致",NA,D3/2020_Ion_Errea_1907.11916.pdf,p.4,"Our calculations for LaD10 in the Fm-3m phase are also in agreement with the experimental point reported"
4,YH6,note,NA,mu_star_6Ry_definition,"库仑项写成 μ*θ(ωc−|ωm|)，即对截止内(|ωm|<ωc=6Ry)所有 Matsubara 频率取同一常数 μ*，截止外为零——盒式锐截止，不是对数 Tolmachev 重整化",ωc=6Ry锐截止,P1/2021_Ivan_A_Troyan_1908.01534.pdf,SI p.S16 Eq.(S12)(S14),"Δn=π/β Σ[λ(iωn−iωm)−μ*θ(ωc−|ωm|)]/√(ωm²Zm²+Δm²)·Δm...We used ωc=6 Ry with the appropriate number of the Matsubara frequencies"
4,YH6,note,NA,standard_range_citation,"未找到专门引用。'standard range'句中唯一的引用编号(31)挂在'isotropic Migdal-Eliashberg equations'上（Eliashberg 1959）",NA,P1/2021_Ivan_A_Troyan_1908.01534.pdf,正文 p.7,"Numerical solution of the isotropic Migdal–Eliashberg equations31 within the standard range of μ*(cutoff frequency is 6 Ry) = 0.15–0.1 yields a Tc = 261-272K"
4,YH6,H,165,formula_match_confirm,"Eq.(S2) 与 P1 推导稿 §3.1 闭式逐项相同",NA,P1/2021_Ivan_A_Troyan_1908.01534.pdf,SI p.S15 Eq.(S2),"β_McM=−dlnTc/dlnM=(1/2)[1−1.04(1+λ)(1+0.62λ)/[λ−μ*(1+0.62λ)]²·μ*²]"
4,"YH6(理论预测D)",D,165,Tc_theory_YDx,"168–176K（μ*=0.15–0.1）或 168K（μ*=0.195）——理论预测值",μ*(6Ry)多值,P1/2021_Ivan_A_Troyan_1908.01534.pdf,SI Table S7,"Tc(YDx), K: 168–176"
```

## 正文要点

**① Drozdov 2019**：LaD10 168 K 相的结构是"假设"（"we can naturally suppose"），Ext. Data 明说未测；α=0.57 是 H 最大 Tc 配 D 最大 Tc，非同压；α=0.62 那对（P4/nmm）两侧 XRD 均确认。Ext. Data Fig. 4b 的 LaD10：约 163–165 K @139、158–160 @145、150 @152 GPa（读图）。正文 168 vs Ext. Data 160 内部不一致。

**② Villa-Cortés 2022**：μ* 逐压力点拟合实验 Tc（Fig. 7：162 GPa 约 0.15–0.25，250 GPa 约 0.10–0.12）。α_el-ph 0.45→0.42（N(0)→N(ε)）@162，0.41→0.39 @250；α_el-el −0.07→−0.09 @162，−0.12→−0.13 @250。N(0) 版总 α 在 250 GPa 高估 24%。逐点值只在图里。

**③ Errea 2020**：Mc/AD 列 μ*=0.1 固定；ME^ani 列用 k 平均静态 RPA 屏蔽库仑，无 μ*；SCDFT 列用能量依赖 μ(ε,ε')，无 μ*。**正文 p.4 明确写"in agreement with the experimentally reported α of 0.46"**——所以 0.46 的上游是 Errea 2020 这句话，不是 Durajski/Villa-Cortés 的独立转引错误；但该句附近没有引用编号。researcher 猜测 Errea 看到的是 Drozdov 预印本早期版本。（**主会话补核**：读了 arXiv:1812.01561v1 全文 p.3，v1 与 Nature 终版一样写 0.57/0.62、LaD10 168 K；所以 Errea 2020 的 0.46 没有可识别的一手来源。）

**④ Troyan 2021**：SI Eq. (S12)/(S14) 库仑核 μ*θ(ω_c−|ω_m|)，ω_c = 6 Ry 盒式截止；无 μ* 与裸 μ 关系式；"standard range 0.10–0.15"无引用。SI Eq. (S2) 与 P1 §3.1 闭式逐项相同。Table S7 理论 Tc(YD6) 168–176 K。
