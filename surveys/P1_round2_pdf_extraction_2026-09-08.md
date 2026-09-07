# P1 二轮支撑：本地 PDF 原文抽取（H/D 数据）+ 公式核对 + 新颖性复核

调研 agent：p1-data-round2（sonnet / researcher）　日期：2026-09-08
方法：Read 本地 PDF（页码 + 英文原句），WebSearch 为辅。原始产出，未手改（两条回传合并）。

---

## 任务 A：数据表（逐条带页码 + 原文引文）

```csv
material,isotope,pressure_GPa,quantity,value,mu_star_convention,source_file,page,quote_or_note
LaH10,H,~170,Tc_exp,250-252,NA,P1/2019_A_P_Drozdov_1812.01561.pdf,p.1摘要/p.4正文,"the pressure dependence... reaching the maximum Tc=250-252K at ~170 GPa"
LaH10,H,150,Tc_exp,249,NA,P1/2019_A_P_Drozdov_1812.01561.pdf,Fig.1图注,"blue curve corresponds to a sample heated under 150 GPa with Tc ~249K"
LaH10,H,NA(P4/nmm competing phase),Tc_exp,215,NA,P1/2019_A_P_Drozdov_1812.01561.pdf,正文+ExtData Fig.5图注,"there are only four distinct, reproducible superconducting transitions in lanthanum hydrides with Tc~250K, ~215K, ~110K, and ~70K"
LaD10.6,D,142,Tc_exp,140,NA,P1/2019_A_P_Drozdov_1812.01561.pdf,正文(Fig.3c附近),"the La+D2 mixture with a Tc~140K at 142 GPa"; P4/nmm结构, 化学计量 LaD10.6
LaD10,D,NA,Tc_exp,"168(正文) vs 160(ExtData Fig.5图注)",NA,P1/2019_A_P_Drozdov_1812.01561.pdf,正文 vs Ext.Data Fig.5,"正文: hydride phases displaying Tc of ~168K and ~140K"; "Ext.Data Fig.5图注: three transitions in lanthanum deuterides with Tc ~160K, 140K and ~70K"——同一篇文章内部160/168小矛盾
LaH10/LaD10,H&D,~170(fcc)/~140(P4nmm),alpha_exp,"0.57(fcc相)/0.62(P4/nmm相)",NA,P1/2019_A_P_Drozdov_1812.01561.pdf,正文,"we obtained the isotope coefficient α of about 0.57 and 0.62 for the superconductivity in the fcc and P4/nmm phases, respectively... α could be explained probably by a soft anharmonism of H/D vibrations"
LaH10/LaD10,note,NA,Tc_exp/alpha_exp,NOT_FOUND,NA,P1/2019_A_P_Drozdov_1812.01561.pdf,全文(正文+Ext.Data)已通读两遍,原文全篇没有任何"Tc(LaD10)=180K"或"α=0.46/0.47"的数字——这是被广泛转引的两个数字,但一手实验论文里查无此数
LaH10,H,129,"lambda_anh,wlog_anh_meV,Tc_ME_anh",3.62;76.4;255.3,μ*=0.1(Tc_Mc/Tc_AD两列),D3/2020_Ion_Errea_1907.11916.pdf,p.11 Extended Data Table I,"Summary of calculated superconducting Tc... anisotropic treatment of Migdal-Eliashberg (Tc_ani^ME)"；λ/ωlog列本身就是SSCHA非谐值,该文未给独立谐波表
LaH10,H,163,"lambda_anh,wlog_anh_meV,Tc_ME_anh",2.67;96.4;242.8,同上,D3/2020_Ion_Errea_1907.11916.pdf,p.11 Ext.Data Table I,同一表格行
LaH10,H,214,"lambda_anh,wlog_anh_meV,Tc_ME_anh",2.06;115.5;237.9,同上,D3/2020_Ion_Errea_1907.11916.pdf,p.11 Ext.Data Table I,同一表格行
LaH10,H,264,"lambda_anh,wlog_anh_meV,Tc_ME_anh",1.73;126.6;216.9,同上,D3/2020_Ion_Errea_1907.11916.pdf,p.11 Ext.Data Table I,同一表格行
LaD10,D,159,"lambda_anh,wlog_anh_meV,Tc_ME_anh",3.14;63.5;180.4,同上,D3/2020_Ion_Errea_1907.11916.pdf,p.11 Ext.Data Table I,同一表格行——注意159GPa与LaH10最近点163GPa相差4GPa,非严格同压
LaD10,D,210,"lambda_anh,wlog_anh_meV,Tc_ME_anh",2.21;81.7;172.9,同上,D3/2020_Ion_Errea_1907.11916.pdf,p.11 Ext.Data Table I,同一表格行
LaD10,D,260,"lambda_anh,wlog_anh_meV,Tc_ME_anh",1.80;92.2;157.9,同上,D3/2020_Ion_Errea_1907.11916.pdf,p.11 Ext.Data Table I,同一表格行
LaH10/LaD10,note,163/159(最近点配对),alpha_theory,"0.43(用Tc_ME_anh) / 0.42(用Tc_AD) / 0.40(用Tc_SCDFT)",同表4种μ*/方法,D3/2020_Ion_Errea_1907.11916.pdf,由p.11表计算,这是我自己用表中数字算的(researcher推断,非原文直接给出的α),ln(242.8/180.4)/ln2=0.43；**与被转引的"0.46/0.47"接近但不完全相等**——不是exact match
LaH10/LaD10,note,~150(Table II均值代表),alpha_theory,0.465(平均),Durajski自己的ME求解;引用ref[23]=Drozdov2019做压力/结构输入,P1/2020_A_P_Durajski_1912.09663.pdf,Table II,"ᾱ=0.465"——这是Durajski自己算的理论α(把H的α²F按ω∝M^-1/2标度近似算D),不是从Drozdov2019直接读出的实验数,但正文把此结果与ref[23]放在一起表述,造成"这是Drozdov实验值"的误读
LaH10/LaD10,note,NA,cite_check,confirmed,NA,P1/2020_A_P_Durajski_1912.09663.pdf,参考文献表,ref[23] = Drozdov et al. Nature 569,528(2019) = arXiv:1812.01561——文号确认无误,但该文不含Durajski引用意图对应的180K/0.46数字
H3S,H,NA,Tc_exp,203,NA(见下),P1/2015_A_P_Drozdov_1506.08190.pdf,p.3,"Tonset=203(1)K"(磁化率测量)；摘要"confirmed by magnetic susceptibility measurements with Tc=203K"
H3S/D3S,H&D,>170(区间拟合,非单点),alpha_exp,0.3,NA,P1/2015_A_P_Drozdov_1506.08190.pdf,p.3,"The BCS theory gives dependence of critical temperature on the atomic mass M as Tc~M^-α where αH~0.5. Comparison Tc's in the pressure range P>170 GPa (Fig. 2c) gives αH~0.3."
H3S/D3S,note,NA,alpha_theory(anharmonic),0.35,NA(引自Errea PRL114,157004(2015)),P1/2015_A_P_Drozdov_1506.08190.pdf,p.4,"Inclusion of zero point vibrations and anharmonicity in the calculations corrected the calculated Tc to ~190K, and the isotope coefficient from α=0.5 to α=0.35 — both in agreement with the present work."
H3S,H,155,Tc_exp,203,NA,P1/2022_S_Villa-Cort_s_2103.12055.pdf,p.1引言,"experimental results from optical reflectivity... Tc=203K under a hydrostatic pressure of 155 GPa"(转引Drozdov2015,非本文一手测量)
LaH10,note,150,alpha_exp(转引声称),0.47,NA,P1/2022_S_Villa-Cort_s_2103.12055.pdf,p.1引言+参考文献表p.11,"For LaH10 a value of α≅0.47 at 150 GPa has been observed[10]"——ref[10]核对为"Drozdov et al., Nature 569,528(2019)"，即arXiv:1812.01561——**与Drozdov2019正文实际所写的0.57/0.62不符**，这是第二篇独立文献重复同一转引错误
H3S,note(N(ε)可变DOS版本),130,alpha_theory,2.37,N(ε)能量依赖(含vHs),P1/2022_S_Villa-Cort_s_2103.12055.pdf,p.1摘要,"in H3S it has a very large value of α≅2.37 at 130 GPa, decreasing it rapidly as pressure increases up to 230 GPa (α≅0.29)"
H3S,note(N(ε)可变DOS版本),230,alpha_theory,0.29,N(ε)能量依赖(含vHs),P1/2022_S_Villa-Cort_s_2103.12055.pdf,p.1摘要,同上
H3S,note,NA,alpha_theory(constant-DOS N(0)版本数值),NOT_EXTRACTED,N(0)=const,P1/2022_S_Villa-Cort_s_2103.12055.pdf,仅读了p.1-2/p.10-12,摘要定性说"对比N(ε)=N(0)=const...单调递减"但具体数值表在结果section(约p.3-9)本轮未读到,需二次补读
H3S,note,NA,el-ph与el-el分量拆分数值,NOT_EXTRACTED(仅有定性符号),NA,P1/2022_S_Villa-Cort_s_2103.12055.pdf,p.1摘要,"one positive coming from the electron-phonon(el-ph) and the other negative from the electron-electron interaction(el-el)"——具体每压力点的两个分量数值同样在结果section未读到
YH6,H,166,Tc_exp,224,NA,P1/2021_Ivan_A_Troyan_1908.01534.pdf,p.5 Fig.3(a)及正文,"Two slightly different superconducting transitions in the YH6 with Tc of 224K(Figure 3a)"; DAC K1
YD6,D,172,Tc_exp,170,NA,P1/2021_Ivan_A_Troyan_1908.01534.pdf,p.5,"The YD6 sample synthesized using deuterated ammonium borane demonstrates superconducting transition at 170 K"; DAC D1
YH6/YD6,H&D,166 vs 172(仅差6GPa),alpha_exp,0.4,NA,P1/2021_Ivan_A_Troyan_1908.01534.pdf,p.5,"which corresponds to the isotope coefficient αexp = 0.4, lower than the BCS theory gives (~0.5)"——验算: ln(224/170)/ln2=0.399，核对通过
YH6,H,165(计算用固定压力),"lambda_har,wlog_har_meV,Tc_har",2.24;80.06(=929K换算);261-272,μ*(6Ry)=0.10-0.15(标准范围),P1/2021_Ivan_A_Troyan_1908.01534.pdf,p.7-8,"Numerical solution of the isotropic ME equations within the standard range of μ*(cutoff 6Ry)=0.15-0.1 yields Tc=261-272K"；929K是原文单位(K),我换算成meV
YH6,H,165,"lambda_anh,wlog_anh_meV,Tc_anh",1.71;114.87(=1333K换算);236-247,同上μ*范围,P1/2021_Ivan_A_Troyan_1908.01534.pdf,p.7-8,"anharmonic α²F(ω) show a decrease in EPC coefficient from 2.24 to 1.71 and an increase of ωlog from 929 to 1333K...resulting Tc is 236-247K"
YH6,H,165,Tc_matched_to_exp,224(拟合值),μ*(6Ry)=0.19-0.22(异常大,非标准范围),P1/2021_Ivan_A_Troyan_1908.01534.pdf,p.8,"To match the experimental data, including isotope coefficient αexp=0.4, unusually high values of the Coulomb pseudopotential μ*(6Ry)=0.19–0.22 are necessary."
YH6,H,165,Tc_SCDFT_anh,160,μ=0.187(SCDFT无经验μ*的库仑势),P1/2021_Ivan_A_Troyan_1908.01534.pdf,p.8,"Solving the SCDFT gap equation with the anharmonic α²F(ω)...at 165 GPa yields Tc=160K and Coulomb potential μ=0.187"
YD6,note,NA,SSCHA_lambda_wlog_Tc,NOT_FOUND,NA,P1/2021_Ivan_A_Troyan_1908.01534.pdf,SI Tables S6-S9(前次session已核),SI中只有YH6/YH4/YH7的SSCHA/ME表，**没有YD6的SSCHA计算**——本文α_exp=0.4是纯实验值,理论侧只对H做了非谐计算
YH9,H,201,Tc_exp,243,NA,P1/2021_P_P_Kong_1909.10482.pdf,p.1摘要/p.4,"the highest Tc of 243 K at 201 GPa"; "the P6₃/mmc YH9 phase with a Tc of 243 K"
YH6,H,237,Tc_exp,227,NA,P1/2021_P_P_Kong_1909.10482.pdf,p.1摘要/p.4,"we also observed a Tc of 227 K at 237 GPa for the YH6 phase"; sample1
YH6,H,201,Tc_exp,210,NA,P1/2021_P_P_Kong_1909.10482.pdf,p.4,"a Tc of 210K in the Im-3m YH6 phase was observed"(sample2,laser加热前)
YD6,D,170,Tc_exp,170,NA,P1/2021_P_P_Kong_1909.10482.pdf,p.4,"We observed a superconductivity with a Tc of 170K in sample 6 which was synthesized through the pulsed laser heating of YD3 under deuterium atmosphere...at a pressure of 170GPa"
YD9,note,NA,Tc_exp,NOT_FOUND,NA,P1/2021_P_P_Kong_1909.10482.pdf,全文p.1-9通读,**本文没有YD9数据**——同位素替代实验只做了YH6/YD6一对;结构部分只提到YH3/YD3、YH6/YD6
YH6/YD6,note,237(H) vs 170(D)——压差67GPa,alpha_naive,0.42(researcher自算,低置信度),NA,P1/2021_P_P_Kong_1909.10482.pdf,p.4数字组合计算,原文**没有给出α数值**;压力严重不匹配(67GPa),不建议采信
PdH,H,常压,"lambda_har,wlog_har_meV,Tc_har",1.55;25.42(=205cm-1换算);47,μ*=0.085(单带ME,来自ref[18]),D3/2013_Ion_Errea_1305.7123.pdf,p.4 Table I,"TABLE I. Calculated λ, ωlog and Tc values for PH...harmonic approximation strongly overestimates Tc"
PdD,D,常压,"lambda_har,wlog_har_meV,Tc_har",1.55;18.60(=150cm-1);34,同上,D3/2013_Ion_Errea_1305.7123.pdf,p.4 Table I,同一表格行
PdT,T,常压,"lambda_har,wlog_har_meV,Tc_har",1.55;15.50(=125cm-1);30,同上,D3/2013_Ion_Errea_1305.7123.pdf,p.4 Table I,同一表格行
PdH,H,常压,"lambda_anh,wlog_anh_meV,Tc_anh",0.40;50.21(=405cm-1);5.0,同上μ*,D3/2013_Ion_Errea_1305.7123.pdf,p.4 Table I,SSCHA列
PdD,D,常压,"lambda_anh,wlog_anh_meV,Tc_anh",0.46;37.69(=304cm-1);6.5,同上μ*,D3/2013_Ion_Errea_1305.7123.pdf,p.4 Table I,SSCHA列
PdT,T,常压,"lambda_anh,wlog_anh_meV,Tc_anh",0.48;31.86(=257cm-1);6.9,同上μ*,D3/2013_Ion_Errea_1305.7123.pdf,p.4 Table I,SSCHA列
PdH,H,常压,Tc_exp,"9[ref7]/8[ref8]",NA,D3/2013_Ion_Errea_1305.7123.pdf,p.4 Table I,ref7=Schirber&Northrup PRB10,3818(1974); ref8=Schirber/Mintz/Wall SSC52,837(1984)
PdD,D,常压,Tc_exp,"11[ref7]/10[ref8]",NA,D3/2013_Ion_Errea_1305.7123.pdf,p.4 Table I,同上两篇实验源
PdH/PdD,note,常压,"alpha_har/alpha_anh/alpha_exp",0.47 / -0.38 / "-0.29[7],-0.32[8]",μ*=0.085 fixed,D3/2013_Ion_Errea_1305.7123.pdf,p.4 Table I底部,"α_PdH(D): harmonic=0.47, SSCHA=-0.38, Exp=-0.29[7]/-0.32[8]"
PdH/PdT,note,常压,"alpha_har/alpha_anh",0.41 / -0.29,μ*=0.085 fixed,D3/2013_Ion_Errea_1305.7123.pdf,p.4 Table I底部,同上表格(无PdT实验α给出)
```

数据缺口（未完全覆盖）：
1. Villa-Cortés 2022 的结果 section（约 p.3–9）未读，N(0)=const 版本的 α(P) 数值表、el-ph/el-el 分量逐压力点数值未抽取。
2. Kong 2021 只读了 p.1–9。
3. Villa-Cortés 的 μ* "fixed vs cutoff-varying" 约定只从附录 A 方程组间接看到 μ* 以固定值出现，置信度中等。

## 任务 B-1：Carbotte 1990 / Garland 1963 的 μ* 修正公式

没有找到 Carbotte RMP 62, 1027 或 Garland PRL 11 的可访问全文（付费墙）。间接证据链：Troyan 2021 (arXiv:1908.01534) 的 SI 方程 (S2) 与推导稿 §3.1 的闭式 α = ½[1 − 1.04(1+λ)(1+0.62λ)μ*²/(λ − μ*(1+0.62λ))²] 逐项一致，且两者都引用 Allen & Dynes PRB 12, 905 (1975)。**推断**：公式通过 Allen–Dynes 原始文献 + Troyan SI 的独立复现得到间接确认，但 Carbotte 1990 原文逐字表述仍未见到。

## 任务 B-2：Rainer & Culetto (1979) / Bergmann & Rainer (1973)

从 Villa-Cortés 2022 正文直接拿到完整方程链。原文："Rainer and Culetto have shown that ΔTc^el-ph and their corresponding total isotope coefficient α^el-ph can be calculated as"：
- ΔTc^(el-ph) = ∫₀^∞ dω [δTc/δα²F(ω)] Δα²F(ω) ——式(3)
- α^(el-ph) = ∫₀^∞ dω β^(el-ph)(ω) ——式(4)
- β^(el-ph)(ω) ≡ −[δTc/δα²F(ω)]·[Δα²F(ω)/(Tc Δln M)]，"differential isotope coefficient" ——式(5)
- β^(el-ph)(ω) ≡ R(ω)α²F(ω) ——式(6)
- R(ω) = d/dω [ω/(2Tc) · δTc/δα²F(ω)] ——式(7)

ref 52 = Rainer & Culetto PRB 19, 2540 (1979)，精确匹配。Bergmann & Rainer 1973 = ref 46（Z. Phys. 263, 59），存在但对应方程本轮未核实。**Σ_i α_i = ½ 求和规则**：两轮检索都没找到一手陈述，维持"检索摘要级"。

## 任务 C：新颖性复核

维持不变：截至 2026-09-08 没有找到把库仑 μ*(M) + 非谐 + 非绝热顶角 + 能量依赖 DOS 四项同时纳入同一氢化物同位素系数分解的论文。最接近：Mishra/Belli/Zurek/Margine arXiv:2608.14949（非谐+非绝热）；Villa-Cortés 2022（库仑+DOS）。建议正式投稿前做系统性主题词穷举。

## 三个关键矛盾点

**矛盾①：H3S 的 α_H。** 推导稿 §4.1 的 0.44 是两点比值；Drozdov 2015 p.3 原句用 P>170 GPa 全区间比较得 α_H ~ 0.3。改用 0.3 后账目几乎清零，"需要非绝热/DOS 填补差额"的证据基础被削弱。**原文直接陈述。**

**矛盾②：LaH10/LaD10 的 α。** 推导稿与多篇二手文献写 0.46–0.47，但 Drozdov 2019 全文只有一句："we obtained the isotope coefficient α of about **0.57** and **0.62** for the superconductivity in the **fcc** and **P4/nmm** phases, respectively"——没有 0.46/0.47，也没有 Tc(LaD10)=180 K。Durajski 2020 ref[23] 和 Villa-Cortés 2022 ref[10] 都指向 Drozdov 2019 却都转引出 0.46/0.47。Errea 2020 SSCHA 表反算 α≈0.43、Tc(LaD10)=180.4 K——"180 K"很可能是理论值被误传（researcher 假说，中等置信度）。真实实验值 0.57/0.62 **远高于** BCS 理想值 0.5，与 H3S、YH6 的"α<0.5"方向相反。

**矛盾③（次要）：** Drozdov 2019 内部 168 vs 160 K 不一致。

## 补充发现

Kong 2021 **没有 YD9 数据**，同位素实验只做了 YH6/YD6（170 K @170 GPa vs YH6 227 K @237 GPa，压差 67 GPa，作者未算 α）。Troyan 2021 是 YH6 体系唯一"干净"的实验 α 来源。
