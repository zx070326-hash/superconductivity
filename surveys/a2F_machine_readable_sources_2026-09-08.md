# 机器可读 α²F(ω) 数据源调研

调研 agent：a2f-data-hunt（sonnet / researcher）　日期：2026-09-08。原始产出，未手改。WebSearch 配额中途用尽，后半段只用 WebFetch；"verified_exists=no"多为 403 拦截而非确认不存在。

```csv
system,isotope,pressure_GPa,harmonic_or_anharmonic,source_type,url,file_name,format_notes,license,verified_exists,notes
H3S/D3S,H+D,~200,both(likely),dataset,https://doi.org/10.6084/m9.figshare.31136635.v1,unknown(403 on file list),unknown - need login to see files,unknown,partial,"DOI来自arXiv:2608.14949 Data Availability引用；重定向后页面标题是'Nonadiabatic and Anharmonic Effects in High-Pressure H3S and D3S Superconductors'，怀疑对应姊妹论文arXiv:2510.25902；优先级最高但文件列表未核实"
H3S,H,unknown,likely anharmonic(SSCHA iteration),tutorial,https://github.com/SSCHAcode/python-sscha/tree/master/Tutorials/H3S,"all_frequencies.dat, get_phonons.py, matdyn0, matdyn1, matdyn2, matdyn3",matdyn0-3疑似SSCHA自洽迭代不同代的声子谱输出,GPL v3(repo整体),yes,"教学示例，压强/k网格参数不一定对应已发表200GPa结果"
H3S,--,200,figure only,paper_figure,https://arxiv.org/abs/1502.02832,Fig.2 (PRL 114 157004),仅PDF图,N/A,no,需数字化图2曲线
YH6/YH9,H,unknown,anharmonic(SSCHA)+mustar renorm,dataset,https://github.com/dingyucheng-2000/yttrium_hydride_data,"EPC_parameters.mat(Fig.2), Tc_renormalized_mustar.mat(Fig.4), k_grid_convergence.mat(Fig.5), PIMD_RDF.mat+PIMD_centroid_trajectory.mat(Fig.1)",MATLAB .mat；对应arXiv:2507.03383；Zenodo 10.5281/zenodo.17472986,CC-BY 4.0,yes,".mat很可能是导出量(lambda/omega_log/Tc)而非原始alpha2F(omega)曲线"
YH6,H,166,figure only,SI_unconfirmed,https://advanced.onlinelibrary.wiley.com/doi/10.1002/adma.202006832,Supporting Information,Wiley SI页403,unknown,no,对应arXiv:1908.01534
LaH10/LaD10,H,various,harmonic(at least),SI_unconfirmed,https://www.nature.com/articles/s41586-020-1955-z,Source Data (name unknown),Nature页面403,unknown,no,对应arXiv:1907.11916
LaH10,H,unknown,structure relax only(no alpha2F),tutorial,https://github.com/SSCHAcode/python-sscha/tree/master/Tutorials/LaH10,"VariableCellRelaxation.ipynb, dyn1, UPF",只是变胞弛豫教程,GPL v3,yes,"与本任务不相关"
multi-hydride(ambient-pressure-stable only),mixed,~0,harmonic(DFPT),high-throughput,https://archive.materialscloud.org/records/3kbt5-r3n56,"files_description.md, README.txt, batch-a...batch-q.tar.bz2 (17 batches, 48.1GiB)",7000个DFPT计算,CC BY 4.0,yes,"Cerqueira/Sanna/Marques；常压凸包稳定相，大概率不覆盖高压氢化物"
general(≤5-atom primitive cells),mixed,ambient,harmonic,dataset,https://figshare.com/articles/dataset/JARVIS-SuperconDB/21370572,"jarvis_epc_data_figshare_1058.json.zip, jarvis_epc_data_2d.json.zip",列含a2F/wlog/lamb/Tc/press;JSON,unknown,partial,"WebFetch被403拦截"
general,mixed,ambient,harmonic,dataset+code,https://github.com/henniggroup/BETE-NET,"database.json, structures/(CIF), indices/, test_preds/",818个材料alpha2F训练数据,Apache-2.0,yes,"是否含氢化物需下载核实"
H3S/YH6/LaBeH8/LaBH8,mixed,unknown,unknown,code+data(claimed),https://github.com/cheil/IsoME.jl,claimed bundled but not found (only Nb test case visible),Julia;Zenodo 10.5281/zenodo.14967551,unknown,no,"论文声称打包但repo只见Nb"
CaH6/CaD6,H+D,unknown,unknown,none_found,,,,,no,"只有论文图"
PdH/PdD/PdT,H+D+T,ambient,unknown,none_found,,,,,no,"只有论文图"
LaBeH8/LaBeD8,H+D,80,unknown,none_found,,,,,no,"只有论文图"
```

## 主会话后续（2026-09-08）
- Figshare 31136635 经 API 确认单文件 `Data-H3S-D3S-PH-EPW-Anharmonic.tar.xz`（24 MB，CC BY 4.0），已下载解压到 `theory/data/external/mishra_h3s_d3s/`，内容远超预期：含 H3S/D3S 200 GPa 谐波+非谐、160 GPa 非谐的 `.a2f` 与四种处理（fsr/fbw × adia/na）的能隙-温度扫描——可直接做四因子分解，见推导稿 §4.1 v0.6。
- Ding 仓库 4 个小 .mat 已下载；`EPC_parameters.mat` 是 λ 的 9 个 Matsubara 点（DFPT/SPIA/SSCHA），`Tc_renormalized_mustar.mat` 是 Tc 有/无 μ* 重整化随某参数的 100 点扫描；无 α²F(ω)。
- python-sscha H3S 教程未取。
