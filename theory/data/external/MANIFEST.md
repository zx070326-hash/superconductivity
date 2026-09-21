# theory/data/external · 外部数据清单（不入 git，见 .gitignore）

| 目录 | 来源 | 许可 | 内容 | 获取日期 |
|---|---|---|---|---|
| `mishra_h3s_d3s/` | Mishra et al., "Nonadiabatic and Anharmonic Effects in High-Pressure H3S and D3S Superconductors", Figshare doi:10.6084/m9.figshare.31136635.v1 = **arXiv:2510.25902 / Annalen der Physik 2026 的配套数据**（reviewer 逐 Tc 核对一致；本地 PDF `library/pdf/P1/2026_Shashi_B_Mishra_2510.25902.pdf`）。**arXiv:2608.14949 表 I 是另一次独立计算**（ML/SSCHA 声子 λ_anh=1.49、μ*=0.13、量子修正压力 210 GPa），不要与本数据集混用 | CC BY 4.0 | `Data-H3S-D3S-PH-EPW-Anharmonic.tar.xz`（24 MB）：H3S-200GPa、D3S-200GPa（谐波 + 非谐 EPW）、H3S-160GPa、D3S-160GPa（非谐）。每个 EPW 目录含 `*.a2f`（ω[meV], α²F, ∫2α²F/ω）、`*.a2fv`、部分含 `*.a2f_proj`（模投影）、`*.dos`（EPW Fermi 窗 DOS）；子目录 `fsr_adia / fsr_na / fbw_adia / fbw_na`（费米面限制/全带宽 × 绝热/非绝热顶角）含 `*.imag_iso_gap0`（能隙边 vs T，meV，网格 1–2 K；Tc 取最后非零温度，√ 拟合外推 RMS 0.54 K）；`muc_0.12/` 为 μ*=0.12 的重复；`band/H3S.dos` 为 QE 全电子 DOS（E_F=17.619 eV）；`phonon/`、`Anharmonic-phonon/`（ZG 迭代）含声子输入输出 | 2026-09-08 |
| `ding_yh6_yh9/` | Ding, Chen, Shi, arXiv:2507.03383，GitHub `dingyucheng-2000/yttrium_hydride_data`（Zenodo 10.5281/zenodo.17472986） | CC BY 4.0 | `EPC_parameters.mat`：YH6/YH9 的 DFPT / SPIA / SSCHA 的 λ(ν)? 9 个数（Fig. 2）；`Tc_renormalized_mustar.mat`：Tc 随某扫描参数（100 点）有/无 μ* 重整化（Fig. 4）；`k_grid_convergence.mat`：λ 对 k 网格（Fig. 5）。**无 α²F(ω) 原始谱**；`PIMD_*.mat`（83 MB）未下载 | 2026-09-08 |
| `cerqueira_raw/`（原始解压，几十 GB，未清理）+ `cerqueira6912_batchq_sample.csv`（解析后，249 行） | Cerqueira, Sanna, Marques, arXiv:2307.10728，Materials Cloud `archive.materialscloud.org/records/3kbt5-r3n56`（doi:10.24435/materialscloud:qv-bq） | CC-BY-4.0 | 17 个 `batch-{a..q}.tar.bz2`（共 6912 材料，51.6 GB）中**只取了 `batch-q`**（250 材料，最小的一包，1.13 GB——纯因体积选中，非按材料属性筛）。每个材料目录的 `McMillan.dat` 文件头直接印出 `lambda`/`wlog[K,meV]`/`w2av[K,meV]`（真实 ⟨ω²⟩^(1/2)，**不需要**手动积分 `a2F.dos`——更正 `theory/notes/materials/omega2_data_sources_2026-09-17.md` Q2(a) 的旧结论）；`mustar-TcMcMillan-TcAllenDynes` 网格步长 0.02；`Eliashberg.dat` 同网格单列。CSV 249 行汇总，解析脚本 `theory/notes/_scratch_ir/cerqueira6912_parse_batch.py`。**⚠ batch-q 全部 249 材料 λ∈[0.02,0.21]，不含论文 Tc>1K 的 2239 个材料中任何一个**，见 `theory/notes/materials/cerqueira6912_batchq_pruning_population_2026-09-18.md` §3。其余 16 个 batch 未取（本机带宽 0.25–0.55 MB/s，全量预计数十小时） | 2026-09-18 |

## 未取得
- Errea 2015 (H3S) 的谱：只有图（本地 PDF Fig. 2）——但 Mishra 数据集覆盖同一体系同一压力，可替代
- Errea 2020 (LaH10/LaD10) 的谱：Nature Source Data 未确认；**LaH10 是当前靶标，最需要**
- Troyan 2021 (YH6) SI、CaH6、PdH/PdD、LaBeH8：无数值文件
- IsoME.jl 声称打包的 H3S/YH6/LaBeH8 benchmark 谱：repo 里只见 Nb

## 使用注意
- Mishra 的 EPW 谐波 H3S 200 GPa：λ=2.29、ω_log=107 meV（Errea 2015 为 2.64、90.4——不同代码/网格/泛函，λ 差 13%）。两套数据不要混用。
- `.a2f` 文件的 ω 网格从 0.48 meV 起，Fermi 窗 0.2 eV，声子 smearing 2.5 meV。
