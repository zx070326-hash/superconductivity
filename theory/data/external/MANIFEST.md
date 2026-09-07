# theory/data/external · 外部数据清单（不入 git，见 .gitignore）

| 目录 | 来源 | 许可 | 内容 | 获取日期 |
|---|---|---|---|---|
| `mishra_h3s_d3s/` | Mishra et al., "Nonadiabatic and Anharmonic Effects in High-Pressure H3S and D3S Superconductors", Figshare doi:10.6084/m9.figshare.31136635.v1（arXiv:2510.25902 / Annalen der Physik 2026 的配套数据；arXiv:2608.14949 的 Data Availability 也引用它） | CC BY 4.0 | `Data-H3S-D3S-PH-EPW-Anharmonic.tar.xz`（24 MB）：H3S-200GPa、D3S-200GPa（谐波 + 非谐 EPW）、H3S-160GPa、D3S-160GPa（非谐）。每个 EPW 目录含 `*.a2f`（ω[meV], α²F, ∫2α²F/ω）、`*.a2fv`、部分含 `*.a2f_proj`（模投影）、`*.dos`（EPW Fermi 窗 DOS）；子目录 `fsr_adia / fsr_na / fbw_adia / fbw_na`（费米面限制/全带宽 × 绝热/非绝热顶角）含 `*.imag_iso_gap0_T` 能隙随温度文件；`muc_0.12/` 为 μ*=0.12 的重复；`band/H3S.dos` 为 QE 全电子 DOS（E_F=17.619 eV）；`phonon/`、`Anharmonic-phonon/`（ZG 迭代）含声子输入输出 | 2026-09-08 |
| `ding_yh6_yh9/` | Ding, Chen, Shi, arXiv:2507.03383，GitHub `dingyucheng-2000/yttrium_hydride_data`（Zenodo 10.5281/zenodo.17472986） | CC BY 4.0 | `EPC_parameters.mat`：YH6/YH9 的 DFPT / SPIA / SSCHA 的 λ(ν)? 9 个数（Fig. 2）；`Tc_renormalized_mustar.mat`：Tc 随某扫描参数（100 点）有/无 μ* 重整化（Fig. 4）；`k_grid_convergence.mat`：λ 对 k 网格（Fig. 5）。**无 α²F(ω) 原始谱**；`PIMD_*.mat`（83 MB）未下载 | 2026-09-08 |

## 未取得
- Errea 2015 (H3S) 的谱：只有图（本地 PDF Fig. 2）——但 Mishra 数据集覆盖同一体系同一压力，可替代
- Errea 2020 (LaH10/LaD10) 的谱：Nature Source Data 未确认；**LaH10 是当前靶标，最需要**
- Troyan 2021 (YH6) SI、CaH6、PdH/PdD、LaBeH8：无数值文件
- IsoME.jl 声称打包的 H3S/YH6/LaBeH8 benchmark 谱：repo 里只见 Nb

## 使用注意
- Mishra 的 EPW 谐波 H3S 200 GPa：λ=2.29、ω_log=107 meV（Errea 2015 为 2.64、90.4——不同代码/网格/泛函，λ 差 13%）。两套数据不要混用。
- `.a2f` 文件的 ω 网格从 0.48 meV 起，Fermi 窗 0.2 eV，声子 smearing 2.5 meV。
