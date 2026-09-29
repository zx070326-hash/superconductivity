#!/usr/bin/env python3
"""
Cerqueira 6912 批处理：Eliashberg 求解 + Carbotte-Leavens 面积界验证

修正记录（2026-09-21）：
- 改用 dos6（σ=0.03 Ry，作者选定的展宽）
- 修正两处量纲 bug（area 和 moment_bound 的单位转换）
- 所有积分在 Ry 单位下进行，与用户复核脚本一致

输出列：
- material: 材料名
- lambda: 电声耦合常数（无量纲）
- omega2_mev: 二阶矩 √⟨ω²⟩ (meV)
- area_ry: ∫α²F dω (Ry 单位)
- M1_ry: ∫α²F·ω dω (Ry 单位)
- Tc_real_K: 真实谱 Tc (K)
- Tc_einstein_K: Einstein 谱 Tc (K)
- Tc_ceiling_K: 矩界 C_∞√(2M₁) (K)
- Tc_carbotte_K: Carbotte-Leavens 面积界 c(0)·A (K)
- Tc_moment_K: 独立验证的矩界（同 Tc_ceiling_K）
- f_deficit: 𝔣 = Tc/(C_∞√λ·ω̄₂)
- f_ratio: 𝔣/𝔣_E (定理 2 判据)
- carbotte_tighter: 面积界是否更紧
- exceeds_carbotte: Tc 是否超出面积界
- exceeds_moment: Tc 是否超出矩界
"""
import os
import sys
import csv
import importlib.util
from pathlib import Path

# 导入求解器
spec = importlib.util.spec_from_file_location(
    'solver',
    os.path.join(os.path.dirname(__file__), '_cerqueira_eliashberg_solver_2026-09-21.py')
)
solver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(solver)

def collect_materials(base_dir, batches):
    """
    收集所有材料目录
    """
    materials = []
    for batch in batches:
        batch_dir = os.path.join(base_dir, batch)
        if not os.path.exists(batch_dir):
            print(f"WARNING: {batch_dir} not found")
            continue

        for entry in os.listdir(batch_dir):
            mat_dir = os.path.join(batch_dir, entry)
            if os.path.isdir(mat_dir):
                materials.append(mat_dir)

    return materials

def main():
    if len(sys.argv) > 1:
        base_dir = sys.argv[1]
    else:
        base_dir = r'D:\physics\superconductivity\theory\data\external\cerqueira_raw\extracted'

    if len(sys.argv) > 2:
        batches = sys.argv[2].split(',')
    else:
        batches = ['batch-q']  # Default to batch-q for local testing

    materials = collect_materials(base_dir, batches)
    print(f"Found {len(materials)} materials in batches: {batches}")

    # 输出 CSV
    output_file = 'theory/data/external/cerqueira_eliashberg_results.csv'

    fieldnames = [
        'material', 'lambda', 'omega2_mev', 'area_ry', 'M1_ry',
        'Tc_real_K', 'Tc_einstein_K', 'Tc_ceiling_K', 'Tc_carbotte_K', 'Tc_moment_K',
        'f_deficit', 'f_ratio', 'f_carbotte',
        'carbotte_tighter', 'exceeds_carbotte', 'exceeds_ceiling', 'exceeds_moment'
    ]

    count_success = 0
    count_fail = 0
    n_carbotte_tighter = 0
    n_exceeds_carbotte = 0
    n_exceeds_moment = 0

    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i, mat_dir in enumerate(materials, 1):
            mat_name = os.path.basename(mat_dir)

            if i % 50 == 0:
                print(f"Progress: {i}/{len(materials)} ({count_success} success, {count_fail} fail)")

            result = solver.process_material(mat_dir, smearing_idx=6, N=400)

            if result is None:
                count_fail += 1
                continue

            row = {
                'material': result['material'],
                'lambda': f"{result['lambda']:.6f}",
                'omega2_mev': f"{result['omega2_mev']:.3f}",
                'area_ry': f"{result['area_ry']:.6f}",
                'M1_ry': f"{result['M1_ry']:.6f}",
                'Tc_real_K': f"{result['Tc_real_K']:.4f}",
                'Tc_einstein_K': f"{result['Tc_einstein_K']:.4f}",
                'Tc_ceiling_K': f"{result['Tc_ceiling_K']:.4f}",
                'Tc_carbotte_K': f"{result['Tc_carbotte_K']:.4f}",
                'Tc_moment_K': f"{result['Tc_moment_K']:.4f}",
                'f_deficit': f"{result['f_deficit']:.6f}",
                'f_ratio': f"{result['f_ratio']:.6f}",
                'f_carbotte': f"{result['f_carbotte']:.6f}",
                'carbotte_tighter': result['carbotte_tighter'],
                'exceeds_carbotte': result['exceeds_carbotte'],
                'exceeds_ceiling': result['exceeds_ceiling'],
                'exceeds_moment': result['exceeds_moment']
            }

            writer.writerow(row)
            count_success += 1

            if result['carbotte_tighter']:
                n_carbotte_tighter += 1
            if result['exceeds_carbotte']:
                n_exceeds_carbotte += 1
            if result['exceeds_moment']:
                n_exceeds_moment += 1

    print(f"\nDone: {count_success} success, {count_fail} fail")
    print(f"Output: {output_file}")
    print(f"\nBound verification (dos6, σ=0.03 Ry):")
    print(f"  Carbotte tighter than moment: {n_carbotte_tighter}/{count_success}")
    print(f"  Tc exceeds Carbotte bound (MUST be 0): {n_exceeds_carbotte}/{count_success}")
    print(f"  Tc exceeds moment bound (MUST be 0): {n_exceeds_moment}/{count_success}")
    print(f"  Tc exceeds ceiling (MUST be 0): {n_exceeds_ceiling}/{count_success}")

    if violations:
        print(f"\nVIOLATIONS of Theorem 2 ({len(violations)}):")
        for mat, ratio in violations[:10]:
            print(f"  {mat}: f_ratio={ratio:.6f}")
        if len(violations) > 10:
            print(f"  ... and {len(violations)-10} more")
    else:
        print("\nNo violations of Theorem 2 found")

if __name__ == '__main__':
    main()
