"""rev3: consistency of the mutilde=0.33 [5,20] infimum 0.163538 with X1 v2's mutilde=0.26 infimum
0.350203: the alpha side is identical (same grid, same Na banding, same weighted L), so
inf g(0.33) = inf g(0.26) - [mu2(0.26) - mu2(0.33)] on the binding cell if the shift is ~constant."""
import sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
from mpmath import mp
from p4d_coulomb_certify import mu2_pencil_cell
from p4_interval_extension import ext_grid
mp.dps = 30
g = ext_grid(5.0, 20.0)
shifts = []
for i in range(len(g)-1):
    a, b = g[i], g[i+1]
    m26 = float(mu2_pencil_cell(a, b, 0.26)[0]); m33 = float(mu2_pencil_cell(a, b, 0.33)[0])
    shifts.append((m26 - m33, a, b))
print(f"mu2(0.26)-mu2(0.33) over the 121 cells: min {min(shifts)[0]:.6f} at [{min(shifts)[1]:.4f},{min(shifts)[2]:.4f}], max {max(shifts)[0]:.6f}")
print(f"0.350203 - shift range -> [{0.350203-max(shifts)[0]:.6f}, {0.350203-min(shifts)[0]:.6f}]  vs logged inf g(0.33) = 0.163538")
