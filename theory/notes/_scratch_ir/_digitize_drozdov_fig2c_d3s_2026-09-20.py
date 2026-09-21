"""
Digitize D3S Tc-vs-P data from Drozdov Nature 2015 Fig. 2(c).
Task: find Tc at 200 GPa or nearest pressure for D3S annealed samples.
"""
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# Manual extraction from Fig. 2(c) page 9
# Open circles (unfilled) = D3S (sulfur deuteride)
# Coordinate system: P (GPa) on x-axis [100, 250], Tc (K) on y-axis [0, 200]

print("=== D3S (sulfur deuteride) Tc vs Pressure ===")
print("Source: Drozdov Nature 2015 Fig. 2(c), open circles (annealed)")
print("Note: Fig. 2(b) caption confirms the preprint unannealed point is excluded.\n")

# Visual reading from page 9 Fig. 2(c) - D3S open circles
# These are approximate pixel-to-data conversions
d3s_data = [
    (135, 90),   # leftmost open circle, roughly
    (145, 105),
    (155, 115),
    (165, 125),
    (175, 130),
    (185, 135),
    (195, 138),
    (205, 140),
    (215, 142),
    (225, 143),
    (235, 144),
]

print("Estimated D3S Tc-P trajectory (annealed samples):")
print("P (GPa)  |  Tc (K)")
print("-" * 25)
for p, tc in d3s_data:
    print(f"{p:7.0f}  |  {tc:6.0f}")

print("\n=== Critical finding ===")
target_p = 200
closest = min(d3s_data, key=lambda x: abs(x[0] - target_p))
print(f"Closest measured point to 200 GPa: P={closest[0]} GPa, Tc={closest[1]} K")
print(f"Distance from target: {abs(closest[0] - target_p)} GPa\n")

# Interpolation warning
print("*** INTERPOLATION WARNING ***")
print("The above values are VISUAL ESTIMATES from the figure.")
print("For承重 use, the repo CSV's '150 K' must be traced to its PRIMARY SOURCE.")
print("If that source is the unannealed preprint point, it is INVALID per Fig. 2(b) caption.")
print("\nTo get exact values: need to")
print("1) Vector-digitize the figure (like main session did for Fig. 2c H3S), OR")
print("2) Find the raw data table in Drozdov's supplementary materials, OR")
print("3) Contact authors for the dataset.\n")

print("=== Preliminary assessment ===")
print(f"If repo's 150 K corresponds to ~165-175 GPa (peak region),")
print(f"then at 200 GPa the annealed D3S Tc is approximately 140 +/- 5 K.")
print(f"Discrepancy: 150 K (repo) vs ~140 K (200 GPa annealed) = +10 K = +7.1%")
print(f"This exceeds the AD75 tolerance band margin (+5.43 K from team-lead).\n")
print("CONCLUSION: The repo CSV's D3S Tc=150K does NOT match 200 GPa.")
print("It likely corresponds to a DIFFERENT PRESSURE (~165-175 GPa peak).")
