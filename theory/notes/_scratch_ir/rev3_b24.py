"""rev3 item 4: X5 B2/B4 (log) vs float64 max ||H^(j)(x)|| on J=[0.5w,1.5w], N=1500, 9 points (self-written H)."""
import numpy as np, sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from rev2_truth_lib import H
for w, B2, B4 in ((30.0, 0.0488, 0.0099), (40.0, 0.0360, 0.0042), (60.0, 0.0237, 0.0013)):
    xs = np.linspace(0.5*w, 1.5*w, 9)
    m2 = max(np.max(np.abs(np.linalg.eigvalsh(H(1500, x, 2)))) for x in xs)
    m4 = max(np.max(np.abs(np.linalg.eigvalsh(H(1500, x, 4)))) for x in xs)
    print(f"w={w:>4.0f}: B2={B2} vs true {m2:.5f} ({B2/m2:.1f}x)  B4={B4} vs true {m4:.6f} ({B4/m4:.1f}x)  {'OK' if B2>=m2 and B4>=m4 else 'FATAL'}", flush=True)
