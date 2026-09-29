import sys
import numpy as np
sys.path.insert(0, "D:/physics/superconductivity/.claude/worktrees/agent-a8932ec4c14045b28/theory/derivations")
from p4_concavity_q_route import O_ker, kernels, top_eigen

THETA = 0.6349
print("Pin the hole edges and the worst point on a fine grid.")
print("")
rows = []
for varpi in np.arange(1.30, 1.86, 0.01):
    x = float(varpi) ** 2
    N = 300
    ker = kernels(x, KX=8 * N)
    A = O_ker(N, ker["f"])
    B = O_ker(N, ker["f2"])
    h, v, l2 = top_eigen(A)
    E = v @ (O_ker(N, ker["g"]) @ v)
    cF = np.linalg.norm(A @ B - B @ A, "fro")
    num = (2.0 / x) * (cF / np.sqrt(2.0)) ** 2
    g_lo = max(0.42 * h, h - THETA)
    cb = num / g_lo ** 3
    need = (num / E) ** (1.0 / 3.0)
    rows.append((float(varpi), h, E, cb, cb / E, need, need / h, l2 / h))

bad = [r for r in rows if r[3] >= r[2]]
print(f"hole edges: varpi in [{bad[0][0]:.2f}, {bad[-1][0]:.2f}]"
      if bad else "no hole")
worst = max(rows, key=lambda r: r[4])
print(f"worst point: varpi={worst[0]:.2f}  |C|an/E={worst[4]:.4f} "
      f"(short by {worst[4]:.3f}x)")
print(f"  needed gap there: {worst[5]:.4f} = {worst[6]:.4f}*h "
      f"(have 0.42h);  true lam2/h={worst[7]:.2e}")
print("")
print("varpi    h       E       |C|an   ratio   need_g/h")
for r in rows[::4]:
    flag = "  HOLE" if r[3] >= r[2] else ""
    print(f"{r[0]:5.2f} {r[1]:7.4f} {r[2]:.5f} {r[3]:.5f} {r[4]:6.3f} "
          f"{r[6]:8.4f}{flag}")
print("")
print("max needed gap/h over the hole:",
      f"{max(r[6] for r in bad):.4f}" if bad else "n/a")
