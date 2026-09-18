"""rev3 item 5 + item 4 truth side: self-written float64 H (rev2_truth_lib), N=3000.
Claim numbers are hard-coded from x4_run2.log (per-cell al= lines, X5 node rows)."""
import numpy as np, time, sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
from rev2_truth_lib import H, lam12, lam1_2x2
from mpmath import mp
from p4_interchange_uniformity import mu2_cell_lb          # claim side (cheap, exact 2x2)
mp.dps = 30
t0 = time.time()
spots = [(20.1, 20.000, 20.115, 340, 0.454), (35.0, 34.906, 35.107, 580, 0.515),
         (59.9, 59.883, 60.000, 960, 0.504)]
print("== item 5: spot w vs cell alpha_cell (log) and mu_2 (mu2_cell_lb at cell a) ==")
print(f"{'w':>6} {'N':>5} {'lam1':>9} {'lam2':>10} {'al(log)':>8} {'al-lam2':>8} {'mu2cell':>9} {'lam1(H2)':>9} {'lam1H2-mu2':>11}")
for w, a, b, Na, al in spots:
    mu2 = float(mu2_cell_lb(a, b)[0])
    for N in (1500, 3000):
        l1, l2 = lam12(N, w)
        t = lam1_2x2(w)
        ok = (al >= l2) and (mu2 <= t)
        print(f"{w:>6.2f} {N:>5} {l1:>9.5f} {l2:>10.6f} {al:>8.3f} {al-l2:>8.4f} {mu2:>9.5f} {t:>9.5f} {t-mu2:>11.2e} {'OK' if ok else 'FATAL'}", flush=True)
print("\n== item 4: X5 nodes, mu (log) <= lam1 and alpha (log) >= lam2 at N=3000 ==")
nodes = [(30.0, 3.94648, 0.42465, 960), (40.0, 4.22923, 0.43859, 1280), (60.0, 4.62865, 0.45956, 1920)]
for w, mu, al, Nv in nodes:
    for N in (1500, 3000):
        l1, l2 = lam12(N, w)
        note = "mu>lam1 (allowed only if N<Nv=%d)" % Nv if mu > l1 else "ok"
        print(f"  w={w:>4.0f} N={N} lam1={l1:.6f} lam2={l2:.6f}  mu-lam1={mu-l1:+.2e} ({note})  alpha-lam2={al-l2:+.4f}", flush=True)
print(f"({time.time()-t0:.0f}s)")
