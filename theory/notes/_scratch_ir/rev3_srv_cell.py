"""rev3: server spot re-run of two X4 cells through the committed cell_job (kna=16 as in x4_run2.log)."""
import sys, time, json
sys.path.insert(0, "/root/p4cert")
import p4_interval_extension_hi as X
X.CFG["kna"] = 16.0
w = float(sys.argv[1]); which = sys.argv[2]
if which == "hi":
    g = X.grid_on(20.0, 60.0, 2)
else:
    g = X.grid_on(5.0, 20.0, 1)
i = max(j for j in range(len(g) - 1) if g[j] <= w + 1e-12)
a, b = g[i], g[i + 1]
t0 = time.time()
r = X.cell_job((a, b, 16.0))
r["floors"] = {str(k): v for k, v in r["floors"].items()}
print(json.dumps(r), flush=True)
print(f"done {time.time()-t0:.0f}s", flush=True)
