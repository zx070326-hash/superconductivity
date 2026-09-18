"""rev3: server re-run of node_job(w) with kna=16, full-precision mu/alpha/mu2/B2/B4."""
import sys, time
sys.path.insert(0, "/root/p4cert")
import p4_interval_extension_hi as X
X.CFG["kna"] = 16.0
w = float(sys.argv[1]); t0 = time.time()
r = X.node_job(w)
print({k: (repr(v) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
print(f"done {time.time()-t0:.0f}s", flush=True)
