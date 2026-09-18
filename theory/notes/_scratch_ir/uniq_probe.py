# Uniqueness strategy probe.  Danskin: d/dt k(t) = <v_t, sum_w  w H'(w t) P  v_t>.  If for every v in a cone C
# containing all Perron vectors the scalar map w -> <v,H(w)v> is nondecreasing, uniqueness follows for ALL P.
# Test: (a) is <v,H'(w)v> >= 0 for v = Perron vector of H(w0), all (w,w0)?  (b) for Perron vectors of 2-atom P?
#       (c) for random positive nonincreasing v?  (d) are Perron vectors positive & nonincreasing in n?
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "derivations"))
from p4_interchange_uniformity import H_mat
N = 800
ws = [0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0]
Hp = {w: H_mat(N, w, 1) for w in ws}
perron = {}
for w0 in ws:
    e, V = np.linalg.eigh(H_mat(N, w0)); v = V[:, -1]; v *= np.sign(v[0]); perron[w0] = v
print("(d) Perron vectors: min entry, #increases in n (v_{n+1}>v_n), #increases of v*sqrt(2n+1)")
n = np.arange(N)
for w0 in ws:
    v = perron[w0]; x = v*np.sqrt(2*n+1)
    print(f"   w0={w0:5.1f}: min v={v.min():+.2e}  inc(v)={int(np.sum(np.diff(v)>1e-14))}  inc(v*sqrt)={int(np.sum(np.diff(x)>1e-14))}  first inc idx={np.argmax(np.diff(x)>1e-14) if np.any(np.diff(x)>1e-14) else -1}")
print("(a) min over w of <v_{w0}, H'(w) v_{w0}> / (||H'(w)||):")
worst = 0
for w0 in ws:
    row = []
    for w in ws:
        q = perron[w0] @ Hp[w] @ perron[w0]; nrm = np.max(np.abs(np.linalg.eigvalsh(Hp[w])))
        row.append(q/nrm); worst = min(worst, q)
    print(f"   w0={w0:5.1f}: " + " ".join(f"{r:+.3f}" for r in row))
print("   worst raw value:", worst)
print("(b) 2-atom P Perron vectors, w-derivative of k(t):  <v, w1 H'(w1) p + w2 H'(w2)(1-p) v>")
rng = np.random.default_rng(1); bad = 0; vals = []
for _ in range(60):
    w1, w2 = np.exp(rng.uniform(np.log(0.1), np.log(40), 2)); p = rng.uniform(0.05, 0.95)
    K = p*H_mat(N, w1) + (1-p)*H_mat(N, w2); e, V = np.linalg.eigh(K); v = V[:, -1]
    M = p*w1*H_mat(N, w1, 1) + (1-p)*w2*H_mat(N, w2, 1); d = v@M@v; vals.append(d/e[-1])
    if d <= 0: bad += 1
print(f"   negatives: {bad}/60 ; min (t k'(t))/k = {min(vals):.4f}, max = {max(vals):.4f}")
print("(c) random positive nonincreasing v: min <v,H'(w)v>/(||v||^2 ||H'||)")
for w in (0.3, 1.0, 5.0, 20.0):
    Hw = Hp[w]; nrm = np.max(np.abs(np.linalg.eigvalsh(Hw))); mn = 1e9
    for _ in range(300):
        inc = rng.exponential(size=N) * np.exp(-rng.uniform(0, 3)*np.arange(N)/rng.uniform(5, N))
        v = np.cumsum(inc[::-1])[::-1]; v /= np.linalg.norm(v); mn = min(mn, v@Hw@v/nrm)
    # adversarial: v = e_0 + ... single-step sequences 1_{n<L}
    for L in (1, 2, 5, 20, 100, 400):
        v = np.zeros(N); v[:L] = 1; v /= np.linalg.norm(v); mn = min(mn, v@Hw@v/nrm)
    print(f"   w={w:5.1f}: min = {mn:+.4f}")
