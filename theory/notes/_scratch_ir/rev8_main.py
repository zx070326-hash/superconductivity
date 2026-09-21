"""rev8: numerical cross-checks of the main argument (P4_uniqueness_strategy.md §3.A). Own construction."""
import numpy as np, time, sys
from scipy.linalg import expm, eigh
N = int(sys.argv[1]) if len(sys.argv) > 1 else 600
rng = np.random.default_rng(8)
t0 = time.time()
n = np.arange(N); u = (2.0 * n + 1.0) ** -0.5
A = np.abs(n[:, None] - n[None, :]); S = n[:, None] + n[None, :] + 1
UU = np.outer(u, u); offmask = (A > 0)

def O_of_f(f):
    """f array indexed k=0..2N (f[0] ignored)."""
    M = f[A] * UU * offmask + f[S] * UU
    cf = np.cumsum(f[:N + 1]) - f[0]
    M[n, n] -= 2.0 * u ** 2 * cf[:N]
    return M

def kernel_P(ws, ps, T):
    """F_T(k) = sum_i p_i [[k]](w_i / (2 pi T)),  [[k]](x) = x^2/(k^2+x^2), k = 0..2N."""
    k = np.arange(2 * N + 1, dtype=float)
    x = ws / (2 * np.pi * T)
    return sum(p * x_ ** 2 / (k ** 2 + x_ ** 2) for p, x_ in zip(ps, x))

def in_cone(v, tol=1e-12):
    v = v / np.sign(v[np.argmax(np.abs(v))])
    return bool(np.all(v >= -tol) and np.all(np.diff(v) <= tol)), v

def random_P():
    m = rng.integers(1, 12)
    ws = np.exp(rng.uniform(np.log(0.05), np.log(50.0), size=m)); ps = rng.dirichlet(np.ones(m))
    return ws, ps

# ---- (5a) monotonicity lambda_1(K_T) > lambda_1(K_T') for T < T', many random P
viol = 0; mind = np.inf; tests = 0
for trial in range(40):
    ws, ps = random_P()
    Ts = np.sort(np.exp(rng.uniform(np.log(0.003), np.log(30.0), size=5)))
    lam = [eigh(O_of_f(kernel_P(ws, ps, T)), eigvals_only=True)[-1] for T in Ts]
    for i in range(4):
        tests += 1; dlt = lam[i] - lam[i + 1]; mind = min(mind, dlt)
        if dlt <= 0: viol += 1
print(f"[5a] N={N}: {tests} pairs T<T' of random multi-atom P: violations of lam1(T)>lam1(T') = {viol}; min diff = {mind:.3e}")
# extreme cases: very small T (kernel ~ 1), very large T, single atom, two atoms far apart
for ws, ps, Ts in [ (np.array([1.0]), np.array([1.0]), [1e-4, 2e-4, 1e-3, 1.0, 100.0]),
                    (np.array([0.01, 100.0]), np.array([0.5, 0.5]), [0.001, 0.01, 0.1, 1, 10, 1000]),
                    (np.array([0.01, 100.0]), np.array([0.999, 0.001]), [0.001, 0.01, 0.1, 1, 10, 1000]) ]:
    lam = [eigh(O_of_f(kernel_P(ws, ps, T)), eigvals_only=True)[-1] for T in Ts]
    print(f"[5a] extreme P ws={ws} ps={ps}: lam1 over T={Ts}: {np.array(lam)}  strictly decreasing: {np.all(np.diff(lam) < 0)}")

# ---- (5b) top eigenvector in C for >= 50 random P (incl. extreme)
bad = 0; worst = 0.0; cnt = 0
for trial in range(50):
    ws, ps = random_P(); T = float(np.exp(rng.uniform(np.log(0.003), np.log(30.0))))
    M = O_of_f(kernel_P(ws, ps, T)); ev, V = eigh(M); v = V[:, -1]
    ok, v = in_cone(v, 1e-10); cnt += 1
    if not ok:
        bad += 1; worst = max(worst, max(-v.min(), np.diff(v).max()))
    # gap check for context
for ws, ps, T in [(np.array([1.0]), np.array([1.0]), 1e-4), (np.array([1.0]), np.array([1.0]), 100.0),
                  (np.array([0.01, 100.0]), np.array([0.999, 0.001]), 0.5), (np.array([0.01, 100.0]), np.array([0.001, 0.999]), 0.5)]:
    M = O_of_f(kernel_P(ws, ps, T)); ev, V = eigh(M); ok, v = in_cone(V[:, -1], 1e-10); cnt += 1
    if not ok: bad += 1
    print(f"[5b] extreme ws={ws} ps={ps} T={T}: top eigvec in C: {ok}; lam1={ev[-1]:.4e}, lam2={ev[-2]:.4e}, min v={v.min():.2e}, max diff(v)={np.diff(v).max():.2e}")
print(f"[5b] top eigenvector in C (nonneg, nonincreasing, tol 1e-10): {cnt - bad}/{cnt} pass; worst violation {worst:.2e}")

# ---- (5c) semigroup cone invariance e^{tK} v in C for random v in C, several t (finite N); plus a NON-monotone kernel control
def rand_cone_vec():
    a = rng.exponential(size=N) * (rng.random(N) < 0.3); a[rng.integers(0, N)] += 1.0
    v = np.cumsum(a[::-1])[::-1]             # v_n = sum_{L>n} a_L
    return v / np.linalg.norm(v)
fails = 0; tot = 0; worst = 0.0
for trial in range(6):
    ws, ps = random_P(); T = float(np.exp(rng.uniform(np.log(0.01), np.log(10.0))))
    M = O_of_f(kernel_P(ws, ps, T))
    for t in [0.1, 1.0, 5.0, 25.0]:
        E = expm(t * M)
        for _ in range(5):
            v = rand_cone_vec(); w = E @ v; w = w / np.linalg.norm(w); tot += 1
            ok, _ = in_cone(w, 1e-11)
            if not ok: fails += 1; worst = max(worst, max(-w.min(), np.diff(w).max()))
print(f"[5c] e^(tK) v in C for random v in C: {tot - fails}/{tot} pass (tol 1e-11), worst violation {worst:.2e}")
# control: non-monotone kernel f = delta_3 (Lemma B fails) -> does e^{tO} leave the cone?
f = np.zeros(2 * N + 1); f[3] = 1.0; E = expm(1.0 * O_of_f(f)); v = np.ones(N) / np.sqrt(N); v[:2] = 0; v = np.zeros(N); v[0] = 1.0
w = E @ v; print(f"[5c] control delta_3 kernel, e^(O) e_1: min = {w.min():.3e}, max diff = {np.diff(w).max():.3e}  (cone left iff either > 0)")
# Metzler check in step basis directly for the P-averaged kernel (sum of steps): coefficient matrix B
ws, ps = random_P(); T = 0.7; M = O_of_f(kernel_P(ws, ps, T))
Y = np.cumsum(M, axis=1)                     # Y[n, L'-1] = (M e_L')_n
B = np.zeros((N, N)); B[:-1, :] = Y[:-1, :] - Y[1:, :]; B[-1, :] = Y[-1, :]
offB = B - np.diag(np.diag(B))
print(f"[5c] step-basis coefficient matrix of K_T (random P): min off-diagonal = {offB.min():+.3e}  (Metzler iff >= 0); last row min = {B[-1, :-1].min():+.3e}")

# ---- (5d) <v, O[f] v> >= 0 for random v in C and random f >= 0 (non-monotone f allowed)
mn = np.inf
for trial in range(300):
    f = rng.exponential(size=2 * N + 1) * (rng.random(2 * N + 1) < 0.5) / (1 + np.arange(2 * N + 1)) ** rng.uniform(0, 3)
    f[0] = 0.0; M = O_of_f(f); v = rand_cone_vec(); mn = min(mn, v @ M @ v)
print(f"[5d] min <v,O[f]v> over 300 random (v in C, f>=0): {mn:+.3e}  (>= 0 claimed)")
# and a NON-cone v with f = delta_1 to show the sign is not trivial
f = np.zeros(2 * N + 1); f[1] = 1.0; M = O_of_f(f); v = np.zeros(N); v[5] = 1.0
print(f"[5d] control: <e_5(unit coord, not in C), O[delta_1] e_5> = {v @ M @ v:+.4f}")
print(f"total {time.time()-t0:.1f}s")
