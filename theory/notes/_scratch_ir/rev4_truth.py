"""rev4 reviewer truth-side code.  Own builder from the KAY II entry formulas
(header of p4_interchange_uniformity.py); imports NOTHING from the derivations."""
import numpy as np, sys, time
np.set_printoptions(precision=10)

def build(N, f):
    """O[f] = -L1[f] + L2[f] + L3[f] on the N-head, f: int array k>=1 -> float."""
    n = np.arange(N)
    kk = np.arange(1, 2 * N + 2)
    ft = np.zeros(2 * N + 2); ft[1:] = f(kk.astype(float))
    M = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            v = ft[i + j + 1]
            if i != j:
                v += ft[abs(i - j)]
            M[i, j] = v / np.sqrt((2 * i + 1) * (2 * j + 1))
    # diagonal L1
    for i in range(N):
        M[i, i] -= 2.0 * ft[1:i + 1].sum() / (2 * i + 1)
    return M

def build_fast(N, f):
    n = np.arange(N)
    kk = np.arange(1, 2 * N + 2)
    ft = np.zeros(2 * N + 2); ft[1:] = f(kk.astype(float))
    sq = np.sqrt(2.0 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :]); S = n[:, None] + n[None, :] + 1
    M = ft[S] + np.where(D > 0, ft[D], 0.0)
    M = M / (sq[:, None] * sq[None, :])
    cum = np.concatenate(([0.0], np.cumsum(2.0 * ft[1:N])))
    M[n, n] -= cum / (2.0 * n + 1.0)
    return M

fG = lambda g: (lambda k: k ** (-g))
fH = lambda w: (lambda k: w * w / (k * k + w * w))
fK = lambda w: (lambda k: 1.0 / (k * k * (k * k + w * w)))
fHp = lambda w: (lambda k: 2 * k * k * w / (k * k + w * w) ** 2)   # d/dw [[k]]

# ---- 0. loop builder == fast builder
for N in (5, 40):
    d = np.max(np.abs(build(N, fH(0.7)) - build_fast(N, fH(0.7))))
    print(f"builder loop vs fast N={N}: {d:.2e}")

# ---- 2(a) identity (I): H(w)/w^2 = G(2) - w^2 K(w)
print("\n[2a] identity (I) at N=300")
G2 = build_fast(300, fG(2))
for w in (0.05, 0.3, 1.0, 2.4, 6.0, 30.0):
    lhs = build_fast(300, fH(w)) / w ** 2
    rhs = G2 - w * w * build_fast(300, fK(w))
    print(f"  w={w:5}: max|H/w^2 - (G2 - w^2 K)| = {np.max(np.abs(lhs-rhs)):.2e}")

# ---- 2(d) G(2) spectrum, residual of the 200-head trial vector
print("\n[2c/2d] G(2) top spectrum and residual of vhat (200-head eigvec)")
for N in (2000, 3000):
    t0 = time.time()
    G = build_fast(N, fG(2))
    ev, V = np.linalg.eigh(G)
    print(f"  N={N}: lam1={ev[-1]:.12f} lam2={ev[-2]:.8e} lam3={ev[-3]:.4e} lam_min={ev[0]:.8f}  ({time.time()-t0:.1f}s)")
    vh = np.linalg.eigh(G[:200, :200])[1][:, -1]
    vh = vh * np.sign(vh[0]); assert vh.min() > 0
    vfull = np.zeros(N); vfull[:200] = vh
    lamhat = vfull @ G @ vfull
    r = G @ vfull - lamhat * vfull
    print(f"       lamhat={lamhat:.15f}  ||r_full||(N={N})={np.linalg.norm(r):.6e}  ||r_head||={np.linalg.norm(r[:200]):.3e}")
    v0 = V[:, -1] * np.sign(V[0, -1])
    print(f"       sin angle(vhat,v0) = {np.sqrt(max(0,1-(vfull@v0)**2)):.3e};  min v0 = {v0.min():.2e}")
    if N == 2000:
        v0_2000 = v0; G2000 = G; ev2000 = ev

# ---- 2(f) kappa, Lam, ehat at N=1500 (true v0 of the 1500-head) and at N=2000
print("\n[2f] kappa, Lam, ehat from own float64 (v0 = top eigvec of G(2) head)")
for N in (1500, 2000):
    G = build_fast(N, fG(2)) if N != 2000 else G2000
    v0 = np.linalg.eigh(G)[1][:, -1]; v0 *= np.sign(v0[0])
    for w in (0.1, 0.5, 1.0, 1.7, 2.4, 4.0, 6.0):
        K = build_fast(N, fK(w))
        kap = v0 @ K @ v0
        Kv = K @ v0
        eh = np.sqrt(max(0.0, Kv @ Kv - kap * kap))
        evK = np.linalg.eigvalsh(K)
        Lam = -evK[0]
        # s=1 structural constraint  w^2 (Lam + 2 ehat)
        print(f"  N={N} w={w:4}: kappa={kap:.9f}  Lam={Lam:.9f}  ehat={eh:.9f}  ||K||={max(abs(evK[0]),abs(evK[-1])):.6f}  w^2(Lam+2eh)={w*w*(Lam+2*eh):.5f}")

# ---- 3. Segment III: ||H(w)|| vs C S1(w), and per-row exact sums vs analytic
print("\n[3] Segment III")
C = 1/3 + 2/np.sqrt(3) + 1/np.sqrt(15)
def S1(w): return (np.pi * w / np.tanh(np.pi * w) - 1) / 2
def phi(n): return 2/(2*n+1) + np.sqrt(2*n-1)/(n*np.sqrt(2*n+1)) + 1/np.sqrt((2*n+1)*(2*n+3)) + 1/np.sqrt(2*n+1)
print(f"  C = {C:.9f}; row0 const = {1+1/np.sqrt(3):.7f}; phi(2) = {phi(2):.7f}; wbar_1 = {C*np.pi/(2*1.31714):.7f}")
print(f"  (note: using [[1]]+[[2]] <= S1 in row 1 would give {1/3+1/np.sqrt(3)+1/3+1/np.sqrt(15):.5f}; then C would be the row-0 constant {1+1/np.sqrt(3):.5f})")
for w in (2.0, 2.08, 3.0, 5.0, 10.0, 30.0):
    N = 3000
    H = build_fast(N, fH(w))
    ev = np.linalg.eigvalsh(H)
    nrm = max(abs(ev[0]), abs(ev[-1]))
    rows = np.abs(H).sum(axis=1)
    b1 = w*w/(1+w*w); b2 = w*w/(4+w*w); s1 = S1(w)
    an0 = 2*s1/np.sqrt(3) + b1*(1-1/np.sqrt(3))
    an1 = b1*(1/3+1/np.sqrt(3)) + b2*(1/np.sqrt(3)-1/3) + s1*(1/3+1/np.sqrt(15))
    ok = rows[0] <= an0 and rows[1] <= an1 and all(rows[n] <= s1*phi(n) for n in range(2, 400))
    print(f"  w={w:5}: lam1={ev[-1]:.6f} ||H||={nrm:.6f}  C*S1={C*s1:.4f}  C*pi*w/2={C*np.pi*w/2:.4f}  1.31714 w^2={1.31714*w*w:.4f}  "
          f"S1<=pi w/2:{s1 <= np.pi*w/2}  rows: R0={rows[0]:.4f}<={an0:.4f} R1={rows[1]:.4f}<={an1:.4f} maxR_n>=2={rows[2:].max():.4f} rowbounds_ok={ok}")

# ---- 4. n=1 sharper constant sanity: is R_1 really the max row for the 'S1' form?
# ---- 5. Stage M spot: finite-difference h' at N=1500
print("\n[5] Stage M spot: central-difference h'(w) at N=1500, and Hellmann-Feynman <v,H'v>")
for w in (1.5, 5.0, 20.0, 29.9):
    N = 1500; d = 1e-3
    lp = np.linalg.eigvalsh(build_fast(N, fH(w + d)))[-1]
    lm = np.linalg.eigvalsh(build_fast(N, fH(w - d)))[-1]
    fd = (lp - lm) / (2 * d)
    H = build_fast(N, fH(w)); e, V = np.linalg.eigh(H); v = V[:, -1]
    hf = v @ build_fast(N, fHp(w)) @ v
    print(f"  w={w:5}: FD h'={fd:.6f}  HF <v,H'v>={hf:.6f}  lam1={e[-1]:.6f} gap={e[-1]-e[-2]:.4f}")
