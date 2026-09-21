"""rev5_truth.py -- adversarial review of p4_extremal_concavity.py.
OWN float64 builder written from the KAY II entry formulas in the header of
p4_interchange_uniformity.py (nothing imported for the truth side)."""
import sys, os, time
import numpy as np

def brk(k, w, order=0):
    k2 = k * k; d = k2 + w * w
    if order == 0: return w * w / d
    if order == 1: return 2 * k2 * w / d ** 2
    if order == 2: return 2 * k2 * (k2 - 3 * w * w) / d ** 3
    raise ValueError

def opK(N, f):
    """O[f] = -L1[f] + L2[f] + L3[f] on the N x N head; f indexed 0..>=2N (f[0] unused)."""
    n = np.arange(N); sq = np.sqrt(2.0 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :]); S = n[:, None] + n[None, :] + 1
    L2 = np.where(D > 0, f[D], 0.0) / (sq[:, None] * sq[None, :])
    L3 = f[S] / (sq[:, None] * sq[None, :])
    cs = np.concatenate(([0.0], np.cumsum(2.0 * f[1:N])))      # cs[n] = sum_{k=1}^n 2 f(k)
    L1 = np.diag(cs / (2.0 * n + 1.0))
    return -L1 + L2 + L3

def ktab(N, fun):
    k = np.arange(0, 2 * N + 2, dtype=float); f = np.zeros_like(k); f[1:] = fun(k[1:]); return f

def Hb(N, w, order=0):
    return opK(N, ktab(N, lambda k: brk(k, w, order)))

def Phi(N, x): return opK(N, ktab(N, lambda k: k * k / (k * k + x) ** 3))
def Psi(N, x): return opK(N, ktab(N, lambda k: 1.0 / (k * k + x) ** 2))

def top(A):
    e, V = np.linalg.eigh(A); return e, V[:, -1] * np.sign(V[0, -1])

def truth(N, w):
    """(lam1, lam2, g, F, E, v, e, V) for the truncated H_N(w)."""
    x = w * w
    e, V = np.linalg.eigh(Hb(N, w)); v = V[:, -1] * np.sign(V[0, -1])
    F = float(v @ Phi(N, x) @ v)
    y = Psi(N, x) @ v; y = y - float(v @ y) * v; E = float(np.linalg.norm(y))
    return e[-1], e[-2], e[-1] - e[-2], F, E, v, e, V

def fd_hpp(N, w, rel=1e-3):
    x = w * w; dx = rel * x
    f = lambda xx: np.linalg.eigvalsh(Hb(N, np.sqrt(xx)))[-1]
    return (f(x + dx) - 2 * f(x) + f(x - dx)) / dx ** 2

if __name__ == "__main__":
    t0 = time.time()
    sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
    import p4_interchange_uniformity as piu   # comparison ONLY

    print("[1] own builder vs piu.H_matrix / H_mat (verbatim T5)")
    for N in (30,):
        for w in (0.05, 0.3, 2.4, 20.0):
            for od in (0, 1, 2):
                d1 = np.max(np.abs(Hb(N, w, od) - piu.H_matrix(N, w, od)))
                print(f"   N={N} w={w} order={od}: max|diff| = {d1:.2e}")

    print("[2] identities (I1),(I2) elementwise, own builder, N=300")
    N = 300
    for w in (0.05, 0.3, 1.0, 7.0, 20.0, 44.0):
        x = w * w
        Axx = (Hb(N, w, 2) - Hb(N, w, 1) / w) / (4 * x)
        P = Phi(N, x)
        d1 = np.max(np.abs(Axx + 2 * P)) / max(1.0, np.max(np.abs(P)))
        Ax = Hb(N, w, 1) / (2 * w)
        Q = Psi(N, x)
        rhs = Hb(N, w) / x - x * Q
        d2 = np.max(np.abs(Ax - rhs)) / max(1.0, np.max(np.abs(rhs)))
        # diagonal -L1 term separately
        dd1 = np.max(np.abs(np.diag(Axx) + 2 * np.diag(P)))
        dd2 = np.max(np.abs(np.diag(Ax) - np.diag(rhs)))
        # scalar kernel identity
        k = np.arange(1, 2000, dtype=float)
        sc = np.max(np.abs(brk(k, w, 2) - brk(k, w, 1) / w + 8 * x * k * k / (k * k + x) ** 3))
        print(f"   w={w:5.2f}: (I1) rel {d1:.2e} (diag abs {dd1:.2e}) | (I2) rel {d2:.2e} (diag abs {dd2:.2e}) | kernel id {sc:.2e}")

    print("[2'] (I3) with the truncated eigenvector, N=1500: ||(I-P)A_x v + x (I-P) Psi v||")
    for w in (0.3, 1.0, 10.0, 44.0):
        N = 1500; x = w * w
        lam1, lam2, g, F, E, v, e, V = truth(N, w)
        Ax = Hb(N, w, 1) / (2 * w)
        a = Ax @ v; a -= float(v @ a) * v
        b = Psi(N, x) @ v; b -= float(v @ b) * v
        print(f"   w={w}: ||resid|| = {np.linalg.norm(a + x * b):.2e}  (||(I-P)A_x v|| = {np.linalg.norm(a):.3e} = x*E, E={E:.3e})")

    print("[3] structure coefficients c_k = T_k + U_k - P_k from OWN code; sign; identity with the quadratic form")
    def support_size(w): return max(140, int(10 * float(w)) // 2 * 2)
    for p in (0.4, 1.0, 3.0, 10.0, 30.0, 44.1):
        Ns = support_size(p)
        e, vh = top(Hb(Ns, p)); vh /= np.linalg.norm(vh)
        u = vh / np.sqrt(2 * np.arange(Ns) + 1)
        K = 2 * Ns
        T = np.zeros(K); U = np.zeros(K); Pk = np.zeros(K)
        for k in range(1, Ns): T[k] = 2 * np.dot(u[:Ns - k], u[k:])
        for s in range(1, 2 * Ns):
            lo = max(0, s - Ns); hi = min(s - 1, Ns - 1)
            if hi >= lo: U[s] = np.dot(u[lo:hi + 1], u[s - 1 - hi:s - lo][::-1])
        cum = np.cumsum((u * u)[::-1])[::-1]     # cum[n] = sum_{m>=n} u_m^2
        for k in range(1, Ns): Pk[k] = 2 * cum[k]
        c = T + U - Pk; c[0] = 0
        # identity with <vh, O[f] vh> for f = phi, psi at wbar = 1.1p and the H kernel
        ok = []
        for fun in (lambda k: k * k / (k * k + (1.1 * p) ** 2) ** 3, lambda k: 1 / (k * k + (1.1 * p) ** 2) ** 2, lambda k: brk(k, p, 0)):
            f = ktab(Ns, fun)
            ok.append(abs(np.dot(c, f[:K]) - vh @ opK(Ns, f) @ vh))
        print(f"   p={p:5.1f} Ns={Ns}: min c_k = {c.min():+.3e} at k={c.argmin()}, max c_k = {c.max():.3e}, #neg = {(c < 0).sum()}, "
              f"min c_k (k<Ns) = {c[1:Ns].min():+.3e}, qform identity err = {max(ok):.1e}, min vhat = {vh.min():.2e}")

    print("[4] (P2)+(I1)+(I3) vs finite-difference H''(x), own builder, N=1500; F, E, g, lam2")
    print(f"   {'w':>6} {'N':>5} {'lam1':>12} {'lam2':>11} {'g':>10} {'F':>11} {'E':>11} {'formula':>13} {'fd(1e-3)':>13} {'fd(3e-3)':>13} {'ratio pos/F':>11}")
    TR = {}
    for w in (0.05, 0.3, 1.0, 2.0, 5.0, 10.0, 20.0, 44.0):
        for N in (1500, 2500) if w in (0.05, 44.0, 1.0) else (1500,):
            x = w * w
            lam1, lam2, g, F, E, v, e, V = truth(N, w)
            cj = V.T @ (Psi(N, x) @ v)
            pos = x * x * float(np.sum(cj[:-1] ** 2 / (e[-1] - e[:-1])))
            formula = -2 * F + 2 * pos
            fd1 = fd_hpp(N, w, 1e-3); fd2 = fd_hpp(N, w, 3e-3)
            TR[(w, N)] = dict(lam1=lam1, lam2=lam2, g=g, F=F, E=E, hpp=formula, pos_true=pos, pos_g=x * x * E * E / g)
            print(f"   {w:6.2f} {N:5d} {lam1:12.6e} {lam2:11.4e} {g:10.4e} {F:11.5e} {E:11.5e} {formula:13.6e} {fd1:13.6e} {fd2:13.6e} {pos/F:11.4f}")
    np.save(r"D:\physics\superconductivity\theory\notes\_scratch_ir\rev5_truth_TR.npy", TR, allow_pickle=True)
    print(f"[done] {time.time()-t0:.1f}s")
