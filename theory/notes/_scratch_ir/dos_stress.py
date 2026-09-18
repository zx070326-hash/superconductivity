import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np

# From-scratch, no repo imports. All f in Bochner class P (Einstein-mode mixtures -> automatically
# positive-definite / in P). All DOS shape factors g(.) kept strictly > 0 by construction (no negative
# tau_n ever -- required by team-lead discipline).

def make_f(nmax, modes):
    f = np.zeros(nmax)
    kk = np.arange(nmax)
    for w, vp in modes:
        f += w * vp**2 / (vp**2 + kk**2)
    return f

def build_P(f, N, u):
    n_idx = np.arange(N)
    absdiff = np.abs(n_idx[:, None] - n_idx[None, :])
    L2 = f[absdiff] * np.outer(u, u)
    np.fill_diagonal(L2, 0.0)
    s = n_idx[:, None] + n_idx[None, :] + 1
    L3 = f[np.clip(s, 0, len(f)-1)] * np.outer(u, u)
    return L2 + L3

def build_L1_ansatzA(f, N, u2):
    # outer weight u2_n multiplies an UNWEIGHTED partial sum (same weight family used in D_u for P)
    diag = np.zeros(N)
    csum = 0.0
    for n in range(N):
        if n >= 1:
            csum += f[n]
        diag[n] = u2[n]*2.0*csum
    return diag

def build_L1_ansatzB(f, N, u2_kin, g):
    # kinematic outer weight (2n+1)^-1, DOS shape g_k placed INSIDE the sum on the source index k
    diag = np.zeros(N)
    csum = 0.0
    for n in range(N):
        if n >= 1:
            csum += g[n]*f[n]
        diag[n] = u2_kin[n]*2.0*csum
    return diag

def top2(M):
    w = np.linalg.eigvalsh(M)
    return w[-1], w[-2]

def gshape(kind, n, **kw):
    n = np.asarray(n, dtype=float)
    if kind == "flat":
        return np.ones_like(n)
    if kind == "logvhs":   # log van Hove peak AT E_F: enhances weight at LOW n (near E_F)
        Lam = kw.get("Lam", 50.0); kap = kw.get("kap", 0.5)
        # physical cutoff: log divergence only holds out to the bandwidth scale Lam; beyond that
        # DOS returns to a (different, but strictly positive) background value g_bg
        g_bg = kw.get("g_bg", 0.3)
        return np.maximum(1.0 + kap*np.log(Lam/(n+1.0)), g_bg)
    if kind == "lorentz_peak":  # narrow peak AT E_F: enhances low n
        A = kw.get("A", 3.0); w0 = kw.get("w0", 1.5)
        return 1.0 + A*w0**2/(w0**2+n**2)
    if kind == "lorentz_dip":  # pseudogap-like SUPPRESSION at E_F: depresses low n (stress test)
        A = kw.get("A", 0.9); w0 = kw.get("w0", 1.5)
        return 1.0 - A*w0**2/(w0**2+n**2)   # stays >0 since A<1
    if kind == "step":
        nc = kw.get("nc", 5); g0 = kw.get("g0", 0.3)
        return np.where(n < nc, 1.0, g0)
    if kind == "step_up":
        nc = kw.get("nc", 5); g1 = kw.get("g1", 3.0)
        return np.where(n < nc, g1, 1.0)
    raise ValueError(kind)

def eval_case(N, modes, gkind, ansatz, **kw):
    nmax = 2*N+2
    f = make_f(nmax, modes)
    n = np.arange(N)
    u2_kin = 1.0/(2*n+1)
    g = gshape(gkind, n, **kw)
    if ansatz == "A":
        u2 = g*u2_kin
        u = np.sqrt(u2)
        P = build_P(f, N, u)
        L1d = build_L1_ansatzA(f, N, u2)
    elif ansatz == "B":
        u = np.sqrt(u2_kin)
        P = build_P(f, N, u)
        L1d = build_L1_ansatzB(f, N, u2_kin, g)
    else:
        raise ValueError
    M = -np.diag(L1d) + P   # M0 = -L1[f]+P[f], single attractive mode (Eliashberg baseline)
    l1, l2 = top2(M)
    return l1, l2, (l1-l2)/l1

if __name__ == "__main__":
    modes = [(1.0, 1.0)]  # varpi=1 attractive Einstein mode, physical-ish ϖ

    print("=== Ansatz A (shared weight u_n^2=g*(2n+1)^-1 in BOTH L1 and P) ===")
    print("Theorem I' predicts: immunity preserved for ANY g>0 (rank-one leak at n=0 unchanged).")
    for gk, kw in [("flat", {}), ("logvhs", {}), ("lorentz_peak", {}), ("lorentz_dip", {}),
                   ("step", {}), ("step_up", {})]:
        rows = []
        for N in (300, 600):
            l1, l2, gr = eval_case(N, modes, gk, "A", **kw)
            rows.append((N, l1, l2, gr))
        (N1, l1a, l2a, gra), (N2, l1b, l2b, grb) = rows
        print(f"{gk:14s} N={N1}: l1={l1a:.5f} l2={l2a:.5f} gap/l1={gra:.5f}  |  "
              f"N={N2}: l1={l1b:.5f} l2={l2b:.5f} gap/l1={grb:.5f}")

    print()
    print("=== Ansatz B stress test (g INSIDE L1 sum only, u kinematic in P) ===")
    print("Testing whether low-n suppression (pseudogap-like g) can spread the leak beyond n=0.")
    for gk, kw in [("flat", {}), ("logvhs", {}), ("lorentz_peak", {}),
                   ("lorentz_dip", {"A": 0.9, "w0": 1.5}),
                   ("lorentz_dip", {"A": 0.99, "w0": 0.8}),
                   ("lorentz_dip", {"A": 0.999, "w0": 0.5}),
                   ("step", {"nc": 3, "g0": 0.1}),
                   ("step", {"nc": 1, "g0": 0.01})]:
        rows = []
        for N in (300, 600):
            l1, l2, gr = eval_case(N, modes, gk, "B", **kw)
            rows.append((N, l1, l2, gr))
        (N1, l1a, l2a, gra), (N2, l1b, l2b, grb) = rows
        print(f"{gk:14s} {str(kw):28s} N={N1}: l1={l1a:.5f} l2={l2a:.5f} gap/l1={gra:.5f}  |  "
              f"N={N2}: l1={l1b:.5f} l2={l2b:.5f} gap/l1={grb:.5f}")
