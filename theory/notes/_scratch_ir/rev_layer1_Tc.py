"""REVIEWER (audit item i + key-number reproduction):
 R1. independent Rayleigh (own u, own mpmath quotient, Nv=140, Nc=0..140) -> delta(1) at mt=0.1/0.13
     vs the certified 0.360987 / 0.282156; own 2x2 + nodal alpha -> g_alt(1) vs 0.349124 / 0.271540.
 R2. collapse point: wbar=0.3, mt=0.1: lam1(M), gap at N=400 and 800, scan Nc=0..30 for the worst.
 L1. N_c(T) is a step function of T: is f(T)=lam1(M(wbar(T), mt, N_c(T))) monotone in T?  Both
     conventions floor(x) and floor(x+1/2) (|omega_n|<omega_c), r=omega_c/omega_E in {3,10}.
     Count crossings of 1/lambda -> is Tc unique?  Then Tc(mu*) monotone (Layer 1)?
 L2. symmetric 2|u><u| vs Marsiglio one-sided 2|v><u| kernel: mu* needed for Tc at wbar=1."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from mpmath import mp
mp.dps = 30
import rev_common as rc
from p4_interchange_uniformity import alpha_certificate_H
t0 = time.time()

def my_rayleigh(w, mt, Nc, Nv, A, Am):
    Mf = A - mt * rc.Kc(Nv, Nc)
    v = np.linalg.eigh(Mf)[1][:, -1]
    vm = [mp.mpf(float(x)) for x in v]
    num = mp.fsum(vm[i] * mp.fsum(Am[i][j] * vm[j] for j in range(Nv)) for i in range(Nv))
    k = min(int(Nc), Nv)
    uv = mp.fsum(vm[n] / mp.sqrt(2 * n + 1) for n in range(k))
    num -= 2 * mp.mpf(mt) * uv * uv
    return num / mp.fsum(x * x for x in vm) - mp.mpf("1e-12")

print("R1. independent reproduction at wbar=1")
w = 1.0; Nv = 140
A = rc.H_mat(Nv, w); Am = [[mp.mpf(float(A[i, j])) for j in range(Nv)] for i in range(Nv)]
alpha = alpha_certificate_H(w)[0]
for mt, ref_d, ref_g in ((0.0, 0.657276, 0.610122), (0.1, 0.360987, 0.349124), (0.13, 0.282156, 0.271540)):
    best = min((my_rayleigh(w, mt, Nc, Nv, A, Am), Nc) for Nc in range(0, Nv + 1))
    delta = float(best[0] - alpha)
    E = np.linalg.eigvalsh(rc.H_mat(400, w) - mt * rc.Kc(400, best[1]))
    mu2 = min(rc.lam1_2x2(w, mt, c) for c in (0, 1, 2))
    galt = mu2 - float(alpha)
    print(f"  mt={mt}: delta={delta:.6f} (cert {ref_d}) worstNc={best[1]} | true gap@400 at that Nc={E[-1]-E[-2]:.6f} lam1={E[-1]:.6f} | "
          f"g_alt={galt:.6f} (cert {ref_g}); alpha={float(alpha):.6f}")

print("\nR2. collapse point wbar=0.3")
for mt in (0.1, 0.13):
    for N in (400, 800):
        H = rc.H_mat(N, 0.3)
        sc = [(np.linalg.eigvalsh(H - mt * rc.Kc(N, Nc))[-2:], Nc) for Nc in range(0, 31)]
        (E, Nc) = min(sc, key=lambda t: t[0][1] - t[0][0])
        print(f"  mt={mt} N={N}: worst Nc={Nc} lam1={E[1]:.4e} lam2={E[0]:.4e} gap={E[1]-E[0]:.4e}; "
              f"lam1 by Nc: " + " ".join(f"{n}:{e[1]:.4f}" for e, n in sc[:8]))
    H = rc.H_mat(400, 0.3)
    for r in (3, 10, 30):
        Nc = int(np.floor(r * 0.3 + 0.5))
        E = np.linalg.eigvalsh(H - mt * rc.Kc(400, Nc))
        print(f"     physical window r=w_c/w_E={r}: N_c={Nc} lam1={E[-1]:.4e} gap={E[-1]-E[-2]:.4e}")

print("\nL1. f(wbar) = lam1(H(wbar) - mt K_c(N_c(wbar))), N_c = window(r*wbar): monotone in wbar?")
N = 200
ws = np.linspace(0.3, 5.0, 471)
for conv, name in ((lambda x: int(np.floor(x)), "floor(x)"), (lambda x: int(np.floor(x + 0.5)), "floor(x+1/2)")):
    for r in (3, 10):
        for mt in (0.1, 0.13):
            f = np.array([np.linalg.eigvalsh(rc.H_mat(N, w) - mt * rc.Kc(N, conv(r * w)))[-1] for w in ws])
            d = np.diff(f); ndown = int(np.sum(d < -1e-12))
            jumps = [(ws[i], d[i]) for i in range(len(d)) if d[i] < -1e-12]
            cross = {}
            for lam in (0.5, 0.7, 1.0, 2.0):
                s = np.sign(f - 1.0 / lam); cross[lam] = int(np.sum(s[1:] != s[:-1]))
            big = min(jumps, key=lambda j: j[1]) if jumps else (float('nan'), 0.0)
            print(f"  {name:>12} r={r:>2} mt={mt}: #decreasing steps={ndown:3d} (largest drop {big[1]:+.2e} at wbar={big[0]:.3f}); crossings of 1/lam: {cross}")
print("  Layer 1 with Tc := sup{T: f>=1/lambda} (= smallest wbar with f>=1/lam), r=10, floor(x+1/2):")
for lam in (0.7, 1.0, 2.0):
    out = []
    for mus in (0.0, 0.05, 0.1, 0.13, 0.2):
        mt = mus / lam
        f = np.array([np.linalg.eigvalsh(rc.H_mat(N, w) - mt * rc.Kc(N, int(np.floor(10 * w + 0.5))))[-1] for w in ws])
        idx = np.where(f >= 1.0 / lam)[0]
        out.append((mus, ws[idx[0]] if len(idx) else np.nan))
    Tc = [1 / (2 * np.pi * o[1]) for o in out]
    mono = all(Tc[i] >= Tc[i + 1] - 1e-12 for i in range(len(Tc) - 1))
    print(f"    lam={lam}: wbar_c(mu*)={[f'{o[1]:.3f}' for o in out]} -> Tc/w_E={[f'{t:.4f}' for t in Tc]} nonincreasing={mono}")

print("\nL2. symmetric vs one-sided kernel: mu* that puts the top eigenvalue at 1 (wbar=1; lambda calibrated so mu*_sym=0.1)")
N = 400; w = 1.0
H = rc.H_mat(N, w)
lam_min = 1 / np.linalg.eigvalsh(H)[-1]
for r in (5, 10, 20, 50):
    Nc = int(np.floor(r * w + 0.5))
    u = np.zeros(N); u[:Nc] = 1 / np.sqrt(2 * np.arange(Nc) + 1.0)
    vfull = 1 / np.sqrt(2 * np.arange(N) + 1.0)
    def mus_sym(lam): return 1 / (2 * u @ np.linalg.solve(lam * H - np.eye(N), u))
    def mus_one(lam): return 1 / (2 * u @ np.linalg.solve(lam * H - np.eye(N), vfull))
    # mu*_sym(lambda) increases with lambda above lam_min; bisect for 0.1
    lo, hi = lam_min * (1 + 1e-9), 5.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mus_sym(mid) < 0.1: lo = mid
        else: hi = mid
        if hi - lo < 1e-13: break
    lam = hi
    print(f"  r={r:>2} N_c={Nc:>2}: lambda={lam:.6f}  mu*_sym={mus_sym(lam):.6f}  mu*_onesided={mus_one(lam):.6f}  rel diff={(mus_one(lam)-mus_sym(lam))/mus_sym(lam):+.3e}")
print(f"\ndone ({time.time()-t0:.0f}s)")
