"""rev5_cells.py -- certified per-cell bounds (from the script's own machinery) vs OWN-builder truth.
Any certified UB on H'' below the truth, F_lo > F, E_hi < E, g_lo > g, or s < true sin-angle is FATAL."""
import sys, time, numpy as np
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
from rev5_truth import Hb, Phi, Psi, fd_hpp, truth
import p4_extremal_concavity as pec
from p4_interchange_uniformity import alpha_certificate_H, mu2_cell_lb
from p4_interval_extension import schur_full_seg_w0
from mpmath import mp
t0 = time.time()
NT = 1500
BAD = []
def flag(cond, msg):
    if not cond:
        BAD.append(msg); print("   !!! VIOLATION:", msg)

def truth_at(w):
    lam1, lam2, g, F, E, v, e, V = truth(NT, w)
    return dict(lam1=lam1, lam2=lam2, g=g, F=F, E=E, v=v, hpp=fd_hpp(NT, w, 1e-3))

print("[A] Stage S first cell (0,0.3] and the point certificates for w<=0.3")
ub, info = pec.cell_S(0.0, 0.3)
print(f"   cell (0,0.3]: hpp_ub={ub:.6f} F_lo={info['F_lo']:.6f} E_hi={info['E_hi']:.6f} g~_lo={info['g_lo']:.6f} s={info['s']:.4e} pos={info['pos']:.6f}")
vS = np.array([float(t) for t in pec.SDATA['tr'].v])
for w in (0.01, 0.05, 0.1, 0.2, 0.3):
    T = truth_at(w)
    vv = np.zeros(NT); vv[:len(vS)] = vS
    sin_true = np.sqrt(max(0.0, 1 - float(vv @ T['v']) ** 2))
    gt = T['g'] / (w * w)
    print(f"   w={w:4.2f}: H''(fd)={T['hpp']:+.6f}  F={T['F']:.6f}  E={T['E']:.6f}  g~={gt:.6f}  sin={sin_true:.3e}  |  x E^2/g~ = {w*w*T['E']**2/gt:.6f}")
    flag(ub >= T['hpp'], f"S cell: ub {ub} < truth {T['hpp']} at w={w}")
    flag(info['F_lo'] <= T['F'], f"S cell F_lo at w={w}"); flag(info['E_hi'] >= T['E'], f"S cell E_hi at w={w}")
    flag(info['g_lo'] <= gt, f"S cell g_lo at w={w}"); flag(info['s'] >= sin_true, f"S cell s at w={w}")

print("[B] point certificates (pec._point_bound) vs truth (N=1500 fd), H'' itself (not 4w^2 H'')")
print(f"   {'w':>6} {'cert UB':>14} {'truth fd':>14} {'UB/truth':>9} | {'F_lo':>11} {'F':>11} | {'E_hi':>11} {'E':>11} | {'g_lo':>9} {'g':>9} | {'s':>9} {'sin':>9}")
for w in (0.05, 0.3, 1.0, 2.0, 5.0, 10.0, 20.0, 44.0):
    ub, info = pec._point_bound(w)
    T = truth_at(w)
    g_t = T['g'] / (w * w) if info['regime'] == 'S' else T['g']
    sin_true = float('nan')
    if info['regime'] == 'L':
        tr = pec.trial_at(info['p'], info['Ns']); vh = np.zeros(NT); vh[:tr.Nv] = [float(t) for t in tr.v]
        sin_true = np.sqrt(max(0.0, 1 - float(vh @ T['v']) ** 2))
    print(f"   {w:6.2f} {ub:14.6e} {T['hpp']:14.6e} {ub/T['hpp']:9.4f} | {info['F_lo']:11.5e} {T['F']:11.5e} | {info['E_hi']:11.5e} {T['E']:11.5e} | {info['g_lo']:9.4f} {g_t:9.4f} | {info['s']:9.3e} {sin_true:9.3e}")
    flag(ub >= T['hpp'], f"point w={w}: cert {ub} < truth {T['hpp']}")
    flag(info['F_lo'] <= T['F'], f"point F_lo w={w}"); flag(info['E_hi'] >= T['E'], f"point E_hi w={w}"); flag(info['g_lo'] <= g_t, f"point g_lo w={w}")
    if info['regime'] == 'L': flag(info['s'] >= sin_true, f"point s w={w}")

print("[C] REAL ladder cells: replicate run_L on the coarse cell containing w*, then check every fine cell at a, p, b")
coarse = pec.coarse_grid(0.3, 45.0)
print(f"   coarse grid: {len(coarse)} nodes, last = {coarse[-1]}")
for wstar in (0.35, 1.0, 10.0, 44.1):
    i = max(j for j in range(len(coarse) - 1) if coarse[j] <= wstar)
    A, B = coarse[i], coarse[i + 1]
    Na = min(pec.na_band(A), pec.na_band(B))
    _, okA, cA, _, _ = alpha_certificate_H(A, Na=Na); _, okB, cB, _, _ = alpha_certificate_H(B, Na=Na)
    assert okA and okB
    L1 = schur_full_seg_w0(A, B, 1)
    alpha = pec.alpha_cell(A, B, cA, cB, L1, Na)
    mu2 = mu2_cell_lb(A, B)[0]
    k = 1; got = None
    for _ in range(pec.MAX_SPLIT + 1):
        sub = np.linspace(A, B, k + 1); acc = []; bad = False
        for j in range(k):
            ub, info = pec.cell_L(float(sub[j]), float(sub[j + 1]), alpha, L1, mu2)
            if ub is None or ub >= 0 or info['ratio'] > pec.TARGET: bad = True; break
            acc.append(info)
        if not bad: got = acc; break
        k *= 2
    print(f"   coarse [{A:.5f},{B:.5f}] Na={Na} cA={float(cA):.4e} cB={float(cB):.4e} alpha={float(alpha):.4f} L1={L1:.4f} mu2={float(mu2):.4f} -> {k} fine cells")
    # true lam2 on the coarse cell endpoints vs alpha
    for w in (A, B):
        T = truth_at(w); flag(float(alpha) >= T['lam2'], f"alpha < lam2 at w={w}")
        print(f"      w={w:.5f}: lam2(true)={T['lam2']:.4e} <= alpha={float(alpha):.4f};  lam1={T['lam1']:.5f} >= mu2={float(mu2):.4f}")
    sel = got if len(got) <= 4 else [got[0], got[len(got)//2], got[-1]]
    for info in sel:
        a, b, p = info['a'], info['b'], info['p']
        tr = pec.trial_at(p, info['Ns']); vh = np.zeros(NT); vh[:tr.Nv] = [float(t) for t in tr.v]
        for w in (a, p, b):
            T = truth_at(w)
            sin_true = np.sqrt(max(0.0, 1 - float(vh @ T['v']) ** 2))
            print(f"      fine [{a:.5f},{b:.5f}] @w={w:.5f}: UB={info['hpp_ub']:.5e} vs H''={T['hpp']:.5e} ({info['hpp_ub']/T['hpp']:.3f}) | F_lo {info['F_lo']:.4e} <= F {T['F']:.4e} | E_hi {info['E_hi']:.4e} >= E {T['E']:.4e} | g_lo {info['g_lo']:.4f} <= g {T['g']:.4f} | s {info['s']:.3e} >= sin {sin_true:.3e} | s_p {info['s_p']:.3e}")
            flag(info['hpp_ub'] >= T['hpp'], f"fine cell [{a},{b}] w={w} UB<truth")
            flag(info['F_lo'] <= T['F'], f"fine F_lo [{a},{b}] w={w}"); flag(info['E_hi'] >= T['E'], f"fine E_hi [{a},{b}] w={w}")
            flag(info['g_lo'] <= T['g'], f"fine g_lo [{a},{b}] w={w}"); flag(info['s'] >= sin_true, f"fine s [{a},{b}] w={w}")
    print(f"   [{time.time()-t0:.0f}s]", flush=True)
print(f"\nVIOLATIONS: {len(BAD)}"); [print("  ", b) for b in BAD]
print(f"[done] {time.time()-t0:.1f}s")
