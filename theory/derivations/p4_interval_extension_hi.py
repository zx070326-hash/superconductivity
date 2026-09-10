r"""
P4 . INTERVAL EXTENSION, PART 2 (server run):  [0.2, 20]  -->  [0.2, W_HI]
      plus the mutilde = 0.33 node (mu* = 0.13 at lambda = 0.4: the last
      uncovered corner of the physics anchor X3 of `p4_interval_extension.py`).

RUN (server, 44 cores; the logged run `_scratch_ir/x4_run2.log` used exactly):
    OMP_NUM_THREADS=1 python -u p4_interval_extension_hi.py X4 X4L X5 X6 \
        --whi 60 --refine 2 --kna 16 --workers 40
(kna=16 is now also the CFG default -- reviewer `interval-hi-review` 2026-09-11
found the committed default was 14 while the log and the notes quote 16; the
script prints the full CFG at start so the provenance is in every log.)

Everything is imported unchanged from the certified modules; the two design
knobs, both documented in the X0-X3 file and its review
(`notes/2026-09-10_interval_extension_review.md`), are
    Na(b)  = head size of the alpha certificate, ~ K_NA * b  (tail envelopes
             ||B||+||D|| at fixed Na grow ~ wbar; Poincare floor mu_2 -> 1.66);
    REFINE = grid refinement factor over the base geometric ratio (the
             Lipschitz slack L*h/2 with h/wbar fixed grows ~ 0.007 wbar even with
             the weighted-row-0 L, so beyond ~20 the grid must get finer).
Each cell is certified INDEPENDENTLY (both endpoints' Cholesky at the cell's
own Na, envelopes at the same Na, L = weighted-row-0 segment Schur bound), so
the work is embarrassingly parallel; no quantity is shared between cells.

STAGES
  X4  cells on [20, W_HI]; certified PREFIX from 20 upward for every mutilde in
      MUTILDES (contiguity with the [.,20] results of X1); plus the [5,20] cells
      for mutilde = 0.33 alone (X1 did not include that node).
  X4L the base grid [0.3, 5] (246 points, Na = 70) for mutilde = 0.33 alone:
      certified SUFFIX ending at 5 -> wbar_min(0.33), by the Stage-E recipe of
      `p4d_coulomb_certify` (base L there; here the weighted-row-0 L, and the
      base-L value is printed alongside for comparison with the published
      0.5391 / 0.6427 / 1.1686).
  X5  nodes {30, 40, 60}: delta, g_alt, B2/B4, K_*, s_max, bound-(6) finite-N
      consistency -- same protocol as X2.
  X6  physics anchor (EVIDENCE): the lambda = 0.4 row of X3 re-tagged with the
      new intervals, and the covered-entry count 60 -> ?

Everything entering a bound is mpmath dps=30; float64 only as trial vectors /
evidence, exactly as in the base files.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import sys
import time
import argparse
from multiprocessing import Pool
import numpy as np
from mpmath import mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_interchange_uniformity as piu                               # noqa: E402
from p4_interchange_uniformity import (H_mat, ROUND_PAD, alpha_certificate_H,  # noqa: E402
                                       mu2_cell_lb, rayleigh_mu, cert_sup_on_J,
                                       kappa_c3, frak_k_of_P, _fixed_shapes,
                                       SIGMA0, ETA_MAX)
from p4d_coulomb_certify import mu2_pencil_cell, K_c, M_mat              # noqa: E402
from p4_interval_extension import (schur_full_seg_w0, alpha_cell_ub_Na,   # noqa: E402
                                   BASE_RATIO)

mp.dps = 30
PASS = [0]
_LAST = [0]


def chk(name, ok, detail=""):
    if not ok:
        raise AssertionError(f"CHECK FAILED: {name} {detail}")
    PASS[0] += 1


def _summary(tag, note=""):
    k = PASS[0] - _LAST[0]
    _LAST[0] = PASS[0]
    print(f"  [{tag}] {k} checks passed. {note}", flush=True)


# ---------------------------------------------------------------------------
W_LO = 20.0
CFG = dict(whi=60.0, refine=2, kna=16.0, nodes=(30.0, 40.0, 60.0))
MUTILDES = (0.0, 0.1, 0.13, 0.26, 0.33)
# certified lower ends quoted from p4d_coulomb_certify Stage E (0.1/0.13/0.26,
# computed there with the BASE L) and p4_interval_extension X0 (0.0).  The 0.33
# lower end is computed here in X4L with the WEIGHTED-row-0 L (1.5666; the base
# L would give 1.6418, also printed).  Both conventions are rigorous; the
# published 0.1/0.13/0.26 values are kept as published rather than re-tightened
# (X4L shows the weighted L would move 0.26 from 1.1686 to 1.1416).
WMIN_QUOTED = {0.0: 0.2, 0.1: 0.5391, 0.13: 0.6427, 0.26: 1.1686}


def head_size(b, kna=None):
    kna = CFG["kna"] if kna is None else kna
    if b <= 5.0 + 1e-12:
        return 70
    if b <= 10.0 + 1e-12:
        return 140
    if b <= 20.0 + 1e-12:
        return 280
    return int(np.ceil(kna * b / 20.0) * 20)


def grid_on(lo, hi, refine):
    r = BASE_RATIO ** (1.0 / refine)
    g = [lo]
    while g[-1] * r < hi - 1e-12:
        g.append(g[-1] * r)
    g.append(hi)
    return [float(x) for x in g]


def cell_job(args):
    """One cell, fully self-contained: returns everything X4 needs."""
    a, b, kna = args
    Na = head_size(b, kna)
    t0 = time.time()
    _, okA, cA, _, _ = alpha_certificate_H(a, Na=Na)
    _, okB, cB, _, _ = alpha_certificate_H(b, Na=Na)
    if not (okA and okB):
        return dict(a=a, b=b, Na=Na, ok=False, t=time.time() - t0)
    Lw = schur_full_seg_w0(a, b, 1)
    Lb = piu.schur_full_seg(a, b, 1)
    al = alpha_cell_ub_Na(a, b, cA, cB, Lw, Na)
    floors = {}
    for mt in MUTILDES:
        if mt == 0.0:
            m, k1, k2, k3 = mu2_cell_lb(a, b)
            case = 0
        else:
            m, k1, k2, k3, case = mu2_pencil_cell(a, b, mt)
        floors[mt] = (float(m), bool(k1 >= 0 and k2 >= 0 and k3 >= 0), case)
    # evidence (float64, never in a bound): L dominates ||H'||, alpha_cell >= lam2
    hp = max(float(np.max(np.abs(np.linalg.eigvalsh(H_mat(1200, wq, 1)))))
             for wq in (a, 0.5 * (a + b), b))
    E = np.linalg.eigvalsh(H_mat(1500, 0.5 * (a + b)))
    return dict(a=a, b=b, Na=Na, ok=True, cA=float(cA), cB=float(cB), Lw=Lw, Lb=Lb,
                slack=Lw * (b - a) / 2, al=float(al), floors=floors,
                Hp=hp, lam1=float(E[-1]), lam2=float(E[-2]), t=time.time() - t0)


RESULT = {}


def run_X4(pool):
    whi, refine = CFG["whi"], CFG["refine"]
    print(f"=== STAGE X4: cells on [20, {whi:g}], refine={refine} (ratio "
          f"{BASE_RATIO**(1/refine):.5f}), Na = ceil({CFG['kna']:g} b / 20)*20 ===", flush=True)
    t0 = time.time()
    grid = grid_on(W_LO, whi, refine)
    cells = [(grid[i], grid[i + 1], CFG["kna"]) for i in range(len(grid) - 1)]
    # the [5,20] cells for mutilde=0.33 (base ratio, base banding): X1 lacked that node
    g520 = grid_on(5.0, 20.0, 1)
    cells520 = [(g520[i], g520[i + 1], CFG["kna"]) for i in range(len(g520) - 1)]
    print(f"  {len(cells)} cells on [20,{whi:g}] + {len(cells520)} cells on [5,20] "
          f"(for mutilde=0.33); {pool._processes} workers", flush=True)
    # heavy cells first for load balance
    order = sorted(range(len(cells)), key=lambda i: -cells[i][1])
    res_hi = [None] * len(cells)
    done = 0
    for i, r in zip(order, pool.imap(cell_job, [cells[i] for i in order])):
        res_hi[i] = r
        done += 1
        if r["ok"]:
            gtxt = " ".join(f"{mt:g}:{r['floors'][mt][0]-r['al']:+.3f}{'' if r['floors'][mt][1] else '!'}" for mt in MUTILDES)
            print(f"    cell {done:3d}/{len(cells)} [{r['a']:.3f},{r['b']:.3f}] Na={r['Na']} {r['t']:.0f}s "
                  f"slack={r['slack']:.3f} al={r['al']:.3f} g: {gtxt}", flush=True)
        else:
            print(f"    cell {done:3d}/{len(cells)} [{r['a']:.3f},{r['b']:.3f}] Na={r['Na']} CHOLESKY FAILED", flush=True)
    res_520 = pool.map(cell_job, cells520)
    print(f"  alpha side done ({time.time()-t0:.0f}s); per-cell time "
          f"{min(r['t'] for r in res_hi):.0f}-{max(r['t'] for r in res_hi):.0f}s; "
          f"Na in [{min(r['Na'] for r in res_hi)}, {max(r['Na'] for r in res_hi)}]", flush=True)
    for r in res_hi + res_520:
        chk("Cholesky-deflation certificate SPD (both endpoints)", r["ok"],
            f"cell [{r['a']:.4f},{r['b']:.4f}] Na={r['Na']}")
        chk("weighted-row-0 L <= base L", r["Lw"] <= r["Lb"] + 1e-15)
        chk("evidence: L_w0 >= ||H'(w)|| (N=1200)", r["Lw"] >= r["Hp"],
            f"cell [{r['a']:.4f},{r['b']:.4f}] L={r['Lw']:.4f} ||H'||={r['Hp']:.4f}")
        chk("evidence: alpha_cell >= true lam_2 (N=1500, midpoint)", r["al"] >= r["lam2"] - 1e-9,
            f"cell [{r['a']:.4f},{r['b']:.4f}] al={r['al']:.4f} lam2={r['lam2']:.5f}")
        chk("evidence: mu_2 <= true lam_1 (N=1500, midpoint)", r["floors"][0.0][0] <= r["lam1"] + 1e-9)
    print(f"  slack L*h/2 in [{min(r['slack'] for r in res_hi):.4f}, {max(r['slack'] for r in res_hi):.4f}], "
          f"||B||+||D||+c ~ alpha_cell-slack in [{min(r['al']-r['slack'] for r in res_hi):.4f}, "
          f"{max(r['al']-r['slack'] for r in res_hi):.4f}]", flush=True)
    print(f"  {'mutilde':>8}  {'certified above 20':>20}  {'cells':>9}  {'inf g_alt':>10}  "
          f"{'worst cell':>22}  {'min mu_2':>9}  {'first drop reason':>18}", flush=True)
    for mt in MUTILDES:
        n_ok = 0
        for r in res_hi:
            g = r["floors"][mt][0] - r["al"]
            if g > 0 and r["floors"][mt][1]:
                n_ok += 1
            else:
                break
        chk("nonempty certified prefix above 20", n_ok > 0, f"mt={mt}")
        kept = res_hi[:n_ok]
        gs = [r["floors"][mt][0] - r["al"] for r in kept]
        wmax = kept[-1]["b"]
        for r, g in zip(kept, gs):
            chk("cell g_alt > 0 on the certified prefix", g > 0)
            chk("cell k-signs hold on the certified prefix", r["floors"][mt][1])
        reason = "-"
        if n_ok < len(res_hi):
            r = res_hi[n_ok]
            g = r["floors"][mt][0] - r["al"]
            reason = ("k-sign" if not r["floors"][mt][1] else f"g={g:+.4f}")
        RESULT[("hi", mt)] = (wmax, min(gs), n_ok, len(res_hi))
        w = kept[int(np.argmin(gs))]
        print(f"  {mt:>8.3f}  [20, {wmax:>8.4f}]{'':>6}  {n_ok:>4}/{len(res_hi):<4}  {min(gs):>10.6f}  "
              f"[{w['a']:.4f},{w['b']:.4f}]  {min(r['floors'][mt][0] for r in kept):>9.5f}  {reason:>18}", flush=True)
    # mutilde = 0.33 on [5,20]: certified SUFFIX ending at 20 (contiguity with above)
    mt = 0.33
    n_ok = 0
    for r in reversed(res_520):
        g = r["floors"][mt][0] - r["al"]
        if g > 0 and r["floors"][mt][1]:
            n_ok += 1
        else:
            break
    kept = res_520[len(res_520) - n_ok:]
    if n_ok:
        gs = [r["floors"][mt][0] - r["al"] for r in kept]
        RESULT[("520", mt)] = (kept[0]["a"], min(gs), n_ok, len(res_520))
        print(f"  mutilde=0.33 on [5,20]: certified suffix [{kept[0]['a']:.4f}, 20], "
              f"{n_ok}/{len(res_520)} cells, inf g_alt {min(gs):.6f}", flush=True)
    else:
        RESULT[("520", mt)] = (None, None, 0, len(res_520))
        print("  mutilde=0.33 on [5,20]: NO cell certified (honest)", flush=True)
    # evidence on the pencil floors of retained hi cells (all window cases)
    for r in res_hi[::10] + [res_hi[-1]]:
        for mt in MUTILDES[1:]:
            for wq in (r["a"], 0.5 * (r["a"] + r["b"]), r["b"]):
                for Nc in (0, 1, 2, 9, 300):
                    v = float(np.linalg.eigvalsh(M_mat(2, wq, 0, mt, Nc))[-1])
                    chk("evidence: pencil mu_2 bound <= true 2x2 top eigenvalue",
                        v >= r["floors"][mt][0] - 1e-13, f"mt={mt} w={wq:.4f} Nc={Nc}")
    # true gap does not collapse on the extension (evidence): min over sampled windows
    for r in res_hi[::20] + [res_hi[-1]]:
        w = 0.5 * (r["a"] + r["b"])
        A = H_mat(1200, w)
        worst = np.inf
        for Nc in (0, 1, 2, 3, int(3 * w), int(10 * w), 1200):
            E = np.linalg.eigvalsh(A - 0.33 * K_c(1200, Nc))
            worst = min(worst, E[-1] - E[-2])
        chk("evidence: true gap(M, mutilde=0.33) > 0 at sampled windows", worst > 0, f"w={w:.3f}")
    _summary("Stage X4", f"({time.time()-t0:.1f}s)")


def run_X4L(pool):
    print("\n=== STAGE X4L: base grid [0.3,5] for mutilde=0.33 -> wbar_min(0.33) ===", flush=True)
    t0 = time.time()
    grid = [float(w) for w in piu.wbar_grid()]
    cells = [(grid[i], grid[i + 1], CFG["kna"]) for i in range(len(grid) - 1)]
    res = pool.map(cell_job, cells)
    for r in res:
        chk("Cholesky-deflation certificate SPD (both endpoints)", r["ok"])
        chk("weighted-row-0 L <= base L", r["Lw"] <= r["Lb"] + 1e-15)
        chk("evidence: alpha_cell >= true lam_2 (N=1500, midpoint)", r["al"] >= r["lam2"] - 1e-9)
    chk("base grid uses the base head size", all(r["Na"] == 70 for r in res))
    for mt in (0.26, 0.33):
        for tag in ("weighted L", "base L"):
            gs = []
            for r in res:
                al = r["al"] if tag == "weighted L" else r["al"] + (r["Lb"] - r["Lw"]) * (r["b"] - r["a"]) / 2
                gs.append((r["floors"][mt][0] - al, r["floors"][mt][1]))
            start = None
            for i in range(len(gs)):
                if all(g > 0 and ks for g, ks in gs[i:]):
                    start = i
                    break
            chk("a nonempty certified suffix exists on [0.3,5]", start is not None, f"mt={mt} {tag}")
            wmin = res[start]["a"]
            ginf = min(g for g, _ in gs[start:])
            if mt == 0.26 and tag == "base L":
                chk("mutilde=0.26 with the base L reproduces the published wbar_min=1.1686",
                    abs(wmin - 1.1686) < 1e-3, f"{wmin:.4f}")
            if mt == 0.33 and tag == "weighted L":
                RESULT[("lo", 0.33)] = (wmin, ginf, len(res) - start)
            print(f"  mutilde={mt:<5} {tag:>10}: wbar_min={wmin:.4f}  cells {len(res)-start}/{len(res)}  "
                  f"inf g_alt={ginf:.6f}", flush=True)
    _summary("Stage X4L", f"({time.time()-t0:.1f}s)")


NODE = {}


def node_job(w):
    Na = head_size(w)
    Nv = 2 * Na
    t0 = time.time()
    mu = rayleigh_mu(w, Nv=Nv)
    alpha, ok, c, nB, nD = alpha_certificate_H(w, Na=Na)
    mu2 = mu2_cell_lb(w, w)[0]
    B2 = cert_sup_on_J(w, 2)
    B4 = cert_sup_on_J(w, 4)
    # evidence truth at N >= 2 Nv: the Rayleigh trial vector has support Nv, so
    # its quotient may EXCEED lam_1 of any truncation smaller than Nv (Cauchy
    # interlacing) while still being <= lam_1 of the true operator.
    Nt = max(2000, 2 * Nv)
    E = np.linalg.eigvalsh(H_mat(Nt, w))
    return dict(w=w, Na=Na, Nv=Nv, ok=ok, mu=float(mu), alpha=float(alpha), mu2=float(mu2),
                B2=B2, B4=B4, lam1=float(E[-1]), lam2=float(E[-2]), Nt=Nt, t=time.time() - t0)


def run_X5(pool):
    nodes = CFG["nodes"]
    print(f"\n=== STAGE X5: nodes {nodes} -- delta, g_alt, B2/B4, K_*, s_max, bound (6) ===", flush=True)
    t0 = time.time()
    res = pool.map(node_job, nodes)
    print(f"  {'wbar':>6}  {'Na':>5}  {'Nv':>5}  {'mu<=lam1':>9}  {'alpha>=lam2':>11}  {'delta':>8}  "
          f"{'mu_2':>8}  {'g_alt':>8}  {'B2^inf':>8}  {'B4^inf':>8}  {'K_*':>10}  {'s_max':>7}", flush=True)
    for r in res:
        w = r["w"]
        chk("node Cholesky SPD", r["ok"], f"w={w}")
        delta = r["mu"] - r["alpha"]
        galt = r["mu2"] - r["alpha"]
        chk("node delta > 0", delta > 0, f"w={w}")
        chk("node g_alt > 0", galt > 0, f"w={w}")
        chk("evidence: mu <= lam1 (N = max(2000, 2 Nv))", r["mu"] <= r["lam1"] + 1e-9,
            f"w={w} mu={r['mu']:.6f} lam1(N={r['Nt']})={r['lam1']:.6f}")
        chk("evidence: alpha >= lam2 (N = max(2000, 2 Nv))", r["alpha"] >= r["lam2"] - 1e-9, f"w={w}")
        B2, B4 = r["B2"], r["B4"]
        s_max = min(SIGMA0, 0.95 * np.sqrt(galt / (2.0 * w**2 * B2)))
        K = w**4 * B4 * ETA_MAX**2 / 24.0 + 0.5 * w**4 * B2**2 / galt
        chk("uniform-denominator condition at s_max", w**2 * B2 * s_max**2 < galt / 2.0, f"w={w}")
        NODE[w] = dict(galt=galt, K=K, s_max=s_max)
        print(f"  {w:>6.1f}  {r['Na']:>5}  {r['Nv']:>5}  {r['mu']:>9.5f}  {r['alpha']:>11.5f}  {delta:>8.5f}  "
              f"{r['mu2']:>8.5f}  {galt:>8.5f}  {B2:>8.4f}  {B4:>8.4f}  {K:>10.4e}  {s_max:>7.4f}", flush=True)
    print(f"\n  finite-N consistency of bound (6) (N=400 proxy, sigma <= s_max):", flush=True)
    worst_ratio = 0.0
    for w in nodes:
        kap, c3, hN = kappa_c3(400, w)
        K = NODE[w]["K"]
        sigmas = tuple(s for s in (0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25) if s <= NODE[w]["s_max"])
        for name, eta, p in _fixed_shapes():
            g3 = float(p @ eta**3)
            per = []
            for s in sigmas:
                exact = frak_k_of_P(400, w, s * eta, p)
                pred = hN * (1.0 + 0.5 * kap * s**2 + c3 * g3 * s**3)
                per.append(abs(exact - pred) / s**4)
            r_ = max(per) / K
            worst_ratio = max(worst_ratio, r_)
            chk("finite-N remainder <= K_* sigma^4 at N=400", r_ <= 1.0, f"w={w} {name}")
            chk("|rem|/sigma^4 does not diverge as sigma->0", per[0] / per[-1] < 30.0, f"w={w} {name}")
            chk("|rem|/sigma^4 settles as sigma->0", abs(per[0] - per[1]) / per[1] < 0.25, f"w={w} {name}")
            print(f"    {w:>6.1f}  {name:>9}  max|rem|/s^4={max(per):.4e}  K_*={K:.4e}  ratio={r_:.4f}", flush=True)
    chk("K_* not vacuous at the new nodes", worst_ratio > 1e-4)
    # N-independence with N up to 1200 (top eigenvector spreads ~ wbar)
    name, eta, p = [s for s in _fixed_shapes() if s[0] == "3pk-asym"][0]
    g3 = float(p @ eta**3)
    for w in nodes:
        vals = []
        for N in (400, 800, 1200):
            kap, c3, hN = kappa_c3(N, w)
            exact = frak_k_of_P(N, w, 0.01 * eta, p)
            vals.append(abs(exact - hN * (1.0 + 0.5 * kap * 1e-4 + c3 * g3 * 1e-6)) / 1e-8)
        rel = (max(vals) - min(vals)) / max(vals)
        chk("|rem|/sigma^4 N-independent for N in [400,1200]", rel < 1e-2, f"w={w} spread={rel:.2e}")
        print(f"    N-independence w={w:g}: " + "  ".join(f"{v:.6f}" for v in vals) + f"  spread={rel:.1e}", flush=True)
    _summary("Stage X5", f"({time.time()-t0:.1f}s)")


def interval_for(mt):
    key = min((m for m in MUTILDES if m >= mt - 1e-12), default=None)
    if key is None:
        return None
    hi = RESULT.get(("hi", key), (20.0,))[0]
    if key == 0.33:
        lo = RESULT.get(("lo", key), (None,))[0]
        if lo is None:
            lo = RESULT.get(("520", key), (None,))[0]
        if lo is None:
            lo = 20.0
        return lo, hi
    return WMIN_QUOTED[key], hi


def wbar_c_job(args):
    lam, mus, r = args
    mt = mus / lam
    ws = np.concatenate([np.arange(0.2, 5, 0.02), np.arange(5, 80.0, 0.1)])
    N = 800
    for w in ws:
        Nc = int(np.floor(r * w + 0.5))
        f = np.linalg.eigvalsh(H_mat(N, w) - mt * K_c(N, Nc))[-1]
        if f >= 1.0 / lam:
            return (lam, mus, r, float(w))
    return (lam, mus, r, np.nan)


def run_X6(pool):
    print("\n=== STAGE X6: physics anchor (EVIDENCE) -- coverage with the new intervals ===", flush=True)
    t0 = time.time()
    lams = (0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0)
    muss = (0.0, 0.05, 0.10, 0.13)
    jobs = [(lam, mus, r) for r in (10, 3) for lam in lams for mus in muss]
    res = {(l, m, r): w for (l, m, r, w) in pool.map(wbar_c_job, jobs)}
    cov_new, cov_20, cov_5 = 0, 0, 0
    for r in (10, 3):
        print(f"\n  r = w_c/w_E = {r}", flush=True)
        for lam in lams:
            cells = []
            for mus in muss:
                mt = mus / lam
                wc = res[(lam, mus, r)]
                iv = interval_for(mt)
                if np.isnan(wc) or iv is None:
                    tag = "--"
                else:
                    tag = "IN" if iv[0] <= wc <= iv[1] else "OUT"
                    cov_new += tag == "IN"
                    if mt <= 0.26 + 1e-12:
                        cov_20 += WMIN_QUOTED[min(m for m in MUTILDES if m >= mt - 1e-12)] <= wc <= 20.0
                        cov_5 += WMIN_QUOTED[min(m for m in MUTILDES if m >= mt - 1e-12)] <= wc <= 5.0
                cells.append(f"{wc:6.2f} {tag:>3} mt={mt:4.2f}")
            print(f"  {lam:>5.1f} | " + " | ".join(f"{c:>16}" for c in cells), flush=True)
    print(f"\n  ==> entries with wbar_c inside the certified interval: {cov_new} now, "
          f"{cov_20} with upper end 20 (X1), {cov_5} with upper end 5 (base); of 64", flush=True)
    chk("coverage does not shrink", cov_new >= cov_20)
    _summary("Stage X6", f"({time.time()-t0:.1f}s)")


def run_verdict():
    print("\n=== VERDICT ===", flush=True)
    for mt in MUTILDES:
        iv = interval_for(mt)
        wmax, ginf, n_ok, n = RESULT[("hi", mt)]
        conv = "lower end: X4L, weighted L" if mt == 0.33 else "lower end: published (base L)"
        print(f"  mutilde={mt:<5}: certified [{iv[0]:.4f}, {iv[1]:.4f}]; above 20: "
              f"{n_ok}/{n} cells, inf g_alt {ginf:.6f}  ({conv})", flush=True)
    for w in sorted(NODE):
        print(f"  node {w:g}: g_alt={NODE[w]['galt']:.5f} K_*={NODE[w]['K']:.4e} s_max={NODE[w]['s_max']:.4f}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stages", nargs="*", default=["X4", "X4L", "X5", "X6"])
    ap.add_argument("--whi", type=float, default=CFG["whi"])
    ap.add_argument("--refine", type=int, default=CFG["refine"])
    ap.add_argument("--kna", type=float, default=CFG["kna"])
    ap.add_argument("--workers", type=int, default=min(40, os.cpu_count() or 1))
    args = ap.parse_args()
    CFG.update(whi=args.whi, refine=args.refine, kna=args.kna)
    want = [s.upper() for s in args.stages]
    print(f"CFG: {CFG}  workers={args.workers}  stages={want}", flush=True)
    T0 = time.time()
    with Pool(args.workers) as pool:
        if "X4" in want:
            run_X4(pool)
        if "X4L" in want:
            run_X4L(pool)
        if "X5" in want:
            run_X5(pool)
        if "X6" in want:
            run_X6(pool)
    if "X4" in want:
        run_verdict()
    elif "X4L" in want:
        print(f"  wbar_min(0.33) = {RESULT[('lo', 0.33)][0]:.4f}")
    print(f"\nALL {PASS[0]} checks passed  ({time.time()-T0:.1f}s).", flush=True)
