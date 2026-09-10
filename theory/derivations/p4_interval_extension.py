r"""
P4 . INTERVAL EXTENSION of the interchange-uniformity theorem:
      [0.3, 5]  -->  [0.2, 20]   (mu* = 0)
      [wbar_min(mutilde), 5]  -->  [wbar_min(mutilde), 20]   (mutilde = 0.1, 0.13)
      [1.1686, 5]  -->  [1.1686, wbar_max(0.26)]              (mutilde = 0.26, lambda<1 corner)

RUN:  python theory/derivations/p4_interval_extension.py          (all stages)
      python theory/derivations/p4_interval_extension.py X1 X2    (a subset)

WHY (closure map 2026-09-10, `notes/2026-09-10_P4_closure_map.md`):
  Sec.5-(1): T5 Part 1 certifies its finite-N gap on [0.2,5] but block 7 lifts it
             to the true operator only on [0.3,5]  -> bookkeeping gap [0.2,0.3).
  Sec.3-(a): the lambda<1 hydride corner has its Tc solution at wbar_c > 5
             (lambda=0.6, mu*=0.13, w_c=10 w_E: 7.1; lambda=0.5: 13.5), so the
             rigorous Tc bound of P4-D says nothing there unless the interval
             theorem is pushed above 5.

WHAT THIS FILE DOES -- NOTHING NEW CONCEPTUALLY.  Every ingredient is imported
UNCHANGED from the two certified modules:
  * `p4_interchange_uniformity` (1668 checks): Cholesky-deflation head bound c,
    Schur segment envelopes ||B||,||D||, closed-form segment Lipschitz constant
    L = schur_full_seg, the exact 2x2 mu_2 with the k-sign monotonicity test,
    Rayleigh mu, cert_sup_on_J for B2/B4;
  * `p4d_coulomb_certify` (6534 checks): mu2_pencil_cell (the 2x2 bound for the
    pencil M = H - mutilde K_c, min over the three window cases, k-signs over
    all three cases).
The ONLY new degree of freedom is the HEAD SIZE Na of the alpha certificate.
The tail envelopes ||B||,||D|| at fixed Na grow ~ wbar (the diagonal of -H1 is
~ wbar for n << wbar), while the Poincare/Cauchy floor mu_2 = lam_1(H^(2))
saturates at ~1.66 as wbar -> inf.  With the base Na = 70 the closure
g_alt = mu_2 - alpha_cell therefore FAILS beyond wbar ~ 11.7 (probe
`_scratch_ir/probe_ends.py`: 47 of 121 cells on [5,20] negative), and the
failure is a CERTIFICATE artefact, not an operator property (true gap at
wbar=20 is 3.55, true lam_2 ~ 0.0024).  Scaling Na ~ 14 wbar restores it:
    Na = 70  for b <= 5   (the base grid, untouched)
    Na = 140 for b <= 10
    Na = 280 for b <= 20
The cost is the mpmath Cholesky at dps=30 (~ Na^3): 0.1 s at 70, 4 s at 280,
14 s at 400.  The extension stops at 20 for cost; nothing collapses there (true
gap 3.55 at wbar=20, nodal g_alt(30) = 1.17 at Na=420), but see the L note
below for what would ALSO have to change beyond ~20.

LIPSCHITZ CONSTANT (reviewer `interval-ext-review`, 2026-09-11, load-bearing
DOWNGRADE, adopted here): the base `schur_full_seg` bounds row 0 of the
envelope matrix by `row_env0` = 2*S1, i.e. it DROPS the (2m+1)^{-1/2} column
weight on row 0 -- rigorous, but 1.6x (wbar~0.2) to 5x (wbar~20) looser than
summing row 0 explicitly with the weight, which is what the n>=1 loop of the
same function already does.  Row 0 is the maximising row on every cell, so
this slack IS the Lipschitz constant.  `schur_full_seg_w0` below is the base
function with the loop started at n=0 (the n>=1 derivation of each term holds
verbatim at n=0: H1 diagonal is empty, H2/H3 rows are the same psi(|n-m|),
psi(n+m+1) sums, the integral tail bound only needs psi decreasing).  It is
asserted <= the base L on every cell.  Consequences (reviewer, confirmed by
the X1 rerun): the mutilde=0.26 prefix, which ended at 16.3172 under the base
L because 61% of alpha_cell there was L*h/2, now closes on the whole [5,20];
and at wbar=30 the geometric cell closes with margin instead of by 0.07.  The
base module is NOT modified (its 1668 checks / Gate A stay untouched).

Everything that enters a bound is mpmath (dps=30) on the same discipline as the
base files; float64 eigvalsh appears ONLY as a trial-vector / evidence device.
The base grids and check counts of the two imported modules are NOT touched
(`p4d_coulomb_certify` Gate A hard-codes the 1668 baseline), which is why this
is a separate file rather than a Stage 3C inside the base module.

STAGES
  X0  lower extension [0.2, 0.3], 36 geometric cells at the base ratio:
        mu*=0: g_alt > 0 and k-signs on every cell  (closes closure-map Sec.5-(1));
        mutilde > 0: reported honestly -- every cell FAILS (k-signs), i.e. the
        collapsed prefix of Stage E stays collapsed, wbar_min(mutilde) unchanged.
  X1  upper extension [5, 20], 121 geometric cells at the base ratio, banded Na:
        mutilde in {0, 0.1, 0.13}: all cells certified, inf g_alt reported;
        mutilde = 0.26: the certified PREFIX [5, wbar_max] is located (it ends
        near 16.3, where the Coulomb-shifted 2x2 floor drops below alpha_cell);
        float64 soundness evidence on sampled cells (alpha_cell >= true lam_2 at
        N=1500; mu_2 bound <= true 2x2 top eigenvalue for every window case).
  X2  new NODES wbar in {7, 10, 15, 20}: delta = mu - alpha (Rayleigh at Nv = 2Na),
        g_alt, B2^inf / B4^inf (cert_sup_on_J), K_*(wbar), s_max(wbar), and the
        N=400/800 finite-N consistency + sigma^4-rate check of bound (6) exactly as
        Stage 5 of the base file -- so theorem (6) is USABLE at the new nodes,
        not merely "the gap does not collapse".
  X3  physics anchor (EVIDENCE, no certificate): wbar_c(lambda, mu*) from
        lam_1(H - mutilde K_c(N_c)) = 1/lambda, N_c = floor(r wbar + 1/2),
        r = w_c/w_E in {10, 3}, N = 600; each entry is tagged with whether it
        lies inside the interval certified for its mutilde = mu*/lambda.
"""
import sys
import os
import time
import numpy as np
from mpmath import mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_interchange_uniformity as piu                               # noqa: E402
from p4_interchange_uniformity import (H_mat, ROUND_PAD, alpha_certificate_H,  # noqa: E402
                                       schur_full_seg, schur_coupling,
                                       schur_tailblock, spec_seg, mu2_cell_lb,
                                       rayleigh_mu, cert_sup_on_J, kappa_c3,
                                       frak_k_of_P, _fixed_shapes, SIGMA0,
                                       ETA_MAX, WBARS)
from p4d_coulomb_certify import mu2_pencil_cell, K_c, M_mat              # noqa: E402

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
    print(f"  [{tag}] {k} checks passed. {note}")


def schur_full_seg_w0(a, b, order, Mrows=200, Mcols=8000):
    """`piu.schur_full_seg` with row 0 summed EXPLICITLY with its column weights
    (loop from n=0) instead of the unweighted `row_env0` = 2*S1.  Same envelope
    matrix, same Schur test, same integral tails -- strictly a subset of the
    base function's own derivation, hence rigorous; never larger than it."""
    spec = spec_seg(a, b, order)
    C, A, p = float(spec[0]), float(spec[1]), spec[2]
    _, int_tail = piu.envelope(spec)
    m = np.arange(Mcols, dtype=np.float64)
    invm = 1.0 / np.sqrt(2.0 * m + 1.0)
    kk = np.arange(Mrows + Mcols + 2, dtype=np.float64)
    psi_k = C / (kk**2 + A**2) ** p
    cum = np.concatenate(([0.0], np.cumsum(psi_k[1:])))
    inv_sqrt_last = 1.0 / np.sqrt(2.0 * Mcols + 1.0)
    best = 0.0
    for n in range(0, Mrows):
        nrm_n = 1.0 / np.sqrt(2.0 * n + 1.0)
        r = 2.0 * cum[n] / (2.0 * n + 1.0)                   # H1 diagonal (empty at n=0)
        d = np.abs(m - n)
        h2 = psi_k[d.astype(np.int64)] * invm
        h2[n] = 0.0
        h3 = psi_k[(n + m + 1).astype(np.int64)] * invm
        tail = (float(int_tail(Mcols - n - 1)) + float(int_tail(Mcols + n))) * inv_sqrt_last
        r += nrm_n * (float(h2.sum()) + float(h3.sum()) + tail)
        best = max(best, r)
    return max(best, float(piu.row_env(Mrows, spec)))


L_RATIO = {}      # (a, b) -> (L_w0, L_base) bookkeeping for the report


def L_seg(a, b):
    Lw = schur_full_seg_w0(a, b, 1)
    Lb = schur_full_seg(a, b, 1)
    chk("weighted-row-0 L <= base L", Lw <= Lb + 1e-15, f"cell [{a},{b}] {Lw} vs {Lb}")
    L_RATIO[(a, b)] = (Lw, Lb)
    return Lw


def _L_evidence(cells, tag):
    """float64 EVIDENCE (not in any bound): L_w0 must dominate ||H'(w)|| at N=1500
    on sampled cells; and report how much tighter it is than the base L."""
    for c in cells[::12] + [cells[-1]]:
        a, b = c[0], c[1]
        for wq in (a, 0.5 * (a + b), b):
            hp = float(np.max(np.abs(np.linalg.eigvalsh(H_mat(1500, wq, 1)))))
            chk("evidence: L_w0 >= ||H'(w)|| at N=1500", c[3] >= hp,
                f"w={wq:.4f} L={c[3]:.4f} ||H'||={hp:.4f}")
    ratios = [L_RATIO[(c[0], c[1])][1] / L_RATIO[(c[0], c[1])][0] for c in cells]
    print(f"  [{tag}] L = weighted-row-0 Schur segment bound; base L / L_w0 in "
          f"[{min(ratios):.2f}, {max(ratios):.2f}] on these cells (reviewer: 1.6x .. 5x)")


# ---------------------------------------------------------------------------
# Grid: SAME geometric ratio as the certified base grid (246 points on [0.3,5]),
# so h/wbar is the constant the base file justified against L = sup||H'||.
# ---------------------------------------------------------------------------
BASE_RATIO = (5.0 / 0.3) ** (1.0 / 245)
LO_EXT, LO_BASE = 0.2, 0.3
HI_BASE, HI_EXT = 5.0, 20.0
MUTILDES = (0.0, 0.1, 0.13, 0.26)
NEW_NODES = (7.0, 10.0, 15.0, 20.0)


def ext_grid(lo, hi):
    g = [lo]
    while g[-1] * BASE_RATIO < hi - 1e-12:
        g.append(g[-1] * BASE_RATIO)
    g.append(hi)
    return [float(x) for x in g]


def head_size(b):
    """Banded head size for the alpha certificate, keyed on the cell's RIGHT end
    (the larger wbar, where the tail envelopes are largest)."""
    if b <= 5.0 + 1e-12:
        return 70
    if b <= 10.0 + 1e-12:
        return 140
    return 280


def alpha_cell_ub_Na(a, b, cA, cB, L, Na):
    """`piu.alpha_cell_ub` with an explicit head size (the base function reads the
    module constant NA_HEAD=70).  Same three terms:
        sup_{[a,b]} lam_2(head) <= max(cA,cB) + L (b-a)/2      (endpoint Cholesky
                                                              + Lipschitz slack)
        + ||B||_seg + ||D||_seg                               (segment envelopes)
        + ROUND_PAD."""
    sp = spec_seg(a, b, 0)
    nB = schur_coupling(Na, b, 0, spec=sp)
    nD = schur_tailblock(Na, b, 0, spec=sp)
    head = max(cA, cB) + mp.mpf(L) * (mp.mpf(b) - mp.mpf(a)) / 2
    return head + nB + nD + ROUND_PAD


def alpha_side(grid):
    """Per-point Cholesky c (banded Na) and per-cell alpha_cell.  mutilde-free."""
    cs, nas = [], []
    for i, w in enumerate(grid):
        Na = head_size(grid[min(i + 1, len(grid) - 1)])   # right neighbour's band
        Na = max(Na, head_size(w))
        _, ok, c, _, _ = alpha_certificate_H(w, Na=Na)
        chk("Cholesky-deflation certificate SPD on the extension grid", ok,
            f"w={w:.4f} Na={Na}")
        cs.append(c)
        nas.append(Na)
    cells = []
    for i in range(len(grid) - 1):
        a, b = grid[i], grid[i + 1]
        Na = min(nas[i], nas[i + 1])      # both endpoint c's were certified at >= Na
        # c is a bound on lam_2 of the Na_i-head; a LARGER head has a larger
        # lam_2 (Cauchy), so a c certified at Na_i is NOT automatically valid at
        # a larger Na.  Using the smaller of the two heads keeps the coupling
        # envelopes consistent with the head actually bounded.  Where the two
        # differ (band boundary) we simply recertify both ends at the smaller Na.
        cA, cB = cs[i], cs[i + 1]
        if nas[i] != nas[i + 1]:
            _, okA, cA, _, _ = alpha_certificate_H(a, Na=Na)
            _, okB, cB, _, _ = alpha_certificate_H(b, Na=Na)
            chk("band-boundary recertification SPD", okA and okB, f"cell [{a},{b}]")
        L = L_seg(a, b)
        al = alpha_cell_ub_Na(a, b, cA, cB, L, Na)
        cells.append((a, b, Na, L, al))
    return cells


def cell_rows(cells, mt):
    rows = []
    for a, b, Na, L, al in cells:
        if mt == 0.0:
            mu2, k1, k2, k3 = mu2_cell_lb(a, b)
            case = 0
        else:
            mu2, k1, k2, k3, case = mu2_pencil_cell(a, b, mt)
        g = float(mu2 - al)
        ks = bool(k1 >= 0 and k2 >= 0 and k3 >= 0)
        rows.append((a, b, Na, L, float(al), float(mu2), g, ks, case))
    return rows


def true_lam12(N, w, mt=0.0, Nc=0):
    E = np.linalg.eigvalsh(M_mat(N, w, 0, mt, Nc))
    return float(E[-1]), float(E[-2])


# ============================================================================
G_LO = {}


def run_X0():
    print("=== STAGE X0: LOWER extension [0.2, 0.3] (T5 domain bookkeeping) ===")
    t0 = time.time()
    grid = ext_grid(LO_EXT, LO_BASE)
    print(f"  grid: {len(grid)} points, {len(grid)-1} geometric cells at the base "
          f"ratio {BASE_RATIO:.5f}; head Na=70 (unchanged)")
    chk("lower grid spans exactly [0.2,0.3]",
        abs(grid[0] - LO_EXT) < 1e-12 and abs(grid[-1] - LO_BASE) < 1e-12)
    cells = alpha_side(grid)
    chk("lower extension uses the base head size only", all(c[2] == 70 for c in cells))
    _L_evidence(cells, "X0")
    for mt in MUTILDES:
        rows = cell_rows(cells, mt)
        good = [r for r in rows if r[6] > 0 and r[7]]
        ginf = min(r[6] for r in rows)
        worst = min(rows, key=lambda r: r[6])
        if mt == 0.0:
            for r in rows:
                chk("mu*=0 lower cell: g_alt > 0", r[6] > 0,
                    f"cell [{r[0]:.4f},{r[1]:.4f}] g={r[6]:.6f}")
                chk("mu*=0 lower cell: k-signs hold", r[7],
                    f"cell [{r[0]:.4f},{r[1]:.4f}]")
            G_LO[mt] = ginf
            # the join with the base theorem: cell [0.2989,0.3] touches 0.3 exactly
            chk("lower extension joins the base interval at 0.3",
                abs(rows[-1][1] - 0.3) < 1e-12)
            print(f"  mu*=0     : ALL {len(rows)} cells certified; "
                  f"inf g_alt = {ginf:.6f} at [{worst[0]:.4f},{worst[1]:.4f}] "
                  f"(mu_2={worst[5]:.5f}, alpha_cell={worst[4]:.5f})")
            # soundness evidence (float64, NOT in the bound)
            for r in rows[::7]:
                for wq in (r[0], 0.5 * (r[0] + r[1]), r[1]):
                    l1, l2 = true_lam12(600, wq)
                    chk("evidence: alpha_cell >= true lam_2 (N=600)", r[4] >= l2 - 1e-9,
                        f"w={wq:.4f}")
                    chk("evidence: mu_2 <= true lam_1 (N=600)", r[5] <= l1 + 1e-9,
                        f"w={wq:.4f}")
        else:
            # HONEST: the pencil's 2x2 floor is collapsed here; no cell should pass.
            chk("mutilde>0: no lower cell is certified (collapsed prefix intact)",
                len(good) == 0, f"mt={mt} good={len(good)}")
            print(f"  mutilde={mt:<5}: 0/{len(rows)} cells certified (inf g="
                  f"{ginf:+.4f}, k-signs fail on all) -> wbar_min(mutilde) UNCHANGED")
    print(f"  ==> mu*=0 theorem domain extends to [0.2, 5]:  "
          f"inf_{{[0.2,0.3]}} g_alt >= {G_LO[0.0]:.6f}  (base [0.3,5]: 0.093225)")
    _summary("Stage X0", f"({time.time()-t0:.1f}s)")


# ============================================================================
G_HI = {}      # mutilde -> (wbar_max, inf g_alt on [5, wbar_max], n cells)
_CELLS_HI = []


def run_X1():
    print("\n=== STAGE X1: UPPER extension [5, 20], banded head size ===")
    t0 = time.time()
    grid = ext_grid(HI_BASE, HI_EXT)
    print(f"  grid: {len(grid)} points, {len(grid)-1} geometric cells at the base "
          f"ratio; Na = 140 for b<=10, 280 for b<=20  (mpmath Cholesky, dps={mp.dps})")
    chk("upper grid spans exactly [5,20]",
        abs(grid[0] - HI_BASE) < 1e-12 and abs(grid[-1] - HI_EXT) < 1e-12)
    cells = alpha_side(grid)
    _CELLS_HI[:] = cells
    _L_evidence(cells, "X1")
    print(f"  alpha side done ({time.time()-t0:.0f}s): "
          f"alpha_cell at [5,{cells[0][1]:.3f}] = {float(cells[0][4]):.4f}, "
          f"at [{cells[-1][0]:.3f},20] = {float(cells[-1][4]):.4f}; "
          f"L_seg in [{min(c[3] for c in cells):.3f}, {max(c[3] for c in cells):.3f}]")
    print(f"  {'mutilde':>8}  {'certified':>18}  {'cells':>8}  {'inf g_alt':>10}  "
          f"{'worst cell':>20}  {'min mu_2':>9}")
    for mt in MUTILDES:
        rows = cell_rows(cells, mt)
        # certified PREFIX from 5 upward (contiguity with the base interval)
        n_ok = 0
        for r in rows:
            if r[6] > 0 and r[7]:
                n_ok += 1
            else:
                break
        chk("a nonempty certified prefix exists above 5", n_ok > 0, f"mt={mt}")
        kept = rows[:n_ok]
        wmax = kept[-1][1]
        ginf = min(r[6] for r in kept)
        worst = min(kept, key=lambda r: r[6])
        G_HI[mt] = (wmax, ginf, n_ok)
        for r in kept:
            chk("cell g_alt > 0 on the certified prefix", r[6] > 0,
                f"mt={mt} cell [{r[0]:.4f},{r[1]:.4f}]")
            chk("cell k-signs hold on the certified prefix", r[7],
                f"mt={mt} cell [{r[0]:.4f},{r[1]:.4f}]")
        if mt <= 0.13:
            chk("mutilde <= 0.13: the WHOLE [5,20] is certified",
                n_ok == len(rows), f"mt={mt} {n_ok}/{len(rows)}")
        else:
            if n_ok < len(rows):
                chk("mutilde=0.26: the first dropped cell fails for a reason",
                    rows[n_ok][6] <= 0 or not rows[n_ok][7], f"n_ok={n_ok}")
            chk("mutilde=0.26 prefix covers the lambda=0.5, mu*=0.13 solution "
                "(EVIDENCE numbers: wbar_c=13.5 at r=10, 15.6 at r=3)", wmax > 15.6,
                f"wmax={wmax:.4f}")
            # what the BASE (unweighted row 0) L would have given on the same cells
            n_base = 0
            for (a, b, Na, L, al) in cells:
                lw, lb = L_RATIO[(a, b)]
                al_b = float(al) + (lb - lw) * (b - a) / 2.0
                mu2, k1, k2, k3, case = mu2_pencil_cell(a, b, mt)
                if float(mu2) - al_b > 0 and k1 >= 0 and k2 >= 0 and k3 >= 0:
                    n_base += 1
                else:
                    break
            print(f"           (with the base unweighted-row-0 L the mutilde=0.26 prefix "
                  f"would end at {cells[n_base-1][1]:.4f} after {n_base} cells -- the "
                  f"reviewer's 16.3172 cutoff; the weighted row 0 is what closes it)")
        print(f"  {mt:>8.3f}  [5, {wmax:>7.4f}]{'':>5}  {n_ok:>3}/{len(rows):<4}  "
              f"{ginf:>10.6f}  [{worst[0]:.4f},{worst[1]:.4f}]  {min(r[5] for r in kept):>9.5f}")
        # soundness EVIDENCE on retained cells (float64; audits, never enters)
        for r in kept[::10]:
            a, b = r[0], r[1]
            for wq in (a, 0.5 * (a + b), b):
                if mt == 0.0:
                    l1, l2 = true_lam12(1500, wq)
                    chk("evidence: alpha_cell >= true lam_2(H) at N=1500",
                        r[4] >= l2 - 1e-9, f"w={wq:.4f} alpha={r[4]:.5f} lam2={l2:.5f}")
                    chk("evidence: mu_2 <= true lam_1(H^(2))",
                        r[5] <= float(np.linalg.eigvalsh(H_mat(2, wq))[-1]) + 1e-12,
                        f"w={wq:.4f}")
                else:
                    for Nc in (0, 1, 2, 9, 300):
                        v = float(np.linalg.eigvalsh(M_mat(2, wq, 0, mt, Nc))[-1])
                        chk("evidence: pencil mu_2 bound <= true 2x2 top eigenvalue",
                            v >= r[5] - 1e-13, f"mt={mt} w={wq:.4f} Nc={Nc}")
    # covering-lemma consistency: wbar_max nonincreasing in mutilde
    for m1, m2 in zip(MUTILDES[:-1], MUTILDES[1:]):
        chk("covering lemma: wbar_max nonincreasing in mutilde",
            G_HI[m1][0] >= G_HI[m2][0] - 1e-12, f"{m1}->{m2}")
    # the base head Na=70 would NOT have closed this (the point of the banding)
    a, b = cells[-1][0], cells[-1][1]
    cA70 = alpha_certificate_H(a, Na=70)[2]
    cB70 = alpha_certificate_H(b, Na=70)[2]
    al70 = float(alpha_cell_ub_Na(a, b, cA70, cB70, cells[-1][3], 70))
    mu2_last = float(mu2_cell_lb(a, b)[0])
    chk("banding is load-bearing: Na=70 fails on the last cell",
        mu2_last - al70 < 0, f"g(Na=70)={mu2_last-al70:.4f}")
    print(f"  banding check at [{a:.3f},20]: Na=70 gives alpha_cell={al70:.3f} > "
          f"mu_2={mu2_last:.3f} (fails); Na=280 gives {float(cells[-1][4]):.3f} (closes)")
    print(f"  ==> mu*=0: inf_{{[5,20]}} g_alt >= {G_HI[0.0][1]:.6f}; combined with the "
          f"base and X0, the theorem domain is [0.2, 20] with inf gap >= "
          f"{min(G_LO.get(0.0, 1.0), 0.093225, G_HI[0.0][1]):.6f}")
    _summary("Stage X1", f"({time.time()-t0:.1f}s)")


# ============================================================================
NODE = {}


def run_X2():
    print("\n=== STAGE X2: new NODES -- delta, g_alt, B2/B4, K_*, s_max, bound (6) ===")
    t0 = time.time()
    print(f"  {'wbar':>6}  {'Na':>4}  {'Nv':>4}  {'mu<=lam1':>9}  {'alpha>=lam2':>11}  "
          f"{'delta':>8}  {'mu_2':>8}  {'g_alt':>8}  {'B2^inf':>8}  {'B4^inf':>8}  "
          f"{'K_*':>10}  {'s_max':>7}")
    for w in NEW_NODES:
        Na = head_size(w)
        Nv = 2 * Na
        mu = rayleigh_mu(w, Nv=Nv)
        alpha, ok, c, nB, nD = alpha_certificate_H(w, Na=Na)
        chk("node Cholesky SPD", ok, f"w={w}")
        delta = float(mu - alpha)
        mu2 = mu2_cell_lb(w, w)[0]
        galt = float(mu2 - alpha)
        chk("node delta > 0", delta > 0, f"w={w} {delta:.6f}")
        chk("node g_alt > 0", galt > 0, f"w={w} {galt:.6f}")
        l1, l2 = true_lam12(1500, w)
        chk("evidence: mu <= lam1 (N=1500)", float(mu) <= l1 + 1e-9, f"w={w}")
        chk("evidence: alpha >= lam2 (N=1500)", float(alpha) >= l2 - 1e-9, f"w={w}")
        B2 = cert_sup_on_J(w, 2)
        B4 = cert_sup_on_J(w, 4)
        chk("B2^inf, B4^inf finite positive", B2 > 0 and B4 > 0 and np.isfinite(B2 + B4))
        s_max = min(SIGMA0, 0.95 * np.sqrt(galt / (2.0 * w**2 * B2)))
        K_tail = w**4 * B4 * ETA_MAX**2 / 24.0
        K_2nd = 0.5 * w**4 * B2**2 / galt
        K = K_tail + K_2nd
        chk("uniform-denominator condition at s_max", w**2 * B2 * s_max**2 < galt / 2.0,
            f"w={w}")
        NODE[w] = dict(Na=Na, mu=float(mu), alpha=float(alpha), delta=delta,
                       mu2=float(mu2), galt=galt, B2=B2, B4=B4, K=K, s_max=s_max)
        print(f"  {w:>6.1f}  {Na:>4}  {Nv:>4}  {float(mu):>9.5f}  {float(alpha):>11.5f}  "
              f"{delta:>8.5f}  {float(mu2):>8.5f}  {galt:>8.5f}  {B2:>8.4f}  {B4:>8.4f}  "
              f"{K:>10.4e}  {s_max:>7.4f}")
    # the certified interval infimum must sit below every nodal g_alt in its cell
    if G_HI:
        for w in NEW_NODES:
            chk("interval infimum <= nodal g_alt", G_HI[0.0][1] <= NODE[w]["galt"] + 1e-12,
                f"w={w}")

    print(f"\n  finite-N consistency of bound (6) at the new nodes (Stage-5 protocol of "
          f"the base file; N=400 proxy, sigma <= s_max):")
    print(f"    {'wbar':>6}  {'shape':>9}  {'max|rem|/s^4':>13}  {'K_*':>10}  "
          f"{'ratio':>7}  {'r(.01)/r(max)':>13}  {'settle':>7}")
    worst_ratio = 0.0
    for w in NEW_NODES:
        kap, c3, hN = kappa_c3(400, w)
        K = NODE[w]["K"]
        sigmas = tuple(s for s in (0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25)
                       if s <= NODE[w]["s_max"])
        for name, eta, p in _fixed_shapes():
            g3 = float(p @ eta**3)
            per = []
            for s in sigmas:
                exact = frak_k_of_P(400, w, s * eta, p)
                pred = hN * (1.0 + 0.5 * kap * s**2 + c3 * g3 * s**3)
                per.append(abs(exact - pred) / s**4)
            worst = max(per)
            r = worst / K
            worst_ratio = max(worst_ratio, r)
            spread = per[0] / per[-1]
            settle = abs(per[0] - per[1]) / per[1]
            chk("finite-N remainder <= K_* sigma^4 at N=400", r <= 1.0,
                f"w={w} {name}: {worst:.4e} > {K:.4e}")
            chk("|rem|/sigma^4 does not diverge as sigma->0", spread < 30.0,
                f"w={w} {name}: {spread:.2f}")
            chk("|rem|/sigma^4 settles as sigma->0", settle < 0.25,
                f"w={w} {name}: {settle:.3f}")
            print(f"    {w:>6.1f}  {name:>9}  {worst:>13.4e}  {K:>10.4e}  {r:>7.4f}  "
                  f"{spread:>13.2f}  {settle:>7.4f}")
    chk("K_* not vacuous at the new nodes", worst_ratio > 1e-4, f"{worst_ratio:.2e}")
    print(f"  worst |rem|/(K_* sigma^4) over the new nodes = {worst_ratio:.4f}")
    # N-independence of the measured remainder (the proxy justification), with
    # N pushed to 800 because the top eigenvector spreads ~ wbar
    print("\n  N-independence of |rem|/sigma^4 (3pk-asym, sigma=0.01), N in {200,400,800}:")
    name, eta, p = [s for s in _fixed_shapes() if s[0] == "3pk-asym"][0]
    g3 = float(p @ eta**3)
    for w in NEW_NODES:
        vals = []
        for N in (200, 400, 800):
            kap, c3, hN = kappa_c3(N, w)
            exact = frak_k_of_P(N, w, 0.01 * eta, p)
            pred = hN * (1.0 + 0.5 * kap * 1e-4 + c3 * g3 * 1e-6)
            vals.append(abs(exact - pred) / 1e-8)
        rel = (max(vals) - min(vals)) / max(vals)
        chk("|rem|/sigma^4 N-independent for N in [200,800]", rel < 1e-2,
            f"w={w} spread={rel:.2e}")
        print(f"    {w:>6.1f}  " + "  ".join(f"{v:>12.6f}" for v in vals) + f"   spread={rel:.1e}")
    _summary("Stage X2", f"({time.time()-t0:.1f}s)")


# ============================================================================
def certified_interval(mt):
    """[wbar_min, wbar_max] certified for this mutilde by base + X0 + X1.
    wbar_min: 0.2 (mu*=0, X0) or the Stage-E value of `p4d_coulomb_certify`
    (0.5391 / 0.6427 / 1.1686, quoted -- NOT recomputed here; those are the
    published certified values and X0 confirms nothing below 0.3 is added)."""
    wmin = {0.0: LO_EXT, 0.1: 0.5391, 0.13: 0.6427, 0.26: 1.1686}
    key = min(m for m in MUTILDES if m >= mt - 1e-12) if mt <= 0.26 else None
    if key is None:
        return None
    if G_HI:
        return wmin[key], G_HI[key][0]
    # X3 run standalone: quote the X1 result (full run uses the live value)
    return wmin[key], X1_QUOTED_WMAX[key]


X1_QUOTED_WMAX = {0.0: 20.0, 0.1: 20.0, 0.13: 20.0, 0.26: 20.0}   # X1 v2 (weighted row 0), 2026-09-11


def run_X3():
    print("\n=== STAGE X3: physics anchor (EVIDENCE) -- wbar_c(lambda, mu*) vs the "
          "certified interval ===")
    t0 = time.time()
    N = 600
    ws = np.concatenate([np.arange(0.2, 5, 0.02), np.arange(5, 60.0, 0.1)])
    print("  Tc := sup{T : lam_1(M) >= 1/lambda}  <=>  wbar_c = smallest wbar with "
          "lam_1(H - mutilde K_c(N_c(wbar))) >= 1/lambda,  N_c = floor(r wbar + 1/2)")
    print("  tag: IN = inside the certified interval for mutilde=mu*/lambda; "
          "OUT = outside; -- = mutilde > 0.26 (no node) or no Tc up to wbar=60")
    if not G_HI:
        print("  (standalone run: wbar_max quoted from X1 -- "
              f"{X1_QUOTED_WMAX}; the default full run uses the live X1 values)")
    covered_now, covered_before = 0, 0
    for r in (10, 3):
        print(f"\n  r = w_c/w_E = {r}")
        print(f"  {'lambda':>7} | " + " | ".join(f"{('mu*='+str(m)):>16}" for m in (0.0, 0.05, 0.10, 0.13)))
        for lam in (0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0):
            cells = []
            for mus in (0.0, 0.05, 0.10, 0.13):
                mt = mus / lam
                wc = np.nan
                for w in ws:
                    Nc = int(np.floor(r * w + 0.5))
                    f = np.linalg.eigvalsh(H_mat(N, w) - mt * K_c(N, Nc))[-1]
                    if f >= 1.0 / lam:
                        wc = float(w)
                        break
                iv = certified_interval(mt)
                if np.isnan(wc) or iv is None:
                    tag = "--"
                else:
                    tag = "IN" if iv[0] <= wc <= iv[1] else "OUT"
                    if tag == "IN":
                        covered_now += 1
                    if iv[0] <= wc <= HI_BASE:
                        covered_before += 1
                cells.append(f"{wc:6.2f} {tag:>3} mt={mt:4.2f}")
            print(f"  {lam:>7.1f} | " + " | ".join(f"{c:>16}" for c in cells))
    chk("extension strictly enlarges physical coverage", covered_now > covered_before,
        f"{covered_now} vs {covered_before}")
    print(f"\n  ==> (lambda, mu*) entries with wbar_c inside the certified interval: "
          f"{covered_now} now vs {covered_before} with the old upper end 5 "
          f"(out of 64 entries, r in {{10,3}})")
    _summary("Stage X3", f"({time.time()-t0:.1f}s)")


def run_verdict():
    print("\n=== VERDICT ===")
    g_all = min(G_LO[0.0], 0.093225, G_HI[0.0][1])
    where = ("[0.2,0.3] (X0)" if g_all == G_LO[0.0] else
             "[5,20] (X1)" if g_all == G_HI[0.0][1] else "[0.3,5] (base, 0.093225)")
    print(f"  mu*=0        : theorem domain [0.2, 20], inf_{{N>=2}} gap_N >= {g_all:.6f} "
          f"on the whole interval; binding segment: {where}.  Per segment: "
          f"[0.2,0.3] {G_LO[0.0]:.6f} | [0.3,5] 0.093225 (base) | [5,20] {G_HI[0.0][1]:.6f}")
    for mt in (0.1, 0.13, 0.26):
        iv = certified_interval(mt)
        print(f"  mutilde={mt:<5}: certified [{iv[0]:.4f}, {iv[1]:.4f}]  "
              f"(inf g_alt on [5,wbar_max] = {G_HI[mt][1]:.6f})")
    print("  new nodes with K_* and s_max: " + ", ".join(f"{w:g}" for w in NEW_NODES))
    print("  NOT done / honest edges: extension stops at 20 for COST (Na^3 Cholesky); "
          "beyond ~20 the binding term is the Lipschitz slack L*h/2 (L saturates ~1.2 "
          "with the weighted row 0, h/wbar fixed => slack ~ 0.007 wbar vs mu_2 -> 1.66), "
          "so a finer grid or a sharper L is needed well before Na cost bites; "
          "lambda<0.5 with mu*=0.13 (mutilde>0.26) has no node at all.")


if __name__ == "__main__":
    want = [s.upper() for s in sys.argv[1:]] or ["X0", "X1", "X2", "X3"]
    T0 = time.time()
    if "X0" in want:
        run_X0()
    if "X1" in want:
        run_X1()
    if "X2" in want:
        run_X2()
    if "X3" in want:
        run_X3()
    if set(want) >= {"X0", "X1"}:
        run_verdict()
    print(f"\nALL {PASS[0]} checks passed  ({time.time()-T0:.1f}s).")
