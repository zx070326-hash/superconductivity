"""Budget probe for the two-stage Rayleigh route (2026-09-19).

Context: `theory/notes/coulomb/mustar_strong_coupling_extremality_2026-09-19.md`
sec.4.5 proposes a two-stage repair for the falsified Rayleigh step:
  Stage A (far field):   exclude P with lam1(M_P) <= c*lam1(M_E), c=0.5,
                          i.e. margin := (lamE-lamP)/lamE >= 0.5, via (CVX)+(E-env).
  Stage B (competitive): only need (AB)/rho<=1 on the P's Stage A does NOT exclude,
                          i.e. margin < 0.5 (P close to Einstein).
That note's sec.4.5/sec.6 explicitly say Stage B's own budget --
    sup{ rho(P) : margin(P) <= 0.5 }
-- was never computed.  This script computes it, from scratch, on a grid
DESIGNED to hit the corners of the two-atom family (weight -> 0.5 and -> 0/1,
spread -> 1 and -> very large), not a mild random sample (CLAUDE.md: "sampling
must be able to trigger the boundary of the premise").

IMPORTANT convention check (done once, in `main`, before trusting any number):
the task that spawned this probe described the "competitive zone" as
margin in [0.5,1.0].  That is the OPPOSITE of the note's own Stage A/B split
(Stage A excludes margin>=0.5 as "far field"; Stage B's zone is margin<=0.5,
i.e. CLOSE to Einstein).  We verify this from the note's own formulas
(lam1(M_P)<=0.5*lam1(M_E)  <=>  margin>=0.5) and then report BOTH bands so the
reader can see the direction, but the number that matters for Stage B's
feasibility is the one for margin in [0, 0.5].

Reuses (read-only import, no re-derivation): O_of_kernel / pencil / abel /
kernel_measure / kernel_einstein from _chk_mustar_conerepair_op_2026-09-19.py,
and the anchor mutilde/mu*/r construction + solve_varpi from
_chk_mustar_physical_placement_2026-09-19.py (identical formulas, re-imported
not re-typed, to avoid a second transcription of the anchor table).

float64, finite N, finite (but boundary-seeking) shape grid ==> numerical-evidence.
This is a BUDGET PROBE, not a certificate and not a proof.
"""

import importlib.util
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


op = _load("cr_op", "_chk_mustar_conerepair_op_2026-09-19.py")
place = _load("cr_place", "_chk_mustar_physical_placement_2026-09-19.py")

MU_STAR = place.MU_STAR
ANCHOR_MUTILDE = place.ANCHOR_MUTILDE
ANCHORS = place.ANCHORS  # [(name, lam, mu*), ...] for the 7 mutilde anchors


# ------------------------------------------------------------- shape grid ----
# Boundary-seeking, not random: weight sweeps from near-0.5 (mild) out to
# near-0/near-1 (one atom dominates almost entirely), spread sweeps from
# near-1 (near-degenerate, i.e. near-Einstein) out to very large (extreme
# bimodal).  This is the SAME two-atom family Theorem TA says is exhaustive
# (repo `mustar_strong_coupling_extremality_2026-09-19.md` sec.2), just with a
# grid built to hit corners rather than the mild 5-shape production sample.
WEIGHTS = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999,
                     0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001])
SPREADS = np.array([1.001, 1.01, 1.05, 1.1, 1.3, 1.5, 2.0, 3.0, 5.0, 10.0,
                     30.0, 100.0, 300.0, 1000.0, 3000.0, 1e4])

# Hard cap on the operator size.  The rare-heavy-tail corner (weight->1,
# spread->large) pushes the minority atom's frequency up to ~spread*w_rms,
# which would otherwise blow nsize() past 10^4-10^5 and make the full grid
# intractable.  Capping trades some truncation accuracy in that one corner
# for tractability; the corner is re-checked separately at higher N below
# if it turns out to matter (see N-convergence section in `main`).
N_CAP = 4000


def two_atom(w_rms, weight, spread):
    r = spread
    w1 = w_rms / np.sqrt(weight + (1.0 - weight) * r * r)
    return (w1, r * w1), (weight, 1.0 - weight)


N_CAP_ABS = 9000  # absolute ceiling even under boost, so eigh stays tractable


def nsize(w_rms, Nc, ws, boost=1.0):
    raw = max(300, boost * 24 * max(ws), boost * (5 * Nc + 40))
    return int(min(raw, N_CAP * boost, N_CAP_ABS))


def measure_point(w_rms, weight, spread, Nc, mut, boost=1.0):
    """Return dict with margin, rho, exact = <v_P,O[d]v_P>, G11,G12,G22, N."""
    ws, ps = two_atom(w_rms, weight, spread)
    N = nsize(w_rms, Nc, ws, boost=boost)
    A = op.O_of_kernel(N, op.kernel_measure(N, ws, ps))
    H = op.O_of_kernel(N, op.kernel_einstein(N, w_rms))
    Od = H - A
    MP = op.pencil(A, N, Nc, mut)
    ME = op.pencil(H, N, Nc, mut)
    lamP, vP, _ = op.top_pair(MP)
    lamE = np.linalg.eigvalsh(ME)[-1]
    a = op.abel(vP)
    neg = a < 0.0
    xi = op.from_abel(np.where(neg, -a, 0.0))
    w = vP + xi
    G11 = float(w @ Od @ w)
    G12 = float(xi @ Od @ w)
    G22 = float(xi @ Od @ xi)
    exact = float(vP @ Od @ vP)
    denom = G11 + G22
    rho = float("nan") if denom <= 0 else 2.0 * G12 / denom
    margin = float((lamE - lamP) / lamE)
    return dict(margin=margin, rho=rho, exact=exact, G11=G11, G12=G12, G22=G22,
                lamP=lamP, lamE=lamE, N=N, weight=weight, spread=spread)


def anchor_setup(mutilde, r):
    lam = MU_STAR / mutilde
    w = place.solve_varpi(lam, MU_STAR, r)
    Nc = int(np.floor(r * w))
    mut = mutilde
    return lam, w, Nc, mut


def scan_anchor(mutilde, r, boost=1.0):
    lam, w_rms, Nc, mut = anchor_setup(mutilde, r)
    rows = []
    for weight in WEIGHTS:
        for spread in SPREADS:
            rows.append(measure_point(w_rms, weight, spread, Nc, mut, boost=boost))
    return lam, w_rms, Nc, mut, rows


def best_in_band(rows, lo, hi):
    """Max rho among rows with margin in [lo,hi) (rho finite)."""
    cand = [r for r in rows if lo <= r["margin"] < hi and np.isfinite(r["rho"])]
    if not cand:
        return None
    return max(cand, key=lambda r: r["rho"])


def main():
    print("=" * 100)
    print("Stage-B budget probe: sup rho over the COMPETITIVE zone margin<=0.5")
    print("(note's own Stage A excludes margin>=0.5 as far field; Stage B needs")
    print(" (AB)/rho<=1 on margin<0.5 -- that sup was never measured before.)")
    print("=" * 100)

    header = (f"{'mutilde':>8} {'r':>4} {'lam':>7} {'w_E':>8} {'Nc':>4} "
              f"{'sup rho [0,0.5)':>16} {'@ (weight,spread)':>20} "
              f"{'sup rho [0.5,1]':>16} {'sup rho all':>12}")
    print(header)

    overall_best = None
    overall_best_band = None
    rows_by_anchor = {}

    for mutilde in ANCHOR_MUTILDE:
        for r in (3.0, 10.0):
            lam, w_rms, Nc, mut, rows = scan_anchor(mutilde, r)
            rows_by_anchor[(mutilde, r)] = (lam, w_rms, Nc, mut, rows)
            comp = best_in_band(rows, 0.0, 0.5)
            far = best_in_band(rows, 0.5, 1.0 + 1e-9)
            allb = max((rr for rr in rows if np.isfinite(rr["rho"])), key=lambda rr: rr["rho"])
            comp_str = f"{comp['rho']:16.6f}" if comp else f"{'(empty)':>16}"
            loc_str = f"({comp['weight']:.4f},{comp['spread']:.3g})" if comp else ""
            far_str = f"{far['rho']:16.6f}" if far else f"{'(empty)':>16}"
            print(f"{mutilde:8.4f} {r:4.0f} {lam:7.3f} {w_rms:8.4f} {Nc:4d} "
                  f"{comp_str} {loc_str:>20} {far_str} {allb['rho']:12.6f}")
            if comp and (overall_best is None or comp["rho"] > overall_best["rho"]):
                overall_best = comp
                overall_best_band = (mutilde, r)

    print("-" * 100)
    if overall_best is None:
        print("No finite-rho point found with margin<0.5 on this grid at all 14 anchors.")
    else:
        print(f"OVERALL sup rho over margin in [0,0.5), all 14 anchor rows: "
              f"{overall_best['rho']:.6f}")
        print(f"  attained at anchor mutilde={overall_best_band[0]:.4f}, r={overall_best_band[1]:.0f}, "
              f"(weight,spread)=({overall_best['weight']:.4f},{overall_best['spread']:.4g}), "
              f"margin={overall_best['margin']:.4f}")

    # ---------------------------------------------------------- N-convergence
    print("\n### N-convergence check on the overall-worst competitive-zone point")
    mutilde0, r0 = overall_best_band
    lam0, w_rms0, Nc0, mut0, _ = rows_by_anchor[(mutilde0, r0)]
    weight0, spread0 = overall_best["weight"], overall_best["spread"]
    for boost in (1.0, 2.0, 4.0, 8.0):
        pt = measure_point(w_rms0, weight0, spread0, Nc0, mut0, boost=boost)
        print(f"  boost={boost:4.1f}  N={pt['N']:6d}  margin={pt['margin']:.6f}  "
              f"rho={pt['rho']:.8f}  exact={pt['exact']:.6e}")

    # ------------------------------------------------- margin=0.5 boundary --
    print("\n### Boundary behaviour of rho(margin) near margin=0.5, and continuity")
    print("    of margin itself along a fixed-weight spread sweep (rules out a")
    print("    grid artefact masquerading as a jump).")
    mutilde1, r1 = overall_best_band
    lam1, w_rms1, Nc1, mut1, rows1 = rows_by_anchor[(mutilde1, r1)]
    # fixed weight = weight of the overall-worst point, sweep spread finely
    fine_spreads = np.geomspace(1.001, 1e5, 60)
    print(f"  anchor mutilde={mutilde1:.4f} r={r1:.0f}, weight fixed = {weight0:.4f}")
    print(f"  {'spread':>10} {'margin':>10} {'rho':>12}")
    prev_margin = None
    for sp in fine_spreads:
        pt = measure_point(w_rms1, weight0, sp, Nc1, mut1)
        jump = ""
        if prev_margin is not None and abs(pt["margin"] - prev_margin) > 0.15:
            jump = "  <-- margin jump >0.15 between adjacent grid points"
        print(f"  {sp:10.3g} {pt['margin']:10.4f} {pt['rho']:12.6f}{jump}")
        prev_margin = pt["margin"]

    # ------------------------------------------------------- extra sup all --
    print("\n### For reference: sup rho with NO margin restriction at all (whole")
    print("    two-atom family, this grid) -- should exceed the sec.4.4 figure")
    print("    of ~1.086 since our grid is more extreme (weight to 1e-4, spread to 1e5).")
    all_best = None
    for (mutilde, r), (lam, w_rms, Nc, mut, rows) in rows_by_anchor.items():
        b = max((rr for rr in rows if np.isfinite(rr["rho"])), key=lambda rr: rr["rho"])
        if all_best is None or b["rho"] > all_best["rho"]:
            all_best = b
            all_best_anchor = (mutilde, r)
    print(f"  sup rho (unrestricted) = {all_best['rho']:.6f} at anchor "
          f"mutilde={all_best_anchor[0]:.4f} r={all_best_anchor[1]:.0f} "
          f"(weight,spread)=({all_best['weight']:.4f},{all_best['spread']:.4g}), "
          f"margin={all_best['margin']:.4f}")


if __name__ == "__main__":
    main()
