"""Locate the weak-coupling extremality-failure boundary lambda_c(mu*) by bisection.

First run (mustar_physical_line.py, /zx/mustar_line/run.log) bracketed it:
  mu*_rms=0.2023: no at 0.52, YES at 0.50
  mu*_rms=0.15  : no at 0.360, YES at 0.330
  mu*_rms=0.13  : no at 0.312, YES at 0.286
i.e. lambda_c ~ 2.2-2.5 mu*.  Bisect s = lambda/mu* on each line to 1% and add
mu* = 0.10 and 0.18 to fit the law lambda_c = c(mu*) mu*.

Cost note: varpi_E blows up as lambda falls (27855 at lambda=0.26, mu*=0.13), so
bisection is capped at 12 steps and any row needing varpi_E > 6e4 is reported
rather than pursued.  Strictness: numerical-evidence (float64), Nfac=10 during
bisection with an Nfac=20 confirmation on the final bracket.
"""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
from mustar_physical_line import row

R = 5.0
LNR = math.log(R)


def bisect(mus, s_lo, s_hi, steps=12, Nfac=10):
    """s = lambda/mu*; s_lo must violate (YES), s_hi must not."""
    mu = mus / (1.0 - mus * LNR)
    print(f"\n=== mu*_rms={mus} (bare mu={mu:.4f})  bisect s=lambda/mu* in "
          f"[{s_lo},{s_hi}] ===", flush=True)
    hint = None
    for i in range(steps):
        s = 0.5 * (s_lo + s_hi)
        res = row(round(s * mus, 5), mu, R, Nfac, vp_hint=hint,
                  log=lambda t: print(f"  s={s:.4f}" + t, flush=True))
        if res is None:
            print(f"  s={s:.4f}: no bracket, stop", flush=True)
            break
        hint = res["vpE"]
        if res["vpE"] > 6e4:
            print(f"  s={s:.4f}: varpi_E={res['vpE']:.0f} > 6e4, stop", flush=True)
            break
        if res["viol"]:
            s_lo = s
        else:
            s_hi = s
        if s_hi - s_lo < 0.01 * s_lo:
            break
    print(f"  -> lambda_c/mu* in [{s_lo:.4f}, {s_hi:.4f}]  "
          f"lambda_c ~ {0.5*(s_lo+s_hi)*mus:.4f}", flush=True)
    # confirm the final bracket at Nfac=20
    for s, tag in ((s_lo, "YES-side"), (s_hi, "no-side")):
        row(round(s * mus, 5), mu, R, 20, vp_hint=hint,
            log=lambda t: print(f"  CONFIRM {tag} s={s:.4f}" + t, flush=True))
    return s_lo, s_hi


def main():
    out = {}
    for mus, lo, hi in ((0.2023, 2.40, 2.60), (0.15, 2.20, 2.40),
                        (0.13, 2.20, 2.40), (0.18, 2.20, 2.55), (0.10, 2.10, 2.50)):
        out[mus] = bisect(mus, lo, hi)
    print("\n=== summary: lambda_c / mu*_rms ===", flush=True)
    for mus, (lo, hi) in sorted(out.items()):
        print(f"  mu*={mus:7.4f}   s_c in [{lo:.4f}, {hi:.4f}]   "
              f"lambda_c ~ {0.5*(lo+hi)*mus:.4f}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
