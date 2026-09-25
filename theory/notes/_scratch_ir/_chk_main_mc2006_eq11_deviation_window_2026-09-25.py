#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Moussa-Cohen 2006 (PRB 74, 094520) eq.(11) -- quantify the authors' own stated
deviation over the lambda window where WE consume it (hydrides, lambda ~ 2-3.3).

Consumption point: theory/notes/ceiling/prescriptive_protocol_and_moussa_cohen_priority_2026-09-24.md
                   section 2.2-(3'') table + conclusions (a)/(b).

Verbatim from p.2 right column (confirmed by TWO instruments, both generated
2026-09-25 by the main session:
  _scratch_ir/_mc2006_p02_mine_utf8.txt  <- pdftotext -enc UTF-8 -f 2 -l 2
  _scratch_ir/_mc2006_mine_p-02.png      <- pdftoppm -r 300 -f 2 -l 2 -png ):

  eq.(10)  f(lam<w^2>/w_max^2, mu*) w_max  <=  Tc  <=  f(lam,mu*) sqrt(<w^2>)
  eq.(11)  f(lam,mu*) = 0.69 exp(-(1+lam)/(lam-mu*)) sqrt((1+0.52 lam)/(1+5.6 mu*))
           "for lam > mu*, otherwise f = 0."
  "...can be empirically fit to"
  "The coefficients were chosen based on a least squares fit on the interval
   0.4 <= lam <= 20 and 0 <= mu* <= 0.2."
  "The average root-mean-squared deviation of f is 0.01 and the max deviation is 0.03."

RIGOUR LEVEL of what this script produces: numerical-evidence (mpmath dps=30
point evaluation, NOT interval arithmetic).  The type judgement "a two-sided fit
is not an upper bound" is analytic and does not depend on anything computed here.

DIRECTION AUDIT (load-bearing -- see CLAUDE.md one-sided bound rule):
  The 0.03 is stated WITHOUT a sign.  That establishes "no direction guarantee";
  it does NOT establish that eq.(11) is ever actually below the true bound.
  Gate B below checks the sign of every point we can actually compare, and finds
  eq.(11) high everywhere -> loose, not false.  Do not upgrade (a) into a
  falsification claim.

Run:  OMP_NUM_THREADS=1 python _chk_main_mc2006_eq11_deviation_window_2026-09-25.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")

from mpmath import mp, mpf, exp, sqrt

mp.dps = 30

MAX_DEV = mpf('0.03')          # authors' stated max deviation of f
RMS_DEV = mpf('0.01')          # authors' stated rms deviation of f
FIT_LAM = (mpf('0.4'), mpf('20'))
FIT_MU = (mpf('0'), mpf('0.2'))

# our analytic constant: C_inf = (1/2pi) sqrt(g(2))
C_INF_OURS = mpf('0.1827262477')


def f_MC(lam, mu):
    """Moussa-Cohen eq.(11)."""
    lam, mu = mpf(lam), mpf(mu)
    if lam <= mu:
        return mpf(0)
    return (mpf('0.69') * exp(-(1 + lam) / (lam - mu))
            * sqrt((1 + mpf('0.52') * lam) / (1 + mpf('5.6') * mu)))


def C_MC(lam, mu):
    """Algebraic bridge: Tc <= f sqrt(<w^2>) and M2 = lam <w^2>
    => Tc <= (f/sqrt(lam)) sqrt(M2), i.e. C_MC = f/sqrt(lam), same object
    and same information class as our C(lam,mu*)."""
    return f_MC(lam, mu) / sqrt(mpf(lam))


def C_MC_asymptote(mu):
    """lam -> infinity limit of C_MC.  exp(-(1+lam)/(lam-mu)) -> e^{-1},
    sqrt((1+0.52lam)/(1+5.6mu)) / sqrt(lam) -> sqrt(0.52/(1+5.6mu))."""
    return mpf('0.69') * exp(-1) * sqrt(mpf('0.52') / (1 + mpf('5.6') * mpf(mu)))


def main():
    print("=" * 74)
    print("MC2006 eq.(11): stated max deviation 0.03 (absolute, in f, UNSIGNED)")
    print("=" * 74)

    print("\n[TABLE] 0.03 as a fraction of f, across the window we consume\n")
    print("  %-8s %-6s %-13s %-10s %s" % ("lam", "mu*", "f(lam,mu*)", "0.03/f", "C_MC=f/sqrt(lam)"))
    rows = []
    for mu in ['0', '0.10', '0.13']:
        for lam in ['0.4', '1', '2', '2.5', '3.2953', '4', '10', '20']:
            f = f_MC(lam, mu)
            rel = MAX_DEV / f if f > 0 else mpf('inf')
            rows.append((lam, mu, f, rel))
            print("  %-8s %-6s %-13s %-10s %s"
                  % (lam, mu, mp.nstr(f, 6), mp.nstr(rel * 100, 4) + "%",
                     mp.nstr(C_MC(lam, mu), 8)))
        print()

    # ---- GATE A: self-consistency of our reading of eq.(11) against eq.(7) ----
    # If the square root in eq.(11) had been read in the wrong place, this would
    # not land on the 0.183 they print in eq.(7).
    asym0 = C_MC_asymptote(0)
    print("[GATE A] eq.(11) lam->inf at mu*=0 must reproduce eq.(7)'s 0.183")
    print("         0.69*e^-1*sqrt(0.52) = %s" % mp.nstr(asym0, 10))
    print("         eq.(7) printed        = 0.183")
    print("         our C_inf (analytic)  = %s" % mp.nstr(C_INF_OURS, 10))
    print("         MC_asym / C_inf       = %s" % mp.nstr(asym0 / C_INF_OURS, 10))
    gate_a = abs(asym0 - mpf('0.183')) < mpf('0.0005')
    print("         -> %s" % ("PASS" if gate_a else "FAIL"))

    # ---- GATE B: direction audit ----
    # Everywhere we can actually compare, is eq.(11) above or below the truth?
    # The only point where we have an independent ANALYTIC value is lam->inf,
    # mu*=0, where the true sup is C_inf.
    print("\n[GATE B] direction where an independent analytic value exists")
    above = asym0 > C_INF_OURS
    print("         lam->inf, mu*=0: C_MC = %s vs C_inf = %s"
          % (mp.nstr(asym0, 10), mp.nstr(C_INF_OURS, 10)))
    print("         eq.(11) sits %s the true sup, by %s%%"
          % ("ABOVE" if above else "BELOW",
             mp.nstr((asym0 / C_INF_OURS - 1) * 100, 4)))
    print("         -> eq.(11) is LOOSE here, not FALSE.")
    print("         => conclusion (a) must stay 'no direction guarantee';")
    print("            it may NOT be upgraded to 'falsified'.")

    # ---- GATE C: vacuity at the fit interval's own lower endpoint ----
    print("\n[GATE C] is the assertable band f +/- 0.03 vacuous anywhere in")
    print("         the authors' own fit rectangle?")
    for mu in ['0', '0.10', '0.13', '0.2']:
        lam = FIT_LAM[0]           # 0.4, the lower endpoint they state
        f = f_MC(lam, mu)
        print("         lam=0.4, mu*=%-5s f=%-12s f-0.03=%-13s %s"
              % (mu, mp.nstr(f, 6), mp.nstr(f - MAX_DEV, 6),
                 "VACUOUS (band includes <=0)" if f - MAX_DEV <= 0 else "ok"))

    # ---- GATE D: sanity of the printed fit interval endpoints ----
    print("\n[GATE D] fit interval as printed (two instruments agree):")
    print("         0.4 <= lam <= 20   and   0 <= mu* <= 0.2")
    print("         (NOT '0 <= mu* < 0.2'; the strict '<' in the paper belongs to")
    print("          a different sentence: 'exact analytic results for 0 <= lam < 0.4')")

    print("\n[RMS] for reference, rms deviation %s is %s%% of f at lam=3.2953, mu*=0.13"
          % (mp.nstr(RMS_DEV, 3),
             mp.nstr(RMS_DEV / f_MC('3.2953', '0.13') * 100, 4)))
    print("\nALL DONE. Grade of everything printed above: numerical-evidence.")


if __name__ == "__main__":
    main()
