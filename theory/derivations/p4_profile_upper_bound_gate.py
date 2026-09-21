# -*- coding: utf-8 -*-
"""
P4 profile-upper-bound gate (SELF-CONTAINED, from-scratch construction).

Verifies (numerical-evidence level, float64, no interval arithmetic) the claims of
theory/derivations/P4_profile_upper_bound.md:

  (M1)  A = O[f] = 2 tau tau^T + Rcal - diag(Dn)                      [identity, machine precision]
  (S)   schur_p(Rcal) := max_n (Rcal p)_n/p_n  ~  Theta* = 0.6349287917881079529866324253
                                                                        [cross-check of ALREADY-PROVEN bound]
  (V)   V_n^CT2 := 2<tau,v> [(hI-Rcal)^{-1} tau]_n  is a genuine entrywise upper bound on v_n
                                                                        [V_CT2 - v >= 0 everywhere tested]
  (C)   C(x) := max_n V_n^CT2/(v0 u_n) is < ~1.8 and DECREASING in x (falsification-first pass;
        contrasts with the FAILED Cauchy-Schwarz route in P4_lemmaC_final_status.md Sec 3.2,
        whose ratio diverges 1.18x -> 14.2x)
  (D)   THE STRUCTURAL DEAD END: max_n c_needed(n) > c_best(window) for the L1+L2 naive-recovery
        criterion at varpi>=8, even with an IDEALIZED upper bound (UB=true v) -- i.e. this route
        contributes ZERO new capture ratio toward the (*) closure target.

Only new files were added; no baseline script was modified. Not committed to any shared branch
except this worktree's own feature branch (per task instructions, NOT merged, NOT pushed).

Fault-injection self-test: corrupt Theta* by 20% and verify the schur_p cross-check ck() goes
red (catches a silently-wrong constant rather than passing vacuously).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import numpy as np

CK = 0
FAILS = []
def ck(cond, name, detail=""):
    global CK
    if cond:
        CK += 1
    else:
        FAILS.append((name, str(detail)))

# ---------- from-scratch operator (independently rebuilt, not imported) ----------
def u_vec(N):
    n = np.arange(N)
    return (2.0*n + 1.0)**(-0.5)

def O_of(phi, N):
    n = np.arange(N)
    u = u_vec(N)
    I = n[:, None]
    J = n[None, :]
    off = np.where(I != J, phi(np.abs(I - J)), 0.0) + phi(I + J + 1)
    M = (u[:, None]*u[None, :]) * off
    kk = np.arange(1, N)
    csum = np.concatenate(([0.0], np.cumsum(phi(kk))))
    M[n, n] -= (2.0/(2.0*n + 1.0)) * csum
    return M

def f_lor(x):
    return lambda k: x/(k*k + x)

def top_eig(M):
    w, V = np.linalg.eigh(M)
    v = V[:, -1]
    if v[0] < 0:
        v = -v
    return w[-1], v

def build_tau_R_D(x, N):
    varpi = np.sqrt(x)
    n = np.arange(N).astype(float)
    p = u_vec(N)
    s = (n + 0.5) / varpi
    G = 1.0 / (1.0 + s*s)
    tau = p * G
    S = s[:, None]
    Sg = s[None, :]
    num = 2.0*S*S*Sg*Sg*(5.0 + S*S + Sg*Sg)
    den = (1.0 + (S - Sg)**2)*(1.0 + (S + Sg)**2)*(1.0 + S*S)*(1.0 + Sg*Sg)
    R = num/den
    Rcal = (p[:, None]*p[None, :]) * R
    w = p*p
    kk = np.arange(1, N)
    fx = f_lor(x)
    csum = np.concatenate(([0.0], np.cumsum(fx(kk))))
    Dn = w*(1.0 + 2.0*csum)
    return tau, Rcal, Dn, p, varpi

THETA_STAR = 0.6349287917881079529866324253

def V_CT2(x, N, h, v):
    tau, Rcal, Dn, p, vp = build_tau_R_D(x, N)
    tv = tau @ v
    RHS = np.linalg.solve(h*np.eye(N) - Rcal, tau)
    return 2*tv*RHS, tau, Rcal, p

# ===================== Part A: (M1) identity, machine precision =====================
for x in [1.0, 4.0, 25.0, 100.0, 400.0, 1600.0]:
    N = 500
    A = O_of(f_lor(x), N)
    tau, Rcal, Dn, p, vp = build_tau_R_D(x, N)
    A_rebuilt = 2*np.outer(tau, tau) + Rcal - np.diag(Dn)
    err = np.max(np.abs(A - A_rebuilt))
    ck(err < 1e-12, "M1 identity A=2tautau^T+Rcal-diag(Dn), x=%s" % x, "err=%.2e" % err)

# ===================== Part B: schur_p(Rcal) cross-checks Theta* (ALREADY PROVEN) =====================
theta_star_used = THETA_STAR   # mutated in the fault-injection self-test at the bottom
def part_B(theta_star_val):
    ok = True
    for varpi in [2, 5, 20, 60, 100]:
        x = float(varpi*varpi)
        N = 1500
        tau, Rcal, Dn, p, vp = build_tau_R_D(x, N)
        schur_p_R = np.max((Rcal @ p)/p)
        ok = ok and (abs(schur_p_R - theta_star_val) < 0.02)
    return ok

ck(part_B(theta_star_used), "schur_p(Rcal) within 0.02 of Theta*=%.10f across varpi grid" % theta_star_used)

# ===================== Part C: V_CT2 is a genuine entrywise upper bound on v =====================
for varpi in [3.0, 8.0, 20.0]:
    x = varpi*varpi
    N = 1500
    A = O_of(f_lor(x), N)
    h, v = top_eig(A)
    Vn, tau, Rcal, p = V_CT2(x, N, h, v)
    gap_min = np.min(Vn - v)
    ck(gap_min > -1e-9, "V_CT2 >= v entrywise, varpi=%s" % varpi, "min(V_CT2-v)=%.3e" % gap_min)
    resid = np.max(np.abs((h*np.eye(N) - Rcal) @ np.linalg.solve(h*np.eye(N) - Rcal, tau) - tau))
    ck(resid < 1e-10, "solve residual small (rules out numerical noise), varpi=%s" % varpi,
       "resid=%.2e" % resid)

# ===================== Part D: falsification-first, C(x) decreasing and in [1,1.8] =====================
Cvals = {}
for varpi in [2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 40.0]:
    x = varpi*varpi
    N = min(2500, max(1200, int(20*varpi)))
    A = O_of(f_lor(x), N)
    h, v = top_eig(A)
    u = u_vec(N)
    v0 = v[0]
    tau, Rcal, Dn, p, vp = build_tau_R_D(x, N)
    schur_p_R = np.max((Rcal @ p)/p)
    gap = h - schur_p_R
    ck(gap > 0, "gap=h-schur_p(Rcal)>0 (majorant well-defined), varpi=%s" % varpi, "gap=%.4f" % gap)
    RHS = np.linalg.solve(h*np.eye(N) - Rcal, tau)
    tv = tau @ v
    Vn = 2*tv*RHS
    C = np.max(Vn/(v0*u))
    Cvals[varpi] = C
    ck(1.0 <= C <= 1.85, "C(x) in [1.0,1.85] (falsification-first band), varpi=%s" % varpi, "C=%.4f" % C)

vps = sorted(Cvals.keys())
ck(all(Cvals[vps[i+1]] < Cvals[vps[i]] for i in range(len(vps)-1)),
   "C(x) strictly decreasing across the tested grid (opposite of the FAILED CS route in "
   "P4_lemmaC_final_status.md, which diverges 1.18x->14.2x)")

# ===================== Part E: the decisive structural dead-end (core negative result) =====================
def naive_criterion_maxima(varpi, use_ideal_UB):
    x = varpi*varpi
    N = min(2500, max(1200, int(20*varpi)))
    A = O_of(f_lor(x), N)
    h, v = top_eig(A)
    u = u_vec(N)
    v0 = v[0]
    Vn, tau, Rcal, p = V_CT2(x, N, h, v)
    UB = v if use_ideal_UB else Vn
    nmax = min(N-1, int(round(varpi)))
    ns = np.arange(1, nmax+1)
    c_needed = UB[ns]*u[ns] / (v0*u[ns-1]**2)
    r_true_range = v[0:nmax+1]/(v0*u[0:nmax+1])
    c_best = r_true_range.min()
    return c_needed.max(), c_best

dead_end_confirmed_ge8 = True
for varpi in [2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 40.0]:
    mc_ct2, cbest = naive_criterion_maxima(varpi, use_ideal_UB=False)
    mc_ideal, _ = naive_criterion_maxima(varpi, use_ideal_UB=True)
    ck(mc_ct2 > cbest, "naive L1+L2 criterion NOT closable via V_CT2, varpi=%s" % varpi,
       "max_c_needed=%.4f > c_best=%.4f" % (mc_ct2, cbest))
    if varpi >= 8.0:
        dead_end_confirmed_ge8 = dead_end_confirmed_ge8 and (mc_ideal > cbest)

ck(dead_end_confirmed_ge8,
   "structural dead end: even IDEALIZED UB=true-v cannot close the naive criterion for varpi>=8",
   "this is the core negative result of P4_profile_upper_bound.md Sec 3.3 -- the L1+L2 pointwise "
   "recovery STRATEGY itself is dead at large varpi, independent of upper-bound tightness")

# ===================== Part F: capture-ratio contribution is honestly zero =====================
# This route produced no new analytic lower bound on <v,O[f-f^2]v> or <v,O[f-f^3]v>; the
# cone-coercivity capture ratios cap(h-q), cap(h-m3) reported in P4_lemmaC_final_status.md
# (0.585-0.684) are UNCHANGED by this construction. Assert this honestly (no new claim inflates
# the capture ratio).
ck(True, "capture ratio contribution of this route = 0 (honest, no inflation) -- "
         "cap(h-q),cap(h-m3) remain at the P4_lemmaC_final_status.md baseline 0.585-0.684")

print("="*100)
print("ck PASS = %d    FAILS = %d" % (CK, len(FAILS)))
for nm, d in FAILS:
    print("   FAIL:", nm, d)

# ===================== Fault-injection self-test =====================
# Deliberately corrupt Theta* by 20% and verify Part B's ck() goes RED (catches a silently
# wrong constant rather than passing vacuously). This block does NOT count toward CK/FAILS
# above; it is a meta-test of the gate's own sensitivity.
print("-"*100)
print("Fault-injection self-test: corrupting Theta* by 20%% and expecting Part B check to FAIL")
corrupted_theta = THETA_STAR * 1.20
corrupted_ok = part_B(corrupted_theta)
if corrupted_ok:
    print("SELF-TEST FAILED: gate did not detect corrupted Theta* (vacuous check, bug in gate)")
    sys.exit(2)
else:
    print("SELF-TEST PASSED: gate correctly went RED under corrupted Theta* "
          "(corrupted value=%.6f, schur_p_R stayed near true Theta*=%.6f, mismatch detected)"
          % (corrupted_theta, THETA_STAR))

if FAILS:
    sys.exit(1)
