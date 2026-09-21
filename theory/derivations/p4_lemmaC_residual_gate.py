"""
P4 - Lemma C: gate hardening the SINGLE residual inequality and three negative
closure attempts, on branch deep-lemmaC-analytic.

Builds its own O_of/f_lor/top_eig (does not import eigh-trust from any other
branch's gate) and imports ONLY the read-only, already-merged
theory/derivations/p4_concentration_bound.py for u_vec/O_of/f_lor/top_eig reuse
(these are elementary, independently re-derivable one-liners; no eigenvalue
CLAIM is imported, only the operator-construction utilities).

Run: OMP_NUM_THREADS=1 python theory/derivations/p4_lemmaC_residual_gate.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p4_concentration_bound import u_vec, O_of, f_lor, top_eig  # noqa: E402

PASS = 0
FAIL = 0


def ck(cond, name, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
    else:
        FAIL += 1
        print(f"  FAIL: {name}  {detail}")


OMEGAS = [2, 3, 5, 8, 12, 20, 30, 40]
N = 1200


# ---------------------------------------------------------------------------
# Part A. The single residual, in the form fixed by P4_lemmaC_var_bound.md
# (commit c2af5d2, worktree agent-ae6c83eb3d884027f): the deficit
# h-q = <v, O[f-f^2] v> and h-m3 = <v, O[f-f^3] v> must be captured to >=0.9
# by any analytic lower bound, or (SUFF)/(QF) do not close.  Reproduce the
# cone-coercivity capture ratio kappa_L3[phi,v]/<v,O[phi]v> for phi=f-f^2
# INDEPENDENTLY (own O_of/top_eig build, not the g2 branch's numbers).
# ---------------------------------------------------------------------------
def kappa_L3(phi_vals_at_S, w):
    # phi_vals_at_S: precomputed phi(n+m+1) matrix; w = u*v
    return float((phi_vals_at_S * np.outer(w, w)).sum())


print("=" * 100)
print("Part A: independent reproduction of the deficit capture ratio (f-f^2, f-f^3)")
print(f"{'om':>5}{'x':>8}{'h-q':>10}{'kL3[f-f2]':>12}{'cap':>8}{'h-m3':>10}{'kL3[f-f3]':>12}{'cap':>8}")
caps_hq = []
caps_hm3 = []
for om in OMEGAS:
    x = om * om
    A = O_of(f_lor(x), N)
    h, v = top_eig(A)
    f = f_lor(x)
    B = O_of(lambda k: f(k) ** 2, N)
    F3 = O_of(lambda k: f(k) ** 3, N)
    q = float(v @ (B @ v))
    m3 = float(v @ (F3 @ v))
    hq = h - q
    hm3 = h - m3
    u = u_vec(N)
    w = u * v
    idx = np.arange(N)
    S = idx[:, None] + idx[None, :] + 1
    phi_m = f(S) - f(S) ** 2
    phi_m3 = f(S) - f(S) ** 3
    kL_hq = kappa_L3(phi_m, w)
    kL_hm3 = kappa_L3(phi_m3, w)
    cap_hq = kL_hq / hq
    cap_hm3 = kL_hm3 / hm3
    caps_hq.append(cap_hq)
    caps_hm3.append(cap_hm3)
    print(f"{om:5.1f}{x:8.1f}{hq:10.4f}{kL_hq:12.4f}{cap_hq:8.3f}{hm3:10.4f}{kL_hm3:12.4f}{cap_hm3:8.3f}")

# The var_bound note reports cap(h-q), cap(h-m3) in [0.58,0.73]-ish (decreasing with x
# for the g_E kernel; here for phi_M=f-f^2/f-f^3 the qualitative claim is the same:
# capture is bounded well below 0.9 for all tested x, so cone-coercivity alone is
# insufficient -- this is the load-bearing claim of Part A.
ck(all(c < 0.90 for c in caps_hq), "cap(h-q) < 0.90 for all tested omega (cone-coercivity insufficient)",
   detail=str([round(c, 3) for c in caps_hq]))
ck(all(c < 0.90 for c in caps_hm3), "cap(h-m3) < 0.90 for all tested omega (cone-coercivity insufficient)",
   detail=str([round(c, 3) for c in caps_hm3]))
ck(max(caps_hq) - min(caps_hq) > 0.05, "cap(h-q) is not a flat constant (varies with x)",
   detail=f"range={min(caps_hq):.3f}..{max(caps_hq):.3f}")

# ---------------------------------------------------------------------------
# Part B. Negative closure attempt #1: restricted summation (concentration
# branch's route to kappa_prof ~ 1/x). Does restricting kappa_L3 to the
# n,m <= a*sqrt(x) box (where profile bound v_n >= c*v0*u_n holds) give a
# useful fraction of the full kappa_L3?  Check the ratio
# kappa_L3_restricted / kappa_L3_full.  From concentration §7, the E-kernel
# (g_E = 2xk^2/(k^2+x)^3) retains ~77% (a=1) / ~40% (a=1/2) of the full sum.
# Here for phi_M=f-f^2, check whether the retention is similarly useful.
# ---------------------------------------------------------------------------
print("=" * 100)
print("Part B: restricted summation (n,m <= a*sqrt(x)) retention ratio")
print(f"{'om':>5}{'x':>8}{'kL3_full':>10}{'kL3_a1':>10}{'kL3_a.5':>10}{'a1/full':>9}{'a.5/full':>10}")
retentions_a1 = []
retentions_a5 = []
for om in OMEGAS:
    x = om * om
    A = O_of(f_lor(x), N)
    h, v = top_eig(A)
    f = f_lor(x)
    u = u_vec(N)
    w = u * v
    idx = np.arange(N)
    S = idx[:, None] + idx[None, :] + 1
    phi_m = f(S) - f(S) ** 2
    kL_full = kappa_L3(phi_m, w)
    sqx = np.sqrt(x)
    n1 = int(round(sqx))
    n5 = int(round(0.5 * sqx))
    kL_a1 = kappa_L3(phi_m[:n1 + 1, :n1 + 1], w[:n1 + 1])
    kL_a5 = kappa_L3(phi_m[:n5 + 1, :n5 + 1], w[:n5 + 1])
    ret_a1 = kL_a1 / kL_full
    ret_a5 = kL_a5 / kL_full
    retentions_a1.append(ret_a1)
    retentions_a5.append(ret_a5)
    print(f"{om:5.1f}{x:8.1f}{kL_full:10.4f}{kL_a1:10.4f}{kL_a5:10.4f}{ret_a1:9.4f}{ret_a5:10.4f}")

# The concentration route works for g_E because its restricted sum retains a
# LARGE fraction (>=0.77 for a=1).  For phi_M=f-f^2, the retention is ~75-88%
# (a=1) and drops with x; a=1/2 gives only ~40-72% and also drops.  Combined
# with the already-insufficient capture ratios (~0.58-0.73), the product
# c^2 * retention * capture is even worse.  This route does not close (★).
ck(all(0.7 < r < 0.9 for r in retentions_a1), "retention(a=1) in [0.7,0.9] range (drops with x)",
   detail=str([round(r, 3) for r in retentions_a1]))
ck(max(retentions_a1) - retentions_a1[-1] > 0.05, "retention(a=1) decreasing with x (not improving)",
   detail=f"om=2: {retentions_a1[0]:.3f}, om=40: {retentions_a1[-1]:.3f}")

# ---------------------------------------------------------------------------
# Part C. Negative closure attempt #2: kappa_L3 for the UNRESTRICTED phi=f
# kernel (pointwise dominates phi_M=f-f^2, so kappa_L3[f] >= kappa_L3[f-f^2]
# automatically).  Does kappa_L3[f] come close enough to h-q to be useful?
# Measure the ratio kappa_L3[f] / (h-q).  A ratio near 1 would give a useful
# analytic lower bound on h-q; a ratio >>1 means kappa_L3[f] captures a
# DIFFERENT quantity (not h-q) and is not helpful for the deficit.
# ---------------------------------------------------------------------------
print("=" * 100)
print("Part C: kappa_L3[f] vs h-q (f pointwise dominates f-f^2)")
print(f"{'om':>5}{'x':>8}{'h-q':>10}{'kL3[f]':>10}{'kL3[f]/hq':>12}")
ratios_f = []
for om in OMEGAS:
    x = om * om
    A = O_of(f_lor(x), N)
    h, v = top_eig(A)
    f = f_lor(x)
    B = O_of(lambda k: f(k) ** 2, N)
    q = float(v @ (B @ v))
    hq = h - q
    u = u_vec(N)
    w = u * v
    idx = np.arange(N)
    S = idx[:, None] + idx[None, :] + 1
    phi_f = f(S)
    kL_f = kappa_L3(phi_f, w)
    ratio_f = kL_f / hq
    ratios_f.append(ratio_f)
    print(f"{om:5.1f}{x:8.1f}{hq:10.4f}{kL_f:10.4f}{ratio_f:12.3f}")

# kappa_L3[f] / (h-q) is ~2-4, meaning kappa_L3[f] captures a quantity
# roughly DOUBLE h-q, not h-q itself.  This is because kappa_L3[f] (by the
# cone-coercivity rung identity, P4_cone_coercivity.md §1) equals the L3
# part of <v,O[f]v>=h, which includes BOTH the deficit h-q AND the q term's
# L3 contribution.  So kappa_L3[f] is not a lower bound on h-q alone; it
# lower-bounds a sum that contains h-q plus other positive terms.  This route
# does not give a useful h-q lower bound.
ck(all(1.5 < r < 6.0 for r in ratios_f), "kL3[f]/(h-q) in [1.5,6] (captures sum, not deficit alone)",
   detail=str([round(r, 3) for r in ratios_f]))
ck(ratios_f[-1] > ratios_f[0] + 1.0, "kL3[f]/(h-q) grows with x (does not stabilize near 1)",
   detail=f"om=2: {ratios_f[0]:.3f}, om=40: {ratios_f[-1]:.3f}")

# ---------------------------------------------------------------------------
# Part D. Negative closure attempt #3: the L1+L2 rung remainder.  The full
# rung identity (P4_cone_coercivity Lemma 1) is <v,O[phi]v> = kappa_L3[phi,v]
# + L12[phi,v], where L12 := sum_k phi(k) * 2 * sum_{n>=k} w_n(w_{n-k}-w_n).
# Since w is nonincreasing, w_{n-k}-w_n <= 0, so L12 <= 0 (cone-coercivity
# discards it to get a lower bound).  Can we lower-bound L12 using the SAME
# profile bound w_n >= c*v0*u_n^2 that gives kappa_L3?  Naive attempt: bound
# w_n from below by c*v0*u_n^2 and w_{n-k} from above by v0*u_{n-k} (since
# v_{n-k}<=v_0), giving w_{n-k}-w_n >= v0*u_{n-k} - c*v0*u_n^2.  Check
# whether this naive bound is ever positive (if always negative, L12's naive
# bound goes the wrong way for a lower bound on <v,O[phi]v>).
# ---------------------------------------------------------------------------
print("=" * 100)
print("Part D: L1+L2 rung remainder — naive profile-bound attempt")
print(f"{'om':>5}{'x':>8}{'frac_pos':>10}  (fraction of n<=sqrt(x) where naive bound > 0)")
for om in OMEGAS:
    x = om * om
    A = O_of(f_lor(x), N)
    h, v = top_eig(A)
    u = u_vec(N)
    v0 = v[0]
    c = 0.20  # concentration c_tail(a=1/2, x>=25) >= 0.19
    n_arr = np.arange(1, N)
    # naive lower bound on w_{n-1} - w_n (k=1 case): c*v0*u_{n-1}^2 - v0*u_n
    # (need w_{n-1} LOWER and w_n UPPER, so roles are c*u^2 for n-1, plain v0*u for n)
    naive_lb = c * v0 * u[np.maximum(n_arr - 1, 0)] ** 2 - v0 * u[n_arr]
    sqx_int = int(np.sqrt(x))
    frac_pos = np.mean(naive_lb[:sqx_int] > 0) if sqx_int > 0 else 0.0
    print(f"{om:5.1f}{x:8.1f}{frac_pos:10.3f}")

# The naive bound is NEVER positive in the relevant range n<=sqrt(x), so the
# L1+L2 remainder cannot be salvaged this way.  The asymptotic reason: u_n ~
# n^{-1/2} decays slower than u_n^2 ~ n^{-1}, so for large n, u_n > c*u_{n-1}^2
# always, making the bound negative.  Recovering L1+L2 would need a profile
# UPPER bound (not just the LOWER bound in hand), which is a genuinely new
# ingredient not provided by the concentration branch.
ck(all(0.0 == fp for fp in [np.mean((c * v0 * u[np.maximum(np.arange(1, N) - 1, 0)] ** 2
                                      - v0 * u[np.arange(1, N)])[:int(np.sqrt(om * om))] > 0)
                             for om in OMEGAS]), "naive L12 bound always <=0 in n<=sqrt(x) (vacuous)",
   detail="(all fractions are 0.000)")

print("=" * 100)
print(f"ck PASS = {PASS}    FAIL = {FAIL}")
if FAIL > 0:
    print("   *** GATE FAILED ***")
    sys.exit(1)
print("   All checks passed.")

