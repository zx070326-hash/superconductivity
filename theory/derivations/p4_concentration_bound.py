# -*- coding: utf-8 -*-
"""
P4 concentration-bound certificate probe (SELF-CONTAINED).

Verifies, from a from-scratch construction of the kernel operator O[.]:
  (S1)  O[1] = 2 u u^T - I
  (star)   (h+1) v = 2 s u - O[gtil] v           (per-component resolvent identity)
  (star2)  v = 2 s ((h+1) I + O[gtil])^{-1} u    (resolvent representation)
  a priori bounds  (P)  w_n <= (2 s / h) u_n^2
                   (P') w_n <= (2 s / h^2) u_n^2 a_n
                   (T)  Sum_{m>=M} w_m <= 3 Ftot s / (h (2M+1))
  Ftot closed bound  Ftot <= 0.502 * pi * sqrt(x)   (x >= 1)
  Rayleigh lower bound  h >= h_lo := <psi,A psi>/||psi||^2 ,  psi = u 1_{[0,L)}
  s-independence of the epsilon bound
  the MAIN rigorous profile bound  v_n >= c(a,x) v_0 u_n  on n <= a sqrt(x),
      with c computed from h_lo alone (no knowledge of v needed).

Notation:  u_n = (2n+1)^{-1/2};  f(k;x)=x/(k^2+x);  gtil(k;x)=k^2/(k^2+x)=1-f;
           A=O[f], h=lambda1(A), v=top unit eigenvector (v0>0), s=<u,v>, w=u*v.
"""
import numpy as np

CK = 0; FAILS = []
def ck(cond, name, detail=""):
    global CK
    if cond: CK += 1
    else: FAILS.append((name, str(detail)))

# ---------- from-scratch operator ----------
def u_vec(N):
    n = np.arange(N); return (2.0*n+1.0)**(-0.5)

def O_of(phi, N):
    n = np.arange(N); u = u_vec(N)
    I = n[:,None]; J = n[None,:]
    off = np.where(I != J, phi(np.abs(I-J)), 0.0) + phi(I+J+1)
    M = (u[:,None]*u[None,:]) * off
    kk = np.arange(1, N)
    csum = np.concatenate(([0.0], np.cumsum(phi(kk))))   # csum[n]=sum_{k=1}^n phi(k)
    M[n, n] -= (2.0/(2.0*n+1.0)) * csum
    return M

def f_lor(x):  return lambda k: x/(k*k + x)
def gtil(x):   return lambda k: (k*k)/(k*k + x)
def one_k(k):  return np.ones_like(np.asarray(k, dtype=float))

def top_eig(M):
    w, V = np.linalg.eigh(M); v = V[:,-1]
    if v[0] < 0: v = -v
    return w[-1], v

# ---------- (S1) ----------
N0 = 600; u0 = u_vec(N0)
S1_err = np.max(np.abs(O_of(one_k, N0) - (2*np.outer(u0,u0) - np.eye(N0))))
ck(S1_err < 1e-12, "S1: O[1]=2uu^T-I", "err=%.2e" % S1_err)

# ---------- identities and A = O[1]-O[gtil] ----------
for x in [0.01,0.1,1.0,4.0,9.0,25.0,100.0,400.0,1600.0]:
    N = 900; u = u_vec(N)
    A  = O_of(f_lor(x), N); Og = O_of(gtil(x), N); O1 = O_of(one_k, N)
    ck(np.max(np.abs(A-(O1-Og))) < 1e-12, "A=O[1]-O[gtil] x=%s" % x)
    h, v = top_eig(A); s = u @ v
    ck(np.max(np.abs((h+1)*v - (2*s*u - Og@v))) < 1e-10, "(star) x=%s" % x)
    vv = 2*s*np.linalg.solve((h+1)*np.eye(N)+Og, u)
    ck(np.max(np.abs(vv-v)) < 1e-9, "(star2) x=%s" % x)
    ck(np.max(np.abs(Og@v - (2*s*u-(h+1)*v))) < 1e-10, "corr-id x=%s" % x)
    ck(v[0] <= 2*s/(h+1) + 1e-12, "v0<=2s/(h+1) x=%s" % x)

# ---------- Ftot closed bound ----------
for x in [1,4,25,100,400,1600]:
    Nbig = 400000; kk = np.arange(1,Nbig)
    Ftot_trunc = np.sum(x/(kk*kk+x))              # undercounts by tail <= x/(Nbig-1)
    tail_ub = 2.0*x/Nbig                          # Sum_{k>=Nbig} x/(k^2+x) < x/(Nbig-1) < 2x/Nbig
    closed = (np.pi*np.sqrt(x)/2)/np.tanh(np.pi*np.sqrt(x)) - 0.5   # EXACT value of the series
    # closed form matches the truncated sum to within the analytic truncation tail:
    ck(Ftot_trunc <= closed + 1e-9 and closed <= Ftot_trunc + tail_ub + 1e-9,
       "Ftot closed form x=%s" % x, "%s vs %s (tail<=%.1e)" % (Ftot_trunc,closed,tail_ub))
    ck(closed <= 0.502*np.pi*np.sqrt(x), "Ftot<=0.502 pi sqx x=%s" % x)

# ---------- a priori bounds ----------
def build_kernels(x, N):
    u = u_vec(N); n = np.arange(N).astype(float); u2 = u**2
    I = n[:,None]; J = n[None,:]
    f = f_lor(x); gt = gtil(x)
    Foff = np.where(I!=J, f(np.abs(I-J)), 0.0); Fl3 = f(I+J+1)
    Goff = np.where(I!=J, gt(np.abs(I-J)), 0.0); Gl3 = gt(I+J+1)
    return u, u2, Foff, Fl3, Goff, Gl3

for x in [4,25,100,400,1600]:
    N = 1600; u, u2, Foff, Fl3, Goff, Gl3 = build_kernels(x, N)
    A = O_of(f_lor(x), N); h, v = top_eig(A); s = u@v; w = u*v
    kk = np.arange(1,N); Ftot = np.sum(x/(kk*kk+x))
    ck(np.max(w - (2*s/h)*u2) <= 1e-6, "(P) w<=(2s/h)u^2 x=%s" % x, np.max(w-(2*s/h)*u2))
    a = Foff@u2 + Fl3@u2
    ck(np.max(w - (2*s/h**2)*u2*a) <= 1e-6, "(P') w<=(2s/h^2)u^2 a x=%s" % x, np.max(w-(2*s/h**2)*u2*a))
    Ttrue = np.array([np.sum(w[M:]) for M in range(1,N,50)])
    Tbnd  = np.array([3*Ftot*s/(h*(2*M+1)) for M in range(1,N,50)])
    ck(np.all(Ttrue <= Tbnd + 1e-9), "(T) tail-mass x=%s" % x)

# ---------- Rayleigh h_lo, s-independence, MAIN bound ----------
def h_lower(x, N):
    A = O_of(f_lor(x), N); u = u_vec(N); best = -1e9
    for L in [2,4,8,16,32,64,128,256,512]:
        if L > N: break
        psi = u.copy(); psi[L:] = 0.0
        best = max(best, (psi@(A@psi))/(psi@psi))
    return best

def eps_bound_sfree(x, N, hval):
    u, u2, Foff, Fl3, Goff, Gl3 = build_kernels(x, N)
    a = Foff@u2 + Fl3@u2
    coeff = u2 * np.minimum(1.0, a/hval)          # = B(m)/(2 s / hval); s cancels in eps
    return (Goff@coeff + Gl3@coeff)/hval

def c_of(x, a_range, N, hval):
    # TRUNCATED c: 1 - max Psi_n^{trunc,N}.  Since a truncated Psi UNDERESTIMATES the true
    # (infinite) Psi (all terms >= 0), this is an UPPER bound on the true valid c -- NOT a
    # rigorous lower bound.  Use c_tail(...) below for a genuine lower bound.
    eb = eps_bound_sfree(x, N, hval); sqx = np.sqrt(x)
    nt = int(round(a_range*sqx))
    if nt < 1: return float('nan')
    return 1.0 - eb[:nt+1].max()

# ---------- rigorous m-tail of Psi (decay lemma) => tail-inclusive LOWER bound on c ----------
# Decay lemma (analytic, see P4_concentration_bound.md 3.1):  for all m>=1,
#   a_m <= a_ub(x,m) := f(ceil(m/2))*S(floor(m/2)) + 2 Ftot/(m+1) + f(m+1)*S(m) + Ftot/(2m+3),
#   S(K):=sum_{j=0}^K 1/(2j+1) <= 1 + 0.5 ln(2K+1);   Ftot <= 0.502 pi sqrt(x)   [ (F) ].
# Then, since 0<=g̃<1 (bracket<=2) and min(1,a_m/h)<=a_m/h,
#   Psi_n^{tail,>=N} <= (2/h^2) sum_{m>=N} u_m^2 a_ub(x,m) =: tail_bnd(x,N,h)   (n-independent).
# Hence c_tail := 1 - (max_n Psi_n^{trunc,N} + tail_bnd) is a rigorous LOWER bound on c.
def Ftot_ub(x): return 0.502*np.pi*np.sqrt(x)
def a_ub(x, m):
    m = np.asarray(m, float); F = Ftot_ub(x); f = f_lor(x)
    S = lambda K: 1.0 + 0.5*np.log(2.0*K + 1.0)
    return (f(np.ceil(m/2.0))*S(np.floor(m/2.0)) + 2.0*F/(m+1.0)
            + f(m+1.0)*S(m) + F/(2.0*m+3.0))
def tail_bnd(x, N, hlo, Mbig=200000):
    m = np.arange(N, Mbig+1, dtype=float)
    part = np.sum((1.0/(2.0*m+1.0))*a_ub(x, m))
    F = Ftot_ub(x); Mb = float(Mbig)
    R = 1.5*F/Mb + 1.5*x*(1.0+0.5*np.log(2.0*Mb+1.0))/(Mb*Mb)   # analytic remainder m>Mbig
    return (2.0/hlo**2)*(part + R)
def c_tail(x, a_range, N, hval):
    eb = eps_bound_sfree(x, N, hval); sqx = np.sqrt(x)
    nt = int(round(a_range*sqx))
    if nt < 1: return float('nan')
    return 1.0 - (eb[:nt+1].max() + tail_bnd(x, N, hval))

for x in [25,100,400]:
    N = 1500; u, u2, Foff, Fl3, Goff, Gl3 = build_kernels(x, N)
    A = O_of(f_lor(x), N); h, v = top_eig(A); s = u@v
    a = Foff@u2 + Fl3@u2
    B = np.minimum((2*s/h)*u2, (2*s/h**2)*u2*a)
    eb_s = (Goff@B + Gl3@B)/(2*s)
    eb_free = eps_bound_sfree(x, N, h)
    ck(np.max(np.abs(eb_s-eb_free)) < 1e-12, "s-independence x=%s" % x)

for x in [16,100,400,1600]:
    N = 1600; A = O_of(f_lor(x), N); h, _ = top_eig(A); hlo = h_lower(x, N)
    ck(hlo <= h + 1e-9, "h_lo<=h x=%s" % x)
    ck(hlo >= 0.9*h,   "h_lo>=0.9h x=%s" % x, "hlo/h=%.3f" % (hlo/h))

for x in [25,100,400,1600]:
    N = 1600; u, u2, Foff, Fl3, Goff, Gl3 = build_kernels(x, N)
    A = O_of(f_lor(x), N); h, v = top_eig(A); s = u@v; hlo = h_lower(x, N)
    eps_true = (O_of(gtil(x), N)@v)/(2*s*u)
    eb = eps_bound_sfree(x, N, hlo)
    idx = int(2*np.sqrt(x))
    ck(np.all(eps_true[:idx+1] <= eb[:idx+1] + 1e-9), "eps_true<=eps_bound(h_lo) x=%s" % x)

# ---------- decay lemma check: a_m <= a_ub(x,m) (validates 3.1) ----------
for x in [4.0, 25.0, 100.0, 400.0]:
    N = 2000; u2 = u_vec(N)**2
    _, _, Foff, Fl3, _, _ = build_kernels(x, N)
    a_exact = Foff@u2 + Fl3@u2
    aub = a_ub(x, np.arange(1, N))
    ck(np.all(a_exact[1:] <= aub + 1e-12), "decay lemma a_m<=a_ub x=%s" % x,
       "max ratio %.3f" % np.max(a_exact[1:]/aub))

# ---------- memoized c-pack (one kernel build per x): (hlo, cT.5, ci.5, cT1, ci1) ----------
_cpack = {}
def c_pack(x, N=2000):
    key = (round(x, 6), N)
    if key in _cpack: return _cpack[key]
    u, u2, Foff, Fl3, Goff, Gl3 = build_kernels(x, N)
    A = O_of(f_lor(x), N); best = -1e9
    for L in [2,4,8,16,32,64,128,256,512,1024]:
        if L > N: break
        psi = u.copy(); psi[L:] = 0.0; best = max(best, (psi@(A@psi))/(psi@psi))
    hlo = best
    a = Foff@u2 + Fl3@u2; coeff = u2*np.minimum(1.0, a/hlo)
    Psi = (Goff@coeff + Gl3@coeff)/hlo
    tb = tail_bnd(x, N, hlo); sqx = np.sqrt(x)
    def cc(ar, inc):
        nt = int(round(ar*sqx))
        if nt < 1: return float('nan')
        return 1.0 - (Psi[:nt+1].max() + (tb if inc else 0.0))
    val = (hlo, cc(0.5,False), cc(0.5,True), cc(1.0,False), cc(1.0,True))
    _cpack[key] = val; return val
def cti(x, a): return c_pack(x)[2] if a == 0.5 else c_pack(x)[4]   # tail-incl c
def ctr(x, a): return c_pack(x)[1] if a == 0.5 else c_pack(x)[3]   # truncated c

# ---------- N-consistency of the tail-inclusive c (tol ~0.006, was 0.05) ----------
for x in [64.0, 400.0, 1600.0]:
    c12 = c_pack(x, 1200)[2]; c20 = c_pack(x, 2000)[2]
    ck(abs(c12 - c20) < 0.006, "N-consistency c_tail(.5) x=%s" % x, "|d|=%.4f" % abs(c12-c20))

# ---------- tail-MAGNITUDE hardening: tail_bnd dominates the TRUE Psi-tail (EXACT a_m) ----------
# The decay-lemma check above validates a_m <= a_ub only for m<N (the truncated region); the
# SUBTRACTED tail (m>=N) is protected here.  Both guards use the EXACT a_m -- computed by direct
# 1D summation over j, INDEPENDENT of the a_ub closed form -- so a regression that corrupts a_ub
# or tail_bnd in the tail cannot pass unseen:
#   (guard-B, pointwise) a_ub(x,m) >= a_m for sampled m in [N,4N]  (decay lemma, TAIL region).
#       Sharp: fires the instant a_ub ceases to dominate a_m (the ONLY soundness-breaking direction).
#   (guard-A, summed)    tail_bnd(x,N,hlo) >= (2/hlo^2) * sum_{m=N}^{4N} u_m^2 a_m .
#       Protects the assembled magnitude (sum + remainder + prefactor), not just the per-term formula.
# The j-sum is truncated at Jmax (all terms >=0 => UNDERESTIMATES a_m), so the RHS is a rigorous
# LOWER proxy for the true infinite tail; both inequalities are therefore SOUND domination checks.
# (Note: a mere shrink of a_ub that STAYS above a_m is not a soundness bug and correctly stays green;
#  the committed a_ub is ~2.6x loose, so guard-B fires only once a_ub drops below the true a_m.)
def a_exact_m(x, mvals, Jmax):
    f = f_lor(x); jj = np.arange(Jmax, dtype=float); u2j = 1.0/(2.0*jj+1.0)
    out = np.empty(len(mvals))
    for i, m in enumerate(mvals):
        off = f(np.abs(m - jj)).copy()
        if m < Jmax: off[int(m)] = 0.0                 # drop the j=m self term
        out[i] = off @ u2j + f(m + jj + 1.0) @ u2j
    return out

for x in [25.0, 400.0]:
    N = 2000; Mwin = 4*N; hlo = c_pack(x, N)[0]        # exact hlo used by tail_bnd in c_pack
    ms = np.arange(N, Mwin+1, 100)
    if ms[-1] != Mwin: ms = np.append(ms, Mwin)
    ae = a_exact_m(x, ms.astype(float), Mwin + 20000)  # generous j-truncation (underestimate)
    aub_s = a_ub(x, ms.astype(float))
    # (guard-B) exact a_m <= a_ub in the tail region
    ck(np.all(ae <= aub_s + 1e-15), "tail decay-lemma a_m<=a_ub (m in [N,4N]) x=%s" % x,
       "min a_ub/a_m=%.3f" % np.min(aub_s/ae))
    # (guard-A) right-Riemann LOWER bound on sum_{m=N}^{4N} u_m^2 a_m: g(m):=u_m^2 a_m is
    # non-increasing (u_m^2 and a_m both positive non-increasing), so on each block (e_{i-1},e_i]
    # the sum >= width * g(e_i) (right-edge value).  Assert the sampled non-increase as a witness.
    g = (1.0/(2.0*ms + 1.0)) * ae
    ck(np.all(np.diff(g) < 0), "g=u^2 a_m non-increasing on [N,4N] x=%s" % x)   # right-Riemann witness
    rr = (2.0/hlo**2) * np.sum(np.diff(ms) * g[1:])
    tb = tail_bnd(x, N, hlo)
    ck(tb >= rr, "tail_bnd dominates true tail (exact a_m) x=%s" % x, "tail_bnd/rr=%.2f" % (tb/rr))

# ---------- MAIN RESULT: TAIL-INCLUSIVE rigorous lower bound c_tail(a,x) ----------
# c_trunc is the truncated UPPER bound; c_tail = c_trunc - tail_bnd is the rigorous LOWER bound.
print("="*100)
print("MAIN RESULT: rigorous profile constant  v_n >= c(a,x) v_0 u_n  on n<=a sqrt(x)")
print("  c_tail = TAIL-INCLUSIVE strict lower bound (via decay-lemma Psi tail);  "
      "c_trunc = truncated UPPER bound")
print("%6s%9s%7s | %9s%9s | %9s%9s" %
      ("omega","x","h_lo","cT(.5)","c_tail.5","cT(1)","c_tail1"))
for om in [1,2,3,5,8,12,20,30,40]:
    x = float(om*om); hlo, ct5, ci5, ct1, ci1 = c_pack(x)
    ck(np.isnan(ci5) or ci5 <= ct5 + 1e-12, "tail<=trunc c(.5) om=%s" % om)  # truncation=UPPER bound
    print("%6.1f%9.2f%7.3f | %9.4f%9.4f | %9.4f%9.4f" % (om, x, hlo, ct5, ci5, ct1, ci1))

# ---------- corrected tail-inclusive constants (the only values labeled 严格下界) ----------
ck(cti(25.0,   0.5) >= 0.19, "c_tail(.5,x=25) >= 0.19",   "%.4f"%cti(25.0,0.5))
ck(cti(144.0,  0.5) >= 0.38, "c_tail(.5,x=144) >= 0.38",  "%.4f"%cti(144.0,0.5))
ck(cti(400.0,  0.5) >= 0.46, "c_tail(.5,x=400) >= 0.46",  "%.4f"%cti(400.0,0.5))
ck(cti(1600.0, 0.5) >= 0.54, "c_tail(.5,x=1600) >= 0.54", "%.4f"%cti(1600.0,0.5))
ck(cti(144.0,  1.0) >  0.0,  "c_tail(1,x=144) > 0",       "%.4f"%cti(144.0,1.0))
ck(cti(400.0,  1.0) >= 0.10, "c_tail(1,x=400) >= 0.10",   "%.4f"%cti(400.0,1.0))
ck(cti(1600.0, 1.0) >= 0.18, "c_tail(1,x=1600) >= 0.18",  "%.4f"%cti(1600.0,1.0))

# ---------- SAWTOOTH structure (c is NOT monotone in x; 4.3) ----------
# tooth-boundary subsequence x=(k/a)^2 is where index n=k enters the max -> local dip;
# these left-ends are the tooth minima and are non-decreasing in k.
ck(cti(36.0, 0.5) < cti(25.0, 0.5), "sawtooth: c_tail(.5,36) < c_tail(.5,25)")   # non-monotone
ck(all(cti(float(4*k*k), 0.5) >= 0.19 for k in range(3, 11)),
   "sawtooth inf: c_tail(.5, 4k^2) >= 0.19 for k=3..10")                          # inf over x>=25
cb = [cti(float(k*k), 1.0) for k in [11,12,13,14,20]]
ck(all(cb[i+1] > cb[i] - 1e-4 for i in range(len(cb)-1)), "a=1 boundary subseq non-decreasing")
ck(cti(100.0, 1.0) < 0.0, "a=1: c_tail(1,x=100) < 0")
ck(cti(121.0, 1.0) < 0.0, "a=1: c_tail(1,x=121) < 0 (sawtooth dip)")
ck(cti(144.0, 1.0) > 0.0, "a=1: c_tail(1,x=144) > 0")

print("="*100)
print("ck PASS = %d    FAILS = %d" % (CK, len(FAILS)))
for nm, d in FAILS: print("   FAIL:", nm, d)
