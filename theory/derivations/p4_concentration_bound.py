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
    eb = eps_bound_sfree(x, N, hval); sqx = np.sqrt(x)
    nt = int(round(a_range*sqx))
    if nt < 1: return float('nan')
    return 1.0 - eb[:nt+1].max()

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

# ---------- N-consistency + true-profile ----------
print("="*96)
print("MAIN RESULT: rigorous profile constant c(a,x)  [ v_n >= c v_0 u_n on n<=a sqrt(x) ], via h_lo")
print("%6s%9s%7s | %10s%10s | %9s%9s | %13s%12s" %
      ("omega","x","h_lo","c(.5)N400","c(.5)N700","c(1)N400","c(1)N700","true c(1)@sqx","min r_n<=sqx"))
for om in [0.5,1,2,3,5,8,12,20,30,40]:
    x = om*om; row = {}
    for N in [400,700]:
        u = u_vec(N); A = O_of(f_lor(x), N); h, v = top_eig(A); s = u@v
        hlo = h_lower(x, N)
        c05 = c_of(x, 0.5, N, hlo); c10 = c_of(x, 1.0, N, hlo)
        sqx = np.sqrt(x); nt = max(1, int(round(sqx)))
        eps_true = (O_of(gtil(x), N)@v)/(2*s*u)
        true_c1 = 1 - eps_true[:nt+1].max() if nt < N else float('nan')
        r = v/(v[0]*u)
        row[N] = (c05, c10, true_c1, r[:nt+1].min())
    ck((not np.isnan(row[400][0])) and abs(row[400][0]-row[700][0]) < 0.05 or np.isnan(row[400][0]),
       "N-consistency c(.5) om=%s" % om)
    print("%6.1f%9.2f%7.3f | %10.3f%10.3f | %9.3f%9.3f | %13.3f%12.3f" %
          (om, x, hlo, row[400][0], row[700][0], row[400][1], row[700][1], row[700][2], row[700][3]))

print("="*96)
print("ck PASS = %d    FAILS = %d" % (CK, len(FAILS)))
for nm, d in FAILS: print("   FAIL:", nm, d)
