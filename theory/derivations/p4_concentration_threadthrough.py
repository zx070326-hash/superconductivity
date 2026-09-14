# -*- coding: utf-8 -*-
"""
P4 concentration bound -- THREAD-THROUGH certificate (SELF-CONTAINED, separate stage).

Companion to `p4_concentration_bound.py`.  Documents/certifies that the profile bound
   v_n >= c(a,x) v_0 u_n   on n <= a*sqrt(x)   (proved in P4_concentration_bound.md)
threads through the cone-coercivity machinery (`P4_cone_coercivity.md`, branch
`g2-cone-coercivity`, section 7.4 open gap #4 and Prop 4 / eq (8)) to deliver the
TRUE order  kappa_prof ~ 1/x  in the large-x regime, replacing that doc's own
kappa_prof ~ 1/(x ln x).  This is a ONE-WAY citation: the g2 branch is not modified.

THREE distinct kernels (never conflate):
   fk(k) = x/(k^2+x)               Lorentz kernel; A=O[fk], h=lambda1(A), v=top unit evec
   gt(k) = k^2/(k^2+x) = 1-fk       resolvent/Psi kernel (this file's c(a,x))
   gE(k) = 2 x k^2/(k^2+x)^3        cone-coercivity E-kernel: E=<v,O[gE]v>,
                                    Phi(x)=sum_{n,m} gE(n+m+1) u_n^2 u_m^2

Reduction certified here (all terms >= 0 by cone algebra, so restriction is legit):
   E(x) = <v,O[gE]v> >= kappa_L3[gE,v] = sum gE(n+m+1) w_n w_m           (their Thm 1)
        >= sum_{n,m<=a sqrt x} gE(n+m+1) w_n w_m
        >= c(a,x)^2 v_0^2 * Phi_{<=a sqrt x}(x)      using w_n = u_n v_n >= c v_0 u_n^2.
Log-cancellation: Phi ~ (ln x)/x and v_0^2 ~ 1/ln x  =>  c^2 v_0^2 Phi ~ const/x.

Run:  OMP_NUM_THREADS=1 python p4_concentration_threadthrough.py
"""
import numpy as np

CK = 0; FAILS = []
def ck(cond, name, detail=""):
    global CK
    if cond: CK += 1
    else: FAILS.append((name, str(detail)))

# ---------- from-scratch operator (identical convention to base cert) ----------
def u_vec(N):
    n = np.arange(N); return (2.0*n+1.0)**(-0.5)

def O_of(phi, N):
    n = np.arange(N); u = u_vec(N)
    I = n[:,None]; J = n[None,:]
    off = np.where(I != J, phi(np.abs(I-J)), 0.0) + phi(I+J+1)
    M = (u[:,None]*u[None,:]) * off
    kk = np.arange(1, N)
    csum = np.concatenate(([0.0], np.cumsum(phi(kk))))
    M[n, n] -= (2.0/(2.0*n+1.0)) * csum
    return M

def f_lor(x): return lambda k: x/(k*k + x)
def gtil(x):  return lambda k: (k*k)/(k*k + x)
def gE(x):    return lambda k: 2.0*x*(k*k)/((k*k + x)**3)
def one_k(k): return np.ones_like(np.asarray(k, dtype=float))

def top_eig(M):
    w, V = np.linalg.eigh(M); v = V[:,-1]
    if v[0] < 0: v = -v
    return w[-1], v

def build_kernels(x, N):
    u = u_vec(N); n = np.arange(N).astype(float); u2 = u**2
    I = n[:,None]; J = n[None,:]
    f = f_lor(x); gt = gtil(x)
    Foff = np.where(I!=J, f(np.abs(I-J)), 0.0); Fl3 = f(I+J+1)
    Goff = np.where(I!=J, gt(np.abs(I-J)), 0.0); Gl3 = gt(I+J+1)
    return u, u2, Foff, Fl3, Goff, Gl3

def eps_bound_sfree(x, N, hval):
    u, u2, Foff, Fl3, Goff, Gl3 = build_kernels(x, N)
    a = Foff@u2 + Fl3@u2
    coeff = u2 * np.minimum(1.0, a/hval)          # = B(m)/(2 s / hval); s cancels
    return (Goff@coeff + Gl3@coeff)/hval

def c_of(x, a_range, N, hval):   # TRUNCATED c = UPPER bound on true c (see base cert 4.1)
    eb = eps_bound_sfree(x, N, hval); sqx = np.sqrt(x)
    nt = int(round(a_range*sqx))
    if nt < 1: return float('nan')
    return 1.0 - eb[:nt+1].max()

# --- rigorous Psi m-tail (decay lemma) => TAIL-INCLUSIVE strict lower bound on c ---
def Ftot_ub(x): return 0.502*np.pi*np.sqrt(x)
def a_ub(x, m):
    m = np.asarray(m, float); F = Ftot_ub(x); f = f_lor(x)
    S = lambda K: 1.0 + 0.5*np.log(2.0*K + 1.0)
    return (f(np.ceil(m/2.0))*S(np.floor(m/2.0)) + 2.0*F/(m+1.0)
            + f(m+1.0)*S(m) + F/(2.0*m+3.0))
def tail_bnd(x, N, hlo, Mbig=200000):
    m = np.arange(N, Mbig+1, dtype=float); part = np.sum((1.0/(2.0*m+1.0))*a_ub(x, m))
    F = Ftot_ub(x); Mb = float(Mbig)
    R = 1.5*F/Mb + 1.5*x*(1.0+0.5*np.log(2.0*Mb+1.0))/(Mb*Mb)
    return (2.0/hlo**2)*(part + R)
def c_tail(x, a_range, Nc=2000):   # tail-inclusive rigorous LOWER bound on c(a,x); no eigh
    u = u_vec(Nc); A = O_of(f_lor(x), Nc); best = -1e9
    for L in [2,4,8,16,32,64,128,256,512,1024]:
        if L > Nc: break
        psi = u.copy(); psi[L:] = 0.0; best = max(best, (psi@(A@psi))/(psi@psi))
    eb = eps_bound_sfree(x, Nc, best); sqx = np.sqrt(x)
    nt = int(round(a_range*sqx))
    if nt < 1: return float('nan')
    return 1.0 - (eb[:nt+1].max() + tail_bnd(x, Nc, best))

def h_lower_from_A(A, N):
    u = u_vec(N); best = -1e9
    for L in [2,4,8,16,32,64,128,256,512,1024]:
        if L > N: break
        psi = u.copy(); psi[L:] = 0.0
        best = max(best, (psi@(A@psi))/(psi@psi))
    return best

def Phi_partial(x, N, a):
    """sum_{n,m<=a*sqrt x} gE(n+m+1) u_n^2 u_m^2 ; a=np.inf -> full."""
    u2 = u_vec(N)**2; g = gE(x)
    cut = N if np.isinf(a) else max(1, min(int(round(a*np.sqrt(x)))+1, N))
    idx = np.arange(cut); I = idx[:,None]; J = idx[None,:]
    G = g(I+J+1); w = u2[:cut]
    return float(w @ G @ w)

# ---------- (0) builder sanity: S1 (real check, not a tautology) ----------
N0 = 600; u0 = u_vec(N0)
ck(np.max(np.abs(O_of(one_k,N0) - (2*np.outer(u0,u0)-np.eye(N0)))) < 1e-12, "S1 builder ok")

# ---------- grid ----------
GRID = [25.0, 64.0, 144.0, 400.0, 900.0, 1600.0, 2500.0, 3600.0]
def Nof(x): return int(max(600, 34*np.sqrt(x)))

# cache the load-bearing quantities per x
h_, hlo_, v0_, c_half_, c_one_, Phi_half_, Phi_one_, Phi_full_ = {},{},{},{},{},{},{},{}
for x in GRID:
    N = Nof(x)
    A = O_of(f_lor(x), N); h, v = top_eig(A); v0 = v[0]
    hlo = h_lower_from_A(A, N)
    h_[x]=h; hlo_[x]=hlo; v0_[x]=v0
    c_half_[x] = c_tail(x, 0.5)          # TAIL-INCLUSIVE rigorous lower bound (Nc=2000)
    c_one_[x]  = c_tail(x, 1.0)
    Phi_half_[x] = Phi_partial(x, N, 0.5)
    Phi_one_[x]  = Phi_partial(x, N, 1.0)
    Phi_full_[x] = Phi_partial(x, N, np.inf)
    ck(hlo <= h + 1e-9 and hlo >= 0.9*h, "h_lo in [0.9h,h] x=%g" % x, "hlo/h=%.3f"%(hlo/h))

# ---------- (1) THRESHOLD ϖ>=5 : c(a=1/2) > 0 across x>=25, and >=0.20 at x=25 ----------
ck(c_half_[25.0] >= 0.19, "c_tail(a=1/2, x=25) >= 0.19 (threshold wpi>=5)", "c=%.4f"%c_half_[25.0])
ck(all(c_half_[x] > 0 for x in GRID), "c_tail(a=1/2) > 0 for all x>=25")
ck(all(c_half_[GRID[i+1]] > c_half_[GRID[i]] - 1e-3 for i in range(len(GRID)-1)),
   "c_tail(a=1/2) nondecreasing on GRID (coarse subseq; c is NOT monotone in x, see base 4.3)")

# ---------- (2) THRESHOLD ϖ>=12 : c(a=1) turns positive at x=144, negative below ----------
ck(c_one_[144.0] > 0, "c_tail(a=1, x=144) > 0 (threshold wpi>=12)", "c=%.4f"%c_one_[144.0])
ck(c_one_[64.0]  < 0, "c_tail(a=1, x=64) < 0 (below wpi=12 threshold)", "c=%.4f"%c_one_[64.0])
ck(all(c_one_[x] > 0 for x in GRID if x >= 144.0), "c_tail(a=1) > 0 for all x>=144")

# ---------- (3) LOG-CANCELLATION: x*v0^2*Phi is FLAT (order 1/x), not decaying ----------
K1 = {x: x*v0_[x]**2*Phi_one_[x]  for x in GRID}     # a=1 kernel constant
Kh = {x: x*v0_[x]**2*Phi_half_[x] for x in GRID}     # a=1/2
vals1 = [K1[x] for x in GRID]; valsh = [Kh[x] for x in GRID]
ck(max(vals1)/min(vals1) < 1.12, "x*v0^2*Phi(a=1) flat (max/min<1.12)",
   "range [%.4f, %.4f]" % (min(vals1), max(vals1)))
ck(all(0.35 <= v <= 0.40 for v in vals1), "x*v0^2*Phi(a=1) ~ 0.378", "%s"%[round(v,4) for v in vals1])
ck(max(valsh)/min(valsh) < 1.25, "x*v0^2*Phi(a=1/2) flat", "range [%.4f, %.4f]"%(min(valsh),max(valsh)))
ck(all(0.17 <= v <= 0.25 for v in valsh), "x*v0^2*Phi(a=1/2) ~ 0.20", "%s"%[round(v,4) for v in valsh])

# ---------- (4) the two log ingredients separately ----------
# v0^2 ~ 1/ln x : v0^2 * ln x rises and stays O(1) in [1.5, 4.5]
vlnx = [v0_[x]**2*np.log(x) for x in GRID]
ck(all(1.5 <= q <= 4.5 for q in vlnx), "v0^2 * ln x is O(1) (v0^2 ~ 1/ln x)", "%s"%[round(q,3) for q in vlnx])
ck(vlnx[-1] > vlnx[0], "v0^2 * ln x rising toward O(1) plateau")
# Phi ~ 0.125 (ln x)/x : the slope d(x*Phi)/d(ln x) is ~0.125 (leading coeff), robust to
# the subleading additive constant (the raw ratio x*Phi/ln x drifts, the slope does not).
xP = [x*Phi_full_[x] for x in GRID]; lx = [np.log(x) for x in GRID]
slopes = [(xP[i+1]-xP[i])/(lx[i+1]-lx[i]) for i in range(len(GRID)-1)]
ck(all(0.11 <= sl <= 0.14 for sl in slopes),
   "d(x*Phi_full)/d(ln x) ~ 0.125 (Phi ~ (ln x)/x)", "%s"%[round(sl,4) for sl in slopes])
ck(all(xP[i+1] > xP[i] for i in range(len(GRID)-1)), "x*Phi_full rising (Phi decays slower than 1/x)")

# ---------- (5) the actual reduction lower bound holds numerically (needs O[gE]) ----------
for x in [144.0, 400.0, 1600.0]:
    N = Nof(x); A = O_of(f_lor(x), N); h, v = top_eig(A); v0 = v[0]
    OgE = O_of(gE(x), N); E = float(v @ (OgE @ v))
    u = u_vec(N); w = u*v
    # kappa_L3 full via w (Theorem 1 layer):
    idx = np.arange(N); g = gE(x)
    G = g(idx[:,None]+idx[None,:]+1)
    kL3_full = float(w @ G @ w)
    ck(kL3_full <= E + 1e-9, "Thm1: kappa_L3 <= E x=%g" % x, "kL3/E=%.4f"%(kL3_full/E))
    # proven reduction: E >= c^2 v0^2 Phi_{<=a sqrt x}
    lb_half = max(c_half_[x],0.0)**2 * v0**2 * Phi_half_[x]
    lb_one  = max(c_one_[x],0.0)**2  * v0**2 * Phi_one_[x]
    ck(E >= lb_half - 1e-12, "E >= c(1/2)^2 v0^2 Phi_{<=sqx/2} x=%g" % x, "E=%.5f lb=%.5f"%(E,lb_half))
    ck(E >= lb_one  - 1e-12, "E >= c(1)^2 v0^2 Phi_{<=sqx} x=%g" % x, "E=%.5f lb=%.5f"%(E,lb_one))
    if abs(x-1600.0) < 1e-9:
        # cross-validation against cone-coercivity probe: kappa_L3,full/E = 0.582 (their [0.582,0.727])
        ck(0.55 <= kL3_full/E <= 0.60, "x-validate kL3/E(x=1600) ~ 0.582", "%.4f"%(kL3_full/E))
        # the improved thread-through/E is a positive constant (order-exact), unlike ~1/(x ln x)
        ck(lb_half/E > 0.05, "thread-through/E is a positive constant (a=1/2)", "%.4f"%(lb_half/E))

# ---------- (6) N-consistency of the flatness constant ----------
for x in [1600.0]:
    for N in [Nof(x), int(1.6*Nof(x))]:
        A = O_of(f_lor(x), N); _, v = top_eig(A); v0 = v[0]
        val = x*v0**2*Phi_partial(x, N, 1.0)
        c_one_[('N',N)] = val
    va = c_one_[('N',Nof(1600.0))]; vb = c_one_[('N',int(1.6*Nof(1600.0)))]
    ck(abs(va-vb)/va < 0.01, "N-consistency x*v0^2*Phi(a=1) x=1600", "%.5f vs %.5f"%(va,vb))

# ---------- report ----------
print("="*92)
print("THREAD-THROUGH into cone-coercivity gap #4  (E-kernel gE, Phi_{<=a sqrt x})")
print("%4s %6s %6s %7s %7s %7s %8s %8s %9s %9s" %
      ("wpi","x","v0^2","c_1/2","c_1","xPhi_f","xv0^2Phi1","xv0^2Ph.5","v0^2lnx","xPhi/lnx"))
for x in GRID:
    print("%4.0f %6.0f %6.4f %7.3f %7.3f %7.4f %8.4f %8.4f %9.4f %9.4f" %
          (np.sqrt(x), x, v0_[x]**2, c_half_[x], c_one_[x], x*Phi_full_[x],
           x*v0_[x]**2*Phi_one_[x], x*v0_[x]**2*Phi_half_[x],
           v0_[x]**2*np.log(x), x*Phi_full_[x]/np.log(x)))
print("="*92)
print("ck PASS = %d    FAILS = %d" % (CK, len(FAILS)))
for nm, d in FAILS: print("   FAIL:", nm, d)
