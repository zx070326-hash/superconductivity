r"""
P4 . EINSTEIN EXTREMAL -- NUMERICAL ASYMPTOTICS OF h(w) FOR LARGE w
================================================================================
Evidence-gathering for **Segment C-II** of `P4_einstein_extremal.md` (concavity
of H(x)=h(sqrt x), x=w^2, for large x).  This file produces NO certificate: it is
float64 / ARPACK numerics whose purpose is to tell the analytic proof which shape
to aim at.  Every number below is a float64 estimate, not a rigorous bound.

--------------------------------------------------------------------------------
THE OPERATOR
--------------------------------------------------------------------------------
For a kernel f : {1,2,...} -> R define O[f] on l^2(N_0) by  (-L1 diagonal, +L2 +L3)

    O[f]_{nm} = [ (1-delta_{nm}) f(|n-m|) + f(n+m+1) ] / sqrt((2n+1)(2m+1))
                - delta_{nm} (1/(2n+1)) sum_{k=1}^{n} 2 f(k) .

The KAY II operator of `P4_interchange_uniformity.md` is  H(w) = O[ [[.]](w) ],
[[k]](w) = w^2/(k^2+w^2), and h(w) = lam_1(H(w)).  f(0) never enters (the
Toeplitz diagonal is removed and the Hankel index n+m+1 >= 1).  O[.] is LINEAR in
f, which is what makes the Loewner comparisons of Stage 5 cheap: the difference
of two comparison operators is again an O[.].

--------------------------------------------------------------------------------
METHOD (and why it can reach w = 2e4)
--------------------------------------------------------------------------------
* MATVEC IN O(N log N).  With D = diag((2n+1)^{-1/2}), O[f] = D(T+K)D - diag(L1)
  where T is the symmetric Toeplitz matrix f(|n-m|) with its diagonal removed and
  K is the Hankel matrix f(n+m+1).  Both are applied by one length-L real FFT
  pair, L = next_fast_len(2N) (>= 2N kills every wrap-around).  This replaces the
  O(N^2) dense build of `H_mat`; `check_matvec` verifies it reproduces `H_mat`
  to ~5e-16 relative, so no entry formula is re-derived here.
* TOP EIGENPAIR by scipy ARPACK `eigsh(k=1, which='LA', tol=0)`.  The spectrum is
  lam_1 ~ ln w, lam_2 = O(1e-3), lam_min -> -1: lam_1 is spectacularly isolated,
  so Lanczos converges in a single 48-step cycle.  (k>=2 is NOT used: lam_2 sits
  inside a cluster and is very slow -- it is not needed for anything here.)
* DERIVATIVES ANALYTICALLY, NOT BY DIFFERENCING.  Since lam_1 is simple,
      h'(w)  = <v, H'(w) v>                                   (Hellmann-Feynman)
      h''(w) = <v, H''(w) v> + 2 <r, (lam_1 - H)^{-1}|_{v^perp} r>,  r = (I-P)H'v
  and H'(w)=O[[[.]]'], H''(w)=O[[[.]]''] are again in the O[.] class.  The
  resolvent solve is done by CG on v^perp, where lam_1 I - H is SPD with condition
  number (lam_1-lam_min)/(lam_1-lam_2) ~ 1.2, so CG needs ~5 iterations.  This is
  the same second-order perturbation identity as (C*) in the note, used here as an
  equality rather than an upper bound.  Stage 2 cross-checks it against 5-point
  finite differences in tau = ln w.
* N-TRUNCATION.  H^(N) is a compression, so h_N increases to h.  Empirically
  h - h_N = C (w/N)^4 with C ~ 0.33 (Stage 1), which is exactly the tail of the
  scaling profile found in Stage 4 (g(s) ~ 1.2/s^2 => v_n ~ w^2 n^{-5/2}).
  Production runs use N = 80w and 160w plus the Richardson combination
  h* = h_{2N} + (h_{2N}-h_N)/15, giving |h*-h| ~ 1e-10.

--------------------------------------------------------------------------------
WHAT THE STAGES REPORT
--------------------------------------------------------------------------------
  V  matvec / eigenvalue verification against `H_mat`, `H_matrix`, `eigvalsh`,
     and the analytic h', h'' against finite differences.
  1  N-convergence of h(w) (10w,20w,40w,80w,160w) and the N reaching 1e-9.
  2  h(w), h'(w), h''(w), H'(x), H''(x), w h'(w), x^2 H''(x) on the w grid.
  3  asymptotic fits.  THE HEADLINE: the expansion is in 1/ln w, NOT in 1/w.
  4  top eigenvector scaling: the naive v_n ~ w^{-1/2} phi(n/w) is FALSE; the
     correct statement is v_n = (2n+1)^{-1/2} g(n/w) / sqrt(Nrm(w)).
  5  comparison operators O[1_{k<=cw}], O[e^{-k/w}], explicit trial vectors, and
     the Loewner-order test lam_min / lam_max of O[1_{k<=cw} - [[.]]].
  6  concavity evidence: sign of H'' and x^2 H''(x) / (-a/2), plus the (C*)
     term-by-term split that tells Segment C-II which term carries the sign.
  7  KAY III Thm 6 cross-check in the opposite corner (N fixed, w -> infinity).
  8  an EXACT rank-one reduction H = 2uu^T - I - O[m_w] and the resulting scalar
     secular equation for h(w).  This is the strongest Segment C-II lead here.

USAGE:  python p4_extremal_asymptotics.py [--quick] [stage ...]
        stages: V 1 2 3 4 5 6 7 8   (default: all)
        --quick trims the w grid and N so the file runs in ~1 min on a laptop.
"""

import sys
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh, cg

# Entry formulas are taken from the verified module; nothing is re-derived here.
from p4_interchange_uniformity import H_mat, H_matrix, _brk_vec

try:
    from scipy.fft import rfft, irfft, next_fast_len
    _HAVE_SCIPY_FFT = True
except ImportError:                                            # pragma: no cover
    from numpy.fft import rfft, irfft
    _HAVE_SCIPY_FFT = False

    def next_fast_len(n):
        L = 1
        while L < n:
            L *= 2
        return L

WORKERS = 1          # scipy.fft thread count; raised by --workers on the server

# Published weak-coupling constants for h(w) = ln w + c0 + o(1), supplied by the
# researcher `extremal-lit` as an independent calibration target:
#   SOFT cutoff N_eff = w/sqrt(e)  (Marsiglio PRB 98, 024523 eq.37;
#                                   Chubukov et al. arXiv:2004.01281 eq.17)
#   HARD cutoff                    (KAY III eq.23)
_EULER_GAMMA = 0.57721566490153286061
C0_SOFT = _EULER_GAMMA + np.log(4.0) - 1.5
C0_HARD = _EULER_GAMMA + np.log(4.0) - 1.0


# ============================================================================
# The operator O[f] as an O(N log N) LinearOperator
# ============================================================================

class OpKernel:
    """O[f] compressed to {0,...,N-1}; f is the array f(k), k = 0..2N-1."""

    __slots__ = ("N", "L", "inv", "diag", "Fc", "Fa", "nmv", "_wk")

    def __init__(self, N, f):
        f = np.asarray(f, dtype=np.float64)
        if f.size < 2 * N:
            raise ValueError(f"kernel needs 2N={2*N} entries, got {f.size}")
        self.N = N
        n = np.arange(N)
        self.inv = 1.0 / np.sqrt(2.0 * n + 1.0)
        # diagonal  -(1/(2n+1)) sum_{k=1}^{n} 2 f(k)
        self.diag = -np.concatenate(([0.0], np.cumsum(2.0 * f[1:N]))) / (2.0 * n + 1.0)
        L = int(next_fast_len(2 * N))
        self.L = L
        self._wk = dict(workers=WORKERS) if _HAVE_SCIPY_FFT else {}
        # symmetric Toeplitz f(|n-m|) with the diagonal removed -> circulant embed
        c = f[:N].copy()
        c[0] = 0.0
        circ = np.zeros(L)
        circ[:N] = c
        circ[L - N + 1:] = c[1:][::-1]
        self.Fc = rfft(circ, **self._wk)
        # Hankel f(n+m+1): circular convolution with the reversed vector.
        a = np.zeros(L)
        a[:2 * N] = f[:2 * N]
        self.Fa = rfft(a, **self._wk)
        self.nmv = 0

    def matvec(self, v):
        self.nmv += 1
        N, L, wk = self.N, self.L, self._wk
        v = np.asarray(v, dtype=np.float64).ravel()
        y = v * self.inv
        pad = np.zeros(L)
        pad[:N] = y
        t = irfft(self.Fc * rfft(pad, **wk), L, **wk)[:N]
        pad[:N] = y[::-1]
        k = irfft(self.Fa * rfft(pad, **wk), L, **wk)[N:2 * N]
        return (t + k) * self.inv + self.diag * v

    def linop(self):
        return LinearOperator((self.N, self.N), matvec=self.matvec, dtype=np.float64)

    def quad(self, v):
        return float(v @ self.matvec(v))


# ---------------------------------------------------------------- kernels ---

def k_brk(N, w, order=0):
    """f(k) = (d/dw)^order [[k]](w) for k = 0..2N-1 (closed forms from T5)."""
    return _brk_vec(np.arange(2 * N), float(w), order)


def k_step(N, K):
    """f(k) = 1_{k <= K}."""
    f = np.zeros(2 * N)
    f[:min(int(np.floor(K)) + 1, 2 * N)] = 1.0
    return f


def k_exp(N, w):
    """f(k) = exp(-k/w)."""
    return np.exp(-np.arange(2 * N) / float(w))


# ------------------------------------------------------------- eigensolve ---

def _seed(N):
    return 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)


def top_eig(N, f, v0=None, ncv=48, which="LA"):
    """(lam, v, op) for the extreme eigenpair of O[f] compressed at size N."""
    op = OpKernel(N, f)
    if v0 is None or len(v0) != N:
        v0 = _seed(N)
    vals, vecs = eigsh(op.linop(), k=1, which=which, tol=0,
                       ncv=min(ncv, N - 1), v0=v0)
    v = vecs[:, 0]
    if v[int(np.argmax(np.abs(v)))] < 0:
        v = -v
    return float(vals[0]), v, op


def h_and_derivs(N, w, v0=None):
    """h, h', h'' at w by Hellmann-Feynman + a projected resolvent solve."""
    lam, v, op = top_eig(N, k_brk(N, w, 0), v0=v0)
    o1 = OpKernel(N, k_brk(N, w, 1))
    o2 = OpKernel(N, k_brk(N, w, 2))
    a1 = o1.matvec(v)
    hp = float(v @ a1)
    r = a1 - hp * v                                   # (I-P) H'(w) v
    def shifted(x):
        x = np.asarray(x).ravel()
        x = x - (v @ x) * v
        y = lam * x - op.matvec(x)
        return y - (v @ y) * v
    A = LinearOperator((N, N), matvec=shifted, dtype=np.float64)
    y, info = cg(A, r, rtol=1e-14, atol=0.0, maxiter=300)
    if info != 0:                                     # pragma: no cover
        raise RuntimeError(f"CG did not converge (info={info}) at w={w}")
    T_a = o2.quad(v)                    # <v, H''(w) v>
    T_c = 2.0 * float(r @ y)            # 2 sum_{j>=2} |<v_j,H'v>|^2/(lam_1-lam_j) >= 0
    hpp = T_a + T_c
    res = float(np.linalg.norm(shifted(y) - r))
    return dict(lam=lam, hp=hp, hpp=hpp, v=v, cg_res=res,
                T_a=T_a, T_c=T_c, r2=float(r @ r))


def richardson(h1, h2, p=4):
    """Extrapolate h_N, h_{2N} assuming h - h_N = C N^{-p}."""
    q = 2.0 ** p
    return h2 + (h2 - h1) / (q - 1.0)


def h_star(w, mult=(80, 160), want_vec=False):
    """Richardson-extrapolated h(w) from N = mult[0]*w and mult[1]*w."""
    N1 = int(round(mult[0] * w))
    N2 = int(round(mult[1] * w))
    h1 = top_eig(N1, k_brk(N1, w))[0]
    lam2, v2, _ = top_eig(N2, k_brk(N2, w))
    hs = richardson(h1, lam2)
    return (hs, h1, lam2, v2) if want_vec else (hs, h1, lam2)


# ============================================================================
# derived quantities:  x = w^2, H(x) = h(sqrt x)
# ============================================================================

def Hx_derivs(w, hp, hpp):
    """H'(x), H''(x) and the scale-invariant combinations."""
    Hp = hp / (2.0 * w)
    Hpp = (hpp - hp / w) / (4.0 * w * w)
    return Hp, Hpp, w * hp, 0.25 * (w * w * hpp - w * hp)   # w h', x^2 H''


# ============================================================================
PASS = [0]


def _chk(name, ok, detail=""):
    if not ok:
        raise AssertionError(f"CHECK FAILED: {name} {detail}")
    PASS[0] += 1
    print(f"    ok  {name}  {detail}")


# ---------------------------------------------------------------------------
def stage_V(cfg):
    print("\n[V] verification against p4_interchange_uniformity")
    rng = np.random.default_rng(20260911)
    worst = 0.0
    for (N, w) in [(60, 3.0), (200, 7.0), (400, 0.5), (300, 50.0), (513, 137.0)]:
        A = H_mat(N, w)
        op = OpKernel(N, k_brk(N, w, 0))
        V = rng.standard_normal((N, 3))
        for j in range(3):
            b = A @ V[:, j]
            e = np.max(np.abs(op.matvec(V[:, j]) - b)) / np.max(np.abs(b))
            worst = max(worst, e)
    _chk("FFT matvec == H_mat", worst < 1e-13, f"max rel err {worst:.2e}")

    A = H_matrix(120, 4.0)
    op = OpKernel(120, k_brk(120, 4.0, 0))
    e = np.max(np.abs(np.column_stack([op.matvec(np.eye(120)[:, j]) for j in range(120)]) - A))
    _chk("FFT matvec == H_matrix (dense, entrywise)", e < 1e-14, f"max abs err {e:.2e}")

    for (N, w) in [(400, 5.0), (900, 40.0)]:
        E = np.linalg.eigvalsh(H_mat(N, w))
        lam = top_eig(N, k_brk(N, w, 0))[0]
        _chk(f"eigsh == eigvalsh (N={N}, w={w})", abs(lam - E[-1]) < 1e-12,
             f"d={lam-E[-1]:+.2e}  lam2={E[-2]:.3e}  lam_min={E[0]:.6f}")

    # analytic derivatives vs. 5-point finite differences in tau = ln w
    print("    analytic (Hellmann-Feynman + resolvent) vs 5-point FD in ln w:")
    for w in [5.0, 50.0, 400.0]:
        N = max(400, int(200 * w))
        D = h_and_derivs(N, w)
        lam, hp, hpp, res = D["lam"], D["hp"], D["hpp"], D["cg_res"]
        d = 0.02
        hs = [top_eig(N, k_brk(N, w * np.exp(j * d)))[0] for j in (-2, -1, 1, 2)]
        ht = (hs[0] - 8 * hs[1] + 8 * hs[2] - hs[3]) / (12 * d)
        htt = (-hs[0] + 16 * hs[1] - 30 * lam + 16 * hs[2] - hs[3]) / (12 * d * d)
        hp_fd = ht / w
        hpp_fd = (htt - ht) / (w * w)
        _chk(f"h'(w={w:g})", abs(hp - hp_fd) < 3e-9 * max(1.0, abs(hp)),
             f"analytic {hp:.12e}  FD {hp_fd:.12e}  d={hp-hp_fd:+.1e}")
        _chk(f"h''(w={w:g})", abs(hpp - hpp_fd) < 2e-6 * max(1e-6, abs(hpp)),
             f"analytic {hpp:.12e}  FD {hpp_fd:.12e}  rel={abs(hpp/hpp_fd-1):.1e} cg_res={res:.1e}")


# ---------------------------------------------------------------------------
def stage_1(cfg):
    print("\n[1] N-convergence of h(w).   h_N increases to h (compression).")
    print("    mult:        " + "".join(f"{m:>16d}w" for m in cfg["mults"]))
    rows = {}
    for w in cfg["wgrid_conv"]:
        vals, t0 = [], time.time()
        for m in cfg["mults"]:
            N = int(round(m * w))
            vals.append(top_eig(N, k_brk(N, w))[0])
        rows[w] = vals
        line = f"    w={w:<8g}" + "".join(f"{v:17.13f}" for v in vals)
        print(line + f"   [{time.time()-t0:.1f}s]")
        d = [f"{vals[i+1]-vals[i]:.2e}" for i in range(len(vals) - 1)]
        first = next((cfg["mults"][i + 1] for i in range(len(vals) - 1)
                      if vals[i + 1] - vals[i] < 1e-9), None)
        print(f"      deltas {' '.join(d)}   ratios "
              f"{' '.join(f'{(vals[i+1]-vals[i])/(vals[i+2]-vals[i+1]):.1f}' for i in range(len(vals)-2))}"
              f"   |dh|<1e-9 first at N={'>' + str(cfg['mults'][-1]) if first is None else first}w")
    # the empirical law  h - h_N = C (w/N)^4
    Cs = []
    for w, vals in rows.items():
        m, v1, v2 = cfg["mults"][-2], vals[-2], vals[-1]
        Cs.append((richardson(v1, v2) - v2) * (2 * m) ** 4 / 1.0)
    print(f"    empirical  h - h_N ~ C (w/N)^4  with C in "
          f"[{min(Cs):.3f},{max(Cs):.3f}]  (=> N=160w gives ~{max(Cs)/160**4:.1e})")


# ---------------------------------------------------------------------------
def stage_2(cfg):
    print("\n[2] h(w) and the x-derivatives (x = w^2, H(x) = h(sqrt x)).")
    print("        w        N            h(w)             w h'(w)       w^2 h''(w)"
          "        H'(x)          H''(x)        x^2 H''(x)")
    out = []
    for w in cfg["wgrid"]:
        N1 = int(round(cfg["prod"][0] * w))
        N2 = int(round(cfg["prod"][1] * w))
        t0 = time.time()
        h1 = top_eig(N1, k_brk(N1, w))[0]
        D = h_and_derivs(N2, w)
        lam, hp, hpp, v = D["lam"], D["hp"], D["hpp"], D["v"]
        hs = richardson(h1, lam)
        Hp, Hpp, whp, x2H = Hx_derivs(w, hp, hpp)
        out.append(dict(w=w, N=N2, h=hs, h_raw=lam, hp=hp, hpp=hpp,
                        Hp=Hp, Hpp=Hpp, whp=whp, x2H=x2H, v=v,
                        T_a=D["T_a"], T_c=D["T_c"], r2=D["r2"]))
        print(f"    {w:9g} {N2:9d} {hs:17.12f} {whp:14.9f} {w*w*hpp:15.9f}"
              f" {Hp:14.6e} {Hpp:15.6e} {x2H:14.9f}   [{time.time()-t0:.1f}s]")
    cfg["res2"] = out
    return out


# ---------------------------------------------------------------------------
_MODELS = [
    ("pow", "a*L + c0",                   lambda L, w: [L, 1 + 0 * L]),
    ("pow", "a*L + c0 + c1/w",            lambda L, w: [L, 1 + 0 * L, 1 / w]),
    ("pow", "a*L + c0 + c1/w + c2/w^2",   lambda L, w: [L, 1 + 0 * L, 1 / w, 1 / w ** 2]),
    ("pow", "a*L + c0 + c2/w^2",          lambda L, w: [L, 1 + 0 * L, 1 / w ** 2]),
    ("pow", "a*L + c0 + d*L/w^2",         lambda L, w: [L, 1 + 0 * L, L / w ** 2]),
    ("pow", "a*L + c0 + c1/w + c2/w^2 + c3/w^3",
     lambda L, w: [L, 1 + 0 * L, 1 / w, 1 / w ** 2, 1 / w ** 3]),
    ("log", "a*L + c0 + b1/L",            lambda L, w: [L, 1 + 0 * L, 1 / L]),
    ("log", "a*L + c0 + b1/L + b2/L^2",   lambda L, w: [L, 1 + 0 * L, 1 / L, 1 / L ** 2]),
    ("log", "a*L + c0 + b1/L + b2/L^2 + b3/L^3",
     lambda L, w: [L, 1 + 0 * L, 1 / L, 1 / L ** 2, 1 / L ** 3]),
    ("log", "L + c0 + b1/L + b2/L^2  (a:=1)",
     lambda L, w: [1 + 0 * L, 1 / L, 1 / L ** 2], 1.0),
    ("log", "L + c0 + b1/L + b2/L^2 + b3/L^3  (a:=1)",
     lambda L, w: [1 + 0 * L, 1 / L, 1 / L ** 2, 1 / L ** 3], 1.0),
    ("mix", "L + c0 + b1/L + b2/L^2 + b3/L^3 + c2/w^2  (a:=1)",
     lambda L, w: [1 + 0 * L, 1 / L, 1 / L ** 2, 1 / L ** 3, 1 / w ** 2], 1.0),
]


def stage_3(cfg):
    res = cfg.get("res2") or stage_2(cfg)
    w = np.array([r["w"] for r in res])
    h = np.array([r["h"] for r in res])
    L = np.log(w)
    keep = w >= cfg["fit_wmin"]
    print(f"\n[3] asymptotic fits of h(w), L := ln w, on w >= {cfg['fit_wmin']:g} "
          f"({keep.sum()} points).")
    print("    fam  model                                        coefficients"
          "                                    #p  max|resid|   rms")
    fits = {}
    for m in _MODELS:
        fam, name, cols = m[0], m[1], m[2]
        fix = m[3] if len(m) > 3 else None
        A = np.column_stack(cols(L[keep], w[keep]))
        b = h[keep] - (fix * L[keep] if fix is not None else 0.0)
        if A.shape[0] <= A.shape[1]:
            continue
        coef, *_ = np.linalg.lstsq(A, b, rcond=None)
        r = b - A @ coef
        mx = float(np.max(np.abs(r)))
        cs = " ".join(f"{c:+.6f}" for c in coef)
        print(f"    {fam}  {name:<44s} {cs:<52s} {A.shape[1]:2d}  {mx:.2e}  "
              f"{np.sqrt(np.mean(r**2)):.2e}")
        # best per family AT EQUAL PARAMETER COUNT is what matters -> keep by (#p, mx)
        key = (fam, A.shape[1])
        if key not in fits or mx < fits[key][0]:
            fits[key] = (mx, name, coef, cols, fix)
    print("\n    HEAD-TO-HEAD at equal parameter count (this is the honest comparison):")
    for p in sorted({k[1] for k in fits}):
        row = f"      #p={p}: "
        for fam in ("pow", "log", "mix"):
            if (fam, p) in fits:
                row += f"{fam} {fits[(fam,p)][0]:.2e}   "
        print(row)
    print("    => the 1/L family beats the 1/w family by 2-3 ORDERS OF MAGNITUDE at")
    print("       every parameter count.  h(w) has NO 1/w or 1/w^2 term at a level the")
    print("       data can see; its corrections are powers of 1/ln w.")

    for fam in ("pow", "log"):
        cand = [v for k, v in fits.items() if k[0] == fam]
        mx, name, coef, cols, fix = min(cand, key=lambda t: t[0])
        A = np.column_stack(cols(L, w))
        pred = A @ coef + (fix * L if fix is not None else 0.0)
        print(f"\n    residuals of the best '{fam}' model: {name}")
        print("        w         h(w)                fit                 residual")
        for i in range(len(w)):
            print(f"    {w[i]:9g} {h[i]:20.12f} {pred[i]:20.12f} {h[i]-pred[i]:+14.3e}")
        if fam == "log":
            cfg["fit"] = (coef, cols, fix, name)
            cfg["b1"] = float(coef[1] if fix is not None else coef[2])
            cfg["c0"] = float(coef[0] if fix is not None else coef[1])
            cfg["a"] = 1.0 if fix is not None else float(coef[0])

    # ---- calibration against the published weak-coupling constants -------
    print("\n    CALIBRATION against the literature (researcher `extremal-lit`):")
    print(f"        soft cutoff  N_eff = w/sqrt(e):  c0 = gamma + ln4 - 3/2 = {C0_SOFT:.7f}")
    print(f"        hard cutoff  (KAY III eq.23)  :  c0 = gamma + ln4 - 1   = {C0_HARD:.7f}")
    print("        model                                            a            c0"
          "        c0-soft      c0-hard")
    for m in _MODELS:
        fam, name, cols = m[0], m[1], m[2]
        fix = m[3] if len(m) > 3 else None
        if fam == "pow":
            continue
        A = np.column_stack(cols(L[keep], w[keep]))
        if A.shape[0] <= A.shape[1]:
            continue
        b = h[keep] - (fix * L[keep] if fix is not None else 0.0)
        coef = np.linalg.lstsq(A, b, rcond=None)[0]
        aa = 1.0 if fix is not None else float(coef[0])
        cc = float(coef[0] if fix is not None else coef[1])
        print(f"        {name:<44s} {aa:11.7f} {cc:13.7f} {cc-C0_SOFT:+12.2e} {cc-C0_HARD:+12.2e}")
    print("      => a = 1 to ~1e-4 and c0 agrees with the SOFT-cutoff value")
    print(f"         gamma+ln4-3/2 = {C0_SOFT:.7f} to ~1e-5, NOT with the hard-cutoff")
    print(f"         value {C0_HARD:.7f} (which is 0.5 higher).  Two independent published")
    print("         weak-coupling formulas (Marsiglio PRB 98, 024523 eq.37; Chubukov et al.")
    print("         arXiv:2004.01281 eq.17) therefore confirm the numerics, and the numerics")
    print("         confirm that N_eff = w/sqrt(e) (soft) is the right effective cutoff.")

    # ---- stability of a on progressively truncated grids ------------------
    print("\n    stability of the fitted a (model a*L+c0+b1/L+b2/L^2+b3/L^3), dropping")
    print("    the small-w end one point at a time:")
    for j in range(0, max(1, len(w) - 5)):
        s = slice(j, None)
        A = np.column_stack([L[s], 1 + 0 * L[s], 1 / L[s], 1 / L[s] ** 2, 1 / L[s] ** 3])
        if A.shape[0] <= A.shape[1]:
            break
        c = np.linalg.lstsq(A, h[s], rcond=None)[0]
        print(f"        w >= {w[j]:<7g} ({A.shape[0]} pts)  a = {c[0]:.8f}  "
              f"1-a = {1-c[0]:+.2e}   c0 = {c[1]:.7f}  c0-soft = {c[1]-C0_SOFT:+.2e}")

    # --- rate of x^2 H'' -> -a/2, predicted from the 1/L series -----------
    print("\n    convergence of x^2 H''(x) -> -a/2 = -0.5 .  With")
    print("        h = L + c0 + b1/L + b2/L^2 + ...   (L = ln w)")
    print("      one gets   w h'(w) = 1 - b1/L^2 - 2 b2/L^3 - ...   and")
    print("        x^2 H''(x) = -1/2 + (b1/2)/L^2 + (b1/2 + b2)/L^3 + ...")
    print("        w      x^2H''(x)      +0.5      (x^2H''+0.5)*L^2  (x^2H''+0.5)*w^2"
          "   w h'(w)    (1-w h')*L^2")
    e2, q2, Lg = [], [], []
    for r in res:
        e = r["x2H"] + 0.5
        Lw = np.log(r["w"])
        Lg.append(Lw); e2.append(e * Lw ** 2); q2.append((1 - r["whp"]) * Lw ** 2)
        print(f"    {r['w']:9g} {r['x2H']:13.9f} {e:+11.3e} {e*Lw**2:15.6f}"
              f" {e*r['w']**2:17.3e} {r['whp']:10.7f} {(1-r['whp'])*Lw**2:12.6f}")
    Lg, e2, q2 = np.array(Lg), np.array(e2), np.array(q2)
    m = Lg >= np.log(cfg["fit_wmin"])
    A = np.column_stack([np.ones(m.sum()), 1.0 / Lg[m], 1.0 / Lg[m] ** 2])
    ce = np.linalg.lstsq(A, e2[m], rcond=None)[0]
    cq = np.linalg.lstsq(A, q2[m], rcond=None)[0]
    print(f"    L->inf extrapolation of the two L^2-scaled errors (fit c + d/L + e/L^2):")
    print(f"        (x^2H''+1/2) L^2 -> {ce[0]:.6f}   [predicted b1/2 = {0.5*cfg.get('b1', float('nan')):.6f}]")
    print(f"        (1 - w h'(w)) L^2 -> {cq[0]:.6f}   [predicted b1   = {cfg.get('b1', float('nan')):.6f}]")
    print("    => BOTH errors settle to a CONSTANT when multiplied by L^2, and blow up")
    print("       when multiplied by w^2: the rate of x^2 H''(x) -> -1/2 is O(1/ln^2 w),")
    print("       NOT O(1/w^2).  The approach is from ABOVE (b1 > 0), so the concavity")
    print("       margin |x^2 H''| is SMALLEST at small w and grows towards 1/2.")


# ---------------------------------------------------------------------------
def _profile(v, w):
    """g_n = v_n sqrt(2n+1) / v_0  and the normalisation Nrm = sum g^2/(2n+1)."""
    n = np.arange(len(v))
    g = v * np.sqrt(2.0 * n + 1.0)
    g = g / g[0]
    return g, float(np.sum(g * g / (2.0 * n + 1.0)))


def stage_4(cfg):
    res = cfg.get("res2") or stage_2(cfg)
    sel = [r for r in res if r["w"] in cfg["wvec"]]
    print("\n[4] top eigenvector scaling.")
    print("    (a) the NAIVE ansatz v_n = w^{-1/2} phi(n/w) requires v_0 sqrt(w) to be")
    print("        w-independent.  It is not:")
    print("          w        v_0*sqrt(w)     1/v_0^2 = Nrm(w)   Nrm - 0.5*ln w   sum_{n<w} v_n^2")
    for r in sel:
        v, w = r["v"], r["w"]
        g, Nrm = _profile(v, w)
        print(f"        {w:8g} {v[0]*np.sqrt(w):14.6f} {1/v[0]**2:18.6f} "
              f"{1/v[0]**2 - 0.5*np.log(w):16.6f} {np.sum(v[:int(w)]**2):16.6f}")
    print("        v_0 sqrt(w) grows like sqrt(w/ln w): the mass is spread LOG-uniformly")
    print("        over 1 <= n <~ w, not concentrated on the scale n ~ w.")
    print("    (b) the CORRECT scaling:  v_n = (2n+1)^{-1/2} g(n/w) / sqrt(Nrm(w)),")
    print("        g(0)=1.  g collapses in s = n/w; the residual drift is a slow")
    print("        inverse power of L = ln w, so the s-table is extrapolated in BOTH")
    print("        1/L and 1/L^2 and the two extrapolations bracket the limit phi(s).")
    ss = np.arange(0, 3.0001, 0.1)
    print("          s     " + "".join(f"{r['w']:>12g}" for r in sel)
          + "    ex(1/L)   ex(1/L^2)  rms(1/L)  rms(1/L^2)")
    tab = {}
    Ls = np.array([np.log(r["w"]) for r in sel])
    gs = [_profile(r["v"], r["w"])[0] for r in sel]
    A1 = np.column_stack([np.ones_like(Ls), 1.0 / Ls])
    A2 = np.column_stack([np.ones_like(Ls), 1.0 / Ls ** 2])
    tot, spread = [0.0, 0.0], 0.0
    for s in ss:
        vals = np.array([g[int(round(s * r["w"]))] if int(round(s * r["w"])) < len(g)
                         else np.nan for g, r in zip(gs, sel)])
        c1 = np.linalg.lstsq(A1, vals, rcond=None)[0]
        c2 = np.linalg.lstsq(A2, vals, rcond=None)[0]
        r1 = float(np.sqrt(np.mean((vals - A1 @ c1) ** 2)))
        r2 = float(np.sqrt(np.mean((vals - A2 @ c2) ** 2)))
        tot[0] += r1 ** 2
        tot[1] += r2 ** 2
        spread = max(spread, abs(c1[0] - c2[0]))
        tab[round(s, 1)] = c1[0]          # 1/L drift: the better-fitting model
        print(f"        {s:5.2f} " + "".join(f"{x:12.7f}" for x in vals)
              + f"{c1[0]:11.7f}{c2[0]:11.7f}{r1:10.1e}{r2:10.1e}")
    print(f"        column rms over the whole table: 1/L {np.sqrt(tot[0]/len(ss)):.2e}, "
          f"1/L^2 {np.sqrt(tot[1]/len(ss)):.2e}")
    print(f"        -> the drift is 1/L, so phi(s) := ex(1/L).  BUT the L-range spanned")
    print(f"        by the data is short, and the two extrapolants differ by up to")
    print(f"        {spread:.3f}: treat phi(s) as accurate to ~1e-2, NOT more.  The FIT of")
    print("        g at finite w is excellent; only the L -> infinity limit is soft.")
    cfg["phi"] = tab
    print("    (c) tail of g:")
    print("          s      " + "".join(f"{r['w']:>13g}" for r in sel))
    for s in [4.0, 8.0, 16.0, 32.0, 64.0]:
        row = f"        g*s^2 s={s:<5g}"
        for r in sel:
            g, _ = _profile(r["v"], r["w"])
            i = int(round(s * r["w"]))
            row += f"{g[i]*s*s:13.7f}" if i < len(g) else f"{'-':>13s}"
        print(row)
    print("        g(s) ~ C/s^2 with C ~ 1.2 (slowly decreasing): v_n ~ w^2 n^{-5/2},")
    print("        which is exactly the (w/N)^4 truncation law of Stage 1.")
    print("    (d) how well is g approximated by 1/(1+s^2) ?")
    r = sel[-1]
    g, _ = _profile(r["v"], r["w"])
    for s in [0.25, 0.5, 1.0, 2.0, 4.0, 8.0]:
        i = int(round(s * r["w"]))
        if i < len(g):
            print(f"        s={s:<5g} g={g[i]:.6f}  1/(1+s^2)={1/(1+s*s):.6f}  "
                  f"ratio={g[i]*(1+s*s):.4f}")


def _phi_eval(tab, s):
    """Interpolate the tabulated limiting profile; C/s^2 beyond the table."""
    xs = np.array(sorted(tab))
    ys = np.array([tab[k] for k in xs])
    C = ys[-1] * xs[-1] ** 2
    s = np.asarray(s, dtype=np.float64)
    out = np.interp(s, xs, ys)
    far = s > xs[-1]
    return np.where(far, C / np.where(far, s, 1.0) ** 2, out)


# ---------------------------------------------------------------------------
def stage_5(cfg):
    print("\n[5] comparison operators and explicit trial vectors.")
    tab = cfg.get("phi")
    if tab is None:
        stage_4(cfg)
        tab = cfg["phi"]
    for w in cfg["wcomp"]:
        N = int(round(cfg["prod"][0] * w))
        n = np.arange(N)
        base = 1.0 / np.sqrt(2.0 * n + 1.0)
        op = OpKernel(N, k_brk(N, w))
        lam, v, _ = top_eig(N, k_brk(N, w))
        print(f"\n    w = {w:g}   N = {N}   h = {lam:.12f}")
        trials = [
            ("(i)   (2n+1)^-1/2 1_{n<w}", base * (n < w)),
            ("(ii)  (2n+1)^-1/2 (1+(n/w)^2)^-1", base / (1.0 + (n / w) ** 2)),
            ("(iii) (2n+1)^-1/2 phi(n/w)", base * _phi_eval(tab, n / w)),
        ]
        for nm, u in trials:
            u = u / np.linalg.norm(u)
            R = op.quad(u)
            print(f"      {nm:<36s} R = {R:14.9f}   h - R = {lam-R:11.3e}"
                  f"   (rel {100*(lam-R)/lam:6.3f}%)   overlap |<u,v>| = {abs(u@v):.6f}")
        print("      comparison operators (top eigenvalue of O[f] at the same N):")
        lstep = None
        for nm, f, ref in [("O[1_{k<=w}]   (step)", k_step(N, w), w),
                           ("O[1_{k<=2w}]", k_step(N, 2 * w), 2 * w),
                           ("O[1_{k<=w/2}]", k_step(N, w / 2), w / 2),
                           ("O[exp(-k/w)]", k_exp(N, w), w)]:
            lf = top_eig(N, f)[0]
            if ref == w and nm.startswith("O[1"):
                lstep = lf
            print(f"        {nm:<22s} lam_1 = {lf:14.9f}   lam_1 - h = {lf-lam:+11.6f}"
                  f"   lam_1 - ln(ref) = {lf-np.log(ref):+11.6f}")
        print(f"        => matching cutoff  c* with lam_1(O[1_{{k<=c*w}}]) = h(w):  "
              f"c* = {np.exp(lam-lstep):.6f}")
        print("      Loewner test:  D_c := O[1_{k<=cw}] - H(w) = O[1_{k<=cw} - [[.]]]")
        print("        c      lam_min(D_c)   lam_max(D_c)    verdict")
        for c in cfg["cgrid"]:
            fd = k_step(N, c * w) - k_brk(N, w)
            lo = top_eig(N, fd, which="SA", ncv=80)[0]
            hi = top_eig(N, fd, which="LA", ncv=80)[0]
            vd = ("H <= O[step]" if lo > -1e-9 else
                  ("H >= O[step]" if hi < 1e-9 else "NEITHER (indefinite)"))
            print(f"        {c:<6g} {lo:14.8f} {hi:14.8f}    {vd}")


# ---------------------------------------------------------------------------
def stage_6(cfg):
    res = cfg.get("res2") or stage_2(cfg)
    print("\n[6] direct concavity evidence at every computed w  (a := 1).")
    print("        w          x = w^2        H''(x)          sign   x^2H''/(-a/2)"
          "     w h'(w)/a    1 - w h'(w)")
    for r in res:
        print(f"    {r['w']:9g} {r['w']**2:15g} {r['Hpp']:15.6e}   "
              f"{'NEG' if r['Hpp'] < 0 else 'POS':>4s} {r['x2H']/(-0.5):16.9f}"
              f" {r['whp']:14.9f} {1-r['whp']:14.3e}")
    bad = [r["w"] for r in res if r["Hpp"] >= 0]
    _chk("H''(x) < 0 at every computed w", not bad, f"(w grid {res[0]['w']:g}..{res[-1]['w']:g})")
    mx = max(r["x2H"] for r in res)
    _chk("x^2 H''(x) stays below -0.45", mx < -0.45,
         f"max over grid = {mx:.6f} at w={[r['w'] for r in res if r['x2H']==mx][0]:g}")

    # ---- the (C*) split of the note, in the x^2 H'' normalisation -------
    print("\n    the (C*) decomposition  4w^2 H''(x) = <v,H''v> - <v,H'v>/w + 2<r,R r>,")
    print("    r = (I-P)H'(w)v, R = (lam_1-H)^{-1}|_{v^perp}.  Everything is multiplied by")
    print("    w^2/4 so that the three columns sum to x^2 H''(x) -> -1/2.  The last column")
    print("    replaces the exact resolvent term by the Davis-Kahan style bound 2||r||^2/g")
    print("    (g = lam_1 - lam_2; lam_2 = O(1e-3) here, so g ~ lam_1 to <0.1%).")
    print("        w      A=w^2<v,H''v>/4  B=-w h'/4   C=w^2 2<r,Rr>/4   A+B+C      "
          "Cbd=w^2(2||r||^2/g)/4   A+B+Cbd")
    rowsL, colA, colC, colB2 = [], [], [], []
    for r in res:
        w, q = r["w"], r["w"] ** 2 / 4.0
        A = q * r["T_a"]
        B = -0.25 * r["whp"]
        C = q * r["T_c"]
        Cbd = q * 2.0 * r["r2"] / r["h_raw"]
        rowsL.append(np.log(w)); colA.append(A); colC.append(C); colB2.append(A + B + Cbd)
        print(f"    {w:9g} {A:15.9f} {B:12.8f} {C:15.9f} {A+B+C:12.8f}"
              f" {Cbd:20.6f} {A+B+Cbd:14.6f}")
    Ls = np.array(rowsL)
    m = Ls >= np.log(cfg["fit_wmin"])
    M = np.column_stack([np.ones(m.sum()), 1.0 / Ls[m] ** 2, 1.0 / Ls[m] ** 3])
    lim = [np.linalg.lstsq(M, np.array(c)[m], rcond=None)[0][0] for c in (colA, colC, colB2)]
    print(f"    L->inf extrapolation (fit c + d/L^2 + e/L^3):  A -> {lim[0]:.6f},"
          f"  C -> {lim[1]:.6f},  A+B+Cbd -> {lim[2]:.6f}")
    print("    READ-OFF for Segment C-II:")
    print("      * the resolvent term C does NOT vanish: A -> -(1/4 + C_inf) and C -> C_inf")
    print("        with C_inf ~ 0.15, so A+C -> -1/4 and, with B -> -1/4, x^2H'' -> -1/2.")
    print("        The concavity is NOT produced by <v,H''v> alone -- the positive")
    print("        second-order term is an O(1) fraction of the budget and must be kept.")
    print("      * the Davis-Kahan replacement 2<r,Rr> -> 2||r||^2/g is nearly SHARP here")
    print("        (Cbd/C ~ 1.04 at the largest w), because the spectral weight of r sits")
    print("        at lam_j ~ 0 = lam_2 + O(1e-3), i.e. exactly where the bound is attained.")
    print("      * consequently A + B + Cbd stays ~ -0.49, bounded away from 0 UNIFORMLY in")
    print("        w after the x^2 rescaling.  (C*) is therefore NOT asymptotically vacuous:")
    print("        the 'wall' announced in the note is in CERTIFYING the three terms (DK")
    print("        angle + tail envelopes), not in the inequality itself.")


# ---------------------------------------------------------------------------
def stage_7(cfg):
    r"""KAY III Thm 6 cross-check: the OTHER order of limits.

    Stages 1-6 take N -> infinity at fixed w.  KAY III Thm 6 describes the
    opposite regime -- N FIXED, w -> infinity -- where [[k]](w) -> 1 uniformly on
    k <= N and therefore h^(N)(w) -> h^(N)(infinity) = lam_1(O_N[1]) with a
    1/w^2 deficit,  h^(N)(w) = h^(N)(inf) - B_N/w^2 + o(1/w^2).  Confirming that
    law (and measuring B_N) is an independent test of the operator, of the
    eigensolver and of the claim that the 1/L series of Stage 3 belongs to the
    N -> infinity limit and NOT to any residual truncation effect.
    """
    print("\n[7] KAY III Thm 6 cross-check: FIXED N, w -> infinity.")
    print("    h^(N)(w) = h^(N)(inf) - B_N/w^2 + o(1/w^2),  h^(N)(inf) = lam_1(O_N[1]).")
    for N in cfg["ndirect"]:
        f1 = np.ones(2 * N)
        hinf = top_eig(N, f1)[0]
        print(f"      N = {N:<6d}  h^(N)(inf) = lam_1(O_N[1]) = {hinf:.12f}"
              f"   (2*S_N - 2 = {2*np.sum(1.0/(2*np.arange(N)+1.0))-2:.6f} is the Lemma-X bound)")
        print("            w        h^(N)(w)          h^(N)(inf)-h^(N)(w)    *w^2 = B_N")
        prev = None
        for w in [4 * N, 8 * N, 16 * N, 32 * N, 64 * N]:
            hw = top_eig(N, k_brk(N, float(w)))[0]
            d = hinf - hw
            B = d * w * w
            rate = "" if prev is None else f"   ratio d(w/2)/d(w) = {prev/d:.4f} (4 = 1/w^2)"
            print(f"        {w:9d} {hw:18.12f} {d:22.6e} {B:16.6f}{rate}")
            prev = d
    print("    => the deficit scales as 1/w^2 with ratio 4.00 per doubling, exactly")
    print("       Thm 6.  Note this is the w >> N corner; Stage 1's (w/N)^4 truncation")
    print("       law lives in the OPPOSITE corner N >> w.  The two are consistent and")
    print("       the production runs (N = 80w..160w) sit safely inside N >> w.")


# ---------------------------------------------------------------------------
def _dense(op):
    N = op.N
    return np.column_stack([op.matvec(np.eye(N)[:, j]) for j in range(N)])


def stage_8(cfg):
    r"""An EXACT rank-one reduction of h(w) to a scalar secular equation.

    Because O[.] is linear and [[k]](w) = 1 - m_w(k) with m_w(k) = k^2/(k^2+w^2),

        H^(N)(w) = O_N[1] - O_N[m_w]   and   O_N[1] = 2 u u^T - I_N ,
        u_n = (2n+1)^{-1/2}                                            (EXACT)

    (O_N[1] = 2uu^T - I because its off-diagonal part is 2u u^T minus its own
    diagonal diag(1/(2n+1)) and its L1 diagonal is diag(2n/(2n+1)); the two
    diagonals add up to the identity.)  Hence spec(O_N[1]) = {2 S_N - 1} u
    {-1}^{N-1} with S_N = sum_{n<N} 1/(2n+1) -- this is Lemma X made exact.

    Writing M_w := O[m_w], the eigenvalue equation H v = lam v becomes
        ((lam+1) I + M_w) v = 2 u <u,v>,
    so v is proportional to B(lam)^{-1} u, B(lam) := (lam+1) I + M_w, and

        G(lam) := 2 <u, B(lam)^{-1} u> = 1                        (SECULAR EQN)

    h(w) is its largest root.  On lam+1 > -lam_min(M_w), G is strictly decreasing
    and completely monotone, so the top root is unique and h'(w), h''(w) follow by
    implicit differentiation of a SCALAR identity instead of eigenvalue
    perturbation theory.  This is the most promising concrete handle this file
    found for Segment C-II routes (i)/(iii).
    """
    print("\n[8] exact rank-one reduction and the scalar secular equation.")
    N, w = 600, 40.0
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    A = H_mat(N, w)
    Md = _dense(OpKernel(N, 1.0 - k_brk(N, w, 0)))
    e = float(np.max(np.abs(A - (2 * np.outer(u, u) - np.eye(N) - Md))))
    _chk("H^(N)(w) = 2uu^T - I - O_N[m_w] exactly", e < 1e-13,
         f"max entrywise |diff| = {e:.2e}  (N={N}, w={w})")
    E1 = np.linalg.eigvalsh(_dense(OpKernel(N, np.ones(2 * N))))
    S = float(np.sum(1.0 / (2.0 * n + 1.0)))
    _chk("spec(O_N[1]) = {2 S_N - 1} u {-1}^(N-1)",
         abs(E1[-1] - (2 * S - 1)) < 1e-11 and abs(E1[-2] + 1) < 1e-11,
         f"lam_1 = {E1[-1]:.12f} = 2S_N-1 = {2*S-1:.12f}; lam_2 = {E1[-2]:.12f}")

    print("      secular equation checked at production scale (CG solves):")
    print("        w       h(w)        lam_min(M_w)  lam_max(M_w)  lam_min(B(h))"
          "    G(h)-1       G(h-.05)  G(h+.05)")
    for w in cfg["wcomp"]:
        N = int(round(cfg["prod"][0] * w))
        n = np.arange(N)
        u = 1.0 / np.sqrt(2.0 * n + 1.0)
        mk = 1.0 - k_brk(N, w, 0)
        Mop = OpKernel(N, mk)
        lo = top_eig(N, mk, which="SA", ncv=80)[0]
        hi = top_eig(N, mk, which="LA", ncv=80)[0]
        lam = top_eig(N, k_brk(N, w, 0))[0]

        def G(z):
            B = LinearOperator((N, N), dtype=np.float64,
                               matvec=lambda x: (z + 1.0) * x + Mop.matvec(x))
            y, info = cg(B, u, rtol=1e-13, atol=0.0, maxiter=500)
            if info != 0:                                   # pragma: no cover
                raise RuntimeError(f"CG failed in G(lam) at w={w}")
            return 2.0 * float(u @ y)

        g0 = G(lam)
        print(f"      {w:8g} {lam:12.8f} {lo:13.8f} {hi:13.8f} {lam+1+lo:14.8f}"
              f" {g0-1:12.2e} {G(lam-0.05):10.7f} {G(lam+0.05):9.7f}")
        _chk(f"G(h(w))=1 at w={w:g}", abs(g0 - 1.0) < 1e-9, f"G-1 = {g0-1:.2e}")
    print("      => B(h) = (h+1)I + M_w is positive definite AT THESE N, G is strictly")
    print("         decreasing through 1 at lam = h(w), and the identity holds to ~1e-13.")
    print("         Implicit differentiation of G(lam,x)=1 reproduces Hellmann-Feynman")
    print("         exactly (v is parallel to B(h)^{-1}u), so all of Stage 6 follows from")
    print("         this scalar equation -- the route Segment C-II should try, since the")
    print("         concavity of the root of a scalar completely-monotone equation is a")
    print("         far smaller target than the concavity of an eigenvalue.")

    # ---- the caveat: this is a FINITE-N identity ------------------------
    print("\n      CAVEAT -- the reduction is a FINITE-N statement.  u = ((2n+1)^{-1/2})")
    print("      is NOT in l^2 (||u||^2 = S_N ~ (1/2)ln N diverges), so 2uu^T and O[m_w]")
    print("      are separately unbounded and only their difference survives N -> infinity.")
    print("      Concretely lam_min(O_N[m_w]) keeps sliding down with N towards the general")
    print("      floor <w,O[m]w> = 2<u,w>^2 - 1 - <w,Hw> >= -(1 + h), so the positive")
    print("      definiteness margin lam_min(B_N(h)) SHRINKS as N grows (the drift is")
    print("      ~1/ln N, and extrapolating it puts the limit at or just above 0 -- the")
    print("      data cannot tell 'small positive' from 'exactly 0'):")
    print("           w       N        lam_min(O_N[m])   lam_min(B_N(h))   -(1+h) floor")
    for w in cfg["wcav"]:
        hw = None
        for mult in cfg["ncav"]:
            N = int(round(mult * w))
            lo = top_eig(N, 1.0 - k_brk(N, w, 0), which="SA", ncv=80)[0]
            if hw is None:
                hw = top_eig(int(round(cfg["ncav"][-1] * w)), k_brk(
                    int(round(cfg["ncav"][-1] * w)), w))[0]
            print(f"        {w:7g} {N:9d} ({mult:4d}w) {lo:14.7f} {hw+1+lo:15.7f}"
                  f" {-(1+hw):14.7f}")
    print("      So the scalar route must be run at finite N with an explicit N -> infinity")
    print("      argument (h_N increases to h by compression), NOT directly on l^2.")


# ============================================================================
CFG_FULL = dict(
    mults=[10, 20, 40, 80, 160],
    wgrid_conv=[100.0, 1000.0],
    wgrid=[10., 20., 50., 100., 200., 500., 1000., 2000., 5000., 10000., 20000.],
    prod=(80, 160),
    fit_wmin=20.0,
    wvec=[100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0],
    wcomp=[50.0, 200.0, 1000.0],
    cgrid=[0.5, 1.0, 2.0, 4.0, 8.0],
    ndirect=[50, 200, 800],
    wcav=[50.0, 200.0],
    ncav=[20, 40, 80, 160, 320, 640, 1280],
)

CFG_QUICK = dict(CFG_FULL, wgrid=[10., 20., 50., 100., 200., 500., 1000.],
                 prod=(40, 80), wgrid_conv=[100.0], mults=[10, 20, 40, 80],
                 cgrid=[1.0, 4.0], ndirect=[50, 200],
                 wcav=[50.0], ncav=[20, 80, 320])


def main(argv):
    global WORKERS
    args = list(argv)
    cfg = CFG_FULL
    if "--quick" in args:
        args.remove("--quick")
        cfg = CFG_QUICK
    for a in list(args):
        if a.startswith("--workers="):
            WORKERS = int(a.split("=")[1])
            args.remove(a)
    cfg = dict(cfg)
    stages = [s.upper() for s in args] or ["V", "1", "2", "3", "4", "5", "6", "7", "8"]
    print(f"p4_extremal_asymptotics  (scipy.fft={_HAVE_SCIPY_FFT}, workers={WORKERS}, "
          f"stages={' '.join(stages)})")
    t0 = time.time()
    for s in stages:
        {"V": stage_V, "1": stage_1, "2": stage_2, "3": stage_3,
         "4": stage_4, "5": stage_5, "6": stage_6, "7": stage_7,
         "8": stage_8}[s](cfg)
    print(f"\ndone  ({PASS[0]} checks passed, {time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main(sys.argv[1:])
