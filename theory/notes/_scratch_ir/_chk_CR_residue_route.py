"""[numerical-evidence] A ROOT-FREE, MONOTONICITY-FREE route to C_R: partial fractions.

Phi_s(sigma) = s^2 sigma (5+s^2+sigma^2) / [(1+(s-sigma)^2)(1+(s+sigma)^2)(1+s^2)(1+sigma^2)]

has exactly six simple poles in sigma, all with |Im| = 1:
    sigma = s+i, s-i, -s+i, -s-i, +i, -i          (for s != 0)
so  Phi_s(sigma) = sum_j r_j/(sigma - p_j)  (deg num 3 < deg den 6) and

    Phi_s''(sigma) = sum_j 2 r_j / (sigma - p_j)^3
    V2[Phi_s] = int_0^inf |Phi_s''| <= 2 sum_j |r_j| int_0^inf dsigma/|sigma-p_j|^3
              <= 2 sum_j |r_j| int_{-inf}^{inf} dsigma/((sigma-a_j)^2+1)^{3/2}
              = 4 sum_j |r_j|            (since int (x^2+1)^{-3/2} dx = 2)

=> C_R <= sup_s 4 s sum_j |r_j|.   NO sign changes, NO M/m, NO root isolation,
and the sigma integral is done EXACTLY.  The whole of T1 collapses to a bound on
an explicit algebraic function of s alone.

Checks:
  A. residue reconstruction reproduces Phi_s (and Phi_s'') pointwise
  B. the residue bound dominates the directly-computed V2[Phi_s]  (must, else bug)
  C. sup_s 4 s sum|r_j|  vs the budget 20
  D. closed form of sum |r_j| via sympy, and its large-s / small-s limits

Run: OMP_NUM_THREADS=1 python _chk_CR_residue_route.py
"""
import os, io, contextlib, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
def _load(name, alias):
    spc = importlib.util.spec_from_file_location(alias, os.path.join(HERE, name))
    m = importlib.util.module_from_spec(spc)
    with contextlib.redirect_stdout(io.StringIO()):
        spc.loader.exec_module(m)
    return m
pr = _load("_chk_padR_rowsum.py", "pr")

s_, si_ = sp.symbols("s sigma", positive=True)
PHI = s_**2 * si_ * (5 + s_**2 + si_**2) / (
    (1 + (s_ - si_)**2) * (1 + (s_ + si_)**2) * (1 + s_**2) * (1 + si_**2))
I = sp.I
POLES = [s_ + I, s_ - I, -s_ + I, -s_ - I, I, -I]

print("== D. residues in closed form ==")
res = []
for p in POLES:
    r = sp.simplify(sp.limit(PHI * (si_ - p), si_, p))
    res.append(sp.simplify(r))
for p, r in zip(POLES, res):
    print(f"   pole {p}:  r = {sp.simplify(r)}")

rf = [sp.lambdify(s_, r, "numpy") for r in res]
phi = sp.lambdify((s_, si_), PHI, "numpy")
d2phi = sp.lambdify((s_, si_), sp.diff(PHI, si_, 2), "numpy")

def sum_abs_r(s):
    return float(sum(abs(complex(f(s))) for f in rf))

print("\n== A/B. reconstruction + domination, and C. the budget ==")
print(f"{'s':>10} {'|Phi-sum r/(x-p)|':>19} {'V2 direct':>11} {'4*sum|r|':>10}"
      f" {'ratio':>7} {'s*V2':>9} {'4s*sum|r|':>11}")
worst_rec, best = 0.0, (-1.0, None)
for s in (0.01, 0.05, 0.2, 0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 3.0, 5.0, 10.0,
          30.0, 100.0, 1000.0):
    x = np.array([0.3, 1.7, 4.2, 11.0])
    rec = np.zeros_like(x, dtype=complex)
    for p, r in zip(POLES, res):
        rec += complex(sp.N(r.subs(s_, s))) / (x - complex(sp.N(p.subs(s_, s))))
    err = float(np.abs(rec - phi(s, x)).max())
    worst_rec = max(worst_rec, err)
    V2, _ = pr.V2_Phi(float(s), ngrid=8001)
    B = 4.0 * sum_abs_r(s)
    if s * B > best[0]:
        best = (s * B, s)
    print(f"{s:10.3f} {err:19.2e} {V2:11.6f} {B:10.6f} {B/V2:7.3f}"
          f" {s*V2:9.6f} {s*B:11.6f}")
print(f"\n   worst reconstruction error = {worst_rec:.2e} (must be ~1e-15)")

ss = np.geomspace(1e-3, 1e4, 4001)
vals = np.array([4.0 * s * sum_abs_r(float(s)) for s in ss])
i = int(np.argmax(vals))
print(f"   grid sup of 4 s sum|r_j| = {vals[i]:.6f} at s = {ss[i]:.6f}"
      f"   (budget 20 -> margin {20/vals[i]:.2f}x)")
print(f"   large-s limit: {vals[-1]:.6f}   small-s: {vals[0]:.3e}")
