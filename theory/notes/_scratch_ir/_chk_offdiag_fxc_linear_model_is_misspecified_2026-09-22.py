"""
team-lead, 2026-09-22 (from scratch): audit of the DOCSTRING MODEL in
_chk_offdiag_fxc_budget_for_sign_flip_2026-09-22.py.

That script's docstring says:

    f^xc_{ll} = 1, f^xc_{ll'} = kappa  (l != l')
    => Lambda_ij(kappa) = D_ij + kappa * R_ij      (LINEAR in kappa)
       kappa* = -D_ij/R_ij  is the off-diagonal f^xc strength needed to flip sign

CLAIM UNDER AUDIT: that linear model is NOT what an orbital off-diagonal f^xc
does to Lambda, because Lambda^SF = f^xc chi f^xc is a SANDWICH -- bilinear in
f^xc, not linear. With F(kappa) = I + kappa*K (K off-diagonal, K_ll = 0):

    F chi F = chi + kappa*(K chi + chi K) + kappa^2 * (K chi K)

so Lambda_ij(kappa) is QUADRATIC, and the linear coefficient is NOT the
off-diagonal part of chi (which is what the audited script uses as R_ij).

Two things are checked here, both exactly / deterministically:
 (A) symbolic-exact (Fraction): the audited script's R_ij  !=  the true linear
     coefficient of the sandwich, and a kappa^2 term exists.
 (B) the practical consequence: kappa*_linear (audited script) vs the true
     smallest positive root of the quadratic -- how far apart are they.

No repo data is loaded; this is a statement about the ALGEBRA of the model,
so a small explicit matrix suffices and can be checked by hand.
"""
from fractions import Fraction as F

N = 3

# An explicit symmetric "Im chi^s"-like matrix (values are arbitrary but fixed;
# the point is algebraic, not material-specific). Symmetric, mixed signs off-diag.
chi = [[F(4), F(-1), F(2)],
       [F(-1), F(3), F(-5)],
       [F(2), F(-5), F(6)]]

# Fermi-surface orbital weights for a pair (i,j): nonneg, as in the audited script.
wi = [F(1, 2), F(1, 4), F(1, 4)]
wj = [F(1, 5), F(3, 5), F(1, 5)]

# Uniform off-diagonal knob K: K_ll = 0, K_ll' = 1 for l != l'.
K = [[F(0) if a == b else F(1) for b in range(N)] for a in range(N)]
I = [[F(1) if a == b else F(0) for b in range(N)] for a in range(N)]


def matmul(A, B):
    return [[sum(A[a][c] * B[c][b] for c in range(N)) for b in range(N)] for a in range(N)]


def add(A, B, s=F(1)):
    return [[A[a][b] + s * B[a][b] for b in range(N)] for a in range(N)]


def bil(u, M, v):
    return sum(u[a] * M[a][b] * v[b] for a in range(N) for b in range(N))


# ---- (A) exact expansion of the sandwich -------------------------------------
# Lambda(kappa) = wi . [ (I+kK) chi (I+kK) ] . wj  =  c0 + c1*k + c2*k^2
c0 = bil(wi, chi, wj)
c1 = bil(wi, add(matmul(K, chi), matmul(chi, K)), wj)
c2 = bil(wi, matmul(matmul(K, chi), K), wj)

# ---- what the audited script actually uses -----------------------------------
# D = sum_l wi_l chi_ll wj_l ;  full1 = wi . chi . wj at kappa=1 ; R = full1 - D
D = sum(wi[a] * chi[a][a] * wj[a] for a in range(N))
full1 = bil(wi, chi, wj)
R_script = full1 - D

print("=== (A) exact algebra, Fraction, N=3 ===")
print(f"  sandwich expansion:  c0 = {c0}   c1 = {c1}   c2 = {c2}")
print(f"  audited script:      D  = {D}    R_script = {R_script}")
print(f"  c0 == D ?            {c0 == D}    <- FALSE: script's 'D' drops chi's off-diag at kappa=0")
print(f"  c1 == R_script ?     {c1 == R_script}  <- the linear coefficients differ")
print(f"  quadratic term c2 = {c2}  (nonzero: {c2 != 0}) -- absent from the script's model")
print()

# ---- (B) flip points ----------------------------------------------------------
# script:  Lambda_lin(k) = D + k*R_script  -> k*_lin = -D/R_script
# truth :  smallest positive root of c2 k^2 + c1 k + c0 = 0
print("=== (B) consequence for the reported kappa* ===")
if R_script != 0:
    k_lin = -F(D, 1) / R_script
    print(f"  kappa*_linear (script formula, as a Fraction) = {k_lin} = {float(k_lin):.6g}")
else:
    k_lin = None
    print("  R_script == 0, script reports no finite flip")

disc = c1 * c1 - 4 * c2 * c0
print(f"  true quadratic:  {c2}k^2 + {c1}k + {c0} = 0,  discriminant = {disc} "
      f"({'real roots' if disc >= 0 else 'NO real root -> sign NEVER flips'})")
if disc >= 0 and c2 != 0:
    import math
    r = math.sqrt(float(disc))
    roots = sorted(((-float(c1) - r) / (2 * float(c2)), (-float(c1) + r) / (2 * float(c2))))
    pos = [x for x in roots if x > 0]
    print(f"  real roots = {roots}   smallest positive = {pos[0] if pos else 'none'}")
    if pos and k_lin is not None:
        print(f"  ratio  kappa*_true / kappa*_linear = {pos[0] / float(k_lin):.4g}")

# ---- (C) is the mis-specification DIRECTIONAL? --------------------------------
# One matrix cannot establish a direction. Sweep exact random rationals and count
# the four cells of (linear model says flippable?) x (quadratic says flippable?).
# "flippable" = exists kappa > 0 with Lambda(kappa) of opposite sign to Lambda(0).
import random
rng = random.Random(20260922)


def rand_sym():
    M = [[F(0)] * N for _ in range(N)]
    for a in range(N):
        for b in range(a, N):
            v = F(rng.randint(-9, 9), rng.randint(1, 4))
            M[a][b] = v
            M[b][a] = v
    return M


def rand_w():
    raw = [F(rng.randint(0, 9) + 1, 10) for _ in range(N)]
    s = sum(raw)
    return [x / s for x in raw]


def quad_flips(c0q, c1q, c2q):
    """True iff some kappa>0 makes c2 k^2 + c1 k + c0 have sign opposite to c0."""
    if c0q == 0:
        return False
    if c2q == 0:
        return (c1q != 0) and ((-c0q / c1q) > 0)
    disc_ = c1q * c1q - 4 * c2q * c0q
    if disc_ < 0:
        return False
    import math as _m
    r_ = _m.sqrt(float(disc_))
    rts = [(-float(c1q) - r_) / (2 * float(c2q)), (-float(c1q) + r_) / (2 * float(c2q))]
    return any(x > 0 for x in rts)


cells = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
coef_mismatch = 0
TRIALS = 4000
for _ in range(TRIALS):
    ch = rand_sym()
    ui, uj = rand_w(), rand_w()
    a0 = bil(ui, ch, uj)
    a1 = bil(ui, add(matmul(K, ch), matmul(ch, K)), uj)
    a2 = bil(ui, matmul(matmul(K, ch), K), uj)
    Dq = sum(ui[a] * ch[a][a] * uj[a] for a in range(N))
    Rq = a0 - Dq
    if a1 != Rq or a0 != Dq:
        coef_mismatch += 1
    lin = (Dq != 0) and (Rq != 0) and ((-Dq / Rq) > 0)
    cells[(lin, quad_flips(a0, a1, a2))] += 1

print("=== (C) direction of the mis-specification, 4000 exact random trials ===")
print(f"  trials where (c0,c1) != (D,R_script): {coef_mismatch}/{TRIALS} "
      f"= {100*coef_mismatch/TRIALS:.1f}%")
print(f"  linear FLIPPABLE & quadratic FLIPPABLE     : {cells[(True, True)]}")
print(f"  linear FLIPPABLE & quadratic NOT flippable : {cells[(True, False)]}"
      "   <- FALSE 'route is open'")
print(f"  linear NOT      & quadratic FLIPPABLE      : {cells[(False, True)]}"
      "   <- FALSE 'route is closed'")
print(f"  linear NOT      & quadratic NOT            : {cells[(False, False)]}")
disagree = cells[(True, False)] + cells[(False, True)]
print(f"  total disagreement: {disagree}/{TRIALS} = {100*disagree/TRIALS:.1f}%"
      "  => the two models are NOT interchangeable, and BOTH error directions occur")
print()
print("VERDICT: the docstring's linear model mis-specifies BOTH coefficients and")
print("omits the kappa^2 term entirely, so the printed kappa* = |D/R| is not the")
print("f^xc off-diagonal strength it is labelled as. (Direction: this does NOT")
print("rescue the sign scan -- theorem (S1) already forbids a sandwich from")
print("producing elementwise positivity for non-rank-one f^xc.)")
