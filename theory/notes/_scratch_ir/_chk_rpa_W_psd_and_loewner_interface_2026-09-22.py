"""
team-lead, 2026-09-22.  Two independent re-verifications, both rational-exact
(fractions only, zero float in any load-bearing line):

  (A) chi0-rebuild's claim that Step 3 of
      literature/rpa_screening_psd_derivation_2026-09-21.md is unsound.
      I find TWO independent defects, and then a one-line repair that makes the
      note's own conclusion TRUE (and stronger than it claimed) in the static case.

  (B) rpa-deep-dive's PSD-order counterexample, and what it does / does not close.
      Their 2x2 is correct against (Z-dom) elementwise, but the LIVE consumption
      point (ML3-abs) still holds on it -- via Theorem (U).  And the Loewner
      hypothesis turns out to be EXACTLY condition (i) of Theorem (U) plus a
      free extra, so the route is half of (U), not a dead end.

Run:  OMP_NUM_THREADS=1 python _chk_rpa_W_psd_and_loewner_interface_2026-09-22.py
"""
from fractions import Fraction as F


def tr_det(M):
    return M[0][0] + M[1][1], M[0][0] * M[1][1] - M[0][1] * M[1][0]


def psd_sym(M):
    """PSD test for a SYMMETRIC 2x2 over exact rationals."""
    assert M[0][1] == M[1][0], "psd_sym called on a non-symmetric matrix"
    t, d = tr_det(M)
    return t >= 0 and d >= 0


def sym_part(M):
    o = (M[0][1] + M[1][0]) / 2
    return [[M[0][0], o], [o, M[1][1]]]


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def lincomb(A, B, s):
    return [[A[i][j] + s * B[i][j] for j in range(2)] for i in range(2)]


print("=" * 78)
print("(A) rpa_screening_psd_derivation_2026-09-21.md, Step 3")
print("=" * 78)

# --- defect A1: the note's own lemma (L.115-124) "A>=0, B>=0 => AB >=0" is false.
A = [[F(1), F(0)], [F(0), F(2)]]
B = [[F(1), F(1)], [F(1), F(1)]]
print("A =", A, " symmetric PSD?", psd_sym(A), "(eigs 1,2)")
print("B =", B, " symmetric PSD?", psd_sym(B), "(eigs 0,2)")
AB = matmul(A, B)
print("AB =", AB, " symmetric?", AB[0][1] == AB[1][0])
S = sym_part(AB)
t, d = tr_det(S)
print(f"sym(AB) = {S}  tr={t} det={d}  PSD? {psd_sym(S)}")
# explicit negative vector: <x,ABx> = (a+b)(a+2b)
a, b = F(-3, 2), F(1)
q = (a + b) * (a + 2 * b)
print(f"explicit x=({a},{b}):  <x,ABx> = (a+b)(a+2b) = {q}  -> NEGATIVE" if q < 0 else "??")
print("=> DEFECT A1 (rational-exact): the product of two PSD matrices need not be PSD;")
print("   the note's 'positivity transfer' lemma and its proof sketch are both false.\n")

# --- defect A2: even granting a transfer lemma, the premise fails.
#     The note's own L.65-67 gives Im chi0(w>0) <= 0, and L.93 gives chi0(q,0) < 0.
#     Terms of W = v + v chi0 v + v chi0 v chi0 v + ... therefore ALTERNATE.
v, c = F(3), F(2)          # v > 0,  chi0 = -c < 0   (scalar, static)
terms = [v * (-c * v) ** n for n in range(5)]
print("scalar static, v =", v, " chi0 =", -c)
print("Dyson terms v*(chi0 v)^n, n=0..4:", terms)
print("signs:", [("+" if t > 0 else "-") for t in terms])
W_closed = v / (1 - (-c) * v)
print(f"closed form W = v/(1 - v*chi0) = {W_closed} > 0  -- the CONCLUSION is true,")
print("=> DEFECT A2 (rational-exact): 'every term PSD => sum PSD' has a FALSE premise;")
print("   the series is alternating.  Conclusion survives, the argument does not.\n")

# --- the repair: W = (v^{-1} - chi0)^{-1}, static, multiband, one line.
#     v >> 0 and chi0 <= 0 (repo theorem: static Re chi0 is NSD)  =>  v^{-1}-chi0 >> 0  =>  W >> 0.
#     Numerical sanity on a 2x2 with a genuinely off-diagonal, indefinite-looking chi0.
vm = [[F(4), F(1)], [F(1), F(3)]]
chi = [[F(-5), F(2)], [F(2), F(-9)]]        # NSD: tr<0, det = 45-4 = 41 > 0
print("repair check:  v =", vm, "PSD?", psd_sym(vm), " chi0 =", chi, "NSD?", psd_sym([[-chi[0][0], -chi[0][1]], [-chi[1][0], -chi[1][1]]]))
detv = vm[0][0] * vm[1][1] - vm[0][1] * vm[1][0]
vinv = [[vm[1][1] / detv, -vm[0][1] / detv], [-vm[1][0] / detv, vm[0][0] / detv]]
Mrep = lincomb(vinv, chi, -1)               # v^{-1} - chi0
t, d = tr_det(Mrep)
print(f"v^-1 - chi0 = {Mrep}  tr={t} det={d}  PD? {t > 0 and d > 0}")
dM = Mrep[0][0] * Mrep[1][1] - Mrep[0][1] * Mrep[1][0]
W = [[Mrep[1][1] / dM, -Mrep[0][1] / dM], [-Mrep[1][0] / dM, Mrep[0][0] / dM]]
t, d = tr_det(W)
print(f"W = (v^-1-chi0)^-1 = {W}  tr={t} det={d}  PD? {t > 0 and d > 0}")
print("=> REPAIR (analytic, one line, holds for the FULL multiband matrix, static only):")
print("   v>>0, chi0<=0  =>  v^-1 - chi0 >> 0  =>  W = (v^-1-chi0)^-1 >> 0.")
print("   This does NOT give Im W(w>0) >= 0, which is the quantity the Coulomb-layer")
print("   task actually needed -- that stays open.\n")

print("=" * 78)
print("(B) rpa-deep-dive's PSD-order counterexample vs the LIVE consumption point")
print("=" * 78)

LZ = [[F(2), F(0)], [F(0), F(2)]]
LD = [[F(0), F(-19, 10)], [F(-19, 10), F(0)]]
print("Lam^Z =", LZ, "  Lam^Delta =", LD)
print("  Lam^Z - Lam^Delta PSD?", psd_sym(lincomb(LZ, LD, -1)),
      " Lam^Z + Lam^Delta PSD?", psd_sym(lincomb(LZ, LD, 1)),
      "  => Loewner order HOLDS")
zdom = all(abs(LD[i][j]) <= LZ[i][j] for i in range(2) for j in range(2))
print("  elementwise (Z-dom) |Lam^D_ij| <= Lam^Z_ij ?", zdom,
      " => their FAIL of (Loewner => Z-dom) is CORRECT")

# the physical split implied by their matrices
Lph = [[(LZ[i][j] + LD[i][j]) / 2 for j in range(2)] for i in range(2)]
Lsf = [[(LZ[i][j] - LD[i][j]) / 2 for j in range(2)] for i in range(2)]
print("\n  implied  Lam^ph =", Lph, "   Lam^sf =", Lsf)
print("  Theorem (U) (i)  Lam^sf PSD                 :", psd_sym(Lsf))
print("  Theorem (U) (ii) offdiag(Lam^ph+Lam^sf) >= 0:", LZ[0][1] >= 0 and LZ[1][0] >= 0)

w = [F(1), F(1)]
R = [(LZ[i][0] * w[0] + LZ[i][1] * w[1]) / w[i] for i in range(2)]
L = [[(R[i] if i == j else 0) - LD[i][j] for j in range(2)] for i in range(2)]
t, d = tr_det(L)
print(f"  LIVE consumption point (ML3-abs):  L = diag(R^Z) - Lam^D = {L}")
print(f"     tr={t} det={d}  L PSD? {psd_sym(L)}   <-- HOLDS on their own counterexample")

print("""
=> Reading.  Their counterexample kills the implication (Loewner => elementwise
   (Z-dom)); that implication is indeed false.  It does NOT close the route,
   because the live consumption point L>=0 survives on the very same matrices
   via Theorem (U).  And the Loewner hypothesis is, identically,

       -Lam^Z <= Lam^D <= Lam^Z   <=>   Lam^Z - Lam^D = 2 Lam^sf >= 0
                                   and  Lam^Z + Lam^D = 2 Lam^ph >= 0,

   i.e. EXACTLY "(i) of Theorem (U), plus Lam^ph >= 0 for free".  The PSD-order
   route therefore delivers half of (U) and contributes nothing to (ii); the
   whole gap is (ii), the off-diagonal SIGN of Lam^Z = Lam^ph + Lam^sf.
   [analytic, one line -- no computation needed for the equivalence itself]
""")
