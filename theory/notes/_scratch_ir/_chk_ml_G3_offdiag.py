"""[numerical-evidence] Does G3 actually cover the OFF-DIAGONAL blocks under (Z-dom)?

(G'-3) needs  N_ij := ||O[Lam^Delta_ij]|| <= c * g2 * (moment)_ij / kappa.
G3 gives ||O[f]|| <= 3 lam1(O[f]) ONLY for f >= 0 non-increasing.
Under (Z-dom) the pairing kernel is  Lam^D_ij(k) = int 2w a2F^D_ij(w) dw/(w^2+(2 pi T k)^2)
with a2F^D = a2F^ph - a2F^sf SIGNED, so Lam^D_ij(k) can change sign IN k.
Then G3's hypothesis fails and lam1(O[f]) can even be negative -- "3 lam1" is not a bound.

Repair (mine): split on the SPECTRAL FUNCTION, not on k.
  a2F^D = a2F^ph - a2F^sf,  Lam^D = Lam^ph - Lam^sf,  each nonneg AND non-increasing in k.
  ||O[Lam^D]|| <= ||O[Lam^ph]|| + ||O[Lam^sf]||
               <= 3 lam1(O[Lam^ph]) + 3 lam1(O[Lam^sf])       [G3 twice]
               <= 3 (g2/kappa) (M2^ph + M2^sf)_ij = 3 (g2/kappa) M2^Z_ij.
So the off-diagonal input is the Z-CHANNEL moment M2^Z, not |M2^D|.
Under (Z-dom) at the spectral level a2F^ph, a2F^sf >= 0 this is automatic.

Test below, kappa = (2 pi T)^2 = 1, Einstein-like modes:
  a2F(w) = sum_r A_r delta(w - w_r)  ==>  Lam(k) = sum_r 2 A_r w_r/(w_r^2+k^2),  M2 = sum_r 2 A_r w_r.
DIRECTION: lam1(P_N O P_N) <= lam1(O) and ||P_N O P_N|| <= ||O||, so a violation seen at
finite N of an UPPER bound claim is only suggestive; here the quantity that must FAIL
(the naive bound) is compared against a LOWER estimate of ||O||, so "naive bound violated"
is sound, while "repaired bound holds" is supported but not proven by this script.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import importlib.util
import numpy as np

spec = importlib.util.spec_from_file_location(
    "B", os.path.join(os.path.dirname(os.path.abspath(__file__)), "_chk_ml_build.py"))
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)

rng = np.random.default_rng(31415)
N = 300
ks = np.arange(0, 2 * N + 2, dtype=float)
G2 = 1.318140212132


def lam_of(modes):
    """modes = [(A, w)]; returns (Lam(k) array, M2)."""
    f = np.zeros_like(ks)
    m2 = 0.0
    for A, w in modes:
        f += 2.0 * A * w / (w ** 2 + ks ** 2)
        m2 += 2.0 * A * w        # M2 = int 2w a2F dw = 2 A w  (Einstein: = lam w^2)
    return f, m2


def norm_and_lam1(f):
    ev = np.linalg.eigvalsh(B.scalar_O(f, N))
    return max(abs(ev[0]), ev[-1]), ev[-1]


print(f"{'case':>28} {'signchg':>8} {'||O||':>9} {'lam1':>9} "
      f"{'3g2|M2^D|':>10} {'3g2 M2^Z':>10} {'naive?':>7} {'rep?':>5}")
bad_naive = bad_rep = n = 0
worst_naive = 0.0
for _ in range(400):
    wph, wsf = rng.uniform(0.2, 4.0), rng.uniform(0.2, 4.0)
    Aph, Asf = rng.uniform(0.05, 1.5), rng.uniform(0.05, 1.5)
    fph, m2ph = lam_of([(Aph, wph)])
    fsf, m2sf = lam_of([(Asf, wsf)])
    fD = fph - fsf
    m2D = m2ph - m2sf
    nrm, l1 = norm_and_lam1(fD)
    naive = 3 * G2 * abs(m2D)
    rep = 3 * G2 * (m2ph + m2sf)
    sc = bool((fD[1:N] > 0).any() and (fD[1:N] < 0).any())
    n += 1
    if nrm > naive + 1e-12:
        bad_naive += 1
        worst_naive = max(worst_naive, nrm / naive if naive > 0 else np.inf)
    if nrm > rep + 1e-12:
        bad_rep += 1
    if bad_naive <= 3 and nrm > naive:
        print(f"{'ph(A,w)=(%.2f,%.2f) sf=(%.2f,%.2f)' % (Aph, wph, Asf, wsf):>28} "
              f"{str(sc):>8} {nrm:9.4f} {l1:9.4f} {naive:10.4f} {rep:10.4f} "
              f"{'VIOL':>7} {'ok' if nrm <= rep else 'VIOL':>5}")
print(f"\n{n} draws: naive bound 3*g2*|M2^D| violated {bad_naive} times "
      f"(worst factor {worst_naive:.2f}x); repaired bound 3*g2*M2^Z violated {bad_rep} times")

# how often does lam1(O[Lam^D]) go NEGATIVE (G3's hypothesis lam1>0 and f>=0 both dead)?
neg = 0
for _ in range(400):
    wph, wsf = rng.uniform(0.2, 4.0), rng.uniform(0.2, 4.0)
    Aph, Asf = rng.uniform(0.05, 1.5), rng.uniform(0.05, 1.5)
    fD = lam_of([(Aph, wph)])[0] - lam_of([(Asf, wsf)])[0]
    if norm_and_lam1(fD)[1] < 0:
        neg += 1
print(f"lam1(O[Lam^D]) < 0 in {neg}/400 draws  ==> '3 lam1' is not even a positive bound there")

# sanity: G3 itself on each nonneg piece
bad = 0
for _ in range(300):
    A, w = rng.uniform(0.05, 1.5), rng.uniform(0.2, 4.0)
    f, m2 = lam_of([(A, w)])
    nrm, l1 = norm_and_lam1(f)
    if nrm > 3 * l1 + 1e-10 or l1 > G2 * m2 + 1e-10:
        bad += 1
print(f"G3 + lemma M on each nonneg piece: {bad}/300 violations")
