"""Does the |M2| ceiling survive on the PHYSICAL operator (Z channel included)?

Level: numerical-evidence (float64, finite N).  Direction bookkeeping is stated per test.

Background.  (G'-abs) as stated by team-lead proves
    lam1(OO_pair) <= g2 * rho(|M2|) / kappa           [pairing-only operator]
and then writes "==> Tc <= C_inf sqrt(rho(|M2|))".  That last arrow needs
    lam1(OO) <= lam1(OO_pair)                          [(ML3)]
whose proof is the Dirichlet identity for  L(k) = diag(R(k)) - Lt(k),  and that identity
gives L >= 0 only because Lt was assumed ELEMENTWISE NONNEGATIVE.  With signed Lt the
Laplacian is not PSD and (ML3) has no proof.  Two readings must be separated:

  (R1) ONE signed matrix drives both channels:  R_i = sum_j Lt_ij (can go negative).
  (R2) s+- physics: the Z channel carries  LtZ = Lam^ph + Lam^sf >= 0 elementwise,
       the PAIRING channel carries  LtD = Lam^ph - Lam^sf  (signed), so that
            |LtD_ij| <= LtZ_ij      (elementwise DOMINANCE).

Claim I derive here (analytic, checked below):
  under dominance |LtD| <= LtZ elementwise, the Dirichlet form of
  L(k) := diag(R^Z(k)) - LtD(k)  is >= 0, hence (ML3) survives sign changes.
  Proof: |x_i||x_j| <= (1/2)(w_j/w_i x_i^2 + w_i/w_j x_j^2), so
    sum_ij LtD_ij x_i x_j <= sum_ij LtZ_ij |x_i||x_j| <= sum_i (sum_j LtZ_ij w_j/w_i) x_i^2
                           = sum_i R^Z_i x_i^2 .
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.linalg import eigvalsh

import _chk_ml_build as B

rng = np.random.default_rng(20260918_2)
G2 = 1.318140212132   # repo value; my own O reproduces it (stage0, N=1600)
FAIL = 0


def report(tag, ok, detail=""):
    global FAIL
    if not ok:
        FAIL += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {tag}  {detail}")


def multiband_O2(lamZ, lamD, N, w):
    """Z channel driven by lamZ row sums, pairing driven by lamD.  Both (K,B,B) symmetric."""
    lamZ = np.asarray(lamZ, float)
    lamD = np.asarray(lamD, float)
    Bb = lamD.shape[1]
    w = np.asarray(w, float)
    u = B.u_vec(N)
    n = np.arange(N)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    R = np.einsum("kij,j->ki", lamZ, w) / w[None, :]
    uu = np.outer(u, u)
    OO = np.zeros((N, Bb, N, Bb))
    pair = np.where(D == 0, 0.0, 1.0)[:, :, None, None] * lamD[D] + lamD[S]
    OO += uu[:, None, :, None] * np.transpose(pair, (0, 2, 1, 3))
    for i in range(Bb):
        for j in range(Bb):
            OO[n, i, n, j] += u ** 2 * lamD[0][i, j]
    csumR = np.concatenate([np.zeros((1, Bb)), np.cumsum(R[1:N], axis=0)])
    for i in range(Bb):
        OO[n, i, n, i] -= u ** 2 * (R[0][i] + 2.0 * csumR[:, i])
    return OO.reshape(N * Bb, N * Bb)


def rho(M):
    return float(np.max(np.abs(eigvalsh(M))))


# self-check: multiband_O2 with lamZ == lamD must equal the stage0-validated builder
N0 = 40
Lt = np.array([[1.0, 0.4], [0.4, 0.7]])
lam = B.lam_factorized(Lt, 1.3, N0)
d = np.max(np.abs(multiband_O2(lam, lam, N0, [1.0, 2.2]) - B.multiband_O(lam, N0, w=[1.0, 2.2])))
report("builder self-check: O2(lam,lam) == validated multiband_O", d < 1e-14, f"maxdiff={d:.2e}")

print()
print("=" * 74)
print("(R1) ONE signed matrix for both channels: does the |M2| ceiling hold?")
print("  direction: lam1(P_N O P_N) <= lam1(O)  ==> numerator is a LOWER estimate")
print("             ==> ratio > 1 at finite N is a TRUE violation (FAIL is sound)")
print("=" * 74)
N = 100
worst = 0.0
worst_cfg = None
n_viol = 0
trials = 1500
for t in range(trials):
    Bb = 2
    a_, b_ = rng.uniform(0.05, 3.0, size=2)
    c_ = -rng.uniform(0.02, 3.0)                 # NEGATIVE interband entry
    Lt = np.array([[a_, c_], [c_, b_]])
    varpi = float(rng.uniform(0.1, 8.0))
    w = np.array([1.0, float(rng.uniform(0.2, 5.0))])
    lam = B.lam_factorized(Lt, varpi, N)
    OO = B.multiband_O(lam, N, w=w)
    ceil = G2 * varpi ** 2 * rho(np.abs(Lt))
    r = B.lam1(OO) / ceil
    if r > 1.0:
        n_viol += 1
    if r > worst:
        worst, worst_cfg = r, (a_, b_, c_, varpi, w[1])
print(f"  worst ratio = {worst:.4f}   violations = {n_viol}/{trials}")
print(f"  worst config (a,b,c,varpi,w2) = "
      f"({worst_cfg[0]:.3f}, {worst_cfg[1]:.3f}, {worst_cfg[2]:.3f}, "
      f"{worst_cfg[3]:.3f}, {worst_cfg[4]:.3f})")
report("(R1) signed single matrix VIOLATES the |M2| ceiling", n_viol > 0,
       f"{n_viol} violations, worst {worst:.3f}x")

# N-convergence of the worst (R1) config: a violation must not be a truncation artifact
a_, b_, c_, varpi, w2 = worst_cfg
print("  N-convergence of the worst violating config:")
for NN in (100, 200, 400, 800):
    lam = B.lam_factorized(np.array([[a_, c_], [c_, b_]]), varpi, NN)
    OO = B.multiband_O(lam, NN, w=[1.0, w2])
    print(f"    N={NN:4d}  ratio={B.lam1(OO) / (G2 * varpi**2 * rho(np.abs([[a_, c_], [c_, b_]]))):.6f}")

print()
print("=" * 74)
print("(R2) s+- structure: LtZ = ph+sf >= 0, LtD = ph-sf signed, |LtD| <= LtZ")
print("=" * 74)
worst2 = 0.0
worst2_cfg = None
nv2 = 0
trials2 = 1500
for t in range(trials2):
    Bb = int(rng.integers(2, 5))
    ph = rng.uniform(0.0, 1.5, size=(Bb, Bb))
    ph = 0.5 * (ph + ph.T)
    sf = rng.uniform(0.0, 1.5, size=(Bb, Bb))
    sf = 0.5 * (sf + sf.T)
    LtZ, LtD = ph + sf, ph - sf
    varpi = float(rng.uniform(0.1, 8.0))
    w = np.concatenate([[1.0], rng.uniform(0.2, 5.0, size=Bb - 1)])
    lamZ = B.lam_factorized(LtZ, varpi, N)
    lamD = B.lam_factorized(LtD, varpi, N)
    OO = multiband_O2(lamZ, lamD, N, w)
    ceil = G2 * varpi ** 2 * rho(np.abs(LtD))
    r = B.lam1(OO) / ceil
    if r > 1.0:
        nv2 += 1
    if r > worst2:
        worst2, worst2_cfg = r, (Bb, varpi)
report("(R2) dominance: |M2| ceiling holds", nv2 == 0,
       f"{nv2}/{trials2} violations, max ratio = {worst2:.9f}")

# same scan but compared against the SHARPER factorized ceiling h(varpi)*lmax(|LtD|)
worst2b = 0.0
nv2b = 0
for t in range(400):
    Bb = int(rng.integers(2, 5))
    ph = rng.uniform(0.0, 1.5, size=(Bb, Bb)); ph = 0.5 * (ph + ph.T)
    sf = rng.uniform(0.0, 1.5, size=(Bb, Bb)); sf = 0.5 * (sf + sf.T)
    LtZ, LtD = ph + sf, ph - sf
    varpi = float(rng.uniform(0.1, 8.0))
    w = np.concatenate([[1.0], rng.uniform(0.2, 5.0, size=Bb - 1)])
    OO = multiband_O2(B.lam_factorized(LtZ, varpi, N), B.lam_factorized(LtD, varpi, N), N, w)
    r = B.lam1(OO) / (B.h_of(varpi, N) * eigvalsh(np.abs(LtD))[-1])
    nv2b += r > 1.0 + 1e-9
    worst2b = max(worst2b, r)
report("(R2) dominance: even the SHARP factorized ceiling h*lmax(|LtD|) holds",
       nv2b == 0, f"{nv2b}/400 violations, max ratio = {worst2b:.9f}")

print()
print("=" * 74)
print("(ML3-abs) Laplacian PSD under elementwise dominance |LtD| <= LtZ")
print("=" * 74)
mn = np.inf
mn_ctrl = 0.0
for t in range(20000):
    Bb = int(rng.integers(2, 7))
    ph = rng.uniform(0.0, 2.0, size=(Bb, Bb)); ph = 0.5 * (ph + ph.T)
    sf = rng.uniform(0.0, 2.0, size=(Bb, Bb)); sf = 0.5 * (sf + sf.T)
    LtZ, LtD = ph + sf, ph - sf
    w = rng.uniform(0.1, 6.0, size=Bb)
    R = (LtZ @ w) / w
    # symmetrized Laplacian in the x-variables (see docstring)
    L = np.diag(R) - LtD
    mn = min(mn, eigvalsh(0.5 * (L + L.T))[0] / max(1.0, np.max(np.abs(L))))
    # control: break dominance by inflating one off-diagonal of |LtD| past LtZ
    LtD2 = LtD.copy()
    i, j = 0, 1
    LtD2[i, j] = LtD2[j, i] = LtZ[i, j] * rng.uniform(1.5, 6.0)
    L2 = np.diag(R) - LtD2
    mn_ctrl = min(mn_ctrl, eigvalsh(0.5 * (L2 + L2.T))[0] / max(1.0, np.max(np.abs(L2))))
report("dominance ==> L PSD (20000 draws, B=2..6, random w)", mn > -1e-12,
       f"min normalized eig = {mn:.3e}")
report("control: breaking dominance destroys PSD", mn_ctrl < -1e-3,
       f"min normalized eig = {mn_ctrl:.4f}")

print()
print("=" * 74)
print(f"TOTAL FAILURES: {FAIL}")
print("=" * 74)
