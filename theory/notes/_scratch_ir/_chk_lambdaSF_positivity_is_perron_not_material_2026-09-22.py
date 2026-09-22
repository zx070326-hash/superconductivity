"""
Is the 90% elementwise-positive rate of Im Lambda^SF_{kk'} INFORMATIVE?

essenberger's corrected scan (misc/fermi_surface_lambdaSF_offdiagonal_fxc_sensitivity
_2026-09-22.md sec 3') reports 89.3-96.7% of off-diagonal Im Lambda^SF > 0 and reads it
as "direction agrees with Eq.(37)".  Theorem (S1) (meta/fxc_sandwich_cannot_give_
elementwise_positivity_2026-09-22.md) measured 86.78% for a RANDOM entrywise-nonneg F
sandwiching a RANDOM PSD chi.  Those are suspiciously close -- but they are different
statistics (S1 = "all N^2 entries nonneg in a trial"; theirs = "fraction of entries").
This script computes the matched null baselines and asks three questions.

Their construction, verbatim from the script:
    Lambda_ij = 2 * (f*w_i) @ chi(k_i - k_j, omega) @ (f*w_j),   w_l = |U_k[l,n]|^2 >= 0
So Lambda_ij is a BILINEAR form in two NONNEGATIVE PROBABILITY VECTORS (sum_l w_l = 1).

Q1 (null A, structural): for random PSD M and random probability vectors u,v, what
    fraction of u^T M v are positive?  If ~90%, their number carries no material info.
Q2 (null B, material-anchored, sharper): keep their REAL Im chi^s(q,omega) but REPLACE
    the real orbital weights by random probability vectors.  If the positive fraction
    survives, all the content is in chi's spectral structure -- which is already visible
    in the orbital basis -- and NONE is in the Fermi-surface band structure.
Q3 (f^xc leverage / DOF): with orbital-diagonal f^xc and f_l = +-1, the l == l' terms
    carry f_l^2 = +1 and are INVARIANT under correlation_flip.  How much of Lambda_ij
    lives in that invariant part?  If it dominates, "robust across the two f^xc branches"
    is a degrees-of-freedom artifact, not evidence of robustness.

Strictness: numerical-evidence.
"""
import numpy as np, importlib.util, os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
def load(f):
    s = importlib.util.spec_from_file_location(f, os.path.join(SCRIPT_DIR, f))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m

rng = np.random.default_rng(20260922)
N = 5

def rand_prob(n):
    x = rng.random(n); return x / x.sum()

print("="*74)
print("Q1  NULL A: random PSD M, random probability vectors u,v -> sign of u^T M v")
print("="*74)
for label, gen in [("M = Wishart PSD (generic)", lambda: (lambda G: G@G.T)(rng.standard_normal((N,N)))),
                   ("M = PSD, strong Perron gap", lambda: (lambda w: np.outer(w,w)*5 + (lambda G: G@G.T)(rng.standard_normal((N,N)))*0.2)(abs(rng.standard_normal(N))))]:
    pos = tot = 0
    for _ in range(4000):
        M = gen(); ws = [rand_prob(N) for _ in range(10)]
        for i in range(10):
            for j in range(i+1, 10):
                v = ws[i] @ M @ ws[j]; tot += 1; pos += (v > 0)
    print(f"  {label:30s}: {pos}/{tot} = {100*pos/tot:5.1f}% positive")
print("  (their measured off-diagonal positive fraction: 89.3-96.7%)")

print()
print("="*74)
print("Q2  NULL B (material-anchored): REAL Im chi^s(q,omega), RANDOMIZED orbital weights")
print("="*74)
gate = load("_chk_gauge_invariance_gate_dynamic_2026-09-22.py")
ham  = load("_graser1003_hamiltonian_2026-09-22.py")
H_func = ham.build_H_graser
mu = float(np.load(os.path.join(SCRIPT_DIR, "_graser_mu_value.npy")))
print(f"  mu = {mu:.4f} eV   (same value the deliverable used)")

def chi0(q, w):   # sign-corrected, same wrapper the deliverable uses
    return -gate.chi0_dynamic(q, w, H_func, mu, T=0.02, eta=0.02, k_grid=16, norb=N, rephase=False)
def dress(c, U=0.3):
    return c @ np.linalg.inv(np.eye(N) - (U*np.eye(N)) @ c)

# reproduce their FS point set exactly (k_grid=36, fs_tol=0.03, 25 pts, seed 20260922)
ks = np.linspace(-np.pi, np.pi, 36, endpoint=False)
pts = []
for kx in ks:
    for ky in ks:
        e, U = np.linalg.eigh(H_func(kx, ky))
        for n in range(N):
            if abs(e[n]-mu) < 0.03:
                pts.append((kx, ky, np.abs(U[:, n])**2))
print(f"  FS points found: {len(pts)}  (deliverable reported 73)")
sub = np.random.default_rng(20260922).choice(len(pts), size=25, replace=False)
pts = [pts[i] for i in sub]

for omega in (0.05, 0.2):
    cache = {}
    def chiq(qx, qy):
        k = (round(qx,6), round(qy,6))
        if k not in cache: cache[k] = dress(chi0((qx,qy), omega))
        return cache[k]
    real_pos = real_tot = 0
    rand_pos = rand_tot = 0
    rw = [rand_prob(N) for _ in pts]          # randomized weights, same count
    gaps, coh = [], 0
    for i in range(len(pts)):
        for j in range(i+1, len(pts)):
            M = chiq(pts[i][0]-pts[j][0], pts[i][1]-pts[j][1]).imag
            real_tot += 1; real_pos += (pts[i][2] @ M @ pts[j][2] > 0)
            rand_tot += 1; rand_pos += (rw[i] @ M @ rw[j] > 0)
    for (kx,ky,_) in pts[:12]:
        M = chiq(kx-pts[0][0], ky-pts[0][1]).imag
        ev, V = np.linalg.eigh((M+M.T)/2)
        k = int(np.argmax(abs(ev)))
        gaps.append(abs(ev[k])/max(abs(np.delete(ev,k)).max(), 1e-30))
        v = V[:,k]; coh += ((v>=0).all() or (v<=0).all())
    print(f"  omega={omega}:  REAL weights {real_pos}/{real_tot} = {100*real_pos/real_tot:5.1f}%"
          f"   |  RANDOM weights {rand_pos}/{rand_tot} = {100*rand_pos/rand_tot:5.1f}%")
    print(f"           Im chi^s Perron gap median {np.median(gaps):7.2f},"
          f"  top eigvec sign-coherent {coh}/12")

print()
print("="*74)
print("Q3  f^xc LEVERAGE: the l==l' part carries f_l^2=+1 and CANNOT be flipped")
print("="*74)
f_uni  = np.array([1.,1.,1.,1.,1.])
f_flip = np.array([1.,1.,-1.,1.,-1.])
for omega in (0.05, 0.2):
    cache = {}
    def chiq(qx, qy):
        k = (round(qx,6), round(qy,6))
        if k not in cache: cache[k] = dress(chi0((qx,qy), omega))
        return cache[k]
    frac, flipped, tot = [], 0, 0
    for i in range(len(pts)):
        for j in range(i+1, len(pts)):
            M = chiq(pts[i][0]-pts[j][0], pts[i][1]-pts[j][1]).imag
            wi, wj = pts[i][2], pts[j][2]
            D = float(np.sum(wi*np.diag(M)*wj))                 # l==l', INVARIANT
            full_u = wi @ M @ wj
            full_f = (f_flip*wi) @ M @ (f_flip*wj)
            frac.append(abs(D)/max(abs(full_u), 1e-30))
            tot += 1; flipped += (np.sign(full_u) != np.sign(full_f))
    print(f"  omega={omega}:  median |invariant l==l' part| / |total| = {np.median(frac):.3f}"
          f"   ;  pairs whose SIGN changes under correlation_flip: {flipped}/{tot}"
          f" = {100*flipped/tot:.1f}%")
print()
print("If the invariant part dominates and few signs move, then 'robust across the two")
print("f^xc branches' measures an algebraic invariance, not physical robustness.")

print()
print("="*74)
print("Q4  PERMUTATION NULL (fair): keep the REAL weight vectors and the REAL q set,")
print("    only break the pairing between orbital character and momentum.")
print("    This preserves the concentration of w exactly (Dirichlet null in Q2 does not).")
print("="*74)
prng = np.random.default_rng(7)
for omega in (0.05, 0.2):
    cache = {}
    def chiq(qx, qy):
        k = (round(qx,6), round(qy,6))
        if k not in cache: cache[k] = dress(chi0((qx,qy), omega))
        return cache[k]
    Ms = {}
    for i in range(len(pts)):
        for j in range(i+1, len(pts)):
            Ms[(i,j)] = chiq(pts[i][0]-pts[j][0], pts[i][1]-pts[j][1]).imag
    W = [p[2] for p in pts]
    obs = sum((W[i] @ Ms[(i,j)] @ W[j] > 0) for (i,j) in Ms)
    null = []
    for _ in range(2000):
        pi = prng.permutation(len(pts))
        null.append(sum((W[pi[i]] @ Ms[(i,j)] @ W[pi[j]] > 0) for (i,j) in Ms))
    null = np.array(null)
    p = (null <= obs).mean()
    print(f"  omega={omega}:  observed {obs}/300 = {100*obs/300:5.1f}%   null median "
          f"{np.median(null):.0f}/300 = {100*np.median(null)/300:5.1f}%  "
          f"[2.5-97.5%: {np.percentile(null,2.5):.0f}-{np.percentile(null,97.5):.0f}]"
          f"   P(null <= obs) = {p:.4f}")

print()
print("="*74)
print("Q5  WHAT ACTUALLY DRIVES IT: sign of diag(Im chi^s)")
print("="*74)
for omega in (0.05, 0.2):
    cache = {}
    dpos = dtot = 0; offneg = offtot = 0
    for i in range(0, len(pts), 2):
        for j in range(i+1, len(pts), 3):
            M = dress(chi0((pts[i][0]-pts[j][0], pts[i][1]-pts[j][1]), omega)).imag
            d = np.diag(M); dpos += (d > 0).sum(); dtot += N
            off = M[~np.eye(N, dtype=bool)]; offneg += (off < 0).sum(); offtot += off.size
    print(f"  omega={omega}:  diag(Im chi^s) > 0 : {dpos}/{dtot} = {100*dpos/dtot:5.1f}%"
          f"   |  off-diag < 0 : {offneg}/{offtot} = {100*offneg/offtot:5.1f}%")
print()
print("w nonneg + diag(Im chi^s) sign-definite => the l==l' term (99.8% of the amplitude,")
print("Q3) has that same fixed sign for EVERY pair.  Elementwise positivity of Lambda^SF is")
print("inherited from a diagonal sign convention, not from Fermi-surface physics.")

print()
print("="*74)
print("Q6  THEOREM (D) (analytic), verified numerically")
print("="*74)
print("""  Let f^xc be orbital-diagonal and REAL, f^xc_{ll'} = f_l delta_{ll'}, and let the
  orbital weights w(i) = |U_{k_i}[:,n_i]|^2 be (manifestly) nonnegative.  Then the
  l == l' contribution to Im Lambda^SF_ij is

      D_ij = sum_l  f_l^2 * w_l(i) * w_l(j) * Im chi^s_ll(q_ij)

  in which f_l^2 >= 0 and w_l(i) w_l(j) >= 0 TERMWISE.  Hence if diag(Im chi^s) is
  sign-definite with sign sigma, D_ij has sign sigma for EVERY pair (i,j) and for
  EVERY real choice of f_l -- including every sign pattern.  QED (analytic)

  Consequence: an f^xc sensitivity sweep that only flips SIGNS of an orbital-diagonal
  kernel has provably zero leverage on D.  It can only move the l != l' remainder.""")
for omega in (0.05, 0.2):
    cache = {}
    def chiq(qx, qy):
        k = (round(qx,6), round(qy,6))
        if k not in cache: cache[k] = dress(chi0((qx,qy), omega))
        return cache[k]
    tests = {"uniform (+1)": f_uni, "correlation_flip": f_flip,
             "random signs": np.sign(rng.standard_normal(N)),
             "random magnitudes": rng.random(N)*3}
    for lbl, fv in tests.items():
        neg = tot = 0
        for i in range(len(pts)):
            for j in range(i+1, len(pts)):
                M = chiq(pts[i][0]-pts[j][0], pts[i][1]-pts[j][1]).imag
                D = float(np.sum(fv**2 * pts[i][2] * np.diag(M) * pts[j][2]))
                tot += 1; neg += (D < 0)
        print(f"  omega={omega}  f_l = {lbl:18s}: D_ij < 0 in {neg}/{tot}  (theory: 0)")
