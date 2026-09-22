"""
team-lead, 2026-09-22, from scratch: WHERE does the omega-integration sit in the
multiband chain, and what does that mean for testing Theorem (U) condition (ii)?

My own note ceiling/unsplit_LZ_route_relaxes_eq37_2026-09-22.md section 4 states:

  "(ii) requires the per-omega (per-k) Lambda(omega), while the literature prints
   the omega-integrated lambda_ij = 2 int Lambda_ij(omega)/omega domega.
   => testing (ii) with printed lambda_ij can only give a sound FAIL, never a PASS."

AUDIT: per multiband_ceiling_lift_2026-09-18.md L.44 the object entering (ML1) is

  lambda_ij(k) = int 2 omega a2F_ij(omega) / (omega^2 + (2 pi T k)^2) domega

i.e. the omega-integration is ALREADY INSIDE Lambda-tilde(k); k is the Matsubara
index, not a frequency. So condition (ii) is a statement about the INTEGRATED
lambda_ij(k) for each k -- NOT about a2F_ij(omega) being pointwise nonneg.
Pointwise nonneg is (ML3), a SUFFICIENT condition, not the requirement.

Four things checked here, all from scratch, no repo imports:
 (A) k=0 identification: lambda_ij(0) == 2 int a2F_ij/omega domega == the
     quantity the literature prints as "lambda_ij".
 (B) large-k asymptotics: lambda_ij(k) * (2 pi T k)^2 -> M2_ij = int 2 omega a2F.
     So the k->infinity end of (ii) is governed by M2_ij, also a printed-ish
     quantity (M2 = lambda <omega^2>).
 (C) THE HONEST LIMIT: exhibit a2F with lambda^Z_ij(0) > 0 but lambda^Z_ij(k) < 0
     for some k >= 1. If this exists, a PASS at k=0 does NOT imply (ii).
 (D) conversely, is a FAIL at k=0 sound? (k=0 is in the (ML1) Laplacian sum,
     L.70: "- Sum_{k>=0} Z_k (x) L(k)", so yes k=0 is required -- check the
     index range claim is what the note says.)
"""
import numpy as np

TWO_PI = 2.0 * np.pi


def lam_of_k(w, a2f, T, k):
    """lambda(k) = int 2 w a2F(w) / (w^2 + (2 pi T k)^2) dw, trapezoid on given grid."""
    a = TWO_PI * T * k
    return np.trapezoid(2.0 * w * a2f / (w ** 2 + a ** 2), w)


def lam0_direct(w, a2f):
    """the printed definition: lambda = 2 int a2F(w)/w dw"""
    return np.trapezoid(2.0 * a2f / w, w)


def M2_direct(w, a2f):
    """M2 = int 2 w a2F dw   (repo definition, multiband_ceiling_lift L.95)"""
    return np.trapezoid(2.0 * w * a2f, w)


# frequency grid (meV-ish units); a2F compactly supported, so no tail issues
w = np.linspace(1e-3, 100.0, 400001)
T = 0.5


def gauss(c, s, A):
    return A * np.exp(-0.5 * ((w - c) / s) ** 2)


print("=== (A) k=0 identification ===")
for name, a2f in [("single peak", gauss(20, 5, 1.0)),
                  ("two peaks", gauss(10, 3, 0.7) + gauss(60, 8, 0.4)),
                  ("sign-changing", gauss(12, 3, 1.0) - gauss(55, 6, 0.55))]:
    l0 = lam_of_k(w, a2f, T, 0)
    ld = lam0_direct(w, a2f)
    print(f"  {name:15s}  lambda(k=0) = {l0: .10f}   2*int a2F/w = {ld: .10f}"
          f"   rel.diff = {abs(l0-ld)/max(abs(ld),1e-30):.2e}")

print()
print("=== (B) large-k asymptotics:  lambda(k)*(2 pi T k)^2  ->  M2 ===")
a2f = gauss(10, 3, 0.7) + gauss(60, 8, 0.4)
M2 = M2_direct(w, a2f)
print(f"  M2 = int 2 w a2F dw = {M2:.8f}")
for k in (10, 100, 1000, 10000, 100000):
    a = TWO_PI * T * k
    print(f"    k={k:7d}   lambda(k)*a^2 = {lam_of_k(w, a2f, T, k)*a*a:.8f}"
          f"   ratio to M2 = {lam_of_k(w, a2f, T, k)*a*a/M2:.8f}")

print()
print("=== (C) does lambda^Z(0) >= 0 imply lambda^Z(k) >= 0 for all k? ===")
# Construct a Z-channel off-diagonal spectral function that is POSITIVE at low
# omega and NEGATIVE at high omega. Low omega is weighted by 1/omega at k=0
# (favouring the positive lobe) but the negative high-omega lobe dominates the
# 2*omega/(omega^2+a^2) weight once a is large.
# Pick the member of the family with the LARGEST RELATIVE violation, so the
# conclusion does not rest on a value close to quadrature noise.
best = None
ks = np.unique(np.round(np.logspace(0, 5, 120)).astype(int))
for A_neg in np.linspace(0.02, 2.0, 400):
    a2f_Z = gauss(6.0, 2.0, 1.0) - gauss(70.0, 8.0, A_neg)
    l0 = lam_of_k(w, a2f_Z, T, 0)
    if l0 <= 0:
        continue
    vals = np.array([lam_of_k(w, a2f_Z, T, int(k)) for k in ks])
    if vals.min() < 0:
        rel = -vals.min() / l0
        if best is None or rel > best[0]:
            best = (rel, A_neg, l0, int(ks[np.argmin(vals)]), vals.min(),
                    M2_direct(w, a2f_Z))

if best:
    rel, A_neg, l0, kbad, lbad, M2Z = best
    print(f"  FOUND: a2F^Z = G(6,2,1.0) - G(70,8,{A_neg:.4f})")
    print(f"    lambda^Z(k=0)    = {l0:+.8f}   <- PRINTED quantity, PASSES (ii) at k=0")
    print(f"    lambda^Z(k={kbad})  = {lbad:+.8f}   <- VIOLATES (ii) at k={kbad}"
          f"   (|viol|/lambda(0) = {rel:.4f})")
    print(f"    M2^Z             = {M2Z:+.8f}   <- large-k end, sign matches the violation")
    # quadrature soundness: is the violation stable under grid refinement?
    print("    grid-refinement check on lambda^Z(k=%d):" % kbad)
    for npts in (50001, 100001, 200001, 400001, 800001):
        ww = np.linspace(1e-3, 100.0, npts)
        aa = (1.0 * np.exp(-0.5 * ((ww - 6.0) / 2.0) ** 2)
              - A_neg * np.exp(-0.5 * ((ww - 70.0) / 8.0) ** 2))
        print(f"      npts={npts:7d}  lambda = {lam_of_k(ww, aa, T, kbad):+.10f}")
    print("  => a PASS at k=0 does NOT imply condition (ii). The interior/large-k")
    print("     members are genuinely independent. (Honest limit of the finding.)")
else:
    print("  no example found in the scanned family -- claim (C) NOT established")

print()
print("=== (D) is a FAIL at k=0 sound? (k=0 must be in the required index set) ===")
a2f_Zneg = gauss(6.0, 2.0, 0.30) - gauss(70.0, 8.0, 1.6)
l0n = lam_of_k(w, a2f_Zneg, T, 0)
assert l0n < 0, f"example is not negative: {l0n}"
print(f"  example a2F^Z = G(6,2,0.30) - G(70,8,1.6):  lambda^Z(0) = {l0n:+.8f} < 0")
print("  (ML1) L.70 sums the Laplacian over k>=0 (Z_0 = diag(u_n^2)), so L(0) IS")
print("  required => lambda^Z_ij(0) < 0 breaks (ii) at k=0 => sound FAIL. ")
print("  Note the PAIR term is Sum_{k>=1}; k=0 enters ONLY via the Laplacian.")

print()
print("=== SUMMARY ===")
print("  * omega-integration sits INSIDE Lambda-tilde(k), BEFORE (ML1).")
print("  * (ii) is a condition on the INTEGRATED lambda_ij(k), one per Matsubara k,")
print("    NOT on a2F_ij(omega) pointwise. Pointwise nonneg = (ML3) = SUFFICIENT only.")
print("  * printed lambda_ij IS the k=0 member exactly => printed data decides the")
print("    k=0 instance OUTRIGHT (both directions), not FAIL-only.")
print("  * but (C) shows k=0 PASS does not extend to k>=1 => the note's overall")
print("    'no PASS from printed data alone' conclusion SURVIVES, for a different")
print("    reason than it gave.")
