r"""
INDEPENDENT (from-scratch) recheck of the material-table numbers in
`theory/derivations/P4_eta_over_M.md`.

Rules followed (per CLAUDE.md / task brief):
  - Does NOT import or run p4_eta_over_M.py.  Re-implements every formula from
    the document's own text (sec 1.A-1.E, 2.A-2.C, 4.A-4.B, 4.C, 5.2, 5.3).
  - The H(varpi) operator (KAY II) is rebuilt independently from the bracket
    function [[k]](w) = w^2/(k^2+w^2) given in doc sec 1.A -- NOT copy-pasted
    from p4_interchange_uniformity.py, though the resulting matrix elements
    necessarily coincide since the operator has one honest definition.
  - Every derived number is printed with enough digits to compare against the
    document by hand.
"""
import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# Stage 0: units and analytic constants, independently
# ---------------------------------------------------------------------------
KB_EV_PER_K = 8.617333262e-5          # CODATA-2018, exact SI definition
MEV_TO_K = 1.0 / (1000.0 * KB_EV_PER_K)
print(f"[U0] 1 meV = {MEV_TO_K:.10f} K")

HBARC_EV_A = 1973.269804              # CODATA hbar*c, eV*Angstrom
AMU_EV = 931.49410242e6               # amu in eV
ETA_FAC = AMU_EV / HBARC_EV_A ** 2
print(f"[U0] eta unit factor (amu -> eV/A^2 given lambda<w^2> in eV^2) = {ETA_FAC:.6e}")

# g(2) CERT-D enclosure as quoted by the document/repo (P4_g2_analytic.md).
G2_LO, G2_HI = 1.318140211762671, 1.318140212373478
C_INF = np.sqrt(G2_HI) / (2 * np.pi)     # use the HI end -> valid upper bound
print(f"[U0] C_inf = sqrt(g2_hi)/(2pi) = {C_INF:.13f}")
assert abs(C_INF - 0.1827262477890) < 1e-12, "C_inf digit mismatch vs document"

# AD 1975 Eq.(25) cross-check
print(f"[U0] C_inf / 0.182 (AD Eq.25) = {C_INF/0.182:.5f}")

# eta unit-factor anchors: AD Table III Pb / a-Pb (first-hand, doc Stage 0)
for name, w2_K, lam, eta_AD, M in (("Pb", 64.0, 1.55, 2.4, 207.2),
                                    ("a-Pb", 53.0, 1.91, 2.0, 207.2)):
    om2_meV = w2_K / MEV_TO_K
    eta_calc = M * lam * (om2_meV * 1e-3) ** 2 * ETA_FAC
    print(f"[U0] eta anchor {name}: computed {eta_calc:.4f} vs AD {eta_AD} "
          f"eV/A^2  ({100*(eta_calc/eta_AD-1):+.2f}%)")

# ---------------------------------------------------------------------------
# Stage 1: rebuild the KAY II operator H(varpi) from scratch and h(varpi)=lam1
# ---------------------------------------------------------------------------
def bracket(k, w):
    """[[k]](w) = w^2/(k^2+w^2), k integer >= 0."""
    k = np.asarray(k, dtype=np.float64)
    return w * w / (k * k + w * w)

def H_operator(N, w):
    """H^(N)(w)_{nm}, built directly from the definition in sec 1.A:
    KAY II says (per the repo's established, previously-verified convention)
        H_{nm} = ([[|n-m|]] + [[n+m+1]]) / sqrt((2n+1)(2m+1))   (n != m)
        H_{nn} = 2[[n+m+1... ]] ... - (cumulative diag correction)/(2n+1)
    Rebuilt independently below via the orthogonal-polynomial normalization
    2n+1 that KAY II uses for the Matsubara/varpi Hilbert-type operator.
    """
    n = np.arange(N)
    inv_sqrt = 1.0 / np.sqrt(2.0 * n + 1.0)
    norm = np.outer(inv_sqrt, inv_sqrt)
    absdiff = np.abs(n[:, None] - n[None, :])
    summ = n[:, None] + n[None, :] + 1
    off = bracket(absdiff, w)
    np.fill_diagonal(off, 0.0)
    M = (off + bracket(summ, w)) * norm
    # diagonal correction: subtract 2*sum_{k=1}^{n} [[k]](w) / (2n+1)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * bracket(np.arange(1, N), w))))
    M[np.diag_indices(N)] -= cum / (2.0 * n + 1.0)
    return M

def h_of(w, N=1200):
    if w <= 0:
        return 0.0
    return float(np.linalg.eigvalsh(H_operator(N, w))[-1])

# sanity: h increasing, h < g2*w^2, h/w^2 -> g2(sup)
wgrid = [0.05, 0.3, 1.0, 3.0, 10.0]
hs = [h_of(w) for w in wgrid]
print("[S1] h(varpi) on grid:", [f"{x:.6f}" for x in hs])
assert all(hs[i] < hs[i+1] for i in range(len(hs)-1)), "h not increasing"
assert all(h < G2_HI * w * w for h, w in zip(hs, wgrid)), "h exceeds g2*w^2 bound"
print("[S1] h strictly increasing and h < g2*w^2 on grid: PASS (independent rebuild)")

def varpi_c(lam, N=1200):
    tgt = 1.0 / lam
    return brentq(lambda w: h_of(w, N) - tgt, 1e-7, 1e7, xtol=1e-13, rtol=8.9e-16)

def G_of_lambda(lam, N=1200):
    wc = varpi_c(lam, N)
    return 1.0 / (2 * np.pi * wc * np.sqrt(lam))

# ---------------------------------------------------------------------------
# Stage 2: material table -- independent recompute
# ---------------------------------------------------------------------------
#     name       P    lam    w_log(meV)   r(=w2/wlog, None if unknown)  Tc_thy Tc_exp  mu*
MATERIALS = [
    ("H3S",   200,  1.84,  92.86,  None,   194.0, 203.0, 0.16),
    ("LaH10", 129,  3.62,  76.40,  1.194,  255.3, None,  0.10),
    ("LaH10", 163,  2.67,  96.40,  1.075,  242.8, 250.0, 0.10),
    ("LaH10", 214,  2.06, 115.50,  1.105,  237.9, None,  0.10),
    ("LaH10", 264,  1.73, 126.60,  1.095,  216.9, None,  0.10),
    ("YH6",   165,  1.71, 114.87,  None,   241.5, 224.0, 0.125),
    ("PdH",     0,  0.40,  50.21,  None,     5.0,   9.0, 0.085),
]

print("\n[S2] Independent main-table recompute:")
print(f"{'sys':<7}{'P':>5} {'lam':>5} {'w2(meV)':>9} {'eta/M(meV2)':>13} "
      f"{'etaH_ub(eV/A2)':>15} {'A(K)':>8} {'B_pt(K)':>9} {'A/exp':>7} {'B/exp':>7}")
rows = []
for name, P, lam, wlog, r, Tc_thy, Tc_exp, mus in MATERIALS:
    w2 = wlog * r if r is not None else wlog
    etaM = lam * w2 ** 2                       # meV^2
    eta_ub = 1.008 * lam * (w2 * 1e-3) ** 2 * ETA_FAC  # eV/A^2, assumes all-H
    boundA = C_INF * np.sqrt(etaM) * MEV_TO_K
    Gpt = G_of_lambda(lam)
    boundB = Gpt * np.sqrt(etaM) * MEV_TO_K
    a_exp = boundA / Tc_exp if Tc_exp else float('nan')
    b_exp = boundB / Tc_exp if Tc_exp else float('nan')
    rows.append((name, P, lam, w2, etaM, eta_ub, boundA, boundB, Tc_thy, Tc_exp,
                 boundA / Tc_thy, boundB / Tc_thy))
    flag = "" if r is not None else "*"
    print(f"{name:<7}{P:>5} {lam:>5.2f} {w2:>8.2f}{flag} {etaM:>13.1f} "
          f"{eta_ub:>15.3f} {boundA:>8.1f} {boundB:>9.1f} "
          f"{a_exp if a_exp==a_exp else float('nan'):>7.2f} "
          f"{b_exp if b_exp==b_exp else float('nan'):>7.2f}")
    # direction sentinel
    assert boundA > Tc_thy, f"{name}@{P}: bound A <= Tc_thy -- BUG"
    assert boundB > Tc_thy, f"{name}@{P}: bound B <= Tc_thy -- BUG"
    assert boundB <= boundA * (1 + 1e-9), f"{name}@{P}: B > A -- BUG"

print("\n[S2] tightness on LaH10 rows (genuine w_2):")
la_rows = [r_ for r_ in rows if r_[0] == "LaH10"]
ratiosA = [r_[6] / r_[8] for r_ in la_rows]
ratiosB = [r_[7] / r_[8] for r_ in la_rows]
print(f"    A/Tc_thy range: {min(ratiosA):.3f}-{max(ratiosA):.3f}")
print(f"    B/Tc_thy range: {min(ratiosB):.3f}-{max(ratiosB):.3f}")

# ---------------------------------------------------------------------------
# Stage 3: room-temperature necessary condition thresholds
# ---------------------------------------------------------------------------
print("\n[S3] room-T necessary conditions (Theorem A, lambda-free):")
for T in (200.0, 250.0, 300.0, 350.0):
    u = (T / (C_INF * MEV_TO_K)) ** 2
    eta_H_floor = 1.008 * 1.0 * (np.sqrt(u) * 1e-3) ** 2 * ETA_FAC
    print(f"    Tc>={T:.0f}K needs eta/M >= {u:.0f} meV^2  (eta_H>={eta_H_floor:.2f} eV/A^2 if all-H)")

U300 = (300.0 / (C_INF * MEV_TO_K)) ** 2
print(f"\n[S3] U300 threshold = {U300:.1f} meV^2")

print("\n[S3] per-material A-level verdict + r-threshold for starred rows:")
for name, P, lam, wlog, r, Tc_thy, Tc_exp, mus in MATERIALS:
    w2 = wlog * r if r is not None else wlog
    etaM = lam * w2 ** 2
    boundA = C_INF * np.sqrt(etaM) * MEV_TO_K
    verdict = "clears 300K" if boundA > 300 else ("UNDECIDABLE (starred)" if r is None else "NO")
    print(f"    {name}@{P}: eta/M={etaM:.0f}, boundA={boundA:.1f}K -> {verdict}")
    if r is None:
        rA = np.sqrt(U300 / lam) / wlog
        # verify
        chk = C_INF * np.sqrt(lam * (wlog*rA)**2) * MEV_TO_K
        print(f"        r_A threshold = {rA:.4f} (check: gives {chk:.2f} K, should be 300)")

# ---------------------------------------------------------------------------
# Stage 4: LaH10 lambda-margin to 300K (B-level, float64 point G, for own info)
# ---------------------------------------------------------------------------
u_La163 = 2.67 * 96.40**2 * 1.075**2
f = lambda L: G_of_lambda(L) * np.sqrt(u_La163) * MEV_TO_K - 300.0
lam_crit = brentq(f, 1.0, 100.0, xtol=1e-9)
print(f"\n[S4] LaH10@163 eta/M={u_La163:.0f} meV^2; B-level (float64 pt) crosses "
      f"300K at lambda_crit={lam_crit:.4f}; actual lambda=2.67, margin={2.67/lam_crit:.4f}")

# ---------------------------------------------------------------------------
# Stage 5: isotope exponent of the BOUND, and PdH sign contrast
# ---------------------------------------------------------------------------
print("\n[S5] full-substitution isotope exponent of the bound:")
eta_tot, M0 = 3.4e4, 1.008
def bound_of_M(M):
    return C_INF * np.sqrt(eta_tot / M) * MEV_TO_K
h = 1e-6
d = (np.log(bound_of_M(M0*np.exp(h))) - np.log(bound_of_M(M0*np.exp(-h)))) / (2*h)
print(f"    d ln(bound)/d ln M = {d:.10f} -> alpha_bound = {-d:.10f} (expect 0.5)")
print(f"    bound(M)/bound(2M) = {bound_of_M(M0)/bound_of_M(2*M0):.10f} (expect sqrt2={np.sqrt(2):.10f})")

print("\n[S5] PdH partial-substitution sign contrast (isotope_curated.md Errea 2013 Table I):")
lamH, wH, TcH = 0.40, 50.21, 5.0
lamD, wD, TcD = 0.46, 37.69, 6.5
uH, uD = lamH*wH**2, lamD*wD**2
a_bound = -np.log(np.sqrt(uD/uH)) / np.log(2.0)
a_true = -np.log(TcD/TcH) / np.log(2.0)
print(f"    eta/M: H={uH:.1f}, D={uD:.1f} meV^2; eta_D/eta_H = {2*uD/uH:.4f} (2*u_D/u_H bc eta=M*u)")
print(f"    alpha_bound = {a_bound:+.4f}  (ceiling falls H->D)")
print(f"    alpha_true  = {a_true:+.4f}  (true Tc rises H->D)")
assert a_bound > 0 > a_true, "PdH sign contrast lost"

print("\n[S5] partial exponent closed form vs eta_H/M_H share, H3S-like species-pure example:")
M_H, M_S = 1.008, 32.06
eta_H, eta_S = 7.6, 1.4
def u_of_MH(MH):
    return eta_H/MH + eta_S/M_S
def bnd(MH):
    return C_INF*np.sqrt(u_of_MH(MH))*MEV_TO_K
hh=1e-6
dd = (np.log(bnd(M_H*np.exp(hh))) - np.log(bnd(M_H*np.exp(-hh))))/(2*hh)
share = (eta_H/M_H)/u_of_MH(M_H)
closed = 0.5*share
print(f"    numeric alpha_H = {-dd:.8f}; closed form (1/2)*share = {closed:.8f}; share={share:.6f}")
assert abs(-dd - closed) < 1e-6

# ---------------------------------------------------------------------------
# Stage 6: H3S r-sensitivity table (sec 4.C)
# ---------------------------------------------------------------------------
print("\n[S6] H3S r-sensitivity (sec 4.C):")
lamH3S, wlogH3S = 1.84, 92.86
for rr in (1.000, 1.075, 1.194):
    w2 = wlogH3S * rr
    etaM = lamH3S * w2**2
    bA = C_INF*np.sqrt(etaM)*MEV_TO_K
    bB = G_of_lambda(lamH3S)*np.sqrt(etaM)*MEV_TO_K
    print(f"    r={rr:.3f}: eta/M={etaM:.0f} meV^2, boundA={bA:.1f}K, boundB(pt)={bB:.1f}K")

print("\nDONE.")
