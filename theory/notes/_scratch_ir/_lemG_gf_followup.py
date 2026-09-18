"""[numerical-evidence] Follow-ups requested / implied by team-lead's review round 2.

A. Lemma g monotonicity at high precision (team-lead warns float64 g(varpi) suffers
   catastrophic cancellation as varpi -> 0, so a float64 "monotone: True" is luck).
   Reported honestly: WHERE float64 actually breaks, and the dps=60 verdict.

B. Budget for gap 1 (cone: top eigenvector of O_pair stays in C) AFTER GF.
   Sharp (cone) route :  Tc <= C_inf sqrt(rho(M2^Z))
   Unconditional GF   :  Tc <= C_inf sqrt(rho(M2''_Z)),
       M2''_Z = diag(M2^ph + c_GF M2^sf) + c_GF * offdiag(M2^ph + M2^sf)
   ==> gap 1 buys exactly sqrt(rho(M2''_Z)/rho(M2^Z)) in [1, sqrt(c_GF)].
   This is the budget-first check CLAUDE.md requires BEFORE taking the assignment.

C. Signed-kernel interface (Z-dom-k): does |Lambda^Delta| stay nonincreasing when
   omega_ph != omega_sf?  (team-lead measured 205/800; independent re-measure here.)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import mpmath as mp

G2_LO, G2_HI = 1.3181093293757127, 1.3181615784372800
G2 = 1.318140212132
PI26 = np.pi ** 2 / 6
C_GF = PI26 / G2_LO                 # 1.247950, the certified-side worst value
rng = np.random.default_rng(20260919)

print("== A. lemma g: float64 vs dps=60 ==")


def g_f64(v):
    return (1 - np.pi * v / np.sinh(np.pi * v)) / v ** 2


def g_mp(v, dps=60):
    with mp.workdps(dps):
        v = mp.mpf(v)
        return (1 - mp.pi * v / mp.sinh(mp.pi * v)) / v ** 2


def g_series_mp(v, dps=60, M=4000):
    with mp.workdps(dps):
        v = mp.mpf(v)
        return 2 * mp.nsum(lambda n: (-1) ** (n + 1) / (v * v + n * n), [1, mp.inf])


for v in (1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 0.4589, 1.0, 5.0):
    a, b, c = g_f64(v), g_mp(v), g_series_mp(v)
    print(f"  v={v:9.2e}  f64={a!r:24s} dps60={mp.nstr(b,12):14s} series={mp.nstr(c,12):14s}"
          f"  |f64-dps60|/dps60 = {float(abs(a - b) / b):.2e}")

for tag, lo, hi, npts in (("log grid 1e-8..1", -8, 0, 4001), ("log grid 1e-3..1", -3, 0, 4001)):
    vs = 10 ** np.linspace(lo, hi, npts)
    mono_f64 = bool(np.all(np.diff(g_f64(vs)) < 0))
    print(f"  {tag}: float64 strictly decreasing? {mono_f64}")
vs = [mp.mpf(10) ** mp.mpf(e) for e in np.linspace(-8, 2, 601)]
gs = [g_mp(v) for v in vs]
print(f"  dps=60 log grid 1e-8..1e2 (601 pts): strictly decreasing? "
      f"{all(gs[i+1] < gs[i] for i in range(len(gs)-1))};  g(1e-8) = {mp.nstr(gs[0], 15)}"
      f"   pi^2/6 = {mp.nstr(mp.pi**2/6, 15)}")
print("  (monotonicity itself is proved analytically in note section 2; this is only a consistency gate)")

print("\n== B. budget for gap 1 (cone) after GF ==")


def rho(X):
    return float(np.abs(np.linalg.eigvals(X)).max())


def offd(X):
    return X - np.diag(np.diag(X))


def gap1_budget(M2ph, M2sf):
    """ratio of the unconditional GF ceiling to the sharp (cone) ceiling."""
    M2Z = M2ph + M2sf
    M2pp = np.diag(np.diag(M2ph + C_GF * M2sf)) + C_GF * offd(M2Z)
    return np.sqrt(rho(M2pp) / rho(M2Z)), rho(M2Z), rho(M2pp)


cases = {
    "MgB2 (Golubov, conv. B), pure phonon":
        (np.array([[1.017, np.sqrt(0.213 * 0.155)], [np.sqrt(0.213 * 0.155), 0.448]]), np.zeros((2, 2))),
    "MgB2 (conv. A), pure phonon":
        (np.array([[1.017, 0.184], [0.184, 0.448]]), np.zeros((2, 2))),
    "all-equal phonon [[1,1],[1,1]]": (np.ones((2, 2)), np.zeros((2, 2))),
    "pure off-diagonal phonon":       (np.array([[0., 1.], [1., 0.]]), np.zeros((2, 2))),
    "pure spin fluctuation (worst)":   (np.zeros((2, 2)), np.ones((2, 2))),
    "LiFeAs-like ph+sf mix":
        (np.array([[0.3, 0.2], [0.2, 0.25]]), np.array([[0.5, 0.9], [0.9, 0.6]])),
}
for tag, (a, b) in cases.items():
    f, r0, r1 = gap1_budget(a, b)
    print(f"  {tag:40s} rho(M2^Z)={r0:.4f} rho(M2''_Z)={r1:.4f}  gap-1 buys x{f:.4f}  ({100*(f-1):.2f}%)")
worst = -np.inf
for _ in range(20000):
    B = int(rng.integers(2, 6))
    a = rng.uniform(0, 1, (B, B)) * (rng.random((B, B)) < 0.8); a = (a + a.T) / 2
    b = rng.uniform(0, 1, (B, B)) * (rng.random((B, B)) < 0.8); b = (b + b.T) / 2
    if rho(a + b) < 1e-9:
        continue
    worst = max(worst, gap1_budget(a, b)[0])
print(f"  20000 random (ph, sf) pairs, B=2..5: worst factor {worst:.6f}"
      f"   (analytic cap sqrt(c_GF) = {np.sqrt(C_GF):.6f})")

print("\n== C. is |Lambda^Delta| nonincreasing when omega_ph != omega_sf? ==")


def lam(k, w, mass):
    return mass * w * w / (k * k + w * w)


bad = 0; worstjump = 0.0
for _ in range(800):
    wph, wsf = 10 ** rng.uniform(-1, 1, 2)
    mph, msf = rng.uniform(0.05, 2.0, 2)
    k = np.arange(1, 400)
    d = np.abs(lam(k, wph, mph) - lam(k, wsf, msf))
    dd = np.diff(d)
    if dd.max() > 1e-14:
        bad += 1
        worstjump = max(worstjump, dd.max())
print(f"  |Lambda^ph - Lambda^sf| fails nonincreasing in {bad}/800 two-frequency draws"
      f"   (worst positive increment {worstjump:.3e})")
print("  ==> the (G'-abs) step ||O[L]|| = ||O[|L|]|| + Lemma M / (G) on |L| needs |L| nonincreasing;")
print("      the triangle-inequality route (Z-dom-k) does not, since it acts on the two nonneg pieces.")
