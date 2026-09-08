r"""
P3 Task A (implementer numerics, 2026-09-09) — closed-form R(lambda) checks + LaH10
decomposition reproduction, per P3_tc_max_pressure.md Sec.3/Sec.5.1/Sec.5.2.

Allen-Dynes closed form, Einstein spectrum (f2=1), Xie 2022 transcription
(P3_tc_max_pressure.md Sec.3, eq (3)):

    Tc_AD = (f1 * omega_log / 1.20) * exp[ -1.04(1+lam) / (lam - mustar*(1+0.62*lam)) ]
    f1 = [ 1 + (lam / (2.46*(1+3.8*mustar)))**1.5 ]**(1/3)
    R_AD(lam, mustar) = Tc_AD / (0.182 * omega_log * sqrt(lam))
                       = f1 / (1.20 * 0.182 * sqrt(lam)) * exp[ -1.04(1+lam)/(lam-mustar*(1+0.62*lam)) ]

Note omega_log cancels out of R_AD -- it is a pure function of (lam, mustar).

Run: python theory/derivations/p3_task_A_R_lambda.py   (from repo root)
"""
import math

MU_STAR = 0.1


def f1(lam, mustar):
    return (1 + (lam / (2.46 * (1 + 3.8 * mustar))) ** 1.5) ** (1 / 3)


def R_AD(lam, mustar, c=0.182):
    expo = -1.04 * (1 + lam) / (lam - mustar * (1 + 0.62 * lam))
    return f1(lam, mustar) / (1.20 * c * math.sqrt(lam)) * math.exp(expo)


def asymptotic_coeff(mustar):
    """Exact lambda->inf limit of Tc_AD/(omega_log*sqrt(lam)), i.e. the true
    strong-coupling coefficient (should reduce to 0.1878 at mustar=0)."""
    return (1 / 1.20) / math.sqrt(2.46 * (1 + 3.8 * mustar)) * math.exp(-1.04 / (1 - 0.62 * mustar))


def loglog_slope(xs, ys):
    """OLS slope of ln(ys) vs ln(xs), all four points (not consecutive differences)."""
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    n = len(xs)
    xbar = sum(lx) / n
    ybar = sum(ly) / n
    num = sum((x - xbar) * (y - ybar) for x, y in zip(lx, ly))
    den = sum((x - xbar) ** 2 for x in lx)
    slope = num / den
    intercept = ybar - slope * xbar
    return slope, intercept


def endpoint_slope(xs, ys):
    """Simple two-point secant slope using only the first and last point
    (ln(y_last/y_first) / ln(x_last/x_first)) -- ignores the middle points entirely."""
    return math.log(ys[-1] / ys[0]) / math.log(xs[-1] / xs[0])


print("=" * 70)
print("Step 1: R_AD(lambda, mu*=0.1) for lambda in [1.5,4] -- monotonicity")
print("=" * 70)
lams = [1.5 + 0.1 * k for k in range(26)]  # 1.5 .. 4.0 step 0.1
Rs = [R_AD(l, MU_STAR) for l in lams]
diffs = [Rs[i + 1] - Rs[i] for i in range(len(Rs) - 1)]
print(f"{'lambda':>8}{'R_AD':>10}")
for l, r in zip(lams, Rs):
    print(f"{l:>8.2f}{r:>10.5f}")
print(f"\nmonotonic increasing (all diffs>0)? {all(d > 0 for d in diffs)}")
print(f"min diff = {min(diffs):.6f}  max diff = {max(diffs):.6f}")

print()
print("=" * 70)
print("Step 2: lambda -> infinity limit of R_AD")
print("=" * 70)
for mustar in (0.0, 0.1):
    for lam_big in (1e3, 1e4, 1e5, 1e6, 1e8):
        r = R_AD(lam_big, mustar)
        print(f"  mu*={mustar:.2f}  lambda={lam_big:>10.0e}  R_AD={r:.6f}")
    A = asymptotic_coeff(mustar)
    print(f"  mu*={mustar:.2f}  analytic asymptotic coefficient A(mu*) = {A:.6f}"
          f"   ->  A/0.182 = {A/0.182:.6f}   (A/0.1878 = {A/0.1878:.6f})")
    print()

print("Document claim: 'lambda->inf (take lambda=1e5) gives R->1.03, because 0.182 vs exact 0.1878'.")
print("This statement is explicitly conditioned on mu*=0 in the document (Sec.3 line: '...~1.03 (mu*=0)').")
print(f"  Numerically confirmed at mu*=0:  R_AD(1e5,0)={R_AD(1e5,0.0):.6f}  (matches ~1.03: "
      f"{'YES' if abs(R_AD(1e5,0.0)-1.03) < 0.01 else 'NO'})")
print(f"  At the task-specified mu*=0.1:   R_AD(1e5,0.1)={R_AD(1e5,0.1):.6f}  "
      f"({'NOT close to 1.03 -- mu*=0 is required for the 1.03 number, this is NOT a mismatch of the doc,'
      f' the doc explicitly scopes it to mu*=0' if abs(R_AD(1e5,0.1)-1.03) > 0.01 else 'matches 1.03'})")

print()
print("=" * 70)
print("Step 3 (decisive test): plug LaH10's 4 lambda values into R_AD closed form,")
print("mu*=0.1 fixed, f2=1 (Einstein). Compare dlnR_AD/dlnP to measured dlnR/dlnP=-0.417.")
print("=" * 70)
# LaH10 four points (P3_tc_max_pressure.md Sec.2 table)
P    = [129, 163, 214, 264]           # GPa
wlog = [886.6, 1118.7, 1340.3, 1469.1]  # K
lam  = [3.62, 2.67, 2.06, 1.73]
TcME = [255.3, 242.8, 237.9, 216.9]   # K

R_AD_vals = [R_AD(l, MU_STAR) for l in lam]
print(f"{'P(GPa)':>8}{'lambda':>8}{'R_AD(lam,mu*=0.1)':>20}")
for p, l, r in zip(P, lam, R_AD_vals):
    print(f"{p:>8}{l:>8.2f}{r:>20.5f}")

slope_RAD, _ = loglog_slope(P, R_AD_vals)
print(f"\nd ln R_AD / d ln P (log-log fit, all 4 points) = {slope_RAD:.4f}")
print(f"Document's closed-form check-4 value: -0.253")
print(f"Measured (real) dlnR/dlnP = -0.417")
print(f"Fraction explained by lambda alone = slope_RAD / measured = {slope_RAD/-0.417:.4f}"
      f"  (document claims ~60.7% = 0.253/0.417 = {0.253/0.417:.4f})")

print()
print("=" * 70)
print("Step 4: reproduce Sec.2 decomposition table exactly (Tc^emp, R, log-log slopes)")
print("=" * 70)
Tcemp = [0.182 * w * math.sqrt(l) for w, l in zip(wlog, lam)]
Rtab  = [tc / te for tc, te in zip(TcME, Tcemp)]

print(f"{'P':>6}{'wlog':>9}{'lambda':>8}{'Tc^emp':>10}{'Tc^ME':>9}{'R':>8}")
for p, w, l, te, tm, r in zip(P, wlog, lam, Tcemp, TcME, Rtab):
    print(f"{p:>6}{w:>9.1f}{l:>8.2f}{te:>10.1f}{tm:>9.1f}{r:>8.3f}")

n_fit, _        = loglog_slope(P, wlog)
m_fit, _        = loglog_slope(P, lam)
dlnR_fit, _     = loglog_slope(P, Rtab)
dlnTcemp_fit, _ = loglog_slope(P, Tcemp)
dlnTc_fit, _    = loglog_slope(P, TcME)

n_ep        = endpoint_slope(P, wlog)
m_ep        = endpoint_slope(P, lam)
dlnR_ep     = endpoint_slope(P, Rtab)
dlnTcemp_ep = endpoint_slope(P, Tcemp)
dlnTc_ep    = endpoint_slope(P, TcME)

print()
print(f"{'quantity':<24}{'4-pt OLS fit':>14}{'endpoint secant':>17}{'document':>10}")
print(f"{'n = dln(wlog)/dlnP':<24}{n_fit:>14.4f}{n_ep:>17.4f}{0.70:>10.2f}")
print(f"{'m = dln(lam)/dlnP':<24}{m_fit:>14.4f}{m_ep:>17.4f}{-1.03:>10.2f}")
print(f"{'dlnR/dlnP':<24}{dlnR_fit:>14.4f}{dlnR_ep:>17.4f}{-0.42:>10.2f}")
print(f"{'dlnTc^emp/dlnP':<24}{dlnTcemp_fit:>14.4f}{dlnTcemp_ep:>17.4f}{0.19:>10.2f}")
print(f"{'dlnTc/dlnP':<24}{dlnTc_fit:>14.4f}{dlnTc_ep:>17.4f}{-0.23:>10.2f}")

print()
print(f"Identity check (4-pt OLS): dlnTc^emp/dlnP + dlnR/dlnP = {dlnTcemp_fit + dlnR_fit:.4f}  "
      f"vs  dlnTc/dlnP = {dlnTc_fit:.4f}   (should match to float precision, OLS is linear in y)")
print(f"Identity check (endpoint):  dlnTc^emp/dlnP + dlnR/dlnP = {dlnTcemp_ep + dlnR_ep:.4f}  "
      f"vs  dlnTc/dlnP = {dlnTc_ep:.4f}   (exact identity, holds for ANY consistent slope definition"
      f" since Tc=Tc^emp*R pointwise)")
print(f"Also cross-check via n + 0.5*m (4-pt OLS) = {n_fit + 0.5*m_fit:.4f}  vs  dlnTc^emp/dlnP = {dlnTcemp_fit:.4f}")
print(f"Also cross-check via n + 0.5*m (endpoint) = {n_ep + 0.5*m_ep:.4f}  vs  dlnTc^emp/dlnP = {dlnTcemp_ep:.4f}")

print()
print("HONEST FINDING: the proper 4-point OLS log-log regression does NOT reproduce the document's")
print("stated single-value slopes as precisely as a simple two-point secant (endpoint 129 vs 264 only,")
print("ignoring the middle two points) does. Compare (fit vs endpoint vs document):")
print(f"  m         : {m_fit:.3f}  vs  {m_ep:.3f}  vs  -1.03   <- endpoint matches, OLS is ~0.007 off")
print(f"  dlnR/dlnP : {dlnR_fit:.3f}  vs  {dlnR_ep:.3f}  vs  -0.42   <- endpoint matches closely, OLS is ~0.02 off")
print(f"  dlnTc/dlnP: {dlnTc_fit:.3f}  vs  {dlnTc_ep:.3f}  vs  -0.23   <- endpoint matches closely, OLS is ~0.02 off")
print("So despite the document's wording 'log-log linear fit over 4 points, not adjacent differences',")
print("the actual reported numbers are best explained as a plain first-to-last secant slope")
print("(ln(y(264)/y(129))/ln(264/129)), not a genuine 4-point least-squares regression.")
print("This is a documentation/methodology ambiguity, not a numerical bug on either side --")
print("flagging it for the main session to adjudicate/clarify in the doc.")
