r"""
P3 §4 stability-bound check on a NON-HYDRIDE spectrum: MgB2.
Independent, self-contained (imports nothing from this workspace).  mpmath dps=25.

WHAT IS TESTED
--------------
The P3 upper bound  lambda_*(sigma, gamma3) = lambda_0 * (1 + kappa/2 * sigma^2
+ c3 * gamma3 * sigma^3)  is a *narrow-spectrum expansion* of the shape-exact
functional  lambda_*[P] = 1 / max_{s>0} Int g(s*omega) P(omega) domega,
with the Semenok kernel g and P(omega) = (2/lambda) * alpha^2F(omega)/omega.

Here we (i) evaluate the shape-EXACT lambda_* by direct integration on a real
first-principles MgB2 alpha^2F, (ii) evaluate the expansion at that spectrum's
(sigma, gamma3), (iii) compare, and (iv) run two narrowing controls to show the
expansion converges onto the exact value once sigma is inside its domain.

DATA PROVENANCE  (read this before quoting any number)
------------------------------------------------------
alpha^2F(omega) embedded below = entry "822", key `Mg1B2_JVASP_19821`, of
    database.json, https://github.com/henniggroup/BETE-NET  (Apache-2.0)
    = training database of Gibson, Hire, Dee, Barrera, Geisler, Hirschfeld,
      Hennig, "Accelerating superconductor discovery through tempered deep
      learning of the electron-phonon spectral function", arXiv:2401.16611.
    Underlying calculation is high-throughput DFPT on the JARVIS structure
    JVASP-19821 (MgB2, P6/mmm).  Downloaded 2026-09-09 via
      curl -sL https://raw.githubusercontent.com/henniggroup/BETE-NET/main/database.json
    Fields used: d["Freq_meV"]["822"] (1025 pts, uniform, 0.2500..102.651 meV)
    and d["a2F"]["822"].  Values below are those a2F entries rounded to 5 dp;
    the grid is regenerated as omega_i = W0 + i*DW (max deviation from the
    stored grid 9.8e-5 meV, i.e. 1e-6 of the range).  NOT a digitized figure.

    HONEST CAVEAT ON THE SPECTRUM ITSELF: this high-throughput DFPT spectrum
    gives lambda = 0.6235, w_log = 65.4 meV, and puts its dominant E2g-derived
    peak at ~76 meV.  The canonical converged MgB2 result (Kong, Dolgov,
    Jepsen, Andersen, PRB 64, 020501(R) (2001)) is lambda ~ 0.87 with the E2g
    peak near 60-70 meV -- the high-throughput q-grid/settings underbind the
    E2g softening.  This matters for the *absolute* lambda of MgB2; it does NOT
    affect the methodology test, because lambda_* depends only on the SHAPE P
    (an overall rescaling of alpha^2F changes lambda and leaves lambda_*
    untouched).  We therefore report lambda/lambda_* for both this spectrum's
    lambda and for the literature lambda = 0.87 on the same shape.
"""
import mpmath as mp

mp.mp.dps = 25

# ---------------------------------------------------------------- kernel ----
A = lambda x: mp.im(mp.polygamma(1, 1j * x))      # Im psi'(ix)
B = lambda x: mp.re(mp.polygamma(2, 1j * x))      # Re psi''(ix)


def g(x):
    x = mp.mpf(x)
    return 6 * x**2 + 12 * x**3 * A(x) + 6 * x**4 * B(x)


def gp(x):
    """g'(x), analytic.  A'=B, B'=-Im psi'''(ix)."""
    x = mp.mpf(x)
    return (12 * x + 36 * x**2 * A(x) + 36 * x**3 * B(x)
            - 6 * x**4 * mp.im(mp.polygamma(3, 1j * x)))


X_MAX = mp.findroot(gp, mp.mpf('0.33'))
G_MAX = g(X_MAX)
LAM0 = 1 / G_MAX
KAPPA = -X_MAX**2 * mp.diff(g, X_MAX, 2) / G_MAX
C3 = -X_MAX**3 * mp.diff(g, X_MAX, 3) / (6 * G_MAX)

REF = dict(x_max='0.32732148', g_max='0.27089315', lam0='3.69149',
           kappa='5.60346', c3='0.89590')


# -------------------------------------------------- discrete measure P ------
def measure(ws, a2f):
    """Return (lambda, [(omega_i, p_i)]) with sum p_i = 1, p ~ trapezoid of a2F/omega."""
    ws = [mp.mpf(w) for w in ws]
    f = [mp.mpf(a) / w for a, w in zip(a2f, ws)]          # alpha^2F/omega
    wt = []
    n = len(ws)
    for i in range(n):
        lo = ws[i] - ws[i - 1] if i > 0 else mp.mpf(0)
        hi = ws[i + 1] - ws[i] if i < n - 1 else mp.mpf(0)
        wt.append((lo + hi) / 2)
    norm = sum(wi * fi for wi, fi in zip(wt, f))           # = lambda/2
    lam = 2 * norm
    pts = [(w, wi * fi / norm) for w, wi, fi in zip(ws, wt, f) if fi > 0]
    tot = sum(p for _, p in pts)
    pts = [(w, p / tot) for w, p in pts]
    return lam, pts


def prune(pts, tol=mp.mpf('1e-12')):
    """Drop weights below tol (pure speed; renormalize)."""
    keep = [(w, p) for w, p in pts if p > tol]
    tot = sum(p for _, p in keep)
    return [(w, p / tot) for w, p in keep]


def moments(pts):
    mu = sum(p * w for w, p in pts)
    c2 = sum(p * (w - mu)**2 for w, p in pts)
    c3m = sum(p * (w - mu)**3 for w, p in pts)
    return mu, mp.sqrt(c2) / mu, c3m / c2**mp.mpf('1.5')


def lam_star_exact(pts, mu):
    """1 / max_s G(s),  G(s) = Int g(s*omega) P(domega).  Root of G'(s)=0.

    Solved in the dimensionless variable u = s*<omega> (y = omega/<omega>) so
    that the residual is O(1) and findroot's tolerance is meaningful.
    """
    ys = [(w / mu, p) for w, p in pts]
    Gp_u = lambda u: sum(p * y * gp(mp.mpf(u) * y) for y, p in ys)
    u0 = mp.findroot(Gp_u, X_MAX, solver='secant', tol=mp.mpf('1e-22'))
    G = sum(p * g(u0 * y) for y, p in ys)
    return 1 / G, u0 / mu


def lam_star_expand(sigma, gam3):
    return LAM0 * (1 + KAPPA / 2 * sigma**2 + C3 * gam3 * sigma**3)


def report(tag, pts, lam=None, lam_lit=None):
    pts = prune(pts)
    mu, sig, gam3 = moments(pts)
    lex, s0 = lam_star_exact(pts, mu)
    lap = lam_star_expand(sig, gam3)
    rel = (lap - lex) / lex
    row = ("%-26s %8.4f %8.4f %9.4f %10.4f %10.4f %9.2f%%"
           % (tag, float(mu), float(sig), float(gam3),
              float(lex), float(lap), float(rel) * 100))
    if lam is not None:
        row += "   lam=%.4f  lam/lam*=%.4f" % (float(lam), float(lam / lex))
    if lam_lit is not None:
        row += "   [lit lam=%.2f -> %.4f]" % (lam_lit, lam_lit / float(lex))
    print(row)
    return dict(mu=mu, sigma=sig, gam3=gam3, exact=lex, expand=lap, rel=rel, s0=s0)


def narrow(pts, t):
    """omega -> mu + t*(omega-mu), weights fixed.  sigma -> t*sigma, gamma3 fixed."""
    mu = sum(p * w for w, p in pts)
    return [(mu + t * (w - mu), p) for w, p in pts]


def window(ws, a2f, lo, hi):
    return zip(*[(w, a) for w, a in zip(ws, a2f) if lo <= w <= hi])


# ---------------------------- embedded alpha^2F(omega) data (see header) ----
A2F_STR = """
    0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
    0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00001 0.00001 0.00001 0.00001 0.00001
    0.00001 0.00001 0.00001 0.00001 0.00002 0.00002 0.00002 0.00002 0.00002 0.00003 0.00003 0.00003
    0.00003 0.00004 0.00004 0.00004 0.00004 0.00005 0.00005 0.00005 0.00006 0.00006 0.00007 0.00007
    0.00007 0.00008 0.00008 0.00009 0.00009 0.00010 0.00010 0.00011 0.00011 0.00012 0.00013 0.00013
    0.00014 0.00015 0.00015 0.00016 0.00017 0.00018 0.00018 0.00019 0.00020 0.00021 0.00022 0.00023
    0.00024 0.00025 0.00026 0.00027 0.00028 0.00029 0.00030 0.00031 0.00032 0.00033 0.00034 0.00036
    0.00037 0.00038 0.00040 0.00041 0.00042 0.00044 0.00045 0.00047 0.00048 0.00050 0.00051 0.00053
    0.00055 0.00056 0.00058 0.00060 0.00062 0.00063 0.00065 0.00067 0.00069 0.00071 0.00073 0.00075
    0.00077 0.00079 0.00081 0.00084 0.00086 0.00088 0.00090 0.00093 0.00095 0.00097 0.00100 0.00102
    0.00105 0.00108 0.00110 0.00113 0.00116 0.00118 0.00121 0.00124 0.00127 0.00130 0.00133 0.00136
    0.00139 0.00142 0.00145 0.00149 0.00152 0.00155 0.00158 0.00162 0.00165 0.00169 0.00172 0.00176
    0.00180 0.00183 0.00187 0.00191 0.00195 0.00199 0.00203 0.00207 0.00212 0.00217 0.00224 0.00231
    0.00239 0.00247 0.00255 0.00264 0.00272 0.00281 0.00290 0.00298 0.00306 0.00313 0.00320 0.00328
    0.00335 0.00343 0.00350 0.00358 0.00366 0.00375 0.00383 0.00391 0.00400 0.00409 0.00418 0.00427
    0.00436 0.00445 0.00455 0.00464 0.00474 0.00483 0.00493 0.00503 0.00513 0.00523 0.00533 0.00543
    0.00553 0.00563 0.00574 0.00584 0.00594 0.00605 0.00615 0.00625 0.00636 0.00646 0.00657 0.00667
    0.00677 0.00688 0.00698 0.00708 0.00719 0.00729 0.00739 0.00749 0.00760 0.00770 0.00780 0.00790
    0.00799 0.00809 0.00819 0.00829 0.00838 0.00850 0.00864 0.00879 0.00893 0.00907 0.00922 0.00937
    0.00951 0.00966 0.00981 0.00996 0.01011 0.01026 0.01041 0.01056 0.01072 0.01087 0.01103 0.01118
    0.01134 0.01149 0.01165 0.01181 0.01197 0.01213 0.01229 0.01245 0.01261 0.01277 0.01293 0.01309
    0.01326 0.01342 0.01359 0.01375 0.01392 0.01415 0.01448 0.01482 0.01519 0.01557 0.01597 0.01639
    0.01682 0.01727 0.01773 0.01820 0.01869 0.01920 0.01971 0.02024 0.02077 0.02132 0.02187 0.02244
    0.02301 0.02359 0.02417 0.02477 0.02536 0.02596 0.02657 0.02718 0.02783 0.02851 0.02921 0.02993
    0.03067 0.03143 0.03220 0.03299 0.03378 0.03458 0.03538 0.03617 0.03697 0.03777 0.03855 0.03933
    0.04012 0.04092 0.04173 0.04255 0.04338 0.04421 0.04507 0.04610 0.04735 0.04894 0.05114 0.05387
    0.05722 0.06098 0.06506 0.06901 0.07308 0.07711 0.08102 0.08479 0.08840 0.09035 0.09104 0.09156
    0.09192 0.09214 0.09222 0.09217 0.09207 0.09192 0.09163 0.09123 0.09076 0.09023 0.08965 0.08902
    0.09785 0.09406 0.09068 0.08748 0.08520 0.08274 0.08006 0.10445 0.09712 0.09083 0.09163 0.09262
    0.08760 0.08207 0.07648 0.07121 0.05986 0.06263 0.06534 0.06796 0.07047 0.07283 0.07501 0.07698
    0.07872 0.08019 0.08137 0.08223 0.08227 0.08138 0.07981 0.07842 0.07826 0.07821 0.07826 0.07840
    0.07863 0.07893 0.07930 0.07974 0.08029 0.08145 0.08334 0.08592 0.08918 0.09309 0.09764 0.10280
    0.10855 0.11712 0.12562 0.13552 0.14081 0.13534 0.13097 0.12218 0.11586 0.10900 0.10167 0.09397
    0.08599 0.07783 0.06958 0.05424 0.04985 0.04544 0.04094 0.03630 0.03144 0.02629 0.02078 0.01702
    0.01778 0.01855 0.01934 0.02016 0.02099 0.02185 0.02273 0.02363 0.02456 0.02534 0.02604 0.02673
    0.02742 0.02811 0.02879 0.02948 0.03016 0.03084 0.03153 0.03221 0.03289 0.03358 0.03427 0.03496
    0.03565 0.03634 0.03704 0.03775 0.03846 0.03917 0.03989 0.04062 0.04136 0.04210 0.04285 0.04361
    0.04437 0.04515 0.04594 0.04711 0.04849 0.04997 0.05156 0.05326 0.05507 0.05701 0.05908 0.06128
    0.06362 0.06610 0.06872 0.07150 0.07443 0.07753 0.08079 0.08422 0.08783 0.09162 0.09560 0.09922
    0.10091 0.10236 0.10366 0.10493 0.10629 0.10784 0.10971 0.11200 0.11481 0.11777 0.12069 0.12356
    0.12636 0.12909 0.13175 0.13433 0.13681 0.13989 0.14441 0.14886 0.15323 0.15752 0.16174 0.16588
    0.16994 0.17397 0.17798 0.18194 0.18584 0.18968 0.19346 0.19719 0.20086 0.20448 0.20803 0.21153
    0.21497 0.21836 0.22166 0.22487 0.22802 0.23111 0.23419 0.23725 0.24034 0.24335 0.24625 0.24904
    0.25171 0.25425 0.25668 0.25899 0.26116 0.26316 0.26497 0.26656 0.28887 0.28823 0.28769 0.28727
    0.28696 0.28676 0.28669 0.28675 0.28693 0.28725 0.28770 0.29049 0.29336 0.29590 0.29813 0.30005
    0.30166 0.30297 0.30397 0.30467 0.30506 0.30517 0.30498 0.30450 0.30361 0.30229 0.30056 0.29848
    0.29499 0.29029 0.28530 0.28006 0.27462 0.26449 0.25708 0.25199 0.24778 0.24729 0.25054 0.25453
    0.25923 0.26464 0.27073 0.27660 0.27658 0.27657 0.27659 0.28567 0.28474 0.28375 0.28272 0.28166
    0.28057 0.27946 0.27834 0.27722 0.27611 0.27501 0.27394 0.27289 0.27189 0.27093 0.27003 0.26919
    0.26843 0.26774 0.26715 0.26665 0.32503 0.33932 0.35157 0.36190 0.37039 0.37716 0.38231 0.38594
    0.38815 0.38905 0.36650 0.36762 0.36759 0.36653 0.36452 0.36167 0.35807 0.35383 0.34936 0.34735
    0.35130 0.36476 0.39124 0.43427 0.51437 0.32840 0.28263 0.27399 0.26533 0.25633 0.24695 0.23717
    0.22694 0.21624 0.20504 0.19329 0.18096 0.16803 0.15444 0.14019 0.12522 0.10950 0.09301 0.02315
    0.02328 0.02341 0.02353 0.02367 0.02380 0.02393 0.02406 0.02420 0.02433 0.02447 0.02461 0.02475
    0.02489 0.02503 0.02518 0.02532 0.02546 0.02561 0.02576 0.02591 0.02606 0.02621 0.02636 0.02651
    0.02667 0.02682 0.02698 0.02713 0.02729 0.02745 0.02761 0.02777 0.02794 0.02810 0.02827 0.02843
    0.02860 0.02877 0.02894 0.02911 0.02928 0.02945 0.02963 0.02980 0.02998 0.03016 0.03034 0.03052
    0.03070 0.03088 0.03106 0.03125 0.03143 0.03162 0.03181 0.03200 0.03218 0.03238 0.03257 0.03276
    0.03296 0.03315 0.03335 0.03355 0.03374 0.06999 0.05977 0.09823 0.15975 0.24460 0.35230 0.48234
    0.63423 0.42099 0.39363 0.38649 0.37463 0.35802 0.33713 0.32880 0.34173 0.35958 0.38233 0.40996
    0.44245 0.47978 0.52194 0.56890 0.62079 0.67458 0.72797 0.78073 0.83259 0.88311 0.93211 0.97952
    1.02527 1.06930 1.11154 1.15194 1.19040 1.22687 1.26155 1.29454 1.32590 1.35570 1.38402 1.41092
    1.43648 1.46076 1.48383 1.50576 1.52662 1.54649 1.56543 1.58371 1.60226 1.62196 1.64370 1.66837
    1.69686 1.73005 1.76884 1.81411 1.86675 1.92765 1.99786 2.08055 2.18040 2.30209 2.45032 2.62977
    2.84515 3.10114 3.40245 2.18196 1.49703 1.04723 0.94570 0.91886 0.89142 0.86339 0.83477 0.80554
    0.77571 0.74528 0.68682 0.62940 0.57765 0.53159 0.49123 0.45660 0.42751 0.40298 0.39410 0.38554
    0.37739 0.36974 0.36267 0.35628 0.35063 0.34583 0.34196 0.33910 0.33735 0.33677 0.33747 0.33953
    0.34303 0.34807 0.35471 0.36307 0.37321 0.38522 0.39920 0.41522 0.43337 0.45375 0.45622 0.45784
    0.45875 0.45883 0.45801 0.45620 0.45328 0.44919 0.44382 0.43709 0.42890 0.41916 0.40793 0.42001
    0.50952 0.50066 0.43704 0.36225 0.34435 0.33740 0.33096 0.32499 0.31961 0.31484 0.31060 0.30683
    0.30352 0.30062 0.29809 0.29548 0.29300 0.29081 0.28886 0.28715 0.28564 0.28430 0.28312 0.28205
    0.28109 0.28019 0.27933 0.27849 0.27764 0.27676 0.27586 0.27492 0.27396 0.27297 0.27195 0.27090
    0.26982 0.26871 0.26758 0.26641 0.26522 0.26400 0.26275 0.26147 0.26017 0.25883 0.25747 0.25608
    0.25466 0.25321 0.25174 0.25023 0.24870 0.24714 0.24555 0.24394 0.24230 0.24062 0.23892 0.23720
    0.23544 0.23366 0.23185 0.23001 0.22814 0.22625 0.22433 0.22238 0.22040 0.21840 0.21637 0.21431
    0.21222 0.25620 0.24717 0.23842 0.22992 0.22166 0.21364 0.20585 0.19826 0.27568 0.17832 0.26596
    0.33182 0.36182 0.48203 0.59644 0.77964 0.73005 0.67755 0.64987 0.61982 0.58717 0.55171 0.51322
    0.47146 0.42621 0.37725 0.32436 0.26878 0.23295 0.20108 0.31252 0.29017 0.27352 0.26267 0.25308
    0.18475 0.18088 0.17759 0.17486 0.17269 0.17107 0.16998 0.16941 0.16936 0.16980 0.17073 0.17215
    0.17402 0.17636 0.17913 0.18224 0.16674 0.16848 0.17010 0.17160 0.17297 0.17421 0.17533 0.17631
    0.17718 0.17791 0.17852 0.17900 0.17936 0.17959 0.17969 0.17966 0.17951 0.17923 0.17882 0.17828
    0.17483 0.16969 0.16437 0.15889 0.15326 0.14753 0.14170 0.13581 0.12987 0.12392 0.11798 0.11207
    0.10621 0.10043 0.09474 0.08903 0.08336 0.07775 0.07221 0.06677 0.06144 0.05625 0.05122 0.04647
    0.04214 0.03793 0.03386 0.03009 0.02661 0.02343 0.02051 0.01786 0.01545 0.01328 0.01133 0.00959
    0.00804 0.00667 0.00547 0.00443 0.00352 0.00275 0.00209 0.00154 0.00109 0.00073 0.00046 0.00027
    0.00013 0.00005 0.00001 0.00000 0.00000
"""

W0 = mp.mpf('0.2500007289')      # first grid point, meV
DW = mp.mpf('0.1000009317')      # uniform spacing, meV  (1025 pts -> 102.651 meV)

A2F = [mp.mpf(t) for t in A2F_STR.split()]
WS = [W0 + DW * i for i in range(len(A2F))]
assert len(A2F) == 1025

# ------------------------------------------------------------------ run ----
print("kernel constants (recomputed vs. P3 reference values)")
for name, val in (('x_max', X_MAX), ('g_max', G_MAX), ('lambda_0', LAM0),
                  ('kappa', KAPPA), ('c3', C3)):
    ref = mp.mpf(REF[{'lambda_0': 'lam0'}.get(name, name)])
    print("   %-9s %-22s   ref %-12s  dev %.1e"
          % (name, mp.nstr(val, 12), REF[{'lambda_0': 'lam0'}.get(name, name)],
             float(abs(val - ref) / abs(ref))))
print()

lam, pts = measure(WS, A2F)
print("MgB2 (BETE-NET/JARVIS DFPT, arXiv:2401.16611):  lambda = %.4f" % float(lam))
mu0, sig0, g30 = moments(prune(pts))
print("   w_log check: %.2f meV     mean of P <w> = %.2f meV"
      % (float(mp.e**sum(p * mp.log(w) for w, p in prune(pts))), float(mu0)))
print()

hdr = ("%-26s %8s %8s %9s %10s %10s %10s"
       % ("spectrum", "<w>meV", "sigma", "gamma3", "lam*_exact", "lam*_exp", "rel.err"))
print(hdr)
print("-" * len(hdr))
main = report("MgB2 full a2F", pts, lam=lam, lam_lit=0.87)
print()

print("narrowing control  omega -> <w> + t*(omega-<w>)   (sigma scales as t, gamma3 fixed)")
print("   t=1.0 is the 'MgB2 full a2F' row above.")
print(hdr)
print("-" * len(hdr))
for t in ('0.7', '0.5', '0.3', '0.15'):
    report("  t=%s" % t, narrow(prune(pts), mp.mpf(t)))
print()

print("physical control: dominant (E2g-derived) peak only, 70-81 meV")
print(hdr)
print("-" * len(hdr))
ws_p, a_p = window(WS, A2F, mp.mpf(70), mp.mpf(81))
lam_p, pts_p = measure(list(ws_p), list(a_p))
report("MgB2 peak-only", pts_p, lam=lam_p)
print()
print("optimal kernel scale for the full spectrum: s_opt = %.6f meV^-1,"
      "  1/s_opt = %.2f meV" % (float(main['s0']), float(1 / main['s0'])))
print("(its physical reading depends on the x convention used in P3 sec.4; not asserted here)")
