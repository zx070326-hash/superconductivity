r"""
P3 sec.4 stability-bound check on NON-HYDRIDE spectra: Nb, MgB2, Pb.
Self-contained (imports nothing from this workspace).  mpmath dps=25.
Round 2: adds Nb as the clean kappa test point and a reconciliation column
against Semenok's own published non-hydride lambda_* (arXiv:2407.12922 Tab. S1).

WHAT IS TESTED
--------------
The P3 bound  lambda_*(sigma, gamma3) = lambda_0 (1 + kappa/2 sigma^2
+ c3 gamma3 sigma^3)  is the narrow-spectrum expansion of the shape-exact
functional  lambda_*[P] = 1 / max_{s>0} Int g(s omega) P(omega) domega,
g = Semenok kernel,  P(omega) = (2/lambda) alpha^2F(omega)/omega.

lambda_* is a pure SHAPE functional: rescaling alpha^2F by a constant changes
lambda and leaves lambda_* invariant.  So the coupling used for a material and
the spectrum shape used for it are independent inputs, and we keep them
separate everywhere below (this is the "same-source lambda" trap: we never let
a DFT lambda silently replace a literature lambda inside xi = lambda/lambda_*).

  (1) shape-exact lambda_* by direct integration on real first-principles a2F
  (2) the expansion evaluated at that spectrum's (sigma, gamma3)
  (3) xi_ours = lambda_lit / lambda_*^exact   vs   Semenok's printed xi
  (4) narrowing controls: expansion must converge onto the exact value as
      sigma -> 0, with residual scaling as sigma^4 if kappa and c3 are right
  (5) inverse solve: what shape would be needed to hit Semenok's MgB2 lambda_*

DATA PROVENANCE
---------------
S1  BETE-NET `database.json`, https://github.com/henniggroup/BETE-NET
    (Apache-2.0) = training set of Gibson, Hire, Dee, Barrera, Geisler,
    Hirschfeld, Hennig, arXiv:2401.16611.  High-throughput DFPT.
    Downloaded 2026-09-09:
      curl -sL https://raw.githubusercontent.com/henniggroup/BETE-NET/main/database.json
    Entries used (d["Freq_meV"][k], d["a2F"][k]):
      k="822" Mg1B2_JVASP_19821  lambda=0.6235 wlog=65.45 meV  (1025 pts)
      k="305" Nb_mp_75           lambda=1.2325 wlog=14.55 meV  ( 276 pts)
      k="21"  Pb_mp_20483        lambda=1.0469 wlog= 5.65 meV  (  98 pts)
S2  JARVIS-SuperconDB, figshare 21370572, file 38307921 (CC BY 4.0), fetched
    2026-09-09 via  curl -sL https://ndownloader.figshare.com/files/38307921
    (the figshare *web* page 403s; the ndownloader API does not).  Fields
    a2F_original_x [meV] / a2F_original_y, only 50 points per material:
      JVASP-934   Nb    lambda=0.7741 wlog=221.57 K
      JVASP-19821 MgB2  lambda=0.6959 wlog=846.34 K
    Used only as an independent-DFT shape robustness check, not as a primary.

    NOT obtained: the Kong-Dolgov-Jepsen-Andersen PRB 64, 020501(R) (2001)
    MgB2 alpha^2F itself (figure only; the session's WebSearch budget was
    exhausted).  We therefore do NOT claim to have reproduced Semenok's
    xi=0.16 from Kong's spectrum; instead sec.(5) below inverts our own
    functional to state exactly what shape his lambda_*=5.44 requires, and
    sec.(0) compares wlog / wmax to test the sourcing conjecture.

    All spectra here are machine-readable first-principles output.  Nothing is
    digitized from a figure.

SEMENOK TABLE S1 TARGETS (arXiv:2407.12922; lambda_* = lambda/xi, his Eq. 6)
    MgB2  lambda=0.87        xi=0.16  -> lambda_* = 5.44   wlog=680 K  wmax=1624 K
    Nb    lambda=0.82..1.05  xi=0.20  -> lambda_* = 4.10..5.25   wlog=163 K
    a-Pb  lambda=1.91        xi=0.29  -> lambda_* = 6.59
    xi is printed to 2 dp, so each lambda_* carries a +-3% rounding band.
    NOTE a-Pb is AMORPHOUS lead; our Pb spectrum is crystalline fcc Pb, a
    different material.  That row is shown but is not a like-for-like target.
"""
import mpmath as mp

mp.mp.dps = 25
MEV_K = mp.mpf('11.604518')          # 1 meV in K

# ---------------------------------------------------------------- kernel ----
A = lambda x: mp.im(mp.polygamma(1, 1j * x))      # Im psi'(ix)
B = lambda x: mp.re(mp.polygamma(2, 1j * x))      # Re psi''(ix)


def g(x):
    x = mp.mpf(x)
    return 6 * x**2 + 12 * x**3 * A(x) + 6 * x**4 * B(x)


def gp(x):
    """g'(x), analytic.  A' = B,  B' = -Im psi'''(ix)."""
    x = mp.mpf(x)
    return (12 * x + 36 * x**2 * A(x) + 36 * x**3 * B(x)
            - 6 * x**4 * mp.im(mp.polygamma(3, 1j * x)))


X_MAX = mp.findroot(gp, mp.mpf('0.33'))
G_MAX = g(X_MAX)
LAM0 = 1 / G_MAX
KAPPA = -X_MAX**2 * mp.diff(g, X_MAX, 2) / G_MAX
C3 = -X_MAX**3 * mp.diff(g, X_MAX, 3) / (6 * G_MAX)
REF = dict(x_max='0.32732148', g_max='0.27089315', lambda_0='3.69149',
           kappa='5.60346', c3='0.89590')


# -------------------------------------------------- discrete measure P ------
def measure(ws, a2f):
    """(lambda, [(omega_i, p_i)]), sum p = 1, p = trapezoid weight of a2F/omega."""
    ws = [mp.mpf(w) for w in ws]
    f = [mp.mpf(a) / w for a, w in zip(a2f, ws)]
    n = len(ws)
    wt = []
    for i in range(n):
        lo = ws[i] - ws[i - 1] if i > 0 else mp.mpf(0)
        hi = ws[i + 1] - ws[i] if i < n - 1 else mp.mpf(0)
        wt.append((lo + hi) / 2)
    norm = sum(a * b for a, b in zip(wt, f))          # = lambda/2
    pts = [(w, wi * fi / norm) for w, wi, fi in zip(ws, wt, f) if fi > 0]
    tot = sum(p for _, p in pts)
    return 2 * norm, [(w, p / tot) for w, p in pts]


def coalesce(pts, k):
    """Merge k adjacent atoms, conserving weight and weighted mean position.

    Pure speed device.  Verified against k=1 on the 1025-pt MgB2 measure:
    k=4 shifts sigma by 2e-5 rel., gamma3 by 6e-5 rel., lambda_* by 5e-6 rel.
    """
    out = []
    for i in range(0, len(pts), k):
        grp = pts[i:i + k]
        P = sum(p for _, p in grp)
        if P > 0:
            out.append((sum(w * p for w, p in grp) / P, P))
    tot = sum(p for _, p in out)
    return [(w, p / tot) for w, p in out]


def prep(pts, k=2, tol=mp.mpf('1e-12')):
    return coalesce([(w, p) for w, p in pts if p > tol], k)


def moments(pts):
    mu = sum(p * w for w, p in pts)
    c2 = sum(p * (w - mu)**2 for w, p in pts)
    c3m = sum(p * (w - mu)**3 for w, p in pts)
    return mu, mp.sqrt(c2) / mu, c3m / c2**mp.mpf('1.5')


def wlog(pts):
    return mp.e**sum(p * mp.log(w) for w, p in pts)


def lam_star_exact(pts, mu):
    """1/max_s G(s).  Solved in u = s*<omega> so the residual is O(1)."""
    ys = [(w / mu, p) for w, p in pts]
    Gp_u = lambda u: sum(p * y * gp(mp.mpf(u) * y) for y, p in ys)
    u0 = mp.findroot(Gp_u, X_MAX, solver='secant', tol=mp.mpf('1e-22'))
    return 1 / sum(p * g(u0 * y) for y, p in ys), u0 / mu


def lam_star_expand(sigma, gam3):
    return LAM0 * (1 + KAPPA / 2 * sigma**2 + C3 * gam3 * sigma**3)


def narrow_lin(pts, t):
    """omega -> mu + t(omega-mu).  sigma -> t*sigma exactly, gamma3 invariant.
    Only valid for t <= omega_min/(omega_min-mu) i.e. narrowing, not broadening."""
    mu = sum(p * w for w, p in pts)
    return [(mu + t * (w - mu), p) for w, p in pts]


def scale_pow(pts, t):
    """omega -> mu*(omega/mu)^t.  Stays positive for any t>0, so it can BROADEN.
    sigma ~ t*sigma to leading order; gamma3 drifts (reported, not assumed)."""
    mu = sum(p * w for w, p in pts)
    return [(mu * (w / mu)**t, p) for w, p in pts]


# --------------------- embedded alpha^2F(omega) data (provenance in header) --
# BETE_MGB2 -- BETE-NET database.json entry 822 (Mg1B2_JVASP_19821), lambda=0.6235 wlog=65.4484 meV
BETE_MGB2_W = ('0.2500007289', '0.1000009317', 1025)   # (w0, dw, n) uniform grid, meV; max grid dev 9.8e-05
BETE_MGB2_A = """
    8.9e-09 2.45e-08 5.21e-08 9.5e-08 1.569e-07 2.41e-07 3.508e-07 4.898e-07 6.613e-07 8.689e-07 1.1158e-06 1.4056e-06
    1.7416e-06 2.1274e-06 2.5663e-06 3.0617e-06 3.6172e-06 4.236e-06 4.9217e-06 5.6777e-06 6.5073e-06 7.4141e-06 8.4014e-06 9.4727e-06
    1.06314e-05 1.1881e-05 1.32248e-05 1.46663e-05 1.62089e-05 1.78561e-05 1.96112e-05 2.14778e-05 2.34592e-05 2.55588e-05 2.77801e-05 3.01265e-05
    3.26015e-05 3.52084e-05 3.79507e-05 4.08319e-05 4.38553e-05 4.70244e-05 5.03425e-05 5.38133e-05 5.74399e-05 6.1226e-05 6.51749e-05 6.929e-05
    7.35748e-05 7.80327e-05 8.26672e-05 8.74816e-05 9.24793e-05 9.76639e-05 0.000103039 0.000108607 0.000114373 0.000120339 0.000126509 0.000132886
    0.000139475 0.000146277 0.000153297 0.000160538 0.000168004 0.000175698 0.000183623 0.000191783 0.000200181 0.000208821 0.000217705 0.000226839
    0.000236224 0.000245865 0.000255764 0.000265926 0.000276353 0.000287049 0.000298018 0.000309263 0.000320787 0.000332594 0.000344687 0.00035707
    0.000369746 0.000382719 0.000395991 0.000409566 0.000423449 0.000437642 0.000452148 0.000466971 0.000482115 0.000497583 0.000513378 0.000529504
    0.000545964 0.000562762 0.000579901 0.000597384 0.000615215 0.000633398 0.000651935 0.00067083 0.000690088 0.00070971 0.0007297 0.000750063
    0.000770801 0.000791918 0.000813416 0.000835301 0.000857574 0.00088024 0.000903302 0.000926764 0.000950628 0.000974898 0.000999578 0.00102467
    0.00105018 0.00107611 0.00110246 0.00112924 0.00115645 0.00118409 0.00121217 0.00124069 0.00126966 0.00129907 0.00132893 0.00135925
    0.00139002 0.00142126 0.00145296 0.00148513 0.00151777 0.00155088 0.00158447 0.00161854 0.0016531 0.00168815 0.00172369 0.00175972
    0.00179626 0.00183329 0.00187083 0.00190888 0.00194744 0.00198652 0.00202612 0.00206691 0.00211554 0.00217284 0.0022377 0.002309
    0.00238564 0.00246652 0.00255053 0.00263655 0.00272348 0.00281022 0.00289566 0.00297868 0.00305818 0.00313306 0.00320391 0.00327613
    0.00335015 0.00342593 0.00350341 0.00358254 0.00366326 0.00374553 0.00382929 0.00391449 0.00400107 0.00408898 0.00417818 0.0042686
    0.00436019 0.00445291 0.00454669 0.00464149 0.00473726 0.00483393 0.00493146 0.00502979 0.00512888 0.00522866 0.00532909 0.00543011
    0.00553168 0.00563373 0.00573622 0.00583908 0.00594228 0.00604576 0.00614945 0.00625332 0.00635731 0.00646136 0.00656543 0.00666945
    0.00677338 0.00687717 0.00698075 0.00708409 0.00718712 0.00728979 0.00739206 0.00749386 0.00759515 0.00769586 0.00779596 0.00789538
    0.00799408 0.008092 0.00818908 0.00828528 0.00838055 0.00850098 0.00864292 0.00878573 0.00892942 0.00907396 0.00921936 0.00936559
    0.00951264 0.00966051 0.00980917 0.00995862 0.0101089 0.0102598 0.0104116 0.0105641 0.0107173 0.0108712 0.0110258 0.0111811
    0.0113371 0.0114938 0.0116511 0.011809 0.0119676 0.0121268 0.0122866 0.012447 0.012608 0.0127696 0.0129317 0.0130943
    0.0132575 0.0134212 0.0135854 0.0137501 0.0139153 0.0141519 0.0144794 0.0148243 0.015188 0.0155698 0.0159692 0.0163853
    0.0168176 0.0172653 0.0177277 0.0182041 0.0186939 0.0191963 0.0197108 0.0202365 0.0207728 0.021319 0.0218744 0.0224383
    0.0230101 0.023589 0.0241743 0.0247655 0.0253617 0.0259623 0.0265666 0.0271829 0.0278306 0.0285061 0.0292065 0.0299291
    0.0306712 0.0314301 0.0322032 0.0329877 0.0337796 0.0345761 0.0353751 0.0361744 0.0369717 0.0377651 0.0385522 0.039331
    0.0401179 0.0409223 0.041732 0.0425492 0.0433759 0.0442144 0.0450738 0.0461042 0.0473505 0.0489359 0.051143 0.0538724
    0.0572195 0.0609841 0.0650608 0.0690119 0.0730811 0.0771051 0.0810155 0.0847871 0.088401 0.0903496 0.0910397 0.0915604
    0.0919227 0.0921374 0.0922155 0.092168 0.0920728 0.0919208 0.0916345 0.091232 0.0907621 0.0902334 0.0896508 0.0890193
    0.0978502 0.094059 0.0906832 0.0874816 0.0851959 0.0827358 0.0800604 0.104451 0.097123 0.0908313 0.091633 0.0926177
    0.0875962 0.0820747 0.076475 0.0712076 0.0598611 0.0626309 0.0653418 0.0679644 0.0704692 0.0728268 0.0750075 0.076982
    0.0787207 0.0801942 0.081373 0.0822276 0.0822684 0.0813758 0.0798082 0.0784174 0.0782641 0.0782147 0.0782628 0.078402
    0.0786257 0.0789276 0.0793012 0.07974 0.0802905 0.0814533 0.0833368 0.0859191 0.0891784 0.093093 0.0976411 0.102801
    0.108552 0.117122 0.125615 0.135517 0.14081 0.135343 0.13097 0.122183 0.115865 0.108996 0.101667 0.0939692
    0.0859945 0.0778347 0.0695812 0.0542397 0.0498452 0.045436 0.0409438 0.0363003 0.0314373 0.0262864 0.0207793 0.0170222
    0.0177764 0.01855 0.0193435 0.0201575 0.0209926 0.0218494 0.0227285 0.0236305 0.024556 0.0253442 0.0260386 0.026731
    0.0274212 0.0281087 0.0287943 0.0294784 0.0301614 0.0308439 0.0315263 0.0322092 0.0328931 0.0335785 0.0342658 0.0349556
    0.0356483 0.0363445 0.0370447 0.0377493 0.0384589 0.0391739 0.0398949 0.0406223 0.0413567 0.0420986 0.0428483 0.0436066
    0.0443737 0.0451503 0.0459369 0.0471123 0.0484896 0.0499696 0.0515572 0.0532574 0.0550749 0.0570148 0.059082 0.0612812
    0.0636176 0.0660958 0.068721 0.0714979 0.0744314 0.0775266 0.0807882 0.0842213 0.0878306 0.0916211 0.0955978 0.0992172
    0.100913 0.102355 0.103656 0.104929 0.106287 0.107843 0.109711 0.112002 0.114814 0.117774 0.120691 0.123556
    0.126358 0.129093 0.131752 0.134328 0.136811 0.139891 0.144414 0.14886 0.15323 0.157523 0.161739 0.165879
    0.169941 0.173966 0.17798 0.181937 0.185837 0.18968 0.193465 0.197194 0.200865 0.204479 0.208035 0.211534
    0.214975 0.218358 0.221659 0.224872 0.228016 0.231114 0.234186 0.237254 0.240337 0.243349 0.24625 0.249036
    0.251705 0.254253 0.256684 0.258994 0.261159 0.263157 0.264967 0.266565 0.288874 0.288229 0.287692 0.287267
    0.286957 0.286764 0.286694 0.286749 0.286934 0.28725 0.287703 0.290491 0.293356 0.295903 0.298135 0.300054
    0.301664 0.302968 0.303967 0.304665 0.305064 0.305167 0.304977 0.304495 0.303612 0.302289 0.300563 0.298476
    0.294992 0.290291 0.285299 0.280062 0.274624 0.264485 0.257077 0.251992 0.247779 0.247286 0.250538 0.254525
    0.259231 0.264639 0.270731 0.276597 0.276576 0.276572 0.276595 0.285674 0.284739 0.283753 0.282724 0.28166
    0.28057 0.279462 0.278344 0.277225 0.276112 0.275013 0.273938 0.272894 0.271889 0.270931 0.27003 0.269193
    0.268428 0.267743 0.267147 0.266648 0.325032 0.33932 0.351574 0.361897 0.370392 0.377162 0.38231 0.385939
    0.388151 0.389049 0.366497 0.367616 0.367593 0.366528 0.36452 0.361667 0.35807 0.353827 0.349363 0.347346
    0.351304 0.364762 0.391243 0.434272 0.514367 0.328399 0.28263 0.273993 0.26533 0.256327 0.246951 0.237168
    0.226944 0.216245 0.205037 0.193287 0.180961 0.168025 0.154444 0.140187 0.125218 0.109504 0.0930108 0.0231495
    0.0232767 0.0234051 0.0235346 0.0236652 0.0237969 0.0239298 0.0240637 0.0241988 0.0243349 0.0244722 0.0246106 0.0247501
    0.0248907 0.0250325 0.0251753 0.0253193 0.0254643 0.0256105 0.0257578 0.0259063 0.0260558 0.0262065 0.0263582 0.0265111
    0.0266651 0.0268202 0.0269765 0.0271338 0.0272923 0.0274519 0.0276126 0.0277745 0.0279374 0.0281015 0.0282667 0.028433
    0.0286004 0.028769 0.0289387 0.0291095 0.0292814 0.0294544 0.0296286 0.0298039 0.0299803 0.0301578 0.0303365 0.0305162
    0.0306971 0.0308792 0.0310623 0.0312466 0.031432 0.0316185 0.0318062 0.031995 0.0321849 0.0323759 0.0325681 0.0327613
    0.0329558 0.0331513 0.033348 0.0335458 0.0337447 0.0699904 0.0597742 0.0982337 0.159752 0.244603 0.352296 0.482337
    0.634233 0.42099 0.393627 0.386493 0.374627 0.358024 0.337127 0.328801 0.34173 0.359578 0.382327 0.409956
    0.442446 0.47978 0.521936 0.568896 0.620787 0.674579 0.727974 0.780731 0.832594 0.883115 0.93211 0.979515
    1.02527 1.0693 1.11154 1.15194 1.1904 1.22687 1.26155 1.29454 1.3259 1.3557 1.38402 1.41092
    1.43648 1.46076 1.48383 1.50576 1.52662 1.54649 1.56543 1.58371 1.60226 1.62196 1.6437 1.66837
    1.69686 1.73005 1.76884 1.81411 1.86675 1.92765 1.99786 2.08055 2.1804 2.30209 2.45032 2.62977
    2.84515 3.10114 3.40245 2.18196 1.49703 1.04723 0.945702 0.918859 0.891424 0.863393 0.834766 0.805539
    0.775712 0.745281 0.686822 0.629399 0.577648 0.531587 0.491232 0.456603 0.427507 0.402981 0.394097 0.385538
    0.37739 0.36974 0.362673 0.356276 0.350634 0.345834 0.341962 0.339103 0.337345 0.336773 0.337472 0.33953
    0.343032 0.348065 0.354714 0.363065 0.373206 0.38522 0.399196 0.415218 0.433374 0.453748 0.456217 0.457844
    0.458747 0.458834 0.458014 0.456195 0.453284 0.449192 0.443824 0.437091 0.4289 0.41916 0.407933 0.420012
    0.509515 0.500659 0.437038 0.362249 0.344354 0.337404 0.330964 0.324992 0.319613 0.314844 0.310596 0.306832
    0.303517 0.300616 0.298095 0.295484 0.293004 0.290806 0.288863 0.28715 0.285638 0.284303 0.283117 0.282053
    0.281086 0.280187 0.279332 0.278492 0.277642 0.276765 0.275859 0.274924 0.27396 0.272968 0.271947 0.270897
    0.269819 0.268712 0.267577 0.266413 0.26522 0.263999 0.26275 0.261473 0.260166 0.258832 0.25747 0.256079
    0.25466 0.253212 0.251737 0.250233 0.248702 0.247142 0.245554 0.243938 0.242295 0.240623 0.238924 0.237196
    0.235441 0.233658 0.231847 0.230009 0.228143 0.226249 0.224328 0.222378 0.220402 0.218398 0.216366 0.214307
    0.21222 0.256199 0.247174 0.238418 0.229918 0.221664 0.213644 0.205849 0.198265 0.275677 0.178317 0.265963
    0.331824 0.361816 0.482033 0.596436 0.779642 0.730047 0.677549 0.649866 0.619815 0.587172 0.551714 0.513216
    0.471455 0.426208 0.37725 0.324359 0.268778 0.232949 0.201082 0.312525 0.290172 0.273523 0.262668 0.253076
    0.184747 0.180877 0.177586 0.174862 0.172694 0.17107 0.169981 0.169413 0.169358 0.169802 0.170735 0.172146
    0.174024 0.176357 0.179135 0.182243 0.166736 0.168483 0.170104 0.171599 0.172967 0.174209 0.175325 0.176314
    0.177176 0.177912 0.178521 0.179003 0.179358 0.179586 0.179687 0.179661 0.179508 0.179227 0.178819 0.178284
    0.174831 0.169694 0.164371 0.158886 0.153263 0.147526 0.141699 0.135807 0.129874 0.123923 0.11798 0.112068
    0.106211 0.100435 0.094739 0.0890346 0.0833648 0.0777502 0.0722112 0.0667682 0.0614417 0.0562521 0.05122 0.0464662
    0.0421361 0.0379263 0.0338642 0.0300895 0.026615 0.0234274 0.0205135 0.01786 0.0154539 0.0132818 0.0113306 0.009587
    0.00803789 0.00667003 0.00547021 0.00442522 0.00352186 0.00274692 0.00208719 0.00153528 0.00108784 0.000734422 0.000464526 0.000267644
    0.000133268 5.08903e-05 1.0005e-05 0 0
"""

# BETE_NB -- BETE-NET database.json entry 305 (Nb_mp_75), lambda=1.2325 wlog=14.5529 meV
BETE_NB_W = ('0.2500020894', '0.1000006972', 276)   # (w0, dw, n) uniform grid, meV; max grid dev 1.2e-04
BETE_NB_A = """
    6.0361e-06 1.6563e-05 3.52024e-05 6.42722e-05 0.00010609 0.000162974 0.000237242 0.000331212 0.000447202 0.000587529 0.000754511 0.000950466
    0.00117771 0.00143857 0.00173535 0.00207038 0.00244597 0.00286444 0.00332811 0.00383929 0.00440031 0.00501348 0.00568112 0.00640554
    0.00718908 0.00803403 0.00894273 0.00991748 0.0109606 0.0120744 0.0132613 0.0145235 0.0158633 0.0172831 0.0187851 0.0203718
    0.0220454 0.0238082 0.0256626 0.0276108 0.0296553 0.0317983 0.034042 0.036389 0.0388414 0.0936365 0.0905089 0.0862454
    0.0823175 0.0801969 0.0813495 0.0844867 0.0877621 0.0911597 0.0946633 0.0982569 0.101924 0.10565 0.109417 0.11321
    0.117012 0.120808 0.124581 0.128316 0.131996 0.135605 0.139127 0.142546 0.145846 0.149012 0.152026 0.154872
    0.157536 0.16 0.162248 0.160187 0.16166 0.163042 0.164334 0.165534 0.166643 0.16766 0.157867 0.159914
    0.161973 0.164045 0.166129 0.168226 0.170335 0.172456 0.17459 0.176735 0.178893 0.181063 0.183245 0.18544
    0.187646 0.189864 0.192094 0.194337 0.196591 0.198857 0.201134 0.203424 0.205725 0.208038 0.210363 0.212699
    0.215047 0.217407 0.219778 0.22216 0.224554 0.22696 0.229377 0.231805 0.234244 0.236695 0.239157 0.241631
    0.244115 0.246611 0.249118 0.251636 0.25812 0.263398 0.268769 0.274232 0.285428 0.29141 0.29749 0.303669
    0.309945 0.316317 0.322784 0.329346 0.336001 0.343339 0.353501 0.366597 0.382709 0.404173 0.432526 0.468789
    0.51398 0.575576 0.626777 0.679318 0.732356 0.785549 0.838753 0.891823 0.947597 1.00526 1.05457 1.10281
    1.1499 1.19576 1.24029 1.5148 1.42745 1.34503 1.26715 1.19345 1.1236 1.05794 0.996357 0.938593
    0.924408 0.889475 0.885735 0.713143 0.692239 0.673078 0.655657 0.63997 0.626013 0.613783 0.550104 0.546369
    0.542454 0.53836 0.534089 0.529643 0.525022 0.520228 0.515263 0.510129 0.504826 0.499356 0.493721 0.487921
    0.48196 0.460593 0.450169 0.450364 0.459344 0.474251 0.492228 0.510415 0.444127 0.435988 0.430961 0.44248
    0.439788 0.436794 0.433785 0.430883 0.42811 0.42549 0.423046 0.497882 0.442427 0.412731 0.407276 0.418878
    0.434691 0.454336 0.477855 0.505289 0.53668 0.572068 0.607853 0.634925 0.663009 0.691536 0.720153 0.748568
    0.77649 0.803629 0.829693 0.925227 1.13849 1.24434 1.20735 1.22713 1.23635 1.23521 1.22395 1.20278
    1.17192 1.1316 1.17977 0.899734 0.437389 0.323411 0.255182 0.226431 0.210793 0.195889 0.181704 0.168224
    0.155434 0.143319 0.131865 0.121057 0.110881 0.110453 0.0988576 0.088104 0.0781701 0.0679742 0.0558518 0.0451808
    0.0358788 0.0278633 0.0210519 0.015362 0.0107113 0.0070172 0.00419727 0.00216906 0.000850088 0.000157898 0 0
"""

# BETE_PB -- BETE-NET database.json entry 21 (Pb_mp_20483), lambda=1.0469 wlog=5.6491 meV
BETE_PB_W = ('0.2500020894', '0.1000008302', 98)   # (w0, dw, n) uniform grid, meV; max grid dev 1.3e-05
BETE_PB_A = """
    2.19998e-05 6.03675e-05 0.000128303 0.000234254 0.000386669 0.000593995 0.000864681 0.00120717 0.00162992 0.00214138 0.00274998 0.00346418
    0.00429243 0.00524318 0.00632486 0.00754594 0.00891486 0.0104401 0.01213 0.0139931 0.0160379 0.0182727 0.0207061 0.0233464
    0.0262022 0.0353355 0.0542463 0.0737528 0.0938423 0.114502 0.135721 0.157484 0.179781 0.20258 0.228068 0.259656
    0.298349 0.34352 0.385866 0.421548 0.450746 0.473055 0.488071 0.502003 0.524357 0.553436 0.58819 0.55801
    0.567819 0.586221 0.608682 0.631547 1.47555 0.613408 0.613297 0.61213 0.609909 0.460637 0.321935 0.341334
    0.361637 0.382832 0.40491 0.42786 0.451669 0.627294 0.614708 0.679968 0.749182 0.856155 0.119516 0.1351
    0.151555 0.168879 0.187074 0.206137 0.233289 0.266289 0.302369 0.341276 0.382755 0.426554 0.563634 0.792594
    0.995166 1.09508 1.17117 1.07599 0.873134 0.601361 0.366932 0.156149 0.103272 0.0619407 0.0315678 0.0115677
    0.0013543 0
"""

# JARVIS_NB -- JARVIS-SuperconDB JVASP-934, lambda=0.7741 wlog=221.57 K
JARVIS_NB_W = ('0.2758759368', '0.5517509018', 50)   # (w0, dw, n) uniform grid, meV; max grid dev 7.3e-05
JARVIS_NB_A = """
    3.16217e-07 8.53786e-06 3.95271e-05 0.000144209 0.000427237 0.000812609 0.00133804 0.00212012 0.00319081 0.00447088 0.00610325 0.00817591
    0.0106238 0.0136671 0.0172665 0.0214696 0.0264115 0.0323082 0.0392549 0.0471254 0.0565354 0.0674052 0.0800239 0.0950079
    0.112484 0.13305 0.158149 0.187964 0.225577 0.274027 0.339005 0.450379 0.560156 0.624002 0.703042 0.802369
    0.721265 0.82995 1.09988 1.04151 1.61389 1.05482 0.664289 0.477954 0.355712 0.266738 0.191981 0.128865
    0.0627881 0
"""

# JARVIS_MGB2 -- JARVIS-SuperconDB JVASP-19821, lambda=0.6959 wlog=846.34 K
JARVIS_MGB2_W = ('1.0250492106', '2.0500983379', 50)   # (w0, dw, n) uniform grid, meV; max grid dev 6.6e-05
JARVIS_MGB2_A = """
    -1.23397e-07 -3.33173e-06 -3.04123e-05 -7.27809e-05 -0.000107415 -0.000276373 -0.000249201 -0.000543535 -0.000618193 -0.000730914 -0.000790262 -0.000202749
    0.000926336 0.00501254 0.0210105 0.0698743 0.110488 0.0417474 0.0626546 0.0177375 0.00754563 0.00160116 -0.00355863 0.0522098
    0.163513 0.202776 0.199526 0.153224 0.173157 0.196048 0.264444 0.314863 0.189827 0.134874 0.24332 0.561207
    0.842102 0.907636 1.0115 1.023 0.929155 0.884957 0.844509 0.844432 0.842295 1.01621 0.335484 0.0331121
    0.00384081 0
"""


# ------------------------------------------------------------ spectra ------
def grid(spec):
    w0, dw, n = spec
    return [mp.mpf(w0) + mp.mpf(dw) * i for i in range(n)]


SPECTRA = {}
for nm, W, Astr in (("BETE_MGB2", BETE_MGB2_W, BETE_MGB2_A),
                    ("BETE_NB", BETE_NB_W, BETE_NB_A),
                    ("BETE_PB", BETE_PB_W, BETE_PB_A),
                    ("JARVIS_NB", JARVIS_NB_W, JARVIS_NB_A),
                    ("JARVIS_MGB2", JARVIS_MGB2_W, JARVIS_MGB2_A)):
    a = [mp.mpf(t) for t in Astr.split()]
    ws = grid(W)
    assert len(ws) == len(a), (nm, len(ws), len(a))
    SPECTRA[nm] = (ws, a)

SEMENOK = {                       # arXiv:2407.12922 Table S1, lambda_* = lambda/xi
    "MgB2": dict(lam=(0.87, 0.87), xi=0.16, wlog_K=680, wmax_K=1624),
    "Nb":   dict(lam=(0.82, 1.05), xi=0.20, wlog_K=163, wmax_K=None),
    "a-Pb": dict(lam=(1.91, 1.91), xi=0.29, wlog_K=None, wmax_K=None),
}


def sem_band(tag):
    """lambda_* implied by Table S1, incl. the +-0.005 rounding band on xi."""
    s = SEMENOK[tag]
    return (s["lam"][0] / (s["xi"] + 0.005), s["lam"][1] / (s["xi"] - 0.005))


def wmax_of(ws, a, frac=mp.mpf('1e-3')):
    thr = max(a) * frac
    return max(w for w, v in zip(ws, a) if v >= thr)


# =============================================================== 0. setup ===
print("=" * 100)
print("0.  kernel constants, recomputed here vs. the P3 reference values")
for name, val in (('x_max', X_MAX), ('g_max', G_MAX), ('lambda_0', LAM0),
                  ('kappa', KAPPA), ('c3', C3)):
    ref = mp.mpf(REF[name])
    print("    %-9s %-24s ref %-11s  rel.dev %.1e"
          % (name, mp.nstr(val, 12), REF[name], float(abs(val - ref) / abs(ref))))

print()
print("0b. spectrum vs. Semenok Table S1 frequency scales (tests same-source conjecture)")
print("    %-13s %10s %10s %10s   %s" % ("spectrum", "wlog[K]", "wmax[K]", "wmax/wlog", "Semenok wlog / wmax [K]"))
for nm, tag in (("BETE_MGB2", "MgB2"), ("JARVIS_MGB2", "MgB2"),
                ("BETE_NB", "Nb"), ("JARVIS_NB", "Nb"), ("BETE_PB", "a-Pb")):
    ws, a = SPECTRA[nm]
    _, pts = measure(ws, a)
    wl = wlog(prep(pts, 1)) * MEV_K
    wm = wmax_of(ws, a) * MEV_K
    s = SEMENOK[tag]
    print("    %-13s %10.1f %10.1f %10.2f   %s / %s"
          % (nm, float(wl), float(wm), float(wm / wl),
             s["wlog_K"], s["wmax_K"]))

# ============================================= 1. main shape-exact table ====
print()
print("=" * 100)
print("1.  shape-exact lambda_* vs. the P3 expansion, and vs. Semenok's published lambda_*")
print("    (lambda_lit is the coupling Semenok used; lambda_dft is this spectrum's own")
print("     coupling.  lambda_* is shape-only, so only the shape comes from the spectrum.)")
hdr = ("%-13s %8s %8s %8s %9s %10s %10s %8s %8s %8s"
       % ("spectrum", "lam_dft", "<w>meV", "sigma", "gamma3",
          "lam*_exact", "lam*_exp", "rel.err", "xi_ours", "xi_Sem"))
print(hdr)
print("-" * len(hdr))
RES = {}
for nm, tag, k in (("BETE_NB", "Nb", 2), ("JARVIS_NB", "Nb", 1),
                   ("BETE_MGB2", "MgB2", 4), ("JARVIS_MGB2", "MgB2", 1),
                   ("BETE_PB", "a-Pb", 1)):
    ws, a = SPECTRA[nm]
    lam, pts = measure(ws, a)
    q = prep(pts, k)
    mu, sig, gam3 = moments(q)
    lex, s0 = lam_star_exact(q, mu)
    lap = lam_star_expand(sig, gam3)
    s = SEMENOK[tag]
    xi_lo, xi_hi = s["lam"][0] / lex, s["lam"][1] / lex
    xi_s = "%.3f-%.3f" % (float(xi_lo), float(xi_hi)) if xi_lo != xi_hi else "%.3f" % float(xi_lo)
    print("%-13s %8.4f %8.3f %8.4f %9.4f %10.4f %10.4f %7.2f%% %8s %8.2f"
          % (nm, float(lam), float(mu), float(sig), float(gam3),
             float(lex), float(lap), float((lap - lex) / lex) * 100, xi_s, s["xi"]))
    RES[nm] = dict(lam=lam, mu=mu, sigma=sig, gam3=gam3, exact=lex, pts=q, tag=tag)

print()
print("    reconciliation: our lambda_*^exact against Semenok Table S1 lambda_* = lambda/xi")
print("    %-13s %14s %20s %10s" % ("spectrum", "lam*_exact", "Semenok lam* band", "inside?"))
for nm in ("BETE_NB", "JARVIS_NB", "BETE_MGB2", "JARVIS_MGB2", "BETE_PB"):
    lo, hi = sem_band(RES[nm]["tag"])
    L = float(RES[nm]["exact"])
    print("    %-13s %14.4f %20s %10s"
          % (nm, L, "%.2f - %.2f" % (lo, hi), "YES" if lo <= L <= hi else "NO"))

# ================================================ 2. narrowing controls =====
print()
print("=" * 100)
print("2.  narrowing control  omega -> <w> + t(omega-<w>)   [sigma -> t*sigma, gamma3 fixed]")
print("    A wrong kappa leaves an O(t^2) residual, a wrong c3 an O(t^3) one; both would")
print("    make err/t^3 diverge as t->0.  Correct coefficients leave O(t^4).")
h2 = ("%-16s %8s %9s %10s %10s %12s %10s %10s"
      % ("case", "sigma", "gamma3", "lam*_exact", "lam*_exp", "rel.err[ppm]",
         "err/t^3", "err/t^4"))
for nm, k in (("BETE_NB", 2), ("BETE_MGB2", 4)):
    base = RES[nm]
    print()
    print("  %s   (t=1 row is the exact/expansion pair from table 1)" % nm)
    print(h2)
    print("-" * len(h2))
    for t in ('1.0', '0.7', '0.5', '0.3', '0.15'):
        tt = mp.mpf(t)
        q = base["pts"] if t == '1.0' else narrow_lin(base["pts"], tt)
        mu, sig, gam3 = moments(q)
        lex, _ = lam_star_exact(q, mu)
        lap = lam_star_expand(sig, gam3)
        rel = (lap - lex) / lex
        print("  t=%-14s %8.4f %9.4f %10.4f %10.4f %12.2f %10.4f %10.4f"
              % (t, float(sig), float(gam3), float(lex), float(lap),
                 float(rel) * 1e6, float(rel / tt**3), float(rel / tt**4)))

# ============================================ 3. inverse solve for Semenok ==
print()
print("=" * 100)
print("3.  inverse solve: what shape reproduces Semenok's MgB2 lambda_* = 5.44?")
print("    Family: omega -> <w>*(omega/<w>)^t applied to the BETE_MGB2 shape.")
print("    t=1 is the DFT spectrum; t>1 broadens it.")
base = prep(measure(*SPECTRA["BETE_MGB2"])[1], 8)
target = mp.mpf('0.87') / mp.mpf('0.16')


def lam_star_at(t):
    q = scale_pow(base, mp.mpf(t))
    mu, _, _ = moments(q)
    return lam_star_exact(q, mu)[0]


print("    %8s %10s %10s %12s" % ("t", "sigma", "gamma3", "lam*_exact"))
for t in ('1.0', '1.3', '1.6'):
    q = scale_pow(base, mp.mpf(t))
    mu, sig, gam3 = moments(q)
    print("    %8s %10.4f %10.4f %12.4f"
          % (t, float(sig), float(gam3), float(lam_star_exact(q, mu)[0])))
t_star = mp.findroot(lambda t: lam_star_at(t) - target, mp.mpf('1.4'),
                     solver='secant', tol=mp.mpf('1e-16'))
q = scale_pow(base, t_star)
mu, sig, gam3 = moments(q)
print("    solved: t* = %.4f  ->  sigma = %.4f, gamma3 = %.4f, lam*_exact = %.4f (target %.4f)"
      % (float(t_star), float(sig), float(gam3),
         float(lam_star_exact(q, mu)[0]), float(target)))
print("    i.e. Semenok's MgB2 lambda_*=5.44 needs sigma = %.3f, vs %.3f for the DFT spectrum"
      % (float(sig), float(RES["BETE_MGB2"]["sigma"])))
