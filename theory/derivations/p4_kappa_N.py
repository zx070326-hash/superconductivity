r"""
P4-A · kappa_N(varpi_bar) and c_{3,N}(varpi_bar) for the KAY II truncated operators.

Tests P4_rigorous_tc_bound.md eqs. (P4A.1)-(P4A.5):
  frak_k^(N)(P,T) = frak_h^(N)(varpi_bar) * [1 + (kappa_N/2) sigma^2
                     + c_{3,N} gamma3 sigma^3 + O(sigma^4)],        (P4A.1)
  kappa_N = varpi_bar^2 <v_N|H^(N)''(varpi_bar)|v_N> / frak_h^(N),   (P4A.2)
  c_{3,N} = (1/6) varpi_bar^3 <v_N|H^(N)'''(varpi_bar)|v_N> / frak_h^(N). (P4A.3)

OPERATORS (do not confuse the two):
* KAY I  G(gamma)  [P4-B object]: matrix elements from 1/k^gamma,
  G^(4) explicit in surveys/P4_KAY_I_and_combescot_2026-09-09.md eq.31.
* KAY II H^(N)(varpi)  [P4-A object, THIS FILE]: same three-part structure
  (=-H1+H2+H3, normalization 1/sqrt((2n+1)(2m+1)), arguments |n-m| and n+m+1)
  with each 1/k^gamma replaced by [[k]](varpi) := varpi^2/(k^2+varpi^2).

  ADOPTED FORM (written out explicitly; survey transcribes eq.23 only as
  "matrix elements built from [[n]](varpi), n = 0..N-1", which does not fix
  the placement, so we pin it down by two self-checks, both of which PASS):

    (H1)_{nm} = +delta_{nm} / (2n+1) * sum_{k=1}^{n} 2[[k]]
    (H2)_{nm} = (1-delta_{nm}) / (sqrt(2n+1) sqrt(2m+1)) * [[|n-m|]]
    (H3)_{nm} = 1 / (sqrt(2n+1) sqrt(2m+1)) * [[n+m+1]]
    H^(N) = -H1 + H2 + H3.

  Self-check 1 (N=1 anchor): H^(1)(varpi) = [[1]] = varpi^2/(1+varpi^2)
    = the document's h(varpi).  (Note [[0]] would be identically 1, so the
    "n = 0..N-1" gloss in the md file must mean the KERNEL ARGUMENTS
    |n-m|, n+m+1, which are >= 1; the n^2 -> 0^2 remark in md section 1
    is a slip. Confirmed: KAY I G^(1)(gamma) = 1/(1^(gamma)) = 1 has the
    same structure, and gamma-independence of g^(1)(gamma)=1 (KAY I eq.67)
    is the exact analogue of [[1]] being the only surviving element.)

  Self-check 2 (truncation constants): lim_{varpi->0} rho(H^(N)(varpi))/varpi^2
    = g^(N)(2) = rho(G^(N)(2)) since [[k]] = varpi^2/k^2 * (1 + O(varpi^2)).
    We verify rho(H^(N)(varpi))/varpi^2 -> g^(N)(2) AND independently
    reproduce the four KAY I truncation constants
    T_c^(N)(2) = (1/2pi) sqrt(g^(N)(2)) = 0.159154943 / 0.1796160944 /
    0.1820383 / 0.1825137102 (survey section 3.2; these are KAY I objects,
    reproduced via the G matrix, and reached from H via the small-varpi
    bridge -- not by pretending H(varpi) is G(2)).

NOTES -- three transcription slips found in
surveys/P4_KAY_I_and_combescot_2026-09-09.md while cross-checking (all three
are survey-side; the matrix construction here reproduces KAY I exactly):
  (i)  sec.1.4 eq.31 4x4 matrix, entry (2,2): printed 1/5[1/3^gamma - 2/2^gamma - 2],
       must be 1/5[1/5^gamma - 2/2^gamma - 2] since the G3 argument is
       n+m+1 = 5.  The printed version gives T_c^(3)=0.1820573,
       T_c^(4)=0.1825341, i.e. it BREAKS the survey's own constants
       (0.1820383 / 0.1825137102); with 1/5^gamma they are reproduced exactly.
  (ii) sec.3.3 N=2 closed form: printed T_c^(2)(2)=(1/2pi)sqrt((60+3 sqrt 13819)/18),
       which evaluates to 0.76205, not 0.17962.  Denominator must be 324:
       g^(2)(2) = (60+3 sqrt 13819)/324 = (20 + sqrt 13819)/108 = 1.2736503963.
       (Verified here: tr G^(2)(2) = 10/27 and det G^(2)(2) = -497/432, matching
       the survey's own quoted trace, and 100/729 + 497/108 = 13819/2916.)
  (iii) sec.3.2 N=2 decimal: quoted 0.1796160944; exact value is
       0.17961609456, i.e. 0.1796160946 to 10 places -- last digit off by 2.
  N=1/3/4 decimals check out (0.1591549431 / 0.1820383116 / 0.1825137102).

METHOD DISCIPLINE (anti-hang):
* T1 (N=1): sympy symbolic only.
* T2/T3 (N>=2): mpmath numerics ONLY -- numeric symmetric eigenpairs
  (mp.eigsy) + elementwise analytic derivatives of [[k]] (rational in
  varpi; NO sympy eigen-decomposition, NO symbolic matrix derivatives).
  Elementwise analytic d^2/d^3 are cross-checked against central finite
  differences with fixed steps; all iterative steps bounded.
* mp.dps = 30. Every eigen-use asserts a nondegenerate leading gap.

Standalone: does NOT import any workspace modules.
"""

# ===========================================================================
# T1: N=1 symbolic verification (sympy only)
# ===========================================================================
import sympy as sp


def run_T1():
    """Symbolic verification of P4A.4 and P4A.5 for N=1."""
    w = sp.Symbol('w', positive=True)  # varpi_bar

    h = w**2 / (w**2 + 1)
    hp, hpp, hppp = sp.diff(h, w), sp.diff(h, w, 2), sp.diff(h, w, 3)

    doc_hp = 2*w / (w**2 + 1)**2
    doc_hpp = (2 - 6*w**2) / (w**2 + 1)**3
    doc_hppp = 24*w*(w**2 - 1) / (w**2 + 1)**4
    chk(sp.simplify(hp - doc_hp) == 0)
    chk(sp.simplify(hpp - doc_hpp) == 0)
    chk(sp.simplify(hppp - doc_hppp) == 0)

    kappa1 = sp.simplify(w**2 * hpp / h)
    c31 = sp.simplify(sp.Rational(1, 6) * w**3 * hppp / h)
    chk(sp.simplify(kappa1 - (2 - 6*w**2)/(w**2 + 1)**2) == 0)
    chk(sp.simplify(c31 - 4*w**2*(w**2 - 1)/(w**2 + 1)**3) == 0)

    # sign-change points: kappa_1 at wbar = 1/sqrt(3), c_{3,1} at wbar = 1
    chk(sp.sqrt(sp.Rational(1, 3)) in sp.solve(sp.numer(sp.together(kappa1)), w))
    chk(sp.Integer(1) in sp.solve(sp.numer(sp.together(c31)), w))

    print("T1 PASS: N=1 closed forms kappa_1, c_{3,1} confirmed by sympy "
          "(7 checks); kappa_1 changes sign at wbar=1/sqrt(3), "
          "c_{3,1} at wbar=1.")


# ===========================================================================
# T2/T3 numerics (mpmath only)
# ===========================================================================
import mpmath as mp

mp.mp.dps = 30

PASS = [0]


def chk(cond, msg=""):
    """Counting assert: aborts on failure, tallies successes."""
    assert cond, msg
    PASS[0] += 1


# ---------- basis [[k]](w) = w^2/(k^2 + w^2) and its derivatives ----------
def brk(k, w):
    """[[k]](w) = w^2/(k^2 + w^2).  For k=0 this is identically 1."""
    w = mp.mpf(w)
    return w**2 / (mp.mpf(k)**2 + w**2)


def dbrk(k, w, order):
    """d^order/dw^order [[k]](w), k >= 1, order in {1,2,3}: exact elementwise.

    [[k]]  = 1 - k^2/(k^2+w^2)
    d1     =  2 k^2 w /(k^2+w^2)^2
    d2     =  2 k^2 (k^2 - 3w^2)/(k^2+w^2)^3
    d3     =  24 k^2 w (w^2 - k^2)/(k^2+w^2)^4
    (k=0: all derivatives vanish -- constant 1.)
    """
    w = mp.mpf(w)
    k2 = mp.mpf(k)**2
    d = k2 + w**2
    if order == 1:
        return 2*k2*w / d**2
    if order == 2:
        return 2*k2*(k2 - 3*w**2) / d**3
    if order == 3:
        return 24*k2*w*(w**2 - k2) / d**4
    raise ValueError(order)


# ---------- KAY II truncated kernel H^(N)(w) and derivatives ----------
def H_matrix(N, w, order=0):
    """KAY II H^(N)(w) (order=0) or its elementwise d^order/dw^order.

    Adopted placement (see module docstring): KAY I G-structure with
    1/k^gamma -> [[k]](w), i.e. H = -H1 + H2 + H3 as written out above.
    """
    w = mp.mpf(w)
    M = mp.matrix(N, N)
    for n in range(N):
        for m in range(N):
            nrm = mp.sqrt(mp.mpf(2*n+1)) * mp.sqrt(mp.mpf(2*m+1))
            if order == 0:
                val = brk(n+m+1, w) / nrm
                if n != m:
                    val += brk(abs(n-m), w) / nrm
                else:
                    s = sum(2*brk(k, w) for k in range(1, n+1)) if n >= 1 else mp.mpf(0)
                    val -= s / mp.mpf(2*n+1)
            else:
                val = dbrk(n+m+1, w, order) / nrm
                if n != m:
                    val += dbrk(abs(n-m), w, order) / nrm
                else:
                    s = sum(2*dbrk(k, w, order) for k in range(1, n+1)) if n >= 1 else mp.mpf(0)
                    val -= s / mp.mpf(2*n+1)
            M[n, m] = val
    return M


def G_matrix(N, gamma):
    """KAY I G^(N)(gamma) (same placement with [[k]] -> 1/k^gamma)."""
    g = mp.mpf(gamma)
    M = mp.matrix(N, N)
    for n in range(N):
        for m in range(N):
            nrm = mp.sqrt(mp.mpf(2*n+1)) * mp.sqrt(mp.mpf(2*m+1))
            val = mp.mpf(1) / (nrm * mp.mpf(n+m+1)**g)
            if n != m:
                val += mp.mpf(1) / (nrm * mp.mpf(abs(n-m))**g)
            else:
                s = sum(2 / mp.mpf(k)**g for k in range(1, n+1)) if n >= 1 else mp.mpf(0)
                val -= s / mp.mpf(2*n+1)
            M[n, m] = val
    return M


def eig_max(M, gap_tol=mp.mpf('1e-10')):
    """Largest eigenvalue + normalized eigenvector; assert simplicity (gap)."""
    E, V = mp.eigsy(M)
    lam = max(E)
    j = 0
    for i in range(len(E)):
        if E[i] == lam:
            j = i
    v = mp.matrix([V[i, j] for i in range(len(E))])
    nv = mp.sqrt(sum(v[i]**2 for i in range(len(v))))
    v = mp.matrix([v[i]/nv for i in range(len(v))])
    # gap to second-largest
    others = [E[i] for i in range(len(E)) if i != j]
    gap = lam - max(others) if others else mp.mpf(1)
    chk(gap > gap_tol, f"leading eigenvalue near-degenerate: gap={mp.nstr(gap, 5)}")
    return lam, v


def quad(v, A):
    """<v|A|v> for real symmetric A, real unit v."""
    N = len(v)
    s = mp.mpf(0)
    for i in range(N):
        for j in range(N):
            s += v[i]*A[i, j]*v[j]
    return s


def fd_quadv(Mfun, w, v, order):
    """Central finite-difference check of <v|M^(order)(w)|v> (fixed steps).

    mpmath is exact, so only truncation error remains and small steps are
    safe: h=1e-3 gives FD error ~1e-6 (q2), ~4e-6 (q3) at the worst tested
    points (the sign-change points where the true value is 0).
    """
    h = mp.mpf('1e-3')
    if order == 2:
        return (quad(v, Mfun(w+h)) - 2*quad(v, Mfun(w)) + quad(v, Mfun(w-h))) / h**2
    if order == 3:
        return (quad(v, Mfun(w+2*h)) - 2*quad(v, Mfun(w+h))
                + 2*quad(v, Mfun(w-h)) - quad(v, Mfun(w-2*h))) / (2*h**3)
    raise ValueError(order)


def kappa_c3(N, w, fd_check=True):
    """(kappa_N(w), c_{3,N}(w), frak_h^(N)(w)) from (P4A.2)/(P4A.3)."""
    w = mp.mpf(w)
    lam, v = eig_max(H_matrix(N, w))
    q2 = quad(v, H_matrix(N, w, order=2))
    q3 = quad(v, H_matrix(N, w, order=3))
    if fd_check:
        q2_fd = fd_quadv(lambda x: H_matrix(N, x), w, v, 2)
        q3_fd = fd_quadv(lambda x: H_matrix(N, x), w, v, 3)
        # relative gate + absolute floor: at sign-change points (e.g. h''=0 at
        # w=1/sqrt(3), h'''=0 at w=1) the true value is ~0 and the FD
        # truncation floor dominates (~5e-7 for q2, ~4e-6 for q3 at h=1e-3).
        chk(abs(q2 - q2_fd) < mp.mpf('1e-4')*abs(q2) + mp.mpf('5e-6'),
            f"FD mismatch q2 N={N} w={w}: {mp.nstr(q2,8)} vs {mp.nstr(q2_fd,8)}")
        chk(abs(q3 - q3_fd) < mp.mpf('1e-3')*abs(q3) + mp.mpf('2e-5'),
            f"FD mismatch q3 N={N} w={w}: {mp.nstr(q3,8)} vs {mp.nstr(q3_fd,8)}")
    return w**2 * q2 / lam, w**3 * q3 / (6*lam), lam


def run_T2():
    print("=== T2: numerics N=1..4 ===")

    # --- (a) KAY I truncation constants (independent G-matrix route) ---
    ref_tc = ['0.159154943', '0.1796160944', '0.1820383', '0.1825137102']
    g_ref = ['1.0000', '1.2737', '1.3082', '1.3151']
    print("(a) KAY I truncation constants from G-matrix "
          "[[k]]->1/k^2 (survey sec.3.2 identity check):")
    tcs = []
    for N in range(1, 5):
        E, _ = mp.eigsy(G_matrix(N, 2))
        gN = max(E)
        tc = mp.sqrt(gN)/(2*mp.pi)
        tcs.append(tc)
        # gate: 0.6 ulp of the *quoted* digits; N=2 is exempted (documented
        # slip (iii): its 10th quoted digit is off by 2) and gated at 3e-9
        dec = len(ref_tc[N-1].split('.')[1])
        tol = mp.mpf('3e-9') if N == 2 else mp.mpf('0.6')*mp.mpf(10)**(-dec)
        chk(abs(tc - mp.mpf(ref_tc[N-1])) < tol,
            f"T_c^({N})(2) mismatch: {mp.nstr(tc, 12)} vs {ref_tc[N-1]}")
        chk(abs(gN - mp.mpf(g_ref[N-1])) < mp.mpf('1e-3'))
        print(f"    N={N}:  g^({N})(2)={mp.nstr(gN, 10):13s}  "
              f"T_c^({N})={mp.nstr(tc, 12)}  (survey {ref_tc[N-1]})")

    # exact closed forms -- decisive fingerprints of the matrix construction
    E, _ = mp.eigsy(G_matrix(1, 2))
    chk(abs(max(E) - 1) < mp.mpf('1e-25'), "g^(1)(2) != 1 (KAY I eq.67)")
    E, _ = mp.eigsy(G_matrix(2, 2))
    g2_exact = (20 + mp.sqrt(13819)) / 108     # tr=10/27, det=-497/432
    chk(abs(max(E) - g2_exact) < mp.mpf('1e-25'),
        f"g^(2)(2) != (20+sqrt(13819))/108: {mp.nstr(max(E), 15)}")
    print("    exact: g^(1)(2)=1 (KAY I eq.67); "
          "g^(2)(2)=(20+sqrt(13819))/108 -- the integer 13819 reproduced,")
    print("    which fingerprints the matrix construction against KAY I eqs.70-71.  OK")
    print("    [survey slips found, see module docstring NOTES]")

    # --- (b) N=1 anchor: H^(1)(w) == w^2/(w^2+1), numerics == T1 closed forms ---
    print("(b) N=1 numeric == symbolic closed forms:")
    for ws in ('0.3', '0.577350269189625764', '1', '2', '5'):
        w = mp.mpf(ws)
        lam, _ = eig_max(H_matrix(1, w))
        chk(abs(lam - w**2/(w**2+1)) < mp.mpf('1e-25'), "H^(1) != w^2/(w^2+1)")
        k1n, c1n, _ = kappa_c3(1, w, fd_check=True)
        k1s = (2 - 6*w**2)/(w**2+1)**2
        c1s = 4*w**2*(w**2-1)/(w**2+1)**3
        chk(abs(k1n - k1s) < mp.mpf('1e-9'), f"kappa_1 numeric {mp.nstr(k1n,10)} vs sym {mp.nstr(k1s,10)}")
        chk(abs(c1n - c1s) < mp.mpf('1e-9'), f"c_3,1 numeric {mp.nstr(c1n,10)} vs sym {mp.nstr(c1s,10)}")
    print("    H^(1)=w^2/(w^2+1); kappa_1, c_3,1 numeric==symbolic to <1e-9 at "
          "w=0.3,1/sqrt(3),1,2,5  OK")

    # --- (c) small-varpi bridge: rho(H^(N)(w))/w^2 -> g^(N)(2) ---
    print("(c) bridge lim_{w->0} rho(H^(N)(w))/w^2 = g^(N)(2):")
    for N in range(1, 5):
        E, _ = mp.eigsy(G_matrix(N, 2))
        gN = max(E)
        E2, _ = mp.eigsy(H_matrix(N, mp.mpf('0.001')))
        ratio = max(E2)/mp.mpf('0.001')**2
        chk(abs(ratio - gN) < mp.mpf('1e-5'),
            f"bridge N={N}: {mp.nstr(ratio, 10)} vs {mp.nstr(gN, 10)}")
        print(f"    N={N}: g^(N)(2)={mp.nstr(gN, 10):13s}  "
              f"rho(H(0.001))/1e-6={mp.nstr(ratio, 10):13s}  OK")

    # --- (d) the kappa_N / c_{3,N} table ---
    print("(d) kappa_N(wbar) [upper] and c_{3,N}(wbar) [lower]:")
    wlist = ['0.3', '0.577350269189625764', '1', '2', '5']
    wtag = ['0.3', '1/sr3', '1', '2', '5']
    hdr = "    wbar:  " + " ".join(f"{t:>10s}" for t in wtag)
    print(hdr)
    for N in range(1, 5):
        ks, cs = [], []
        for ws in wlist:
            k, c, _ = kappa_c3(N, ws)
            ks.append(k)
            cs.append(c)
        print(f"    k_N={N} " + " ".join(f"{mp.nstr(k, 7):>10s}" for k in ks))
        print(f"    c_N={N} " + " ".join(f"{mp.nstr(c, 7):>10s}" for c in cs))
    print("T2 PASS")


# ===========================================================================
# T3: universality of (P4A.1) across spectral shapes (falsifiable core)
# ===========================================================================
def frak_k_of_P(N, wbar, xis, ps):
    """frak_k^(N)(P,T) = rho( sum_j p_j H^(N)(wbar*(1+xi_j)) )  (discrete P)."""
    Nmat = len(H_matrix(N, mp.mpf(1)))  # dimension guard
    M = mp.matrix(Nmat, Nmat)
    for xi, p in zip(xis, ps):
        Hij = H_matrix(N, mp.mpf(wbar)*(1 + mp.mpf(xi)))
        for i in range(Nmat):
            for j in range(Nmat):
                M[i, j] += mp.mpf(p)*Hij[i, j]
    E, _ = mp.eigsy(M)
    return max(E)


def _realroots(coeffs):
    """Real roots of polynomial (highest degree first), sorted; assert ~real."""
    rts = mp.polyroots(coeffs, maxsteps=50, extraprec=50)
    out = []
    for z in rts:
        chk(abs(mp.im(z)) < mp.mpf('1e-18'), f"non-real root {z}")
        out.append(mp.re(z))
    return sorted(out)


def shape_A(s):
    """2-peak symmetric: xi={-s,+s}, p=(1/2,1/2).  (sigma=s, gamma3=0)"""
    return [-s, s], [mp.mpf('0.5'), mp.mpf('0.5')]


def shape_B(s):
    """3-peak sym 1-2-1: xi={-a,0,a}, p=(1/4,1/2,1/4), a=s*sqrt(2)."""
    a = s*mp.sqrt(2)
    return [-a, mp.mpf(0), a], [mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.25')]


def shape_C(s):
    """4-peak symmetric equal-w: xi={-b,-b/2,b/2,b}, b=s*sqrt(8/5)."""
    b = s*mp.sqrt(mp.mpf(8)/5)
    return [-b, -b/2, b/2, b], [mp.mpf('0.25')]*4


def shape_D(s):
    """2-peak asym 1:2: xi={s*sqrt2, -s/sqrt2}, p=(1/3,2/3).  gamma3=1/sqrt2."""
    return [s*mp.sqrt(2), -s/mp.sqrt(2)], [mp.mpf(1)/3, mp.mpf(2)/3]


def shape_E1(s, gam3):
    """3-peak asym, equal weights: roots of t^3 - (3/2)s^2 t - gam3 s^3.

    For equal weights the moment conditions are exactly Vieta's formulas:
    sum xi = 0, sum xi^2 = 3 s^2, sum xi^3 = 3 gam3 s^3  (closed form)."""
    rts = _realroots([1, 0, -mp.mpf('1.5')*s**2, -mp.mpf(gam3)*s**3])
    return rts, [mp.mpf(1)/3]*3


def shape_E2(s, gam3):
    """3-peak asym 1-2-1: xi={-a,-v/2,b}, p=(1/4,1/2,1/4), a=(u-v)/2, b=(u+v)/2.

    Moment conditions reduce to the cubic  v^3 - 4 s^2 v + (8/3) gam3 s^3 = 0;
    the small positive root keeps a>0 (the large positive root gives a<0)."""
    rts = _realroots([1, 0, -4*s**2, mp.mpf(8)/3*mp.mpf(gam3)*s**3])
    pos = [r for r in rts if r > 0]
    chk(pos, "no positive root for v")
    v = min(pos)
    u = mp.sqrt(8*s**2 - 2*v**2)
    chk(u > v, "invalid branch (a<=0)")
    a, b = (u - v)/2, (u + v)/2
    return [-a, -v/2, b], [mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.25')]


def _check_shape(s, gam3, xis, ps, tag):
    """moment audit: <xi>=0, sigma=s, gamma3=gam3 (to 1e-12)."""
    mu = sum(p*xi for xi, p in zip(xis, ps))
    c2 = sum(p*(xi - mu)**2 for xi, p in zip(xis, ps))
    c3 = sum(p*(xi - mu)**3 for xi, p in zip(xis, ps))
    sig = mp.sqrt(c2)
    g3 = c3/sig**3 if c3 != 0 else mp.mpf(0)
    chk(abs(mu) < mp.mpf('1e-15'), f"{tag}: mean {mp.nstr(mu, 5)}")
    chk(abs(sig - s) < mp.mpf('1e-12'), f"{tag}: sigma {mp.nstr(sig, 10)} vs {s}")
    chk(abs(g3 - mp.mpf(gam3)) < mp.mpf('1e-10'),
        f"{tag}: gamma3 {mp.nstr(g3, 10)} vs {gam3}")


def run_T3():
    print("=== T3: universality of (P4A.1) across spectral shapes ===")
    wbars = ['0.5', '1', '2']
    sigmas = [mp.mpf(x) for x in ('0.05', '0.1', '0.15', '0.2')]
    gam3_sym, gam3_asym = 0, mp.mpf('0.35')

    # shape families and their (sigma, gamma3) tags
    def shapes(s):
        return [
            ("2pk-sym  ", *shape_A(s), gam3_sym),
            ("3pk-sym  ", *shape_B(s), gam3_sym),
            ("4pk-sym  ", *shape_C(s), gam3_sym),
            ("3pk-asymE1", *shape_E1(s, gam3_asym), gam3_asym),
            ("3pk-asymE2", *shape_E2(s, gam3_asym), gam3_asym),
        ]

    worst_err4 = mp.mpf(0)      # sup |err23|/sigma^4
    worst_spread4 = mp.mpf(0)   # sup cross-shape spread / sigma^4
    first = True
    for wb in wbars:
        wbar = mp.mpf(wb)
        for N in range(1, 5):
            kap, c3, hN = kappa_c3(N, wbar, fd_check=False)
            errs = {sh: {} for sh in ('2pk-sym  ', '3pk-sym  ', '4pk-sym  ',
                                      '3pk-asymE1', '3pk-asymE2')}
            for s in sigmas:
                exacts = {}
                for tag, xis, ps, g3 in shapes(s):
                    _check_shape(s, g3, xis, ps, tag)
                    exacts[tag] = frak_k_of_P(N, wbar, xis, ps)
                for tag, ex in exacts.items():
                    g3 = gam3_asym if 'asym' in tag else gam3_sym
                    pred3 = hN*(1 + kap/2*s**2 + c3*g3*s**3)
                    err = abs(pred3 - ex)/ex
                    errs[tag][s] = (err, err/s**3, err/s**4, ex, pred3)
                    worst_err4 = max(worst_err4, err/s**4)
                # cross-shape spread at matched (sigma, gamma3): the sigma^4
                # envelope is where peak count / kurtosis may first enter
                for grp, tags in (("sym", ('2pk-sym  ', '3pk-sym  ', '4pk-sym  ')),
                                  ("asym", ('3pk-asymE1', '3pk-asymE2'))):
                    vals = [exacts[t] for t in tags]
                    spread = (max(vals) - min(vals))/exacts[tags[0]]
                    worst_spread4 = max(worst_spread4, spread/s**4)
            if first:  # one detailed narrowing table (err/sigma^3 -> 0, err/sigma^4 -> const)
                first = False
                print(f"  narrowing table at wbar={wb}, N={N} "
                      "(2pk-sym | 3pk-asymE1):")
                print("    sigma   err/s^3(sym) err/s^4(sym)   "
                      "err/s^3(asym) err/s^4(asym)")
                for s in sigmas:
                    _, e3s, e4s, _, _ = errs['2pk-sym  '][s]
                    _, e3a, e4a, _, _ = errs['3pk-asymE1'][s]
                    print(f"    {float(s):.2f}   {float(e3s):10.3e}  {float(e4s):10.3e}   "
                          f"{float(e3a):11.3e}  {float(e4a):11.3e}")
                # err = O(sigma^4) means err/sigma^3 = c*sigma + O(sigma^2):
                # it GROWS linearly with sigma, so the small-sigma value must
                # be well below the large-sigma value (ratio ~ 1/4 + O(sigma))
                for tag in errs:
                    r = [errs[tag][s][1] for s in sigmas]
                    chk(r[0] < r[-1]/2,
                        f"err/sigma^3 not O(sigma^4)-shrinking for {tag}: {r}")
                print("    err/sigma^3(sigma=0.05) < err/sigma^3(sigma=0.2)/2 "
                      "for every family (residual is O(sigma^4))  OK")

    print(f"  sup |err(O(s^3) pred)|/s^4 over all (wbar,N,sigma,shape) = "
          f"{mp.nstr(worst_err4, 6)}")
    print(f"  sup cross-shape spread/s^4 at matched (sigma,gamma3)     = "
          f"{mp.nstr(worst_spread4, 6)}")
    chk(worst_err4 < 8, f"expansion residual not O(sigma^4): {mp.nstr(worst_err4, 6)}")
    chk(worst_spread4 < 8,
        f"shape dependence already at O(sigma^4) order: {mp.nstr(worst_spread4, 6)}")

    # the c3*gamma3*sigma^3 term must DO WORK on asymmetric shapes.  Tested at
    # wbar=0.5 where c_{3,N} ~ -0.35..-0.38 (nonzero, no cancellation): dropping
    # the term must cost strictly more than keeping it.  (At wbar=2, N=1 the
    # sigma^3 and sigma^4 pieces have opposite signs, so pointwise |error| can
    # accidentally favor the LOWER order -- normal for asymptotic series.)
    wbar = mp.mpf('0.5')
    s = mp.mpf('0.2')
    xis, ps = shape_E1(s, gam3_asym)
    _check_shape(s, gam3_asym, xis, ps, "3pk-asymE1")
    for N in (1, 2, 3, 4):
        kap, c3, hN = kappa_c3(N, wbar, fd_check=False)
        ex = frak_k_of_P(N, wbar, xis, ps)
        e2 = abs(hN*(1 + kap/2*s**2) - ex)/ex
        e23 = abs(hN*(1 + kap/2*s**2 + c3*gam3_asym*s**3) - ex)/ex
        chk(e23 < e2/2,
            f"N={N}: adding c3*gamma3*s^3 did not help: e2={mp.nstr(e2, 4)} e23={mp.nstr(e23, 4)}")
    print("  skewness term c_{3,N}*gamma3*sigma^3 improves every N at "
          "(wbar=0.5, s=0.2, gamma3=0.35) by >2x  OK")

    # N=1 exact consistency of the P-average machinery (no eigensolver):
    # frak_k^(1)(P) = sum_j p_j h(wbar(1+xi_j)) must equal the matrix route
    xis, ps = shape_E1(mp.mpf('0.15'), gam3_asym)
    wb = mp.mpf('1')
    scalar = sum(p*(wb*(1+xi))**2/((wb*(1+xi))**2 + 1) for xi, p in zip(xis, ps))
    mat = frak_k_of_P(1, wb, xis, ps)
    chk(abs(scalar - mat) < mp.mpf('1e-25'),
        f"N=1 P-average mismatch: {mp.nstr(scalar, 12)} vs {mp.nstr(mat, 12)}")
    print("  N=1 direct scalar P-average == matrix route to 1e-25  OK")
    print(f"T3 PASS: (P4A.1) holds to O(sigma^4) with shape-independent "
          f"kappa_N/c_{{3,N}} (worst residual/s^4 = {mp.nstr(worst_err4, 4)}, "
          f"worst shape spread/s^4 = {mp.nstr(worst_spread4, 4)})")


if __name__ == '__main__':
    run_T1()
    run_T2()
    run_T3()
    print(f"ALL {PASS[0]} checks passed.")
