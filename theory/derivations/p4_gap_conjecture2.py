r"""
P4-C 试验田 (AI exploration scratch): anatomy of the KAY 2.034 T_c gap
=====================================================================

Target operator (Kiessling-Altshuler-Yuzbashyan I, arXiv:2409.00533v2), the
gamma-model operator on l^2(N_0) (index n from 0), verified upstream (P4-B):

    G(gamma) = -G1 + G2 + G3        (compact, self-adjoint)
    (G1)_{n,m} = [ (1/(2n+1)) sum_{k=1}^{n} 2/k^gamma ] delta_{n,m}   (diag >=0)
    (G2)_{n,m} = (1-delta_{n,m}) / ( sqrt(2n+1) |n-m|^gamma sqrt(2m+1) )
    (G3)_{n,m} = 1 / ( sqrt(2n+1) (n+m+1)^gamma sqrt(2m+1) )          (>0, PSD)

At gamma=2:  g(2) := rho(G(2)) ~ 1.3181405  (P4-B enclosure width 3.9e-8),
             C_inf = (1/2pi) sqrt(g(2)) = 0.1827262477... (sharp constant).

KAY's PROVEN analytic upper-bound chain (survey Sec.2.3-2.4), gamma=2, eps=0.65:
    rho(G) <= rho(G2+G3)          [loss (1): drop the negative -G1]
           <= rho(G2)+rho(G3)     [loss (2): subadditivity, Lemma 9]
           <= 2*Prop5 + 1         [loss (3): Cauchy-Schwarz + zeta, Props 5&6]
    Prop5 = sqrt( (2^{1+eps}-1) zeta(1+2gamma-eps) zeta(1+eps) ),  eps=min{gamma,0.65}.

This script (HONEST NUMERICAL SCOUTING, not a proof of anything new):

  PROBE A1  anatomize the three losses in g-units -> which one dominates?
            (a prior scratch said loss (3) dominates ~90%, CONTRADICTING the
             survey Sec.2.4 which claims loss (1) is largest -- we adjudicate).

  PROBE A2  a one-parameter geometric trial vector v_n = t^n gives
            * a rigorous Rayleigh LOWER bound  mu = <v|G|v>/<v|v> <= g,
            * a rigorous Kato-Temple UPPER bound g <= mu + eta^2/(mu-alpha),
              where alpha is a rigorous upper bound on the SECOND eigenvalue
              of G.  Point: show the factor ~2.03 is an artifact of the loose
              (2)+(3) chain, not intrinsic.

  PROBE B   (secondary) finite-lambda Conjecture 2 + the "Semenok crossing":
            reported honestly, including what could NOT be reconstructed.

Numerics discipline: numpy float64 + eigvalsh/eigh for all heavy linear
algebra (symmetric); mpmath dps=30 ONLY for final one-matvec certificates and
one Cholesky.  Never generic eig, never full mpmath eigendecomposition.

Run:  python theory/derivations/p4_gap_conjecture2.py
"""
import time
import numpy as np
import mpmath as mp

mp.mp.dps = 30

CHECKS = []


def _check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f"  {detail}" if detail else ""))
    assert ok, f"{name} FAILED: {detail}"


def _stage_summary(k, note=""):
    n_pass = sum(1 for _, ok in CHECKS if ok)
    n_total = len(CHECKS)
    print(f"\nStage {k}: {n_pass}/{n_total} checks passed (cumulative){note}\n")


# ============================================================================
# Operator construction (inlined & self-contained; copied from the verified
# P4-B builder p4_C_infinity_enclosure.py -- NO import).
# ============================================================================

def build_blocks_np(N, gamma=2.0):
    """float64 blocks (-G1 diag, G2, G3) and G = -G1+G2+G3 as N x N arrays."""
    gamma = float(gamma)
    n = np.arange(N)
    sq = np.sqrt(2.0 * n + 1.0)
    g1d = np.zeros(N)
    running = 0.0
    for k in range(N):
        if k >= 1:
            running += 2.0 / (k ** gamma)
        g1d[k] = running / (2 * k + 1)
    idx = np.arange(N)
    Dm = np.abs(idx[:, None] - idx[None, :])
    G2 = np.where(Dm > 0,
                  1.0 / (sq[:, None] * np.where(Dm > 0, Dm, 1) ** gamma * sq[None, :]),
                  0.0)
    Sm = (idx[:, None] + idx[None, :] + 1.0)
    G3 = 1.0 / (sq[:, None] * Sm ** gamma * sq[None, :])
    G1 = np.diag(g1d)
    return G1, G2, G3, (-G1 + G2 + G3)


def build_G_mp(N, gamma=2):
    """Exact mpmath G = -G1 + G2 + G3 (dps=30). O(N^2), seconds for N<=400."""
    gamma = mp.mpf(gamma)
    sq = [mp.sqrt(2 * n + 1) for n in range(N)]
    G1neg = []
    run = mp.mpf(0)
    for n in range(N):
        if n >= 1:
            run += 2 / mp.mpf(n) ** gamma
        G1neg.append(-run / (2 * n + 1))
    G = mp.matrix(N)
    for n in range(N):
        sqn = sq[n]
        for m in range(N):
            if n == m:
                G[n, m] = G1neg[n] + 1 / (sqn * mp.mpf(n + m + 1) ** gamma * sq[m])
            else:
                G[n, m] = (1 / (sqn * mp.mpf(abs(n - m)) ** gamma * sq[m])
                           + 1 / (sqn * mp.mpf(n + m + 1) ** gamma * sq[m]))
    return G


def prop5_rhs(gamma):
    """Prop 5 RHS = sqrt((2^{1+eps}-1) zeta(1+2gamma-eps) zeta(1+eps)), eps=min{gamma,0.65}."""
    gamma = mp.mpf(gamma)
    eps = min(gamma, mp.mpf("0.65"))
    return mp.sqrt((mp.mpf(2) ** (1 + eps) - 1) * mp.zeta(1 + 2 * gamma - eps) * mp.zeta(1 + eps)), eps


# ============================================================================
# Rigorous tail certificates for the SECOND eigenvalue (copied/adapted from
# the P4-B enclosure code; used to bound rho(BB^T) and rho(D_tail)).
# ============================================================================

def tail_d_bound(N):
    """Gershgorin bound d with rho(D_tail) <= d (tail block n,m >= N), gamma=2.

    |G1_nn| <= (pi^2/3)/(2N+1); G2 row <= (pi^2/3)/(2N+1); G3 row <= 1/(N(2N+1)).
    (2*zeta(2) = pi^2/3 exact; sum_{j>=J} j^-2 <= 1/(J-1).)
    """
    zeta2 = mp.pi ** 2 / 3
    return (zeta2 + zeta2) / (2 * N + 1) + mp.mpf(1) / (N * (2 * N + 1))


def beta_bbT_bound(N, w, M_cut):
    """Rigorous beta >= rho(B B^T), B = head(N) x tail coupling, gamma=2.

    Collatz-Wielandt on nonnegative S=BB^T with positive w:
      rho(S) <= max_n (S w)_n / w_n,  (S w)_n = sum_{m>=N} B_{nm} c_m,
      c_m = sum_{n'<N} B_{n'm} w_{n'}.
    Partial float64 sums m in [N, M_cut) (all positive; 1e-9 inflated) plus a
    closed-form remainder for m >= M_cut >= 4N:
      B_{nm} <= 5/(m^2 sqrt(2n+1)); sum_{m>=M} m^-4 <= 1/(3(M-1)^3).
    """
    assert M_cut >= 4 * N, "need M_cut >= 4N for the tail bound"
    inv_sqrt = 1.0 / np.sqrt(2.0 * np.arange(float(N)) + 1.0)
    w_np = np.array([float(x) for x in w])
    Sw = np.zeros(N)
    m = N
    while m < M_cut:
        sqm = np.sqrt(2.0 * m + 1.0)
        mdm = m - np.arange(float(N))
        Bcol = (1.0 / (mdm * mdm) + 1.0 / ((m + np.arange(float(N)) + 1.0) ** 2)) * inv_sqrt / sqm
        Sw += Bcol * float(np.dot(Bcol, w_np))
        m += 1
    M = mp.mpf(M_cut)
    W1 = sum(w)
    tail_max = 25 * W1 / (mp.sqrt(2 * M + 1) * 3 * (M - 1) ** 3)
    beta = mp.mpf(0)
    for n in range(N):
        Swn = mp.mpf(float(Sw[n])) * (1 + mp.mpf("1e-9"))
        beta = max(beta, (Swn + tail_max / mp.sqrt(2 * n + 1)) / w[n])
    return beta * (1 + mp.mpf("1e-9"))


# reference true value (P4-B rigorous enclosure midpoint)
G_TRUE_REF = mp.mpf("1.3181405")
C_INF_REF = mp.mpf("0.1827262477")


# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    T_START = time.time()

    # ------------------------------------------------------------------
    # Stage 0: operator self-tests
    # ------------------------------------------------------------------
    print("=" * 70)
    print("Stage 0: operator construction self-tests")
    print("=" * 70)

    # g^(1)(2) = 1 exactly (KAY I eq.67)
    _, _, _, G1x1 = build_blocks_np(1, 2.0)
    _check("g^(1)(2) = 1 exactly (float64)", abs(G1x1[0, 0] - 1.0) < 1e-14,
           f"got {G1x1[0,0]:.15f}")

    # 4x4 explicit entries (survey Sec.1.4): G[0,0]=1, G[0,1]=(1/2^2+1)/sqrt(3)
    _, _, _, G4 = build_blocks_np(4, 2.0)
    _check("G[0,0]=1 (4x4)", abs(G4[0, 0] - 1.0) < 1e-14, f"got {G4[0,0]}")
    _check("G[0,1]=(1/4+1)/sqrt3", abs(G4[0, 1] - (1.25) / np.sqrt(3)) < 1e-14,
           f"got {G4[0,1]:.12f}")
    _check("G symmetric (4x4)", np.allclose(G4, G4.T), "")

    # g_true from truncations converges up to ~1.3181405 (Theorem 2)
    g_prev = -1.0
    g_seq = {}
    for Nt in [50, 100, 200, 400]:
        _, _, _, GN = build_blocks_np(Nt, 2.0)
        g_seq[Nt] = float(np.linalg.eigvalsh(GN)[-1])
        _check(f"monotone-up g^({Nt}) > g^(prev)", g_seq[Nt] > g_prev - 1e-15,
               f"g^({Nt})={g_seq[Nt]:.9f}")
        g_prev = g_seq[Nt]
    _check("g^(400)(2) matches P4-B reference 1.3181405 (5e-7)",
           abs(g_seq[400] - float(G_TRUE_REF)) < 5e-7,
           f"g^(400)={g_seq[400]:.9f}, ref={float(G_TRUE_REF)}")
    print(f"    g^(N)(2): " + ", ".join(f"N={k}:{v:.9f}" for k, v in g_seq.items()))
    _stage_summary(0)

    # ------------------------------------------------------------------
    # Stage 1 (PROBE A1): anatomy of the three losses, in g-units
    # ------------------------------------------------------------------
    print("=" * 70)
    print("Stage 1 (Probe A1): loss decomposition of the 2.03 gap, in g-units")
    print("=" * 70)

    prop5, eps = prop5_rhs(2.0)
    chain = 2 * prop5 + 1  # = Prop5 + Prop6
    print(f"    eps(2)={float(eps)}, Prop5={mp.nstr(prop5,10)}, "
          f"Prop6=1+Prop5={mp.nstr(1+prop5,10)}, 2*Prop5+1={mp.nstr(chain,10)}")

    # N-convergence table of the three spectral radii
    ladder = {}
    for Nt in [100, 200, 400]:
        G1b, G2b, G3b, GN = build_blocks_np(Nt, 2.0)
        g_true = float(np.linalg.eigvalsh(GN)[-1])
        r23 = float(np.linalg.eigvalsh(G2b + G3b)[-1])
        r2 = float(np.linalg.eigvalsh(G2b)[-1])
        r3 = float(np.linalg.eigvalsh(G3b)[-1])
        ladder[Nt] = (g_true, r23, r2, r3)
        print(f"    N={Nt:4d}: g_true={g_true:.7f}  rho(G2+G3)={r23:.7f}  "
              f"rho(G2)={r2:.7f}  rho(G3)={r3:.7f}  rho(G2)+rho(G3)={r2+r3:.7f}")

    # N-convergence check: values stable between N=200 and N=400 (< 1e-5)
    for j, lab in enumerate(["g_true", "rho(G2+G3)", "rho(G2)", "rho(G3)"]):
        d = abs(ladder[400][j] - ladder[200][j])
        _check(f"N-converged: {lab} (|400-200|<1e-5)", d < 1e-5, f"diff={d:.2e}")

    g_true, r23, r2, r3 = ladder[400]
    r2p3 = r2 + r3
    chain_f = float(chain)
    loss1 = r23 - g_true          # drop -G1
    loss2 = r2p3 - r23            # subadditivity
    loss3 = chain_f - r2p3        # Cauchy-Schwarz + zeta
    total = chain_f - g_true
    s1, s2, s3 = loss1 / total, loss2 / total, loss3 / total
    print()
    print(f"    LADDER (g-units, N=400):")
    print(f"      g_true              = {g_true:.7f}")
    print(f"      after (1) rho(G2+G3)= {r23:.7f}   [loss (1) = {loss1:.7f}]")
    print(f"      after (2) rho2+rho3 = {r2p3:.7f}   [loss (2) = {loss2:.7f}]")
    print(f"      after (3) 2*Prop5+1 = {chain_f:.7f}   [loss (3) = {loss3:.7f}]")
    print(f"      total gap           = {total:.7f}")
    print(f"    SHARES:  (1) drop -G1     = {100*s1:6.2f}%")
    print(f"             (2) subadditivity= {100*s2:6.2f}%")
    print(f"             (3) C-S + zeta   = {100*s3:6.2f}%")

    _check("shares sum to 100%", abs(s1 + s2 + s3 - 1.0) < 1e-12)
    _check("loss ordering (3) >> (2) > (1)", loss3 > loss2 > loss1)
    _check("loss (3) dominates (> 80% of gap)", s3 > 0.80, f"share={100*s3:.2f}%")
    _check("loss (1) is the SMALLEST (< 5%)", s1 < 0.05, f"share={100*s1:.2f}%")

    # T_c-unit cross-check: the overall factor sqrt(chain/g_true)
    tc_factor = float(mp.sqrt(chain / mp.mpf(g_true)))
    print(f"\n    T_c chain factor sqrt((2*Prop5+1)/g_true) = {tc_factor:.5f}")
    print(f"    (= T_c*(2)/C_inf; T_c*(2)=(1/2pi)sqrt(2*Prop5+1)="
          f"{float(mp.sqrt(chain)/(2*mp.pi)):.7f}; survey's 2.032 is vs T_c^(4).)")
    _check("chain factor ~ 2.03 (in [2.02,2.05])", 2.02 < tc_factor < 2.05,
           f"factor={tc_factor:.5f}")

    print()
    print("    *** VERDICT (Probe A1):")
    print(f"        Loss (3) [Cauchy-Schwarz + zeta, Props 5&6] DOMINATES:"
          f" {100*s3:.1f}% of the gap.")
    print(f"        Loss (1) [dropping -G1] is the SMALLEST: {100*s1:.1f}%.")
    print("        => Survey Sec.2.4's claim that loss (1) is largest is WRONG.")
    print("        The prior scratch (1~0.134, 2~0.284, 3~3.69) is CONFIRMED.")
    _stage_summary(1)

    # ------------------------------------------------------------------
    # Stage 2 (PROBE A2): geometric trial vector -> rigorous 2-sided bound
    # ------------------------------------------------------------------
    print("=" * 70)
    print("Stage 2 (Probe A2): geometric trial vector v_n = t^n")
    print("  rigorous lower bound (Rayleigh) + rigorous upper bound (Kato-Temple)")
    print("=" * 70)

    # ---- (i) rigorous alpha >= lambda_2(G), needed for Kato-Temple ----
    # Weyl + operator order:  G = X + Y,  X = A(+)0  (head block A, zero tail),
    #   Y = [[0,B],[B^T,D]].   lambda_2(G) <= lambda_2(X) + ||Y||_op
    #   = lambda_2(A) + (sqrt(rho(BB^T)) + rho(D)) <= lambda_2(A) + sqrt(beta)+d.
    # lambda_2(A) is certified by a Cholesky-deflation:  with P = I - u u^T
    #   (u = numerical top eigenvector), lambda_2(A) <= max_{x _|_ u} <x,Ax>
    #   = lambda_max(P A P);  we certify lambda_max(PAP) <= c by showing
    #   c I - P A P is SPD (mpmath Cholesky).  Valid even for approximate u.
    Na = 200
    print(f"  [alpha certificate] head size Na={Na}")
    G1b, G2b, G3b, A_np = build_blocks_np(Na, 2.0)
    wA, vA = np.linalg.eigh(A_np)
    lam2A_num = float(wA[-2])
    u = vA[:, -1]
    c = mp.mpf(lam2A_num) + mp.mpf("1e-7")   # target upper bound on lambda_2(A)
    t_c0 = time.time()
    G_mp_a = build_G_mp(Na, 2)
    umb = mp.matrix([mp.mpf(float(x)) for x in u])
    Au = G_mp_a * umb
    uAu = sum(umb[i] * Au[i] for i in range(Na))
    M = mp.matrix(Na)
    for i in range(Na):
        for jj in range(Na):
            pap = G_mp_a[i, jj] - umb[i] * Au[jj] - Au[i] * umb[jj] + umb[i] * umb[jj] * uAu
            M[i, jj] = (c if i == jj else mp.mpf(0)) - pap
    chol_ok = True
    try:
        mp.cholesky(M)
    except Exception as e:  # noqa
        chol_ok = False
    t_c1 = time.time()
    _check("Cholesky certifies lambda_2(A_Na) <= c (SPD of cI-PAP)", chol_ok,
           f"c={float(c):.3e}, {t_c1-t_c0:.1f}s")
    lam2A_up = c
    d_tail = tail_d_bound(Na)
    beta = beta_bbT_bound(Na, [mp.mpf(1)] * Na, 4 * Na)
    alpha = lam2A_up + mp.sqrt(beta) + d_tail
    print(f"      lambda_2(A_Na) <= {float(lam2A_up):.4e}   (num {lam2A_num:.4e})")
    print(f"      sqrt(rho(BB^T)) <= {float(mp.sqrt(beta)):.4e}, rho(D) <= {float(d_tail):.4e}")
    print(f"      => alpha (rigorous upper bound on lambda_2(G)) = {float(alpha):.4e}")

    # ---- (ii) optimize t in float64 (fast), then certify at t* in mpmath ----
    NG = 400
    Nv = 120                       # support of v; t^120 < 1e-40 for t<0.45
    _, _, _, G_NG_np = build_blocks_np(NG, 2.0)

    def kt_float(t):
        v = np.zeros(NG)
        v[:Nv] = t ** np.arange(Nv)
        Gv = G_NG_np @ v
        vv = v @ v
        mu = (v @ Gv) / vv
        eta2 = (Gv @ Gv) / vv - mu * mu
        return mu, eta2, mu + eta2 / (mu - float(alpha))

    ts = np.linspace(0.34, 0.46, 121)
    kts = [kt_float(t)[2] for t in ts]
    t_star = float(ts[int(np.argmin(kts))])
    print(f"\n  [trial optimization, float64]  best t* = {t_star:.4f}")

    # ---- (iii) rigorous mu and eta^2 at t* (mpmath, one matvec) ----
    # mu = <v|G|v>/<v|v> EXACT (v supported on n<Nv<=NG): rigorous LOWER bound.
    # eta^2 = ||Gv||^2/<v|v> - mu^2.  ||Gv||^2 = sum_{m>=0}(Gv)_m^2; we sum
    # m<NG exactly and ADD a rigorous tail (m>=NG) so eta^2_up >= eta^2_true,
    # keeping the Kato-Temple value a valid UPPER bound.
    t_m0 = time.time()
    G_mp = build_G_mp(NG, 2)
    v = mp.matrix(NG, 1)
    tt = mp.mpf(t_star)
    for n in range(Nv):
        v[n, 0] = tt ** n
    vv = sum(v[i, 0] ** 2 for i in range(NG))
    Gv = G_mp * v
    xGx = sum(v[i, 0] * Gv[i, 0] for i in range(NG))
    mu = xGx / vv                                  # rigorous lower bound
    normGv2 = sum(Gv[i, 0] ** 2 for i in range(NG))
    # tail: for m>=NG>=4Nv, (Gv)_m <= (5/(m^2 sqrt(2m+1))) S, S = sum t^n/sqrt(2n+1) <= 1/(1-t)
    S = mp.mpf(1) / (1 - tt)
    tail_eta = 25 * S ** 2 / (4 * (mp.mpf(NG) - 1) ** 4)   # sum_{m>=NG} 25 S^2/m^5
    eta2_up = (normGv2 + tail_eta) / vv - mu ** 2
    assert alpha < mu, "Kato-Temple needs alpha < mu"
    kt_upper = mu + eta2_up / (mu - alpha)
    t_m1 = time.time()

    print(f"  [rigorous certificate at t*={t_star:.4f}, NG={NG}, {t_m1-t_m0:.1f}s]")
    print(f"      mu   (rigorous LOWER bound on g) = {mp.nstr(mu, 12)}")
    print(f"      eta^2 (rigorous upper bound)     = {mp.nstr(eta2_up, 6)}  "
          f"(tail added {mp.nstr(tail_eta,3)})")
    print(f"      Kato-Temple UPPER bound on g     = {mp.nstr(kt_upper, 12)}")
    print(f"      => g(2) in [ {mp.nstr(mu,12)} , {mp.nstr(kt_upper,12)} ]")

    width = kt_upper - mu
    C_lo = mp.sqrt(mu) / (2 * mp.pi)
    C_up = mp.sqrt(kt_upper) / (2 * mp.pi)
    tc_up_factor = mp.sqrt(kt_upper / G_TRUE_REF)   # T_c(upper)/C_inf
    print(f"      enclosure width on g = {mp.nstr(width, 4)}")
    print(f"      C_inf in [ {mp.nstr(C_lo,10)} , {mp.nstr(C_up,10)} ]  (true {C_INF_REF})")
    print(f"      T_c upper-bound factor vs C_inf  = {mp.nstr(tc_up_factor, 8)}")
    print(f"      (KAY analytic chain factor was   = {tc_factor:.5f})")

    # self-tests: valid two-sided bracket around the P4-B reference
    _check("mu is a valid LOWER bound (mu <= g_ref + tiny)",
           mu <= G_TRUE_REF + mp.mpf("1e-6"), f"mu={mp.nstr(mu,12)}")
    _check("Kato-Temple is a valid UPPER bound (kt >= g_ref)",
           kt_upper >= G_TRUE_REF - mp.mpf("1e-6"), f"kt={mp.nstr(kt_upper,12)}")
    _check("bracket actually contains the reference g(2)",
           mu <= G_TRUE_REF <= kt_upper + mp.mpf("1e-9"))
    _check("upper bound beats the analytic chain by >100x",
           (chain - G_TRUE_REF) / (kt_upper - G_TRUE_REF) > 100,
           f"chain excess={mp.nstr(chain-G_TRUE_REF,4)}, KT excess={mp.nstr(kt_upper-G_TRUE_REF,4)}")
    _check("T_c upper factor < 1.001 (vs chain 2.03)", tc_up_factor < mp.mpf("1.001"),
           f"factor={mp.nstr(tc_up_factor,8)}")

    print()
    print("    *** VERDICT (Probe A2):")
    print(f"        A rigorous ONE-parameter trial vector v_n=t^n (t*={t_star:.3f})")
    print(f"        + Kato-Temple gives g(2) <= {mp.nstr(kt_upper,10)}, i.e. a T_c")
    print(f"        upper bound only {mp.nstr(tc_up_factor,7)}x C_inf -- vs the KAY")
    print(f"        chain's {tc_factor:.3f}x.  The ~2.03 factor is an ARTIFACT of the")
    print("        subadditivity+zeta steps (2)+(3), not intrinsic to the operator.")
    print("        Auxiliary quantity required: a rigorous upper bound alpha on")
    print("        lambda_2(G) with alpha < mu (here alpha=%.3e << mu)." % float(alpha))
    _stage_summary(2)

    # ------------------------------------------------------------------
    # Stage 3 (PROBE B): finite-lambda Conjecture 2 + Semenok crossing
    # ------------------------------------------------------------------
    print("=" * 70)
    print("Stage 3 (Probe B): Conjecture 2 (finite-lambda) & Semenok crossing")
    print("=" * 70)

    # B1: the ONLY finite-"size" analogue we can compute rigorously from the
    # material at hand is the KAY truncation sequence T_c^(N)(gamma), which
    # (KAY Thm 2, eqs.15-16) converges UPWARD to T_c(gamma).  This is a genuine
    # "approach the limit from below, saturating as the cutoff grows" structure
    # -- but it is finite-N, NOT the finite-lambda (finite coupling) gap
    # equation of Conjecture 2, whose kernel is NOT in the upstream material.
    print("  B1 (finite-lambda Conjecture 2):")
    print("     The exact finite-lambda (finite-coupling) gap-equation kernel is")
    print("     NOT present in the upstream survey; reconstructing it would be a")
    print("     guess, so we do NOT fabricate a 'finite-lambda T_c'.  What we CAN")
    print("     show rigorously is the KAY truncation sequence T_c^(N)(2), which")
    print("     converges UPWARD to T_c(2) (Thm 2): ratio < 1, -> 1 as N grows.")
    Tc_inf = float(mp.sqrt(G_TRUE_REF)) / (2 * np.pi)
    ratios = []
    for Nt in [1, 2, 4, 10, 50, 200, 400]:
        _, _, _, GN = build_blocks_np(max(Nt, 1), 2.0)
        gN = float(np.linalg.eigvalsh(GN)[-1])
        TcN = np.sqrt(gN) / (2 * np.pi)
        ratios.append(TcN / Tc_inf)
        print(f"       N={Nt:4d}: T_c^(N)(2)={TcN:.9f}  ratio T_c^(N)/T_c={TcN/Tc_inf:.7f}")
    _check("B1: truncation ratios strictly increasing", all(np.diff(ratios) > 0))
    _check("B1: all ratios < 1 (bound approached from below)", all(r < 1 + 1e-9 for r in ratios))
    _check("B1: ratio -> 1 (N=400 within 1e-6)", abs(ratios[-1] - 1.0) < 1e-6)

    # B2: the "Semenok crossing 0.31748 vs 0.3175".  Scan T_c(gamma) =
    # (1/2pi) g(gamma)^{1/gamma}.  It is monotone DECREASING in gamma (no
    # interior extremum): -> 1/2pi as gamma->inf, blows up as gamma->0.
    print("\n  B2 (Semenok crossing ~0.31748 vs 0.3175):")
    tc_of_gamma = {}
    for gm in [0.75, 0.85, 1.0, 1.5, 2.0, 3.0, 5.0]:
        _, _, _, GN = build_blocks_np(300, gm)
        gg = float(np.linalg.eigvalsh(GN)[-1])
        tc_of_gamma[gm] = (1 / (2 * np.pi)) * gg ** (1.0 / gm)
        print(f"       gamma={gm:4.2f}: g={gg:.6f}  T_c(gamma)=(1/2pi)g^(1/gamma)={tc_of_gamma[gm]:.6f}")
    # locate gamma with T_c = 0.3175
    lo, hi = 0.75, 1.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        _, _, _, GN = build_blocks_np(300, mid)
        gg = float(np.linalg.eigvalsh(GN)[-1])
        val = (1 / (2 * np.pi)) * gg ** (1.0 / mid)
        if val > 0.3175:
            lo = mid
        else:
            hi = mid
    gamma_cross = 0.5 * (lo + hi)
    print(f"       T_c(gamma)=0.3175 is hit at gamma ~ {gamma_cross:.4f} (interpolated,"
          f" not a distinguished point).")
    print(f"       Note 1/pi = {1/np.pi:.6f} (close to 0.31748 but NOT equal).")
    print("       VERDICT B2: T_c(gamma) from the gamma-model is MONOTONE in gamma;")
    print("       0.31748/0.3175 is not a natural crossing/extremum of it.  Its exact")
    print("       'Semenok' definition is NOT reconstructable from upstream material")
    print("       -> reported honestly, NOT fabricated.")
    _check("B2: T_c(gamma) monotone decreasing on tested grid",
           all(tc_of_gamma[a] > tc_of_gamma[b]
               for a, b in zip([0.75, 0.85, 1.0, 1.5, 2.0, 3.0], [0.85, 1.0, 1.5, 2.0, 3.0, 5.0])))
    _stage_summary(3)

    # ------------------------------------------------------------------
    # Final synthesis
    # ------------------------------------------------------------------
    print("=" * 70)
    print("SYNTHESIS")
    print("=" * 70)
    print(f"  A1  loss shares (g-units): (1) drop -G1 = {100*s1:.1f}%,"
          f" (2) subadd = {100*s2:.1f}%, (3) C-S+zeta = {100*s3:.1f}%")
    print(f"      -> loss (3) dominates; survey Sec.2.4 (claims (1) largest) is WRONG.")
    print(f"  A2  rigorous g(2) in [{mp.nstr(mu,11)}, {mp.nstr(kt_upper,11)}]"
          f" via v_n=t^n + Kato-Temple")
    print(f"      -> T_c upper factor {mp.nstr(tc_up_factor,7)}x vs chain {tc_factor:.3f}x"
          f" (2.03 is an artifact).")
    print("  B   finite-lambda kernel & Semenok crossing NOT reconstructable from")
    print("      upstream material -> reported honestly, not fabricated; finite-N")
    print("      monotone convergence shown as the computable analogue.")
    n_pass = sum(1 for _, ok in CHECKS if ok)
    n_total = len(CHECKS)
    print("=" * 70)
    print(f"FINAL: {n_pass}/{n_total} checks passed   (total {time.time()-T_START:.1f}s)")
    print("=" * 70)
