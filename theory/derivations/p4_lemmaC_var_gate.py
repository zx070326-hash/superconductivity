"""
p4_lemmaC_var_gate.py  —  Gate for the Lemma C variance reduction.

Reproduces every quantitative claim of  theory/derivations/P4_lemmaC_var_bound.md :
the exact identity  Var_v(B) = Var_v(M)  with  M := A - B = O[f - f^2],
the two operator-norm bounds on Var_v(M), the Perron/Schur upper bound on
lambda_max(M) via  M <= M_od + D^+ , the two scalar sufficient reductions
(SUFF and the q-free QF), the certificate-level "ceiling" (the terminal
inequality holds with the ANALYTIC gap lower bound gamma_lo and true spectral
data), and the cone-coercivity capture ratios that localise the single residual.

Run:   OMP_NUM_THREADS=1 python p4_lemmaC_var_gate.py
Single thread, N=400, ~11 grid points; well under 10 min.

Builder + gap machinery are imported dirname-relative (no hardcoded worktree path):
  p4_concavity_q_route.py   -> O_ker, kernels, top_eigen
  p4_concavity_gap2_gates.py-> h_lo_trial, lam2_hi, THETA_STAR
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p4_concavity_q_route import O_ker, kernels, top_eigen
from p4_concavity_gap2_gates import h_lo_trial, lam2_hi, THETA_STAR

# ------------------------------------------------------------------ ck counter
_NP = 0
_NF = 0


def ck(cond, msg):
    global _NP, _NF
    if cond:
        _NP += 1
    else:
        _NF += 1
        print("  FAIL:", msg)
    return cond


# ------------------------------------------------------------------ helpers
def kernel_array(x, KX, which):
    """values phi(k), k=0..KX-1 (phi(0) never used by O_ker)."""
    k = np.arange(KX, dtype=float)
    d = k * k + x
    f = x / d
    if which == "f2":
        a = f * f
    elif which == "f3":
        a = f * f * f
    elif which == "f-f2":
        a = f - f * f
    elif which == "f-f3":
        a = f - f * f * f
    else:
        raise ValueError(which)
    a[0] = 0.0
    return a


def kL3(v, karr, N, KX):
    """cone-coercivity lower bound  sum_{n,m} phi(n+m+1) u_n u_m v_n v_m  (P4_cone_coercivity Thm 1)."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    w = u * v
    S = np.minimum(n[:, None] + n[None, :] + 1, KX - 1)
    return float((karr[S] * np.outer(w, w)).sum())


def wschur(Mabs_dom, psi):
    """Schur-test upper bound on lambda_max of a nonneg symmetric matrix: max_n (Mdom psi)_n / psi_n."""
    return float(((Mabs_dom @ psi) / psi).max())


def perron_schur(Mdom, iters=400):
    """near-Perron weight by power iteration, then Schur bound (= lambda_max(Mdom) in the limit)."""
    N = Mdom.shape[0]
    psi = np.ones(N)
    for _ in range(iters):
        psi = Mdom @ psi
        psi = np.abs(psi) / (np.linalg.norm(psi) + 1e-300) + 1e-300
    return wschur(Mdom, psi)


# ------------------------------------------------------------------ main grid
VARPI = [0.1, 0.3, 0.5, 0.7, 1.0, 1.3, 1.5, 2.0, 3.0, 6.0, 10.0, 20.0, 40.0]
N = 400

print("=" * 78)
print("P4 Lemma C variance-reduction gate   (N=%d, single thread)" % N)
print("=" * 78)
hdr = ("  wpi     VarB    opm/Var  BD/Var  Ldom/L   Su/L   VarM=VarB "
       " capMq  capM3  budget  QFok")
print(hdr)

opm_ratios = []
bd_ratios = []
ldom_ratios = []
su_ratios = []
budgets = []
cap_min = 1.0
small = {}

for w in VARPI:
    x = w * w
    KX = 8 * N
    ker = kernels(x, KX=KX)
    A = O_ker(N, ker["f"])
    B = O_ker(N, ker["f2"])
    B3 = O_ker(N, ker["f3"])
    M = A - B
    h, v, l2 = top_eigen(A)
    q = float(v @ (B @ v))
    m3 = float(v @ (B3 @ v))
    Bv = B @ v
    VarB = float(Bv @ Bv - q * q)

    # ---- (I) exact identity  Var_v(B) = Var_v(M)
    Mv = M @ v
    mm = float(v @ Mv)
    VarM = float(Mv @ Mv - mm * mm)
    ck(abs(VarB - VarM) < 1e-9, f"[{w}] VarB!=VarM ({abs(VarB-VarM):.2e})")
    ck(abs(mm - (h - q)) < 1e-9, f"[{w}] <v,Mv> != h-q")
    ck(VarB > 0, f"[{w}] VarB not positive")

    # ---- (II) operator-norm bounds on Var_v(M)
    eM = np.linalg.eigvalsh(M)
    Lmax, Lmin = float(eM[-1]), float(eM[0])
    ck(Lmax >= abs(Lmin) - 1e-12, f"[{w}] top eigenvalue does not dominate |Lmin|")  # ||M||_2 = Lmax
    hq = h - q
    opm = Lmax ** 2 - hq ** 2                       # op-norm-minus upper bound
    BD = (Lmax - hq) * (hq - Lmin)                  # Bhatia-Davis upper bound
    ck(opm >= VarB - 1e-12, f"[{w}] op-norm-minus not an upper bound")
    ck(BD >= VarB - 1e-12, f"[{w}] Bhatia-Davis not an upper bound")
    ck(BD <= opm + 1e-12, f"[{w}] BD not tighter than opm")
    opm_ratios.append(opm / VarB)
    bd_ratios.append(BD / VarB)

    # ---- (III) analytic upper bound on lambda_max(M):  M <= M_od + D^+
    di = np.diag(M).copy()
    Mdom = M.copy()
    np.fill_diagonal(Mdom, 0.0)
    Mdom = np.abs(Mdom) + np.diag(np.maximum(di, 0.0))   # nonneg symmetric, entrywise >= M on diag & |.| off-diag
    Ldom = float(np.linalg.eigvalsh(Mdom)[-1])
    ck(Ldom >= Lmax - 1e-9, f"[{w}] lambda_max(Mdom) < lambda_max(M) (M<=Mdom violated)")
    u = 1.0 / np.sqrt(2 * np.arange(N) + 1.0)
    Su = wschur(Mdom, u)
    Sp = perron_schur(Mdom)
    ck(Su >= Ldom - 1e-9, f"[{w}] Schur(psi=u) below lambda_max(Mdom)")
    ck(Sp >= Ldom - 1e-6, f"[{w}] Perron-Schur below lambda_max(Mdom)")
    ldom_ratios.append(Ldom / Lmax)
    su_ratios.append(Su / Lmax)

    # ---- (IV) scalar reductions SUFF & QF  (pure algebra, verified numerically)
    gamma_true = h - l2
    # QF: min over q of (h-q)^2 + gamma(q-m3) is gamma(h-m3) - gamma^2/4 at q*=h-gamma/2
    gam = gamma_true
    qstar = h - gam / 2.0
    lhs_parab = (h - qstar) ** 2 + gam * (qstar - m3)
    rhs_parab = gam * (h - m3) - gam ** 2 / 4.0
    QFalg = abs(lhs_parab - rhs_parab) < 1e-9
    ck(QFalg, f"[{w}] QF vertex algebra mismatch")

    # ---- (V) the ceiling: terminal inequality with ANALYTIC gap lower bound gamma_lo
    hlo = h_lo_trial(x, KX)
    l2hi = lam2_hi(x, N, KX)
    gamma_lo = hlo - l2hi
    ck(gamma_lo <= gamma_true + 1e-9, f"[{w}] gamma_lo not <= gamma_true")
    ck(gamma_lo > 0, f"[{w}] gamma_lo not positive")
    qm3 = q - m3                                    # = x E / 2 (Cor 1.1)
    RHS_an = gamma_lo * qm3
    ck(VarB < RHS_an, f"[{w}] CEILING FAILS: VarB !< gamma_lo*(q-m3)")
    budgets.append(RHS_an / VarB)

    # ---- (VI) cone captures (localise the single residual; kL3 <= true, Thm 1)
    cMq = kL3(v, kernel_array(x, KX, "f-f2"), N, KX)
    cM3 = kL3(v, kernel_array(x, KX, "f-f3"), N, KX)
    ck(cMq <= hq + 1e-9, f"[{w}] kL3[f-f2] not <= h-q")
    ck(cM3 <= (h - m3) + 1e-9, f"[{w}] kL3[f-f3] not <= h-m3")
    rMq = cMq / hq
    rM3 = cM3 / (h - m3)
    cap_min = min(cap_min, rMq, rM3)

    # QF closure with the tight (Perron) lambda_max(M) upper bound and TRUE h-m3:
    QF_close = (Sp ** 2) < (gamma_lo * (h - m3) - gamma_lo ** 2 / 4.0)

    if w in (0.1, 0.3):
        small[w] = (VarB, qm3, gamma_lo)

    print(f" {w:5.2f} {VarB:9.2e} {opm/VarB:7.3f} {BD/VarB:6.3f} {Ldom/Lmax:6.4f} "
          f"{Su/Lmax:6.3f} {abs(VarB-VarM):9.1e} {rMq:5.3f} {rM3:5.3f} {RHS_an/VarB:7.2f} {str(QF_close):>5}")

# ------------------------------------------------------------------ summary claims
print("-" * 78)
ck(max(opm_ratios) < 1.70, "op-norm-minus looseness exceeds 1.70x somewhere")
ck(max(bd_ratios) < 1.47, "Bhatia-Davis looseness exceeds 1.47x somewhere")
ck(max(ldom_ratios) < 1.12, "Perron/Mdom bound looseness exceeds 1.12x somewhere")
ck(min(budgets) > 1.0, "ceiling budget dropped to/below 1.0")
ck(min(budgets) > 3.9, "ceiling budget below 3.9 (expected worst ~4.0 at large wpi)")
ck(0.55 < cap_min < 0.75, "cone capture outside expected [0.55,0.75] band")

# small-wpi asymptotic power counting:  q ~ x^2, m3 ~ x^3  =>  q-m3 ~ x^2;  VarB ~ x^4.
# Hence RHS = gamma*(q-m3) ~ x^3  (gamma ~ x) while VarB ~ x^4, so VarB/RHS ~ x -> 0.
w1, w2 = 0.1, 0.3
V1, R1, _ = small[w1]
V2, R2, _ = small[w2]
x1, x2 = w1 * w1, w2 * w2
# VarB/x^4 and (q-m3)/x^2 should each be roughly constant across the two small points
r_V = (V2 / x2 ** 4) / (V1 / x1 ** 4)
r_R = (R2 / x2 ** 2) / (R1 / x1 ** 2)
print(f"small-wpi scaling: VarB/x^4 ratio(0.3/0.1)={r_V:.3f}, (q-m3)/x^2 ratio={r_R:.3f} (both ~1 confirms orders)")
ck(0.5 < r_V < 2.0, "VarB ~ x^4 scaling not confirmed at small wpi")
ck(0.5 < r_R < 2.0, "(q-m3) ~ x^2 scaling not confirmed at small wpi")
ck(budgets[0] > budgets[1] > budgets[2], "ceiling budget not increasing as wpi->0 (asymptotic margin ~1/x)")

print("=" * 78)
print(f"PASS {_NP} checks" + ("" if _NF == 0 else f"   ({_NF} FAILED)"))
print("=" * 78)
sys.exit(1 if _NF else 0)
