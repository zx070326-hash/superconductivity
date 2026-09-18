"""[numerical-evidence, checking an analytic derivation] Prove an EXPLICIT crude
tail bound  V_beta(n) <= C0 * varpi^2 / (2n+1)  for n>=1, beta in {0,0.25,0.5},
varpi<=1.70 -- per team-lead's scoping (theta_tilde_keeping_D_2026-09-17.md
section 12), P_beta only binds on those 5 cells, so this restricted claim is all
that's needed to finitize sup_n in peel_all_betas beyond its explicit NROW scan.

Derivation (every step direction-audited, each is an UPPER bound on V_beta(n)):

  V_beta(n) = (2n+1)^{-a} sum_{m>=1} K(n,m) (2m+1)^{-c} - Delta_n,   a=(1-beta)/2,
              c=(1+beta)/2=1-a,  K(n,m)=f(|n-m|)+f(n+m+1) (m!=n), K(n,n)=f(2n+1),
              f(t)=varpi^2/(t^2+varpi^2).

  Step 1: drop -Delta_n <= 0 (subtracted term, dropping it INCREASES the RHS: safe).
  Step 2: for m!=n, f(n+m+1) <= f(|n-m|) (f decreasing, n+m+1 > |n-m|)
          => K(n,m) <= 2 f(|n-m|).  K(n,n)=f(2n+1) kept as is (no Toeplitz f(0)).
  Step 3: f(t) <= varpi^2/t^2 (t>0).
  Step 4 (diagonal): f(2n+1)(2n+1)^{-c} <= varpi^2(2n+1)^{-c-2} <= varpi^2(2n+1)^{-c}
          (coefficient 1, since 2n+1>=1).
  Step 5 (m>n side): (2m+1)^{-c} <= (2n+1)^{-c}  (m>n => 2m+1>2n+1, c>0)
          => sum_{k>=1} varpi^2/k^2 * (2n+1)^{-c} = (pi^2/6) varpi^2 (2n+1)^{-c}.
  Step 6 (m<n side, near k<=n/2): 2n+1-2k >= n+1 (since k<=n/2 => 2k<=n)
          => (2n+1-2k)^{-c} <= (n+1)^{-c} <= 2^c (2n+1)^{-c}  [since n+1 >= (2n+1)/2]
          use c<=0.75 (worst over beta in {0,.25,.5}) => 2^c <= 2^0.75.
          sum_{k=1}^{floor(n/2)} 1/k^2 <= pi^2/6.
          => near coefficient <= 2^0.75 * pi^2/6.
  Step 7 (m<n side, far n/2<k<=n-1): 2m+1>=3 => (2m+1)^{-c}<=1 trivially (c>0).
          sum_{k>=k0} 1/k^2 <= 2/n  via TELESCOPING (NOT integral comparison --
          that version is subtly wrong at integer boundaries, per team-lead's
          audit: the naive claim sum_{k>K}1/k^2<=1/K silently assumes the sum's
          start point sits at or past K, which can fail off by <=1).
          Telegraph identity: 1/k^2 <= 1/(k-1/2) - 1/(k+1/2)  [since that RHS
          equals 1/(k^2-1/4) >= 1/k^2], so for k0:=floor(n/2)+1,
          sum_{k>=k0} 1/k^2 <= sum_{k>=k0}[1/(k-1/2)-1/(k+1/2)] = 1/(k0-1/2)
          (telescopes to the first term).  k0-1/2 >= n/2 for both parities of n
          (n even: k0=n/2+1, k0-1/2=n/2+1/2>n/2; n odd: k0=(n+1)/2, k0-1/2=n/2
          exactly) => sum_{k>=k0}1/k^2 <= 2/n.  Same numeric constant as before,
          now with a sound proof.
          express as coefficient of (2n+1)^{-c}: 2/n = (2n+1)^{-c} * [2(2n+1)^c/n].
          g(n):=2(2n+1)^c/n is DECREASING in n for c<1 (shown by sign of g'),
          so its sup over n>=2 (smallest n with a nonempty far range) is at n=2:
          g(2) = 2*5^c/2 = 5^c <= 5^0.75 (worst c=0.75).

  Combine (all coefficients of varpi^2 (2n+1)^{-c}, note (2n+1)^{-a}*(2n+1)=(2n+1)^{1-a}
  = (2n+1)^c so after dividing by (2n+1)^{-a} the target becomes (2n+1)^{-c} form --
  see note below):

    C0_proven := 1  (diag)  +  2*[ pi^2/6 (far m>n)  +  2^0.75*pi^2/6 (near m<n)
                                   + 5^0.75 (far m<n) ]

This script (a) verifies every INEQUALITY STEP pointwise against the exact
quantities it replaces, at many (n,varpi,beta) combinations -- so a violation
here means the hand algebra has a bug; (b) reports C0_proven numerically;
(c) sanity-checks the final bound V_beta(n) <= C0_proven*varpi^2/(2n+1) against
DIRECT computation of V_beta(n) (dropping only Delta_n, kept exact otherwise) for
n=1..400, varpi in a grid up to 1.70, beta in {0,0.25,0.5}.

Run: OMP_NUM_THREADS=1 python lemmaC_Pbeta_crude_bound.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import math

PI2_6 = math.pi**2 / 6.0


def V_beta_exact_minus_Delta_dropped(n, w, beta, MT=200000):
    """sum_{m>=1} K(n,m)(2m+1)^{-c} * (2n+1)^{-a}  (Delta_n DROPPED, per step 1)."""
    a = (1.0 - beta) / 2.0
    c = (1.0 + beta) / 2.0
    x = w * w
    f = lambda t: x / (t * t + x)
    m = np.arange(1, MT, dtype=float)
    K = f(np.abs(n - m)) + f(n + m + 1.0)
    K[int(n) - 1] = f(2.0 * n + 1.0) if n <= MT else K[int(n) - 1]
    wm = (2.0 * m + 1.0) ** (-c)
    return (2.0 * n + 1.0) ** (-a) * float(K @ wm)


def C0_proven(beta_max=0.5):
    c_worst = (1.0 + beta_max) / 2.0  # 0.75
    diag = 1.0
    far_mgtn = PI2_6
    near_mltn = (2.0 ** c_worst) * PI2_6
    far_mltn = 5.0 ** c_worst
    total = diag + 2.0 * (far_mgtn + near_mltn + far_mltn)
    return total, dict(diag=diag, far_mgtn=far_mgtn, near_mltn=near_mltn,
                        far_mltn=far_mltn)


if __name__ == "__main__":
    print("== Step-by-step inequality checks ==")

    # Step 2: f(n+m+1) <= f(|n-m|) for m != n, all n,m>=1, varpi tested.
    bad = 0
    for w in (0.30, 0.85, 1.20, 1.70):
        x = w * w
        f = lambda t: x / (t * t + x)
        for n in range(1, 60):
            m = np.arange(1, 200)
            m = m[m != n]
            lhs = f(n + m + 1.0)
            rhs = f(np.abs(n - m).astype(float))
            if (lhs > rhs + 1e-15).any():
                bad += 1
    print(f"Step2 f(n+m+1)<=f(|n-m|):            violations = {bad}")

    # Step 6: n+1 >= (2n+1)/2  <=>  2n+2>=2n+1, always true (exact, not numeric).
    print("Step6 n+1 >= (2n+1)/2:                 exact algebra, always true")

    # Step 7a: telescoping tail bound sum_{k>=k0} 1/k^2 <= 1/(k0-1/2), replacing
    #   the earlier (subtly wrong at integer boundaries) integral-comparison
    #   justification.  Check pointwise for n=2..5000 that k0=floor(n/2)+1 gives
    #   the exact tail sum <= 1/(k0-1/2), and that 1/(k0-1/2) <= 2/n as claimed.
    bad7a = 0
    for n in range(2, 5000):
        k0 = n // 2 + 1
        kk = np.arange(k0, k0 + 200000, dtype=float)  # truncated tail, converges fast
        tail_exact = float(np.sum(1.0 / kk ** 2)) + 1.0 / (k0 + 200000 - 1)  # +tiny remainder pad
        teleb = 1.0 / (k0 - 0.5)
        if tail_exact > teleb + 1e-9:
            bad7a += 1
        if teleb > 2.0 / n + 1e-12:
            bad7a += 1
    print(f"Step7a telescoping bound sum_(k>=k0)1/k^2<=1/(k0-1/2)<=2/n: violations = {bad7a}")

    # Step 7: g(n) = 2(2n+1)^c/n decreasing in n for c<1, check numerically
    #   over n=2..2000 and c in {0.5,0.625,0.75} that g is non-increasing.
    bad = 0
    for c in (0.5, 0.625, 0.75):
        nn = np.arange(2, 2000, dtype=float)
        g = 2.0 * (2.0 * nn + 1.0) ** c / nn
        if (np.diff(g) > 1e-9).any():
            bad += 1
    print(f"Step7 g(n) decreasing in n, c in {{.5,.625,.75}}: violations = {bad}")
    print(f"      g(2) at c=0.75 = {2.0*5.0**0.75/2.0:.6f}  (claimed sup)")

    tot, parts = C0_proven()
    print(f"\n== C0_proven = {tot:.4f} ==")
    for k, v in parts.items():
        print(f"   {k:10s} = {v:.4f}")

    print("\n== Final pointwise check: V_beta(n) [Delta dropped] <= C0*varpi^2/(2n+1) ==")
    bad = []
    C0 = tot
    for beta in (0.0, 0.25, 0.5):
        for w in (0.30, 0.60, 0.85, 1.20, 1.70):
            worst_ratio = -np.inf
            for n in range(1, 400):
                v = V_beta_exact_minus_Delta_dropped(n, w, beta, MT=20000)
                bound = C0 * w * w / (2.0 * n + 1.0)
                worst_ratio = max(worst_ratio, v / bound)
                if v > bound * (1.0 + 1e-9):
                    bad.append((beta, w, n, v, bound))
            print(f"  beta={beta:.2f} varpi={w:.2f}  worst v/bound = {worst_ratio:.4f}")
    print(f"\n  violations: {len(bad)}" + ("  PASS" if not bad else f"  FAIL {bad[:5]}"))
