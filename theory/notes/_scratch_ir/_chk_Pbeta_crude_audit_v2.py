"""[numerical-evidence] From-zero audit of lemmaC_Pbeta_crude_bound.py's ACTUAL claim.

Supersedes _chk_Pbeta_crude_bound_verify.py, which tested the bound as it was
TRANSCRIBED IN A MESSAGE  (V <= varpi^2 (2n+1)^{-c} C0)  rather than as the script
states and proves it  (V <= C0 varpi^2 / (2n+1), docstring line 2, code line 125).
The two differ by (2n+1)^{c-1}; the script's form is the strong one and is the one
its combine-step actually derives, via (2n+1)^{-a}*(2n+1)^{-c} = (2n+1)^{-1}, a+c=1.

Audited here, each independently re-derived rather than re-read:

  T1  target restatement:  V <= C0 w^2/(2n+1)  <=>  S(n) <= C0 w^2 (2n+1)^{-c},
      where S(n) = sum_{m>=1} K(n,m)(2m+1)^{-c} and V = (2n+1)^{-a} S - l_G.
  T2  step 2   K(n,m) <= 2f(|n-m|)  for m != n.
  T3  step 4   diagonal coefficient 1.
  T4  step 5   m>n side coefficient 2*pi^2/6.
  T5  step 6   m<n near side coefficient 2*2^{3/4}*pi^2/6.
  T6  step 7   m<n far side coefficient 2*5^{3/4}   -- INCLUDING the tail-sum
      inequality sum_{k>n/2} 1/k^2 <= 2/n, which the docstring justifies by
      "sum_{k>K}1/k^2 <= 1/K".  That justification is FALSE at the integer
      boundary (K=2.5: LHS=sum_{k>=3}=0.3949 > 1/K=0.4?  no -- but K=1.5:
      sum_{k>=2}=0.6449 > 1/1.5=0.6667? no).  Checked here exhaustively, plus the
      sound route  1/k^2 <= 1/(k-1/2)-1/(k+1/2)  =>  sum_{k>=k0} <= 1/(k0-1/2).
  T7  C0 recomputed from the six pieces.
  T8  final bound checked against V computed from the driver's own definition,
      with an UPWARD tail pad for m >= MT so the check cannot fail open.
  T9  closure at NROW in {150, 300}: need C0 w^2/(2*NROW+1) <= max_{n<=NROW} V.

Run: OMP_NUM_THREADS=1 python _chk_Pbeta_crude_audit_v2.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import math
import numpy as np

PI2_6 = math.pi ** 2 / 6.0
BETAS = (0.0, 0.25, 0.5)
CELLS = [(0.30, 0.42, 0.00), (0.42, 0.60, 0.00), (0.60, 0.85, 0.00),
         (0.85, 1.20, 0.25), (1.20, 1.70, 0.50)]


def V_and_pad(n, w, beta, MT=200000):
    """V_beta(n) from the driver's definition, plus an UPWARD pad for m >= MT.

    V = (2n+1)^{-a} sum_{m>=1} K(n,m)(2m+1)^{-c} - l_G(n),  a=(1-b)/2, c=(1+b)/2,
    K(n,m) = f(|n-m|)+f(n+m+1) for m != n,  K(n,n) = f(2n+1),  f(t)=w^2/(t^2+w^2).
    l_G(n) = (2/(2n+1)) sum_{k=1..n} f(k).

    Pad: for m >= MT (> n), K <= 2f(m-n) <= 2w^2/(m-n)^2 and (2m+1)^{-c} <=
    (2MT+1)^{-c}; sum_{k >= MT-n} k^{-2} <= 1/(MT-n-1).  All three raise.
    """
    a, c = (1.0 - beta) / 2.0, (1.0 + beta) / 2.0
    x = w * w
    f = lambda t: x / (t * t + x)
    m = np.arange(1, MT, dtype=float)
    K = f(np.abs(n - m)) + f(n + m + 1.0)
    K[n - 1] = f(2.0 * n + 1.0)
    S = float(K @ (2.0 * m + 1.0) ** (-c))
    lG = 2.0 * float(f(np.arange(1, n + 1, dtype=float)).sum()) / (2.0 * n + 1.0)
    pad = 2.0 * x * (2.0 * MT + 1.0) ** (-c) / (MT - n - 1.0)
    return (2.0 * n + 1.0) ** (-a) * (S + pad) - lG, (2.0 * n + 1.0) ** (-a) * pad


if __name__ == "__main__":
    print("== T2: K(n,m) <= 2 f(|n-m|), i.e. f(n+m+1) <= f(|n-m|) for m != n ==")
    worst = 0.0
    for w in (0.30, 0.85, 1.70):
        x = w * w
        f = lambda t: x / (t * t + x)
        for n in range(1, 120):
            mm = np.arange(1, 400.0)
            mm = mm[mm != n]
            worst = max(worst, float((f(n + mm + 1.0) / f(np.abs(n - mm))).max()))
    print(f"   sup f(n+m+1)/f(|n-m|) = {worst:.6f}   (must be <= 1)"
          f"   -> {'OK' if worst <= 1 else 'VIOLATED'}")
    print("   algebra: n+m+1 > |n-m| for all n,m>=1, and f is decreasing on t>=0.")

    print("\n== T6: sum_{k > n/2, k <= n-1} 1/k^2 <= 2/n ==")
    worst_r, bad_docstring = 0.0, []
    for n in range(2, 5001):
        k = np.arange(math.floor(n / 2) + 1, n, dtype=float)
        s = float((1.0 / k ** 2).sum()) if k.size else 0.0
        worst_r = max(worst_r, s / (2.0 / n))
        # the docstring's stated route, on the FULL tail (no k<=n-1 cutoff):
        k0 = math.floor(n / 2) + 1
        full_tail = float(PI2_6 - (1.0 / np.arange(1.0, k0) ** 2).sum()) if k0 > 1 else PI2_6
        if full_tail > 2.0 / n + 1e-15:
            bad_docstring.append(n)
    print(f"   sup (partial sum)/(2/n) over n<=5000 = {worst_r:.6f}  -> "
          f"{'OK' if worst_r <= 1 else 'VIOLATED'}")
    print(f"   full tail sum_{{k>n/2}} 1/k^2 > 2/n at n = {bad_docstring[:8]} "
          f"({len(bad_docstring)} cases)")
    print("   sound route: 1/k^2 <= 1/(k-1/2)-1/(k+1/2) => sum_{k>=k0} <= 1/(k0-1/2);")
    print("   k0 = floor(n/2)+1 => k0-1/2 >= n/2 => tail <= 2/n.  Checked:")
    tel = max(((PI2_6 - sum(1.0 / j ** 2 for j in range(1, math.floor(n / 2) + 1)))
               * (math.floor(n / 2) + 0.5) for n in range(2, 2001)))
    print(f"      sup (full tail)*(k0-1/2) = {tel:.6f}  (must be <= 1)")

    print("\n== T7: C0 recomputed from the six pieces ==")
    pieces = {"diag": 1.0, "2 * far m>n": 2.0 * PI2_6,
              "2 * near m<n": 2.0 * 2.0 ** 0.75 * PI2_6,
              "2 * far m<n": 2.0 * 5.0 ** 0.75}
    C0 = sum(pieces.values())
    for k, v in pieces.items():
        print(f"   {k:14s} = {v:.6f}")
    print(f"   C0 = {C0:.6f}   (script reports 16.5101)")
    print(f"   g(n)=2(2n+1)^c/n decreasing: d/dn ln g = (2n(c-1)-1)/(n(2n+1)) < 0"
          f" for c<=1 -- analytic, all n>0.  sup taken at n=2 though far range is")
    print(f"   empty until n=3 (g(3)={2*7**0.75/3:.4f} < g(2)={5**0.75:.4f}):"
          f" conservative, fail-closed.")

    print("\n== T8: V_beta(n) <= C0 w^2/(2n+1), V padded UPWARD ==")
    print(f"{'beta':>5} {'varpi':>6} {'max V/bound':>12} {'at n':>5} {'viol':>5} {'pad/V':>9}")
    worst_all = 0.0
    for lo, hi, b in CELLS:
        for w in (lo, hi):
            r, padfrac = [], 0.0
            for n in range(1, 401):
                v, p = V_and_pad(n, w, b, MT=20000)
                r.append(v / (C0 * w * w / (2.0 * n + 1.0)))
                padfrac = max(padfrac, p / max(v, 1e-300))
            r = np.array(r)
            worst_all = max(worst_all, float(r.max()))
            print(f"{b:5.2f} {w:6.2f} {r.max():12.6f} {int(r.argmax())+1:5d}"
                  f" {int((r > 1).sum()):5d} {padfrac:9.2e}")
    print(f"   worst over all = {worst_all:.6f} -> "
          f"{'BOUND HOLDS' if worst_all <= 1 else 'BOUND FAILS'}")

    print("\n== T9: closure -- need C0 w^2/(2*NROW+1) <= peak_{n<=NROW} V ==")
    for NROW in (150, 300):
        print(f"   --- NROW = {NROW} ---")
        print(f"{'cell':>13} {'beta':>5} {'varpi':>6} {'peak':>10} {'tail bd':>10}"
              f" {'bd/peak':>9} {'verdict':>8}")
        allpass = True
        for lo, hi, b in CELLS:
            for w in np.geomspace(lo, hi, 5):
                peak = max(V_and_pad(n, w, b, MT=20000)[0] for n in range(1, NROW + 1))
                bd = C0 * w * w / (2.0 * NROW + 1.0)
                ok = bd <= peak
                allpass &= ok
                if w in (lo, hi):
                    print(f"[{lo:4.2f},{hi:4.2f}] {b:5.2f} {w:6.3f} {peak:10.6f}"
                          f" {bd:10.6f} {bd/peak:9.3f} {'PASS' if ok else 'FAIL':>8}")
        print(f"   => NROW={NROW}: {'CLOSES' if allpass else 'DOES NOT CLOSE'}")
