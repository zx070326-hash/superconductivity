#!/usr/bin/env python3
r"""
Root-cause probe: Is Tc universality a consequence of Fermi statistics?

Numerically verifies each link of the chain:
  Fermi Matsubara equidistance → harmonic row sum → h ~ ln w
    → H(x) nearly linear → small Jensen gap → spectral-shape insensitivity

And constructs the bosonic counterpart to show where the chain breaks.

All numbers here are float64 / numerical-evidence, NOT certificates.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'theory', 'derivations'))
import numpy as np
from scipy import linalg

# ============================================================
# 0.  Fermi H-matrix builder  (from p4_interchange_uniformity)
# ============================================================
def brk(k, w):
    """[[k]](w) = w^2 / (k^2 + w^2), with [[0]]=1."""
    if k == 0:
        return 1.0
    return w*w / (k*k + w*w)

def H_fermi(N, w):
    """KAY II  H^(N)(w) — Fermi Matsubara, frequencies (2n+1)pi T."""
    M = np.zeros((N, N))
    for n in range(N):
        for m in range(N):
            nrm = np.sqrt((2*n+1)*(2*m+1))
            val = brk(n+m+1, w) / nrm
            if n != m:
                val += brk(abs(n-m), w) / nrm
            else:
                s = sum(2*brk(k, w) for k in range(1, n+1))
                val -= s / (2*n+1)
            M[n, m] = val
    return M

# ============================================================
# 1.  Bosonic analogue: replace (2n+1) → (2n) everywhere
#     Frequencies omega_m = 2m pi T   (m = 0,1,2,...)
#     The m=0 mode has zero frequency → potential divergence.
#     We handle this in two ways:
#       (a) Skip m=0: sum over m=1,2,... (infrared cutoff)
#       (b) Include m=0 with regularization
# ============================================================
def brk_B(k, w):
    """Bosonic kernel: w^2 / ((2k)^2 + w^2) for k>=1."""
    if k == 0:
        return 1.0  # placeholder; zero-mode treated separately
    return w*w / (4*k*k + w*w)

def H_boson_skip0(N, w):
    """Bosonic analogue of H, skipping m=0. Indices n=1,2,...,N.
    Weight: 1/(2n) instead of 1/(2n+1).
    Kernel: w^2/((2k)^2 + w^2) instead of w^2/(k^2+w^2).
    """
    # Use n=1,...,N to avoid zero-frequency mode
    M = np.zeros((N, N))
    for i in range(N):
        n = i + 1
        for j in range(N):
            m = j + 1
            nrm = np.sqrt(4*n*m)  # = 2*sqrt(n*m), weight factor
            # H3 analogue: kernel at n+m
            val = brk_B(n+m, w) / nrm
            # H2 analogue: kernel at |n-m| (nonzero for n!=m)
            if n != m:
                val += brk_B(abs(n-m), w) / nrm
            else:
                # H1 analogue: diagonal subtraction
                s = sum(2*brk_B(k, w) for k in range(1, n))
                val -= s / (2*n)
            M[i, j] = val
    return M

# ============================================================
# 2.  Row sums  s_n  and their asymptotics
# ============================================================
def row_sum_fermi(n, w, N=5000):
    """Sum_m |H_{nm}| for fermi H, m=0..N-1."""
    s = 0.0
    for m in range(N):
        nrm = np.sqrt((2*n+1)*(2*m+1))
        val = brk(n+m+1, w) / nrm
        if n != m:
            val += brk(abs(n-m), w) / nrm
        else:
            sub = sum(2*brk(k, w) for k in range(1, n+1))
            val -= sub / (2*n+1)
        s += abs(val)
    return s

def S1_fermi(w, N=5000):
    """S1 = sum_{k>=1} [[k]](w) = sum w^2/(k^2+w^2)."""
    return sum(w*w/(k*k + w*w) for k in range(1, N+1))

# ============================================================
# 3.  Eigenvalue computations
# ============================================================
def h_fermi(w, N=None):
    """Top eigenvalue of H_fermi^(N)(w)."""
    if N is None:
        N = max(200, int(20*w))
    H = H_fermi(N, w)
    return float(np.linalg.eigvalsh(H)[-1])

def h_boson(w, N=None):
    """Top eigenvalue of H_boson^(N)(w)."""
    if N is None:
        N = max(200, int(20*w))
    H = H_boson_skip0(N, w)
    return float(np.linalg.eigvalsh(H)[-1])

# ============================================================
# MAIN ANALYSIS
# ============================================================
if __name__ == '__main__':
    print("=" * 72)
    print("ROOT-CAUSE PROBE: Fermi statistics and Tc universality")
    print("=" * 72)

    # ----------------------------------------------------------
    # STEP 1→2: Matsubara equidistance → row-sum structure
    # ----------------------------------------------------------
    print("\n--- Step 1→2: Row-sum asymptotics (Fermi vs Boson) ---")
    print(f"{'w':>8} {'S1_fermi':>12} {'ln(w)+c0':>12} {'pi*w/2':>12}")
    c0 = np.euler_gamma + np.log(4) - 1.5  # 0.46351
    for w in [1, 2, 5, 10, 20, 50, 100]:
        s1 = S1_fermi(w, 10000)
        asymp = np.log(w) + c0 + 1  # 2B = ln w + gamma + ln4 - 1/2
        piw2 = np.pi * w / 2
        print(f"{w:8d} {s1:12.6f} {asymp:12.6f} {piw2:12.6f}")

    # Bosonic: S1_B = sum w^2/((2k)^2 + w^2) = sum w^2/(4k^2+w^2)
    print(f"\n{'w':>8} {'S1_boson':>12} {'(pi*w/4)coth(pi*w/2)':>18}")
    for w in [1, 2, 5, 10, 20, 50, 100]:
        s1b = sum(w*w/(4*k*k + w*w) for k in range(1, 10001))
        # Exact: (pi*w/4)*coth(pi*w/2) - 1/2
        exact = (np.pi*w/4)/np.tanh(np.pi*w/2) - 0.5
        print(f"{w:8d} {s1b:12.6f} {exact:12.6f}")

    # Key difference: Fermi S1 ~ (pi*w/2)*coth(pi*w) ~ pi*w/2 (large w)
    # while Boson S1 ~ (pi*w/4)*coth(pi*w/2) ~ pi*w/4 (large w)
    # Both are linear in w for large w.
    # But the DIAGONAL subtraction changes things for the eigenvalue.

    # ----------------------------------------------------------
    # STEP 2→3: h(w) asymptotics — Fermi vs Boson
    # ----------------------------------------------------------
    print("\n--- Step 2→3: Eigenvalue h(w) asymptotics ---")
    print(f"\n{'w':>8} {'h_fermi':>12} {'ln(w)+c0':>12} {'h_fermi-ln(w)':>14}")
    for w in [3, 5, 10, 20, 50, 100, 200]:
        N = max(400, int(30*w))
        hf = h_fermi(w, N)
        lnw = np.log(w) + c0
        print(f"{w:8d} {hf:12.6f} {lnw:12.6f} {hf - np.log(w):14.6f}")

    print(f"\n{'w':>8} {'h_boson':>12} {'h_boson/w':>12} {'h_boson/w^2':>12}")
    for w in [3, 5, 10, 20, 50]:
        N = max(400, int(30*w))
        hb = h_boson(w, N)
        print(f"{w:8d} {hb:12.6f} {hb/w:12.6f} {hb/w**2:12.8f}")

    # ----------------------------------------------------------
    # STEP 3→4: H(x) = h(sqrt(x)) and its curvature
    # ----------------------------------------------------------
    print("\n--- Step 3→4: H(x) curvature (Jensen gap) ---")
    print("Fermi: x^2 * H''(x) should -> -1/2 (nearly linear in log x)")
    dx = 0.01
    print(f"{'w':>8} {'H(x)':>12} {'x^2*H\"':>12} {'H\"/H':>12}")
    for w in [2, 5, 10, 20, 50, 100]:
        N = max(600, int(40*w))
        x = w*w
        h0 = h_fermi(w, N)
        hp = h_fermi(w*(1+dx), N)
        hm = h_fermi(w*(1-dx), N)
        # H(x) = h(sqrt(x)), H'(x) = h'(w)/(2w), H''(x) ≈ ...
        # Use finite difference on h(w) to get H''(x):
        # d^2h/dw^2 from three points
        h2 = (hp - 2*h0 + hm) / (w*dx)**2
        # H''(x) = (h''(w) - h'(w)/w) / (4w^2)
        h1 = (hp - hm) / (2*w*dx)
        Hpp = (h2 - h1/w) / (4*w*w)
        print(f"{w:8d} {h0:12.6f} {x*x*Hpp:12.6f} {Hpp/h0*x*x:12.6f}")

    # ----------------------------------------------------------
    # STEP 4→5: Jensen gap for model spectral distributions
    # ----------------------------------------------------------
    print("\n--- Step 4→5: Jensen gap for two-atom model spectra ---")
    print("If P = (1-a)*delta(w1) + a*delta(w2), Jensen gap = ")
    print("h(w_rms) - integral h(w)P(dw)")
    print("For Fermi (H concave), this gap > 0 -> Tc insensitive to shape")
    print()

    w1_phys = 5.0  # omega_1 / 2pi T
    N = max(800, int(40*max(20, w1_phys)))
    alphas = [0.1, 0.3, 0.5]
    ratios = [1.5, 2.0, 3.0, 5.0]
    print(f"{'r=w2/w1':>8} {'alpha':>6} {'h(w_rms)':>10} {'<h>':>10} {'gap':>10} {'gap/h':>10}")
    for r in ratios:
        for a in alphas:
            w2 = w1_phys * r
            wrms = np.sqrt((1-a)*w1_phys**2 + a*w2**2)
            Nuse = max(N, int(40*w2))
            h_rms = h_fermi(wrms, Nuse)
            h_avg = (1-a)*h_fermi(w1_phys, Nuse) + a*h_fermi(w2, Nuse)
            gap = h_rms - h_avg
            print(f"{r:8.1f} {a:6.1f} {h_rms:10.5f} {h_avg:10.5f} {gap:10.5f} {gap/h_rms:10.5f}")

    # ----------------------------------------------------------
    # Bosonic Jensen gap (for comparison)
    # ----------------------------------------------------------
    print("\n--- Bosonic Jensen gap ---")
    print(f"{'r=w2/w1':>8} {'alpha':>6} {'h_B(w_rms)':>10} {'<h_B>':>10} {'gap':>10} {'gap/h':>10}")
    for r in [1.5, 2.0, 3.0]:
        for a in [0.3, 0.5]:
            w2 = 5.0 * r
            wrms = np.sqrt((1-a)*25.0 + a*w2**2)
            Nuse = max(400, int(30*w2))
            hb_rms = h_boson(wrms, Nuse)
            hb_avg = (1-a)*h_boson(5.0, Nuse) + a*h_boson(w2, Nuse)
            gapb = hb_rms - hb_avg
            print(f"{r:8.1f} {a:6.1f} {hb_rms:10.5f} {hb_avg:10.5f} {gapb:10.5f} {gapb/hb_rms:10.5f}")

    # ----------------------------------------------------------
    # Hypothetical: tau_n with different decay rates
    # ----------------------------------------------------------
    print("\n--- Hypothetical: effect of tau_n decay rate on eigenvalue ---")
    print("Fermi: tau_n ~ 1/sqrt(2n+1) * 1/(1+(n/w)^2) → log-uniform mass")
    print("If tau decays as n^{-s}:")
    print("  s=1/2 (Fermi-like): sum ~ ln N → h ~ ln w")
    print("  s=1   (faster):     sum ~ const → h bounded")
    print("  s=1/4 (slower):     sum ~ N^{1/2} → h ~ w^{1/2} ?")

    # Construct model operator O[f] with weight w_n = n^{-2s} instead of 1/(2n+1)
    # to test how decay rate affects eigenvalue scaling
    def model_operator(N, w, s_decay):
        """Model rank-1-plus-residual with weight n^{-2s} and Lorentzian kernel."""
        M = np.zeros((N, N))
        for n in range(N):
            nn = max(n, 1)  # avoid zero
            for m in range(N):
                mm = max(m, 1)
                wt = 1.0 / (nn**s_decay * mm**s_decay)
                k_nm = w*w / ((n+m+1)**2 + w*w)
                val = k_nm * wt
                M[n, m] = val
        return M

    print(f"\n{'w':>8} {'s=0.5(Fermi)':>14} {'s=1.0':>14} {'s=0.25':>14}")
    for w in [5, 10, 20, 50, 100]:
        N = max(200, int(10*w))
        vals = []
        for s in [0.5, 1.0, 0.25]:
            M = model_operator(N, w, s)
            lam1 = float(np.linalg.eigvalsh(M)[-1])
            vals.append(lam1)
        print(f"{w:8d} {vals[0]:14.5f} {vals[1]:14.5f} {vals[2]:14.5f}")

    # Check scaling
    print("\nScaling analysis: h(200)/h(100)")
    for s in [0.5, 1.0, 0.25]:
        h100 = float(np.linalg.eigvalsh(model_operator(2000, 100, s))[-1])
        h200 = float(np.linalg.eigvalsh(model_operator(4000, 200, s))[-1])
        print(f"  s={s:.2f}: h(200)/h(100) = {h200/h100:.4f}  (ln ratio = {np.log(h200/h100)/np.log(2):.4f})")

    # ----------------------------------------------------------
    # Classical limit: discrete → continuous
    # ----------------------------------------------------------
    print("\n--- Classical limit (T → ∞, i.e. w → 0 for fixed omega) ---")
    print("In the classical limit, Matsubara sum → integral:")
    print("  sum_n f(n) → integral f(n) dn  (spacing = 1)")
    print("  h(w)/w^2 → g(2) as w → 0  (known)")
    print()
    print(f"{'w':>8} {'h(w)':>12} {'h/w^2':>12} {'g(2)=1.318':>12}")
    for w in [0.01, 0.03, 0.1, 0.3, 0.5, 1.0, 2.0]:
        hf = h_fermi(w, 400)
        print(f"{w:8.3f} {hf:12.8f} {hf/w**2:12.6f} {'1.31814':>12}")

    print("\n--- Crossover scale ---")
    print("The crossover from quantum (discrete) to classical (continuous)")
    print("occurs when w ~ 1, i.e., omega_E ~ 2*pi*T")
    print("At T = Tc: w_c = omega_E/(2*pi*Tc)")
    print("  Weak coupling: w_c = omega_E/(2*pi*0.253*omega_E*exp(-1/lambda))")
    print("                    = exp(1/lambda)/(2*pi*0.253) >> 1")
    print("  → always deep quantum regime at Tc")
    print("  Strong coupling (lambda=3): w_c ~ 2-3")
    print("  → quantum discreteness still dominant")

    # ----------------------------------------------------------
    # Summary: Goldilocks argument
    # ----------------------------------------------------------
    print("\n" + "=" * 72)
    print("GOLDILOCKS SUMMARY")
    print("=" * 72)
    print(
"The universality chain:\n"
"\n"
"  FERMI STATISTICS\n"
"    |  (antiperiodic BC -> odd Matsubara: omega_n = (2n+1)*pi*T)\n"
"  EQUIDISTANT MATSUBARA FREQUENCIES\n"
"    |  (row sum ~ Sum 1/(n^2+w^2) with weight 1/(2m+1))\n"
"  HARMONIC SERIES with Lorentz envelope\n"
"    |  (Sum_{n=0}^{~w} 1/(2n+1) ~ (1/2) ln w)\n"
"  h(w) ~ ln w\n"
"    |  (H(x) = h(sqrt(x)) = (1/2) ln x + c0)\n"
"  H(x) NEARLY LINEAR in ln x\n"
"    |  (H''(x) = -1/(2x^2) -> 0 as x -> inf)\n"
"  SMALL JENSEN GAP\n"
"    |  (gap/h ~ Var(x)/(2<x>^2 h) -> 0 as w -> inf)\n"
"  Tc INSENSITIVE TO SPECTRAL SHAPE\n"
"\n"
"KEY LINK: tau_n = (2n+1)^{-1/2} * g(n/w)\n"
"  The 1/(2n+1) weight is FERMI. It gives Sum tau_n^2 ~ (1/2) ln w.\n"
"\n"
"BOSONIC COMPARISON:\n"
"  Replace (2n+1) -> (2n): weight 1/(2n) and kernel w^2/((2n)^2+w^2)\n"
"  -> Sum tau_n^2 ~ (1/4) ln w  (same log growth!)\n"
"  -> h_B(w) also grows logarithmically\n"
"  -> But the ZERO MODE (n=0) diverges: needs IR cutoff\n"
"  -> With cutoff: h_B ~ ln w with different c0\n"
"\n"
"CONCLUSION: The LOG GROWTH is generic to equidistant Matsubara frequencies\n"
"(both fermi and boson). What Fermi statistics SPECIFICALLY provides:\n"
"  1. No zero-mode divergence (omega_0 = pi*T > 0)\n"
"  2. The specific value of c0 (hence Tc prefactor)\n"
"  3. The specific H''(x) profile (hence quantitative Jensen bound)\n"
"\n"
"The truly CRITICAL feature is EQUIDISTANT MATSUBARA FREQUENCIES,\n"
"which is a consequence of PERIODICITY in imaginary time --\n"
"shared by BOTH Fermi and Bose statistics.\n"
    )

    print("DONE.")
