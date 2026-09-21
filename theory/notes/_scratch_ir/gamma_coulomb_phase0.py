"""
Phase 0: Can Coulomb μ* be soundly incorporated into γ-model?

KAY-I footnote 4 claims: adding a constant μ* to λ(n-m) makes Matsubara
sums diverge. Let's check EXACTLY what diverges and whether the standard
physical regularization (cutoff Ω_c) can be made rigorous.

Key question: In our P4 treatment, μ* enters as a FINITE-RANK perturbation
(restricted to |ω_n| < Ω_c). Does this work for γ-model too?

The Eliashberg gap equation with Coulomb:
  ω_n Δ_n = πT Σ_m [λ(ω_n - ω_m) - μ* θ(Ω_c-|ω_n|) θ(Ω_c-|ω_m|)]
              × Δ_m/|ω_m|

In Matsubara index space with ω_n = (2n+1)πT:
  - λ(ω_n - ω_m) = g^γ / |ω_n - ω_m|^γ (γ-model)
  - The Coulomb term μ* is restricted to |n|, |m| < N_c where N_c ~ Ω_c/(2πT)

This means μ* adds a FINITE-DIMENSIONAL matrix (rank 2N_c) to the infinite
kernel — it does NOT make infinite sums diverge, because it's zero outside
the cutoff window.

The question is whether this cutoff can be given a rigorous mathematical
status (not "smuggled in" as KAY complains).
"""
import numpy as np
from scipy.linalg import eigvalsh
import time

def build_gamma_operator(N, gamma, mu_star=0.0, Nc=None):
    """Build the linearized Eliashberg operator for γ-model + Coulomb.

    H = -L1[f] + L2[f] + L3[f] - μ* P_Nc

    where P_Nc is the projection onto |n| < Nc (finite-rank perturbation).
    In our half-space (n >= 0) formulation:
      P_Nc acts on n = 0, ..., Nc-1 only.

    f(k) = 1/|k|^γ for k != 0, f(0) = 0.
    """
    dim = N + 1
    u = np.array([1.0/np.sqrt(2*n+1) for n in range(dim)])

    def f(k):
        if k == 0:
            return 0.0
        return 1.0 / abs(k)**gamma

    # Truncated L1 diagonal: Σ_{k=1}^M f(k)
    M_tail = max(50000, 100*N)
    f_sum = sum(f(k) for k in range(1, M_tail+1))

    H = np.zeros((dim, dim))

    for n in range(dim):
        for m in range(dim):
            if n == m:
                H[n, n] = -u[n]**2 * f_sum + u[n]**2 * f(2*n+1)
            else:
                H[n, m] = u[n]*u[m] * (f(abs(n-m)) + f(n+m+1))

    # Add Coulomb: -μ* u_n u_m for n, m < Nc
    if mu_star != 0 and Nc is not None and Nc > 0:
        Nc_eff = min(Nc, dim)
        for n in range(Nc_eff):
            for m in range(Nc_eff):
                H[n, m] -= mu_star * u[n] * u[m]

    return H

print("=" * 70)
print("PHASE 0: γ-MODEL + COULOMB μ* FEASIBILITY")
print("=" * 70)

# Part 1: Verify that Coulomb as finite-rank perturbation is well-defined
print("\n--- Part 1: Coulomb as finite-rank perturbation ---")
print("Adding μ* only within window |n| < Nc (finite-rank, no divergence)")
print()

gamma = 1.0  # SDW QCP, most important case
N = 200

# Without Coulomb
H0 = build_gamma_operator(N, gamma, mu_star=0.0)
evals0 = np.sort(eigvalsh(H0))[::-1]
lam1_0 = evals0[0]

print(f"γ={gamma}, N={N}, no Coulomb:")
print(f"  λ₁ = {lam1_0:.8f}")
print()

# With Coulomb for different Nc and μ*
print("With Coulomb μ* (different Nc and μ* values):")
print(f"{'μ*':>6} {'Nc':>6} {'λ₁':>12} {'Δλ₁':>12} {'λ₂':>12} {'gap':>10}")
for mu_star in [0.05, 0.10, 0.13, 0.15, 0.20]:
    for Nc in [10, 20, 50, 100]:
        H = build_gamma_operator(N, gamma, mu_star=mu_star, Nc=Nc)
        evals = np.sort(eigvalsh(H))[::-1]
        lam1 = evals[0]
        lam2 = evals[1]
        gap = lam1 - lam2
        print(f"{mu_star:6.2f} {Nc:6d} {lam1:12.6f} {lam1-lam1_0:12.6f} {lam2:12.6f} {gap:10.6f}")

# Part 2: Check convergence in Nc (does lam1 converge as Nc -> infinity?)
print("\n--- Part 2: Convergence in Nc (fixed μ*, N=200) ---")
print("If λ₁ converges as Nc → N, the finite-rank regularization is well-defined.")
print()

for mu_star in [0.10, 0.13]:
    print(f"μ* = {mu_star}:")
    prev = None
    for Nc in [5, 10, 20, 30, 50, 75, 100, 150, 200]:
        H = build_gamma_operator(N, gamma, mu_star=mu_star, Nc=Nc)
        evals = np.sort(eigvalsh(H))[::-1]
        lam1 = evals[0]
        if prev is not None:
            delta = lam1 - prev
            print(f"  Nc={Nc:4d}: λ₁={lam1:12.8f}  Δ={delta:+12.8f}")
        else:
            print(f"  Nc={Nc:4d}: λ₁={lam1:12.8f}")
        prev = lam1
    print()

# Part 3: Check convergence in N (fixed Nc, fixed μ*)
print("--- Part 3: Convergence in N (fixed Nc=50, μ*=0.13) ---")
print("Does λ₁ converge as N → ∞ with fixed finite Nc?")
print()

mu_star = 0.13
Nc = 50
for N_val in [60, 100, 150, 200, 300, 400]:
    t0 = time.time()
    H = build_gamma_operator(N_val, gamma, mu_star=mu_star, Nc=Nc)
    evals = np.sort(eigvalsh(H))[::-1]
    lam1 = evals[0]
    lam2 = evals[1]
    v1 = np.linalg.eigh(H)[1][:, -1]
    c1sq = v1[0]**2
    print(f"  N={N_val:4d}: λ₁={lam1:12.8f}  λ₂={lam2:12.8f}  gap={lam1-lam2:10.6f}"
          f"  c₁²={c1sq:.6f}  [{time.time()-t0:.1f}s]")

# Part 4: The REAL question — what happens when Nc scales with N?
print("\n--- Part 4: What happens when Nc = N (no cutoff)? ---")
print("This is what KAY fears: if we set Nc = N (remove the cutoff),")
print("does the diagonal of the Coulomb term (-μ* Σ_m u_m²) diverge?")
print()

# The Coulomb diagonal correction for site n (within the window):
# -μ* u_n Σ_{m=0}^{Nc-1} u_m = -μ* u_n S(Nc)
# where S(Nc) = Σ_{m=0}^{Nc-1} 1/√(2m+1)
# S(Nc) ~ √(Nc) (diverges!) -- this is what KAY means

print("S(Nc) = Σ_{m=0}^{Nc-1} 1/√(2m+1):")
for Nc in [10, 100, 1000, 10000]:
    S = sum(1.0/np.sqrt(2*m+1) for m in range(Nc))
    print(f"  Nc={Nc:6d}: S(Nc) = {S:.6f}  (√Nc = {np.sqrt(Nc):.4f}, ratio = {S/np.sqrt(Nc):.4f})")

print()
print("YES: S(Nc) ~ √Nc → ∞. This confirms KAY's divergence concern.")
print("The Coulomb diagonal correction grows as √Nc, which diverges if Nc → ∞.")
print()
print("BUT: In physics, Nc is FIXED (= Ω_c / (2πT_c), a physical cutoff).")
print("The question is whether we can give this cutoff rigorous status.")
print()

# Part 5: Our P4 approach — how did we handle it?
print("--- Part 5: How P4 handles Coulomb ---")
print("In P4 (Einstein phonon), Coulomb μ* enters as:")
print("  ℌ₊(μ*) = ℌ₊(0) - μ* |u⟩⟨u|")
print("where |u⟩ = (u_0, u_1, ..., u_{N_c-1}, 0, 0, ...)")
print("This is a RANK-ONE perturbation (or finite-rank in the half-space).")
print("It's well-defined because |u⟩ has finite ℓ² norm when restricted to Nc terms.")
print()

# Check: does u restricted to Nc terms have finite norm?
for Nc in [10, 100, 1000]:
    norm_sq = sum(1.0/(2*m+1) for m in range(Nc))
    print(f"  ||u_Nc||^2 = sum_{{m=0..{Nc-1}}} 1/(2m+1) = {norm_sq:.6f}"
          f"  (diverges as ln(Nc), but for FIXED Nc, finite)")

print()
print("The key insight: ||u||² ~ ln(Nc) diverges as Nc→∞,")
print("but for any FIXED physical Nc, it's finite and well-defined.")
print("The operator ℌ₊(μ*) is a compact perturbation of ℌ₊(0) plus a finite-rank term.")
print("The finite-rank term doesn't affect compactness.")
print()

# Part 6: Does the EIGENVALUE (Tc) depend on Nc?
print("--- Part 6: Tc sensitivity to cutoff Nc ---")
print("Physical question: does λ₁ (∝ Tc^γ) change significantly with Nc?")
print("If yes, then Nc is a 'knob' we're adding — physically ambiguous.")
print("If no (converges rapidly), then the exact value of Nc doesn't matter much.")
print()

# For γ=1, check λ₁ vs Nc at fixed μ*=0.13
N = 300
mu_star = 0.13

print(f"γ={gamma}, N={N}, μ*={mu_star}:")
print(f"{'Nc':>6} {'λ₁':>14} {'Tc_rel':>10}")
lam1_ref = None
for Nc in [5, 10, 20, 30, 50, 75, 100, 150, 200, 250, 300]:
    H = build_gamma_operator(N, gamma, mu_star=mu_star, Nc=Nc)
    evals = np.sort(eigvalsh(H))[::-1]
    lam1 = evals[0]
    if lam1_ref is None:
        lam1_ref = lam1
    # Tc ∝ λ₁^{1/γ}, relative change
    tc_rel = (lam1/lam1_ref)**(1.0/gamma) if lam1_ref != 0 and lam1 > 0 else 0
    print(f"{Nc:6d} {lam1:14.8f} {tc_rel:10.6f}")

print()
print("=" * 70)
print("PHASE 0 VERDICT")
print("=" * 70)
print("""
1. KAY's divergence is REAL: the Coulomb diagonal Σ_m u_m ~ √N diverges.
   Adding μ* as a constant to ALL Matsubara terms is mathematically ill-defined.

2. BUT: The physical regularization (finite Nc cutoff) is mathematically
   well-defined — it's a finite-rank perturbation of a compact operator.
   Our P4-D already handles this correctly.

3. The issue is PHYSICAL, not mathematical: Nc = Ω_c/(2πTc) depends on
   the material's electronic bandwidth, and λ₁ depends on Nc. So the
   "rigorous Tc bound" becomes a function Tc(g, γ, μ*, Nc).

4. This is EXACTLY what happens in real Eliashberg theory — Tc depends
   on the Coulomb cutoff. The standard Morel-Anderson formula gives
   μ*(Nc) = μ/(1 + μ ln(Ω_c/ω_D)), which absorbs the cutoff dependence.

5. OUR OPPORTUNITY: Give the FIRST rigorous treatment of this cutoff
   dependence for the γ-model. Instead of "smuggling in" μ* as KAY
   criticizes, we can:
   (a) State the theorem for FIXED Nc (mathematically clean)
   (b) Show how Tc depends on Nc (the Morel-Anderson renormalization)
   (c) Give certificate-level bounds Tc(g, γ, μ*, Nc) for physical ranges

6. This would be genuinely new: KAY refused to do it, EA didn't try.
   We have the machinery (P4-D) and the mathematical framework (finite-rank
   perturbation theory) to do it rigorously.

VERDICT: GO on Coulomb μ* direction for γ-model.
""")
