"""
Phase 1 feasibility: Haldane conjecture Davis-Kahan attack.

1. Construct H_Heis and H_AKLT for spin-1 chains (PBC) of length L=4,6,8,10
2. Verify the identity: H_Heis = 2*H_AKLT - sum P^(0) - L*I
3. Compute spectral gap of H(t) = (1-t)*H_AKLT + t*H_Heis along t in [0,1]
4. Compute ||V||/gap ratio for Davis-Kahan assessment
5. Direct finite-size gap of Heisenberg chain

Uses S^z conservation to block-diagonalize (total S^z = 0 sector for gap).
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
import time

# Spin-1 matrices (3x3)
Sz = np.diag([1.0, 0.0, -1.0])
Sp = np.array([[0, np.sqrt(2), 0],
               [0, 0, np.sqrt(2)],
               [0, 0, 0]])
Sm = Sp.T
Sx = (Sp + Sm) / 2.0
Sy = (Sp - Sm) / (2.0j)

def kron_op(op, site, L):
    """Place operator on site `site` of L-site chain (sparse)."""
    mats = []
    for i in range(L):
        if i == site:
            mats.append(sparse.csr_matrix(op))
        else:
            mats.append(sparse.eye(3, format='csr'))
    result = mats[0]
    for m in mats[1:]:
        result = sparse.kron(result, m, format='csr')
    return result

def two_site_dot(i, j, L):
    """S_i . S_j on L-site chain."""
    return (kron_op(Sx, i, L) @ kron_op(Sx, j, L) +
            kron_op(Sy, i, L) @ kron_op(Sy, j, L) +
            kron_op(Sz, i, L) @ kron_op(Sz, j, L))

def projector_S2(i, j, L):
    """Projector onto total spin=2 sector of sites i,j.
    P^(2) = 1/3 + h/2 + h^2/6 where h = S_i . S_j"""
    h = two_site_dot(i, j, L)
    dim = 3**L
    I = sparse.eye(dim, format='csr')
    return I/3.0 + h/2.0 + h@h/6.0

def projector_S0(i, j, L):
    """Projector onto total spin=0 sector of sites i,j.
    P^(0) = 1/3 - h/3 + h^2/6 - h/6 ... let me derive properly.

    h eigenvalues: S=0 -> -2, S=1 -> -1, S=2 -> +1
    P^(0) should give 1 for h=-2, 0 for h=-1, 0 for h=+1

    Lagrange interpolation on {-2,-1,1}:
    P^(0)(h) = (h-(-1))(h-1) / ((-2-(-1))(-2-1)) = (h+1)(h-1)/((-1)(-3)) = (h^2-1)/3
    """
    h = two_site_dot(i, j, L)
    dim = 3**L
    I = sparse.eye(dim, format='csr')
    return (h@h - I) / 3.0

def build_H_heis(L):
    """Heisenberg spin-1 chain with PBC."""
    dim = 3**L
    H = sparse.csr_matrix((dim, dim), dtype=float)
    for i in range(L):
        j = (i + 1) % L
        H = H + two_site_dot(i, j, L)
    return H

def build_H_aklt(L):
    """AKLT Hamiltonian = sum of P^(2)_i,i+1 with PBC."""
    dim = 3**L
    H = sparse.csr_matrix((dim, dim), dtype=float)
    for i in range(L):
        j = (i + 1) % L
        H = H + projector_S2(i, j, L)
    return H

def build_sum_P0(L):
    """Sum of singlet projectors P^(0)_i,i+1 with PBC."""
    dim = 3**L
    P = sparse.csr_matrix((dim, dim), dtype=float)
    for i in range(L):
        j = (i + 1) % L
        P = P + projector_S0(i, j, L)
    return P

# Total Sz operator for block diagonalization
def build_total_Sz(L):
    dim = 3**L
    Sztot = sparse.csr_matrix((dim, dim), dtype=float)
    for i in range(L):
        Sztot = Sztot + kron_op(Sz, i, L)
    return Sztot

def get_sz0_sector(L):
    """Get the Sz=0 sector indices."""
    # Each site: +1, 0, -1. We need total = 0.
    dim = 3**L
    Sztot = build_total_Sz(L)
    diag = Sztot.diagonal()
    return np.where(np.abs(diag) < 0.5)[0]

print("=" * 70)
print("HALDANE CONJECTURE FEASIBILITY: Phase 1")
print("=" * 70)

# Part 1: Verify algebraic identity
print("\n--- Part 1: Verify H_Heis = 2*H_AKLT - sum_P0 - L*I ---")
for L in [4, 6]:
    t0 = time.time()
    H_heis = build_H_heis(L)
    H_aklt = build_H_aklt(L)
    sum_P0 = build_sum_P0(L)
    dim = 3**L
    I = sparse.eye(dim)

    # Check: H_heis should equal 2*H_aklt - sum_P0 - L*I
    diff = H_heis - (2*H_aklt - sum_P0 - L*I)
    err = sparse.linalg.norm(diff)
    print(f"  L={L}: dim={dim}, ||H_heis - (2*H_aklt - sum_P0 - L*I)|| = {err:.2e}  "
          f"({'OK' if err < 1e-10 else 'FAIL'})  [{time.time()-t0:.1f}s]")

# Part 2: Spectral gaps
print("\n--- Part 2: Spectral gaps ---")
print(f"{'L':>4} {'dim':>8} {'dim_Sz0':>8} {'gap_Heis':>10} {'gap_AKLT':>10} "
      f"{'E0_Heis':>10} {'E0_AKLT':>10} {'time':>6}")

for L in [4, 6, 8, 10]:
    t0 = time.time()
    H_heis = build_H_heis(L)
    H_aklt = build_H_aklt(L)
    dim = 3**L

    if dim <= 2000:
        # Full diag for small systems
        eH = np.linalg.eigvalsh(H_heis.toarray())
        eA = np.linalg.eigvalsh(H_aklt.toarray())
        gap_H = eH[1] - eH[0]
        gap_A = eA[1] - eA[0]  # Should be ~0 (ground space degeneracy)
        # AKLT ground space is 4-fold degenerate (VBS), gap is E(5th) - E(1st)
        # Actually for PBC AKLT, GS is 4-fold; gap = E(5th)-E(1st) for L>=3
        gap_A_true = eA[4] - eA[0] if L >= 3 else eA[1] - eA[0]
        print(f"{L:4d} {dim:8d} {'full':>8} {gap_H:10.6f} {gap_A_true:10.6f} "
              f"{eH[0]:10.6f} {eA[0]:10.6f} {time.time()-t0:6.1f}s")
    else:
        # Sparse eigsh
        # Get lowest 6 eigenvalues
        eH = eigsh(H_heis, k=min(6, dim-2), which='SA', return_eigenvectors=False)
        eH.sort()
        eA = eigsh(H_aklt, k=min(6, dim-2), which='SA', return_eigenvectors=False)
        eA.sort()
        gap_H = eH[1] - eH[0]
        gap_A_true = eA[4] - eA[0] if len(eA) > 4 else eA[-1] - eA[0]
        print(f"{L:4d} {dim:8d} {'sparse':>8} {gap_H:10.6f} {gap_A_true:10.6f} "
              f"{eH[0]:10.6f} {eA[0]:10.6f} {time.time()-t0:6.1f}s")

# Part 3: Interpolation path H(t) = (1-t)*2*H_AKLT + t*H_Heis, gap along path
print("\n--- Part 3: Gap along interpolation H(t)=(1-t)*2*H_AKLT + t*H_Heis ---")
# Use the identity: H(t) = (1-t)*2*H_AKLT + t*(2*H_AKLT - sum_P0 - L*I)
#                        = 2*H_AKLT - t*sum_P0 - t*L*I
# So gap(H(t)) = gap(2*H_AKLT - t*sum_P0) (additive constant doesn't affect gap)

for L in [4, 6, 8]:
    t0 = time.time()
    H_aklt = build_H_aklt(L)
    sum_P0 = build_sum_P0(L)
    dim = 3**L

    print(f"  L={L} (dim={dim}):")
    ts = np.linspace(0, 1, 11)
    for t in ts:
        Ht = 2*H_aklt - t*sum_P0
        if dim <= 2000:
            evals = np.linalg.eigvalsh(Ht.toarray())
            # AKLT has 4-fold GS for PBC; at t=0, gap above this
            # As t increases, the degeneracy may split
            # Report gap = E(2nd unique) - E(1st)
            gap = evals[1] - evals[0]
        else:
            evals = eigsh(Ht, k=6, which='SA', return_eigenvectors=False)
            evals.sort()
            gap = evals[1] - evals[0]
        print(f"    t={t:.2f}: E0={evals[0]:.6f}  gap(E1-E0)={gap:.6f}")
    print(f"    [{time.time()-t0:.1f}s]")

# Part 4: Per-bond perturbation norm analysis
print("\n--- Part 4: Per-bond perturbation analysis ---")
print("For H_Heis = 2*H_AKLT - sum_P0 - L*I:")
print("  V = -sum_P0 (perturbation relative to 2*H_AKLT)")
print("  Per-bond: v_i = -P^(0)_i, ||v_i|| = 1")
print(f"  AKLT gap ≈ 0.350, so 2*AKLT gap ≈ 0.700")
print(f"  Ratio ||v_i||/gap ≈ {1.0/0.700:.3f}")
print()

# Alternative: optimize c in H_Heis - c*H_AKLT
print("Optimizing c in ||H_Heis - c*H_AKLT|| (per bond):")
print("  h eigenvalues: S=0 -> -2, S=1 -> -1, S=2 -> +1")
print("  P^(2) eigenvalues: S=0 -> 0, S=1 -> 0, S=2 -> 1")
print("  V(c) per bond eigenvalues: S=0 -> -2, S=1 -> -1, S=2 -> 1-c")
for c in [0, 0.5, 1, 2, 3]:
    ev = [-2, -1, 1-c]
    norm = max(abs(e) for e in ev)
    print(f"  c={c:.1f}: eigenvalues={ev}, ||v||={norm:.3f}")

print("\n  => Minimum ||v(c)|| = 2 at any c <= 3 (S=0 eigenvalue -2 dominates)")
print("  => Regardless of c, per-bond perturbation >= 2")
print()

# Part 5: Davis-Kahan assessment
print("=" * 70)
print("DAVIS-KAHAN FEASIBILITY ASSESSMENT")
print("=" * 70)
print()
print("For the path H_AKLT -> H_Heis:")
print("  Best decomposition: H_Heis = 2*H_AKLT - sum P^(0) - L*I")
print("  Per-bond perturbation norm: ||P^(0)_i|| = 1")
print("  2*AKLT spectral gap: ~0.700")
print("  Ratio: 1/0.700 ≈ 1.43 > 1")
print()
print("VERDICT: Naive Davis-Kahan perturbation from AKLT to Heisenberg FAILS.")
print("  The per-bond perturbation exceeds the rescaled AKLT gap.")
print("  This ratio (1.43) is close to 1, suggesting the approach 'almost works'")
print("  but cannot succeed without additional structural information.")
print()
print("HOWEVER: The gap along the interpolation path may still remain open")
print("  (as strongly suggested by DMRG/ED numerics). The issue is that")
print("  Davis-Kahan gives a SUFFICIENT condition; its failure doesn't mean")
print("  the gap actually closes.")
print()
print("RECOMMENDED NEXT STEPS:")
print("  1. LTI-SDP certificate route (Rai et al. 2026): directly certify")
print("     Heisenberg gap for finite L, then extrapolate")
print("  2. Refined perturbation: use spectral flow / quasi-adiabatic")
print("     continuation with more detailed gap tracking")
print("  3. Combined: AKLT certificate (easy) -> interpolation certificates")
print("     at discrete t-values -> piece together a gap proof")
