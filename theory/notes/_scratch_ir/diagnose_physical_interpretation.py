"""Question 5: Physical interpretation — connection to Anderson-Morel / Tolmachev mechanism.

Check if concavity failure correlates with:
  - Sign changes in v (cone exit)
  - Sign changes occurring near n ~ Nc (cutoff edge)
  - Spatial structure of c2 = <w2, M' v> vs uc's truncation boundary

Strictness: numerical-evidence (float64, finite N).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

N = 300
_n = np.arange(N)
_u = 1.0 / np.sqrt(2.0 * _n + 1.0)
_ks = np.arange(2 * N + 1, dtype=float)
_absd = np.abs(_n[:, None] - _n[None, :])
_sumd = _n[:, None] + _n[None, :] + 1
_uu = _u[:, None] * _u[None, :]


def O(f):
    A = (f[_absd] + f[_sumd]) * _uu
    L1 = np.zeros(N)
    L1[1:] = np.cumsum(2.0 * f[1:N])[:N - 1] / (2 * _n[1:] + 1)
    np.fill_diagonal(A, -L1 + f[2 * _n + 1] * _u ** 2)
    return A


def f_ein(x):
    f = x / (_ks ** 2 + x)
    f[0] = 1.0
    return f


def coulomb(mu_t, Nc):
    uc = np.where(_n < (Nc if Nc else N), _u, 0.0)
    return 2.0 * mu_t * np.outer(uc, uc), uc


print("=" * 100)
print("QUESTION 5: Physical interpretation — sign changes and spatial structure")
print("=" * 100)
print()

# Compare two cases: (mu_t=0, Nc=None) vs (mu_t=0.13, Nc=6) at x=0.09
cases = [
    ("Baseline: mu*=0", 0.0, None, 0.09),
    ("Failure case: mu_t=0.13, Nc=6", 0.13, 6, 0.09),
    ("Mild case: mu_t=0.13, Nc=45", 0.13, 45, 0.09),
]

fig, axes = plt.subplots(3, 3, figsize=(15, 12))
fig.suptitle("Spatial structure analysis", fontsize=14)

for row_idx, (label, mu_t, Nc, x) in enumerate(cases):
    s = _ks ** 2
    f = f_ein(x)
    A_ = O(f)
    B_ = O(f ** 2)
    Cm, uc = coulomb(mu_t, Nc) if mu_t > 0 else (np.zeros((N, N)), np.zeros(N))
    M = A_ - Cm
    lam, W = eigh(M)
    v = W[:, -1]
    w2 = W[:, -2]
    gam = lam[-1] - lam[-2]

    Mp = (A_ - B_) / x
    Psi = float(v @ (O(s / (s + x) ** 3) @ v))
    cj = W[:, :-1].T @ (Mp @ v)
    R = float(2.0 * np.sum(cj ** 2 / (lam[-1] - lam[:-1])))
    E = 2 * x * Psi
    ratio = (x * R) / E if E > 1e-12 else 0.0

    sign_changes = np.sum(np.diff(np.sign(v)) != 0)

    # Find sign change locations
    sign_change_locs = np.where(np.diff(np.sign(v)) != 0)[0]

    # Compute spatial structure of M'v
    Mpv = Mp @ v

    print(f"{label}:")
    print(f"  gap = {gam:.4f}, ratio = {ratio:.3f}, sign_changes = {sign_changes}")
    if Nc:
        print(f"  Nc = {Nc}, sign change locations: {sign_change_locs}")
    else:
        print(f"  sign change locations: {sign_change_locs}")

    # Plot v
    ax = axes[row_idx, 0]
    ax.plot(_n, v, 'b-', linewidth=0.8)
    ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    if Nc:
        ax.axvline(Nc, color='r', linestyle='--', linewidth=1.5, label=f'Nc={Nc}')
        ax.legend(fontsize=8)
    for loc in sign_change_locs:
        ax.axvline(loc, color='orange', linestyle=':', linewidth=1, alpha=0.7)
    ax.set_xlim(0, 100)
    ax.set_xlabel('n', fontsize=9)
    ax.set_ylabel('v[n]', fontsize=9)
    ax.set_title(f'{label}\nTop eigenvector v', fontsize=9)
    ax.grid(True, alpha=0.3)

    # Plot w2
    ax = axes[row_idx, 1]
    ax.plot(_n, w2, 'g-', linewidth=0.8)
    ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    if Nc:
        ax.axvline(Nc, color='r', linestyle='--', linewidth=1.5)
    ax.set_xlim(0, 100)
    ax.set_xlabel('n', fontsize=9)
    ax.set_ylabel('w2[n]', fontsize=9)
    ax.set_title('Second eigenvector w2', fontsize=9)
    ax.grid(True, alpha=0.3)

    # Plot M'v
    ax = axes[row_idx, 2]
    ax.plot(_n, Mpv, 'm-', linewidth=0.8)
    ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    if Nc:
        ax.axvline(Nc, color='r', linestyle='--', linewidth=1.5)
    for loc in sign_change_locs:
        ax.axvline(loc, color='orange', linestyle=':', linewidth=1, alpha=0.7)
    ax.set_xlim(0, 100)
    ax.set_xlabel('n', fontsize=9)
    ax.set_ylabel("(M'v)[n]", fontsize=9)
    ax.set_title(f"M'v (drives level repulsion)\nc2 = <w2, M'v> = {cj[-1]:.3e}", fontsize=9)
    ax.grid(True, alpha=0.3)

    print()

plt.tight_layout()
plt.savefig('D:/physics/superconductivity/theory/notes/_scratch_ir/concavity_spatial_structure.png', dpi=150)
print("Saved figure: concavity_spatial_structure.png")
print()

# Detailed sign change analysis
print("=" * 100)
print("Sign change analysis vs Nc")
print("=" * 100)
print()

mu_t = 0.13
x = 0.09
print(f"Fixed: mu_t = {mu_t}, x = {x} (varpi = {np.sqrt(x):.3f})")
print(f"\n{'Nc':>4} {'gap':>8} {'ratio':>7} {'sign_chg':>9} {'locations':>30} {'near Nc?':>10}")

for Nc in [3, 4, 6, 8, 12, 20, 30, 45, 64]:
    s = _ks ** 2
    f = f_ein(x)
    A_ = O(f)
    B_ = O(f ** 2)
    Cm, uc = coulomb(mu_t, Nc)
    M = A_ - Cm
    lam, W = eigh(M)
    v = W[:, -1]
    gam = lam[-1] - lam[-2]

    Mp = (A_ - B_) / x
    Psi = float(v @ (O(s / (s + x) ** 3) @ v))
    cj = W[:, :-1].T @ (Mp @ v)
    R = float(2.0 * np.sum(cj ** 2 / (lam[-1] - lam[:-1])))
    E = 2 * x * Psi
    ratio = (x * R) / E if E > 1e-12 else 0.0

    sign_changes = np.sum(np.diff(np.sign(v)) != 0)
    sign_change_locs = np.where(np.diff(np.sign(v)) != 0)[0]

    # Check if any sign change is within 3 indices of Nc
    near_nc = any(abs(loc - Nc) <= 3 for loc in sign_change_locs) if len(sign_change_locs) > 0 else False

    locs_str = ",".join(map(str, sign_change_locs[:5]))  # show first 5
    if len(sign_change_locs) > 5:
        locs_str += "..."

    print(f"{Nc:>4} {gam:>8.4f} {ratio:>7.3f} {sign_changes:>9} {locs_str:>30} {'YES' if near_nc else 'no':>10}")
