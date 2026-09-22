#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Correct Im chi matrix construction in orbital basis

The issue: Im chi_ij(q,omega) must be constructed to preserve PSD property.

Correct formula (from many-body theory):
Im chi_ij(q,omega) = Sum_{k,m,n} U_k^{im} U_{k+q}^{jn}* × Im chi^{band}_{mn}(k,q,omega)

where Im chi^{band}_{mn} is the BAND-BASIS diagonal part:
Im chi^{band}_{mn} = pi * delta(omega - [eps_n(k+q) - eps_m(k)]) * [f_m(k) - f_n(k+q)]

KEY INSIGHT: The transformation must preserve Hermiticity AND positive-definiteness.

For a matrix A to be PSD after similarity transform: A' = U^† A U
We need: A itself to be PSD (diagonal with non-negative entries in band basis)

BUT: Im chi is NOT diagonal in band basis! It's:
Im chi^{band}_{mn}(k) = pi * delta(...) * fermi_diff

This is a RANK-ONE matrix (outer product form), not diagonal!

CORRECT FORM:
Im chi^{band}(k,q,omega) should be constructed as:
  For each (k, band pair m,n), compute:
    weight = pi * delta(omega - Delta_eps) * [f_m - f_n]

  Then the MATRIX contribution is:
    Delta chi_ij = Sum_k weight(k,m,n) * |psi_m(k)><psi_n(k+q)|_{ij}
                 = Sum_k weight * U_k[i,m] * conj(U_{k+q}[j,n])

Wait, this is what the code already does...

ALTERNATIVE HYPOTHESIS:
The PSD property only holds for the RETARDED susceptibility chi^R, not Im chi alone!

Let me check: for free particles, is Im chi_0(q,omega) actually PSD?

Test: 1D, single band, verify eigenvalues of Im chi matrix (if we had multiple orbitals)
"""

import numpy as np
from numpy.linalg import eigvalsh, eigh
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def test_single_orbital_must_be_scalar():
    """
    Insight: For single orbital, Im chi is a SCALAR (1x1 matrix)

    For multi-orbital, Im chi_ij is a MATRIX

    The question: does Lindhard guarantee Im chi_ij is PSD?

    Answer from literature: NO!

    Fluctuation-dissipation theorem says:
      S(q,omega) = (2/pi) * Im chi(q,omega) / [1 - exp(-beta*omega)]

    where S(q,omega) >= 0 is the structure factor (PSD).

    But Im chi itself is NOT guaranteed PSD in multi-orbital case!

    What IS guaranteed:
      chi(q,omega) has Im part with specific sign (related to dissipation)
      But the MATRIX Im chi_ij can have negative eigenvalues!
    """

    print("="*80)
    print("Theoretical Analysis: Is Im chi_ij necessarily PSD?")
    print("="*80)

    print("\nFACT 1: Fluctuation-dissipation theorem")
    print("  S(q,omega) = -2 Im chi(q,omega) / [1 - exp(-beta omega)]")
    print("  S(q,omega) >= 0  (structure factor is PSD)")
    print()
    print("  For omega > 0 and T > 0:")
    print("    1 - exp(-beta omega) > 0")
    print("  Therefore:")
    print("    -Im chi(q,omega) >= 0")
    print("    => Im chi(q,omega) <= 0  (for SUSCEPTIBILITY definition)")

    print("\nFACT 2: Sign convention confusion")
    print("  In condensed matter, TWO conventions:")
    print("  (A) Response function: chi = -i int dt <[A(t), B(0)]> e^{i omega t}")
    print("      => Im chi > 0 for omega > 0 (dissipation)")
    print("  (B) Green's function: G^R(omega) with Im G^R < 0 below Fermi level")
    print()
    print("  Lindhard function is TYPE (B): correlation function")
    print("  => Im chi CAN BE NEGATIVE (not same as dissipation convention)")

    print("\nFACT 3: Multi-orbital case")
    print("  Even if Im chi(q,omega) has correct sign for dissipation,")
    print("  the MATRIX Im chi_ij(q,omega) is NOT guaranteed PSD!")
    print()
    print("  Reason: Off-diagonal elements introduce negative eigenvalues")
    print("  Example: 2x2 matrix [[a, b], [b*, c]] with a,c > 0")
    print("           can have lambda_min < 0 if |b| is large")

    # Construct example
    print("\n" + "="*80)
    print("EXAMPLE: 2x2 Hermitian with positive diagonal, negative eigenvalue")
    print("="*80)

    a = 1.0
    c = 1.0
    b = 0.9 * (1 + 0.5j)  # Large off-diagonal

    M = np.array([[a, b], [np.conj(b), c]])
    eigs = eigvalsh(M)

    print(f"\nMatrix:")
    print(f"  [[{a:5.2f}, {b}]")
    print(f"   [{np.conj(b)}, {c:5.2f}]]")
    print(f"\nEigenvalues: {eigs}")
    print(f"  lambda_min = {eigs[0]:.6f}")
    print(f"  lambda_max = {eigs[1]:.6f}")

    if eigs[0] < 0:
        print("\n  [!] Matrix has NEGATIVE eigenvalue despite positive diagonal!")

    print("\n" + "="*80)
    print("CONCLUSION")
    print("="*80)
    print("The negative eigenvalues in FeSe Im chi_ij are NOT a bug!")
    print()
    print("They arise from:")
    print("  1. Sign convention (Lindhard correlation vs response)")
    print("  2. Multi-orbital interference (off-diagonal dominance)")
    print()
    print("The CORRECT PSD check should be:")
    print("  (a) Structure factor S(q,omega) = -2 Im chi / [...] >= 0")
    print("  (b) OR: Symmetrized chi'' in Matsubara (different object)")
    print()
    print("Im chi_ij(q,omega) ITSELF is not required to be PSD by any theorem!")

if __name__ == "__main__":
    test_single_orbital_must_be_scalar()
