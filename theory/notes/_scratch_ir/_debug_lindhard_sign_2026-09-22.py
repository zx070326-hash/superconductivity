#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Debug: Lindhard function Im chi sign convention

Analytic derivation check for imaginary part on real axis.

Starting from Matsubara:
chi_0(q, i*omega_n) = (1/beta) Sum_k [f_k - f_{k+q}] / [i*omega_n + eps_k - eps_{k+q}]

Analytic continuation i*omega_n -> omega + i*eta:
chi_0(q, omega+i*eta) = Sum_k [f_k - f_{k+q}] / [omega + i*eta + eps_k - eps_{k+q}]

Define Delta_eps = eps_{k+q} - eps_k, then:
chi_0(q, omega+i*eta) = Sum_k [f_k - f_{k+q}] / [omega - Delta_eps + i*eta]

Taking imaginary part:
Im chi_0 = Sum_k [f_k - f_{k+q}] * Im[1/(omega - Delta_eps + i*eta)]

Using: Im[1/(x + i*eta)] = -eta / (x^2 + eta^2) -> -pi * delta(x) as eta->0

Therefore:
Im chi_0(q, omega) = -pi * Sum_k [f_k - f_{k+q}] * delta(omega - Delta_eps)

KEY SIGN CHECK:
- For omega > 0 and Delta_eps > 0 (excited states):
  If f_k > f_{k+q} (k on Fermi surface, k+q above), contribution is POSITIVE * (-pi) = NEGATIVE

BUT: Fluctuation-dissipation theorem says Im chi MUST be PSD for omega > 0.

RESOLUTION: The standard physics convention uses chi with a DIFFERENT sign:
chi^R(omega) = chi_0(i*omega_n -> omega + i*0^+)

And the FDT form is:
S(q,omega) = (1 - exp(-beta*omega))^{-1} * Im chi(q,omega)

Where S(q,omega) >= 0 is guaranteed.

For T -> 0 and omega > 0:
Im chi(q,omega) / omega = pi * Sum_k delta(eps_k) delta(eps_{k+q} - eps_k - omega) >= 0

This suggests the formula should be:
Im chi = +pi * Sum_k delta(omega - Delta_eps) * [f_k - f_{k+q}]

NOT with a negative sign!

Let me verify this with a simple 1D test case.
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def test_1d_free_electron_gas():
    """
    Test case: 1D free electron gas

    eps_k = k^2 / (2m), take m=1, k in [-pi, pi]
    T -> 0, mu = 0 (half filling)

    Im chi_0(q, omega) should be >= 0 for omega > 0
    """

    print("="*80)
    print("1D Free Electron Gas Test - Sign Convention Check")
    print("="*80)

    # k grid
    N_k = 200
    k_arr = np.linspace(-np.pi, np.pi, N_k, endpoint=False)
    dk = 2*np.pi / N_k

    # Dispersion
    eps_k = k_arr**2 / 2.0

    # Fermi function (T=0.01)
    mu = 0.0
    T = 0.01
    f_k = 1.0 / (1.0 + np.exp((eps_k - mu) / T))

    print(f"\nParameters:")
    print(f"  N_k = {N_k}")
    print(f"  mu = {mu}")
    print(f"  T = {T}")

    # Test q = 0.5*pi, omega = 0.5
    q = 0.5 * np.pi
    omega = 0.5
    eta = 0.02

    print(f"\nTest point: q = {q:.3f}, omega = {omega:.3f}")

    # Compute Im chi with BOTH sign conventions
    im_chi_negative_pi = 0.0
    im_chi_positive_pi = 0.0

    for k in k_arr:
        eps_kq = ((k + q) % (2*np.pi) - np.pi)**2 / 2.0  # Periodic BC

        f_kq_val = 1.0 / (1.0 + np.exp((eps_kq - mu) / T))

        delta_eps = eps_kq - eps_k[np.where(k_arr == k)[0][0]]

        # Lorentzian broadening
        lorentzian = (eta / np.pi) / ((omega - delta_eps)**2 + eta**2)

        fermi_diff = f_k[np.where(k_arr == k)[0][0]] - f_kq_val

        # Two conventions
        im_chi_negative_pi += (-np.pi) * lorentzian * fermi_diff * dk / (2*np.pi)
        im_chi_positive_pi += (+np.pi) * lorentzian * fermi_diff * dk / (2*np.pi)

    print(f"\nResults:")
    print(f"  Im chi (with -pi): {im_chi_negative_pi:+.6f}")
    print(f"  Im chi (with +pi): {im_chi_positive_pi:+.6f}")

    print(f"\nFluctuation-dissipation theorem requires: Im chi > 0 for omega > 0")

    if im_chi_positive_pi > 0:
        print("[OK] +pi convention gives POSITIVE Im chi")
        print("=> Correct formula: Im chi = +pi * delta(...) * [f_k - f_{k+q}]")
    elif im_chi_negative_pi > 0:
        print("[OK] -pi convention gives POSITIVE Im chi")
        print("=> Correct formula: Im chi = -pi * delta(...) * [f_k - f_{k+q}]")
    else:
        print("[ERROR] Both conventions give wrong sign!")

    # Extended test: scan omega
    print("\n" + "="*80)
    print("Extended test: omega scan")
    print("="*80)

    omega_arr = np.linspace(0.01, 2.0, 20)
    im_chi_arr_neg = []
    im_chi_arr_pos = []

    for omega in omega_arr:
        im_neg = 0.0
        im_pos = 0.0

        for k in k_arr:
            eps_kq = ((k + q) % (2*np.pi) - np.pi)**2 / 2.0
            f_kq_val = 1.0 / (1.0 + np.exp((eps_kq - mu) / T))
            delta_eps = eps_kq - eps_k[np.where(k_arr == k)[0][0]]
            lorentzian = (eta / np.pi) / ((omega - delta_eps)**2 + eta**2)
            fermi_diff = f_k[np.where(k_arr == k)[0][0]] - f_kq_val

            im_neg += (-np.pi) * lorentzian * fermi_diff * dk / (2*np.pi)
            im_pos += (+np.pi) * lorentzian * fermi_diff * dk / (2*np.pi)

        im_chi_arr_neg.append(im_neg)
        im_chi_arr_pos.append(im_pos)

    im_chi_arr_neg = np.array(im_chi_arr_neg)
    im_chi_arr_pos = np.array(im_chi_arr_pos)

    print(f"\nomega range: [{omega_arr[0]:.2f}, {omega_arr[-1]:.2f}]")
    print(f"\nWith -pi convention:")
    print(f"  Min Im chi: {im_chi_arr_neg.min():+.6f}")
    print(f"  Max Im chi: {im_chi_arr_neg.max():+.6f}")
    print(f"  Negative values: {np.sum(im_chi_arr_neg < 0)} / {len(omega_arr)}")

    print(f"\nWith +pi convention:")
    print(f"  Min Im chi: {im_chi_arr_pos.min():+.6f}")
    print(f"  Max Im chi: {im_chi_arr_pos.max():+.6f}")
    print(f"  Negative values: {np.sum(im_chi_arr_pos < 0)} / {len(omega_arr)}")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(omega_arr, im_chi_arr_neg, 'o-', label='Im chi (with -pi)', color='red')
    ax.plot(omega_arr, im_chi_arr_pos, 's-', label='Im chi (with +pi)', color='blue')
    ax.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax.set_xlabel('omega')
    ax.set_ylabel('Im chi_0(q, omega)')
    ax.set_title(f'1D Free Electron Gas: q = {q:.2f}')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('D:/physics/superconductivity/theory/notes/_scratch_ir/_debug_lindhard_sign.png', dpi=150)
    print(f"\nPlot saved to _debug_lindhard_sign.png")

    print("\n" + "="*80)
    print("CONCLUSION")
    print("="*80)

    if np.all(im_chi_arr_pos >= -1e-6):
        print("[DEFINITIVE] The correct sign is: Im chi = +pi * delta(...) * [f_k - f_{k+q}]")
        print("Current implementation uses -pi => WRONG SIGN, must flip!")
    elif np.all(im_chi_arr_neg >= -1e-6):
        print("[DEFINITIVE] The correct sign is: Im chi = -pi * delta(...) * [f_k - f_{k+q}]")
        print("Current implementation is CORRECT")
    else:
        print("[INCONCLUSIVE] Neither convention works - deeper issue in formula")

if __name__ == "__main__":
    test_1d_free_electron_gas()
