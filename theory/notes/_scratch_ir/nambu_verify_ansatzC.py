"""
Verify Ansatz C (derived from raw Nambu-Eliashberg xi-integral, see
nambu_weight_derivation_2026-09-16.md) against Lemma 4's a_n >= 0 claim.

Ansatz C:
  L1[g,f]_n = (1/(2n+1)) * 2 * sum_{k=1}^n g(w_k) f(k)      (g on SOURCE index k)
  P[g,f]    = L2+L3 with u_n^2 = g(w_n)/(2n+1)              (same as Ansatz A for P)

Claim: a_n := L1[g,f]_n - rho(f) g(w_n)/(2n+1) >= 0 for n>=1 whenever
g is NON-INCREASING (Fermi-surface-enhanced / van Hove-type DOS shape),
and can go negative for g NON-DECREASING (pseudogap-type).

From scratch, no repo imports. OMP_NUM_THREADS=1.
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np

def bose_kernel(k, omegas, weights):
    # f(k) = sum_i weights_i * omega_i^2/(omega_i^2+k^2), k in units of 2*pi*T
    k = np.asarray(k, dtype=float)
    out = np.zeros_like(k)
    for w, om in zip(weights, omegas):
        out += w * om**2 / (om**2 + k**2)
    return out

def rho_and_fhat_pi(f_at, omegas, weights):
    # rho(f) = f(0) - fhat(pi) = f(0) - sum_i w_i * pi*om_i/sinh(pi*om_i)
    f0 = f_at(np.array([0.0]))[0]
    fhat_pi = sum(w * np.pi*om/np.sinh(np.pi*om) for w, om in zip(weights, omegas))
    return f0 - fhat_pi

def a_n_ansatzC(N, f_at, g_of_n, omegas, weights):
    rho = rho_and_fhat_pi(f_at, omegas, weights)
    a = np.zeros(N)
    fk = f_at(np.arange(0, N+1))  # fk[k] = f(k)
    gk = np.array([g_of_n(k) for k in range(N+1)])
    csum = 0.0
    for n in range(N):
        if n >= 1:
            csum += gk[n]*fk[n]  # add f(n)*g(n) to running sum for k=1..n
        L1n = (2.0*csum)/(2*n+1)
        a[n] = L1n - rho*gk[n]/(2*n+1)
    return a, rho

if __name__ == "__main__":
    omegas = [1.0]
    weights = [1.0]
    f_at = lambda k: bose_kernel(k, omegas, weights)
    N = 400

    profiles = {
        "flat (g=1)": lambda n: 1.0,
        "van-Hove-like DECREASING g: 0.3+0.7/(1+n/20)": lambda n: 0.3+0.7/(1.0+n/20.0),
        "pseudogap-like INCREASING g: 1-0.7/(1+n/20)": lambda n: 1.0-0.7/(1.0+n/20.0),
        "step DOWN at n=5 (g=0.05 for n<5 else 1)": lambda n: (0.05 if n<5 else 1.0),
        "step UP at n=5 (g=1 for n<5 else 0.05)": lambda n: (1.0 if n<5 else 0.05),
    }

    for name, gfun in profiles.items():
        a, rho = a_n_ansatzC(N, f_at, gfun, omegas, weights)
        print(f"{name:55s}  rho={rho:.4f}  min(a[1:])={a[1:].min():+.6e}  a0={a[0]:+.6e}  (a0 always = -rho*g0 by construction)")
