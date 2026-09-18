"""[analytic attempt + numerical-evidence] Reconcile three claims from lam2bound.

(A) The "factor 2 in Phi_s" claim.  Decisive test: Phi_s must satisfy
    int_0^inf Phi_s dsigma = Theta(s), because Theta IS that integral.
    If a spurious 2 sits in the denominator the integral gives Theta/2.

(B) The "D_nn = w_n vs w_n + l_G(n)" claim.  Hypothesis: lam2bound's H was
    built WITHOUT the L1 diagonal -(2/(2n+1))sum_{k<=n} phi(k) = -l_G(n).
    Test both builds and see which diagonal each one forces.

(C) Finitization WITHOUT needing e_n monotone.  We never need delta_n monotone:
    the bound actually used is the ANALYTIC ENVELOPE
        delta_n <= E(s_n) := Thetatil(s_n) + c/s_n ,  c = (D^2/8)(C_R + V2[G]/2)
    and E is decreasing wherever Thetatil is (c/s always is).  So the only
    missing analytic piece is  Thetatil'(s) < 0 for s >= some explicit s1.
    Thetatil'(s) = N(s)/[s^2 (1+s^2)^2] with
        N(s) = -s^3 ln(1+s^2) + s^3 + 3 s^2 arctan(s) - s + arctan(s).
    Using arctan(s) <= pi/2 - 1/s + 1/(3 s^3)  (from arctan x >= x - x^3/3,
    x = 1/s in (0,1]) and arctan(s) <= pi/2:
        N(s) <= -psi(s),  psi(s) := s^3(ln(1+s^2)-1) - (3pi/2)s^2 + 4s - 1/s - pi/2.
    Claim: psi(3) > 0 and psi'(s) > 0 for s >= 3  =>  N(s) < 0 for all s >= 3.

Run: OMP_NUM_THREADS=1 python _chk_reconcile_lam2bound.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import mpmath as mp

mp.mp.dps = 30


def G(y):
    return 1.0 / (1.0 + y * y)


def phi_of(w):
    return lambda k: w * w / (k * k + w * w)


def build_H(w, N, with_L1):
    """Repo operator; with_L1=False drops the -l_G(n) diagonal (the suspected bug)."""
    n = np.arange(N)
    I, J = n[:, None], n[None, :]
    pv = phi_of(w)(np.arange(0, 2 * N + 2, dtype=float))
    T = pv[np.abs(I - J)].copy()
    np.fill_diagonal(T, 0.0)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    M = (T + pv[I + J + 1]) * u[:, None] * u[None, :]
    if with_L1:
        M[n, n] -= 2.0 * np.concatenate(([0.], np.cumsum(pv[1:N]))) / (2 * n + 1.0)
    return M


def pieces(w, N):
    n = np.arange(N, dtype=float)
    s = (n + 0.5) / w
    p = 1.0 / np.sqrt(2.0 * n + 1.0)
    tau = p * G(s)
    Rm = (G(s[:, None] - s[None, :]) + G(s[:, None] + s[None, :])
          - 2.0 * G(s[:, None]) * G(s[None, :])) * p[:, None] * p[None, :]
    lG = 2.0 * np.concatenate(([0.0], np.cumsum(G(np.arange(1, N) / w)))) / (2 * n + 1.0)
    w_n = 1.0 / (2.0 * n + 1.0)
    return tau, Rm, w_n, lG


def theta(s):
    return (np.log1p(s * s) + 2.0 * s * np.arctan(s)) / (2.0 * (1.0 + s * s))


def thetatil_mp(s):
    s = mp.mpf(s)
    return (mp.log(1 + s ** 2) + 2 * s * mp.atan(s)) / (2 * (1 + s ** 2)) - mp.atan(s) / s


def N_of(s):
    s = mp.mpf(s)
    return (-s ** 3 * mp.log(1 + s ** 2) + s ** 3 + 3 * s ** 2 * mp.atan(s)
            - s + mp.atan(s))


def psi(s):
    s = mp.mpf(s)
    return (s ** 3 * (mp.log(1 + s ** 2) - 1) - 3 * mp.pi / 2 * s ** 2
            + 4 * s - 1 / s - mp.pi / 2)


if __name__ == "__main__":
    print("== (A) is there a spurious 2 in Phi_s?  test: int Phi == Theta ==")

    def Phi_mine(s, t):
        q = 1 + s * s + t * t
        return (s * s * t * (5 + s * s + t * t)
                / ((q * q - 4 * s * s * t * t) * (1 + s * s) * (1 + t * t)))

    def Phi_peer(s, t):
        return Phi_mine(s, t) / 2

    def R_repo(s, t):
        return (2 * s * s * t * t * (5 + s * s + t * t)
                / ((1 + (s - t) ** 2) * (1 + (s + t) ** 2) * (1 + s * s) * (1 + t * t)))

    for s in (0.6, 1.4, 2.86906499489):
        Im = float(mp.quad(lambda x: Phi_mine(mp.mpf(s), x), [0, 1, 10, mp.inf]))
        Ip = float(mp.quad(lambda x: Phi_peer(mp.mpf(s), x), [0, 1, 10, mp.inf]))
        t0 = 1.3
        print(f"   s={s:<14g} int Phi_mine={Im:.12f}  int Phi_peer={Ip:.12f}"
              f"  Theta={theta(s):.12f}")
        print(f"{'':>18} R_repo/(2t) at t={t0}: {R_repo(s, t0) / (2 * t0):.12f}"
              f"   Phi_mine={Phi_mine(s, t0):.12f}")

    print("\n== (B) which D_nn closes the identity, with vs without the L1 diagonal ==")
    print(f"{'varpi':>7} {'N':>5} | {'with L1: D=w':>13} {'with L1: D=w+lG':>16}"
          f" | {'no L1: D=w':>11} {'no L1: D=w+lG':>14}")
    for w in (0.30, 1.00, 3.40, 20.0):
        N = 600
        tau, Rm, w_n, lG = pieces(w, N)
        base = 2.0 * np.outer(tau, tau) + Rm
        for tag, H in (("", build_H(w, N, True)),):
            r1 = np.abs(H - (base - np.diag(w_n))).max()
            r2 = np.abs(H - (base - np.diag(w_n + lG))).max()
        H0 = build_H(w, N, False)
        r3 = np.abs(H0 - (base - np.diag(w_n))).max()
        r4 = np.abs(H0 - (base - np.diag(w_n + lG))).max()
        print(f"{w:7.2f} {N:5d} | {r1:13.3e} {r2:16.3e} | {r3:11.3e} {r4:14.3e}")
    print("   (residual = max|H - (2*tau tau^T + R - D)|)")

    print("\n== (C) analytic finitization: N(s) <= -psi(s), and psi>0 for s>=3 ==")
    print(f"{'s':>8} {'d/ds Thetatil':>14} {'N/(s^2(1+s^2)^2)':>18} {'N(s)':>12}"
          f" {'-psi(s)':>12} {'N<=-psi':>8} {'psi>0':>6}")
    ok = True
    for s in (1.0, 2.0, 2.869065, 3.0, 3.5, 5.0, 10.0, 100.0, 1e4):
        d = mp.diff(thetatil_mp, mp.mpf(s))
        Nv, pv_ = N_of(s), psi(s)
        pred = Nv / (mp.mpf(s) ** 2 * (1 + mp.mpf(s) ** 2) ** 2)
        a, b = Nv <= -pv_, pv_ > 0
        ok &= bool(a)
        print(f"{s:8.4g} {float(d):14.6e} {float(pred):18.6e} {float(Nv):12.4f}"
              f" {float(-pv_):12.4f} {str(bool(a)):>8} {str(bool(b)):>6}")

    print("\n   psi and psi' on a fine grid s in [3, 1e6]:")
    ss = np.geomspace(3.0, 1e6, 400)
    psiv = np.array([float(psi(s)) for s in ss])
    dpsi = np.array([float(mp.diff(psi, mp.mpf(float(s)))) for s in ss])
    print(f"      min psi  = {psiv.min():.6g}  at s={ss[psiv.argmin()]:.4f}"
          f"   (psi(3)={float(psi(3)):.6f})")
    print(f"      min psi' = {dpsi.min():.6g}  at s={ss[dpsi.argmin()]:.4f}")
    print(f"      => N(s) < 0 for all s >= 3 :"
          f" {'PASS' if psiv.min() > 0 and dpsi.min() > 0 else 'FAIL'}")

    print("\n   direct check N(s) < 0 on s in [2.87, 1e6] (2000 pts):")
    s2 = np.geomspace(2.8691, 1e6, 2000)
    Nv = np.array([float(N_of(float(s))) for s in s2])
    print(f"      max N = {Nv.max():.6g} at s={s2[Nv.argmax()]:.6f}"
          f"   ({'all negative' if Nv.max() < 0 else 'SIGN PROBLEM'})")
