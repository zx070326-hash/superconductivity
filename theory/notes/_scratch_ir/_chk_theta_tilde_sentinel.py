"""[numerical-evidence] Direction sentinel + stationary-point uniqueness for the
Thetatil chain.

Claim under test (the chain that now carries the 12-cell criterion):

  (1) lam2(H) <= lam1(R - D)                             [repo 1.F.3(c), analytic]
  (2)         <= Stil(w) := sup_n [sum_m w_m Rk(s_n,s_m) - D_nn]   [Lemma Lambda]
  (3)         -> Thetatil* = sup_s [Theta(s) - arctan(s)/s] + discretization pad

Sentinel A: Stil(w) >= lam2(H_N) at every (w, N).  A single failure means the
            derivation is wrong -> abort.
Sentinel B: Stil(w) >= lam1(R_N - D_N) too (the tighter intermediate).
Check  C  : is Thetatil' sign-changing exactly once on (0,inf)?  (findroot only
            returns *a* root; the sup claim needs uniqueness.)
Check  D  : how much of the excess  Stil(w) - Thetatil*  is covered by the
            analytic midpoint pad on D_nn alone,
                pad_D = (D^2/8) * V2[G] / (2 s_n),  V2[G] = int|G''| = 3*sqrt(3)/2
            (the remainder must be covered by the row-sum pad, repo M2 -- which I
             have NOT independently re-verified; that is logged as debt).

Run: OMP_NUM_THREADS=1 python _chk_theta_tilde_sentinel.py
"""
import os, io, contextlib, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("tt", os.path.join(HERE, "_chk_theta_tilde.py"))
tt = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(tt)

spec2 = importlib.util.spec_from_file_location("rp", os.path.join(HERE, "_chk_Rop_identity_peel.py"))
rp = importlib.util.module_from_spec(spec2)
with contextlib.redirect_stdout(io.StringIO()):
    spec2.loader.exec_module(rp)

THETATIL_STAR = 0.0737378188046
V2G = 3.0 * np.sqrt(3.0) / 2.0          # int |G''| dt, exact


def lam2_and_l1RD(w, N):
    H = rp.O(lambda k: w * w / (k * k + w * w), N)
    tau, Rm, D = rp.pieces(w, N)
    ev = np.linalg.eigvalsh(H)
    l1RD = np.linalg.eigvalsh(Rm - np.diag(D))[-1]
    return ev[-2], l1RD, ev[-1]


if __name__ == "__main__":
    mp.mp.dps = 40

    print("== Check C: uniqueness of the Thetatil stationary point ==")
    ss = np.geomspace(1e-3, 1e5, 4001)
    dv = np.array([float(mp.diff(tt.theta_tilde_mp, s)) for s in ss])
    sign_changes = np.where(np.sign(dv[:-1]) != np.sign(dv[1:]))[0]
    print(f"   sign changes of Thetatil' on [1e-3,1e5], 4001 log pts:"
          f" {len(sign_changes)}")
    for i in sign_changes:
        print(f"     between s={ss[i]:.6f} and s={ss[i+1]:.6f}")
    print(f"   Thetatil'(1e-3) = {dv[0]:+.4e}  (rising)")
    print(f"   Thetatil'(1e5)  = {dv[-1]:+.4e}  (falling)")
    fine = np.linspace(2.0, 4.5, 2501)
    fv = np.array([float(tt.theta_tilde_mp(s)) for s in fine])
    print(f"   dense max on [2.0,4.5]: {fv.max():.12f} at s={fine[fv.argmax()]:.6f}")
    print(f"   findroot value        : {THETATIL_STAR:.12f}")

    print("\n== Sentinels A & B: Stil(w) vs lam2(H_N) and lam1(R_N - D_N) ==")
    print(f"{'varpi':>7} {'N':>6} {'lam2(H)':>12} {'lam1(R-D)':>12} {'Stil':>10}"
          f" {'A':>3} {'B':>3}")
    WS = [0.20, 0.30, 0.42, 0.60, 0.85, 1.00, 1.20, 1.70, 2.00, 2.40,
          3.40, 4.80, 6.80, 9.60, 13.60, 20.0, 30.0, 45.0]
    nfail = 0
    STIL = {}
    for w in WS:
        STIL[w] = tt.discrete_sup(w)[0]
    for w in WS:
        St = STIL[w]
        for N in (800, 1600, 3200):
            l2, l1rd, h = lam2_and_l1RD(w, N)
            a, b = St >= l2, St >= l1rd
            nfail += (not a) + (not b)
            if N == 1600 or not (a and b):
                print(f"{w:7.2f} {N:6d} {l2:12.4e} {l1rd:12.4e} {St:10.6f}"
                      f" {'ok' if a else 'FAIL':>3} {'ok' if b else 'FAIL':>3}")
    print(f"\n   total sentinel failures: {nfail}   "
          f"({'PASS' if nfail == 0 else 'ABORT - derivation is wrong'})")

    print("\n== Check D: is the excess covered by pad_D alone? ==")
    print(f"{'varpi':>7} {'Stil':>10} {'excess':>10} {'pad_D':>10} {'covered?':>9}")
    for w in WS:
        St, narg = tt.discrete_sup(w)
        s_arg = (narg + 0.5) / w
        pad_D = (1.0 / (w * w) / 8.0) * V2G / (2.0 * s_arg)
        exc = St - THETATIL_STAR
        print(f"{w:7.2f} {St:10.6f} {exc:+10.6f} {pad_D:10.6f}"
              f" {'yes' if exc <= pad_D else 'no':>9}")
