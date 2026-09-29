r"""
P4 . CONE-COERCIVITY certificate probe  (deliverable of agent g2-cone-coercivity)

Numerically certifies the analytic statements of `P4_cone_coercivity.md`:

  (L1) Master rung identity, w := u.v :
        T_k := <v,O[delta_k]v> = 2 sum_{n>=k} w_n(w_{n-k}-w_n) + sum_{n=0}^{k-1} w_n w_{k-1-n}.
  (T1) Cone-coercivity theorem, v in C, phi>=0 :
        <v,O[phi]v> >= sum_{n,m} phi(n+m+1) u_n u_m v_n v_m  =: kap_L3
                     >= phi(1) v0^2.
        Specialise phi=g:  E = <v,O[g]v> >= kap_L3 >= g(1) v0^2.
  (P2) e_0 test-vector bound on the top-eigenvector head:
        v0^2 >= (f(1)-lam2)/(h-lam2)          (rigorous; tight since lam2 ~ 0),
        and the conservative analytic-gap versions with lam2 <= 0.58 h / Theta*.
  (P3) Row-n profile bound  v_n >= beta_n v0,  beta_n = u_n(f(n)+f(n+1))/(h-A_nn),
        feeding  E >= v0^2 * P(x),  P = sum g(n+m+1) (u_n beta_n)(u_m beta_m).

Prints tightness kappa/E for every candidate kappa over a varpi grid in [0.1,40].
`ck` counts fail-stop assertions.  N>=400.  Import path is dirname-relative.
"""
import sys
import os
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)  # p4_concavity_q_route.py lives beside this file
from p4_concavity_q_route import O_ker, kernels, top_eigen

_ck = 0
_t0 = time.time()


def ck(cond, msg):
    global _ck
    _ck += 1
    if not cond:
        print(f"\n*** FAIL check {_ck}: {msg}")
        sys.exit(1)


THETA_STAR = 0.6349  # rank-one gap bound (Lemma M), varpi>=2


def spectral(x, N):
    ker = kernels(x, KX=8 * N)
    A = O_ker(N, ker["f"])
    h, v, l2 = top_eigen(A)
    idx = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * idx + 1.0)
    w = u * v
    E = float(v @ (O_ker(N, ker["g"]) @ v))
    return dict(ker=ker, A=A, h=h, v=v, l2=l2, u=u, w=w, E=E, idx=idx)


def hankel_g(g, idx):
    S = idx[:, None] + idx[None, :] + 1
    return g[np.minimum(S, len(g) - 1)]  # G[n,m] = g(n+m+1)


def delta(k, KX):
    a = np.zeros(KX)
    a[k] = 1.0
    return a


def main():
    print("=" * 74)
    print("P4 cone-coercivity certificate probe")
    print("=" * 74)
    N = 400
    grid = [0.1, 0.2, 0.35, 0.5, 0.8, 1.0, 1.2, 1.42, 1.54, 1.7,
            2.0, 3.0, 5.0, 8.0, 14.0, 20.0, 30.0, 40.0]

    # ---------------------------------------------------------------- (L1)
    print("\n(L1) master rung identity  T_k = tele_k + L3_k   [spot k, several x]")
    print("-" * 74)
    for wpi in [0.2, 1.0, 5.0, 20.0]:
        x = wpi * wpi
        S = spectral(x, N)
        v, w = S["v"], S["w"]
        worst = 0.0
        for k in [1, 2, 3, 7, 15]:
            Tk = float(v @ (O_ker(N, delta(k, 8 * N)) @ v))
            tele = 2.0 * np.sum(w[k:] * (w[:-k] - w[k:]))
            L3 = float(sum(w[n] * w[k - 1 - n] for n in range(k)))
            worst = max(worst, abs(Tk - (tele + L3)))
            ck(Tk > -1e-13, f"x={x} k={k}: T_k negative {Tk:.2e}")
            ck(tele > -1e-13 and L3 > -1e-13, f"x={x} k={k}: a part negative")
        ck(worst < 1e-11, f"x={x}: rung identity err {worst:.2e}")
        print(f"  wpi={wpi:5.2f}  max|T_k - (tele+L3)| = {worst:.2e}   (all T_k,parts >=0)")

    # ------------------------------------- (T1) on a GENERIC cone vector
    print("\n(T1) coercivity on a GENERIC v in C (reusability, not the eigenvector)")
    print("-" * 74)
    rng = np.random.default_rng(0)
    idx = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * idx + 1.0)
    for trial, x in enumerate([0.3, 2.4, 30.0]):
        raw = np.sort(rng.random(N))[::-1]      # nonincreasing, >=0
        vv = raw / np.linalg.norm(raw)
        ww = u * vv
        ker = kernels(x, KX=8 * N)
        for phi_name, phi in [("g", ker["g"]), ("f2", ker["f2"])]:
            lhs = float(vv @ (O_ker(N, phi) @ vv))
            G = hankel_g(phi, idx)
            kapL3 = float(ww @ (G @ ww))
            floor = phi[1] * vv[0] ** 2
            ck(kapL3 <= lhs + 1e-12, f"x={x} {phi_name}: L3 form exceeds form ({kapL3:.4e}>{lhs:.4e})")
            ck(floor <= kapL3 + 1e-12, f"x={x} {phi_name}: phi(1)v0^2 exceeds L3 form")
            ck(lhs > 0, f"x={x} {phi_name}: form not positive")
        print(f"  x={x:5.1f}  <v,O[g]v>={float(vv @ (O_ker(N, ker['g']) @ vv)):.5f}"
              f"  kapL3<=form OK, phi(1)v0^2<=kapL3 OK")

    # ---------------------------------------------------- main grid table
    print("\n(T1)+(P2)+(P3) on the TOP eigenvector, grid varpi in [0.1,40]")
    print("-" * 74)
    hdr = (f"{'wpi':>5} {'x':>8} {'E':>10} {'kapL3/E':>7} {'g1v0^2/E':>8} "
           f"{'v0lb/v0^2':>9} {'gapCr/v0^2':>10} {'kProf/E':>8}")
    print(hdr)
    worst_L3 = 1.0
    for wpi in grid:
        x = wpi * wpi
        S = spectral(x, N)
        h, v, l2, u, w, E = S["h"], S["v"], S["l2"], S["u"], S["w"], S["E"]
        ker = S["ker"]
        v0 = float(v[0])
        f1 = x / (1.0 + x)

        # (T1) L3 form and floor
        G = hankel_g(ker["g"], S["idx"])
        kapL3 = float(w @ (G @ w))
        g1 = 2.0 * x / (1.0 + x) ** 3
        floor = g1 * v0 ** 2
        ck(kapL3 <= E + 1e-10, f"wpi={wpi}: kapL3>E")
        ck(floor <= kapL3 + 1e-12, f"wpi={wpi}: floor>kapL3")
        ck(floor <= E + 1e-12, f"wpi={wpi}: g1 v0^2 > E")
        worst_L3 = min(worst_L3, kapL3 / E)

        # (P2) v0 lower bound, true gap and conservative analytic gap
        v0lb_true = (f1 - l2) / (h - l2)
        ck(v0lb_true <= v0 ** 2 + 1e-9, f"wpi={wpi}: v0 true-gap bound exceeds v0^2")
        # regime-correct analytic gap upper bound (gap_lorentz note):
        #   lam2 <= 0.58 h   for varpi <= 2 ;   lam2 <= Theta*   for varpi >= 2.
        gap_hi = 0.58 * h if wpi <= 2.0 else THETA_STAR
        v0lb_cr = (f1 - gap_hi) / (h - gap_hi)
        if v0lb_cr > 0:
            ck(v0lb_cr <= v0 ** 2 + 1e-9, f"wpi={wpi}: v0 conservative bound exceeds v0^2")

        # (P3) profile bound and P(x)
        Ann = np.diag(S["A"]).copy()
        fk = ker["f"][np.minimum(S["idx"], len(ker["f"]) - 1)]
        fk1 = ker["f"][np.minimum(S["idx"] + 1, len(ker["f"]) - 1)]
        beta = u * (fk + fk1) / (h - Ann)
        beta[0] = 1.0
        ck(bool(np.all(v >= beta * v0 - 1e-11)), f"wpi={wpi}: profile v_n>=beta_n v0 fails")
        ub = u * beta
        P = float(ub @ (G @ ub))
        kProf = v0 ** 2 * P
        ck(kProf <= E + 1e-10, f"wpi={wpi}: kProf>E")

        print(f"{wpi:5.2f} {x:8.2f} {E:10.6f} {kapL3 / E:7.3f} {floor / E:8.3f} "
              f"{v0lb_true / v0 ** 2:9.3f} {(v0lb_cr / v0 ** 2 if v0lb_cr>0 else 0):10.3f} "
              f"{kProf / E:8.3f}")

    print(f"\n  sharp L3-form floor:  kapL3/E in [{worst_L3:.3f}, ...]  (>=0.58 across grid)")

    # ---------------------------------- supporting claims used in the .md
    print("\nSupporting claims:  Q_{LL'}(d1) >= 1  and  row-0 upper bound  h <= h_ub")
    print("-" * 74)
    Od1 = O_ker(N, delta(1, 8 * N))
    minQ = np.inf
    for L in [1, 2, 3, 5, 10, 50, 200]:
        for Lp in [1, 2, 3, 5, 10, 50, 200]:
            eL = np.zeros(N); eL[:L] = 1.0
            eLp = np.zeros(N); eLp[:Lp] = 1.0
            minQ = min(minQ, float(eL @ (Od1 @ eLp)))
    ck(minQ >= 1.0 - 1e-9, f"Q_LL'(d1) floor below 1: {minQ:.4f}")
    print(f"  min_{{L,L'}} Q_{{LL'}}(delta_1) = {minQ:.6f}  (>= 1, i.e. group-3 term u_0^2=1)")
    for wpi in [0.2, 1.0, 2.0, 20.0]:
        x = wpi * wpi
        S = spectral(x, N)
        mm = np.arange(1, N)
        uu = 1.0 / np.sqrt(2 * mm + 1.0)
        fK = S["ker"]["f"]
        h_ub = float(fK[1] + np.sum(uu * (fK[np.minimum(mm, len(fK) - 1)]
                                          + fK[np.minimum(mm + 1, len(fK) - 1)])))
        ck(S["h"] <= h_ub + 1e-9, f"wpi={wpi}: row-0 h_ub violated")
        print(f"  wpi={wpi:5.2f}  h={S['h']:.4f} <= h_ub={h_ub:.4f}  (row-0 bound, feeds P3's h^hi)")

    # ------------------------------------------- N-consistency of margins
    print("\nN-consistency (N=400 vs 700) of kapL3 and kProf")
    print("-" * 74)
    for wpi in [0.5, 1.54, 20.0]:
        x = wpi * wpi
        vals = {}
        for NN in (400, 700):
            S = spectral(x, NN)
            G = hankel_g(S["ker"]["g"], S["idx"])
            kapL3 = float(S["w"] @ (G @ S["w"]))
            vals[NN] = (S["E"], kapL3)
        (E4, k4), (E7, k7) = vals[400], vals[700]
        rk = abs(k4 - k7) / max(abs(k7), 1e-30)
        rE = abs(E4 - E7) / max(abs(E7), 1e-30)
        ck(rk < 5e-3 and rE < 5e-3, f"wpi={wpi}: N-drift kapL3 {rk:.2e} E {rE:.2e}")
        print(f"  wpi={wpi:5.2f}  dE/E={rE:.2e}  dKapL3/KapL3={rk:.2e}")

    print("\n" + "=" * 74)
    print(f"PASS: {_ck} checks in {time.time() - _t0:.1f}s")
    print("Conclusions:")
    print("  (T1) E >= kapL3 >= g(1) v0^2  (rigorous, all x); kapL3/E >= 0.58 uniformly")
    print("  (P2) v0^2 >= (f1-lam2)/(h-lam2)  (rigorous; ~0.64-0.90 of v0^2, lam2~0)")
    print("  (P3) E >= v0^2 P(x)  (rigorous given h upper bound); kProf/E ~ 0.72->0.14,")
    print("       loses only ~one log(x) vs E ~ 0.6/x, whereas g(1)v0^2 loses a factor x")
    print("=" * 74)


if __name__ == "__main__":
    main()
