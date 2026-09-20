"""Phi_E(lam, mu*, r) MAIN TABLE + N_c dependence + sawtooth + AD comparison + lam->inf.

Phi_E(lam,mu*,r) := Tc_E(lam,mu*,r)/Tc_E(lam,0),  r := omega_c/omega_E,
Nc = floor(r*wbar) is the REALIZED Coulomb window index at the root.

Sections
  1  Phi_E table on lam x mu* x r                      [numerical-evidence]
  2  N_c dependence at FROZEN N_c (smooth branch) + dlnPhi/dlnN_c   [num-ev]
  3  Anderson-Morel collapse test: does Tc depend on (mu*,Nc) only
     through  mu_AM := 1/(1/mu* + ln Nc)?                           [num-ev]
  4  Sawtooth: the staircase Nc=floor(r wbar) makes Phi_E(r) NON-C1;
     jump amplitudes at integer-Nc boundaries; integer vs period-average [num-ev]
  5  AD/McMillan  exp[-1.04(1+lam)/(lam-mu*(1+0.62 lam))]  ratio vs Phi_E [num-ev]
  6  lam -> infinity verdict: does Phi_E -> 1 (mu* "diluted")?      [num-ev + analytic]

LEVEL: numerical-evidence throughout (float64, N=220 truncation) unless a line
says otherwise.  Certified grid points are produced by the separate script
_chk_PhiE_cert_2026-09-21.py.
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib.util as _ilu  # noqa: E402


def _load(stem):
    """Load a sibling scratch module whose filename contains hyphens."""
    spec = _ilu.spec_from_file_location(
        stem.replace("-", "_"), os.path.join(HERE, stem + ".py"))
    mod = _ilu.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


_E = _load("_phiE_engine_2026-09-21")
C_INF, TWO_PI = _E.C_INF, _E.TWO_PI
g_fixedNc, nc_of, tc_over_wE, wstar_sup = (
    _E.g_fixedNc, _E.nc_of, _E.tc_over_wE, _E.wstar_sup)


LAMS = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 5.0]
MUS = [0.05, 0.10, 0.13, 0.15]
RS = [3.0, 5.0, 10.0, 20.0]


# ---------------------------------------------------------------- fixed-Nc root
def wstar_fixedNc(lam, mu_star, Nc, N=220, w_lo=0.02, w_hi=40.0, sub=200):
    """Smallest root of g(wbar; Nc FROZEN) = 1/lam.  No staircase: smooth branch."""
    mu_t = 0.0 if mu_star == 0.0 else mu_star / lam
    gstar = 1.0 / lam
    ws = np.geomspace(w_lo, w_hi, sub)
    prev = None
    for w in ws:
        g = g_fixedNc(w, mu_t, Nc, N)
        if prev is not None and prev[1] < gstar <= g:
            a, b = prev[0], w
            for _ in range(60):
                m = 0.5 * (a + b)
                if g_fixedNc(m, mu_t, Nc, N) < gstar:
                    a = m
                else:
                    b = m
            return 0.5 * (a + b)
        prev = (w, g)
    return None


def tc_fixedNc(lam, mu_star, Nc, N=220):
    w = wstar_fixedNc(lam, mu_star, Nc, N)
    return None if w is None else 1.0 / (TWO_PI * w)


# ------------------------------------------------------------------ AD/McMillan
def ad_ratio(lam, mu_star):
    """Allen-Dynes/McMillan exponent ratio Tc(mu*)/Tc(0) (f1 f2 and prefactor cancel:
    they carry no mu* dependence in the AD formula)."""
    e0 = -1.04 * (1 + lam) / lam
    den = lam - mu_star * (1 + 0.62 * lam)
    if den <= 0:
        return 0.0
    e1 = -1.04 * (1 + lam) / den
    return float(np.exp(e1 - e0))


def sec(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


def main():
    # ---------------------------------------------------------------- SECTION 1
    sec("1  Phi_E(lam, mu*, r) TABLE   [numerical-evidence, float64 N=220]")
    tc0 = {}
    for lam in LAMS:
        t, w, _ = tc_over_wE(lam, 0.0, 1.0)
        tc0[lam] = t
    print("  mu*=0 reference (= f_E(lam) * C_inf * sqrt(lam)):")
    print("   lam    Tc/wE     f_E(lam)")
    for lam in LAMS:
        print(f"  {lam:4.1f}  {tc0[lam]:.6f}   {tc0[lam]/(C_INF*np.sqrt(lam)):.4f}")
    print()
    table = {}
    for r in RS:
        print(f"  --- r = omega_c/omega_E = {r:g} ---")
        print("   lam  " + "".join(f"  mu*={m:<5.2f}(Nc)" for m in MUS))
        for lam in LAMS:
            row = f"  {lam:4.1f} "
            for mus in MUS:
                t, w, nc = tc_over_wE(lam, mus, r)
                if t is None:
                    row += "      ---      "
                    table[(lam, mus, r)] = (None, None)
                    continue
                phi = t / tc0[lam]
                table[(lam, mus, r)] = (phi, nc)
                row += f"  {phi:.4f}({nc:2d})"
            print(row)
        print()

    # ---------------------------------------------------------------- SECTION 2
    sec("2  N_c DEPENDENCE at FROZEN N_c (smooth branch)  [numerical-evidence]")
    print("Phi at frozen Nc: Tc(lam,mu*,Nc)/Tc(lam,0).  This removes the staircase")
    print("so d lnPhi / d lnNc is well defined (the realized-Nc curve is NOT C1).")
    print()
    for lam in (1.0, 2.0, 3.0, 5.0):
        for mus in (0.10, 0.13):
            ncs = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]
            vals = []
            for nc in ncs:
                t = tc_fixedNc(lam, mus, nc)
                vals.append(None if t is None else t / tc0[lam])
            print(f"  lam={lam:4.1f} mu*={mus:.2f}:")
            print("    Nc  " + "".join(f"{n:>8d}" for n in ncs))
            print("    Phi " + "".join(
                ("     ---" if v is None else f"{v:8.4f}") for v in vals))
            # local log-log slopes
            sl = []
            for i in range(1, len(ncs)):
                if vals[i] is None or vals[i - 1] is None:
                    sl.append(None)
                    continue
                sl.append((np.log(vals[i]) - np.log(vals[i - 1])) /
                          (np.log(ncs[i]) - np.log(ncs[i - 1])))
            print("    dlnPhi/dlnNc " + "".join(
                ("    ---" if s is None else f"{s:+7.4f}") for s in sl))
            print()

    # ---------------------------------------------------------------- SECTION 3
    sec("3  ANDERSON-MOREL COLLAPSE: does Tc depend on (mu*,Nc) only via "
        "mu_AM=1/(1/mu*+lnNc)?")
    print("If the Morel-Anderson log is the WHOLE Nc dependence, then Phi at frozen")
    print("Nc should be a function of (lam, mu_AM) alone.  Test: for each lam, pick")
    print("pairs (mu*,Nc) and (mu*',Nc') with equal mu_AM and compare Phi.")
    print()
    for lam in (1.0, 2.0, 3.0):
        print(f"  lam={lam}:")
        print("    mu*     Nc    mu_AM      Phi     |  partner check")
        recs = []
        for mus in (0.05, 0.08, 0.10, 0.13, 0.15, 0.20, 0.30):
            for nc in (2, 4, 8, 16, 32, 64):
                t = tc_fixedNc(lam, mus, nc)
                if t is None:
                    continue
                mam = 1.0 / (1.0 / mus + np.log(nc))
                recs.append((mus, nc, mam, t / tc0[lam]))
        recs.sort(key=lambda z: z[2])
        # bin by mu_AM and report within-bin spread of Phi
        import collections
        bins = collections.defaultdict(list)
        for mus, nc, mam, phi in recs:
            bins[round(np.log(mam) * 8)].append((mus, nc, mam, phi))
        worst = 0.0
        for key in sorted(bins):
            grp = bins[key]
            if len(grp) < 2:
                continue
            phis = [g[3] for g in grp]
            spread = (max(phis) - min(phis)) / np.mean(phis)
            worst = max(worst, spread)
            memb = " ".join(f"({g[0]:.2f},{g[1]:d})" for g in grp)
            print(f"    mu_AM~{np.mean([g[2] for g in grp]):.4f}  Phi in "
                  f"[{min(phis):.4f},{max(phis):.4f}]  spread {spread*100:5.2f}%  {memb}")
        print(f"    => worst within-bin Phi spread at lam={lam}: {worst*100:.2f}%")
        print()

    # ---------------------------------------------------------------- SECTION 4
    sec("4  SAWTOOTH / NON-C1: Phi_E(r) across integer-Nc boundaries")
    print("At the root, Nc=floor(r*wbar*).  As r sweeps, Nc steps by 1 and Phi_E")
    print("JUMPS.  Below: Phi_E on a fine r-grid, jump amplitudes at the steps.")
    print()
    for lam, mus in ((2.0, 0.13), (3.0, 0.10)):
        rs = np.linspace(2.0, 24.0, 441)
        out = []
        for r in rs:
            t, w, nc = tc_over_wE(lam, mus, float(r))
            out.append((r, None if t is None else t / tc0[lam], nc))
        print(f"  lam={lam} mu*={mus}:")
        jumps = []
        for i in range(1, len(out)):
            if out[i][2] != out[i - 1][2] and out[i][1] and out[i - 1][1]:
                jumps.append((out[i][0], out[i - 1][2], out[i][2],
                              out[i][1] - out[i - 1][1]))
        phivals = [o[1] for o in out if o[1] is not None]
        print(f"    Phi_E range over r in [2,24]: [{min(phivals):.4f},{max(phivals):.4f}]")
        print(f"    Nc steps detected: {len(jumps)}")
        amps = [abs(j[3]) for j in jumps]
        if amps:
            print(f"    |jump| in Phi_E: max {max(amps):.4f}  median "
                  f"{np.median(amps):.4f}  (relative max {max(amps)/np.mean(phivals)*100:.2f}%)")
            for j in jumps[:8]:
                print(f"      r={j[0]:6.3f}  Nc {j[1]:2d}->{j[2]:2d}  dPhi={j[3]:+.5f}")
        print()
    print("  READING: the jumps are the r-image of the same discontinuity that makes")
    print("  Tc(T) non-C1 (repo memory: 'coulomb ramp makes Tc non-C1'): the")
    print("  discriminant is the distance of r*wbar* to an integer.  A table indexed")
    print("  by INTEGER Nc is single-valued and exact for the model; a table indexed")
    print("  by r must either quote the realized Nc or report a period average.")

    # ---------------------------------------------------------------- SECTION 5
    sec("5  Phi_E vs the ALLEN-DYNES/McMILLAN mu* dependence  [numerical-evidence]")
    print("AD ratio := exp[-1.04(1+lam)/(lam-mu*(1+0.62lam))] / exp[-1.04(1+lam)/lam]")
    print("(prefactor, f1, f2 carry no mu* in the AD formula, so they cancel exactly).")
    print()
    for mus in (0.10, 0.13):
        print(f"  mu* = {mus}")
        print("   lam   AD-ratio |  Phi_E r=3   r=5    r=10   r=20  |  Phi_E(r=10)/AD")
        for lam in LAMS:
            adr = ad_ratio(lam, mus)
            cells = []
            for r in RS:
                v = table.get((lam, mus, r), (None, None))[0]
                cells.append("  ---  " if v is None else f"{v:.4f} ")
            v10 = table.get((lam, mus, 10.0), (None, None))[0]
            rat = "  ---" if v10 is None else f"{v10/adr:.3f}"
            print(f"  {lam:4.1f}   {adr:.4f}  | " + " ".join(cells) + f" |  {rat}")
        print()

    # ---------------------------------------------------------------- SECTION 6
    sec("6  lam -> INFINITY: is mu* DILUTED by strong coupling?")
    print("ANALYTIC (level: analytic, given the repo's Cor(ii) strong-coupling law):")
    print("  h(wbar) ~ g(2) wbar^2  =>  wbar*(lam) ~ 1/sqrt(lam g(2)) -> 0,")
    print("  Tc/omega_E -> C_inf sqrt(lam).   But Nc = floor(r*wbar*), so once")
    print("  wbar*(lam) < 1/r  the Coulomb window is EMPTY (Nc=0) and K_c = 0")
    print("  IDENTICALLY.  Hence Phi_E(lam,mu*,r) = 1 EXACTLY for all")
    print("  lam > lam_c(r), where lam_c(r) solves  r*wbar*(lam) = 1.")
    print("  Leading estimate: wbar* ~ 1/(2 pi C_inf sqrt(lam) f_E(lam)) =>")
    print("  lam_c(r) ~ [r/(2 pi C_inf)]^2 = (r/1.14809)^2  (times f_E^-2 correction).")
    print()
    for r in (3.0, 5.0, 10.0, 20.0):
        est = (r / (TWO_PI * C_INF)) ** 2
        # locate lam_c numerically: smallest lam with realized Nc == 0
        lo, hi = 1.0, 4000.0
        for _ in range(60):
            m = 0.5 * (lo + hi)
            _, w, nc = tc_over_wE(m, 0.13, r, w_max=40.0)
            if nc is None:
                lo = m
                continue
            if nc == 0:
                hi = m
            else:
                lo = m
        print(f"  r={r:5.1f}: naive estimate lam_c ~ {est:8.2f};  measured lam_c "
              f"(mu*=0.13) = {0.5*(lo+hi):8.2f}")
    print()
    print("  Approach to 1 BELOW lam_c (the physically reachable part):")
    for r in (3.0, 10.0):
        print(f"    r={r:g}, mu*=0.13:")
        print("      lam   " + "".join(f"{x:>9.1f}" for x in
                                       (1, 2, 3, 5, 8, 12, 20, 40, 80)))
        vals = []
        for lam in (1, 2, 3, 5, 8, 12, 20, 40, 80):
            t, w, nc = tc_over_wE(float(lam), 0.13, r, w_max=40.0)
            t0, _, _ = tc_over_wE(float(lam), 0.0, 1.0, w_max=40.0)
            vals.append((None if t is None else t / t0, nc))
        print("      Phi   " + "".join(
            ("      ---" if v[0] is None else f"{v[0]:9.4f}") for v in vals))
        print("      Nc    " + "".join(
            ("      ---" if v[1] is None else f"{v[1]:9d}") for v in vals))
    print()
    print("  CONTRAST with Allen-Dynes: AD ratio -> exp[-1.04*0.62*mu*/(1-0.62mu*)]")
    for mus in (0.10, 0.13, 0.15):
        lim = np.exp(-1.04 * 0.62 * mus / (1 - 0.62 * mus))
        print(f"    mu*={mus}: AD limit = {lim:.4f}  (AD says mu* is NEVER fully "
              f"diluted; Phi_E says it IS, abruptly, at lam_c(r))")


if __name__ == "__main__":
    main()
