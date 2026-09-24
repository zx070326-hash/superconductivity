import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv

KB = 8.617333262e-2          # meV/K
HB = 64.65415134             # meV per sqrt(eV/Å^2/amu)
MUS = 0.13
BASE = os.path.join(r"D:\physics\superconductivity", "theory", "data", "external",
                     "mishra_h3s_d3s", "Data-H3S-D3S-PH-EPW-Anharmonic")

# ---------------------------------------------------------------------
# Machinery reused verbatim from the team lead's
# _chk_main_realistic_shape_300K_requirement_2026-09-23.py.
# ---------------------------------------------------------------------

def load_a2f(path):
    rows = []
    for ln in open(path, errors="replace"):
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        try:
            v = [float(x) for x in s.split()]
        except ValueError:
            continue
        if len(v) >= 2:
            rows.append(v[:2])
    a = np.array(rows)
    w, f = a[:, 0], np.clip(a[:, 1], 0.0, None)
    m = w > 1e-9
    return w[m], f[m]


quad = (lambda y, x: np.trapezoid(y, x)) if hasattr(np, "trapezoid") else (lambda y, x: np.trapz(y, x))


def moments(w, a):
    lam = 2.0 * quad(a / w, w)
    return lam, np.sqrt(2.0 * quad(a * w, w) / lam)


def G_of_spectrum(kmax, T, w, a):
    x2 = (2.0 * np.pi * T * np.arange(kmax, dtype=float)) ** 2
    den = w[None, :] ** 2 + x2[:, None]
    return quad((2.0 * w * a)[None, :] / den, w)


def top_eig_G(T, G, N, mus, wc):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    L1 = np.zeros(N)
    L1[1:] = 2.0 * np.cumsum(G[1:N]) / (2.0 * n[1:] + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    M = (G[np.abs(i - j)] + G[i + j + 1]) * u[i] * u[j] - np.diag(L1)
    M[n, n] -= G[0] * u * u
    if mus > 0.0:
        Nc = int(np.floor(wc / (2.0 * np.pi * T)))
        if Nc > 0:
            uc = np.where(n < Nc, u, 0.0)
            M -= 2.0 * mus * np.outer(uc, uc)
    return eigh(M, eigvals_only=True, subset_by_index=[N - 1, N - 1])[0]


def w2_for_Tc(T0_K, lam, shape, mus, r, N=400):
    T = T0_K * KB
    if shape is None:                                   # Einstein: wE = wbar2
        def ev(w2):
            k = np.arange(2 * N + 2, dtype=float)
            G = lam * w2 ** 2 / (w2 ** 2 + (2.0 * np.pi * T * k) ** 2)
            return top_eig_G(T, G, N, mus, r * w2) - 1.0
    else:
        w0, a0, lam0, w20 = shape
        def ev(w2):
            s, c = w2 / w20, lam / lam0
            G = G_of_spectrum(2 * N + 2, T, s * w0, c * a0)
            return top_eig_G(T, G, N, mus, r * w2) - 1.0
    lo, hi = 1.0, 30000.0
    if ev(hi) < 0:
        return np.nan
    if ev(lo) > 0:
        return lo
    for _ in range(46):
        mid = np.sqrt(lo * hi)
        if ev(mid) < 0:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


def Tc_K(lam, wE, mus, r, N=400):
    # pseudo-root guard (team-lead audit 2026-09-23): once N_c=floor(wc/(2*pi*T)) > N,
    # the Coulomb term saturates and lambda_1(T)=1 can pick up a spurious extra root
    # below the true Tc; the bisection start point then decides which root you land
    # on, independent of N. Floor T so N_c never exceeds N inside this bisection.
    wc = r * wE
    f = lambda T: top_eig_G(T, lam * wE ** 2 / (wE ** 2 + (2.0 * np.pi * T * np.arange(2 * N + 2)) ** 2), N, mus, wc) - 1.0
    lo, hi = max(1e-3, wc / (2.0 * np.pi * (N + 1))), 900.0
    if not (f(lo) > 0 > f(hi)):
        return float("nan")
    for _ in range(44):
        mid = np.sqrt(lo * hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi) / KB


def eta_over_M(lam, wE_meV):
    return lam * wE_meV ** 2 / HB ** 2


# ---------------------------------------------------------------------
# self-check
# ---------------------------------------------------------------------
anchor = Tc_K(1.721, 149.26, MUS, 10.0)
print(f"[self-check] Tc_K(1.721, 149.26, 0.13, 10.0) = {anchor:.4f} K (expect ~246.8 K)")
assert abs(anchor - 246.8) / 246.8 < 0.01

# ---------------------------------------------------------------------
# load all three real alpha^2F shapes
# ---------------------------------------------------------------------
SHAPES = {}
for tag, rel in [("H3S anharm", "H3S-200GPa/EPW_anharmonic/H3S.a2f"),
                  ("H3S harm", "H3S-200GPa/EPW_harmonic/H3S.a2f"),
                  ("D3S anharm", "D3S-200GPa/EPW_anharmonic/D3S.a2f")]:
    p = os.path.join(BASE, rel.replace("/", os.sep))
    w0, a0 = load_a2f(p)
    l0, v0 = moments(w0, a0)
    SHAPES[tag] = (w0, a0, l0, v0)
    print(f"[shape] {tag}: lambda0={l0:.4f}, wbar2_0={v0:.2f} meV")

# ---------------------------------------------------------------------
# CONSISTENCY CHECK (not an independent implementation — same solver
# (w2_for_Tc/G_of_spectrum), same .a2f data as the team lead's own
# _chk_main_realistic_shape_300K_requirement_2026-09-23.py; re-running
# it here only confirms the call is wired correctly, it has zero
# degrees of freedom against the team lead's own numbers).
# ---------------------------------------------------------------------
LAM_TL = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]
TL_TABLE_300K = {   # team lead's numbers, quoted for the printed re-run comparison only
    1.5: (15.682, 25.436, 20.540, 22.068),
    2.0: (11.981, 17.331, 14.757, 15.627),
    2.5: (10.373, 13.911, 12.298, 12.872),
    3.0: (9.507, 12.119, 10.953, 11.405),
    3.5: (8.968, 10.981, 10.121, 10.469),
    4.0: (8.595, 10.244, 9.538, 9.819),
    5.0: (8.103, 9.304, 8.812, 9.011),
}
print("\n[consistency check — SAME solver/data as team lead, NOT an independent implementation]:")
for lam in LAM_TL:
    wE = w2_for_Tc(300.0, lam, None, MUS, 10.0)
    e_ein = eta_over_M(lam, wE)
    wA = w2_for_Tc(300.0, lam, SHAPES["H3S anharm"], MUS, 10.0)
    e_hA = eta_over_M(lam, wA)
    wH = w2_for_Tc(300.0, lam, SHAPES["H3S harm"], MUS, 10.0)
    e_hH = eta_over_M(lam, wH)
    wD = w2_for_Tc(300.0, lam, SHAPES["D3S anharm"], MUS, 10.0)
    e_hD = eta_over_M(lam, wD)
    tl = TL_TABLE_300K[lam]
    print(f"  lam={lam:.1f}  mine=({e_ein:.3f},{e_hA:.3f},{e_hH:.3f},{e_hD:.3f})  "
          f"team-lead=({tl[0]:.3f},{tl[1]:.3f},{tl[2]:.3f},{tl[3]:.3f})")

# ---------------------------------------------------------------------
# Full-range band (lambda = 1.2..6.0): all three real shapes now drawn as
# complete curves, not just the 7 team-lead points, per team-lead request.
# ---------------------------------------------------------------------
lam_band = np.array([1.2, 1.5, 1.72, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0])
curve_ein, curve_h3s_anh, curve_h3s_har, curve_d3s_anh = [], [], [], []
print("\n[300K full-range band] lambda=1.2..6.0:")
for lam in lam_band:
    wE = w2_for_Tc(300.0, lam, None, MUS, 10.0)
    e_ein = eta_over_M(lam, wE)
    wA = w2_for_Tc(300.0, lam, SHAPES["H3S anharm"], MUS, 10.0)
    e_hA = eta_over_M(lam, wA)
    wH = w2_for_Tc(300.0, lam, SHAPES["H3S harm"], MUS, 10.0)
    e_hH = eta_over_M(lam, wH)
    wD = w2_for_Tc(300.0, lam, SHAPES["D3S anharm"], MUS, 10.0)
    e_hD = eta_over_M(lam, wD)
    curve_ein.append(e_ein); curve_h3s_anh.append(e_hA); curve_h3s_har.append(e_hH); curve_d3s_anh.append(e_hD)
    print(f"  lam={lam:.2f}  Einstein={e_ein:.3f}  H3S-anh={e_hA:.3f}  H3S-har={e_hH:.3f}  D3S-anh={e_hD:.3f}")

curve_ein = np.array(curve_ein); curve_h3s_anh = np.array(curve_h3s_anh)
curve_h3s_har = np.array(curve_h3s_har); curve_d3s_anh = np.array(curve_d3s_anh)

# ---------------------------------------------------------------------
# Einstein-family full grid (background dashed contours, unchanged)
# ---------------------------------------------------------------------
lam_grid = np.unique(np.round(np.concatenate([np.linspace(0.8, 2.0, 13), np.linspace(2.2, 6.0, 15)]), 6))
wE_grid = np.geomspace(20.0, 400.0, 30)
rows = []
for lam in lam_grid:
    for wE in wE_grid:
        tc = Tc_K(lam, wE, MUS, 10.0, 400)
        rows.append((lam, wE, eta_over_M(lam, wE), tc))
n_nan = sum(1 for r_ in rows if np.isnan(r_[3]))
csv_path = r"D:\physics\superconductivity\theory\notes\_scratch_ir\_agent_feasibility_grid_2026-09-23.csv"
with open(csv_path, "w", newline="") as fh:
    wtr = csv.writer(fh)
    wtr.writerow(["lambda", "omega_E_meV", "eta_over_M", "Tc_K"])
    for r_ in rows:
        wtr.writerow([f"{r_[0]:.6f}", f"{r_[1]:.6f}", f"{r_[2]:.6f}", "nan" if np.isnan(r_[3]) else f"{r_[3]:.6f}"])
print(f"\n[grid] Einstein: {len(rows)} cells, {n_nan} nan")

# ---------------------------------------------------------------------
# materials: ONLY the two sourced anchors (Nb / MgB2 / Pb dropped —
# team lead confirmed those (lambda, wbar2) values had no provenance
# in-repo; my own "textbook recall" substitutes were equally unsourced)
# ---------------------------------------------------------------------
materials = [
    ("H3S-200GPa (anharm)", 1.721, 149.26, 203),
    ("H3S-200GPa (harm)", 2.293, 129.31, 203),
    ("D3S-200GPa (anharm)", 1.784, 105.43, 150),
]

# ---------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 7.5))

LAM = np.array([[lam for wE in wE_grid] for lam in lam_grid])
WE = np.array([[wE for wE in wE_grid] for lam in lam_grid])
ETA = LAM * WE ** 2 / HB ** 2
lookup = {(round(r_[0], 6), round(r_[1], 6)): r_[3] for r_ in rows}
TC = np.array([[lookup[(round(lam, 6), round(wE, 6))] for wE in wE_grid] for lam in lam_grid])

levels_all = [100, 200, 203, 250, 350]   # 300 dropped: duplicates the bold Einstein curve below (same
                                          # solver, same r/mu*/N -> necessarily coincide, zero d.o.f.)
cs = ax.contour(LAM, ETA, TC, levels=levels_all, colors="#5b7fae", linestyles="dashed", linewidths=1.0, alpha=0.55)
ax.clabel(cs, fmt=lambda v: f"{v:g} K", fontsize=7, inline=True)

# the 300K band: 4 curves (full lambda=1.2..6.0 range) + fill
ax.plot(lam_band, curve_ein, color="#5b7fae", lw=2.4, ls="dashed", marker="s", ms=4,
        label="300K, Einstein — NECESSARY, conditional on unproven $\\mu^*{>}0$ extremality (excl.-only)")
ax.plot(lam_band, curve_h3s_har, color="#c98a2a", lw=1.4, marker="^", ms=3.5,
        label="300K, H3S harmonic shape — SUFFICIENT, no premise (constructive)")
ax.plot(lam_band, curve_d3s_anh, color="#7a3ca0", lw=1.4, marker="v", ms=3.5,
        label="300K, D3S anharmonic shape — SUFFICIENT, no premise (constructive)")
ax.plot(lam_band, curve_h3s_anh, color="#8a1f1f", lw=2.4, marker="o", ms=4,
        label="300K, H3S anharmonic shape — SUFFICIENT, no premise (constructive)")
ax.fill_between(lam_band, curve_ein, curve_h3s_anh, color="#f4c04a", alpha=0.30, zorder=0,
                 label="300K band: Einstein (lower, conditional) to H3S-anharm (upper, sufficient) — interior UNDETERMINED")
_ui = 2  # lambda~1.72: band is widest here (Einstein vs H3S-anharm most separated), well clear of eta/M=9.172 line
_uy = curve_ein[_ui] * (curve_h3s_anh[_ui] / curve_ein[_ui]) ** 0.75  # biased toward the upper curve, away from the horizontal line
ax.text(lam_band[_ui], _uy, "UNDETERMINED",
        fontsize=8.5, color="#8a6d1a", ha="center", va="center", style="italic",
        rotation=0, alpha=0.9,
        bbox=dict(boxstyle="round", fc="#fff6df", ec="none", alpha=0.7, pad=0.2))

# room-T headline: drop-lines from the eta/M=9.172 crossings down to the x-axis.
# Solved directly (bisection on lambda, not interpolated off the 12-point curve) so the
# printed numbers match the textbox/note's 3.30 / 5.20 to 2 decimals.
def _lam_for_eta(target_eta, shape):
    def h(lm):
        w2 = w2_for_Tc(300.0, lm, shape, MUS, 10.0)
        return eta_over_M(lm, w2) - target_eta
    lo, hi = 1.2, 6.0
    if h(lo) * h(hi) > 0:
        return np.nan
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if h(lo) * h(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)

lam_lo = _lam_for_eta(9.172, None)
lam_hi = _lam_for_eta(9.172, SHAPES["H3S anharm"])
for lx, lbl in ((lam_lo, f"{lam_lo:.2f}"), (lam_hi, f"{lam_hi:.2f}")):
    ax.plot([lx, lx], [9.172, 0.11], color="0.4", lw=0.8, ls=(0, (1, 1)), zorder=4)
    ax.text(lx, 0.115, lbl, fontsize=7.5, color="0.25", ha="center", va="bottom")
ax.annotate("room-T band", (np.sqrt(lam_lo * lam_hi), 0.13), fontsize=7.5, color="0.25",
            ha="center", va="bottom", style="italic")

ax.set_yscale("log")
ax.set_xlabel(r"$\lambda$")
ax.set_ylabel(r"$\eta/M$  [eV/Å$^2$/amu]  ($=\lambda\bar\omega_2^2/\hbar^2$)")
ax.set_title(
    "Feasibility map: isotropic Eliashberg, $\\mu^*=0.13$, $\\omega_c=10\\bar\\omega_2$, N=400, finite $\\lambda$ grid\n"
    "(1) dashed = Einstein = NECESSARY, conditional on UNPROVEN $\\mu^*{>}0$ extremality, excl.-only.\n"
    "(2) solid = real $\\alpha^2F$ shapes = SUFFICIENT, constructive, no premise.\n"
    "(3) severity: numerical-evidence; isotropic Eliashberg only; no lattice-stability constraint.",
    fontsize=8.3, loc="left")

offsets = {"H3S-200GPa (anharm)": (-72, 26), "H3S-200GPa (harm)": (10, 22), "D3S-200GPa (anharm)": (14, -34)}
h3s_pts = [m for m in materials if m[0].startswith("H3S")]
for name, lam, wE, tc_exp in materials:
    eta = eta_over_M(lam, wE)
    ax.scatter([lam], [eta], color="black", zorder=6, s=32, marker="D")
    label = name if name.startswith("H3S") else f"{name}\nTc_expt={tc_exp}K"
    ax.annotate(label, (lam, eta), fontsize=6.8,
                textcoords="offset points", xytext=offsets[name],
                arrowprops=(dict(arrowstyle="-", color="0.55", lw=0.5, shrinkA=0, shrinkB=2)
                            if name.startswith("D3S") else None))
# the two H3S points are the SAME material/experiment under two ab initio phonon
# treatments (anharmonic vs harmonic), not two independently-matching materials —
# connect them and give a single shared Tc_expt label instead of two copies.
h3s_x = [m[1] for m in h3s_pts]; h3s_y = [eta_over_M(m[1], m[2]) for m in h3s_pts]
ax.plot(h3s_x, h3s_y, color="black", lw=0.8, ls=":", alpha=0.6, zorder=5)
ax.annotate("shared Tc_expt=203K\n(same H3S sample, anharm vs harm phonons)",
            (np.mean(h3s_x), np.sqrt(h3s_y[0] * h3s_y[1])), fontsize=6.6,
            textcoords="offset points", xytext=(135, 30),
            arrowprops=dict(arrowstyle="-", color="0.5", lw=0.6, shrinkA=0, shrinkB=2))

ax.axhline(9.172, color="gray", ls=":", lw=1)
ax.annotate(r"$\eta/M=9.172$ (H$_3$S anharm.; fixed design" + "\n" + r"constraint, NOT a physical invariant)",
            (4.6, 9.172), fontsize=7.0, color="gray", ha="left", va="top",
            textcoords="offset points", xytext=(0, -30),
            arrowprops=dict(arrowstyle="-", color="0.6", lw=0.6, shrinkA=0, shrinkB=2))

txt = (r"$\lambda$ needed for 300 K at fixed $\eta/M{=}9.172$ (r=6/10/20):" + "\n"
       r"  Einstein (conditional, lower bound): 3.26 / 3.30 / 3.12" + "\n"
       r"  H3S anharm (real, sufficient):        4.96 / 5.20 / 5.00" + "\n"
       r"  H3S harm (real, sufficient):          4.24 / 4.44 / 4.24" + "\n"
       r"  D3S anharm (real, sufficient):        4.54 / 4.75 / 4.54" + "\n"
       r"Shape penalty is NOT constant: H3S-anharm/Einstein ratio (by $\eta/M$)" + "\n"
       r"  runs 1.148 ($\lambda$=5) $\to$ 1.622 ($\lambda$=1.5) — no single $\times$1.24 factor." + "\n"
       r"H3S harm/anharm $\eta/M$ agreement (=9.172) is NOT a physical check: it is" + "\n"
       r"  circular (frozen e-ph matrix elements, only freq. remapped) — only the" + "\n"
       r"  k=+1 moment is conserved (1e-6), neighbours k=-1,0,+2,+3 differ 20-54%." + "\n"
       r"D3S/H3S first-moment ratio IS content-bearing: 0.51712 vs. frozen-matrix-" + "\n"
       r"  element ideal $M_H/M_D$=0.5000 (residual 3.4%); at 60meV cut, high band" + "\n"
       r"  (H-dominated) gives 0.504/conserved, low band (S-mode, ~2% weight) not." + "\n"
       r"Anchors: only H3S/D3S shown (Nb/MgB2/Pb dropped, no in-repo source)." + "\n"
       r"Re-ran the team lead's solver on the same data (consistency check," + "\n"
       r"  not an independent implementation — zero degrees of freedom)."
       )
ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=6.6, va="bottom", ha="left",
        bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9))

ax.legend(loc="upper right", fontsize=7, framealpha=0.9)
fig.tight_layout()
png_path = r"D:\physics\superconductivity\theory\notes\_scratch_ir\_agent_feasibility_map_2026-09-23.png"
fig.savefig(png_path, dpi=150)
print(f"\n[plot] saved to {png_path}")
