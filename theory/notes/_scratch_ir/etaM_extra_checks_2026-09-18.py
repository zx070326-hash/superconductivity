r"""Extra independent checks for P4_eta_over_M.md: sec 2.C table (G(lambda))
and sec 1.E (E1) isotope-drift eta ratios.  Re-derives the H(varpi) operator
from scratch (same as etaM_material_recheck_2026-09-18.py); does not import
p4_eta_over_M.py or p4_interchange_uniformity.py.
"""
import numpy as np
from scipy.optimize import brentq

G2_HI = 1.318140212373478
C_INF = np.sqrt(G2_HI) / (2 * np.pi)

def bracket(k, w):
    k = np.asarray(k, dtype=np.float64)
    return w * w / (k * k + w * w)

def H_operator(N, w):
    n = np.arange(N)
    inv_sqrt = 1.0 / np.sqrt(2.0 * n + 1.0)
    norm = np.outer(inv_sqrt, inv_sqrt)
    absdiff = np.abs(n[:, None] - n[None, :])
    summ = n[:, None] + n[None, :] + 1
    off = bracket(absdiff, w)
    np.fill_diagonal(off, 0.0)
    M = (off + bracket(summ, w)) * norm
    cum = np.concatenate(([0.0], np.cumsum(2.0 * bracket(np.arange(1, N), w))))
    M[np.diag_indices(N)] -= cum / (2.0 * n + 1.0)
    return M

def h_of(w, N=1200):
    if w <= 0:
        return 0.0
    return float(np.linalg.eigvalsh(H_operator(N, w))[-1])

def varpi_c(lam, N=1200):
    tgt = 1.0 / lam
    return brentq(lambda w: h_of(w, N) - tgt, 1e-7, 1e7, xtol=1e-13, rtol=8.9e-16)

def G_of_lambda(lam, N=1200):
    wc = varpi_c(lam, N)
    return 1.0 / (2 * np.pi * wc * np.sqrt(lam)), wc

print("=== sec 2.C table: lambda, varpi_c, G(lambda), G/C_inf ===")
for lam in (0.30, 0.40, 0.60, 1.00, 1.71, 2.06, 2.67, 3.62, 10, 30, 100):
    G, wc = G_of_lambda(lam)
    print(f"    lam={lam:7.2f}  varpi_c={wc:.6f}  G={G:.8f}  G/Cinf={G/C_INF:.6f}")

print("\n=== sec 1.E (E1): isotope eta-drift ratios ===")
# H3S/D3S: SSCHA anharmonic values from isotope_curated.md (Errea 2015 Table I)
lamH, wlogH = 1.84, 92.86
lamD, wlogD = 1.87, 73.3
s = 2.0  # full mass-doubling approximation (H->D), as used generically in the doc
etaD_over_etaH = s * (lamD/lamH) * (wlogD/wlogH)**2
print(f"    H3S(anh): lam H->D {lamH}->{lamD} (+{100*(lamD/lamH-1):.1f}%), "
      f"wlog H->D {wlogH}->{wlogD} (ratio {wlogD/wlogH:.4f})")
print(f"    eta_D/eta_H (s=2, using wlog as proxy for w2) = {etaD_over_etaH:.4f}  "
      f"(doc claims '~1.27')")

# PdH: harmonic vs anharmonic lambda listed in isotope_curated.md
# (this is a DIFFERENT ratio -- lambda_har vs lambda_anh at fixed isotope, NOT
#  an isotope substitution; kept separate from the H->D eta ratio in Stage 3)
lamH_anh, wlogH_anh, TcH = 0.40, 50.21, 5.0
lamD_anh, wlogD_anh, TcD = 0.46, 37.69, 6.5
etaD_over_etaH_pdh = 2.0 * (lamD_anh/lamH_anh) * (wlogD_anh/wlogH_anh)**2
print(f"    PdH/PdD (anh): eta_D/eta_H (s=2, w_log proxy) = {etaD_over_etaH_pdh:.4f}  "
      f"(doc/Stage3 claims 1.296, rounded '1.30')")
