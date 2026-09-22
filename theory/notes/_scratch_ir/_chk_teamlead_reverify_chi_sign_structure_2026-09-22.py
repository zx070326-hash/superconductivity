"""
team-lead independent re-verification (2026-09-22 evening) of the two structural
claims in misc/essenberger_chi_sign_structure_2026-09-22.md, section 11.1.

Written from scratch (does NOT reuse the third-append-block scripts) so that the
two claims have two independent implementations.  Both claims were ALREADY
established in that note's third append block; this is confirmation, not a new
finding -- see section 11 of the note.

(a) The section-2.3 "coherent" orbital projection
        chi0[i,j] += conj(U_k[i,m]) * chi_band_mn * U_kq[j,n]
    is GAUGE DEPENDENT: the U_k phase belongs to band m and the U_kq phase to
    band n, neither paired with its own conjugate.  Re-phasing each eigenvector
    column (a pure gauge transformation for a Hermitian matrix) changes the
    eigenvalues, and even the NUMBER of negative ones.
    Stronger gauge group than the note's +/-1 enumeration: random U(1) per band.

(b) The gauge-invariant 4-U-insertion form (chi0_dynamic in
    _chk_gauge_invariance_gate_dynamic_2026-09-22.py) gives a static Re chi0
    that is negative DEFINITE, and the RPA-dressed chi^s stays negative for
    every U, with min eig(I - U chi0) >= 1 and INCREASING in U -- i.e. raising
    the Hubbard U moves AWAY from the Stoner instability, which is the overall
    sign-flip fingerprint.

Run:  OMP_NUM_THREADS=1 python _chk_teamlead_reverify_chi_sign_structure_2026-09-22.py
"""
import os
import importlib.util
import numpy as np
from scipy.linalg import eigh

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def _load(fname):
    spec = importlib.util.spec_from_file_location(fname, os.path.join(SCRIPT_DIR, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gate = _load("_chk_gauge_invariance_gate_dynamic_2026-09-22.py")
ham = _load("_graser1003_hamiltonian_2026-09-22.py")
H_func = ham.build_H_graser
MU = float(np.load(os.path.join(SCRIPT_DIR, "_graser_mu_value.npy")))

N = 5
T = 0.02
ETA = 0.02
KG = 16


def fermi(e, mu, T):
    return 1.0 / (1.0 + np.exp(np.clip((e - mu) / T, -500, 500)))


def chi0_coherent(q, seed=None):
    """The note's section-2.3 projection.  seed=None -> whatever gauge eigh returns."""
    rng = np.random.default_rng(seed) if seed is not None else None
    qx, qy = q
    out = np.zeros((N, N), dtype=complex)
    ks = np.linspace(-np.pi, np.pi, KG)
    dk = (2 * np.pi / KG) ** 2 / (2 * np.pi) ** 2
    for kx in ks:
        for ky in ks:
            ek, Uk = eigh(H_func(kx, ky))
            ekq, Ukq = eigh(H_func(kx + qx, ky + qy))
            if rng is not None:                      # pure gauge: one phase per band
                Uk = Uk * np.exp(1j * rng.uniform(0, 2 * np.pi, size=N))[None, :]
                Ukq = Ukq * np.exp(1j * rng.uniform(0, 2 * np.pi, size=N))[None, :]
            w = (fermi(ek, MU, T)[:, None] - fermi(ekq, MU, T)[None, :]) \
                / (ek[:, None] - ekq[None, :] + 1j * ETA)
            out += np.einsum('im,mn,jn->ij', np.conj(Uk), w, Ukq) * dk
    return out


def part_a():
    print("=== (a) section-2.3 coherent projection is gauge dependent ===")
    q = (np.pi, 0.0)
    base = chi0_coherent(q)
    print("  q=(pi,0)  eig(Re sym), no rephasing :",
          np.round(np.linalg.eigvalsh((base.real + base.real.T) / 2), 4))
    for s in (1, 2, 3):
        c = chi0_coherent(q, seed=s)
        print(f"  q=(pi,0)  eig(Re sym), gauge seed={s}:",
              np.round(np.linalg.eigvalsh((c.real + c.real.T) / 2), 4))

    def nneg(c):
        return int((np.linalg.eigvalsh((c + c.conj().T) / 2) < 0).sum())

    rng = np.random.default_rng(5)
    qs = [(np.pi, 0.0), (0.0, 0.0), (np.pi, np.pi)] + \
         [(rng.uniform(-np.pi, np.pi), rng.uniform(-np.pi, np.pi)) for _ in range(5)]
    print("  negative-eigenvalue count of the Hermitian part:")
    print("    q                    orig  |  6 random gauges")
    for q in qs:
        row = [nneg(chi0_coherent(q))] + [nneg(chi0_coherent(q, seed=100 + s)) for s in range(6)]
        print(f"    ({q[0]:+.3f},{q[1]:+.3f})      {row[0]}    |  {row[1:]}")
    print("  -> the count is a function of the phase choice, not of the material.\n")


def part_b():
    print("=== (b) gauge-invariant 4-U form: static Re chi0 is negative definite ===")
    rng = np.random.default_rng(2026)
    mx, mn = -1e9, 1e9
    for _ in range(24):
        q = (rng.uniform(-np.pi, np.pi), rng.uniform(-np.pi, np.pi))
        c = gate.chi0_dynamic(q, 0.0, H_func, MU, T=T, eta=ETA, k_grid=KG, norb=N, rephase=False)
        ev = np.linalg.eigvalsh((c.real + c.real.T) / 2)
        mx, mn = max(mx, ev.max()), min(mn, ev.min())
    print(f"  over 24 random q:  max eig = {mx:+.6e}   min eig = {mn:+.6e}")
    print("  -> max eigenvalue strictly negative => Re chi0 negative DEFINITE\n")

    print("=== (b') RPA dressing keeps it negative; Stoner denominator recedes ===")
    I = np.eye(N)
    rng = np.random.default_rng(99)
    for u in (0.3, 0.8, 1.5, 3.0):
        mx, mn, st = -1e9, 1e9, 1e9
        r2 = np.random.default_rng(99)
        for _ in range(12):
            q = (r2.uniform(-np.pi, np.pi), r2.uniform(-np.pi, np.pi))
            c = gate.chi0_dynamic(q, 0.0, H_func, MU, T=T, eta=ETA, k_grid=KG, norb=N, rephase=False)
            X = (c.real + c.real.T) / 2
            st = min(st, np.linalg.eigvalsh(I - u * X).min())
            R = X @ np.linalg.inv(I - u * X)
            ev = np.linalg.eigvalsh((R + R.T) / 2)
            mx, mn = max(mx, ev.max()), min(mn, ev.min())
        print(f"  U={u:4.1f}: chi^s eig in [{mn:+.4e}, {mx:+.4e}]   min eig(I-U*chi0) = {st:.4f}")
    _ = rng
    print("  analytic: chi0<=0, u>0 => lam/(1-u*lam)<=0 and eig(I-u*chi0)>=1.")
    print("  min eig(I-U*chi0) INCREASES with U => moving away from Stoner => overall sign flipped.")


if __name__ == "__main__":
    part_a()
    part_b()
