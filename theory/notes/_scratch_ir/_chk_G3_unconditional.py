"""[numerical-evidence] What does G3 (c=3) cost the ceiling on REAL multiband
matrices?  Degradation factor is 1+(c-1)s, s = off-diagonal Perron share."""
import numpy as np

def share(M):
    ev, V = np.linalg.eigh(M); rho = ev[-1]; w = np.abs(V[:, -1])
    if w.min() <= 1e-12: return rho, np.nan
    off = (M - np.diag(np.diag(M))) @ w / w / rho
    return rho, off.max()

CASES = {
  # MgB2, Golubov JPCM 14,1353 lambda_ij (convention B: lam_ij = N_j K_ij), N_s=0.300 N_p=0.410
  "MgB2 (B-symmetrised)": (np.array([[1.017, 0.213],[0.155, 0.448]]), (0.300, 0.410)),
  "MgB2 (A, as printed)": (np.array([[1.017, 0.184],[0.184, 0.448]]), None),
  # LiFeAs 4-band fitted (Ummarino Table 2, partial -- indicative only)
  "LiFeAs-like strong off": (np.array([[0.90, 1.15],[1.15, 1.50]]), None),
  "equal-coupling worst":   (np.array([[1.0, 1.0],[1.0, 1.0]]), None),
  "fully off-diagonal":     (np.array([[0.0, 1.0],[1.0, 0.0]]), None),
}
print(f"{'case':>24} {'rho':>8} {'s':>7} {'1+2s':>7} {'sqrt':>7}")
for tag, (L, N) in CASES.items():
    if N is not None:                     # symmetrise convention (B): sqrt(Ni/Nj) lam_ij
        Ni = np.array(N); S = np.sqrt(np.outer(Ni, 1/Ni)); L = L * S
        L = (L + L.T) / 2                 # residual asymmetry is rounding
    rho, s = share(np.abs(L))
    print(f"{tag:>24} {rho:8.4f} {s:7.4f} {1+2*s:7.4f} {np.sqrt(1+2*s):7.4f}")
print("\nbound on the factor: s in [0,1] => 1+2s in [1,3], sqrt in [1, 1.7321]")

print("\n=== sharp version: optimise the Gershgorin test vector ===")
print("factor = rho(M2') / rho(M2),  M2' := M2 + 2*offdiag(M2)   (>= no w to choose)")
print(f"{'case':>24} {'rho(M2)':>9} {'rho(M2p)':>9} {'ratio':>7} {'sqrt':>7} {'vs 1+2s':>8}")
for tag, (L, N) in CASES.items():
    if N is not None:
        Ni = np.array(N); L = L * np.sqrt(np.outer(Ni, 1/Ni)); L = (L + L.T)/2
    A = np.abs(L); Ap = A + 2*(A - np.diag(np.diag(A)))
    r0 = np.linalg.eigvalsh(A)[-1]; r1 = np.linalg.eigvalsh(Ap)[-1]
    _, s = share(A)
    print(f"{tag:>24} {r0:9.4f} {r1:9.4f} {r1/r0:7.4f} {np.sqrt(r1/r0):7.4f}"
          f" {np.sqrt(1+2*s):8.4f}")
print("\nalways rho(M2') <= 3 rho(M2) (elementwise M2' <= 3 M2) => factor <= sqrt(3) = 1.7321")
