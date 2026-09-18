"""
Exploratory (non-certificate) probe: does h(w) = lambda_1(O[f_w]) ~ ln(w) + c0
survive when the Lorentz kernel f(k,w) = w^2/(k^2+w^2) is replaced by a
non-Lorentz family?  Kernel-generic operator O[f] is VERBATIM the KAY
construction (L1 diagonal cumulative sum, L2/L3 off-diagonal), only the
scalar function f(k,w) is swapped.  float64 point evaluation only ->
numerical-evidence level, not certificate.

Not committed; scratch per CLAUDE.md.
"""
import numpy as np
import time

def O_matrix(N, f):
    """Generic kernel-generic operator O[f]: f is a vectorized function
    f(k_array, w) -> values, applied to positive integers k>=1.
    Structure copied from p4_interchange_uniformity.H_mat, f swapped in
    for _brk_vec."""
    n = np.arange(N)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)
    A = np.abs(n[:, None] - n[None, :])   # |n-m|, 0 on diagonal (unused there)
    S = n[:, None] + n[None, :] + 1        # n+m+1 >= 1
    off = f(np.where(A == 0, 1, A))        # placeholder value on diag, zeroed next
    np.fill_diagonal(off, 0.0)
    M = (off + f(S)) * nrm
    if N > 1:
        cum = np.concatenate(([0.0], np.cumsum(2.0 * f(np.arange(1, N)))))
    else:
        cum = np.array([0.0])
    M[np.diag_indices(N)] -= cum / (2.0 * n + 1.0)
    return M


def h_of(w, fkernel, N):
    f = lambda k: fkernel(k.astype(np.float64), w)
    M = O_matrix(N, f)
    return float(np.linalg.eigvalsh(M)[-1])


KERNELS = {
    "lorentz (gamma=2, baseline)": lambda k, w: w**2 / (k**2 + w**2),
    "yukawa gapped alpha=2":       lambda k, w: w**2 / ((k + 2.0)**2 + w**2),
    "yukawa gapped alpha=5":       lambda k, w: w**2 / ((k + 5.0)**2 + w**2),
    "power-law gamma=1.5":         lambda k, w: w**1.5 / (k**1.5 + w**1.5),
    "power-law gamma=3":           lambda k, w: w**3.0 / (k**3.0 + w**3.0),
    "power-law gamma=4":           lambda k, w: w**4.0 / (k**4.0 + w**4.0),
    "gaussian cutoff":             lambda k, w: np.exp(-(k / w) ** 2),
    "exponential cutoff":          lambda k, w: np.exp(-k / w),
}

W_VALUES = [5.0, 15.0, 45.0, 135.0]

if __name__ == "__main__":
    print(f"{'kernel':30s} " + " ".join(f"h(w={w:g})".rjust(12) for w in W_VALUES) + "   slope(last two)  c0=h-ln(w) at w=135")
    results = {}
    t0 = time.time()
    for name, fk in KERNELS.items():
        hs = []
        for w in W_VALUES:
            N = min(3000, int(np.ceil(50 * w)))
            hs.append(h_of(w, fk, N))
        results[name] = hs
        # convergence check at largest w: double N
        w_big = W_VALUES[-1]
        N1 = min(3000, int(np.ceil(50 * w_big)))
        N2 = min(4000, int(np.ceil(90 * w_big)))
        h_ref1 = h_of(w_big, fk, N1)
        h_ref2 = h_of(w_big, fk, N2)
        conv = abs(h_ref2 - h_ref1)
        slope = (hs[-1] - hs[-2]) / (np.log(W_VALUES[-1]) - np.log(W_VALUES[-2]))
        c0 = hs[-1] - np.log(W_VALUES[-1])
        print(f"{name:30s} " + " ".join(f"{h:12.5f}" for h in hs)
              + f"   {slope:14.5f}  {c0:8.5f}   (N-conv |dh|={conv:.2e})")
    print(f"\nelapsed {time.time()-t0:.1f}s")

# --- extension: push to larger w for slope-convergence check on the kernels
# with the slowest transient (yukawa gap, slow power-law tail, exponential) ---
if __name__ == "__main__":
    print("\n--- extended w range (adds w=400) for slope convergence ---")
    extra = {
        "lorentz (gamma=2, baseline)": KERNELS["lorentz (gamma=2, baseline)"],
        "yukawa gapped alpha=5": KERNELS["yukawa gapped alpha=5"],
        "power-law gamma=1.5": KERNELS["power-law gamma=1.5"],
        "exponential cutoff": KERNELS["exponential cutoff"],
    }
    W2 = [135.0, 400.0]
    for name, fk in extra.items():
        hs = []
        for w in W2:
            N = min(4500, int(np.ceil(35 * w)))
            hs.append(h_of(w, fk, N))
        slope = (hs[-1] - hs[-2]) / (np.log(W2[-1]) - np.log(W2[-2]))
        print(f"{name:30s} h(135)={hs[0]:.5f} h(400)={hs[1]:.5f}  slope={slope:.5f}")
