import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
rng = np.random.default_rng(20260918)

# ---- Claim 1: the Coulomb block matrix is PSD even with BAND-DEPENDENT cutoffs.
# Blocks  M_ij = mu_ij * u^(i) (u^(j))^T  with u^(i)_n = (2n+1)^{-1/2} 1_{n < Nc_i}
# Note delivery assumes a single common K_c (a clean Kronecker product). Test the weaker structure.
def blockmat(mu, Ncs, N):
    B = len(Ncs)
    n = np.arange(N); u = (2*n+1.0)**-0.5
    U = np.array([np.where(n < Nc, u, 0.0) for Nc in Ncs])      # B x N
    M = np.zeros((B*N, B*N))
    for i in range(B):
        for j in range(B):
            M[i*N:(i+1)*N, j*N:(j+1)*N] = mu[i, j]*np.outer(U[i], U[j])
    return M, U

worst_psd, worst_nonpsd = np.inf, np.inf
for trial in range(4000):
    B = int(rng.integers(2, 6)); N = 40
    Ncs = rng.integers(1, N, size=B)                 # band-dependent cutoffs, deliberately unequal
    A = rng.normal(size=(B, B)); mu = A @ A.T        # PSD by construction
    M, U = blockmat(mu, Ncs, N)
    worst_psd = min(worst_psd, np.linalg.eigvalsh(M)[0])
    # control: elementwise-nonneg but NOT PSD
    mu2 = np.abs(rng.normal(size=(B, B))); mu2 = (mu2+mu2.T)/2
    if np.linalg.eigvalsh(mu2)[0] < -1e-9:
        worst_nonpsd = min(worst_nonpsd, np.linalg.eigvalsh(blockmat(mu2, Ncs, N)[0])[0])
print("Claim 1  band-dependent cutoffs, mu PSD : min eig over 4000 trials =", worst_psd)
print("Claim 1  control (elementwise>=0, non-PSD): min eig            =", worst_nonpsd)

# identity behind it: <v, Mv> = a^T mu a  with a_i = <v_i, u^(i)>   -> check to machine precision
err = 0.0
for _ in range(300):
    B = int(rng.integers(2, 6)); N = 30
    Ncs = rng.integers(1, N, size=B)
    A = rng.normal(size=(B, B)); mu = A @ A.T
    M, U = blockmat(mu, Ncs, N)
    v = rng.normal(size=(B, N))
    a = np.einsum('in,in->i', v, U)
    err = max(err, abs(v.ravel() @ M @ v.ravel() - a @ mu @ a))
print("Claim 1  identity <v,Mv> = a^T mu a : max abs err =", err)

# ---- Claim 2: mu_ij = <rho_i, V rho_j> is a Gram matrix when V is a PSD kernel.
# screened Coulomb V(q)=4pi/(q^2+ks^2) > 0 pointwise in Fourier  => PSD kernel (Bochner)
nq = 64
q = np.fft.fftfreq(nq)*nq
Q = np.sqrt(q[:,None]**2 + q[None,:]**2)
worst = np.inf
for _ in range(500):
    ks = rng.uniform(0.3, 5.0)
    Vq = 4*np.pi/(Q**2 + ks**2)                 # >0 pointwise
    B = int(rng.integers(2, 6))
    rho = np.abs(rng.normal(size=(B, nq, nq)))  # Fermi-surface weights, nonneg
    # mu_ij = sum_k rhohat_i* Vq rhohat_j   (convolution in real space = product in Fourier)
    rh = np.fft.fft2(rho, axes=(1,2))
    mu = np.real(np.einsum('ikl,kl,jkl->ij', rh.conj(), Vq, rh))
    worst = min(worst, np.linalg.eigvalsh((mu+mu.T)/2)[0]/max(1.0, np.abs(mu).max()))
print("Claim 2  mu from screened-Coulomb Gram : min normalised eig over 500 =", worst)

# ---- Claim 3: Anderson-Morel transport mu* = mu (I + L mu)^{-1} preserves PSD, L>0 scalar
worst = np.inf
for _ in range(2000):
    B = int(rng.integers(2, 6)); L = rng.uniform(0.05, 8.0)
    A = rng.normal(size=(B,B)); mu = A @ A.T
    ms = mu @ np.linalg.inv(np.eye(B) + L*mu)
    worst = min(worst, np.linalg.eigvalsh((ms+ms.T)/2)[0]/max(1.0, np.abs(ms).max()))
    assert np.allclose(ms, ms.T, atol=1e-9), "mu* lost symmetry"
print("Claim 3  AM transport mu(I+L mu)^-1    : min normalised eig over 2000 =", worst)
