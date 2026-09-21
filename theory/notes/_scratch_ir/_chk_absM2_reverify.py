"""[numerical-evidence] W1 re-verification of MY OWN correction to (G'-abs):
does the corrected statement (Z-dom premise) survive B=2..5 and kernels other
than Einstein?  A correction can introduce a NEW false claim."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
G2 = 1.318140212132
rng = np.random.default_rng(31337)

def build(LZ, LD, fv, N):
    B = len(LZ); n = np.arange(N); u = 1/np.sqrt(2*n+1.)
    T = fv[np.abs(n[:,None]-n[None,:])].copy(); np.fill_diagonal(T, 0.0)
    P = (T + fv[n[:,None]+n[None,:]+1]) * np.outer(u,u)
    d = 2.0*np.concatenate(([0.], np.cumsum(fv[1:N])))/(2*n+1.)
    rsZ = LZ.sum(axis=1); M = np.zeros((B*N, B*N))
    for i in range(B):
        for j in range(B):
            blk = LD[i,j]*P
            if i == j: blk = blk - np.diag(rsZ[i]*d) - np.diag(u*u*(rsZ[i]-LD[i,i]))
            else:      blk = blk + np.diag(u*u*LD[i,j])
            M[i*N:(i+1)*N, j*N:(j+1)*N] = blk
    return M

N = 80; kk = np.arange(0, 2*N+2, dtype=float)
KERNELS = {}
for w in (0.15, 0.8, 4.0):
    KERNELS[f"einstein w={w}"] = (w**2/(kk**2+w**2), w**2)          # fv, M2/lam per unit
# Debye-like / power-law / two-Einstein mixtures: M2 = int 2w a2F = sum lam_p * wp^2
for w in (0.3, 2.0):
    fv = 0.5*(w**2/(kk**2+w**2)) + 0.5*((3*w)**2/(kk**2+(3*w)**2))
    KERNELS[f"two-mode w={w}"] = (fv, 0.5*w**2 + 0.5*(3*w)**2)
for a in (1.2, 2.5):                                   # f(k) = (1+ (k/a)^2)^-2  (super-Einstein tail)
    fv = 1.0/(1.0+(kk/a)**2)**2
    KERNELS[f"tail^-4 a={a}"] = (fv, None)

print(f"{'kernel':>18} {'B':>2} {'draws':>6} {'viol':>5} {'worst':>9}")
for tag, (fv, m2unit) in KERNELS.items():
    if m2unit is None:
        # M2 = 2 * sum_k f(k) * (2 pi T)^2 ... recover via the scalar identity:
        # lemma M gives lam1(O[f]) <= g2 * M2/kappa with M2/kappa = 2*sum_{k>=1} k*f(k)?  
        # Avoid guessing: calibrate m2unit from the SCALAR case so the test is
        # a pure premise test, not a moment-convention test.
        one = np.array([[1.0]])
        m2unit = np.linalg.eigvalsh(build(one, one, fv, N))[-1]/G2*1.0
        cal = " (calibrated to scalar, ratio<=1 by construction at B=1)"
    else:
        cal = ""
    for B in (2, 3, 5):
        worst, nv, D = -np.inf, 0, 400
        for _ in range(D):
            Z = rng.uniform(0, 2, (B,B)); Z = (Z+Z.T)/2
            Dm = np.clip(Z*rng.uniform(-1.5,1.5,(B,B)), -Z, Z); Dm = (Dm+Dm.T)/2
            Dm = np.clip(Dm, -Z, Z)
            lam1 = np.linalg.eigvalsh(build(Z, Dm, fv, N))[-1]
            rho = np.abs(np.linalg.eigvals(np.abs(Dm))).max()
            r = lam1/(G2*rho*m2unit)
            nv += r > 1.0; worst = max(worst, r)
        print(f"{tag:>18} {B:2d} {D:6d} {nv:5d} {worst:9.4f}{cal if B==2 else ''}")

# control: same kernels, premise DELIBERATELY broken
print("\ncontrol (premise broken: |LD| allowed up to 3x LZ):")
for tag, (fv, m2unit) in list(KERNELS.items())[:3]:
    worst, nv = -np.inf, 0
    for _ in range(400):
        Z = rng.uniform(0,2,(3,3)); Z=(Z+Z.T)/2
        Dm = Z*rng.uniform(-3,3,(3,3)); Dm=(Dm+Dm.T)/2
        lam1 = np.linalg.eigvalsh(build(Z, Dm, fv, N))[-1]
        rho = np.abs(np.linalg.eigvals(np.abs(Dm))).max()
        r = lam1/(G2*rho*m2unit); nv += r>1.0; worst=max(worst,r)
    print(f"{tag:>18}  {nv:4d}/400 violate, worst {worst:9.3f}")
