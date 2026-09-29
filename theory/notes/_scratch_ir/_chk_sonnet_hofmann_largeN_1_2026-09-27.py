import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np

# Independent numerical solve of Hofmann-Chowdhury-Kivelson-Berg (arXiv:2105.09322) Eq. (31):
#
#   Delta(n) = [lam*w0^2/(2*pi*T)^2] * sum_{n'!=n over all integers}
#              [Delta(n') - Delta(n)*(2n'+1)/(2n+1)] / |2n'+1| / (n-n')^2
#
# with the particle-hole symmetry Delta(-n-1) = Delta(n) folded in (standard
# Eliashberg gap-equation symmetry; the equation is written for n>=0 after folding).
# Tc is the temperature at which this linear homogeneous system first has a
# nontrivial solution (dimensionless eigenvalue lam*w0^2/(2 pi Tc)^2 = c(lam)).
#
# We solve, for fixed dimensionless coupling lam and w0=1, for the largest T
# at which the matrix operator has eigenvalue 1, using a Matsubara cutoff N.
# Richardson-extrapolate in N to remove truncation bias (per CLAUDE.md discipline:
# multi-N agreement checks truncation, not rounding).

def build_matrix(lam, w0, T, Nmax):
    # indices n = 0..Nmax-1 (folded, positive branch only)
    n = np.arange(Nmax)
    M = np.zeros((Nmax, Nmax))
    pref = lam * w0**2 / (2*np.pi*T)**2
    for i in range(Nmax):
        row = np.zeros(Nmax)
        ni = n[i]
        # positive-branch n' = m >= 0, m != ni
        for m in range(Nmax):
            if m == ni:
                continue
            num = 1.0  # coefficient of Delta(m) from n'=m term
            denom = abs(2*m+1) * (ni-m)**2
            row[m] += pref * num / denom
            # diagonal contribution -Delta(ni)*(2m+1)/(2ni+1) / |2m+1| /(ni-m)^2
            row[ni] += pref * (-(2*m+1)/(2*ni+1)) / abs(2*m+1) / (ni-m)**2
        # negative-branch n' = -m-1, m>=0 (all m, since n'=-m-1 != ni always for ni>=0,m>=0)
        for m in range(Nmax):
            denom = abs(2*m+1) * (ni+m+1)**2
            row[m] += pref * 1.0 / denom
            row[ni] += pref * ((2*m+1)/(2*ni+1)) / abs(2*m+1) / (ni+m+1)**2
        M[i, :] = row
    return M

def max_eig(lam, w0, T, Nmax):
    M = build_matrix(lam, w0, T, Nmax)
    ev = np.linalg.eigvals(M)
    return np.max(ev.real)

def solve_Tc(lam, w0, Nmax, Tlo=1e-4, Thi=1e4, tol=1e-10):
    # eigenvalue is monotonically decreasing in T; find T where max_eig=1
    lo, hi = Tlo, Thi
    flo = max_eig(lam, w0, lo, Nmax) - 1.0
    fhi = max_eig(lam, w0, hi, Nmax) - 1.0
    assert flo > 0 and fhi < 0, (flo, fhi, "bracket failed")
    for _ in range(80):
        mid = np.sqrt(lo*hi)  # log-bisection since Tc ~ sqrt(lam) w0, wide range
        fm = max_eig(lam, w0, mid, Nmax) - 1.0
        if fm > 0:
            lo = mid
        else:
            hi = mid
        if hi/lo - 1 < tol:
            break
    return np.sqrt(lo*hi)

w0 = 1.0
results = []
for lam in [5, 10, 20, 40, 80]:
    Ns = [40, 60, 90]
    Tcs = []
    for Nmax in Ns:
        Tc = solve_Tc(lam, w0, Nmax)
        Tcs.append(Tc)
    # Richardson extrapolate assuming truncation error ~ 1/Nmax
    # simple 3-point Richardson using last two points (finer pair)
    T2, T3 = Tcs[-2], Tcs[-1]
    N2, N3 = Ns[-2], Ns[-1]
    # error ~ A/N -> extrapolate to N=inf: Tinf = T3 + (T3-T2)*N2/(N3-N2)
    Tinf = T3 + (T3 - T2) * (N2/(N3-N2))
    coeff = Tinf / (w0*np.sqrt(lam))
    results.append((lam, Tcs, Tinf, coeff))
    print(f"lam={lam:6.1f}  Tc(N=40,60,90)={Tcs}  Tc_extrap={Tinf:.6f}  Tc/(w0*sqrt(lam))={coeff:.6f}")

print()
print("Reference: Hofmann et al. quote Tc = 0.1827 w0 sqrt(lam) for 1<<lam<<N (their Eq. 33, citing Allen-Dynes 1975).")
print("Repo-certified C_inf = 0.1827262477 (asymptotic Allen-Dynes strong-coupling constant).")
