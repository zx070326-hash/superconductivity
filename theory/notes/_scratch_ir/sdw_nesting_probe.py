import numpy as np
from scipy import integrate

# Toy nesting susceptibility: two bands xi_e(k)=eps, xi_h(k+Q) = -eps + delta
# chi_0(Q,T) = N0 * integral_{-W}^{W} [f(eps) - f(delta-eps)] / (delta - 2 eps) d eps
# f = Fermi function. N0 set to 1 (per spin per unit energy, flat DOS approx).

def fermi(x, T):
    # numerically stable
    return 1.0/(1.0+np.exp(np.clip(x/T, -700, 700)))

def chi0(T, delta, W=50.0, N0=1.0):
    def integrand(eps):
        denom = delta - 2*eps
        if abs(denom) < 1e-12:
            # l'Hopital-ish: derivative of f at that point ~ -f'(eps)/... just use small shift
            denom = 1e-12
        return N0*(fermi(eps,T) - fermi(delta-eps,T))/denom
    val, err = integrate.quad(integrand, -W, W, limit=400, points=[delta/2.0])
    return val

W = 50.0
print("=== perfect nesting delta=0: expect chi0 ~ N0 ln(1.13 W/T) (Cooper-log form) ===")
for T in [1.0, 0.3, 0.1, 0.03, 0.01]:
    c = chi0(T, 0.0, W=W)
    pred = np.log(1.13*W/T)
    print(f"T={T:8.4f}  chi0={c:9.5f}  ln(1.13W/T)={pred:9.5f}  ratio={c/pred:.4f}")

print()
print("=== finite mismatch delta: log should SATURATE as T->0 (no divergence) ===")
for delta in [0.0, 0.5, 2.0, 5.0, 10.0]:
    print(f" delta={delta}")
    for T in [1.0, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001]:
        c = chi0(T, delta, W=W)
        print(f"   T={T:8.4f}  chi0={c:9.5f}")
