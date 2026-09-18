import mpmath as mp
mp.mp.dps = 30

# Build KAY operator from scratch for power-law kernel f(k) = (varpi)^gamma / (k^2+varpi^2)^(gamma/2)
# i.e. generalize the Einstein/Lorentzian kernel llbracket k rrbracket(varpi)=varpi^2/(k^2+varpi^2)
# to gamma-power: [[k]](varpi) = (varpi^2/(k^2+varpi^2))^(gamma/2)  -- test H(x)=h(sqrt(x)) concavity near x->0 for gamma>2
N = 150
u = [ (2*n+1)**mp.mpf('-0.5') for n in range(N) ]

def build_O(f_vals):
    # f_vals[k] for k=0..2N (need f at |n-m| and n+m+1)
    M = mp.zeros(N,N)
    for n in range(N):
        # L1 diagonal: (2/(2n+1)) * sum_{k=1}^n f(k)
        s = mp.mpf(0)
        for k in range(1,n+1):
            s += f_vals[k]
        M[n,n] = -s*2/(2*n+1)
    for n in range(N):
        for m in range(N):
            if n!=m:
                M[n,m] += f_vals[abs(n-m)]*u[n]*u[m]
            M[n,m] += f_vals[n+m+1]*u[n]*u[m]
    return M

def kernel_gamma(k,varpi,gamma):
    return (varpi**2/(k**2+varpi**2))**(mp.mpf(gamma)/2)

def h_of_varpi(varpi,gamma):
    f_vals = [kernel_gamma(k,varpi,gamma) for k in range(2*N+1)]
    M = build_O(f_vals)
    ev = mp.eigsy(M, eigvals_only=True)
    return max(ev)

def H_of_x(x,gamma):
    varpi = mp.sqrt(x)
    return h_of_varpi(varpi,gamma)

print("=== Check 3: H(x)=h(sqrt(x)) concavity near x->0 for power-law kernel, gamma=3.0 ===")
gamma = 3.0
xs = [mp.mpf('0.02'), mp.mpf('0.05'), mp.mpf('0.1'), mp.mpf('0.2'), mp.mpf('0.35'), mp.mpf('0.5')]
hh = mp.mpf('0.005')
for x in xs:
    Hm = H_of_x(x-hh,gamma)
    H0 = H_of_x(x,gamma)
    Hp = H_of_x(x+hh,gamma)
    H2 = (Hp-2*H0+Hm)/hh**2
    print(f"x={float(x):.3f}  H''(x) approx = {float(H2):.6f}  {'VIOLATION(convex)' if H2>0 else 'concave ok'}")
