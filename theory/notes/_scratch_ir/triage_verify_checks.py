import mpmath as mp
mp.mp.dps = 50

# Check 1: beta_function_rg_flow_section X.2 digamma asymptotic
def S(N):
    return mp.mpf('0.5')*(mp.digamma(N+mp.mpf('0.5')) - mp.digamma(mp.mpf('0.5')))

const = mp.euler + mp.log(4)
print("=== Check 1: S(N) ~ 0.5 ln N + 0.5(gamma+ln4) ===")
for N in [10,50,200,800]:
    val = S(N)
    resid = 2*val - mp.log(N) - const
    print("N=",N,"2S-lnN-const=",resid)
print("gamma+ln4 =", const)
print()

# Check 2: explore_independent_measures Z-channel convexity: f_k(x)=x/(k^2+x) concave check
print("=== Check 2: f_k(x)=x/(k^2+x) concavity (f''<0 for x>0) ===")
def f(k,x):
    return x/(k**2+x)
h = mp.mpf('1e-6')
for k in [1,2,5,20]:
    for x in [0.01, 1, 10, 100]:
        k_=mp.mpf(k); x_=mp.mpf(x)
        f2 = (f(k_,x_+h)-2*f(k_,x_)+f(k_,x_-h))/h**2
        analytic = -2*k_**2/(k_**2+x_)**3
        print(f"k={k} x={x} f''(numeric)={float(f2):.6e} f''(analytic)={float(analytic):.6e}")
