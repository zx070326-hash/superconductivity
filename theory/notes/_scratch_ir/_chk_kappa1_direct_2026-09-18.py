import mpmath as mp
mp.mp.dps = 30


def h1_kernel(x):
    x = mp.mpf(x)
    return x**2/(x**2+1)


u0 = mp.mpf('1.64057683931113599607318111086')
wbar0 = 1/u0
print("wbar0=", wbar0)

for sv in ['0.001', '0.002', '0.005', '0.01', '0.02', '0.05']:
    s = mp.mpf(sv)
    w1, w2 = 1-s, 1+s
    avg = (h1_kernel(w1/u0) + h1_kernel(w2/u0))/2
    h0 = h1_kernel(wbar0)
    rel = (avg-h0)/h0
    coeff = rel/s**2
    print(f"s={float(s):<8} avg_h1={float(avg):.12f}  rel_dev/s^2={float(coeff):.6f}")

# analytic kappa1(wbar0)
kappa1 = (2-6*wbar0**2)/(wbar0**2+1)**2
print("kappa1/2 predicted =", kappa1/2)
