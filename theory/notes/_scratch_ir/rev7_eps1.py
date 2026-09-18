# rev7: eps1 = q(ln(1+Vw) + V/(V+D)) vs the EXACT sum of the per-n (L3) bounds  sum_n w_n eps_f(s_n)  (G<=1 dropped; closed forms via psi / Hurwitz zeta),
# decomposed into: T1 = Lam_f(D/2) * sum w_n/s_n   (NOT present in the eps1 formula),  T2 = q sum_{n<=n1} n/(n+1/2)^2,  T3 = q (V/2D) sum_{n>n1} (n+1/2)^-2
import mpmath as mp
mp.mp.dps = 30
VF = {"phi": mp.mpf("0.96059156149540"), "kap": mp.mpf("1.09092291179862")}
LAM = {"phi": lambda s: (mp.atan(s) + s*(s*s-1)/(1+s*s)**2)/8, "kap": lambda s: (mp.atan(s) - s/(1+s*s))/2}
print(" w | f | T1/q  T2/q  T3/q | (T1+T2+T3)/eps1 | (T2+T3)/eps1 | ln(1+Vw)+V/(V+D)")
for w in (6, 8, mp.mpf("8.4761"), 10, 45, 1000, 1e4, 1e6, 1e9, 1e12):
    w = mp.mpf(w); D = 1/w; q = D*D/8
    for nm in ("phi", "kap"):
        V = VF[nm]; n1 = int(mp.floor(V/(2*D)))          # 2 n D <= V  <=>  n <= V/(2D)
        T1 = LAM[nm](D/2) * (mp.pi**2/(4*D))            # sum_n w_n/s_n = (D/2) sum s_n^-2 = pi^2/(4D)
        # T2 = q * sum_{n=1}^{n1} n/(n+1/2)^2 = q [ (psi(n1+3/2)-psi(3/2)) - (1/2)(zeta(2,3/2)-zeta(2,n1+3/2)) ]
        T2 = q*((mp.digamma(n1+mp.mpf(3)/2) - mp.digamma(mp.mpf(3)/2)) - (mp.zeta(2, mp.mpf(3)/2) - mp.zeta(2, n1+mp.mpf(3)/2))/2)
        T3 = q*(V/(2*D))*mp.zeta(2, n1+mp.mpf(3)/2)
        e1 = q*(mp.log(1+V*w) + V/(V+D))
        print(f" {float(w):8g} {nm} | {float(T1/q):.4f} {float(T2/q):.4f} {float(T3/q):.4f} | {float((T1+T2+T3)/e1):.4f} | {float((T2+T3)/e1):.4f} | {float(e1/q):.4f}")
print("\n(T1/q -> pi^2/12 = %.4f as D->0;  the eps1 bracket ln(1+Vw)+V/(V+D) has no such term)" % float(mp.pi**2/12))
