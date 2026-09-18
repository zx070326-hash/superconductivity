"""Is a gamma>2 power-law kernel even inside Theorem E's hypothesis class?
Class P := { f(k) = int omega^2/(omega^2+k^2) dP(omega), P >= 0 finite }.
Claim: f in P, P != 0  =>  k^2 f(k) -> int omega^2 dP  (monotone convergence),
so f decays EXACTLY like k^-2 (finite 2nd moment) or SLOWER (infinite).
Never faster. Hence gamma>2 is disjoint from P.
"""
import numpy as np
def f_of(meas, k):
    k = np.asarray(k, float); o = np.zeros_like(k)
    for w, om in meas: o += w*om**2/(om**2+k**2)
    return o
ks = np.array([1e1, 1e2, 1e3, 1e4, 1e5, 1e6])
tests = {
  "single Einstein w=1,om=3":      [(1.0, 3.0)],
  "4-peak broad":                  [(0.7,0.3),(1.3,2.0),(0.5,11.0),(2.0,40.0)],
  "heavy tail (approx inf 2nd mom)":[(1.0/om**2.5, om) for om in np.logspace(0,6,400)],
}
for name, m in tests.items():
    m2 = sum(w*om**2 for w, om in m)
    print(f"\n{name}:  int omega^2 dP = {m2:.6g}")
    print("   k        k^2 f(k)        f(k)*k^3 (would ->const if gamma=3)")
    for k in ks:
        fk = f_of(m, np.array([k]))[0]
        print(f"   {k:8.0e}  {k**2*fk:14.6g}  {k**3*fk:16.6g}")
print("\n=> k^2 f(k) saturates at the 2nd moment; k^3 f(k) DIVERGES like k.")
print("=> no positive measure reproduces a k^-3 tail. gamma>2 is OUTSIDE class P.")
