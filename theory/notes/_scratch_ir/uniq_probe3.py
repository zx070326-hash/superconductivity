# (S2) exhaustive: Q^{(c)}_{LL'}(delta_k) >= 0 for all k<=K, L,L'<=Lmax, weights c=1 and c=u.
# (S3) stress: Perron vectors of O[F] in the cone for adversarial P (20-atom log-uniform, extreme widths) and
#      for the Coulomb pencil M = H - mt*Kc (does the cone survive mu*?).
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "derivations"))
from p4_interchange_uniformity import H_mat
Lmax, K = 400, 400
n = np.arange(Lmax); u = 1/np.sqrt(2*n+1)
worst = {1: (np.inf, None), 'u': (np.inf, None)}
for k in range(1, K+1):
    M = np.zeros((Lmax, Lmax))
    # L2 off-diagonal |n-m|=k ; L3 n+m+1=k ; -L1 diagonal 2/(2n+1) for n>=k
    i = np.arange(Lmax-k); M[i, i+k] += u[i]*u[i+k]; M[i+k, i] += u[i]*u[i+k]
    if k-1 < Lmax:
        for nn in range(0, k):
            mm = k-1-nn
            if mm < Lmax: M[nn, mm] += u[nn]*u[mm]
    d = np.arange(k, Lmax); M[d, d] -= 2/(2*d+1)
    for tag, c in ((1, np.ones(Lmax)), ('u', u)):
        W = M*np.outer(c, c); Q = np.cumsum(np.cumsum(W, 0), 1)
        mn = Q.min()
        if mn < worst[tag][0]: worst[tag] = (mn, (k,)+tuple(int(x)+1 for x in np.unravel_index(Q.argmin(), Q.shape)))
print("(S2) min over k<=%d, L,L'<=%d of Q_{LL'}(delta_k):" % (K, Lmax))
for tag in (1, 'u'): print(f"   c={tag}: min = {worst[tag][0]:+.3e} at (k,L,L') = {worst[tag][1]}")
# closed form check at L=L'=infinity for c=u: sum over line1 2*sum_n 1/((2n+1)(2n+2k+1)) + line2 - 2*sum_{n>=k}1/(2n+1)^2
N = 800; nn = np.arange(N)
print("(S3) stress cone membership (v nonincreasing / v*sqrt(2n+1) nonincreasing):")
rng = np.random.default_rng(5); vio1 = viou = 0; T = 60
for t in range(T):
    k = 20; ws = np.exp(rng.uniform(np.log(0.02), np.log(200), k)); p = rng.dirichlet(0.3*np.ones(k))
    Kmat = sum(pi*H_mat(N, wi) for pi, wi in zip(p, ws)); e, V = np.linalg.eigh(Kmat); v = V[:, -1]*np.sign(V[0,-1])
    if v.min() < -1e-12 or np.any(np.diff(v) > 1e-13): vio1 += 1
    if np.any(np.diff(v*np.sqrt(2*nn+1)) > 1e-13): viou += 1
print(f"   20-atom log-uniform P (0.02..200), {T} trials: violations v-cone {vio1}, u-cone {viou}")
def Kc(Nn, Nc):
    z = np.zeros(Nn); kk = min(Nc, Nn); z[:kk] = 1/np.sqrt(2*np.arange(kk)+1); return 2*np.outer(z, z)
vio1 = viou = 0; cnt = 0
for w in (0.7, 1.0, 2.0, 5.0, 20.0):
    for mt in (0.1, 0.26, 0.33):
        for Nc in (1, 3, 10, 40, 200):
            e, V = np.linalg.eigh(H_mat(N, w) - mt*Kc(N, Nc)); v = V[:, -1]*np.sign(V[0,-1]); cnt += 1
            if v.min() < -1e-12 or np.any(np.diff(v) > 1e-13): vio1 += 1
            if np.any(np.diff(v*np.sqrt(2*nn+1)) > 1e-13): viou += 1
print(f"   Coulomb pencil H - mt*Kc, {cnt} cases: violations v-cone {vio1}, u-cone {viou}")
# (S1) sign of the Coulomb-pencil derivative: d/dT lam1(H(w)-mt Kc) = -(1/T)<v, w H'(w) v>  (Kc has no T-dependence at fixed Nc)
