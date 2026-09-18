import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from lemmaC_compact_probe import O_of, top, THETA_STAR

N=600

def build(w, N=N):
    x = w*w
    f = lambda k: x/(k*k+x)
    psi = lambda k: 1.0/(k*k+x)**2
    A = O_of(f,N); Opsi = O_of(psi,N)
    h,v,l2,ev,EV = top(A)
    g = h-l2
    npsi = np.abs(np.linalg.eigvalsh(Opsi)).max()
    return dict(w=w,x=x,N=N,A=A,Opsi=Opsi,h=h,v=v,l2=l2,g=g,npsi=npsi)

def gap_an(w,h): return 0.4*h if w<=2 else h-THETA_STAR

def trial_vector(pr, alpha, c=0.5):
    n=np.arange(pr["N"]); u=1.0/np.sqrt(2.0*n+1.0)
    a=alpha*pr["w"]; t=u*a*a/(a*a+(n+c)**2)
    return t/np.linalg.norm(t)

def sinTheta_bound(pr,t,g_an):
    A=pr["A"]
    rho=t@(A@t); res=np.linalg.norm(A@t-rho*t)
    return min(res/g_an,1.0) if g_an>0 else np.inf

CELLS = [(0.30,0.42,5.0),(0.60,0.85,2.7),(0.85,1.20,2.1),(1.20,1.70,1.8),(2.40,3.40,1.5),(13.6,20.0,1.24)]
print(f"{'cell':>15} {'a':>5} {'w_worst':>8} {'S(g_an)':>9} {'S(g_true)':>10} {'sinTh_true':>10} {'S(sinTh_true)':>14}")
for wa,wb,alpha in CELLS:
    ws = np.geomspace(wa,wb,11)
    worst=None
    for w in ws:
        pr = build(w)
        g_an = gap_an(w,pr["h"])
        t = trial_vector(pr,alpha)
        s_an = sinTheta_bound(pr,t,g_an)
        # true sinTheta: angle between t and true top eigenvector v
        cos = abs(t@pr["v"])
        sin_true = np.sqrt(max(1-cos*cos,0.0))
        g_true = pr["g"]
        npsi = pr["npsi"]
        Opsi_t = pr["Opsi"]@t
        qpsi_t = float(t@Opsi_t)
        e_psi = np.linalg.norm(Opsi_t-qpsi_t*t)
        Kt0 = e_psi
        def Sval(s,g):
            if s<=0: return np.inf
            return np.sqrt(max(1-s*s,0))*Kt0/(2*s*npsi)
        S_gan = Sval(s_an,g_an)
        S_gtrue = Sval(s_an,g_true)
        S_sinTrue = Sval(sin_true,g_true)
        if worst is None or S_gan < worst[1]:
            worst = (w,S_gan,S_gtrue,sin_true,S_sinTrue)
    w_,S1,S2,st,S3 = worst
    print(f"[{wa:5.2f},{wb:5.2f}] {alpha:5.2f} {w_:8.3f} {S1:9.4f} {S2:10.4f} {st:10.4f} {S3:14.4f}")
