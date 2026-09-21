"""Route R' feasibility: M1's Lorentz trial with freed width, tau_n = u_n a^2/(a^2+(n+1/2)^2),
a = alpha*w with alpha piecewise-constant on cells; DK ratio with analytic gap formula
(0.4*rho for w<=2, rho-Theta* for w>=2). Trial quantities in float64 (numerical-evidence)."""
import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
from lemmaC_compact_probe import O_of, top, THETA_STAR
N=600; n=np.arange(N); u=1/np.sqrt(2*n+1)
def ratio(w, alpha):
    x=w*w; a=alpha*w
    t=u*a*a/(a*a+(n+0.5)**2); t/=np.linalg.norm(t)
    A=O_of(lambda k: x/(k*k+x),N)
    rho=t@A@t; res=np.linalg.norm(A@t-rho*t)
    glo=0.4*rho if w<=2 else rho-THETA_STAR
    if glo<=0: return np.inf,res,np.inf,rho
    s=min(res/glo,1.0)
    Op=O_of(lambda k: k*k/(k*k+x)**3,N); Os=O_of(lambda k: 1/(k*k+x)**2,N)
    Ft=t@Op@t; e_phi=np.linalg.norm(Op@t-Ft*t); Kt=Os@t; q=t@Kt; e_psi=np.linalg.norm(Kt-q*t)
    nphi=np.abs(np.linalg.eigvalsh(Op)).max(); npsi=np.abs(np.linalg.eigvalsh(Os)).max()
    F_lo=(1-s*s)*Ft-2*s*e_phi-s*s*nphi; E_hi=e_psi+s*(abs(q)+npsi)
    return (x*x*E_hi**2/(glo*F_lo) if F_lo>0 else np.inf), res, s, rho
cells=[(0.30,0.42,5.0),(0.42,0.60,3.6),(0.60,0.85,2.7),(0.85,1.2,2.1),(1.2,1.7,1.8),(1.7,2.4,1.6),(2.4,3.4,1.5),(3.4,4.8,1.42),(4.8,6.8,1.35),(6.8,9.6,1.30),(9.6,13.6,1.27),(13.6,20,1.24)]
print(f"{'w-cell':>14} {'alpha':>6} | worst ratio   worst res   worst s_bnd   (over 9 sample points)")
for wa,wb,al in cells:
    ws=np.exp(np.linspace(np.log(wa),np.log(wb),9))
    R=[ratio(w,al) for w in ws]
    print(f"[{wa:5.2f},{wb:5.2f}] {al:6.2f} | {max(r[0] for r in R):11.3f} {max(r[1] for r in R):11.4f} {max(r[2] for r in R):12.4f}")
