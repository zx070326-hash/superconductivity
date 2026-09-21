# Main-session INDEPENDENT recomputation of lambda / w_log / wbar2 from the raw EPW .a2f.
# Own parser, own quadrature. Does NOT import the repo Spectrum class.
from pathlib import Path
import numpy as np

D = Path("theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic")
FILES = {("H3S","anh"): D/"H3S-200GPa/EPW_anharmonic/H3S.a2f",
         ("D3S","anh"): D/"D3S-200GPa/EPW_anharmonic/D3S.a2f",
         ("H3S","har"): D/"H3S-200GPa/EPW_harmonic/H3S.a2f",
         ("D3S","har"): D/"D3S-200GPa/EPW_harmonic/D3S.a2f"}

def parse(p):
    w,a=[],[]
    for ln in p.read_text().splitlines()[1:]:
        if "Integrated" in ln: break
        t=ln.split()
        if len(t)<2: continue
        try: x,y=float(t[0]),float(t[1])
        except ValueError: continue
        w.append(x); a.append(y)
    return np.array(w), np.array(a)

def simpson(y,x):
    n=len(x)-1
    if n%2: return np.trapezoid(y,x)          # fall back if odd #intervals
    h=np.diff(x)
    if not np.allclose(h,h[0],rtol=1e-9): return np.trapezoid(y,x)
    return h[0]/3*(y[0]+y[-1]+4*y[1:-1:2].sum()+2*y[2:-1:2].sum())

def moments(w,a,quad):
    m = w>0                                    # a2F/w: drop the w=0 node (a2F(0)=0 anyway)
    ww,aa = w[m],a[m]
    lam   = 2*quad(aa/ww, ww)
    wlog  = np.exp((2/lam)*quad(np.log(ww)*aa/ww, ww))
    wbar2 = np.sqrt((2/lam)*quad(ww*aa, ww))
    return lam,wlog,wbar2

print("grid check + two quadrature schemes + grid halving")
for k,p in FILES.items():
    w,a=parse(p)
    dw=np.diff(w)
    print(f"\n{k[0]}-{k[1]}: n={len(w)}  w in [{w[0]:.4f},{w[-1]:.3f}] meV  "
          f"uniform={np.allclose(dw,dw[0],rtol=1e-9)}  dw={dw[0]:.6f}  a2F(0)={a[0]:.3e}")
    for nm,q in (("trapz",np.trapezoid),("simpson",simpson)):
        lam,wl,w2=moments(w,a,q)
        print(f"   {nm:8s} lambda={lam:.6f}  w_log={wl:.6f}  wbar2={w2:.6f}  w_log/wbar2={wl/w2:.6f}")
    lam,wl,w2=moments(w[::2],a[::2],np.trapezoid)      # halve the grid: discretisation probe
    print(f"   {'trapz/2':8s} lambda={lam:.6f}  w_log={wl:.6f}  wbar2={w2:.6f}   <- grid halved")

print("\n=== claims under test ===")
w,a=parse(FILES[("H3S","anh")]); lam,wl,w2=moments(w,a,np.trapezoid)
print(f"  agent: lambda=1.720958  w_log=109.995646 | mine: {lam:.6f} {wl:.6f}")
print(f"  csv:3  lambda=1.721     w_log=109.995  wbar2=149.26 | mine wbar2={w2:.6f}")
w,a=parse(FILES[("D3S","anh")]); lamD,wlD,w2D=moments(w,a,np.trapezoid)
print(f"  D3S anh: lambda={lamD:.6f} w_log={wlD:.6f} wbar2={w2D:.6f}")
print(f"  repo wbar2(D3S)=105.43 | mine={w2D:.4f}   ratio wbar2 D/H = {w2D/w2:.6f}  (1/sqrt2={2**-0.5:.6f})")
print(f"  ratio w_log D/H (Mishra a2F) = {wlD/wl:.6f}   <- compare Errea 73.3/92.86=0.789361")
print(f"  lambda D/H = {lamD/lam:.6f}  (uniform mass scaling requires exactly 1)")
