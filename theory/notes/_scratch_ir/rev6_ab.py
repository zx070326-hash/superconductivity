# rev6: (a) A,B closed forms, elementary bound, and the "(1/24) int|Xi''|" midpoint claim
import mpmath as mp
mp.mp.dps = 40
gam = mp.euler
def AB_direct(w):
    x = mp.mpf(w)**2
    A = mp.nsum(lambda n: (1/(2*n+1)) * x/(x+(n+mp.mpf(1)/2)**2), [0, mp.inf])
    B = mp.nsum(lambda n: (1/(2*n+1)) * (x/(x+(n+mp.mpf(1)/2)**2))**2, [0, mp.inf])
    return A, B
def AB_closed(w):
    z = mp.mpc(mp.mpf(1)/2, w)
    A = (mp.re(mp.digamma(z)) + gam + 2*mp.log(2))/2
    B = A + mp.mpf(w)/4*mp.im(mp.polygamma(1, z))
    return A, B
def elem(w):
    w = mp.mpf(w)
    return mp.log(w) + gam + mp.log(4) - mp.mpf(3)/2 - 1/(12*w) - 1/(4*w**2) - 1/(16*w**2)
def elem_corrected(w):   # with the Peano-valid constant (1/8) int|Xi''| = 3/(16 w^2)
    w = mp.mpf(w)
    return mp.log(w) + gam + mp.log(4) - mp.mpf(3)/2 - 1/(12*w) - 1/(4*w**2) - 3/(16*w**2)
print(" w      | A_dir-A_cl   B_dir-B_cl | 2B-1            elem bound      2B-1-elem   | w^2*(SumXi-1/2)  | Binet remainder*12w")
for w in [8, 20, 45, 1000]:
    Ad, Bd = AB_direct(w); Ac, Bc = AB_closed(w)
    x = mp.mpf(w)**2
    SXi = mp.nsum(lambda n: (n+mp.mpf(1)/2)*x/((n+mp.mpf(1)/2)**2+x)**2, [0, mp.inf])
    z = mp.mpc(mp.mpf(1)/2, w)
    rem = mp.re(mp.digamma(z)) - (mp.log(x+mp.mpf(1)/4)/2 - 1/(4*x+1))
    print(f"{w:6d} | {mp.nstr(Ad-Ac,3):>10} {mp.nstr(Bd-Bc,3):>10} | {mp.nstr(2*Bc-1,14)} {mp.nstr(elem(w),14)} {mp.nstr(2*Bc-1-elem(w),4):>10} | {mp.nstr(x*(SXi-mp.mpf(1)/2),10):>14} | {mp.nstr(rem*12*w,6)}")
print("claimed midpoint majorant (1/24)V_eta = 1/16 =", 1/16., "; Peano-valid (1/8)V_eta = 3/16 =", 3/16., "; EM leading -f'(0)/24 -> 1/24 =", 1/24.)
# Is 2B-1 >= elem for ALL w >= 1? scan
bad = [w for w in [0.5,1,2,3,5,7.78,10,15,30,100,300,3000] if 2*AB_closed(w)[1]-1 < elem(w)]
print("elem bound violated at:", bad)
# the sign of E = int Xi - Sum Xi(n+1/2): here Sum < int (S positive means Sum-1/2>0)
