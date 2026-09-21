from mpmath import mp, mpf
import mpmath
mp.dps = 60
def N(s): return 2*s + (1-s*s)*mp.atan(s) - s*mp.log(1+s*s)
def Th(s): return (mp.log(1+s*s) + 2*s*mp.atan(s))/(2*(1+s*s))
a,b = mpf(1),mpf(2)
for _ in range(200):
    m=(a+b)/2
    if N(m)>0: a=m
    else: b=m
star=(a+b)/2
true_theta = Th(star)
print("true Theta* (dps60)   =", mp.nstr(true_theta, 30))
# replicate the shipped routine at dps 40
with mp.workdps(40):
    a2,b2=mpf(1),mpf(2)
    for _ in range(150):
        m=(a2+b2)/2
        if N(m)>0: a2=m
        else: b2=m
    val = max(Th(a2),Th(b2)) + (b2-a2)
    margin = float(b2-a2); v = float(val)
print("shipped lo,hi          =", repr(v-margin), repr(v+margin))
print("margin (half-width)    =", margin)
lo,hi = mpf(v)-mpf(margin), mpf(v)+mpf(margin)
print("true >= lo ?", true_theta >= lo, "  true <= hi ?", true_theta <= hi)
print("ENCLOSES?", bool(lo <= true_theta <= hi))
print("excursion true-hi      =", mp.nstr(true_theta-hi, 10))
print("float64 rounding err   =", mp.nstr(mpf(v)-val, 10))
